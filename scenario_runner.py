import os
import json
import numpy as np
import pandas as pd 
from scipy.integrate import solve_ivp
from scipy.sparse import lil_matrix
from scipy.interpolate import interp1d
import sys
import re
import warnings
import time

# Internal Module Imports
from mechanism_solver import MechanismSolver
from easy_parser import EasyParser
from mechanism import Mechanism, Reaction
from load_saphir_environment import load_saphir_environment
import mcm_kinetics

class CallableFloat(float):
    """A bulletproof float that returns 0.0 if accidentally called as a function."""
    def __call__(self, *args, **kwargs):
        return CallableFloat(0.0)

class SafeDict(dict):
    """A bulletproof dictionary for chemical math. If a variable is missing, it safely returns CallableFloat(0.0)."""
    def __missing__(self, key):
        return CallableFloat(0.0)

class ScenarioRunner:
    def __init__(self, config_path):
        print(f">>> [DEBUG] Initializing ScenarioRunner with: {config_path}")
        if not os.path.exists(config_path):
            raise FileNotFoundError(f"Config file not found: {config_path}")
            
        with open(config_path, 'r') as f:
            self.cfg = json.load(f)
            
        self.root = self.cfg['folder']
        self.mechanism = None
        self.solver = None
        self.t_span = None
        self.ctx_data = {}
        self.var_assignments = []
        self.last_progress_print_time = -1
        self.t_start = 0
        self.t_end = 0
        self.eval_count = 0
        self.real_start_time = 0.0

    def run(self):
        try:
            _run_start = time.time()
            print("\n>>> Starting SAPHIR Simulation Pipeline")
            print("-" * 40)
            self._setup_environment()
            self._build_mechanism()
            self._solve()
            out_path = self._export_results()
            _elapsed = time.time() - _run_start
            print("-" * 40)
            print(f">>> All steps completed successfully.  Total time: {_elapsed//60:.0f}m {_elapsed%60:.1f}s")
            self._show_results(out_path)
        except Exception as e:
            print(f"\n[ERROR] [CRITICAL ERROR] Simulation failed!")
            print(f"Error Details: {str(e)}")
            import traceback
            traceback.print_exc()

    def _setup_environment(self):
        enz_file = next(f for f in os.listdir(self.root) if 'kammer' in f.lower() and f.endswith('.enz'))
        enz_path = os.path.join(self.root, enz_file)
        providers, t_min, t_max, headers = load_saphir_environment(enz_path)
        
        self.t_span = np.linspace(t_min, t_max, num=2000)

        for raw_header, interp_func in providers.items():
            self.ctx_data[raw_header.upper()] = interp_func

        # The ENZ 'dilute' column is raw flow in L/min. Mechanism .dil files use
        # CONST(DILUTE) as a rate constant (s⁻¹), which would be ~10⁶× too large.
        # The runner already applies correct dilution analytically (dydt -= y * k_dil),
        # so disable the legacy mechanism dil reactions by setting DILUTE = 0.
        self.ctx_data['DILUTE'] = lambda t: 0.0

        for internal_key, mapping in self.cfg['env_map'].items():
            if mapping['type'] == "SAPHIR File":
                file_header = mapping['value'].lower()
                if file_header not in providers:
                    actual_header = next((h for h in headers if h.lower() == file_header), None)
                    if actual_header: file_header = actual_header
                if file_header in providers:
                    raw_interp = providers[file_header]
                    # For H2O: determine unit scaling once from the full-timeseries maximum,
                    # then build a pre-corrected interpolation so the per-timestep code never
                    # sees the ambiguous 0.02–1.0 grey zone.
                    if internal_key.upper() == 'H2O':
                        t_arr = np.linspace(t_min, t_max, 500)
                        raw_vals = raw_interp(t_arr)
                        peak = float(np.nanmax(np.abs(raw_vals)))
                        if peak > 1.0:
                            # Stored as percentage (e.g. 1.27 for 1.27%) — divide once here
                            scale = 1.0 / 100.0
                        else:
                            # Stored as VMR fraction (e.g. 0.0127)
                            scale = 1.0
                        corrected_vals = np.where(raw_vals < 0.02, 0.0, raw_vals * scale)
                        self.ctx_data[internal_key] = interp1d(t_arr, corrected_vals,
                                                               fill_value="extrapolate")
                        print(f"  [ENV] H2O unit auto-detected: {'percentage /100' if scale < 1 else 'VMR fraction'}, "
                              f"peak={peak:.4g}, noise floor < 0.02 zeroed.")
                    else:
                        self.ctx_data[internal_key] = raw_interp
            else:
                val = float(mapping['value'])
                self.ctx_data[internal_key] = lambda t, v=val: np.full_like(t, v) if isinstance(t, np.ndarray) else v

    def get_M(self, t):
        temp = self.ctx_data['TEMP'](t)
        press = self.ctx_data['PRESS'](t)
        return (press * 100) / (1.3806e-23 * temp) * 1e-6

    def _build_mechanism(self):
        chem_path = os.path.join(self.root, "chemistry")
        parser = EasyParser(base_path=chem_path)
        self.mechanism = Mechanism()
        
        # Sort files to ensure .rate definitions are parsed before .mech files
        mech_files = sorted(self.cfg['mechanism_files'], key=lambda x: ".rate" not in x)
        
        for f_name in mech_files:
            f_path = os.path.join(chem_path, f_name)
            if os.path.exists(f_path):
                m_part = parser.parse(f_path)
                self.mechanism.species.extend(m_part.species)
                self.mechanism.reactions.extend(m_part.reactions)

                # Only use uncommented parameters — the parser already filters invalid names
                for p_name, p_expr in m_part.parameters.items():
                    p_clean = p_expr.replace('@', '**')
                    p_clean = re.sub(r'(?<=\d)[dD](?=[+-]?\d)', 'E', p_clean)
                    p_clean = re.sub(r'J\((\d+)\)', r'J\1', p_clean)
                    p_clean = p_clean.strip()
                    try:
                        self.var_assignments.append((p_name.upper(), compile(p_clean, '<string>', 'eval')))
                    except Exception as e:
                        print(f"[WARNING] Could not compile parameter {p_name}: {e}")
                    
        # Apply blocked_parameters: explicit suppression from GUI Chemistry Overview tab
        blocked_params = {p.upper() for p in self.cfg.get('blocked_parameters', [])}
        if blocked_params:
            before = len(self.var_assignments)
            self.var_assignments = [(n, c) for n, c in self.var_assignments if n not in blocked_params]
            removed = before - len(self.var_assignments)
            if removed:
                print(f"  ... Blocked {removed} parameter(s): {', '.join(sorted(blocked_params))}")

        for j_target, j_expr in self.cfg.get('j_scalings', {}).items():
            j_clean = j_expr.upper().replace('@', '**')
            try: self.var_assignments.append((j_target.upper(), compile(j_clean, '<string>', 'eval')))
            except: pass
            
        all_species = set(self.mechanism.species)
        for spec in self.cfg.get('initials', {}).keys():
            all_species.add(spec.upper())
        for b in self.cfg.get('burst_injections', []):
            if b.get('target'): all_species.add(b['target'].upper())
        for s in self.cfg.get('continuous_sources', []):
            if s.get('target'): all_species.add(s['target'].upper())
            
        for r in self.cfg.get('custom_reactions', []):
            if not r.get('reactants') or not r.get('products') or not r.get('rate'): continue
            parts = re.split(r'\s*\+\s*|\s*->\s*|\s*=\s*', r['reactants'] + " + " + r['products'])
            for p in parts:
                clean_p = p.strip().upper()
                clean_p = re.sub(r'^[\d\.]+\s*\*\s*', '', clean_p) 
                if clean_p and not clean_p[0].isdigit(): all_species.add(clean_p)
            
            reactants = [x.strip().upper() for x in r['reactants'].split('+') if x.strip()]
            products = [x.strip().upper() for x in r['products'].split('+') if x.strip()]
            self.mechanism.reactions.append(Reaction(reactants, products, r['rate']))

        # --- TwoReaction injections: virtual intermediate species + first-order decay reactions ---
        # Each TwoReaction burst gets three internal species:
        #   _INJ1_<uid>  — fast buffer tank  (fraction fast_frac of total)
        #   _INJ2_<uid>  — slow buffer tank  (fraction 1-fast_frac)
        #   _INJA_<uid>  — delivery pool     → TARGET  at k_delivery
        # The source term (filling the tanks during the valve-open window) is
        # applied in _solve() as a burst_eval entry, NOT as a reaction.
        self._tr_specs = {}  # burst_idx → {'uid', 'sp1', 'sp2', 'spa', 'fast_frac'}
        for b_idx, b in enumerate(self.cfg.get('burst_injections', [])):
            if b.get('profile', 'Gaussian') != 'TwoReaction':
                continue
            if not b.get('target') or not b.get('ppb') or b['ppb'] == '':
                continue
            raw_uid = re.sub(r'[^A-Z0-9]', '_', b.get('name', f'INJ{b_idx}').upper())[:12].rstrip('_') or f'INJ{b_idx}'
            uid, suffix = raw_uid, 0
            while any(v['uid'] == uid for v in self._tr_specs.values()):
                suffix += 1; uid = f"{raw_uid[:10]}_{suffix}"
            sp1 = f'_INJ1_{uid}'
            sp2 = f'_INJ2_{uid}'
            spa = f'_INJA_{uid}'
            target_up = b['target'].upper()
            all_species.update([sp1, sp2, spa])
            try:
                k_fast     = float(eval(str(b.get('k_fast',     '0.01667'))))
                k_slow     = float(eval(str(b.get('k_slow',     '0.000556'))))
                k_delivery = float(eval(str(b.get('k_delivery', '0.01667'))))
                fast_frac  = float(eval(str(b.get('fast_frac',  '0.9'))))
            except Exception:
                k_fast, k_slow, k_delivery, fast_frac = 1/60.0, 1/1800.0, 1/60.0, 0.9
            self.mechanism.reactions.append(Reaction([sp1], [spa],       str(k_fast)))
            self.mechanism.reactions.append(Reaction([sp2], [spa],       str(k_slow)))
            self.mechanism.reactions.append(Reaction([spa], [target_up], str(k_delivery)))
            byp_t = b.get('byproduct_target', '').strip().upper()
            byp_f = 0.0
            if b.get('byproduct_frac'):
                try: byp_f = float(eval(str(b['byproduct_frac'])))
                except: pass
            if byp_t and byp_f > 0.0:
                all_species.add(byp_t)
                self.mechanism.reactions.append(Reaction([spa], [byp_t], str(k_delivery * byp_f)))
            self._tr_specs[b_idx] = {'uid': uid, 'sp1': sp1, 'sp2': sp2, 'spa': spa, 'fast_frac': fast_frac}
            print(f"  [TwoReaction] '{b.get('name')}' -> {sp1}/{sp2}/{spa} -> {target_up} "
                  f"(k_fast={k_fast:.4g}  k_slow={k_slow:.4g}  k_del={k_delivery:.4g}  "
                  f"f_fast={fast_frac:.2f})")

        # --- BUG FIX: Only purge explicitly mapped GUI variables, leave O3 alone! ---
        for mapped_key in self.cfg['env_map'].keys():
            all_species.discard(mapped_key.upper())
        for bg in ['RH', 'O2', 'N2', 'M']:
            all_species.discard(bg)

        self.mechanism.species = sorted(list(all_species))

        # Apply reaction exclusions / rate overrides defined in the scenario JSON
        excluded_rxns = self.cfg.get('excluded_reactions', [])
        if excluded_rxns:
            def _canon(species_list):
                result = set()
                for s in species_list:
                    name = s.split('*')[-1].strip().upper()
                    if name and not name.replace('.', '').isdigit():
                        result.add(name)
                return frozenset(result)

            exclude_keys = set()
            rate_overrides = []  # list of (r_canon, p_canon, new_rate_str)
            for e in excluded_rxns:
                r_canon = _canon([t.strip() for t in e['reactants'].split('+') if t.strip()])
                p_canon = _canon([t.strip() for t in e['products'].split('+') if t.strip()])
                if e.get('action', 'Remove') == 'Modify Rate' and e.get('rate', '').strip():
                    rate_overrides.append((r_canon, p_canon, e['rate'].strip()))
                else:
                    exclude_keys.add((r_canon, p_canon))

            before = len(self.mechanism.reactions)
            self.mechanism.reactions = [
                rxn for rxn in self.mechanism.reactions
                if (_canon(rxn.reactants), _canon(rxn.products)) not in exclude_keys
            ]
            removed = before - len(self.mechanism.reactions)
            if removed:
                print(f"  ... Excluded {removed} reaction(s) per scenario config.")

            # Patch rate expressions — safe_expr is compiled afterwards so this takes effect
            modified = 0
            for rxn in self.mechanism.reactions:
                rxn_key = (_canon(rxn.reactants), _canon(rxn.products))
                for r_canon, p_canon, new_rate in rate_overrides:
                    if rxn_key == (r_canon, p_canon):
                        rxn.rate_expression = new_rate
                        modified += 1
            if modified:
                print(f"  ... Modified rate for {modified} reaction(s) per scenario config.")

        self._add_artifacts()

    def _add_artifacts(self):
        if self.cfg['artifacts']['active']:
            scale = eval(str(self.cfg['artifacts']['scale']))
            if self.cfg['artifacts'].get('hono_active', True):
                self.mechanism.reactions.append(Reaction([], ["HONO"], f"({self.cfg['artifacts']['hono_eqn']}) * {scale}"))
            if self.cfg['artifacts'].get('hcho_active', True):
                self.mechanism.reactions.append(Reaction([], ["HCHO"], f"{self.cfg['artifacts']['hcho_eqn']}"))
            if self.cfg['artifacts'].get('y_active', True):
                y_rate = f"({self.cfg['artifacts']['y_eqn']}) * 1e-9 * M * 2.4e-13"
                self.mechanism.reactions.append(Reaction(['OH'], ["HO2"], y_rate))
            if self.cfg['artifacts'].get('wall_active', True):
                # dark_hono expression is a direct production rate in molec cm⁻³ s⁻¹.
                # H2O dependence (if any) must be explicit in the expression itself.
                self.mechanism.reactions.append(Reaction([], ["HONO"], str(self.cfg['artifacts']['dark_hono'])))

    def _solve(self):
        self.solver = MechanismSolver(self.mechanism, self.cfg['env_map'])
        y0 = np.zeros(self.solver.num_species)
        
        saphir_base_date = pd.Timestamp("2000-01-01 00:00:00")

        for spec, ppb in self.cfg.get('initials', {}).items():
            if not ppb: continue
            spec_up = spec.upper()
            if spec_up in self.solver.species_to_idx:
                y0[self.solver.species_to_idx[spec_up]] = (float(ppb) * 1e-9) * self.get_M(self.t_span[0])

        print("\n>>> [DIAGNOSTICS] Testing chemical mechanism evaluation...")
        
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            
            test_ctx = SafeDict()
            test_ctx.update({
                "T": 298.15, "TEMP": 298.15, "M": 2.5e19, "RH": 50.0, 
                "O2": 0.2095 * 2.5e19, "N2": 0.7809 * 2.5e19,
                "np": np, "exp": np.exp, "EXP": np.exp, 
                "LOG10": np.log10, "log10": np.log10,
                "LOG": np.log, "log": np.log,
                "COS": np.cos, "SIN": np.sin, "TAN": np.tan, "SQRT": np.sqrt
            })
            
            for func_name in dir(mcm_kinetics):
                if not func_name.startswith('_'):
                    func = getattr(mcm_kinetics, func_name)
                    if func_name.startswith('KMT'):
                        test_ctx[func_name.upper()] = lambda *args, f=func: f(*args, test_ctx['TEMP'], test_ctx['M'])
                    else:
                        test_ctx[func_name.upper()] = func  # COMPRC loaded as direct callable

            # Seed test PHTFKT callable using SZA=45 test conditions
            _test_cosx = max(np.cos(45.0 * 0.0174533), 0.01)
            test_ctx['COSX']     = _test_cosx
            test_ctx['LNO2']     = 1.165e-2
            test_ctx['MNO2']     = 0.244
            test_ctx['NNO2']     = -0.267
            test_ctx['PHTFKT']   = lambda a, b, c, _c=_test_cosx: a * _c**b * np.exp(c / _c)
            test_ctx['PHTFKT_CS']= test_ctx['PHTFKT']
            test_ctx['AMAX']     = np.maximum
            test_ctx['AMIN']     = np.minimum

            for k in self.ctx_data.keys(): test_ctx[k.upper()] = 1.0
            for s in self.solver.species: test_ctx[s] = 1.0

            failed_vars = []
            for _ in range(2):
                for var_name, comp_expr in self.var_assignments:
                    try:
                        test_ctx[var_name] = eval(comp_expr, globals(), test_ctx)
                    except Exception as e:
                        if var_name not in [v[0] for v in failed_vars]:
                            failed_vars.append((var_name, str(e)))

            if failed_vars:
                print(f"[WARNING] {len(failed_vars)} variables failed to evaluate and will default to 0.0!")

            failed_rxns = []
            mapped_keys_upper = [k.upper() for k in self.cfg['env_map'].keys()]
            bulk_gases = ['O2', 'N2', 'M', 'RH', 'H2O'] 
            
            for rxn in self.mechanism.reactions:
                expr = getattr(rxn, 'rate_expression', getattr(rxn, 'rate_str', '0.0')).upper()
                expr = expr.replace('@', '**') 
                expr = re.sub(r'(?<=\d)[dD](?=[+-]?\d)', 'E', expr)
                expr = re.sub(r'J\((\d+)\)', r'J\1', expr)
                try: 
                    rxn.safe_expr = compile(expr, '<string>', 'eval')
                    eval(rxn.safe_expr, globals(), test_ctx)
                except Exception as e: 
                    failed_rxns.append((expr, str(e)))
                    rxn.safe_expr = compile('0.0', '<string>', 'eval')

                rxn.reactants_parsed = []
                rxn.reactants_env = []  
                for r in rxn.reactants:
                    coeff, name = 1.0, r.strip().upper()
                    if '*' in name: 
                        coeff, name = float(name.split('*', 1)[0]), name.split('*', 1)[1]
                        
                    if name in self.solver.species_to_idx:
                        rxn.reactants_parsed.append({'idx': self.solver.species_to_idx[name], 'coeff': coeff})
                    elif name in mapped_keys_upper or name in bulk_gases:
                        rxn.reactants_env.append({'name': name, 'coeff': coeff})
                        
                rxn.products_parsed = []
                for p in rxn.products:
                    coeff, name = 1.0, p.strip().upper()
                    if '*' in name: 
                        coeff, name = float(name.split('*', 1)[0]), name.split('*', 1)[1]
                    if name in self.solver.species_to_idx:
                        rxn.products_parsed.append({'idx': self.solver.species_to_idx[name], 'coeff': coeff})

            if failed_rxns:
                print(f"[WARNING] {len(failed_rxns)} reactions failed to evaluate and will be DEAD!")
            else:
                print(">>> [DIAGNOSTICS] All chemical equations evaluated successfully!")

        cont_sources = []
        for s in self.cfg.get('continuous_sources', []):
            if not s.get('target') or not s.get('rate'): continue
            target = s['target'].upper()
            if target in self.solver.species_to_idx:
                cont_sources.append({'idx': self.solver.species_to_idx[target], 'code': compile(s['rate'], '<string>', 'eval')})

        burst_eval = []
        min_dur = 3600.0
        for b_idx, b in enumerate(self.cfg.get('burst_injections', [])):
            if not b.get('target') or not b.get('ppb') or b['ppb'] == "": continue
            target  = b['target'].upper()
            profile = b.get('profile', 'Gaussian')

            if profile == 'TwoReaction':
                # Source term feeds _INJ1 and _INJ2 during the valve-open window.
                # Decay reactions (_INJ1→_INJA→TARGET) are already in mechanism.
                tr = self._tr_specs.get(b_idx)
                if tr is None: continue
                sp1, sp2 = tr['sp1'], tr['sp2']
                if sp1 not in self.solver.species_to_idx: continue
                target_utc  = pd.to_datetime(b['start_utc'], dayfirst=True)
                start_s     = (target_utc - saphir_base_date).total_seconds()
                dur         = float(eval(str(b.get('duration_s', 60))))
                min_dur     = min(min_dur, dur)
                total_conc  = float(b['ppb']) * 1e-9 * self.get_M(start_s)
                fast_frac   = tr['fast_frac']
                rate_per_s  = total_conc / max(dur, 1.0)
                burst_eval.append({
                    'profile': 'tworeaction',
                    'idx1': self.solver.species_to_idx[sp1],
                    'idx2': self.solver.species_to_idx[sp2],
                    'start_s': start_s, 'dur': dur,
                    'rate1': rate_per_s * fast_frac,
                    'rate2': rate_per_s * (1.0 - fast_frac),
                })
            elif target in self.solver.species_to_idx:
                target_utc = pd.to_datetime(b['start_utc'], dayfirst=True)
                start_s = (target_utc - saphir_base_date).total_seconds()
                dur = float(eval(str(b.get('duration_s', 60))))
                min_dur = min(min_dur, dur)
                total_conc = float(b['ppb']) * 1e-9 * self.get_M(start_s)  # molec/cm³
                try:
                    t90_frac = max(0.01, min(0.99, float(b.get('t90_frac', 0.5))))
                except (ValueError, TypeError):
                    t90_frac = 0.5
                # Optional byproduct (e.g. NO2 from NO injection)
                byp_t = b.get('byproduct_target', '').strip().upper()
                byp_f = 0.0
                if byp_t and b.get('byproduct_frac'):
                    try:
                        byp_f = float(eval(str(b['byproduct_frac'])))
                    except Exception:
                        pass
                byp_idx = self.solver.species_to_idx.get(byp_t) if byp_t and byp_f > 0 else None
                if profile == 'Saturating':
                    lam = np.log(10.0) / (t90_frac * dur)
                    norm = (1.0 - np.exp(-lam * dur)) / lam
                    coeff = total_conc / norm
                    burst_eval.append({'profile': 'saturating', 'idx': self.solver.species_to_idx[target],
                                       'start_s': start_s, 'dur': dur, 'lam': lam, 'coeff': coeff,
                                       'byp_idx': byp_idx, 'byp_frac': byp_f})
                else:  # Gaussian (default)
                    # t90_frac > 0.5: use it to set σ (z₉₀ = 1.2816); otherwise fall back to dur/6
                    if t90_frac > 0.5:
                        sigma = dur * (t90_frac - 0.5) / 1.2816
                    else:
                        sigma = dur / 6.0
                    mu = start_s + (dur / 2.0)
                    coeff = total_conc / (sigma * np.sqrt(2 * np.pi))
                    burst_eval.append({'profile': 'gaussian', 'idx': self.solver.species_to_idx[target],
                                       'mu': mu, 'sigma': sigma, 'coeff': coeff,
                                       'byp_idx': byp_idx, 'byp_frac': byp_f})

        # --- H2O injection impurities (CO, NO2, …) tied to H2O injection window ---
        h2o_imps = []
        h2o_imp_defs = self.cfg.get('h2o_injection_impurities', [])
        if h2o_imp_defs:
            h2o_inj = next(
                (b for b in self.cfg.get('burst_injections', [])
                 if b.get('target', '').upper() == 'H2O'
                 or b.get('name', '').upper() in ('H2O', 'WATER')),
                None
            )
            if h2o_inj and h2o_inj.get('start_utc'):
                h2o_start = (pd.to_datetime(h2o_inj['start_utc'], dayfirst=True)
                             - saphir_base_date).total_seconds()
                h2o_dur   = float(eval(str(h2o_inj.get('duration_s', 600))))
                M_h2o     = self.get_M(h2o_start)
                for imp in h2o_imp_defs:
                    tgt = imp.get('target', '').upper()
                    if not tgt or tgt not in self.solver.species_to_idx:
                        continue
                    try:
                        rate_ppb_per_min = float(eval(str(imp.get('rate_ppb_per_min', '0'))))
                    except Exception:
                        continue
                    if rate_ppb_per_min <= 0:
                        continue
                    h2o_imps.append({
                        'idx':     self.solver.species_to_idx[tgt],
                        'start_s': h2o_start,
                        'dur':     h2o_dur,
                        'rate':    rate_ppb_per_min * 1e-9 * M_h2o / 60.0,
                    })
                    print(f"  [H2O Impurity] {tgt}: {rate_ppb_per_min} ppb/min during H2O injection")

        def system_ode(t, y):
            dydt = np.zeros(self.solver.num_species)
            M_now = self.get_M(t)
            T_now = self.ctx_data['TEMP'](t)

            ctx = SafeDict()

            # 1. Baseline Environment Data (Loads everything from .enz)
            active_constraints = set()
            for k, v_func in self.ctx_data.items():
                try: 
                    ctx[k.upper()] = float(v_func(t))
                    if k.upper() in self.cfg['env_map']:
                        active_constraints.add(k.upper())
                except Exception:
                    pass

            # --- THE WATER AUTO-BALANCER ---
            if 'RH' in ctx and 'H2O' not in ctx:
                p_sat = np.exp(21.36469 - 5339.66 / T_now)
                press = float(self.ctx_data['PRESS'](t)) if 'PRESS' in self.ctx_data else 1013.25
                ctx['H2O'] = (ctx['RH'] / 100.0) * (p_sat / press) 
                
            if 'H2O' in ctx:
                # ctx['H2O'] is already a corrected VMR fraction (scaled + noise-floored at setup).
                h2o_vmr = ctx['H2O']
                ctx['H2O'] = h2o_vmr * M_now
                # Back-compute RH (%) from H2O VMR so artifact expressions work correctly
                if 'RH' not in ctx:
                    if h2o_vmr > 0:
                        p_sat = np.exp(21.36469 - 5339.66 / T_now)
                        press = float(self.ctx_data['PRESS'](t)) if 'PRESS' in self.ctx_data else 1013.25
                        ctx['RH'] = h2o_vmr * press / p_sat * 100.0
                    else:
                        ctx['RH'] = 0.0
            else:
                # No H2O channel at all — use 1% fallback so kinetics don't crash
                ctx['H2O'] = 0.01 * M_now

            ctx.update({
                "T": T_now, "TEMP": T_now, "M": M_now,
                "O2": 0.2095 * M_now,
                "N2": 0.7809 * M_now,
                "np": np, "exp": np.exp, "EXP": np.exp,
                "LOG10": np.log10, "log10": np.log10,
                "LOG": np.log, "log": np.log,
                "COS": np.cos, "SIN": np.sin, "TAN": np.tan, "SQRT": np.sqrt
            })

            # Load mcm_kinetics: COMPRC gets loaded as a direct callable (no KMT wrapping).
            # KMT lambdas are loaded but will be overwritten by COMPRC-based var_assignments
            # which produce correct scalar values at the current T and M.
            for func_name in dir(mcm_kinetics):
                if not func_name.startswith('_'):
                    func = getattr(mcm_kinetics, func_name)
                    if func_name.startswith('KMT'):
                        ctx[func_name.upper()] = lambda *args, f=func: f(*args, T_now, M_now)
                    else:
                        ctx[func_name.upper()] = func  # COMPRC loaded here as a true callable

            for i, name in enumerate(self.solver.species):
                ctx[name] = float(y[i])

            # Re-apply constrained env_map species AFTER the species loop so that
            # measured values (e.g. H2O, NO, NO2) always override the ODE state in
            # rate expressions.  The ODE may track these as species internally, but
            # for chemistry evaluation the measurement must win.
            for k in active_constraints:
                try:
                    ctx[k] = float(self.ctx_data[k.lower() if k.lower() in self.ctx_data else k](t))
                except Exception:
                    pass
            # Re-run the water block so H2O molec/cm³ and RH are consistent with
            # the restored measured H2O VMR.
            if 'H2O' in active_constraints:
                # ctx['H2O'] re-read from the pre-corrected interpolation — already a VMR fraction.
                h2o_vmr = ctx['H2O']
                ctx['H2O'] = h2o_vmr * M_now
                if h2o_vmr > 0:
                    p_sat = np.exp(21.36469 - 5339.66 / T_now)
                    press = float(self.ctx_data['PRESS'](t)) if 'PRESS' in self.ctx_data else 1013.25
                    ctx['RH'] = h2o_vmr * press / p_sat * 100.0
                else:
                    ctx['RH'] = 0.0

            # --- PASS 1: Evaluate constants and COMPRC-based KMT rates ---
            # This computes COSX, LNO2, MNO2, NNO2, JSCALE, KMT01-KMT14, KFPAN, KBPAN, etc.
            # PHTFKT-based J-values will fail here (PHTFKT not yet callable) — that is expected.
            for var_name, comp_expr in self.var_assignments:
                if var_name in active_constraints:
                    continue
                try:
                    ctx[var_name] = eval(comp_expr, globals(), ctx)
                except:
                    ctx[var_name] = 0.0

            # --- BUILD PHTFKT PYTHON CALLABLE using COSX and measured JNO2 from pass 1 ---
            _cosx = max(float(ctx.get('COSX', np.cos(float(ctx.get('SZ', 45.0)) * 0.0174533))), 0.01)
            _jno2 = float(ctx.get('JNO2', 0.0))
            _lno2 = float(ctx.get('LNO2', 1.165e-2))
            _mno2 = float(ctx.get('MNO2', 0.244))
            _nnno2 = float(ctx.get('NNO2', -0.267))
            _hv   = float(ctx.get('HV', 0.0))

            def _phtfkt(a, b, c, _c=_cosx, _j=_jno2, _l=_lno2, _m=_mno2, _n=_nnno2, _h=_hv):
                if _h == 0 or _c <= 0:
                    return 0.0
                cs_j = a * _c**b * np.exp(c / _c)
                if _j > 0:
                    jno2_cs = _l * _c**_m * np.exp(_n / _c)
                    scale = min(_j / jno2_cs, 1.0) if jno2_cs > 0 else 1.0
                else:
                    scale = 1.0
                return cs_j * scale

            ctx['PHTFKT']    = _phtfkt
            ctx['PHTFKT_CS'] = lambda a, b, c, _c=_cosx: a * _c**b * np.exp(c / _c)
            ctx['AMAX']      = np.maximum
            ctx['AMIN']      = np.minimum

            # --- PASS 2: Evaluate J-values using the real PHTFKT callable ---
            for var_name, comp_expr in self.var_assignments:
                if var_name in active_constraints:
                    continue
                try:
                    ctx[var_name] = eval(comp_expr, globals(), ctx)
                except:
                    ctx[var_name] = 0.0

            # --- PHOTOLYSIS AUTO-LINKER (runs AFTER var_assignments as fallback) ---
            # Measured j-values always beat JSCALE/PHTFKT approximations.
            # J7 and J3 are defined in anorg.mcm331.phot.rate as JSCALE*l*HV (non-zero),
            # so the ==0 guard would prevent the override — they must be unconditional.
            # MCM J-numbering: J1=JO1D (O3->O1D), J4=JNO2, J7=JHONO, J3=JH2O2
            if ctx.get('JO1D',   0.0) > 0: ctx['J1']  = ctx['JO1D']
            if ctx.get('JNO2',   0.0) > 0: ctx['J4']  = ctx['JNO2']
            if ctx.get('JHONO',  0.0) > 0: ctx['J7']  = ctx['JHONO']   # unconditional: beats JSCALE approx
            if ctx.get('JH2O2',  0.0) > 0: ctx['J3']  = ctx['JH2O2']   # unconditional: beats JSCALE approx
            if ctx.get('J11', 0.0) == 0.0 and ctx.get('JHCHO_R',0.0) > 0: ctx['J11'] = ctx['JHCHO_R']
            if ctx.get('J12', 0.0) == 0.0 and ctx.get('JHCHO_M',0.0) > 0: ctx['J12'] = ctx['JHCHO_M']

            for rxn in self.mechanism.reactions:
                try: 
                    k = eval(rxn.safe_expr, globals(), ctx)
                except:
                    k = 0.0

                if callable(k):
                    k = 0.0  # KMT lambda was not overwritten by var_assignments — config error
                if k != 0.0:
                    v = k
                    for r in rxn.reactants_parsed: v *= (y[r['idx']] ** r['coeff'])
                    for r in rxn.reactants_env: v *= (ctx[r['name']] ** r['coeff']) 
                    
                    for r in rxn.reactants_parsed: dydt[r['idx']] -= (v * r['coeff'])
                    for p in rxn.products_parsed: dydt[p['idx']] += (v * p['coeff'])

            k_dil = (self.ctx_data['FLOW'](t) + 0.6) / (self.cfg['chamber_vol'] * 3600)
            dydt -= y * k_dil

            for s in cont_sources:
                dydt[s['idx']] += eval(s['code'], globals(), ctx)

            for b in burst_eval:
                if b['profile'] == 'tworeaction':
                    if b['start_s'] <= t <= b['start_s'] + b['dur']:
                        dydt[b['idx1']] += b['rate1']
                        dydt[b['idx2']] += b['rate2']
                elif b['profile'] == 'saturating':
                    dt = t - b['start_s']
                    if 0.0 <= dt <= b['dur']:
                        rate = b['coeff'] * np.exp(-b['lam'] * dt)
                        dydt[b['idx']] += rate
                        if b['byp_idx'] is not None:
                            dydt[b['byp_idx']] += rate * b['byp_frac']
                else:  # gaussian
                    rate = b['coeff'] * np.exp(-0.5 * ((t - b['mu']) / b['sigma'])**2)
                    dydt[b['idx']] += rate
                    if b['byp_idx'] is not None:
                        dydt[b['byp_idx']] += rate * b['byp_frac']

            for imp in h2o_imps:
                if imp['start_s'] <= t <= imp['start_s'] + imp['dur']:
                    dydt[imp['idx']] += imp['rate']
            
            self.eval_count += 1
            if self.eval_count % 500 == 0 or t > self.last_progress_print_time + 300:
                progress = max(0.0, min((t - self.t_start) / (self.t_end - self.t_start), 1.0))
                
                eta_str = ""
                if progress > 0:
                    elapsed = time.time() - self.real_start_time
                    total_est = elapsed / progress
                    rem_secs = int(total_est - elapsed)
                    eta_str = f" | ETA: {rem_secs // 60}m {rem_secs % 60}s"
                
                print(f"  ... ODE Eval Step: {self.eval_count} | Sim Time: {t-self.t_start:.1f}s / {self.t_end-self.t_start:.1f}s | Fraction Done: {progress:.3f} ({progress*100:.1f}%){eta_str}")
                if t > self.last_progress_print_time + 300:
                    self.last_progress_print_time = t

            return dydt

        solver_method = 'Radau'
        step_limit = min(min_dur / 2.0, 60.0)

        # --- BUILD JACOBIAN SPARSITY PATTERN ---
        # Uses compiled code object .co_names to detect which species each
        # var_assignment expression depends on (e.g. RO2 = sum of 6 peroxy radicals).
        # For each reaction the non-zero Jacobian columns are:
        #   direct reactant species  +  any species-dependent vars in the rate expression
        # Rows are all affected species (reactants + products).
        # Diagonal added for the dilution term (dydt_i -= y_i * k_dil).
        _sp_idx = self.solver.species_to_idx
        N = self.solver.num_species
        sparsity = lil_matrix((N, N), dtype=np.int8)

        # diagonal: dilution
        for i in range(N):
            sparsity[i, i] = 1

        # var_assignment expressions that depend on species concentrations (e.g. RO2 sum)
        var_species_deps = {}
        for var_name, comp_expr in self.var_assignments:
            dep = frozenset(_sp_idx[n] for n in comp_expr.co_names if n in _sp_idx)
            if dep:
                var_species_deps[var_name] = dep

        for rxn in self.mechanism.reactions:
            reactant_idxs = frozenset(r['idx'] for r in rxn.reactants_parsed)
            affected_idxs = reactant_idxs | frozenset(p['idx'] for p in rxn.products_parsed)
            if not affected_idxs:
                continue  # zeroth-order source with no ODE species on either side

            # Extra species the rate constant depends on via var_assignments (e.g. RO2)
            rate_extra = frozenset(
                idx
                for vname in rxn.safe_expr.co_names
                if vname in var_species_deps
                for idx in var_species_deps[vname]
            )
            all_rate_deps = reactant_idxs | rate_extra
            if not all_rate_deps:
                continue  # pure parameter-driven rate (T/M/J only), no species dependence

            for ai in affected_idxs:
                for rj in all_rate_deps:
                    sparsity[ai, rj] = 1

        sparsity_csr = sparsity.tocsr()
        nnz = sparsity_csr.nnz
        density = nnz / (N * N) * 100
        print(f">>> [SPARSE] Jacobian: {N}×{N}, {nnz} non-zeros ({density:.1f}% dense) vs {N*N} dense")

        self.t_start = self.t_span[0]
        self.t_end = self.t_span[-1]
        self.last_progress_print_time = self.t_start
        self.eval_count = 0
        self.real_start_time = time.time()

        print(f">>> [DEBUG] Solving for {len(self.mechanism.species)} species using {solver_method} (Max Step: {step_limit}s)...")
        self.sol = solve_ivp(
            system_ode,
            [self.t_span[0], self.t_span[-1]],
            y0,
            method=solver_method,
            t_eval=self.t_span,
            max_step=step_limit,
            rtol=1e-3,
            atol=1e4,
            jac_sparsity=sparsity_csr
        )
        if not self.sol.success: print(f"\n[WARNING] ODE Solver stopped early!\nReason: {self.sol.message}")

    def _export_results(self):
        print("\n>>> [DATA] Exporting results to CSV...")
        
        # Transpose sol.y so rows are time steps and columns are species.
        # Clip negatives here (solver internal arithmetic can produce small negatives near 0).
        df = pd.DataFrame(np.maximum(self.sol.y, 0.0).T, columns=self.solver.species, index=self.sol.t)
        df.index.name = "Time_s"
        
        # Calculate Air Density (M) for every time step to allow for ppb conversions later
        df['M'] = [self.get_M(t) for t in self.sol.t]
        
        # Merge SAPHIR Interpolated Data into the Results CSV for direct comparison
        for k, v_func in self.ctx_data.items():
            if callable(v_func):
                try:
                    saphir_vals = v_func(self.sol.t)
                    if isinstance(saphir_vals, (float, int)):
                        saphir_vals = np.full_like(self.sol.t, saphir_vals)
                    df[f"SAPHIR_{k}"] = saphir_vals
                except Exception:
                    pass

        # Evaluate diagnostic var_assignments (e.g. RO2 = sum of RO2 radicals) that are
        # not ODE species but are defined as expressions over species that ARE in df.
        _derived_ctx = {col: df[col].values for col in df.columns}
        _derived_ctx['np'] = np
        for var_name, compiled_expr in self.var_assignments:
            if var_name in df.columns:
                continue  # already a species column
            try:
                result = eval(compiled_expr, globals(), _derived_ctx)
                if isinstance(result, np.ndarray) and result.shape == (len(df),):
                    df[var_name] = np.maximum(result, 0.0)
            except Exception:
                pass

        # --- O1D pseudo-steady-state correction ------------------------------------
        # O1D has a lifetime of ~1 ns (quenched mainly by N2 and O2).  The ODE
        # solver (atol=1e4) cannot resolve it accurately — atol is ~10⁶× the
        # actual PSS value (~0.01 molec/cm³), so the integrated column carries
        # numerical noise / an initial-condition spike.  Replace with the
        # analytical PSS:  [O1D]_pss = J1·[O3] / (k_N2·[N2] + k_O2·[O2] + k_H2O·[H2O])
        if 'O1D' in df.columns:
            _T   = (df['SAPHIR_TEMP'].values if 'SAPHIR_TEMP' in df.columns else
                    df['SAPHIR_T'].values     if 'SAPHIR_T'    in df.columns else
                    np.full(len(df), 298.15))
            _M   = df['M'].values
            _HV  = df['SAPHIR_HV'].values  if 'SAPHIR_HV'  in df.columns else np.ones(len(df))
            _o3  = df['O3'].values          if 'O3'         in df.columns else np.zeros(len(df))
            # J1: prefer direct SAPHIR JO1D measurement, fall back to J1 diagnostic column
            _j1  = (df['SAPHIR_JO1D'].values if 'SAPHIR_JO1D' in df.columns else
                    df['J1'].values            if 'J1'          in df.columns else
                    np.zeros(len(df)))
            _j1  = np.where(_HV > 0.5, _j1, 0.0)   # zero when lights are off
            # H2O: SAPHIR_H2O is stored as VMR fraction → convert to molec/cm³
            _h2o = (df['SAPHIR_H2O'].values * _M) if 'SAPHIR_H2O' in df.columns else np.zeros(len(df))
            _kq  = (2.15e-11 * np.exp(110.0 / _T) * 0.7809 * _M +
                    3.30e-11 * np.exp( 55.0 / _T) * 0.2095 * _M +
                    2.20e-10 * _h2o)
            df['O1D'] = np.where(_kq > 0, _j1 * _o3 / np.maximum(_kq, 1e-30), 0.0)
            print("  [O1D] Replaced ODE output with pseudo-steady-state (J1·O3/k_quench)")

        # Y background: diagnostic column showing background VOC concentration (molec/cm³)
        if self.cfg.get('artifacts', {}).get('y_active', False):
            try:
                y_eqn = str(self.cfg['artifacts']['y_eqn'])
                df['Y'] = eval(y_eqn, {"__builtins__": None}, {"M": df['M'].values}) * 1e-9 * df['M'].values
            except Exception:
                pass

        out_path = os.path.join(self.root, "simulation_results.csv")
        df.to_csv(out_path)
        
        print(f"  ... Saved to {out_path}")
        return out_path

    def _show_results(self, out_path):
        from results_viewer import show_results
        show_results(out_path)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python scenario_runner.py <path_to_scenario_config.json>")
        sys.exit(1)
        
    config_path = sys.argv[1]
    runner = ScenarioRunner(config_path)
    runner.run()