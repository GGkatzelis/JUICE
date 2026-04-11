import os
import json
import numpy as np
import pandas as pd 
from scipy.integrate import solve_ivp
from scipy.sparse import coo_matrix
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
        # Look for ENZ file: prefer 'kammer*.enz', fall back to any '*.enz'
        enz_candidates = [f for f in os.listdir(self.root) if f.lower().endswith('.enz')]
        enz_file = next((f for f in enz_candidates if 'kammer' in f.lower()), None)
        if enz_file is None:
            enz_file = enz_candidates[0] if enz_candidates else None
        if enz_file is None:
            raise FileNotFoundError(f"No .enz file found in {self.root}")
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
        temp_fn = self.ctx_data.get('TEMP') or self.ctx_data.get('T')
        press_fn = self.ctx_data.get('PRESS') or self.ctx_data.get('P')
        if temp_fn is None:
            raise KeyError("No temperature source found in ctx_data (expected 'TEMP' or 'T')")
        if press_fn is None:
            raise KeyError("No pressure source found in ctx_data (expected 'PRESS' or 'P')")
        temp = temp_fn(t)
        press = press_fn(t)
        # Use the EASY/FACSIMILE DENSITY convention: M = P / (T * 1.379e-19)
        # where P is in hPa and the constant 1.379e-19 = kB in hPa·cm³/K.
        # This matches the EASY DENSITY=FUNCTION[(%1)/((%2)*1.379E-19)] macro
        # used in every MCM .incl file. Using the exact kB (1.380649e-23 J/K)
        # instead would introduce a 0.116% offset in M that propagates into
        # every M-dependent rate expression.
        return press / (temp * 1.379e-19)

    def _build_mechanism(self):
        chem_path = os.path.join(self.root, "chemistry")
        parser = EasyParser(base_path=chem_path)
        self.mechanism = Mechanism()
        
        # Sort files to ensure .rate definitions are parsed before .mech files
        mech_files = sorted(self.cfg['mechanism_files'], key=lambda x: ".rate" not in x)
        
        for f_name in mech_files:
            f_path = os.path.join(chem_path, f_name)
            if not os.path.exists(f_path):
                print(f"[WARNING] Mechanism file not found (skipped): {f_path}")
                continue
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
                    
        # Post-process all accumulated reactions using the fully-populated rates_map.
        # When files are parsed in sequence, a later SAPHIR override file may update
        # rates_map entries for reactions already materialised in an earlier .mech file.
        # Re-applying the final rates_map here ensures those overrides take effect.
        def _rxn_canon(r):
            reac_str = '+'.join(str(x) for x in r.reactants)
            prod_str  = '+'.join(str(x) for x in r.products)
            return EasyParser._canon_rxn_key(reac_str, prod_str)

        for rxn in self.mechanism.reactions:
            reac_str = '+'.join(str(x) for x in rxn.reactants)
            prod_str  = '+'.join(str(x) for x in rxn.products)
            raw_key   = reac_str.replace(' ', '').upper() + '-->' + prod_str.replace(' ', '').upper()
            canon_key = _rxn_canon(rxn)
            final_rate = parser.rates_map.get(raw_key) or parser.rates_map.get(canon_key)
            if final_rate is None:
                # FACSIMILE rate-inheritance: bare reaction with no explicit rate
                # inherits the -/->-disabled rate for the same reactant pair.
                reac_key = frozenset(r.upper() for r in rxn.reactants)
                final_rate = parser._reactant_fallback_rates.get(reac_key)
            if final_rate is not None:
                rxn.rate_expression = final_rate

        # Remove reactions disabled by -/-> across all parsed files.
        # These are typically FACSIMILE override lines that replace an MCM reaction
        # with a SAPHIR-specific version (e.g. NO3+hv→NO2+O disabled, replaced by
        # NO3+hv→NO2+O3).  The parser may have materialised the original reaction in
        # an earlier file before encountering the -/-> in a later include file.
        if parser.disabled_rxn_keys:
            before_d = len(self.mechanism.reactions)
            self.mechanism.reactions = [
                r for r in self.mechanism.reactions
                if _rxn_canon(r) not in parser.disabled_rxn_keys
            ]
            removed_d = before_d - len(self.mechanism.reactions)
            if removed_d:
                print(f"  ... Removed {removed_d} reaction(s) disabled by -/-> notation")

        # Deduplicate reactions by canon key: when the same reaction (same reactants+
        # products) appears in multiple files (e.g. MCM phot.mech AND saphir override),
        # keep only the LAST occurrence so the later SAPHIR file's version wins.
        # This mirrors the single-file _rxn_map deduplication in EasyParser.parse().
        seen_canon = {}
        for i, rxn in enumerate(self.mechanism.reactions):
            seen_canon[_rxn_canon(rxn)] = i          # last index wins
        dedup_indices = set(seen_canon.values())
        before_dd = len(self.mechanism.reactions)
        self.mechanism.reactions = [
            r for i, r in enumerate(self.mechanism.reactions)
            if i in dedup_indices
        ]
        removed_dd = before_dd - len(self.mechanism.reactions)
        if removed_dd:
            print(f"  ... Deduplicated {removed_dd} reaction(s) (multi-file overrides)")

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

        # Detect whether the imported mechanism contains explicit per-species
        # dilution reactions ([SPECIES + DIL → ]). If so, the analytical
        # dilution path must be skipped to avoid double-counting (the explicit
        # reactions and the analytical decay use the same DILUTE rate).
        self._has_explicit_dil = any(
            any(r.strip().upper().split('*')[-1].strip() == 'DIL' for r in rxn.reactants)
            for rxn in self.mechanism.reactions
        )
        if self._has_explicit_dil:
            print(f"  [Dilution] Detected explicit +DIL reactions in mechanism — "
                  f"analytical dilution disabled to avoid double-counting")

        # IDL-mode computed parameters (DILUTE, DIL, …) — evaluated each ODE step,
        # overriding the ctx_data DILUTE=0 stub so explicit +DIL reactions use the
        # correct rate constant.
        for p_name, p_expr in self.cfg.get('idl_parameters', {}).items():
            p_clean = p_expr.upper().replace('@', '**')
            p_clean = re.sub(r'(?<=\d)[dD](?=[+-]?\d)', 'E', p_clean)
            try:
                self.var_assignments.append((p_name.upper(), compile(p_clean, '<string>', 'eval')))
                print(f"  [IDL param] {p_name} = {p_expr}")
            except Exception as e:
                print(f"[WARNING] Could not compile IDL parameter {p_name}: {e}")
            
        all_species = set(self.mechanism.species)
        for spec in self.cfg.get('initials', {}).keys():
            all_species.add(spec.upper())
        for b in self.cfg.get('burst_injections', []):
            if b.get('target'): all_species.add(b['target'].upper())
        for s in self.cfg.get('continuous_sources', []):
            if s.get('target'): all_species.add(s['target'].upper())
            
        for r in self.cfg.get('custom_reactions', []):
            # Allow empty products (pure decay) or empty reactants (zeroth-order source); only rate is required.
            if not r.get('rate'): continue
            reacts_str = r.get('reactants', '') or ''
            prods_str  = r.get('products',  '') or ''
            parts = re.split(r'\s*\+\s*|\s*->\s*|\s*=\s*', reacts_str + " + " + prods_str)
            for p in parts:
                clean_p = p.strip().upper()
                clean_p = re.sub(r'^[\d\.]+\s*\*\s*', '', clean_p)
                if clean_p and not clean_p[0].isdigit(): all_species.add(clean_p)

            reactants = [x.strip().upper() for x in reacts_str.split('+') if x.strip()]
            products  = [x.strip().upper() for x in prods_str.split('+')  if x.strip()]
            self.mechanism.reactions.append(Reaction(reactants, products, r['rate']))

        # --- IDL chamber_source_reactions: SHONO/SHCHO (zeroth-order) & H2O dark source ---
        # SHONO/SHCHO are FACSIMILE CONST pseudo-species: their values are already folded
        # into the rate expression, so reactants are empty (zeroth-order production).
        # H2O --> HONO has H2O as a real env reactant (first-order, from env_map).
        mapped_up_set = {k.upper() for k in self.cfg.get('env_map', {})}
        for r in self.cfg.get('chamber_source_reactions', []):
            if not r.get('products') or not r.get('rate'):
                continue
            reacts_str = r.get('reactants', '')
            reactants = [x.strip().upper() for x in reacts_str.split('+') if x.strip()] if reacts_str.strip() else []
            products = [x.strip().upper() for x in r['products'].split('+') if x.strip()]
            # Only add to all_species if NOT an env-map entry (H2O etc. are env, not ODE)
            for s in reactants:
                clean_s = re.sub(r'^[\d\.]+\s*\*\s*', '', s.strip())
                if clean_s and not clean_s[0].isdigit() and clean_s not in mapped_up_set and clean_s not in self.ctx_data:
                    all_species.add(clean_s)
            for s in products:
                clean_s = re.sub(r'^[\d\.]+\s*\*\s*', '', s.strip())
                if clean_s and not clean_s[0].isdigit():
                    all_species.add(clean_s)
            self.mechanism.reactions.append(Reaction(reactants, products, r['rate']))

        # --- IDL injection_reactions: switch-gated zeroth-order sources + chain reactions ---
        for r in self.cfg.get('injection_reactions', []):
            if not r.get('products') or not r.get('rate'):
                continue
            reacts_str = r.get('reactants', '')
            reactants = [x.strip().upper() for x in reacts_str.split('+') if x.strip()] if reacts_str.strip() else []
            products_str = r['products']
            products = [x.strip().upper() for x in products_str.split('+') if x.strip()]
            for s in reactants:
                clean_s = re.sub(r'^[\d\.]+\s*\*\s*', '', s.strip())
                # Skip ENZ switch columns (they are constrained inputs in ctx_data,
                # not ODE species; the gating is already in the rate expression).
                if clean_s and not clean_s[0].isdigit() and clean_s not in self.ctx_data:
                    all_species.add(clean_s)
            for s in products:
                clean_s = re.sub(r'^[\d\.]+\s*\*\s*', '', s.strip())
                if clean_s and not clean_s[0].isdigit():
                    all_species.add(clean_s)
            self.mechanism.reactions.append(Reaction(reactants, products, r['rate']))

        # --- IDL dilution reactions: explicit per-species decay at rate DILUTE --------------------
        # These come from k[SPECIES+DIL-->]=CONST(DILUTE) lines in the EASY master file.
        # Used only when analytical_dilution=False (general .dil.mech/.dil.rate commented out).
        for r in self.cfg.get('idl_dilution_reactions', []):
            if not r.get('reactants') or not r.get('rate'):
                continue
            reacts_str = r['reactants']
            reactants  = [x.strip().upper() for x in reacts_str.split('+') if x.strip()]
            prods_str  = r.get('products', '')
            products   = [x.strip().upper() for x in prods_str.split('+') if x.strip()] if prods_str.strip() else []
            for s in reactants:
                clean_s = re.sub(r'^[\d\.]+\s*\*\s*', '', s)
                if clean_s and not clean_s[0].isdigit() and clean_s not in self.ctx_data:
                    all_species.add(clean_s)
            for s in products:
                clean_s = re.sub(r'^[\d\.]+\s*\*\s*', '', s)
                if clean_s and not clean_s[0].isdigit():
                    all_species.add(clean_s)
            self.mechanism.reactions.append(Reaction(reactants, products, r['rate']))
        n_dil = len(self.cfg.get('idl_dilution_reactions', []))
        if n_dil:
            print(f"  [IDL Dilution] Added {n_dil} per-species dilution reaction(s) from EASY file")

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
            # IDL bursts are informational summaries — the actual delivery is
            # handled by the [QSIG]→TARGET injection_reactions chain. Skipping
            # avoids double-injection.
            if b.get('source') == 'IDL':
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

        # --- PSS species: O1D lifetime ~1 ns, far below ODE atol=1e4.
        # Remove from ODE state; compute analytically in ctx so reactions using
        # O1D as a reactant (e.g. O1D+H2O→2OH) see the correct PSS concentration.
        self._pss_species = set()
        if 'O1D' in all_species:
            all_species.discard('O1D')
            self._pss_species.add('O1D')
            _o1d_pss = (
                "np.maximum(JO1D, 0.0) * np.maximum(O3, 0.0) / "
                "np.maximum("
                "2.0e-11 * N2 * np.exp(130.0 / T) + "
                "3.20e-11 * O2 * np.exp(67.0 / T) + "
                "2.14e-10 * H2O, 1e-30)"
            )
            try:
                self.var_assignments.append(('O1D', compile(_o1d_pss, '<O1D_pss>', 'eval')))
                print("  [PSS] O1D removed from ODE state → J1·O3/k_quench analytical expression added")
            except Exception as e:
                print(f"  [WARNING] Could not compile O1D PSS expression: {e}")

        # --- O atom: kept in ODE (NOT PSS) ---
        # FACSIMILE also integrates O in the ODE, so both models share the same numerical
        # treatment. Applying PSS in JUICE but not FACSIMILE makes the O→O3 production
        # rate diverge and breaks O3 agreement (tested: PSS caused O3 slope ±50% swings).
        # The ODE O value may be stiff (~76× PSS analytically) but agrees with FACSIMILE
        # exactly because both run the same reactions at the same numerical rates.
        # Post-processing (_export_results) replaces O in the output CSV with the PSS
        # formula for a physically meaningful reported value without affecting chemistry.

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
            bulk_gases = ['O2', 'N2', 'M', 'RH', 'H2O'] + list(getattr(self, '_pss_species', set()))
            
            for rxn in self.mechanism.reactions:
                expr = getattr(rxn, 'rate_expression', getattr(rxn, 'rate_str', '0.0')).upper()
                expr = expr.replace('@', '**') 
                expr = re.sub(r'(?<=\d)[dD](?=[+-]?\d)', 'E', expr)
                expr = re.sub(r'J\((\d+)\)', r'J\1', expr)
                try:
                    rxn.safe_expr = compile(expr, '<string>', 'eval')
                    eval(rxn.safe_expr, globals(), test_ctx)
                    rxn.rate_str_clean = expr
                except Exception as e:
                    failed_rxns.append((expr, str(e)))
                    rxn.safe_expr = compile('0.0', '<string>', 'eval')
                    rxn.rate_str_clean = '0.0'

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

        # ── ODE HOT-PATH PRECOMPILATION ──────────────────────────────────────────────
        # Pre-collect mcm_kinetics functions once — avoids dir()+getattr() every ODE call
        _kmt_funcs   = [(n, getattr(mcm_kinetics, n)) for n in dir(mcm_kinetics)
                        if not n.startswith('_') and n.upper().startswith('KMT')]
        _other_funcs = [(n, getattr(mcm_kinetics, n)) for n in dir(mcm_kinetics)
                        if not n.startswith('_') and not n.upper().startswith('KMT')]

        # Cache species name list for fast ctx update (avoids repeated attribute lookup)
        _species_names = list(self.solver.species)

        # Generate a single compiled function that evaluates all reaction rate constants
        # in one call.  Unique names are unpacked from ctx into locals once (one dict
        # lookup per name, not per reaction), then all expressions run as plain Python
        # bytecode — eliminates N_reactions × eval() dispatch overhead per ODE step.
        _rate_co_names = set()
        for rxn in self.mechanism.reactions:
            _rate_co_names.update(rxn.safe_expr.co_names)

        _fn_lines = ["def _compute_rates(_ctx, _k):"]
        for name in sorted(_rate_co_names):
            _fn_lines.append(f"    {name} = _ctx['{name}']")
        for i, rxn in enumerate(self.mechanism.reactions):
            _fn_lines.append(f"    _k[{i}] = {rxn.rate_str_clean}")
        _fn_src = "\n".join(_fn_lines)

        _batch_ns = {
            '__builtins__': {'abs': abs, 'max': max, 'min': min, 'float': float},
            'exp': np.exp, 'EXP': np.exp,
            'log': np.log,  'LOG': np.log,
            'log10': np.log10, 'LOG10': np.log10,
            'cos': np.cos,  'sin': np.sin,
            'sqrt': np.sqrt, 'SQRT': np.sqrt,
            'np': np,
        }
        exec(_fn_src, _batch_ns)
        _compute_rates = _batch_ns['_compute_rates']
        _k_buf = np.zeros(len(self.mechanism.reactions))
        _v_buf = np.zeros(len(self.mechanism.reactions))  # reaction velocity buffer for Jacobian
        print(f">>> [ODE] Precompiled batch rate function: {len(self.mechanism.reactions)} reactions, "
              f"{len(_rate_co_names)} unique rate vars")

        # Store for use in _export_kinetics (called after _solve completes)
        self._compute_rates_fn = _compute_rates
        self._kmt_funcs_cached = _kmt_funcs
        self._other_funcs_cached = _other_funcs

        # PSS species set — captured in system_ode closure for O1D (and future additions)
        _pss_species_set = getattr(self, '_pss_species', set())

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
            # IDL bursts are informational summaries — the [QSIG]→TARGET chain
            # in injection_reactions is the actual delivery. Skip to avoid 2x.
            if b.get('source') == 'IDL':
                continue
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
            y = np.maximum(y, 0.0)  # clip: BDF can produce slightly negative y; prevents NaN in k*y^n
            t_abs = t + _t_offset   # convert relative solver time → absolute interpolation time
            dydt = np.zeros(self.solver.num_species)
            M_now = self.get_M(t_abs)
            _temp_fn = self.ctx_data.get('TEMP') or self.ctx_data.get('T')
            T_now = _temp_fn(t_abs)

            ctx = SafeDict()

            # 1. Baseline Environment Data (Loads everything from .enz)
            active_constraints = set()
            for k, v_func in self.ctx_data.items():
                try:
                    ctx[k.upper()] = float(v_func(t_abs))
                    if k.upper() in self.cfg['env_map']:
                        active_constraints.add(k.upper())
                except Exception:
                    pass

            # --- THE WATER AUTO-BALANCER ---
            if 'RH' in ctx and 'H2O' not in ctx:
                p_sat = np.exp(21.36469 - 5339.66 / T_now)
                press = float(self.ctx_data['PRESS'](t_abs)) if 'PRESS' in self.ctx_data else 1013.25
                ctx['H2O'] = (ctx['RH'] / 100.0) * (p_sat / press)

            if 'H2O' in ctx:
                # ctx['H2O'] is already a corrected VMR fraction (scaled + noise-floored at setup).
                h2o_vmr = ctx['H2O']
                ctx['H2O'] = h2o_vmr * M_now
                # Back-compute RH (%) from H2O VMR so artifact expressions work correctly
                if 'RH' not in ctx:
                    if h2o_vmr > 0:
                        p_sat = np.exp(21.36469 - 5339.66 / T_now)
                        press = float(self.ctx_data['PRESS'](t_abs)) if 'PRESS' in self.ctx_data else 1013.25
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
            for func_name, func in _kmt_funcs:
                ctx[func_name.upper()] = lambda *args, f=func, _T=T_now, _M=M_now: f(*args, _T, _M)
            for func_name, func in _other_funcs:
                ctx[func_name.upper()] = func

            ctx.update(zip(_species_names, y.tolist()))

            # Re-apply constrained env_map species AFTER the species loop so that
            # measured values (e.g. H2O, NO, NO2) always override the ODE state in
            # rate expressions.  The ODE may track these as species internally, but
            # for chemistry evaluation the measurement must win.
            for k in active_constraints:
                try:
                    ctx[k] = float(self.ctx_data[k.lower() if k.lower() in self.ctx_data else k](t_abs))
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
                    press = float(self.ctx_data['PRESS'](t_abs)) if 'PRESS' in self.ctx_data else 1013.25
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
            # J7, J3, J11, J12 are defined in anorg.mcm331.phot.rate as JSCALE*l*HV
            # (non-zero during day), so the ==0 guard would prevent the override —
            # all measured J overrides must be unconditional.
            # MCM J-numbering: J1=JO1D, J4=JNO2, J7=JHONO, J3=JH2O2, J11=JHCHO_R, J12=JHCHO_M
            if ctx.get('JO1D',   0.0) > 0: ctx['J1']  = ctx['JO1D']
            if ctx.get('JNO2',   0.0) > 0: ctx['J4']  = ctx['JNO2']
            if ctx.get('JHONO',  0.0) > 0: ctx['J7']  = ctx['JHONO']   # unconditional: beats JSCALE approx
            if ctx.get('JH2O2',  0.0) > 0: ctx['J3']  = ctx['JH2O2']   # unconditional: beats JSCALE approx
            if ctx.get('JHCHO_R',0.0) > 0: ctx['J11'] = ctx['JHCHO_R'] # unconditional: measured > PHTFKT by ~35%
            if ctx.get('JHCHO_M',0.0) > 0: ctx['J12'] = ctx['JHCHO_M'] # unconditional: measured > PHTFKT by ~27%

            # --- Recompute O1D PSS using final J1 (post-autolinker) ---
            if 'O1D' in _pss_species_set:
                _j1_now  = float(ctx.get('J1',  0.0))
                _o3_now  = float(ctx.get('O3',  0.0))
                _t_now   = float(ctx.get('T',   298.15))
                _n2_now  = float(ctx.get('N2',  0.7809 * M_now))
                _o2_now  = float(ctx.get('O2',  0.2095 * M_now))
                _h2o_now = float(ctx.get('H2O', 0.0))
                _kq_now  = (2.0e-11 * _n2_now * np.exp(130.0 / _t_now) +
                            3.20e-11 * _o2_now * np.exp(67.0  / _t_now) +
                            2.14e-10 * _h2o_now)
                ctx['O1D'] = max(_j1_now, 0.0) * max(_o3_now, 0.0) / max(_kq_now, 1e-30)

            # O atom is kept in the ODE state (not PSS) — see note in _build_mechanism.

            try:
                _compute_rates(ctx, _k_buf)
            except Exception:
                # Fallback: individual eval() per reaction (handles edge cases gracefully)
                for i, rxn in enumerate(self.mechanism.reactions):
                    try: _k_buf[i] = eval(rxn.safe_expr, globals(), ctx)
                    except: _k_buf[i] = 0.0
            np.nan_to_num(_k_buf, nan=0.0, posinf=0.0, neginf=0.0, copy=False)

            for i, rxn in enumerate(self.mechanism.reactions):
                k = float(_k_buf[i])
                if k != 0.0:
                    v = k
                    for r in rxn.reactants_parsed: v *= (y[r['idx']] ** r['coeff'])
                    for r in rxn.reactants_env: v *= (ctx[r['name']] ** r['coeff'])
                    for r in rxn.reactants_parsed: dydt[r['idx']] -= (v * r['coeff'])
                    for p in rxn.products_parsed: dydt[p['idx']] += (v * p['coeff'])

            # Analytical dilution: applied only when explicit +DIL reactions are
            # NOT already handling it. We detect this at setup via _has_explicit_dil.
            if self.cfg.get('analytical_dilution', True) and not self._has_explicit_dil:
                _flow = self.ctx_data['FLOW'](t_abs) if 'FLOW' in self.ctx_data else 0.0
                k_dil = (_flow + 0.6) / (self.cfg['chamber_vol'] * 3600)
                dydt -= y * k_dil

            # Thermal expansion is handled by instantaneous y rescaling at
            # segment boundaries (see below), not as a continuous ODE term.
            # With piecewise-constant T, dM/dt is zero within each segment
            # and a delta function at the boundary — a continuous correction
            # cannot capture this correctly.

            for s in cont_sources:
                dydt[s['idx']] += eval(s['code'], globals(), ctx)

            for b in burst_eval:
                if b['profile'] == 'tworeaction':
                    if b['start_s'] <= t_abs <= b['start_s'] + b['dur']:
                        dydt[b['idx1']] += b['rate1']
                        dydt[b['idx2']] += b['rate2']
                elif b['profile'] == 'saturating':
                    dt = t_abs - b['start_s']
                    if 0.0 <= dt <= b['dur']:
                        rate = b['coeff'] * np.exp(-b['lam'] * dt)
                        dydt[b['idx']] += rate
                        if b['byp_idx'] is not None:
                            dydt[b['byp_idx']] += rate * b['byp_frac']
                else:  # gaussian
                    rate = b['coeff'] * np.exp(-0.5 * ((t_abs - b['mu']) / b['sigma'])**2)
                    dydt[b['idx']] += rate
                    if b['byp_idx'] is not None:
                        dydt[b['byp_idx']] += rate * b['byp_frac']

            for imp in h2o_imps:
                if imp['start_s'] <= t_abs <= imp['start_s'] + imp['dur']:
                    dydt[imp['idx']] += imp['rate']

            self.eval_count += 1
            if self.eval_count % 500 == 0 or t > self.last_progress_print_time + 300:
                progress = max(0.0, min(t / _t_duration, 1.0))

                eta_str = ""
                if progress > 0:
                    elapsed = time.time() - self.real_start_time
                    total_est = elapsed / progress
                    rem_secs = int(total_est - elapsed)
                    eta_str = f" | ETA: {rem_secs // 60}m {rem_secs % 60}s"

                print(f"  ... ODE Eval Step: {self.eval_count} | Sim Time: {t:.1f}s / {_t_duration:.1f}s | Fraction Done: {progress:.3f} ({progress*100:.1f}%){eta_str}")
                if t > self.last_progress_print_time + 300:
                    self.last_progress_print_time = t

            return dydt

        solver_method = 'BDF'
        step_limit = min(min_dur / 2.0, 60.0)
        N = self.solver.num_species

        # --- PRE-COMPUTE ANALYTICAL JACOBIAN STRUCTURE ---
        # For each reaction i with ODE reactants r1..rN and affected species a1..aM:
        #   J[aj, rj] += (sign_aj * coeff_rj) * v_i / y[rj]
        # where v_i = k_i * prod(y_ode^c) * prod(env^c)  and  dv/dy[rj] = coeff_rj * v_i / y[rj]
        # Replaces expensive finite-difference (~200 ODE evals/Jacobian → 1 algebraic pass).
        _jac_rows_l, _jac_cols_l, _jac_factors_l, _jac_rxn_idxs_l = [], [], [], []
        for rxn_idx, rxn in enumerate(self.mechanism.reactions):
            if not rxn.reactants_parsed:
                continue  # zero-order or env-only reactants: dv/dy = 0
            # affected[i] = (species_idx, signed_stoich): reactants lose, products gain
            affected = [(r['idx'], -r['coeff']) for r in rxn.reactants_parsed] + \
                       [(p['idx'], +p['coeff']) for p in rxn.products_parsed]
            for rj in rxn.reactants_parsed:
                for ai, sign_ca in affected:
                    _jac_rows_l.append(ai)
                    _jac_cols_l.append(rj['idx'])
                    _jac_factors_l.append(sign_ca * rj['coeff'])
                    _jac_rxn_idxs_l.append(rxn_idx)
        _jac_rows     = np.array(_jac_rows_l,      dtype=np.int32)
        _jac_cols     = np.array(_jac_cols_l,      dtype=np.int32)
        _jac_factors  = np.array(_jac_factors_l,   dtype=np.float64)
        _jac_rxn_idxs = np.array(_jac_rxn_idxs_l, dtype=np.int32)
        print(f">>> [JAC] Analytical Jacobian: {N}×{N}, {len(_jac_rows)} terms "
              f"across {len(self.mechanism.reactions)} reactions")

        def analytical_jac(t, y_raw):
            """Analytical Jacobian: exact d(dydt)/dy computed from reaction structure.
            Much faster than finite-difference: O(reactions) vs O(N_colors * ODE_cost)."""
            y = np.maximum(y_raw, 0.0)
            t_abs = t + _t_offset   # convert relative solver time → absolute interpolation time
            M_now = self.get_M(t_abs)
            _temp_fn_j = self.ctx_data.get('TEMP') or self.ctx_data.get('T')
            T_now = _temp_fn_j(t_abs)
            ctx_j = SafeDict()
            # 1. Load all measured/environmental inputs
            for k, v_func in self.ctx_data.items():
                try:
                    ctx_j[k.upper()] = float(v_func(t_abs))
                except Exception:
                    pass
            # 2. H2O VMR → molec/cm³ (mirrors system_ode logic)
            if 'RH' in ctx_j and 'H2O' not in ctx_j:
                p_sat = np.exp(21.36469 - 5339.66 / T_now)
                press = float(self.ctx_data['PRESS'](t_abs)) if 'PRESS' in self.ctx_data else 1013.25
                ctx_j['H2O'] = (ctx_j['RH'] / 100.0) * (p_sat / press)
            if 'H2O' in ctx_j:
                h2o_vmr = ctx_j['H2O']
                ctx_j['H2O'] = h2o_vmr * M_now
                if 'RH' not in ctx_j:
                    if h2o_vmr > 0:
                        p_sat = np.exp(21.36469 - 5339.66 / T_now)
                        press = float(self.ctx_data['PRESS'](t_abs)) if 'PRESS' in self.ctx_data else 1013.25
                        ctx_j['RH'] = h2o_vmr * press / p_sat * 100.0
                    else:
                        ctx_j['RH'] = 0.0
            else:
                ctx_j['H2O'] = 0.01 * M_now
            # 3. Base physics variables + math functions
            ctx_j.update({
                "T": T_now, "TEMP": T_now, "M": M_now,
                "O2": 0.2095 * M_now, "N2": 0.7809 * M_now,
                "np": np, "exp": np.exp, "EXP": np.exp,
                "LOG10": np.log10, "log10": np.log10,
                "LOG": np.log, "log": np.log,
                "COS": np.cos, "SIN": np.sin, "TAN": np.tan, "SQRT": np.sqrt,
            })
            # 4. KMT functions
            for func_name, func in _kmt_funcs:
                ctx_j[func_name.upper()] = lambda *args, f=func, _T=T_now, _M=M_now: f(*args, _T, _M)
            for func_name, func in _other_funcs:
                ctx_j[func_name.upper()] = func
            # 5. ODE species concentrations
            ctx_j.update(zip(_species_names, y.tolist()))
            # 6. Pass 1: COMPRC/KMT var_assignments (COSX, KMT01-14, etc.)
            for var_name, comp_expr in self.var_assignments:
                try:
                    ctx_j[var_name] = eval(comp_expr, globals(), ctx_j)
                except Exception:
                    ctx_j[var_name] = 0.0
            # 7. Build PHTFKT callable (requires COSX from pass 1)
            _cosx  = max(float(ctx_j.get('COSX', np.cos(float(ctx_j.get('SZ', 45.0)) * 0.0174533))), 0.01)
            _jno2  = float(ctx_j.get('JNO2',  0.0))
            _lno2  = float(ctx_j.get('LNO2',  1.165e-2))
            _mno2  = float(ctx_j.get('MNO2',  0.244))
            _nnno2 = float(ctx_j.get('NNO2',  -0.267))
            _hv    = float(ctx_j.get('HV',    0.0))
            def _phtfkt_j(a, b, c, _c=_cosx, _j=_jno2, _l=_lno2, _m=_mno2, _n=_nnno2, _h=_hv):
                if _h == 0 or _c <= 0: return 0.0
                cs_j = a * _c**b * np.exp(c / _c)
                if _j > 0:
                    jno2_cs = _l * _c**_m * np.exp(_n / _c)
                    scale = min(_j / jno2_cs, 1.0) if jno2_cs > 0 else 1.0
                else:
                    scale = 1.0
                return cs_j * scale
            ctx_j['PHTFKT']    = _phtfkt_j
            ctx_j['PHTFKT_CS'] = lambda a, b, c, _c=_cosx: a * _c**b * np.exp(c / _c)
            ctx_j['AMAX']      = np.maximum
            ctx_j['AMIN']      = np.minimum
            # 8. Pass 2: J-value var_assignments (uses PHTFKT from step 7)
            for var_name, comp_expr in self.var_assignments:
                try:
                    ctx_j[var_name] = eval(comp_expr, globals(), ctx_j)
                except Exception:
                    ctx_j[var_name] = 0.0
            # 9. Measured J-value override
            if ctx_j.get('JO1D',   0.0) > 0: ctx_j['J1']  = ctx_j['JO1D']
            if ctx_j.get('JNO2',   0.0) > 0: ctx_j['J4']  = ctx_j['JNO2']
            if ctx_j.get('JHONO',  0.0) > 0: ctx_j['J7']  = ctx_j['JHONO']
            if ctx_j.get('JH2O2',  0.0) > 0: ctx_j['J3']  = ctx_j['JH2O2']
            if ctx_j.get('JHCHO_R', 0.0) > 0: ctx_j['J11'] = ctx_j['JHCHO_R']  # unconditional
            if ctx_j.get('JHCHO_M', 0.0) > 0: ctx_j['J12'] = ctx_j['JHCHO_M']  # unconditional
            # 9b. Recompute O1D PSS using final J1 (Jacobian must match system_ode)
            if 'O1D' in _pss_species_set:
                _j1j = float(ctx_j.get('J1', 0.0))
                _o3j = float(ctx_j.get('O3', 0.0))
                _tj  = float(ctx_j.get('T', 298.15))
                _n2j = float(ctx_j.get('N2', 0.7809 * M_now))
                _o2j = float(ctx_j.get('O2', 0.2095 * M_now))
                _hwj = float(ctx_j.get('H2O', 0.0))
                _kqj = (2.0e-11 * _n2j * np.exp(130.0 / _tj) +
                        3.20e-11 * _o2j * np.exp(67.0  / _tj) +
                        2.14e-10 * _hwj)
                ctx_j['O1D'] = max(_j1j, 0.0) * max(_o3j, 0.0) / max(_kqj, 1e-30)

            # 9c. O atom is kept in the ODE state (not PSS) — see note in _build_mechanism.

            # 10. Rate constants via precompiled batch function
            try:
                _compute_rates(ctx_j, _k_buf)
            except Exception:
                for i, rxn in enumerate(self.mechanism.reactions):
                    try: _k_buf[i] = eval(rxn.safe_expr, globals(), ctx_j)
                    except: _k_buf[i] = 0.0
            np.nan_to_num(_k_buf, nan=0.0, posinf=0.0, neginf=0.0, copy=False)
            # 11. Reaction velocities: v[i] = k[i] * prod(y_ode^c) * prod(env^c)
            np.copyto(_v_buf, _k_buf)
            for i, rxn in enumerate(self.mechanism.reactions):
                if _v_buf[i] == 0.0: continue
                for r in rxn.reactants_parsed:
                    _v_buf[i] *= y[r['idx']] ** r['coeff']
                for r in rxn.reactants_env:
                    env_v = float(ctx_j.get(r['name'], 0.0))
                    _v_buf[i] *= env_v ** r['coeff']
            # 12. Assemble sparse Jacobian
            # J[ai, rj] += factor * v[rxn] / y[rj]
            # where factor = sign_affected * coeff_reactant_rj (precomputed)
            y_col    = np.maximum(y[_jac_cols], 1e-40)
            jac_data = _jac_factors * _v_buf[_jac_rxn_idxs] / y_col
            # Dilution diagonal: d(dydt[i])/d(y[i]) -= k_dil (only in analytical mode)
            diag_rate = 0.0
            if self.cfg.get('analytical_dilution', True) and not self._has_explicit_dil:
                _flow    = self.ctx_data['FLOW'](t_abs) if 'FLOW' in self.ctx_data else 0.0
                diag_rate -= (_flow + 0.6) / (self.cfg['chamber_vol'] * 3600)
            if diag_rate != 0.0:
                dil_idx  = np.arange(N, dtype=np.int32)
                dil_data = np.full(N, diag_rate)
                rows = np.concatenate([_jac_rows, dil_idx])
                cols = np.concatenate([_jac_cols, dil_idx])
                d    = np.concatenate([jac_data, dil_data])
            else:
                rows, cols, d = _jac_rows, _jac_cols, jac_data
            # coo_matrix sums duplicate (row,col) entries — correct for multi-reaction contributions
            return coo_matrix((d, (rows, cols)), shape=(N, N)).tocsc()

        self.t_start = self.t_span[0]
        self.t_end = self.t_span[-1]
        _t_offset   = self.t_start          # solver works in relative time (t_rel = t_abs - offset)
        _t_duration = self.t_end - self.t_start  # total integration span in seconds
        self.last_progress_print_time = 0.0  # relative-time tracking
        self.eval_count = 0
        self.real_start_time = time.time()

        # --- DETECT DISCONTINUITY EVENTS IN ALL ENZ COLUMNS ---
        # All ENZ columns use piecewise-constant interpolation (kind='previous').
        # BDF polynomial history becomes inconsistent at step changes. We detect
        # value jumps across ALL columns and split integration at each one.
        event_times = []
        n_events_by_col = {}

        def _find_step_changes(func, label):
            """Scan a piecewise-constant signal for any value jumps."""
            probe_t = np.linspace(self.t_start, self.t_end, 5000)
            vals = np.array([float(func(t)) for t in probe_t])
            transitions = []
            for i in range(1, len(probe_t)):
                v0, v1 = vals[i - 1], vals[i]
                if v0 == v1:
                    continue
                rel = abs(v1 - v0) / max(abs(v0), abs(v1), 1e-30)
                if rel > 1e-6 or abs(v1 - v0) > 1e-10:
                    lo, hi = probe_t[i - 1], probe_t[i]
                    for _ in range(10):
                        mid = 0.5 * (lo + hi)
                        if float(func(lo)) != float(func(mid)):
                            hi = mid
                        else:
                            lo = mid
                    transitions.append(0.5 * (lo + hi))
            return transitions

        for key, func in self.ctx_data.items():
            if key.upper() == 'TIME':
                continue
            try:
                trans = _find_step_changes(func, key)
                if trans:
                    event_times.extend(trans)
                    n_events_by_col[key] = len(trans)
            except Exception:
                pass

        # Deduplicate events within 0.5s
        if event_times:
            event_times = sorted(set(event_times))
            deduped = [event_times[0]]
            for t in event_times[1:]:
                if t - deduped[-1] > 0.5:
                    deduped.append(t)
            event_times = deduped

        if event_times:
            summary = ', '.join(f'{k}:{v}' for k, v in sorted(n_events_by_col.items()))
            print(f">>> [SOLVER] {len(event_times)} discontinuity event(s) detected "
                  f"({summary}) — splitting integration for clean BDF restarts.")
        else:
            print(f">>> [SOLVER] No discontinuity events detected.")

        # Build list of (t0, t1) segment boundaries
        boundaries = sorted(set([self.t_start] + event_times + [self.t_end]))
        # Remove duplicate / too-close boundaries
        boundaries = [boundaries[0]] + [b for i, b in enumerate(boundaries[1:], 1)
                                         if b - boundaries[i - 1] > 1.0]

        print(f">>> [DEBUG] Solving {len(self.mechanism.species)} species with {solver_method} "
              f"in {len(boundaries)-1} segment(s) (max_step={step_limit}s)...")

        # Solve each segment, passing the final state as y0 for the next
        seg_t_list, seg_y_list = [], []
        y_seg = y0.copy()
        for seg_idx, (t0, t1) in enumerate(zip(boundaries[:-1], boundaries[1:])):
            # Extract t_eval points that fall within this segment (absolute time)
            mask = (self.t_span >= t0) & (self.t_span <= t1)
            t_eval_seg = self.t_span[mask]
            # Always include the segment endpoints so the final y is exact
            if len(t_eval_seg) == 0 or t_eval_seg[0] > t0 + 1.0:
                t_eval_seg = np.concatenate([[t0], t_eval_seg])
            if t_eval_seg[-1] < t1 - 1.0:
                t_eval_seg = np.concatenate([t_eval_seg, [t1]])
            t_eval_seg = np.unique(t_eval_seg)

            # Shift to relative time before passing to solver (avoids float64 precision loss
            # when t_abs ~ 7.7e8 s and the required step ~ 1e-4 s: eps*t_abs ~ 1.7e-7 s)
            t0_rel      = t0 - _t_offset
            t1_rel      = t1 - _t_offset
            t_eval_rel  = t_eval_seg - _t_offset

            seg_sol = solve_ivp(
                system_ode,
                [t0_rel, t1_rel],
                y_seg,
                method=solver_method,
                t_eval=t_eval_rel,
                max_step=step_limit,
                rtol=1e-2,
                atol=1e4,
                jac=analytical_jac,
            )

            if not seg_sol.success:
                print(f"\n[WARNING] Segment {seg_idx+1}/{len(boundaries)-1} "
                      f"({t0:.0f}–{t1:.0f}s) stopped early: {seg_sol.message}")

            # Collect results (skip first point of subsequent segments to avoid duplicates)
            # Restore absolute time for the output arrays
            skip = 1 if seg_idx > 0 and len(seg_t_list) > 0 else 0
            seg_t_list.append(seg_sol.t[skip:] + _t_offset)
            seg_y_list.append(seg_sol.y[:, skip:])
            y_seg = seg_sol.y[:, -1]

            # Thermal expansion: rescale concentrations at segment boundary
            # to conserve mixing ratios when M changes due to T step.
            # With piecewise-constant T, M is constant within each segment
            # but may jump at the boundary. FACSIMILE does this implicitly.
            if seg_idx < len(boundaries) - 2:
                t_next = boundaries[seg_idx + 1]
                M_end  = self.get_M(t1)        # M at end of this segment
                M_next = self.get_M(t_next + 0.01)  # M at start of next segment
                if M_end > 0 and abs(M_next / M_end - 1.0) > 1e-10:
                    y_seg = y_seg * (M_next / M_end)

        # Stitch segments into a single solution object
        all_t = np.concatenate(seg_t_list)
        all_y = np.concatenate(seg_y_list, axis=1)

        # Wrap in a namespace that _export_results expects (sol.t and sol.y)
        class _SolResult:
            pass
        self.sol = _SolResult()
        self.sol.t = all_t
        self.sol.y = all_y
        self.sol.success = True

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

        # --- O atom export: disabled (was post-process PSS overwrite) ---------------
        # FACSIMILE keeps O in the ODE state and reports values O(~10²) molec/cm³,
        # which is ~100× smaller than the analytical PSS J4·NO2/k_O+O2+M (~10⁴).
        # The PSS overwrite produced spurious +10880% bias in the diagnostic.
        # Per commit 5da582f, JUICE keeps O in the ODE — leave the raw ODE state in
        # the export so comparisons reflect what the solver actually computes.
        if False and 'O' in df.columns:
            _T   = (df['SAPHIR_TEMP'].values if 'SAPHIR_TEMP' in df.columns else
                    df['SAPHIR_T'].values     if 'SAPHIR_T'    in df.columns else
                    np.full(len(df), 298.15))
            _M   = df['M'].values
            _HV  = df['SAPHIR_HV'].values if 'SAPHIR_HV' in df.columns else np.ones(len(df))
            _no2 = df['NO2'].values        if 'NO2'       in df.columns else np.zeros(len(df))
            _o3  = df['O3'].values         if 'O3'        in df.columns else np.zeros(len(df))
            _no3 = df['NO3'].values        if 'NO3'       in df.columns else np.zeros(len(df))
            _j4  = (df['SAPHIR_JNO2'].values if 'SAPHIR_JNO2' in df.columns else
                    df['J4'].values           if 'J4'          in df.columns else
                    np.zeros(len(df)))
            _j2  = df['J2'].values if 'J2' in df.columns else np.zeros(len(df))   # O3 + hv → O
            _j6  = df['J6'].values if 'J6' in df.columns else np.zeros(len(df))   # NO3 + hv → NO2 + O
            _j4  = np.where(_HV > 0.5, _j4, 0.0)
            _j2  = np.where(_HV > 0.5, _j2, 0.0)
            _j6  = np.where(_HV > 0.5, _j6, 0.0)
            _N2  = 0.7809 * _M
            _O2  = 0.2095 * _M
            _ks  = (5.6e-34 * _N2 + 6.0e-34 * _O2) * ((_T / 300.0) ** -2.6) * _O2
            _o_prod = (np.maximum(_j4, 0.0) * np.maximum(_no2, 0.0) +
                       np.maximum(_j2, 0.0) * np.maximum(_o3,  0.0) +
                       np.maximum(_j6, 0.0) * np.maximum(_no3, 0.0))
            df['O'] = np.where(_ks > 0, _o_prod / np.maximum(_ks, 1e-30), 0.0)
            print("  [O]   Replaced ODE output with pseudo-steady-state (J4·NO2+J2·O3+J6·NO3)/k_O+O2+M")

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

        # Export reaction kinetics (k-values + turnovers) to separate files
        self._export_kinetics(df)

        return out_path

    def _export_kinetics(self, df):
        """Re-evaluate all reaction rate constants and turnovers at every output timestep.

        Design for speed:
        - All env time series pre-extracted as numpy arrays (one vectorized call each).
        - Species matrix pre-extracted as numpy (n_t × n_species) — no pandas per-step.
        - Rate-constant loop: only unavoidable per-timestep work (scalar ctx + compiled eval).
        - Turnover computation: fully vectorized in numpy AFTER the k loop — no Python
          loop over (timestep × reaction).

        Saves:
        - simulation_kinetics.npz  : compressed float32 arrays (n_t × n_rxn)
            keys: 'time_s', 'k_values', 'turnovers' (molec/cm³/s)
        - simulation_reactions.csv : human-readable reaction manifest
        M is already in simulation_results.csv.
        """
        print("\n>>> [DATA] Computing reaction kinetics (k-values + turnovers)...")
        n_t   = len(self.sol.t)
        n_rxn = len(self.mechanism.reactions)

        _species_names = list(self.solver.species)
        _compute_rates = self._compute_rates_fn
        _kmt_funcs     = self._kmt_funcs_cached
        _other_funcs   = self._other_funcs_cached

        # ── PRE-EXTRACT: species concentrations as contiguous numpy matrix ─────────
        # df[col].iloc[t_i] in a loop is ~100× slower than numpy indexing.
        # Only keep species that are actually ODE columns in df.
        _sp_cols    = [s for s in _species_names if s in df.columns]
        _sp_missing = [s for s in _species_names if s not in df.columns]
        y_matrix = np.maximum(df[_sp_cols].values.astype(np.float64), 0.0)  # (n_t, len(_sp_cols))
        _sp_idx_map = {s: i for i, s in enumerate(_sp_cols)}  # species → column in y_matrix

        # ── PRE-EXTRACT: env time series by calling interpolators vectorized ───────
        env_arrs = {}   # key.upper() → numpy array shape (n_t,)
        for k_env, v_func in self.ctx_data.items():
            try:
                arr = v_func(self.sol.t)
                if not isinstance(arr, np.ndarray):
                    arr = np.full(n_t, float(arr))
                env_arrs[k_env.upper()] = arr.astype(np.float64)
            except Exception:
                pass

        # Pre-extract M and T arrays (used every step)
        M_arr = df['M'].values.astype(np.float64) if 'M' in df.columns else \
                np.array([self.get_M(t) for t in self.sol.t])
        T_arr = env_arrs.get('TEMP', np.full(n_t, 298.15))

        # H2O VMR → molec/cm³ array (pre-converted once)
        if 'H2O' in env_arrs:
            h2o_arr = env_arrs['H2O'] * M_arr          # VMR fraction × M
        else:
            h2o_arr = 0.01 * M_arr                     # 1% fallback

        # ── RATE-CONSTANT LOOP (per-timestep, scalar ctx) ─────────────────────────
        # This cannot be fully vectorized because the compiled rate expressions
        # contain scalar math (exp, KMT functions, etc.).  Everything else is
        # optimised: numpy indexing instead of pandas, pre-extracted arrays instead
        # of repeated interpolation calls.
        k_arr  = np.zeros((n_t, n_rxn), dtype=np.float64)
        _k_loc = np.zeros(n_rxn)

        for t_i in range(n_t):
            t_abs = self.sol.t[t_i]
            M_now = M_arr[t_i]
            T_now = T_arr[t_i]

            ctx = SafeDict()
            # Load env from pre-extracted arrays (fast numpy scalar)
            for k_up, arr in env_arrs.items():
                ctx[k_up] = float(arr[t_i])
            # H2O already converted
            ctx['H2O'] = float(h2o_arr[t_i])
            ctx.update({"T": T_now, "TEMP": T_now, "M": M_now,
                        "O2": 0.2095 * M_now, "N2": 0.7809 * M_now,
                        "np": np, "exp": np.exp, "EXP": np.exp,
                        "LOG10": np.log10, "log10": np.log10,
                        "LOG": np.log, "log": np.log, "SQRT": np.sqrt,
                        "COS": np.cos, "SIN": np.sin, "TAN": np.tan,
                        "AMAX": np.maximum, "AMIN": np.minimum})
            for fn, func in _kmt_funcs:
                ctx[fn.upper()] = lambda *args, f=func, _T=T_now, _M=M_now: f(*args, _T, _M)
            for fn, func in _other_funcs:
                ctx[fn.upper()] = func
            # Species from numpy (no pandas .iloc)
            for j, s in enumerate(_sp_cols):
                ctx[s] = float(y_matrix[t_i, j])
            # Pass 1: constants + KMT rates
            for var_name, comp_expr in self.var_assignments:
                try: ctx[var_name] = eval(comp_expr, globals(), ctx)
                except: ctx[var_name] = 0.0
            # PHTFKT callable
            _cosx = max(float(ctx.get('COSX', np.cos(float(ctx.get('SZ', 45.0)) * 0.0174533))), 0.01)
            _jno2 = float(ctx.get('JNO2', 0.0))
            _lno2 = float(ctx.get('LNO2', 1.165e-2))
            _mno2 = float(ctx.get('MNO2', 0.244))
            _nnno2= float(ctx.get('NNO2', -0.267))
            _hv   = float(ctx.get('HV', 0.0))
            def _phtfkt(a, b, c, _c=_cosx, _j=_jno2, _l=_lno2, _m=_mno2, _n=_nnno2, _h=_hv):
                if _h == 0 or _c <= 0: return 0.0
                cs_j = a * _c**b * np.exp(c / _c)
                if _j > 0:
                    jcs = _l * _c**_m * np.exp(_n / _c)
                    return cs_j * min(_j / jcs, 1.0) if jcs > 0 else cs_j
                return 0.0
            ctx['PHTFKT']    = _phtfkt
            ctx['PHTFKT_CS'] = lambda a, b, c, _c=_cosx: a * _c**b * np.exp(c / _c)
            # Pass 2: J-values using PHTFKT
            for var_name, comp_expr in self.var_assignments:
                try: ctx[var_name] = eval(comp_expr, globals(), ctx)
                except: pass
            # Measured J overrides
            if ctx.get('JO1D', 0) > 0: ctx['J1'] = ctx['JO1D']
            if ctx.get('JNO2', 0) > 0: ctx['J4'] = ctx['JNO2']
            if ctx.get('JHONO', 0) > 0: ctx['J7'] = ctx['JHONO']
            if ctx.get('JH2O2', 0) > 0: ctx['J3'] = ctx['JH2O2']
            if ctx.get('JHCHO_R', 0) > 0: ctx['J11'] = ctx['JHCHO_R']  # unconditional
            if ctx.get('JHCHO_M', 0) > 0: ctx['J12'] = ctx['JHCHO_M']  # unconditional

            try:
                _compute_rates(ctx, _k_loc)
            except Exception:
                for i, rxn in enumerate(self.mechanism.reactions):
                    try: _k_loc[i] = eval(rxn.safe_expr, globals(), ctx)
                    except: _k_loc[i] = 0.0
            np.nan_to_num(_k_loc, nan=0.0, posinf=0.0, neginf=0.0, copy=False)
            k_arr[t_i] = _k_loc

            if (t_i + 1) % 400 == 0 or t_i == n_t - 1:
                print(f"  ... k-values {t_i + 1}/{n_t}")

        # ── VECTORISED TURNOVER COMPUTATION ───────────────────────────────────────
        # v[t,i] = k[t,i] * prod_j( y[t,j]^c_j ) * prod_e( env[t,e]^c_e )
        # All multiplications are numpy array ops over the time axis — no Python loop
        # over timesteps.  Only the outer loop over reactions (n_rxn) remains.
        print("  ... computing turnovers (vectorised)...")
        v_arr = k_arr.copy()   # shape (n_t, n_rxn), start with k

        for i, rxn in enumerate(self.mechanism.reactions):
            # ODE species reactants
            for r in rxn.reactants_parsed:
                col_i = _sp_idx_map.get(_species_names[r['idx']])
                if col_i is None:
                    v_arr[:, i] = 0.0
                    break
                y_col = y_matrix[:, col_i]           # shape (n_t,)
                if r['coeff'] == 1.0:
                    v_arr[:, i] *= y_col
                else:
                    v_arr[:, i] *= y_col ** r['coeff']
            # Env reactants (e.g. H2O, M)
            for r in rxn.reactants_env:
                env_col = env_arrs.get(r['name'])
                if env_col is None:
                    # Special cases: M, O2, N2 not in env_arrs but derivable
                    if r['name'] == 'M':
                        env_col = M_arr
                    elif r['name'] == 'O2':
                        env_col = 0.2095 * M_arr
                    elif r['name'] == 'N2':
                        env_col = 0.7809 * M_arr
                    elif r['name'] == 'H2O':
                        env_col = h2o_arr
                    else:
                        v_arr[:, i] = 0.0
                        break
                if r['coeff'] == 1.0:
                    v_arr[:, i] *= np.maximum(env_col, 0.0)
                else:
                    v_arr[:, i] *= np.maximum(env_col, 0.0) ** r['coeff']

        np.nan_to_num(v_arr, nan=0.0, posinf=0.0, neginf=0.0, copy=False)

        # ── BUILD REACTION MANIFEST ───────────────────────────────────────────────
        rxn_labels    = []
        rows_manifest = []
        for i, rxn in enumerate(self.mechanism.reactions):
            r_str  = ' + '.join(rxn.reactants) if rxn.reactants else '(none)'
            p_str  = ' + '.join(rxn.products)  if rxn.products  else '(none)'
            rate_e = getattr(rxn, 'rate_expression', getattr(rxn, 'rate_str', '?'))
            label  = f"R{i:05d}"
            rxn_labels.append(label)
            rows_manifest.append({"index": i, "label": label,
                                   "reactants": r_str, "products": p_str,
                                   "rate_expression": rate_e[:120]})

        # ── SAVE ─────────────────────────────────────────────────────────────────
        npz_path = os.path.join(self.root, "simulation_kinetics.npz")
        np.savez_compressed(npz_path,
                            time_s=self.sol.t.astype(np.float64),
                            k_values=k_arr.astype(np.float32),
                            turnovers=v_arr.astype(np.float32),
                            reaction_labels=np.array(rxn_labels))
        print(f"  ... Kinetics saved to {npz_path}  "
              f"({n_rxn} reactions × {n_t} timesteps, float32)")

        man_path = os.path.join(self.root, "simulation_reactions.csv")
        pd.DataFrame(rows_manifest).to_csv(man_path, index=False)
        print(f"  ... Reaction manifest saved to {man_path}")

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