import os
import json
import numpy as np
import pandas as pd 
from scipy.integrate import solve_ivp
import sys
import re
import warnings

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
        self.debug_counter = 0
        self.debug_counter = 0

    def run(self):
        try:
            print("\n>>> Starting SAPHIR Simulation Pipeline")
            print("-" * 40)
            self._setup_environment()
            self._build_mechanism()
            self._solve()
            out_path = self._export_results()
            print("-" * 40)
            print(">>> All steps completed successfully.")
            # self._show_results(out_path)
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

        for internal_key, mapping in self.cfg['env_map'].items():
            if mapping['type'] == "SAPHIR File":
                file_header = mapping['value'].lower()
                if file_header in providers:
                    self.ctx_data[internal_key] = providers[file_header]
                else:
                    actual_header = next((h for h in headers if h.lower() == file_header), None)
                    if actual_header: self.ctx_data[internal_key] = providers[actual_header]
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
                
                for p_name, p_expr in m_part.parameters.items():
                    p_clean = p_expr.replace('@', '**')
                    p_clean = re.sub(r'(?<=\d)[dD](?=[+-]?\d)', 'E', p_clean)
                    p_clean = re.sub(r'J\((\d+)\)', r'J\1', p_clean)
                    try: self.var_assignments.append((p_name.upper(), compile(p_clean, '<string>', 'eval')))
                    except: pass
                    
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

        for mapped_key in self.cfg['env_map'].keys():
            all_species.discard(mapped_key.upper())
        for bg in ['RH', 'O2', 'N2', 'M']:
            all_species.discard(bg)

        self.mechanism.species = sorted(list(all_species))
        self._add_artifacts()

    def _add_artifacts(self):
        if self.cfg['artifacts']['active']:
            scale = eval(str(self.cfg['artifacts']['scale']))
            if self.cfg['artifacts'].get('hono_active', True):
                self.mechanism.reactions.append(Reaction([], ["HONO"], f"({self.cfg['artifacts']['hono_eqn']}) * {scale}"))
            if self.cfg['artifacts'].get('hcho_active', True):
                self.mechanism.reactions.append(Reaction([], ["HCHO"], f"{self.cfg['artifacts']['hcho_eqn']}"))
            if self.cfg['artifacts'].get('y_active', True):
                if 'Y' not in self.mechanism.species: self.mechanism.species.append('Y')
                y_rate = f"({self.cfg['artifacts']['y_eqn']}) * 1e-9 * M * 2.4e-13" 
                self.mechanism.reactions.append(Reaction(['OH'], ["HO2"], y_rate))
            if self.cfg['artifacts'].get('wall_active', True):
                self.mechanism.reactions.append(Reaction(["H2O"], ["HONO"], str(self.cfg['artifacts']['dark_hono'])))

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
                        test_ctx[func_name.upper()] = func
            
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
                    failed_rxns.append(f"'{expr}' ({e})")
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
                print(f"[WARNING] {len(failed_rxns)} reactions failed to evaluate and will be DEAD! Details: {failed_rxns}")
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
        for b in self.cfg.get('burst_injections', []):
            if not b.get('target') or not b.get('ppb') or b['ppb'] == "": continue
            target = b['target'].upper()
            if target in self.solver.species_to_idx:
                target_utc = pd.to_datetime(b['start_utc'], dayfirst=True)
                start_s = (target_utc - saphir_base_date).total_seconds()
                dur = float(b.get('duration_s', 60))
                min_dur = min(min_dur, dur)
                mu = start_s + (dur / 2.0)
                sigma = dur / 6.0 
                # Use M from the start of the simulation to prevent extrapolation errors
                coeff = (float(b['ppb']) * 1e-9 * self.get_M(self.t_span[0])) / (sigma * np.sqrt(2 * np.pi))
                burst_eval.append({'idx': self.solver.species_to_idx[target], 'mu': mu, 'sigma': sigma, 'coeff': coeff})

        def system_ode(t, y):
            y = np.maximum(y, 0.0) 
            dydt = np.zeros(self.solver.num_species)
            M_now = self.get_M(t)
            T_now = self.ctx_data['TEMP'](t)

            ctx = SafeDict()
            co2_idx = self.solver.species_to_idx.get('CO2', -1)

            # 1. Baseline Environment Data (Loads everything from .enz)
            active_constraints = set()
            for k, v_func in self.ctx_data.items():
                try: 
                    ctx[k.upper()] = float(v_func(t))
                    if k.upper() in self.cfg['env_map']:
                        active_constraints.add(k.upper())
                except Exception:
                    pass

            # --- CORRECTED WATER AND RH LOGIC ---
            if 'H2O' in self.ctx_data:
                h2o_ppb = float(self.ctx_data['H2O'](t))
                ctx['H2O'] = (h2o_ppb * 1e-9) * M_now
                if 'PRESS' in self.ctx_data:
                    press_hpa = float(self.ctx_data['PRESS'](t))
                    # Magnus formula for saturation vapor pressure in hPa
                    p_sat_hpa = 6.112 * np.exp((17.67 * (T_now - 273.15)) / (T_now - 29.65))
                    p_h2o_hpa = (h2o_ppb * 1e-9) * press_hpa
                    ctx['RH'] = (p_h2o_hpa / p_sat_hpa) * 100.0
                else:
                    ctx['RH'] = 50.0 # Fallback
            else:
                ctx['H2O'] = 0.01 * M_now # 1% Fallback
                ctx['RH'] = 50.0

            # --- THE PHOTOLYSIS AUTO-LINKER ---
            if 'J41' not in ctx or ctx['J41'] == 0.0:
                if 'JO1D' in ctx: ctx['J41'] = ctx['JO1D']
            if 'J1' not in ctx or ctx['J1'] == 0.0:
                if 'JNO2' in ctx: ctx['J1'] = ctx['JNO2']
            if 'J3' not in ctx or ctx['J3'] == 0.0:
                if 'JHONO' in ctx: ctx['J3'] = ctx['JHONO']

            # 2. Math & Atmospheric Bulk Gases
            ctx.update({
                "T": T_now, "TEMP": T_now, "M": M_now, 
                "O2": 0.2095 * M_now, "N2": 0.7809 * M_now, "np": np, 
                "exp": np.exp, "EXP": np.exp, "LOG10": np.log10, "log10": np.log10,
                "LOG": np.log, "log": np.log, "COS": np.cos, "SIN": np.sin, "TAN": np.tan, "SQRT": np.sqrt
            })
            for func_name in dir(mcm_kinetics):
                if not func_name.startswith('_'):
                    func = getattr(mcm_kinetics, func_name)
                    if func_name.startswith('KMT'):
                        ctx[func_name.upper()] = lambda *args, f=func: f(*args, T_now, M_now)
                    else:
                        ctx[func_name.upper()] = func

            # 3. Inject Simulated Species
            for i, name in enumerate(self.solver.species):
                ctx[name] = float(y[i])

            # 4. Variable Assignments
            for _ in range(2):
                for var_name, comp_expr in self.var_assignments:
                    if var_name in active_constraints: continue 
                    try: ctx[var_name] = eval(comp_expr, globals(), ctx)
                    except: ctx[var_name] = 0.0

            # 5. Chemistry Kinetics
            for rxn in self.mechanism.reactions:
                try: 
                    k = eval(rxn.safe_expr, globals(), ctx)
                except:
                    k = 0.0

                if k > 0.0:
                    v = k
                    for r in rxn.reactants_parsed: v *= (y[r['idx']] ** r['coeff'])
                    for r in rxn.reactants_env: 
                        if r['name'] in ctx: v *= (ctx[r['name']] ** r['coeff'])
                        else: v = 0.0; break
                    
                    if v > 0.0:
                        dydt_before = dydt[co2_idx] if co2_idx != -1 else 0
                        
                        for r in rxn.reactants_parsed: dydt[r['idx']] -= (v * r['coeff'])
                        for p in rxn.products_parsed: dydt[p['idx']] += (v * p['coeff'])
                        
                        dydt_after = dydt[co2_idx] if co2_idx != -1 else 0
                        if self.debug_counter < 5 and abs(dydt_before - dydt_after) > 1e-10:
                            print(f"\n!!! ROGUE REACTION altering CO2 gradient !!!")
                            print(f"  Reaction: {' + '.join(rxn.reactants)} -> {' + '.join(rxn.products)}")
                            print(f"  Rate Expression: {getattr(rxn, 'rate_expression', 'N/A')}")
                            print(f"  Velocity v = {v:.4e}")

            # 6. Physical Chamber Dynamics
            k_dil = (self.ctx_data['FLOW'](t) + 0.6) / (self.cfg['chamber_vol'] * 3600)
            dydt -= y * k_dil

            for s in cont_sources:
                dydt[s['idx']] += eval(s['code'], globals(), ctx)

            for b in burst_eval:
                dydt[b['idx']] += b['coeff'] * np.exp(-0.5 * ((t - b['mu']) / b['sigma'])**2)
            
            self.debug_counter += 1
            return dydt

        step_limit = min(min_dur / 4.0, 15.0)
        print(f">>> [DEBUG] Solving for {len(self.mechanism.species)} species using BDF (Max Step: {step_limit}s)...")
        k_dil_initial = (self.ctx_data['FLOW'](self.t_span[0]) + 0.6) / (self.cfg['chamber_vol'] * 3600)
        print(f">>> [DEBUG] Initial Dilution Rate (k_dil) = {k_dil_initial:.4e} s-1")
        self.sol = solve_ivp(system_ode, [self.t_span[0], self.t_span[-1]], y0, method='BDF', t_eval=self.t_span, max_step=step_limit)
        if not self.sol.success: print(f"\n⚠️ [WARNING] ODE Solver stopped early!\nReason: {self.sol.message}")

    def _export_results(self):
        output_path = os.path.join(self.root, "model_output.csv")
        out_molec = pd.DataFrame(self.sol.y.T, index=self.sol.t, columns=self.solver.species)
        M_array = np.array([self.get_M(t) for t in self.sol.t])
        out_ppb = out_molec.div(M_array * 1e-9, axis=0) 
        
        # --- THE RADICAL VISIBILITY EXPORTER ---
        for rad in ['OH', 'O1D']:
            if rad in out_molec.columns: out_ppb[f"{rad}_molec"] = out_molec[rad] # Molecules/cm3
        for rad in ['HO2', 'CH3O2']:
            if rad in out_molec.columns: out_ppb[f"{rad}_ppt"] = out_molec[rad] / (M_array * 1e-12) # ppt
        
        for k, func in self.ctx_data.items():
            try: out_ppb[f"ENV_{k}"] = [float(func(t_val)) for t_val in self.sol.t]
            except: pass
            
        out_ppb["INJECTION_EVENT"] = ""
        saphir_base_date = pd.Timestamp("2000-01-01 00:00:00")
        actual_start_utc = saphir_base_date + pd.to_timedelta(self.t_span[0], unit='s')
        col_idx = out_ppb.columns.get_loc("INJECTION_EVENT")
        
        for b in self.cfg.get('burst_injections', []):
            if not b.get('target') or not b.get('ppb') or b['ppb'] == "": continue
            user_dt = pd.to_datetime(b['start_utc'], dayfirst=True)
            target_utc = user_dt.replace(year=actual_start_utc.year, month=actual_start_utc.month, day=actual_start_utc.day)
            start_s = (target_utc - saphir_base_date).total_seconds()
            
            idx = np.argmin(np.abs(self.sol.t - start_s))
            existing = out_ppb.iat[idx, col_idx]
            evt_str = f"Injected {b['ppb']}ppb {b['target']} @ {b['start_utc']}"
            out_ppb.iat[idx, col_idx] = (existing + " | " + evt_str) if existing else evt_str

        out_ppb.to_csv(output_path)
        print(f">>> [DEBUG] Results (including ENV observations) saved to: {output_path}")
        return output_path

    def _show_results(self, output_path):
        print(">>> [DEBUG] Launching Interactive Results Viewer...")
        import tkinter as tk
        from tkinter import ttk
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
        import matplotlib.pyplot as plt

        df = pd.read_csv(output_path, index_col=0)
        time_hours = (df.index - df.index[0]) / 3600.0
        species_list = [c for c in df.columns if c != "INJECTION_EVENT"]

        root = tk.Tk()
        root.title("SAPHIR Simulation Results Viewer")
        root.geometry("1000x600")

        left_frame = ttk.Frame(root, width=200, padding=10)
        left_frame.pack(side=tk.LEFT, fill=tk.Y)

        ttk.Label(left_frame, text="Search Species:").pack(side=tk.TOP, anchor='w')
        search_var = tk.StringVar()
        search_entry = ttk.Entry(left_frame, textvariable=search_var)
        search_entry.pack(side=tk.TOP, fill=tk.X, pady=(0, 5))

        btn_frame = ttk.Frame(left_frame)
        btn_frame.pack(side=tk.TOP, fill=tk.X, pady=(0, 5))
        ttk.Label(btn_frame, text="Actions (Or Double-Click list):").pack(side=tk.TOP, anchor='w')
        ttk.Button(btn_frame, text="Replace Plot (Left Axis)", command=lambda: plot_action('replace')).pack(side=tk.TOP, fill=tk.X, pady=1)
        ttk.Button(btn_frame, text="+ Append Left Axis", command=lambda: plot_action('add_left')).pack(side=tk.TOP, fill=tk.X, pady=1)
        ttk.Button(btn_frame, text="+ Append Right Axis", command=lambda: plot_action('add_right')).pack(side=tk.TOP, fill=tk.X, pady=1)
        ttk.Button(btn_frame, text="Clear Plot", command=lambda: plot_action('clear')).pack(side=tk.TOP, fill=tk.X, pady=1)

        list_frame = ttk.Frame(left_frame)
        list_frame.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        scrollbar = ttk.Scrollbar(list_frame)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        listbox = tk.Listbox(list_frame, yscrollcommand=scrollbar.set, exportselection=False)
        listbox.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.config(command=listbox.yview)

        def update_list(*args):
            search_term = search_var.get().lower()
            listbox.delete(0, tk.END)
            for sp in species_list:
                if search_term in sp.lower(): listbox.insert(tk.END, sp)

        search_var.trace_add("write", update_list)
        update_list()

        right_frame = ttk.Frame(root)
        right_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)
        fig, ax = plt.subplots(figsize=(6, 4), dpi=100)
        canvas = FigureCanvasTkAgg(fig, master=right_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        toolbar = NavigationToolbar2Tk(canvas, right_frame)
        toolbar.update()

        plot_state = {'left': [], 'right': []}

        def redraw():
            ax.clear()
            if hasattr(ax, 'right_ax'):
                ax.right_ax.remove()
                del ax.right_ax

            lines = []
            labels = []
            colors = plt.cm.tab10.colors
            c_idx = 0

            if plot_state['left']:
                for sp in plot_state['left']:
                    l, = ax.plot(time_hours, df[sp], label=sp, color=colors[c_idx % 10], linewidth=2)
                    lines.append(l)
                    labels.append(sp)
                    c_idx += 1
                ax.set_ylabel("Measured Value" if any(s.startswith("ENV_") for s in plot_state['left']) else "Concentration")
            else:
                ax.set_ylabel("")

            if plot_state['right']:
                ax.right_ax = ax.twinx()
                for sp in plot_state['right']:
                    l, = ax.right_ax.plot(time_hours, df[sp], label=sp+' (R)', color=colors[c_idx % 10], linewidth=2, linestyle='--')
                    lines.append(l)
                    labels.append(sp+' (R)')
                    c_idx += 1
                ax.right_ax.set_ylabel("Measured Value" if any(s.startswith("ENV_") for s in plot_state['right']) else "Concentration")

            ax.set_xlabel("Time (Hours Since Start)")
            ax.grid(True, linestyle='--', alpha=0.7)

            if lines: ax.legend(lines, labels, loc='best')

            if "INJECTION_EVENT" in df.columns:
                inj_mask = df["INJECTION_EVENT"].fillna("").astype(str) != ""
                trans = ax.get_xaxis_transform()
                for inj_time_s, event_str in df.loc[inj_mask, "INJECTION_EVENT"].items():
                    h = (inj_time_s - df.index[0]) / 3600.0
                    ax.axvline(x=h, color='k', linestyle=':', alpha=0.6)
                    tag = str(event_str).split(' @ ')[0] 
                    ax.text(h, 1.02, f" {tag}", transform=trans, rotation=45, verticalalignment='bottom', horizontalalignment='left', alpha=0.8, fontsize=9)
                    
            ax.set_title("Simulation Results", fontsize=12, fontweight='bold', pad=60)
            fig.tight_layout()
            canvas.draw()

        def plot_action(action):
            selection = listbox.curselection()
            species = listbox.get(selection[0]) if selection else None
            
            if action == 'clear': plot_state['left'] = []; plot_state['right'] = []
            elif not species: return
            elif action == 'replace': plot_state['left'] = [species]; plot_state['right'] = []
            elif action == 'add_left' and species not in plot_state['left']: plot_state['left'].append(species)
            elif action == 'add_right' and species not in plot_state['right']: plot_state['right'].append(species)
            redraw()

        listbox.bind('<Double-1>', lambda e: plot_action('replace'))
        default_sp = "C2H5OH" if "C2H5OH" in species_list else (species_list[0] if species_list else None)
        if default_sp:
            try:
                idx = listbox.get(0, tk.END).index(default_sp)
                listbox.selection_set(idx)
                listbox.see(idx)
                plot_action('replace')
            except ValueError: pass

        root.mainloop()

if __name__ == "__main__":
    config_file = sys.argv[1] if len(sys.argv) > 1 else "scenario_config.json"
    try:
        runner = ScenarioRunner(config_file)
        runner.run()
    except Exception as e:
        print(f"[ERROR] [STARTUP ERROR] Could not initialize runner: {e}")