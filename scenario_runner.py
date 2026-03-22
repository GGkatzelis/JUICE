import os
import json
import numpy as np
import pandas as pd 
from scipy.integrate import solve_ivp
import sys
import re

# Internal Module Imports
from mechanism_solver import MechanismSolver
from easy_parser import EasyParser
from mechanism import Mechanism, Reaction
from load_saphir_environment import load_saphir_environment
import mcm_kinetics

class SafeDict(dict):
    """A bulletproof dictionary for chemical math. If a variable is missing, it safely returns 0.0."""
    def __missing__(self, key):
        return 0.0

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

    def run(self):
        try:
            print("\n🚀 Starting SAPHIR Simulation Pipeline")
            print("-" * 40)
            self._setup_environment()
            self._build_mechanism()
            self._solve()
            out_path = self._export_results()
            print("-" * 40)
            print("✨ All steps completed successfully.")
            self._show_results(out_path)
        except Exception as e:
            print(f"\n❌ [CRITICAL ERROR] Simulation failed!")
            print(f"Error Details: {str(e)}")
            import traceback
            traceback.print_exc()

    def _setup_environment(self):
        enz_file = next(f for f in os.listdir(self.root) if 'kammer' in f.lower() and f.endswith('.enz'))
        enz_path = os.path.join(self.root, enz_file)
        providers, t_min, t_max, headers = load_saphir_environment(enz_path)
        
        self.t_span = np.linspace(t_min, t_max, num=2000)

        # 1. AUTO-LOAD ALL RAW VARIABLES (Crucial for measured JHONO, JHCHO, etc.)
        for raw_header, interp_func in providers.items():
            self.ctx_data[raw_header.upper()] = interp_func

        # 2. APPLY GUI MAPPINGS
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
        
        # Parse standard MCM definitions
        for f_name in self.cfg['mechanism_files']:
            f_path = os.path.join(chem_path, f_name)
            if os.path.exists(f_path):
                m_part = parser.parse(f_path)
                self.mechanism.species.extend(m_part.species)
                self.mechanism.reactions.extend(m_part.reactions)
                
                # --- Map Variable Assignments & Parameters ---
                for p_name, p_expr in m_part.parameters.items():
                    p_clean = p_expr.replace('@', '**')
                    p_clean = re.sub(r'(?<=\d)[dD](?=[+-]?\d)', 'E', p_clean)
                    p_clean = re.sub(r'J\((\d+)\)', r'J\1', p_clean)
                    try: self.var_assignments.append((p_name.upper(), compile(p_clean, '<string>', 'eval')))
                    except: pass
                    
        # Inject Derived J-Scalings
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
            
        # Compile Custom Reactions
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
        actual_start_utc = saphir_base_date + pd.to_timedelta(self.t_span[0], unit='s')

        for spec, ppb in self.cfg.get('initials', {}).items():
            if not ppb: continue
            spec_up = spec.upper()
            if spec_up in self.solver.species_to_idx:
                y0[self.solver.species_to_idx[spec_up]] = (float(ppb) * 1e-9) * self.get_M(self.t_span[0])

        # --- Diagnostic Pre-Check ---
        print("\n>>> [DIAGNOSTICS] Testing chemical mechanism evaluation...")
        test_ctx = SafeDict()
        test_ctx.update({"T": 298.15, "TEMP": 298.15, "M": 2.5e19, "RH": 50.0, "np": np, "exp": np.exp, "EXP": np.exp, "LOG10": np.log10, "log10": np.log10})
        test_ctx['KMT06'] = lambda A, B, C, D: mcm_kinetics.KMT06(A, B, C, D, test_ctx['TEMP'], test_ctx['M'])
        test_ctx['KMT07'] = lambda A, B, C, D, E, F: mcm_kinetics.KMT07(A, B, C, D, E, F, test_ctx['TEMP'], test_ctx['M'])
        
        for k in self.ctx_data.keys(): test_ctx[k.upper()] = 1.0
        for s in self.solver.species: test_ctx[s] = 1.0
        for _ in range(2):
            for var_name, comp_expr in self.var_assignments:
                try: test_ctx[var_name] = eval(comp_expr, {"__builtins__": None}, test_ctx)
                except: pass

        # Pre-compile Reactions
        failed_rxns = []
        for rxn in self.mechanism.reactions:
            expr = getattr(rxn, 'rate_expression', getattr(rxn, 'rate_str', '0.0')).upper()
            expr = expr.replace('@', '**') 
            expr = re.sub(r'(?<=\d)[dD](?=[+-]?\d)', 'E', expr)
            expr = re.sub(r'J\((\d+)\)', r'J\1', expr)
            try: 
                rxn.safe_expr = compile(expr, '<string>', 'eval')
                eval(rxn.safe_expr, {"__builtins__": None}, test_ctx)
            except Exception as e: 
                failed_rxns.append((expr, str(e)))
                rxn.safe_expr = compile('0.0', '<string>', 'eval')

            rxn.reactants_parsed = []
            for r in rxn.reactants:
                coeff, name = 1.0, r.strip().upper()
                if '*' in name: 
                    coeff, name = float(name.split('*', 1)[0]), name.split('*', 1)[1]
                if name in self.solver.species_to_idx:
                    rxn.reactants_parsed.append({'idx': self.solver.species_to_idx[name], 'coeff': coeff})
                    
            rxn.products_parsed = []
            for p in rxn.products:
                coeff, name = 1.0, p.strip().upper()
                if '*' in name: 
                    coeff, name = float(name.split('*', 1)[0]), name.split('*', 1)[1]
                if name in self.solver.species_to_idx:
                    rxn.products_parsed.append({'idx': self.solver.species_to_idx[name], 'coeff': coeff})

        if failed_rxns:
            print(f"⚠️ [WARNING] {len(failed_rxns)} reactions failed to evaluate and will be DEAD (k=0.0)!")
            for expr, err in failed_rxns[:5]: print(f"   -> Eq: {expr} | Error: {err}")
            if len(failed_rxns) > 5: print(f"   -> ...and {len(failed_rxns)-5} more.")
        else:
            print(">>> [DIAGNOSTICS] All chemical equations evaluated successfully!")

        # Pre-compile Continuous Sources
        cont_sources = []
        for s in self.cfg.get('continuous_sources', []):
            if not s.get('target') or not s.get('rate'): continue
            target = s['target'].upper()
            if target in self.solver.species_to_idx:
                cont_sources.append({'idx': self.solver.species_to_idx[target], 'code': compile(s['rate'], '<string>', 'eval')})

        # Process Burst Injections
        burst_eval = []
        min_dur = 3600.0
        for b in self.cfg.get('burst_injections', []):
            if not b.get('target') or not b.get('ppb') or b['ppb'] == "": continue
            target = b['target'].upper()
            if target in self.solver.species_to_idx:
                user_dt = pd.to_datetime(b['start_utc'], dayfirst=True)
                target_utc = user_dt.replace(year=actual_start_utc.year, month=actual_start_utc.month, day=actual_start_utc.day)
                start_s = (target_utc - saphir_base_date).total_seconds()
                
                dur = float(b.get('duration_s', 60))
                min_dur = min(min_dur, dur)
                
                mu = start_s + (dur / 2.0)
                sigma = dur / 6.0 
                coeff = (float(b['ppb']) * 1e-9 * self.get_M(start_s)) / (sigma * np.sqrt(2 * np.pi))
                burst_eval.append({'idx': self.solver.species_to_idx[target], 'mu': mu, 'sigma': sigma, 'coeff': coeff})

        def system_ode(t, y):
            dydt = np.zeros(self.solver.num_species)
            M_now = self.get_M(t)
            T_now = self.ctx_data['TEMP'](t)

            ctx = SafeDict()
            ctx.update({
                "T": T_now, "TEMP": T_now, "M": M_now, 
                "np": np, "exp": np.exp, "EXP": np.exp, 
                "LOG10": np.log10, "log10": np.log10,
                "KMT06": lambda A, B, C, D: mcm_kinetics.KMT06(A, B, C, D, T_now, M_now),
                "KMT07": lambda A, B, C, D, E, F: mcm_kinetics.KMT07(A, B, C, D, E, F, T_now, M_now)
            })
            
            for func_name in dir(mcm_kinetics):
                if not func_name.startswith('_') and func_name not in ['KMT06', 'KMT07']: 
                    ctx[func_name.upper()] = getattr(mcm_kinetics, func_name)

            for k, v_func in self.ctx_data.items():
                try: ctx[k.upper()] = float(v_func(t))
                except: ctx[k.upper()] = 0.0

            if 'H2O' in self.ctx_data and 'PRESS' in self.ctx_data:
                p_sat = np.exp(21.36469 - 5339.66 / T_now)
                ctx['RH'] = (self.ctx_data['H2O'](t) * 100) / (p_sat / self.ctx_data['PRESS'](t) * M_now)

            # 1. Inject Species into Math Context
            for i, name in enumerate(self.solver.species):
                ctx[name] = float(y[i])

            # 2. Evaluate Variable Assignments (2 passes to cascade RO2 sums & J-Values)
            for _ in range(2):
                for var_name, comp_expr in self.var_assignments:
                    try: ctx[var_name] = eval(comp_expr, {"__builtins__": None}, ctx)
                    except: pass

            # 3. Chemistry Processing
            for rxn in self.mechanism.reactions:
                try: k = eval(rxn.safe_expr, {"__builtins__": None}, ctx)
                except: k = 0.0

                if k != 0.0:
                    v = k
                    for r in rxn.reactants_parsed: v *= (y[r['idx']] ** r['coeff'])
                    for r in rxn.reactants_parsed: dydt[r['idx']] -= (v * r['coeff'])
                    for p in rxn.products_parsed: dydt[p['idx']] += (v * p['coeff'])

            # 4. Dilution
            k_dil = (self.ctx_data['FLOW'](t) + 0.6) / (self.cfg['chamber_vol'] * 3600)
            dydt -= y * k_dil

            for s in cont_sources:
                dydt[s['idx']] += eval(s['code'], ctx)

            for b in burst_eval:
                dydt[b['idx']] += b['coeff'] * np.exp(-0.5 * ((t - b['mu']) / b['sigma'])**2)
                
            return dydt

        step_limit = min(min_dur / 2.0, 60.0)
        print(f">>> [DEBUG] Solving for {len(self.mechanism.species)} species using BDF (Max Step: {step_limit}s)...")
        self.sol = solve_ivp(system_ode, [self.t_span[0], self.t_span[-1]], y0, method='BDF', t_eval=self.t_span, max_step=step_limit)
        if not self.sol.success: print(f"\n⚠️ [WARNING] ODE Solver stopped early!\nReason: {self.sol.message}")

    def _export_results(self):
        output_path = os.path.join(self.root, "model_output.csv")
        out_molec = pd.DataFrame(self.sol.y.T, index=self.sol.t, columns=self.solver.species)
        M_array = np.array([self.get_M(t) for t in self.sol.t])
        out_ppb = out_molec.div(M_array * 1e-9, axis=0) # Export back into easy PPB
        
        # --- NEW: Append Kammer Data & Mapped Environmental Variables ---
        for k, func in self.ctx_data.items():
            try: out_ppb[f"ENV_{k}"] = [float(func(t_val)) for t_val in self.sol.t]
            except: pass
            
        # --- NEW: Record Injection Events as a Timeseries ---
        out_ppb["INJECTION_EVENT"] = ""
        saphir_base_date = pd.Timestamp("2000-01-01 00:00:00")
        actual_start_utc = saphir_base_date + pd.to_timedelta(self.t_span[0], unit='s')
        col_idx = out_ppb.columns.get_loc("INJECTION_EVENT")
        
        for b in self.cfg.get('burst_injections', []):
            if not b.get('target') or not b.get('ppb') or b['ppb'] == "": continue
            user_dt = pd.to_datetime(b['start_utc'], dayfirst=True)
            target_utc = user_dt.replace(year=actual_start_utc.year, month=actual_start_utc.month, day=actual_start_utc.day)
            start_s = (target_utc - saphir_base_date).total_seconds()
            
            # Find the closest time step in the ODE output to attach the string
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
        # Hide the text event column from the plotting selector
        species_list = [c for c in df.columns if c != "INJECTION_EVENT"]

        root = tk.Tk()
        root.title("SAPHIR Simulation Results Viewer")
        root.geometry("1000x600")

        # Left sidebar for list and search
        left_frame = ttk.Frame(root, width=200, padding=10)
        left_frame.pack(side=tk.LEFT, fill=tk.Y)

        ttk.Label(left_frame, text="Search Species:").pack(side=tk.TOP, anchor='w')
        search_var = tk.StringVar()
        search_entry = ttk.Entry(left_frame, textvariable=search_var)
        search_entry.pack(side=tk.TOP, fill=tk.X, pady=(0, 5))

        # Plotting Controls
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

        # Right area for plotting
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

            # Left Axis
            if plot_state['left']:
                for sp in plot_state['left']:
                    l, = ax.plot(time_hours, df[sp], label=sp, color=colors[c_idx % 10], linewidth=2)
                    lines.append(l)
                    labels.append(sp)
                    c_idx += 1
                ax.set_ylabel("Measured Value" if any(s.startswith("ENV_") for s in plot_state['left']) else "Concentration (ppb) [Left]")
            else:
                ax.set_ylabel("")

            # Right Axis
            if plot_state['right']:
                ax.right_ax = ax.twinx()
                for sp in plot_state['right']:
                    l, = ax.right_ax.plot(time_hours, df[sp], label=sp+' (R)', color=colors[c_idx % 10], linewidth=2, linestyle='--')
                    lines.append(l)
                    labels.append(sp+' (R)')
                    c_idx += 1
                ax.right_ax.set_ylabel("Measured Value" if any(s.startswith("ENV_") for s in plot_state['right']) else "Concentration (ppb) [Right]")

            ax.set_xlabel("Time (Hours Since Start)")
            ax.grid(True, linestyle='--', alpha=0.7)

            if lines:
                ax.legend(lines, labels, loc='best')

            # Injections
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
            
            if action == 'clear':
                plot_state['left'] = []; plot_state['right'] = []
            elif not species: return
            elif action == 'replace':
                plot_state['left'] = [species]; plot_state['right'] = []
            elif action == 'add_left' and species not in plot_state['left']: plot_state['left'].append(species)
            elif action == 'add_right' and species not in plot_state['right']: plot_state['right'].append(species)
            redraw()

        listbox.bind('<Double-1>', lambda e: plot_action('replace'))

        # Set default selection to Ethanol if it exists
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
        print(f"❌ [STARTUP ERROR] Could not initialize runner: {e}")