import os
import json
import numpy as np
import pandas as pd 
from scipy.integrate import solve_ivp
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
            print("\n>>> Starting SAPHIR Simulation Pipeline")
            print("-" * 40)
            self._setup_environment()
            self._build_mechanism()
            self._solve()
            out_path = self._export_results()
            print("-" * 40)
            print(">>> All steps completed successfully.")
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
        
        active_parameters = self.cfg.get('active_parameters')
        
        for f_name in mech_files:
            f_path = os.path.join(chem_path, f_name)
            if os.path.exists(f_path):
                m_part = parser.parse(f_path)
                self.mechanism.species.extend(m_part.species)
                self.mechanism.reactions.extend(m_part.reactions)
                
                # Process both active and commented parameters from the parser
                all_params = {**m_part.parameters, **m_part.commented_parameters}

                for p_name, p_expr in all_params.items():
                    # If active_parameters list exists, only include parameters from that list
                    if active_parameters is not None and p_name.upper() not in active_parameters:
                        continue

                    p_clean = p_expr.replace('@', '**')
                    p_clean = re.sub(r'(?<=\d)[dD](?=[+-]?\d)', 'E', p_clean)
                    p_clean = re.sub(r'J\((\d+)\)', r'J\1', p_clean)
                    try:
                        self.var_assignments.append((p_name.upper(), compile(p_clean, '<string>', 'eval')))
                    except Exception as e:
                        print(f"[WARNING] Could not compile parameter {p_name}: {e}")
                    
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

        # --- BUG FIX: Only purge explicitly mapped GUI variables, leave O3 alone! ---
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
        actual_start_utc = saphir_base_date + pd.to_timedelta(self.t_span[0], unit='s')

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
                coeff = (float(b['ppb']) * 1e-9 * self.get_M(start_s)) / (sigma * np.sqrt(2 * np.pi))
                burst_eval.append({'idx': self.solver.species_to_idx[target], 'mu': mu, 'sigma': sigma, 'coeff': coeff})

        def system_ode(t, y):
            y = np.maximum(y, 0.0) 
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
                
            if 'H2O' in ctx and ctx['H2O'] > 0:
                if ctx['H2O'] > 1.0: ctx['H2O'] = ctx['H2O'] / 100.0 # Fix percentage data
                ctx['H2O'] = ctx['H2O'] * M_now
            else:
                ctx['H2O'] = 0.01 * M_now # Safe 1% Fallback

            # --- THE PHOTOLYSIS AUTO-LINKER ---
            if 'J41' not in ctx or ctx['J41'] == 0.0:
                if 'JO1D' in ctx: ctx['J41'] = ctx['JO1D']
            if 'J1' not in ctx or ctx['J1'] == 0.0:
                if 'JNO2' in ctx: ctx['J1'] = ctx['JNO2']
            if 'J3' not in ctx or ctx['J3'] == 0.0:
                if 'JHONO' in ctx: ctx['J3'] = ctx['JHONO']

            ctx.update({
                "T": T_now, "TEMP": T_now, "M": M_now, 
                "O2": 0.2095 * M_now, 
                "N2": 0.7809 * M_now, 
                "np": np, "exp": np.exp, "EXP": np.exp, 
                "LOG10": np.log10, "log10": np.log10,
                "LOG": np.log, "log": np.log,
                "COS": np.cos, "SIN": np.sin, "TAN": np.tan, "SQRT": np.sqrt
            })
            
            for func_name in dir(mcm_kinetics):
                if not func_name.startswith('_'):
                    func = getattr(mcm_kinetics, func_name)
                    if func_name.startswith('KMT'):
                        ctx[func_name.upper()] = lambda *args, f=func: f(*args, T_now, M_now)
                    else:
                        ctx[func_name.upper()] = func

            for i, name in enumerate(self.solver.species):
                ctx[name] = float(y[i])

            for _ in range(2):
                for var_name, comp_expr in self.var_assignments:
                    if var_name in active_constraints:
                        continue 
                    try: 
                        ctx[var_name] = eval(comp_expr, globals(), ctx)
                    except:
                        ctx[var_name] = 0.0

            for rxn in self.mechanism.reactions:
                try: 
                    k = eval(rxn.safe_expr, globals(), ctx)
                except:
                    k = 0.0

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
                dydt[b['idx']] += b['coeff'] * np.exp(-0.5 * ((t - b['mu']) / b['sigma'])**2)
            
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

        solver_method = 'BDF'
        step_limit = min(min_dur / 2.0, 60.0)
        
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
            atol=1e4
        )
        if not self.sol.success: print(f"\n[WARNING] ODE Solver stopped early!\nReason: {self.sol.message}")

    def _export_results(self):
        print("\n>>> [DATA] Exporting results to CSV...")
        
        # Transpose sol.y so rows are time steps and columns are species
        df = pd.DataFrame(self.sol.y.T, columns=self.solver.species, index=self.sol.t)
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

        out_path = os.path.join(self.root, "simulation_results.csv")
        df.to_csv(out_path)
        
        print(f"  ... Saved to {out_path}")
        return out_path

    def _show_results(self, out_path):
        try:
            from PyQt6.QtWidgets import QApplication
            has_pyqt = True
        except ImportError:
            has_pyqt = False
            print("[WARNING] PyQt6 is not installed. Falling back to simple Matplotlib visualization.")
            
        if has_pyqt:
            print(">>> [VISUALIZATION] Launching Interactive PyQt Viewer...")
            df = pd.read_csv(out_path, index_col=0)
            app = QApplication.instance()
            if app is None:
                app = QApplication(sys.argv)
            viewer = ResultsViewer(df)
            viewer.show()
            app.exec()
            return
            
        # Fallback to simple generic plot if PyQt is missing
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            return
        
        df = pd.read_csv(out_path, index_col=0)
        t_hours = (df.index - df.index[0]) / 3600.0
        
        fig, ax1 = plt.subplots(figsize=(12, 6))
        ax1.set_xlabel('Time (Hours)')
        ax1.set_ylabel('Mixing Ratio (ppb)', color='black')
        
        macro_species = ['O3', 'NO', 'NO2', 'HONO', 'C2H5OH', 'HCHO']
        for spec in macro_species:
            if spec in df.columns:
                ppb_vals = (df[spec] / df['M']) * 1e9
                ax1.plot(t_hours, ppb_vals, label=f"{spec} (ppb)")
                
        ax1.grid(True, linestyle='--', alpha=0.6)
        
        radical_species = ['OH', 'HO2', 'RO2']
        if any(r in df.columns for r in radical_species):
            ax2 = ax1.twinx()
            ax2.set_ylabel('Radicals (molec/cm³)', color='red')
            for r in radical_species:
                if r in df.columns:
                    ax2.plot(t_hours, df[r], label=f"{r} (molec/cm³)", linestyle='--')
            
        fig.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05), ncol=5)
        plt.title('SAPHIR Scenario Simulation')
        plt.tight_layout()
        plt.show()

try:
    from PyQt6.QtWidgets import (QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, 
                                 QListWidget, QPushButton, QLabel, QComboBox, 
                                 QSplitter, QLineEdit)
    from PyQt6.QtCore import Qt
    from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
    from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
    from matplotlib.figure import Figure
except ImportError:
    pass

class ResultsViewer(QMainWindow):
    def __init__(self, df):
        super().__init__()
        self.df = df.copy()
        self.setWindowTitle("Interactive SAPHIR/MCM Results Viewer")
        self.setGeometry(50, 50, 1500, 850)

        # Convert index to Hours for X-Axis
        self.df.index = (self.df.index - self.df.index[0]) / 3600.0
        
        central = QWidget()
        self.setCentralWidget(central)
        main_layout = QHBoxLayout(central)
        splitter = QSplitter(Qt.Orientation.Horizontal)
        main_layout.addWidget(splitter)
        
        # --- LEFT PANEL (Controls) ---
        left_panel = QWidget()
        left_layout = QVBoxLayout(left_panel)
        
        self.search_box = QLineEdit()
        self.search_box.setPlaceholderText("Search variable...")
        self.search_box.textChanged.connect(self.filter_list)
        left_layout.addWidget(self.search_box)
        
        self.var_list = QListWidget()
        self.var_list.setSelectionMode(QListWidget.SelectionMode.ExtendedSelection)
        self.var_list.itemDoubleClicked.connect(self.add_selected)
        self.all_vars = sorted([c for c in self.df.columns if c != 'M'])
        self.var_list.addItems(self.all_vars)
        left_layout.addWidget(self.var_list)
        
        control_layout = QVBoxLayout()
        control_layout.addWidget(QLabel("Target Axis:"))
        self.axis_combo = QComboBox()
        self.axis_combo.addItems(["Left Axis (Main)", "Right Axis (Main)", "Left Axis (Offset)", "Right Axis (Offset)"])
        control_layout.addWidget(self.axis_combo)
        
        control_layout.addWidget(QLabel("Conversion (MCM Data):"))
        self.unit_combo = QComboBox()
        self.unit_combo.addItems(["Convert to ppb", "Native (molec/cm³ or SAPHIR unit)"])
        control_layout.addWidget(self.unit_combo)
        
        btn_add = QPushButton("Plot Selected (Or Double-Click List)")
        btn_add.clicked.connect(self.add_selected)
        control_layout.addWidget(btn_add)
        
        btn_clear = QPushButton("Clear Canvas")
        btn_clear.clicked.connect(self.clear_plot)
        control_layout.addWidget(btn_clear)
        
        left_layout.addLayout(control_layout)
        splitter.addWidget(left_panel)
        
        # --- RIGHT PANEL (Canvas) ---
        right_panel = QWidget()
        right_layout = QVBoxLayout(right_panel)
        
        self.fig = Figure(figsize=(10, 6))
        self.canvas = FigureCanvas(self.fig)
        self.toolbar = NavigationToolbar(self.canvas, self)
        
        right_layout.addWidget(self.toolbar)
        right_layout.addWidget(self.canvas)
        splitter.addWidget(right_panel)
        splitter.setSizes([300, 1200])
        
        self.init_axes()
        self._plot_defaults()

    def filter_list(self, text):
        self.var_list.clear()
        self.var_list.addItems([v for v in self.all_vars if text.lower() in v.lower()])

    def init_axes(self):
        self.fig.clear()
        self.ax1 = self.fig.add_subplot(111)
        self.ax1.set_xlabel("Time (Hours)", fontweight='bold')
        self.ax1.grid(True, linestyle='--', alpha=0.4)
        self.ax2 = self.ax1.twinx()
        self.ax3 = self.ax1.twinx()
        self.ax3.spines['right'].set_position(('outward', 60))
        self.ax4 = self.ax1.twinx()
        self.ax4.spines['left'].set_position(('outward', 60))
        self.ax4.yaxis.set_ticks_position('left')
        self.ax4.yaxis.set_label_position('left')
        
        self.axes_map = {
            "Left Axis (Main)": self.ax1, "Right Axis (Main)": self.ax2,
            "Right Axis (Offset)": self.ax3, "Left Axis (Offset)": self.ax4
        }
        self.fig.subplots_adjust(left=0.2, right=0.8, top=0.85)

    def add_selected(self):
        selected = [item.text() for item in self.var_list.selectedItems()]
        if not selected and self.var_list.currentItem():
            selected = [self.var_list.currentItem().text()]
            
        ax_name = self.axis_combo.currentText()
        unit = self.unit_combo.currentText()
        
        for var in selected:
            self.plot_variable(var, ax_name, unit)
        self._update_legends()
        self.canvas.draw()

    def clear_plot(self):
        self.init_axes()
        self.canvas.draw()

    def plot_variable(self, var_name, axis_name, unit):
        ax = self.axes_map[axis_name]
        data = self.df[var_name].copy()
        
        if "ppb" in unit and not var_name.startswith("SAPHIR_") and 'M' in self.df.columns:
            data = (data / self.df['M']) * 1e9
            label = f"{var_name} (ppb)"
        else:
            label = var_name
            
        ax.plot(self.df.index, data, label=label, linewidth=2.0)
        if not ax.get_ylabel(): ax.set_ylabel(axis_name, fontweight='bold')

    def _update_legends(self):
        lines, labels = [], []
        for ax in self.axes_map.values():
            h, l = ax.get_legend_handles_labels()
            lines.extend(h); labels.extend(l)
        if self.fig.legends: self.fig.legends[0].remove()
        if lines: self.fig.legend(lines, labels, loc='upper center', bbox_to_anchor=(0.5, 0.98), ncol=6)

    def _plot_defaults(self):
        if 'C2H5OH' in self.df.columns: self.plot_variable('C2H5OH', "Left Axis (Main)", "Convert to ppb")
        if 'NO' in self.df.columns: self.plot_variable('NO', "Right Axis (Main)", "Convert to ppb")
        if 'NO2' in self.df.columns: self.plot_variable('NO2', "Right Axis (Main)", "Convert to ppb")
        if 'O3' in self.df.columns: self.plot_variable('O3', "Right Axis (Offset)", "Convert to ppb")
        if 'OH' in self.df.columns: self.plot_variable('OH', "Left Axis (Offset)", "Native")
        if 'HO2' in self.df.columns: self.plot_variable('HO2', "Left Axis (Offset)", "Native")
        self._update_legends()
        self.canvas.draw()

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python scenario_runner.py <path_to_scenario_config.json>")
        sys.exit(1)
        
    config_path = sys.argv[1]
    runner = ScenarioRunner(config_path)
    runner.run()