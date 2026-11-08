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
import datetime

try:
    from version import __version__ as JUICE_VERSION
except ImportError:
    JUICE_VERSION = "unknown"

# ── Numba @njit stoichiometry kernel ────────────────────────────────────────
# Defined at module level so Numba's cache=True persists across _solve() calls.
# Replaces the per-reaction Python loop in system_ode: computes reaction
# velocities from k_buf + y + env_vals, then scatters dydt contributions.
try:
    import numba as _numba

    @_numba.njit(cache=True)
    def _stoich_eval_numba(k_buf, y, env_vals,
                           ode_r_rxn, ode_r_sp, ode_r_coeff,
                           env_r_rxn, env_r_name_idx, env_r_coeff,
                           stoich_rxn, stoich_sp, stoich_sc,
                           v_buf, dydt_out):
        """In-place kernel: fills dydt_out (pre-zeroed) and uses v_buf as scratch.
        No internal array allocations — caller owns both buffers."""
        for i in range(len(k_buf)):
            v_buf[i] = k_buf[i]
        for j in range(len(ode_r_rxn)):
            i = ode_r_rxn[j]
            if v_buf[i] != 0.0:
                v_buf[i] *= y[ode_r_sp[j]] ** ode_r_coeff[j]
        for j in range(len(env_r_rxn)):
            i = env_r_rxn[j]
            if v_buf[i] != 0.0:
                ev = env_vals[env_r_name_idx[j]]
                v_buf[i] *= ev ** env_r_coeff[j]
        for j in range(len(stoich_rxn)):
            dydt_out[stoich_sp[j]] += v_buf[stoich_rxn[j]] * stoich_sc[j]

    _NUMBA_AVAILABLE = True
except ImportError:
    _NUMBA_AVAILABLE = False

    def _stoich_eval_numba(*_a, **_kw):  # noqa: F811
        raise RuntimeError("numba not installed")

_FAC_AMP = "&" * 44   # separator line before data (≥8 & required)
_FAC_TIME_COL = "TIME / SEC SINCE 1.1.2000 00:00:00"


def _col_to_fac_header(col: str, col_idx: int, file_type: str) -> str:
    """Map a DataFrame column name to a FACSIMILE-style COLUMN N= header string."""
    if col_idx == 1:  # first column is always TIME
        return _FAC_TIME_COL
    if file_type == 'kv':
        m = re.match(r'^([kv])_R(\d+):(.+)$', col)
        if m:
            return f"{m.group(1)}<{int(m.group(2))},0> !U{m.group(3)}!N"
        return f"{col}<0>"
    if file_type == 'pd':
        if col.startswith('P_'):
            return f"p<0> {col[2:]}"
        if col.startswith('D_'):
            return f"d<0> {col[2:]}"
        return f"{col}<0>"
    return f"{col}<0>"   # cs: every column gets <0>


def _write_csv(path: str, df: pd.DataFrame, comment: str = "",
               file_type: str = 'cs') -> None:
    """Write a FACSIMILE-compatible output file (COLUMN headers + &&&& + space-sep data).

    Format:
      # comment
      COLUMN 1=TIME / SEC SINCE 1.1.2000 00:00:00
      COLUMN 2=OH<0>          (cs)  or  k<1,0> !U...!N  (kv)  or  p<0> SP  (pd)
      ...
      NUMBER OF COLUMNS=N
      &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
       value value ...         ← space-separated data rows
    """
    all_cols = [df.index.name or 'TIME_s'] + list(df.columns)
    idx_vals = df.index.values.astype(float)
    col_arrs = [df[c].values.astype(float) for c in df.columns]

    with open(path, 'w', newline='', encoding='utf-8') as fh:
        if comment:
            fh.write(f"# {comment}\n")
        for i, col in enumerate(all_cols, 1):
            fh.write(f"COLUMN {i}={_col_to_fac_header(col, i, file_type)}\n")
        fh.write(f"NUMBER OF COLUMNS={len(all_cols)}\n")
        fh.write(_FAC_AMP + "\n")
        # Time column with full float64 precision (.15E = 16 sig figs) so that
        # BDF boundary times like 773812899.9999847 are never rounded to the same
        # value as the exact boundary 773812900.0 (relative diff ~2e-14 needs >=14 sig figs).
        # Formatting is pure Python (no C fast path exists for '%.6E'); large files
        # (kv: 27 M cells) are formatted in parallel row chunks, byte-identical.
        data = np.column_stack(col_arrs) if col_arrs else np.empty((len(idx_vals), 0))
        n_cells = data.size
        n_workers = max(1, min(8, (os.cpu_count() or 2) - 1))
        if n_cells >= 3_000_000 and n_workers > 1:
            from concurrent.futures import ProcessPoolExecutor
            n_chunks = n_workers * 2
            bounds = np.linspace(0, len(idx_vals), n_chunks + 1).astype(int)
            jobs = [(idx_vals[a:b], data[a:b]) for a, b in zip(bounds[:-1], bounds[1:]) if b > a]
            try:
                with ProcessPoolExecutor(max_workers=n_workers) as ex:
                    chunks = list(ex.map(_format_rows, jobs))
            except Exception as _e:          # e.g. BrokenProcessPool (frozen/odd main)
                print(f"  [EXPORT] parallel formatting unavailable ({type(_e).__name__}); serial")
                chunks = [_format_rows(j) for j in jobs]
            for chunk in chunks:
                fh.write(chunk)
        else:
            fh.write(_format_rows((idx_vals, data)))


def _format_rows(job):
    """Format (times, 2-D data) into FACSIMILE data lines: time .15E, values .6E."""
    idx_vals, data = job
    out = []
    for t_i in range(len(idx_vals)):
        row = f'{idx_vals[t_i]:.15E}'
        for v in data[t_i]:
            row += f' {v:.6E}'
        out.append(row)
    return '\n'.join(out) + ('\n' if out else '')


def _juice_header(scenario_stem: str, backend: str = "python") -> str:
    """Return the standard JUICE file header comment."""
    ts = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    return (f"Calculated by JUICE v{JUICE_VERSION}  "
            f"Backend: {backend}  "
            f"Scenario: {scenario_stem}  "
            f"Date: {ts}")

def _load_step_col_patterns():
    """Return step-column patterns: hardcoded defaults + TOML additions.

    Layer 1 (always active): Q*  HV
    Layer 2 (from FR_GUI.toml [solver] step_column_patterns_extra): added on top
    Layer 3 (future): EASY-script comment annotations — see
             docs/plan_easy_step_col_annotation.md
    """
    defaults = ['Q*', 'HV']
    toml_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "FR_GUI.toml")
    try:
        import tomllib
        with open(toml_path, 'rb') as fh:
            cfg = tomllib.load(fh)
        extra = cfg.get('solver', {}).get('step_column_patterns_extra', [])
        if isinstance(extra, str):
            import ast as _ast
            try: extra = _ast.literal_eval(extra)
            except Exception: extra = []
        # Merge defaults + extras, preserving order and removing duplicates
        seen = set(p.upper() for p in defaults)
        combined = list(defaults)
        for p in extra:
            if p.upper() not in seen:
                combined.append(p)
                seen.add(p.upper())
        return combined
    except Exception:
        return defaults

# Internal Module Imports
from mechanism_solver import MechanismSolver
from easy_parser import EasyParser
from mechanism import Mechanism, Reaction
from load_saphir_environment import load_saphir_environment
# mcm_kinetics is intentionally NOT imported into the execution path.
# All KMT/COMPRC/PHTFKT functions are inline-expanded by easy_parser.expand_macros
# at parse time — the EASY script is the single authoritative source.

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
            
        self.root = self.cfg.get('folder', '') or os.path.dirname(os.path.abspath(config_path))
        # Derive output file stem from the JSON filename (e.g. scene_03.json → "scene_03")
        self._stem = os.path.splitext(os.path.basename(config_path))[0]
        self.mechanism = None
        self.solver = None
        self.t_span = None
        self.ctx_data = {}
        self.var_assignments = []
        self._ctx_data_protected = set()   # raw ENZ aliases (e.g. H2O_IN) — never overwritten by var_assignments
        self._idl_params_upper  = set()    # idl_parameter names — allowed to override env_map
        self.last_progress_print_time = -1
        self.t_start = 0
        self.t_end = 0
        self.eval_count = 0
        self.real_start_time = 0.0
        self._o1d_prod_rxns = []   # populated in _solve() after mechanism is built
        self._o1d_loss_rxns = []
        self._sp_by_idx = {}       # idx → species name

    def _pss_o1d(self, ctx):
        """Compute [O1D] dynamically from mechanism production/loss reactions.

        Works with scalar (ODE/snap) and numpy-array (export) contexts.
        Production = Σ(rate × Π[other reactants]) for reactions producing O1D.
        k_eff      = Σ(rate × Π[other reactants]) per unit [O1D] for loss reactions.
        [O1D]      = production / k_eff
        """
        prod  = 0.0
        k_eff = 0.0
        for rxn in self._o1d_prod_rxns:
            try:
                v = eval(rxn.safe_expr, globals(), ctx)
                for r in rxn.reactants_parsed:   # ODE species: look up by idx
                    sp = self._sp_by_idx.get(r['idx'], '')
                    v = v * np.maximum(ctx.get(sp, 0.0), 0.0) ** r['coeff']
                for r in rxn.reactants_env:       # ENV/PSS species: look up by name
                    v = v * np.maximum(ctx.get(r['name'], 0.0), 0.0) ** r['coeff']
                prod = prod + np.maximum(v, 0.0)
            except Exception:
                pass
        for rxn in self._o1d_loss_rxns:
            try:
                v = eval(rxn.safe_expr, globals(), ctx)
                for r in rxn.reactants_parsed:
                    sp = self._sp_by_idx.get(r['idx'], '')
                    v = v * np.maximum(ctx.get(sp, 0.0), 0.0) ** r['coeff']
                for r in rxn.reactants_env:
                    if r['name'].upper() != 'O1D':   # exclude O1D itself (k_eff side)
                        v = v * np.maximum(ctx.get(r['name'], 0.0), 0.0) ** r['coeff']
                k_eff = k_eff + np.maximum(v, 0.0)
            except Exception:
                pass
        return prod / np.maximum(k_eff, 1e-30)

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
            if 'O1D' in getattr(self, '_pss_species', set()):
                print(">>> Note: O1D was not integrated by the ODE solver — its concentration "
                      "was computed analytically at each step as a pseudo-steady-state "
                      "(PSS: [O1D] = J1·[O3] / k_quench), matching FACSIMILE behaviour.")
            for _note in getattr(self, '_autolink_notes', []):
                print(f">>> Note: {_note}")
            self._show_results(out_path)
        except Exception as e:
            print(f"\n[ERROR] [CRITICAL ERROR] Simulation failed!")
            print(f"Error Details: {str(e)}")
            import traceback
            traceback.print_exc()

    def _setup_environment(self):
        # Prefer an absolute enz_path stored in the config (written by FR_GUI),
        # then fall back to scanning the scenario folder for any *.enz file.
        cfg_enz = self.cfg.get('enz_path', '')
        if cfg_enz and not os.path.isabs(cfg_enz):
            cfg_enz = os.path.join(self.root, cfg_enz)
        if cfg_enz and os.path.isfile(cfg_enz):
            enz_path = cfg_enz
        else:
            enz_candidates = [f for f in os.listdir(self.root) if f.lower().endswith('.enz')]
            enz_file = next((f for f in enz_candidates if 'kammer' in f.lower()), None)
            if enz_file is None:
                enz_file = enz_candidates[0] if enz_candidates else None
            if enz_file is None:
                raise FileNotFoundError(f"No .enz file found in {self.root}")
            enz_path = os.path.join(self.root, enz_file)
        interp_mode = self.cfg.get('interpolation', 'linear')
        step_patterns = _load_step_col_patterns()
        providers, t_min, t_max, headers, enz_times, step_change_times = \
            load_saphir_environment(enz_path, interpolation=interp_mode,
                                    step_col_patterns=step_patterns)
        self._enz_times        = enz_times
        self._step_change_times = step_change_times
        print(f"  [ENZ] Interpolation mode: {interp_mode}  |  "
              f"step patterns: {step_patterns}  |  "
              f"{len(step_change_times)} restart point(s)")

        # Build output time grid: OUTSTEP regular grid + every ENZ row, merged & deduped
        outstep = self.cfg.get('outstep', None)
        if outstep and outstep > 0:
            grid = np.arange(t_min, t_max + outstep * 0.5, outstep)
            t_span = np.union1d(grid, enz_times)
            t_span = t_span[(t_span >= t_min) & (t_span <= t_max)]
            print(f"  [OUTSTEP] {outstep}s grid + {len(enz_times)} ENZ rows = {len(t_span)} output points")
        else:
            t_span = np.union1d(np.linspace(t_min, t_max, num=2000), enz_times)
        self.t_span = t_span

        # Load all ENZ columns under their raw (uppercased) names.
        # In ELOG mode this is sufficient: the EASY script uses INPUT(col) to
        # read them and CONST(formula) to derive everything else.
        for raw_header, interp_func in providers.items():
            self.ctx_data[raw_header.upper()] = interp_func
        # Remember raw ENZ column names for INPUT-alias auto-detection later.
        self._enz_providers_up = {h.upper() for h in providers}

        _using_elog = bool(self.cfg.get('elog_path', ''))

        if _using_elog:
            # ELOG mode: only Fixed Constant env_map entries go into ctx_data.
            # SAPHIR File entries are already present as raw ENZ column names above.
            # No renaming, no overrides, no alias syncing — the EASY script handles all of that.
            for internal_key, mapping in self.cfg['env_map'].items():
                if mapping['type'] != "SAPHIR File":
                    val = float(mapping['value'])
                    self.ctx_data[internal_key] = (
                        lambda t, v=val: np.full_like(t, v) if isinstance(t, np.ndarray) else v)
        else:
            # Non-ELOG (legacy mechanism_files) mode:
            # Process env_map fully: rename ENZ columns to internal names, force
            # DILUTE=0 to disable legacy .dil rate-constant abuse, sync T/TEMP aliases.
            self.ctx_data['DILUTE'] = lambda t: 0.0

            for internal_key, mapping in self.cfg['env_map'].items():
                if mapping['type'] == "SAPHIR File":
                    file_header = mapping['value'].lower()
                    if file_header not in providers:
                        actual_header = next((h for h in headers if h.lower() == file_header), None)
                        if actual_header:
                            file_header = actual_header
                    if file_header in providers:
                        self.ctx_data[internal_key] = providers[file_header]
                else:
                    val = float(mapping['value'])
                    self.ctx_data[internal_key] = (
                        lambda t, v=val: np.full_like(t, v) if isinstance(t, np.ndarray) else v)

            # Sync T↔TEMP and P↔PRESS so a Fixed Constant set under one name
            # also wins over any ENZ-sourced value stored under the other name.
            for const_key, enz_alias in (('TEMP', 'T'), ('T', 'TEMP'),
                                         ('PRESS', 'P'), ('P', 'PRESS')):
                cfg_entry = self.cfg.get('env_map', {}).get(const_key, {})
                if cfg_entry.get('type') not in ('SAPHIR File', '') and const_key in self.ctx_data:
                    self.ctx_data[enz_alias] = self.ctx_data[const_key]

    def get_M(self, t):
        temp_fn = self.ctx_data.get('TEMP') or self.ctx_data.get('T')
        press_fn = self.ctx_data.get('PRESS') or self.ctx_data.get('P')
        # Fallback: read Fixed Constant directly from env_map (handles JSONs where
        # idl_parser did not yet write T/PRESS into env_map)
        if temp_fn is None:
            for _k in ('TEMP', 'T'):
                _m = self.cfg.get('env_map', {}).get(_k)
                if _m and _m.get('type') != 'SAPHIR File':
                    _v = float(_m['value'])
                    temp_fn = lambda t, v=_v: v
                    self.ctx_data[_k] = temp_fn
                    break
        if press_fn is None:
            for _k in ('PRESS', 'P'):
                _m = self.cfg.get('env_map', {}).get(_k)
                if _m and _m.get('type') != 'SAPHIR File':
                    _v = float(_m['value'])
                    press_fn = lambda t, v=_v: v
                    self.ctx_data[_k] = press_fn
                    break
        if temp_fn is None:
            print("[WARNING] No temperature source found — defaulting to 295 K. "
                  "Add T=CONST(295) to the EASY file or set TEMP in env_map.")
            temp_fn = lambda t: 295.0
            self.ctx_data['TEMP'] = temp_fn
        if press_fn is None:
            print("[WARNING] No pressure source found — defaulting to 1013.25 hPa. "
                  "Add P=CONST(1013) to the EASY file or set PRESS in env_map.")
            press_fn = lambda t: 1013.25
            self.ctx_data['PRESS'] = press_fn
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
        # If a JUICE elog exists (generated by FR_GUI preprocessor), use it as the
        # single authoritative flat input — no separate mechanism files needed.
        # This matches FACSIMILE's preprocessor: everything inline, last wins.
        elog_path = self.cfg.get('elog_path', '')
        if elog_path and not os.path.isabs(elog_path):
            elog_path = os.path.join(self.root, elog_path)
        self._va_src = getattr(self, '_va_src', {})     # NAME → source text of its compiled expression
        if elog_path and os.path.isfile(elog_path):
            print(f"  [MECHANISM] Using JUICE elog: {os.path.basename(elog_path)}")
            parser = EasyParser(base_path=os.path.dirname(elog_path))
            self.mechanism = Mechanism()
            m_part = parser.parse(elog_path)
            self.mechanism.species.extend(m_part.species)
            self.mechanism.reactions.extend(m_part.reactions)
            self.mechanism.functions.update(m_part.functions)
            # Phase 1/2/3 structured entries — copy from parsed result
            self.mechanism.function_entries.extend(m_part.function_entries)
            self.mechanism.k_entries.extend(m_part.k_entries)
            self.mechanism.s_entries.extend(m_part.s_entries)
            for p_name, p_expr in m_part.parameters.items():
                p_clean = p_expr.replace('@', '**')
                p_clean = re.sub(r'(?<=\d)[dD](?=[+-]?\d)', 'E', p_clean)
                p_clean = p_clean.strip()
                try:
                    self._va_src[p_name.upper()] = p_clean
                    self.var_assignments.append(
                        (p_name.upper(), compile(p_clean, '<string>', 'eval')))
                except Exception as e:
                    print(f"[WARNING] Could not compile parameter {p_name}: {e}")
            # Store INITIAL[] directives from elog — used instead of cfg['initials']
            self._elog_initial_conditions = m_part.initial_conditions

            # CONST override of ENZ-sourced variables: if the elog ends with e.g.
            # HV=CONST(0) or DIL=CONST(1) after HV=INPUT(hv), "last wins" means
            # the CONST value should override the ENZ input.  Detect this by checking
            # which physical env variables appear in m_part.parameters (now that
            # easy_parser no longer excludes them in elog mode), and update env_map
            # to "Fixed Constant" so ctx_data no longer overrides them.
            _phys_overridden = {}
            for p_name, p_expr in m_part.parameters.items():
                p_up = p_name.upper()
                if p_up in self.cfg.get('env_map', {}):
                    # Evaluate the expression to get the constant value
                    p_clean = p_expr.replace('@', '**')
                    p_clean = re.sub(r'(?<=\d)[dD](?=[+-]?\d)', 'E', p_clean).strip()
                    try:
                        val = float(eval(compile(p_clean, '<s>', 'eval'), {}, {}))
                        _phys_overridden[p_up] = val
                        print(f"  [OVERRIDE] {p_up} = CONST({val}) overrides ENZ input")
                    except Exception:
                        pass  # expression depends on other vars — keep ENZ source
            if _phys_overridden:
                cfg_env = self.cfg.get('env_map', {})
                for name, val in _phys_overridden.items():
                    cfg_env[name] = {'type': 'Fixed Constant', 'value': str(val)}
                    # Also remove from ctx_data so ENZ interpolator is no longer used
                    for k in list(self.ctx_data.keys()):
                        if k.upper() == name:
                            import scipy.interpolate
                            self.ctx_data[k] = lambda t, _v=float(val): _v

            self._finalize_mechanism(parser)
            return

        # Legacy path: parse mechanism files individually (no elog available)
        cfg_chem = self.cfg.get('chemistry_folder', '')
        chem_path = cfg_chem if (cfg_chem and os.path.isdir(cfg_chem)) \
                    else os.path.join(self.root, "chemistry")
        toml_include = ''
        try:
            import tomllib as _tl
            _tp = os.path.join(os.path.dirname(os.path.abspath(__file__)), "FR_GUI.toml")
            with open(_tp, 'rb') as _fh:
                toml_include = _tl.load(_fh).get('paths', {}).get('include_dir', '')
        except Exception:
            pass
        parser = EasyParser(
            base_path=chem_path,
            include_dir=toml_include if (toml_include and os.path.isdir(toml_include)) else None)
        self.mechanism = Mechanism()
        self._mechanism_s_names = set()   # s_entry names from mechanism_files (for cs export)

        mech_files = list(self.cfg['mechanism_files'])

        for f_name in mech_files:
            # Each file may reside in chem_path or include_dir — EasyParser._find_file handles both
            f_path = f_name if os.path.isabs(f_name) else os.path.join(chem_path, f_name)
            if not os.path.exists(f_path):
                # Try absolute path stored by idl_parser (resolved at parse time)
                if os.path.isabs(f_name) and os.path.exists(f_name):
                    f_path = f_name
                else:
                    print(f"[WARNING] Mechanism file not found (skipped): {f_path}")
                    continue
            m_part = parser.parse(f_path)
            self.mechanism.species.extend(m_part.species)
            self.mechanism.reactions.extend(m_part.reactions)
            for _e in m_part.s_entries:
                if not getattr(_e, 'is_synthetic', False):
                    self._mechanism_s_names.add(_e.name.upper())

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
                    
        self._finalize_mechanism(parser)

    def _finalize_mechanism(self, parser):
        """Common post-processing after all mechanism files (or elog) have been parsed."""
        # Post-process all accumulated reactions using the fully-populated rates_map.
        # When files are parsed in sequence, a later SAPHIR override file may update
        # rates_map entries for reactions already materialised in an earlier .mech file.
        # Re-applying the final rates_map here ensures those overrides take effect.
        def _rxn_canon(r):
            reac_str = '+'.join(str(x) for x in r.reactants)
            prod_str  = '+'.join(str(x) for x in r.products)
            return EasyParser._canon_rxn_key(reac_str, prod_str)

        for rxn in self.mechanism.reactions:
            # rates_map is keyed exclusively on canonical form — single lookup.
            canon_key = _rxn_canon(rxn)
            final_rate = parser.rates_map.get(canon_key)
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

        # IDL-mode computed parameters (FLOW, DILUTE, DIL, SZ, H2O, …) — evaluated
        # each ODE step. Prepended so chemistry-file parameters (e.g. COSX =
        # COS(SZ * 0.0174533) in j_clear_sky.incl) see the master-easy scenario
        # values. If SZ/COSX are evaluated AFTER the J-values that reference them,
        # the Js see stale zero values (Bug 5, April 2026).
        # Precondition (handled in idl_parser): source_consts that duplicate
        # env_map keys are NOT promoted here — their formulas may reference ENZ
        # column aliases that aren't in ctx and would silently zero the value.
        _using_elog = bool(self.cfg.get('elog_path', ''))

        # Prepend idl_parameters to var_assignments.
        # In elog mode these come from the EASY script (parsed via elog by idl_parser)
        # and set physical variables (H2O, COSX, T, P, etc.) early in var_assignments
        # so that dependent parameters (KMT06 uses H2O, J-values use COSX) see
        # correct values in the first evaluation pass.
        # idl_parameters come from the JSON (written by idl_parser) and may contain
        # unexpanded FUNCTION calls like COMPRC(...).  Apply expand_macros here using
        # the parser's function table so the JSON path gets the same inline-expansion
        # as the elog path.
        from easy_parser import expand_macros as _expand_macros
        _idl_params = list(self.cfg.get('idl_parameters', {}).items())
        # idl_parameters are computed values (e.g. H2O = 1e-2*M*H2O_IN) that must
        # override raw ENZ-sourced values even when the name is in env_map/active_s.
        _idl_params_upper = {k.upper() for k, _ in _idl_params}
        self._idl_params_upper = _idl_params_upper   # expose for _solve / _export_results
        for p_name, p_expr in reversed(_idl_params):
            p_clean = p_expr.upper().replace('@', '**')
            p_clean = re.sub(r'(?<=\d)[dD](?=[+-]?\d)', 'E', p_clean)
            p_clean = _expand_macros(p_clean, parser.functions)   # expand COMPRC, PHTFKT, …
            p_clean = p_clean.replace('@', '**')   # re-apply: macro bodies may reintroduce @
            p_clean = re.sub(r'(?<=\d)[dD](?=[+-]?\d)', 'E', p_clean)   # re-apply D→E
            p_clean = p_clean.strip()
            try:
                self._va_src[p_name.upper()] = p_clean
                self.var_assignments.insert(0, (p_name.upper(), compile(p_clean, '<string>', 'eval')))
            except Exception as e:
                print(f"[WARNING] Could not compile IDL parameter {p_name}: {e}")

        all_species = set(self.mechanism.species)
        for b in self.cfg.get('burst_injections', []):
            if b.get('target'): all_species.add(b['target'].upper())
        for s in self.cfg.get('continuous_sources', []):
            if s.get('target'): all_species.add(s['target'].upper())

        if _using_elog and hasattr(self, '_elog_initial_conditions'):
            # Species named in INITIAL[X] are ODE species even if not in the mechanism
            # (e.g. H2O defined via INITIAL[H2O] instead of as a parameter formula).
            for spec in self._elog_initial_conditions.keys():
                if spec != '*':
                    all_species.add(spec.upper())
        elif not _using_elog:
            for spec in self.cfg.get('initials', {}).keys():
                all_species.add(spec.upper())

            for r in self.cfg.get('custom_reactions', []):
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

            mapped_up_set = {k.upper() for k in self.cfg.get('env_map', {})}
            for r in self.cfg.get('chamber_source_reactions', []):
                if not r.get('products') or not r.get('rate'): continue
                reacts_str = r.get('reactants', '')
                reactants = [x.strip().upper() for x in reacts_str.split('+') if x.strip()] if reacts_str.strip() else []
                products = [x.strip().upper() for x in r['products'].split('+') if x.strip()]
                for s in reactants:
                    clean_s = re.sub(r'^[\d\.]+\s*\*\s*', '', s.strip())
                    if clean_s and not clean_s[0].isdigit() and clean_s not in mapped_up_set and clean_s not in self.ctx_data:
                        all_species.add(clean_s)
                for s in products:
                    clean_s = re.sub(r'^[\d\.]+\s*\*\s*', '', s.strip())
                    if clean_s and not clean_s[0].isdigit(): all_species.add(clean_s)
                self.mechanism.reactions.append(Reaction(reactants, products, r['rate']))

            for r in self.cfg.get('injection_reactions', []):
                if not r.get('products') or not r.get('rate'): continue
                reacts_str = r.get('reactants', '')
                reactants = [x.strip().upper() for x in reacts_str.split('+') if x.strip()] if reacts_str.strip() else []
                products = [x.strip().upper() for x in r['products'].split('+') if x.strip()]
                for s in reactants:
                    clean_s = re.sub(r'^[\d\.]+\s*\*\s*', '', s.strip())
                    if clean_s and not clean_s[0].isdigit() and clean_s not in self.ctx_data:
                        all_species.add(clean_s)
                for s in products:
                    clean_s = re.sub(r'^[\d\.]+\s*\*\s*', '', s.strip())
                    if clean_s and not clean_s[0].isdigit(): all_species.add(clean_s)
                self.mechanism.reactions.append(Reaction(reactants, products, r['rate']))

            for r in self.cfg.get('idl_dilution_reactions', []):
                if not r.get('reactants') or not r.get('rate'): continue
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
                    if clean_s and not clean_s[0].isdigit(): all_species.add(clean_s)
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

        # Determine which "auto" behaviours are needed because the EASY script
        # did not define the quantity explicitly.
        _idl_params  = set(self.cfg.get('idl_parameters', {}).keys())
        _env_map     = self.cfg.get('env_map', {})
        _j_scalings  = set(self.cfg.get('j_scalings', {}).keys())

        # H2O must be defined in the EASY script (CONST, INPUT, or formula).
        # No automatic unit conversion — if ENZ stores H2O in %, the EASY script
        # must convert explicitly (e.g. H2O = CONST(INPUT(h2o)/100 * M)).
        self._h2o_from_env_map = False   # no longer used for conversion
        self._h2o_from_easy    = 'H2O' in _env_map or 'H2O' in _idl_params

        self._autolink_notes = []
        if not self._h2o_from_easy:
            self._autolink_notes.append(
                "H2O auto-balancer was active: H2O is not defined in the EASY script — "
                "it was computed from RH (ENZ) at every ODE step. "
                "Add 'H2O = ...' to the EASY script to make the run fully self-contained.")
        if 'O1D' in all_species:
            all_species.discard('O1D')
            self._pss_species.add('O1D')
            # Dynamic PSS: pre-identify which reactions produce and consume O1D.
            # This works regardless of how the photolysis reaction is written
            # (O3+HV→O1D, O3→O1D, etc.) and uses the actual EASY-script rate.
            self._o1d_prod_rxns = []   # reactions that have O1D as a product
            self._o1d_loss_rxns = []   # reactions that consume O1D
            print("  [PSS] O1D removed from ODE — PSS computed dynamically from mechanism reactions")


        # --- O atom: kept in ODE (not PSS) ---
        # FACSIMILE also integrates O in the ODE.  Both models share the same
        # numerical treatment (O << atol=1e4, so both are inaccurate by the same
        # amount).  Switching JUICE to PSS while FACSIMILE stays in ODE would
        # cause systematic O3 disagreement (~20%) because FACSIMILE's elevated
        # numerical O drives faster O→O3 conversion.  The 16× O discrepancy in
        # the comparison reflects a shared numerical artefact, not a chemistry bug.

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

        # ── VALIDATION: no unresolved function names in var_assignments ──────────────
        # If expand_macros in easy_parser failed to inline-expand a FUNCTION call
        # (e.g. COMPRC, KMT01, PHTFKT), its name would appear as a co_name in the
        # compiled expression and evaluate to SafeDict(0.0) at runtime — silently
        # wrong.  Catch this here and STOP so the user can fix the EASY script.
        _va_names    = {n.upper() for n, _ in self.var_assignments}
        _sp_names    = {s.upper() for s in all_species}
        _env_names   = {k.upper() for k in self.cfg.get('env_map', {}).keys()}
        _ctx_names   = {k.upper() for k in getattr(self, 'ctx_data', {}).keys()}
        _math_names  = {'AMAX', 'AMIN', 'EXP', 'LOG', 'LOG10', 'SQRT', 'ABS',
                        'COS', 'SIN', 'TAN', 'NP', 'MAX', 'MIN', 'TRUE', 'FALSE',
                        'NONE', 'INT', 'FLOAT', 'STR', 'ROUND', 'POW',
                        # numpy attributes accessed via np.xxx (show up in co_names)
                        'MAXIMUM', 'MINIMUM', 'WHERE', 'POWER', 'INTERP',
                        'NANMAX', 'NANMIN', 'ISNAN', 'ISINF', 'ZEROS', 'ONES',
                        'FULL', 'ARRAY', 'ASARRAY', 'RAVEL', 'RADIANS',
                        'CONCATENATE', 'LINSPACE', 'ARANGE'}
        # Standard EASY/atmospheric variables that are always available in ctx
        # (T, M via get_M; HV, JNO2 from ENZ; COSX via var_assignments or SNZ).
        # RESERVED in easy_parser so they don't appear in m_part.parameters, but
        # they ARE correctly set in ctx at runtime — not a validation problem.
        _phys_names  = {'M', 'T', 'TEMP', 'P', 'PRESS', 'N2', 'O2', 'H2O',
                        'HV', 'RH', 'DIL', 'DILUTE', 'FLOW', 'SZ', 'SZA',
                        'COSX', 'JNO2', 'LNO2', 'MNO2', 'NNO2', 'JSCALE',
                        'JO3', 'JNO3', 'JNO2_1', 'DENSITY'}
        _known = _va_names | _sp_names | _env_names | _ctx_names | _math_names | _phys_names
        _unresolved = []
        _va_src = getattr(self, '_va_src', {})
        _added_species = {}
        for va_name, code in self.var_assignments:
            for co_name in code.co_names:
                co_up = co_name.upper()
                if co_up in _known or co_up.startswith('_') or co_up in _added_species:
                    continue
                src = _va_src.get(va_name.upper(), '')
                if src and not re.search(rf'\b{re.escape(co_name)}\s*\(', src, re.I):
                    # A plain name, not a call: FACSIMILE takes any such name in an expression as
                    # a variable (a species with no reactions, holding INITIAL[*]).  Do the same,
                    # so a sum such as RO2 that names a species the mechanism does not have (a
                    # truncated include: 'COHM2CO' for COHM2CO3) runs as it did in FACSIMILE.
                    _added_species[co_up] = va_name
                    continue
                _unresolved.append((va_name, co_name))
        if _added_species:
            for co_up, va_name in _added_species.items():
                print(f"  [WARNING] '{co_up}' in {va_name} is in no reaction: kept as a species without "
                      f"chemistry, as FACSIMILE does (a misspelt species in the include file?)")
            self.mechanism.species = sorted(set(self.mechanism.species) | set(_added_species))
            all_species.update(_added_species)
        if _unresolved:
            lines = '\n'.join(f"  {va} uses '{co}'" for va, co in _unresolved[:10])
            raise RuntimeError(
                f"Unresolved name(s) in compiled EASY parameters — "
                f"expand_macros may have missed a FUNCTION definition:\n{lines}\n"
                f"Fix: ensure the FUNCTION is defined in the EASY script before "
                f"it is used, or add it to the elog.")

        # ── Auto-register INPUT alias variables ──────────────────────────────────
        # Any var_assignment whose compiled expression has exactly ONE co_name and
        # that co_name is a raw ENZ column (present in _enz_providers_up) is an
        # INPUT alias — it just reads a raw ENZ column under a different name.
        # Examples:  H2O_IN = H2O,  T_IN = T,  HV_IN = HV,  FLOW = DILUTE
        # These must be frozen in the 3-pass evaluation loop so that downstream
        # formulas (DILUTE = f(FLOW), H2O = 1e-2*M*H2O_IN) always see the raw
        # ENZ value, not a value that was computed in a previous pass.
        _enz_up = getattr(self, '_enz_providers_up', set())
        for vn, ce in self.var_assignments:
            co = set(getattr(ce, 'co_names', []))
            if len(co) == 1:
                src = next(iter(co))
                if (src in _enz_up and vn not in self._ctx_data_protected
                        and not getattr(ce, 'co_consts', ())):
                    # Pure INPUT alias: vn = src (no constants, no operations).
                    # Formulas like OH = 3e6*HV have co_consts=(3e6,) → not an alias.
                    if vn not in self.ctx_data:
                        self.ctx_data[vn] = self.ctx_data[src]
                    self._ctx_data_protected.add(vn)

        # ── Phase 3: build k/s/c/enz/function vectors ────────────────────────
        # These run alongside the existing var_assignments for now.
        # Phase 3 wires them into the ODE loops; Phase 4 removes the old code.
        if self.cfg.get('elog_path'):
            from mechanism import SEntry, EnzEntry
            from evaluator import update_enz_vec as _uenz, update_s_vec as _us
            # enz_vec: one entry per raw ENZ column (original case)
            self._enz_vec = [EnzEntry(name=col) for col in self._enz_providers_up
                             if col.lower() in {h.lower() for h in
                                                getattr(self, '_enz_providers_up', set())}]
            # rebuild properly from ctx_data keys that are ENZ providers
            _prov_names = {k for k in self.ctx_data
                           if k.upper() in self._enz_providers_up}
            self._enz_vec = [EnzEntry(name=n) for n in _prov_names]
            # Store the providers dict for update_enz_vec
            self._enz_providers = {n: self.ctx_data[n] for n in _prov_names}

            # s_vec: all entries from mech.s_entries (populated by easy_parser)
            self._s_vec = [SEntry(name=e.name, definition=e.definition,
                                  is_synthetic=getattr(e, 'is_synthetic', False))
                           for e in self.mechanism.s_entries]

            # function_entries from mechanism
            self._function_entries = list(self.mechanism.function_entries)

            # Pre-compute evaluation order once (avoid topological_sort at every ODE call)
            from evaluator import topological_sort as _tsort, precompile_s_vec as _pcs
            self._s_vec_ordered = _tsort(self._s_vec, self._function_entries)
            _pcs(self._s_vec_ordered, self._function_entries)  # pre-compile code objects

            # Build a single precompiled batch function for s_vec evaluation.
            # Replaces N eval() calls in update_s_vec with one Python function call.
            # Same principle as _compute_rates for reaction rates.
            _svec_fn_lines = ["def _eval_s_vec_batch(_ctx):"]
            _svec_fn_lines.append("    _c = _ctx  # shorthand")
            for _e in self._s_vec_ordered:
                _n = _e.name
                if getattr(_e, '_input_col', None) is not None:
                    _col = _e._input_col
                    _svec_fn_lines.append(
                        f"    {_n} = _c.get('{_col}', _c.get('{_col.lower()}', 0.0)); _c['{_n}'] = {_n}")
                elif getattr(_e, '_expr_str', None) is not None:
                    _svec_fn_lines.append(
                        f"    try:\n"
                        f"        {_n} = float({_e._expr_str}); _c['{_n}'] = {_n}\n"
                        f"    except: pass")
            if len(_svec_fn_lines) == 2:
                _svec_fn_lines.append("    pass")
            _svec_src = "\n".join(_svec_fn_lines)
            _svec_ns = {
                '__builtins__': {'abs': abs, 'max': max, 'min': min,
                                 'float': float, 'int': int},
                'exp': np.exp, 'EXP': np.exp, 'log': np.log, 'LOG': np.log,
                'log10': np.log10, 'LOG10': np.log10,
                'cos': np.cos, 'sin': np.sin, 'sqrt': np.sqrt, 'SQRT': np.sqrt,
                'np': np,
            }
            try:
                exec(_svec_src, _svec_ns)
                self._eval_s_vec_batch = _svec_ns['_eval_s_vec_batch']
            except Exception as _exc:
                self._eval_s_vec_batch = None

            # ── @njit version of the s_vec batch (array-based, no dict) ──────
            # Build index map: name → position in a float64 vals array.
            # vals[0..N_enz-1] = ENZ column values (filled by caller before call).
            # vals[N_enz..]    = s_vec computed values (written by @njit).
            self._eval_s_vec_njit = None
            self._svec_enz_names  = []   # ENZ column names in vals order
            self._svec_s_names    = []   # s_vec output names in vals order
            if _NUMBA_AVAILABLE and self._eval_s_vec_batch is not None:
                import re as _re2
                _MATH_WORDS = {'EXP','LOG','LOG10','SQRT','SIN','COS','TAN','ABS',
                               'AMAX','AMIN','MAX','MIN','NP','FLOAT','INT'}
                _enz_names_ord = [_e.name.upper() for _e in self._enz_vec]
                _svec_names_ord = []
                _seen_n = set(_enz_names_ord)
                for _e in self._s_vec_ordered:
                    nm = _e.name.upper()
                    if nm not in _seen_n:
                        _svec_names_ord.append(nm)
                        _seen_n.add(nm)
                _all_names = _enz_names_ord + _svec_names_ord
                _n2i = {nm: i for i, nm in enumerate(_all_names)}
                _N_enz = len(_enz_names_ord)
                self._svec_enz_names = _enz_names_ord
                self._svec_s_names   = _svec_names_ord
                self._svec_n2i       = _n2i
                self._svec_N_enz     = _N_enz

                def _replace_vars(expr_s):
                    out = expr_s
                    for _nm, _idx in sorted(_n2i.items(), key=lambda x: -len(x[0])):
                        if _nm in _MATH_WORDS: continue
                        out = _re2.sub(rf'\b{_re2.escape(_nm)}\b',
                                       f'_v[{_idx}]', out)
                    return out

                _nj_lines = ["def _s_vec_njit_fn(_v):"]
                _nj_lines.append("    import math as _m")
                for _e in self._s_vec_ordered:
                    _nm  = _e.name.upper()
                    _idx = _n2i.get(_nm, -1)
                    if _idx < 0: continue
                    _ic  = getattr(_e, '_input_col', None)
                    _es  = getattr(_e, '_expr_str', None)
                    if _ic is not None:
                        _col_idx = _n2i.get(_ic.upper(), -1)
                        if _col_idx >= 0:
                            _nj_lines.append(f"    _v[{_idx}] = _v[{_col_idx}]")
                    elif _es is not None:
                        _expr_r = _replace_vars(_es)
                        _expr_r = _re2.sub(r'\bEXP\b',   '_m.exp',   _expr_r)
                        _expr_r = _re2.sub(r'\bLOG10\b', '_m.log10', _expr_r)
                        _expr_r = _re2.sub(r'\bLOG\b',   '_m.log',   _expr_r)
                        _expr_r = _re2.sub(r'\bSQRT\b',  '_m.sqrt',  _expr_r)
                        _expr_r = _re2.sub(r'\bCOS\b',   '_m.cos',   _expr_r)
                        _expr_r = _re2.sub(r'\bSIN\b',   '_m.sin',   _expr_r)
                        _expr_r = _re2.sub(r'\bAMAX\b\s*\(', 'max(', _expr_r)
                        _expr_r = _re2.sub(r'\bAMIN\b\s*\(', 'min(', _expr_r)
                        _nj_lines.append(
                            f"    try: _v[{_idx}] = float({_expr_r})\n    except: pass")
                if len(_nj_lines) == 2: _nj_lines.append("    pass")
                _nj_src = "\n".join(_nj_lines)
                _nj_ns = {'__builtins__': {'abs': abs, 'max': max, 'min': min,
                                           'float': float, 'int': int}}
                try:
                    exec(_nj_src, _nj_ns)
                    _py_fn = _nj_ns['_s_vec_njit_fn']
                    self._eval_s_vec_njit = _numba.njit(cache=True)(_py_fn)
                    # Warm up JIT with a dummy call
                    _dummy = np.zeros(len(_all_names), dtype=np.float64)
                    self._eval_s_vec_njit(_dummy)
                    _nj_label = "njit+batch"
                except Exception as _exc2:
                    self._eval_s_vec_njit = None
                    _nj_label = f"batch only (njit failed: {_exc2})"
            else:
                _nj_label = "batch only (numba unavailable)"

            _batch_ok = "OK" if self._eval_s_vec_batch else "FAILED"
            print(f"  [Phase3] enz_vec: {len(self._enz_vec)}  "
                  f"s_vec: {len(self._s_vec)}  "
                  f"function_entries: {len(self._function_entries)}  "
                  f"(batch={_batch_ok}, {_nj_label})")
        else:
            self._enz_vec = []
            self._s_vec = []
            self._function_entries = []

        self._add_artifacts()

    def _add_artifacts(self):
        # When loading from an EASY/ELOG file, all chamber chemistry (Y+OH,
        # HONO wall source, HCHO source, …) must be defined in the EASY script
        # itself.  Adding extra reactions here would silently override or
        # duplicate what the EASY script specifies.  If a reaction is missing,
        # fix the EASY script — not JUICE.
        if self.cfg.get('elog_path'):
            return   # EASY/ELOG mode: no artifact backdoors

        if not self.cfg.get('artifacts') or not self.cfg['artifacts'].get('active', False):
            return
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
                self.mechanism.reactions.append(Reaction([], ["HONO"], str(self.cfg['artifacts']['dark_hono'])))

    def _solve(self):
        self.solver = MechanismSolver(self.mechanism, self.cfg['env_map'])
        y0 = np.zeros(self.solver.num_species)
        # Species held to a measurement (idl_parser 'held_species': X = CONST(<expression of an
        # INPUT>)).  Only ODE species count; the value is set from the settings snapshot at the
        # start of every ENZ row and kept constant within the row, as FACSIMILE does for a
        # species that is also a parameter.  {species index: (name, compiled expression)}
        self._held = {}
        _upper_idx = {s.upper(): i for s, i in self.solver.species_to_idx.items()}
        for _hn, _hx in (self.cfg.get('held_species') or {}).items():
            if _hn.upper() in _upper_idx:
                self._held[_upper_idx[_hn.upper()]] = (self.solver.species[_upper_idx[_hn.upper()]],
                                                       compile(str(_hx), '<held>', 'eval'))
        if self._held:
            print(f"  [HELD] {len(self._held)} species held to the ENZ: "
                  + ", ".join(n for n, _c in self._held.values()))

        saphir_base_date = pd.Timestamp("2000-01-01 00:00:00")

        M0 = self.get_M(self.t_span[0])
        _temp_fn = self.ctx_data.get('TEMP') or self.ctx_data.get('T')
        T0 = float(_temp_fn(self.t_span[0])) if _temp_fn else 298.15

        _using_elog = bool(self.cfg.get('elog_path', ''))
        if _using_elog and hasattr(self, '_elog_initial_conditions'):
            # Evaluate INITIAL[spec] = CONST(expr) directly from elog parse.
            # Expressions are already in molec/cm³ (e.g. M*1e-9*10 = 10 ppb).
            _elog_ctx = {'M': M0, 'T': T0, 'EXP': np.exp, 'LOG': np.log,
                         'LOG10': np.log10, 'SQRT': np.sqrt}
            elog_ics = self._elog_initial_conditions

            floor_expr = elog_ics.get('*', '0.')
            try:
                floor_val = float(eval(floor_expr.replace('@', '**'), {}, _elog_ctx))
            except Exception:
                floor_val = 0.0
            if floor_val > 0:
                for idx in range(self.solver.num_species):
                    y0[idx] = floor_val
                print(f"  [INIT] INITIAL[*] = {floor_val:.2e} molec/cm³ applied to all {self.solver.num_species} species")

            for spec, expr in elog_ics.items():
                if spec == '*': continue
                spec_up = spec.upper()
                if spec_up in self.solver.species_to_idx:
                    try:
                        val = float(eval(expr.replace('@', '**'), {}, _elog_ctx))
                        y0[self.solver.species_to_idx[spec_up]] = val
                        print(f"  [INIT] {spec_up} = {val:.3e} molec/cm³  ({val/M0*1e9:.4g} ppb)")
                    except Exception as e:
                        print(f"[WARNING] Could not evaluate INITIAL[{spec}] = {expr}: {e}")
        else:
            initials = self.cfg.get('initials', {})
            floor_ppb = float(initials.get('*', 0.0))
            if floor_ppb > 0:
                for idx in range(self.solver.num_species):
                    y0[idx] = floor_ppb * 1e-9 * M0
                print(f"  [INIT] INITIAL[*] floor applied: {floor_ppb:.2e} ppb to all {self.solver.num_species} species")
            for spec, ppb in initials.items():
                if spec == '*' or not ppb: continue
                spec_up = spec.upper()
                if spec_up in self.solver.species_to_idx:
                    y0[self.solver.species_to_idx[spec_up]] = (float(ppb) * 1e-9) * M0

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
            
            # KMT/COMPRC are now inline-expanded parameters — no mcm_kinetics injection

            _test_cosx = max(np.cos(45.0 * 0.0174533), 0.01)
            test_ctx['COSX']     = _test_cosx
            test_ctx['LNO2']     = 1.165e-2
            test_ctx['MNO2']     = 0.244
            test_ctx['NNO2']     = -0.267
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
                for vn, err in failed_vars[:10]:
                    print(f"  {vn}: {err[:80]}")

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
                    else:
                        # All non-ODE reactants (DIL, HV, Q, env vars, …) go to
                        # reactants_env so their ctx value multiplies the rate.
                        # ctx.get(name, 0.0) ensures missing names contribute 0.
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
        # All rate functions (KMT01–16, COMPRC, PHTFKT, …) are inline-expanded at parse
        # time by easy_parser.expand_macros.  No separate KMT callable injection needed.

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
        if len(_fn_lines) == 1:   # no reactions — guard against empty function body
            _fn_lines.append("    pass")
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

        # Store compiled rate function for use in _export_kinetics
        self._compute_rates_fn = _compute_rates

        # ── VECTORISED STOICHIOMETRY ENGINE ─────────────────────────────────────────
        # Pre-build flat index arrays so the per-reaction Python loop in system_ode
        # becomes a handful of NumPy operations. For 17,000 reactions this cuts each
        # RHS call from ~5ms to ~0.5ms.
        #
        # For each reaction we need:
        #   v[i] = k[i] * prod(y[r_idx]^r_coeff) * prod(ctx[r_name]^r_coeff)
        #   dydt[r_idx] -= v[i] * r_coeff   (for each ODE reactant)
        #   dydt[p_idx] += v[i] * p_coeff   (for each ODE product)
        #
        # We pre-compute:
        #   _ode_reactant arrays: (rxn_idx, species_idx, coeff) for ODE species in reactants
        #   _env_reactant arrays: (rxn_idx, env_name, coeff) for env variables in reactants
        #   _stoich arrays: (rxn_idx, species_idx, signed_coeff) for all dydt contributions

        _ode_r_rxn, _ode_r_sp, _ode_r_coeff = [], [], []
        _env_r_rxn, _env_r_names, _env_r_coeff = [], [], []
        _stoich_rxn, _stoich_sp, _stoich_sign_coeff = [], [], []

        for i, rxn in enumerate(self.mechanism.reactions):
            for r in rxn.reactants_parsed:
                _ode_r_rxn.append(i)
                _ode_r_sp.append(r['idx'])
                _ode_r_coeff.append(r['coeff'])
            for r in rxn.reactants_env:
                _env_r_rxn.append(i)
                _env_r_names.append(r['name'])
                _env_r_coeff.append(r['coeff'])
            # Stoichiometry: reactants consumed (negative), products produced (positive)
            for r in rxn.reactants_parsed:
                _stoich_rxn.append(i)
                _stoich_sp.append(r['idx'])
                _stoich_sign_coeff.append(-r['coeff'])
            for p in rxn.products_parsed:
                _stoich_rxn.append(i)
                _stoich_sp.append(p['idx'])
                _stoich_sign_coeff.append(p['coeff'])

        _ode_r_rxn   = np.array(_ode_r_rxn, dtype=np.int32)
        _ode_r_sp    = np.array(_ode_r_sp, dtype=np.int32)
        _ode_r_coeff = np.array(_ode_r_coeff, dtype=np.float64)
        _env_r_rxn   = np.array(_env_r_rxn, dtype=np.int32)
        _env_r_coeff = np.array(_env_r_coeff, dtype=np.float64)
        _stoich_rxn  = np.array(_stoich_rxn, dtype=np.int32)
        _stoich_sp   = np.array(_stoich_sp, dtype=np.int32)
        _stoich_sc   = np.array(_stoich_sign_coeff, dtype=np.float64)
        # Pre-group env reactants by unique name for batch lookup
        _env_r_unique_names = list(dict.fromkeys(_env_r_names))  # preserve order, dedup
        _env_r_name_idx = np.array([_env_r_unique_names.index(n) for n in _env_r_names], dtype=np.int32)
        n_rxn = len(self.mechanism.reactions)

        print(f">>> [ODE] Vectorised stoichiometry: {len(_ode_r_rxn)} ODE reactant terms, "
              f"{len(_env_r_rxn)} env reactant terms, {len(_stoich_rxn)} stoichiometry terms")

        # PSS species set — captured in system_ode closure for O1D (and future additions)
        _pss_species_set = getattr(self, '_pss_species', set())

        # Identify O1D production/loss reactions now that reactants_parsed/env are set.
        # Use raw rxn.products (not products_parsed) since PSS species are not in ODE.
        if 'O1D' in _pss_species_set and hasattr(self, '_o1d_prod_rxns'):
            for rxn in self.mechanism.reactions:
                prods_up = {p.split('*')[-1].strip().upper() for p in rxn.products}
                reacts_env_up = {r['name'].upper() for r in rxn.reactants_env}
                if 'O1D' in prods_up:
                    self._o1d_prod_rxns.append(rxn)
                elif 'O1D' in reacts_env_up:
                    self._o1d_loss_rxns.append(rxn)
            print(f"  [PSS] O1D: {len(self._o1d_prod_rxns)} production, "
                  f"{len(self._o1d_loss_rxns)} loss reactions")

        # ── Dynamic PSS helper ─────────────────────────────────────────────────
        # Store index→name map on self so _pss_o1d (bound method) can use it.
        self._sp_by_idx = {idx: sp for sp, idx in self.solver.species_to_idx.items()}
        _pss_o1d = self._pss_o1d   # alias for closures inside _solve

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

        # ── FACSIMILE architecture: compute all settings once per output step ──────
        # s-values (KMT rates, J-values, DILUTE, COSX, …) are computed ONCE at the
        # start of each output step from the current c-vector (ODE state) and the
        # current s-vector (ENZ inputs).  They remain fixed while the ODE integrates
        # that step — exactly as FACSIMILE does.  Only c-species (H2O, OH, …) are
        # updated at every RHS call from the ODE state vector y.
        _s_snap: dict = {}          # settings snapshot, rebuilt before each step
        _snap_active: set = set()   # active_constraints for the current snapshot

        def _make_snap(y_state, t_abs):
            """Build a settings snapshot from the current ODE state y_state and ENZ at t_abs."""
            M_s = self.get_M(t_abs)
            T_fn_s = self.ctx_data.get('TEMP') or self.ctx_data.get('T')
            T_s = float(T_fn_s(t_abs))
            ctx_s = SafeDict()
            active_s = set()
            _env_map_up  = {k.upper() for k in self.cfg['env_map']}
            for k, v_func in self.ctx_data.items():
                try:
                    ctx_s[k.upper()] = float(v_func(t_abs))
                    if k.upper() in _env_map_up or k.upper() in self._ctx_data_protected:
                        active_s.add(k.upper())
                except Exception:
                    pass
            for a, b in (('TEMP','T'),('T','TEMP'),('PRESS','P'),('P','PRESS')):
                if a in active_s and b in ctx_s:
                    active_s.add(b)
            ctx_s.update({"T": T_s, "TEMP": T_s, "M": M_s,
                          "np": np, "exp": np.exp, "EXP": np.exp,
                          "LOG10": np.log10, "log10": np.log10,
                          "LOG": np.log, "log": np.log,
                          "COS": np.cos, "SIN": np.sin, "TAN": np.tan,
                          "SQRT": np.sqrt, "AMAX": np.maximum, "AMIN": np.minimum})
            # Load current ODE state (c-values feed into settings like KMT06=f(H2O))
            for i, sp in enumerate(_species_names):
                ctx_s[sp] = float(y_state[i])
            # Held species: their value for this row, from the ENZ (ctx_s holds every column) and M
            for _hi, (_hname, _hcode) in self._held.items():
                try:
                    _hv = max(float(eval(_hcode, globals(), ctx_s)), 0.0)
                except Exception:
                    continue
                y_state[_hi] = _hv
                ctx_s[_hname] = _hv
            # Re-apply constrained species
            for k in active_s:
                try:
                    ctx_s[k] = float(self.ctx_data[k.lower() if k.lower() in self.ctx_data else k](t_abs))
                except Exception:
                    pass
            # Evaluate settings via new vector system (ELOG) or old 3-pass (non-ELOG)
            if self._s_vec:
                from evaluator import update_enz_vec as _uenz, update_s_vec as _us, build_ctx as _bctx
                from mechanism import CEntry
                _uenz(self._enz_vec, self._enz_providers, t_abs)
                _c_vec = [CEntry(name=sp, value=float(y_state[i]))
                          for i, sp in enumerate(_species_names)]
                _us(self._s_vec, self._function_entries, _c_vec, self._enz_vec,
                    _ordered=self._s_vec_ordered)
                ctx_s.update(_bctx(self._s_vec, _c_vec, self._enz_vec))
            else:
                _snap_env_guard = active_s - self._idl_params_upper if not self.cfg.get('elog_path') else set()
                for _p in range(3):
                    for vn, ce in self.var_assignments:
                        if vn in self._ctx_data_protected:
                            continue
                        if vn in _snap_env_guard:
                            continue
                        try: ctx_s[vn] = eval(ce, globals(), ctx_s)
                        except:
                            if _p == 0: ctx_s[vn] = 0.0
            # O1D PSS — dynamic from mechanism reactions
            if 'O1D' in _pss_species_set:
                ctx_s['O1D'] = float(_pss_o1d(ctx_s))
            _s_snap.clear()
            _s_snap.update(ctx_s)
            _snap_active.clear()
            _snap_active.update(active_s)
            if _k_row_cacheable:
                # ctx_s holds the frozen settings + species at row start; since no
                # rate expression references a species, the result is row-constant.
                try:
                    _compute_rates(ctx_s, _k_row_buf)
                    np.nan_to_num(_k_row_buf, nan=0.0, posinf=0.0, neginf=0.0, copy=False)
                    _k_row_valid[0] = True
                except Exception:
                    _k_row_valid[0] = False
            # Optional observer (f0am_export.py): sees every settings snapshot exactly as the
            # solver uses it.  Never set during normal runs.
            _hook = getattr(self, 'snap_hook', None)
            if _hook is not None:
                _hook(t_abs, ctx_s, y_state, active_s)

        # Set of actual ODE species (c-values) — these update from y during ODE;
        # everything else (settings: T, HV, KMT*, J-values, …) stays frozen from snap.
        _species_set_upper = {s.upper() for s in _species_names}


        # ODE backend — read from scenario JSON (set by GUI dropdown).
        _ode_backend = self.cfg.get('ode_backend', 'python')
        _use_numba   = ('numba' in _ode_backend) and _NUMBA_AVAILABLE
        _numba_missing = ('numba' in _ode_backend) and not _NUMBA_AVAILABLE
        self._ode_backend = _ode_backend if not _numba_missing else \
                            _ode_backend.replace('numba', 'python')
        if _numba_missing:
            print(">>> [ODE] WARNING: numba not installed — kernel falls back to python")
        _backend_label = _ode_backend + (" (python kernel)" if _numba_missing else "")
        print(f">>> [ODE] Backend: {_backend_label}")

        # Pre-allocated buffers — avoids array allocation on every ODE call.
        _env_vals_buf    = np.zeros(len(_env_r_unique_names), dtype=np.float64)
        _v_buf_stoich    = np.zeros(len(self.mechanism.reactions), dtype=np.float64)
        _dydt_stoich_buf = np.zeros(self.solver.num_species, dtype=np.float64)
        # vals buffer for @njit s_vec evaluation
        _svec_njit_fn  = getattr(self, '_eval_s_vec_njit', None)
        _svec_n2i      = getattr(self, '_svec_n2i', {})
        _svec_enz_nms  = getattr(self, '_svec_enz_names', [])
        _svec_s_nms    = getattr(self, '_svec_s_names', [])
        _svec_N_enz    = getattr(self, '_svec_N_enz', 0)
        _n_svals       = len(_svec_enz_nms) + len(_svec_s_nms)
        _svals_buf     = np.zeros(_n_svals, dtype=np.float64) if _n_svals > 0 else None

        # Profiling accumulators — time (s) spent in each section of system_ode.
        _prof = {'snap_read': 0.0, 'ctx_data': 0.0, 'elog': 0.0,
                 'rates': 0.0, 'stoich': 0.0, 'other': 0.0, 'calls': 0}

        # Cache for t-dependent environment (perstep_snap=False only).
        # T, HV, ENZ columns and k-values depend only on t, not on y.
        # Newton iterations within one solver step share the same t → reuse.
        _env_cache_t = [float('nan')]  # mutable container for closure
        _env_cache    = {}             # cached environment dict (t-dependent only)

        # Vectorised ENZ evaluator: replace 28 individual scipy calls with one
        # numpy searchsorted + fancy-index operation per cache miss.
        _venz_names  = []   # column names (uppercase)
        _venz_x      = []   # breakpoint arrays
        _venz_y      = []   # value arrays
        _venz_kind   = []   # 'prev' or 'lin'
        _venz_scalar = []   # fallback callables for non-standard providers
        for _k, _fn in self.ctx_data.items():
            if hasattr(_fn, 'x') and hasattr(_fn, 'y'):
                _venz_names.append(_k.upper())
                _venz_x.append(np.asarray(_fn.x, dtype=np.float64))
                _venz_y.append(np.asarray(_fn.y, dtype=np.float64))
                _venz_kind.append('prev' if getattr(_fn, '_kind', '') == 'previous' else 'lin')
                _venz_scalar.append(None)
            else:
                _venz_names.append(_k.upper())
                _venz_x.append(None)
                _venz_y.append(None)
                _venz_kind.append('call')
                _venz_scalar.append(_fn)

        def _eval_ctx_data_fast(t_q):
            """Evaluate all ctx_data providers at t_q using vectorised numpy."""
            result = {}
            for _i, _name in enumerate(_venz_names):
                if _venz_kind[_i] == 'call':
                    try: result[_name] = float(_venz_scalar[_i](t_q))
                    except: pass
                else:
                    _xs = _venz_x[_i]
                    _ys = _venz_y[_i]
                    _idx = np.searchsorted(_xs, t_q, side='right') - 1
                    _idx = max(0, min(_idx, len(_ys) - 1))
                    if _venz_kind[_i] == 'prev':
                        result[_name] = float(_ys[_idx])
                    else:
                        if _idx < len(_xs) - 1:
                            _frac = (t_q - _xs[_idx]) / (_xs[_idx+1] - _xs[_idx])
                            result[_name] = float(_ys[_idx] + _frac * (_ys[_idx+1] - _ys[_idx]))
                        else:
                            result[_name] = float(_ys[_idx])
            return result

        def system_ode(t, y):
            _t0 = time.perf_counter()
            y = np.maximum(y, 0.0)
            t_abs = t + _t_offset
            dydt = np.zeros(self.solver.num_species)
            _prof['calls'] += 1

            ctx = SafeDict()
            if _perstep_snap:
                # FACSIMILE-exact: ALL settings frozen from step start (snap).
                # Only c-species update from ODE y.
                _t1 = time.perf_counter()
                ctx.update(_s_snap)
                ctx.update(zip(_species_names, y.tolist()))
                for k in _snap_active:
                    if k.upper() not in _species_set_upper:
                        continue
                    try:
                        ctx[k] = float(self.ctx_data[k.lower() if k.lower() in self.ctx_data else k](t_abs))
                    except Exception:
                        pass
                _prof['snap_read'] += time.perf_counter() - _t1
            else:
                # Hybrid: settings evaluated live. Cache t-dependent env so Newton
                # iterations at the same t reuse it without recomputing.
                if t_abs != _env_cache_t[0]:
                    _t1 = time.perf_counter()
                    M_now = self.get_M(t_abs)
                    _ec = _eval_ctx_data_fast(t_abs)
                    _temp_val = _ec.get('TEMP', _ec.get('T', 298.15))
                    T_now = float(_temp_val) if _temp_val else 298.15
                    _ec.update({"T": T_now, "TEMP": T_now, "M": M_now,
                                "np": np, "EXP": np.exp, "LOG": np.log,
                                "LOG10": np.log10, "AMAX": np.maximum, "AMIN": np.minimum,
                                "ABS": np.abs, "SQRT": np.sqrt,
                                "COS": np.cos, "SIN": np.sin, "TAN": np.tan})
                    _prof['ctx_data'] += time.perf_counter() - _t1
                    _t1 = time.perf_counter()
                    if self._s_vec:
                        from evaluator import update_enz_vec as _uenz, build_ctx as _bctx
                        from mechanism import CEntry
                        _uenz(self._enz_vec, self._enz_providers, t_abs)
                        if _svec_njit_fn is not None and _svals_buf is not None:
                            # @njit path: fill ENZ values, call compiled fn, read back
                            for _si, _sn in enumerate(_svec_enz_nms):
                                for _ev in self._enz_vec:
                                    if _ev.name.upper() == _sn:
                                        _svals_buf[_si] = _ev.value
                                        break
                                else:
                                    _svals_buf[_si] = _ec.get(_sn, 0.0)
                            _svec_njit_fn(_svals_buf)
                            for _si, _sn in enumerate(_svec_s_nms):
                                _ec[_sn] = _svals_buf[_svec_N_enz + _si]
                        elif getattr(self, '_eval_s_vec_batch', None) is not None:
                            _ctx_batch = SafeDict(_ec)
                            for _e in self._enz_vec:
                                _ctx_batch[_e.name.upper()] = _e.value
                                _ctx_batch[_e.name.lower()] = _e.value
                            self._eval_s_vec_batch(_ctx_batch)
                            _ec.update(_ctx_batch)
                        else:
                            from evaluator import update_s_vec as _us
                            _c_vec_h = [CEntry(name=sp, value=0.0)
                                        for sp in self.solver.species]
                            _us(self._s_vec, self._function_entries, _c_vec_h, self._enz_vec)
                            _ec.update(_bctx(self._s_vec, _c_vec_h, self._enz_vec))
                    else:
                        _ctx_tmp = SafeDict(_ec)
                        for _p in range(3):
                            for vn, ce in self.var_assignments:
                                if vn in self._ctx_data_protected:
                                    continue
                                try: _ec[vn] = eval(ce, globals(), _ctx_tmp)
                                except:
                                    if _p == 0: _ec[vn] = 0.0
                    _prof['elog'] += time.perf_counter() - _t1
                    _env_cache.clear()
                    _env_cache.update(_ec)
                    _env_cache_t[0] = t_abs
                ctx.update(_env_cache)
                ctx.update(zip(_species_names, y.tolist()))

            if 'H2O' in _snap_active and 'H2O' in _species_set_upper:
                h2o_vmr = ctx['H2O']
                M_now_s = float(_s_snap.get('M', self.get_M(t_abs)))
                T_now_s = float(_s_snap.get('T', _s_snap.get('TEMP', 298.15)))
                ctx['H2O'] = h2o_vmr * M_now_s
                if h2o_vmr > 0:
                    p_sat = np.exp(21.36469 - 5339.66 / T_now_s)
                    press = float(_s_snap.get('PRESS', _s_snap.get('P', 1013.25)))
                    ctx['RH'] = h2o_vmr * press / p_sat * 100.0
                else:
                    ctx['RH'] = 0.0

            # O1D PSS: use measured JO1D (k[O3+HV→O1D]) not theoretical J1.
            if 'O1D' in _pss_species_set:
                ctx['O1D'] = float(_pss_o1d(ctx))

            _t1 = time.perf_counter()
            if _k_row_cacheable and _k_row_valid[0]:
                np.copyto(_k_buf, _k_row_buf)
            else:
                try:
                    _compute_rates(ctx, _k_buf)
                except Exception:
                    for i, rxn in enumerate(self.mechanism.reactions):
                        try: _k_buf[i] = eval(rxn.safe_expr, globals(), ctx)
                        except: _k_buf[i] = 0.0
                np.nan_to_num(_k_buf, nan=0.0, posinf=0.0, neginf=0.0, copy=False)
            _prof['rates'] += time.perf_counter() - _t1

            _t1 = time.perf_counter()
            if _use_numba:
                for _ni, _nn in enumerate(_env_r_unique_names):
                    _env_vals_buf[_ni] = ctx[_nn]
                _dydt_stoich_buf[:] = 0.0
                _stoich_eval_numba(
                    _k_buf, y, _env_vals_buf,
                    _ode_r_rxn, _ode_r_sp, _ode_r_coeff,
                    _env_r_rxn, _env_r_name_idx, _env_r_coeff,
                    _stoich_rxn, _stoich_sp, _stoich_sc,
                    _v_buf_stoich, _dydt_stoich_buf,
                )
                dydt += _dydt_stoich_buf
            else:
                for i, rxn in enumerate(self.mechanism.reactions):
                    k = float(_k_buf[i])
                    if k != 0.0:
                        v = k
                        for r in rxn.reactants_parsed: v *= (y[r['idx']] ** r['coeff'])
                        for r in rxn.reactants_env: v *= (ctx.get(r['name'], 0.0) ** r['coeff'])
                        for r in rxn.reactants_parsed: dydt[r['idx']] -= (v * r['coeff'])
                        for p in rxn.products_parsed: dydt[p['idx']] += (v * p['coeff'])
            _prof['stoich'] += time.perf_counter() - _t1
            _t1 = time.perf_counter()   # checkpoint: start of "other" section

            # Analytical dilution: applied only when explicit +DIL reactions are
            # NOT already handling it. We detect this at setup via _has_explicit_dil.
            if self.cfg.get('analytical_dilution', True) and not self._has_explicit_dil:
                _dil_switch = float(ctx.get('DIL', 1.0))   # DIL=0 disables dilution
                if _dil_switch != 0.0:
                    _flow = self.ctx_data['FLOW'](t_abs) if 'FLOW' in self.ctx_data else 0.0
                    k_dil = (_flow + 0.6) / (self.cfg['chamber_vol'] * 3600)
                    dydt -= y * k_dil * _dil_switch

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

            if _held_idx.size:
                dydt[_held_idx] = 0.0

            _prof['other'] += time.perf_counter() - _t1

            self.eval_count += 1
            if self.eval_count % 500 == 0 or t > self.last_progress_print_time + 300:
                progress = max(0.0, min(t / _t_duration, 1.0))

                eta_str = ""
                if progress > 0:
                    elapsed = time.time() - self.real_start_time
                    total_est = elapsed / progress
                    rem_secs = int(total_est - elapsed)
                    eta_str = f" | ETA: {rem_secs // 60}m {rem_secs % 60}s"

                _n = max(_prof['calls'], 1)
                _pt = (_prof['snap_read'] + _prof['ctx_data'] + _prof['elog']
                       + _prof['rates'] + _prof['stoich'] + _prof['other'])
                _prof_str = (f" | snap={_prof['snap_read']/_n*1e3:.2f} "
                             f"ctx={_prof['ctx_data']/_n*1e3:.2f} "
                             f"elog={_prof['elog']/_n*1e3:.2f} "
                             f"rates={_prof['rates']/_n*1e3:.2f} "
                             f"stoich={_prof['stoich']/_n*1e3:.2f} ms/call")
                print(f"  ... ODE Eval Step: {self.eval_count} | Sim Time: {t:.1f}s / {_t_duration:.1f}s | Fraction Done: {progress:.3f} ({progress*100:.1f}%){eta_str}{_prof_str}", flush=True)
                if t > self.last_progress_print_time + 300:
                    self.last_progress_print_time = t

            return dydt

        # Solver method and Jacobian mode from scenario JSON (set by GUI dropdowns).
        solver_method = solver_method_cont = self.cfg.get('solver_method', 'BDF')
        _jac_mode = self.cfg.get('jac_mode', 'analytic')
        # Solver settings: FACS[HMAX] from EASY takes precedence over TOML defaults
        _hmax_easy = self.cfg.get('hmax', None)
        try:
            import tomllib as _tl
            _tp = getattr(self, '_toml_path', None) or \
                  os.path.join(os.path.dirname(os.path.abspath(__file__)), "FR_GUI.toml")
            with open(_tp, 'rb') as _fh:
                _solver_toml = _tl.load(_fh).get('solver', {})
        except Exception:
            _solver_toml = {}
        # Solver profile (fast | precise): a JSON written by scenario_prep carries rtol, perstep_warm and
        # solver_profile.  An older JSON without perstep_warm gets the TOML's profile, so rtol and the
        # restart mode always come as a validated pair (rtol 1e-6 with warm continuation diverges).
        if 'perstep_warm' not in self.cfg:
            try:
                from scenario_prep import resolve_profile as _resolve_profile
                _prof = _resolve_profile(_solver_toml)
                if 'rtol' in self.cfg and float(self.cfg['rtol']) != _prof['rtol']:
                    _prof['perstep_warm'] = float(self.cfg['rtol']) > 1e-5    # keep the JSON's rtol, pair it safely
                    _prof['solver_profile'] = 'custom'
                    _prof['rtol'] = float(self.cfg['rtol'])
                self.cfg.update(_prof)
            except Exception as _e:
                print(f"  [SOLVER] could not resolve the solver profile: {_e}")
        print(f"  [SOLVER] profile={self.cfg.get('solver_profile', '?')}")
        step_limit = float(_hmax_easy) if _hmax_easy is not None \
                     else float(_solver_toml.get('hmax', 60.0))
        # rtol from JSON (set by GUI) takes precedence over TOML default
        _rtol_json = self.cfg.get('rtol', None)
        if _rtol_json is not None:
            _rtol = float(_rtol_json)
        else:
            _rc = _solver_toml.get('rtol_choices', [1e-2])
            if isinstance(_rc, str):
                import ast as _ast
                try: _rc = _ast.literal_eval(_rc)
                except Exception: _rc = [1e-2]
            _rtol = float(_rc[0]) if _rc else 1e-2
        # perstep_snap: False → FACSIMILE-konform (default): BDF restarts only at ENZ
        #                       boundaries where settings actually change (~7 restarts).
        #                       Matches FACSIMILE's architecture: k-values updated only at
        #                       ENZ boundaries, not within a segment. ~15s for scene_03.
        #               True  → Sub-step mode: settings frozen + BDF restart at EVERY
        #                       10-s output step (2520 restarts). Not what FACSIMILE does.
        _perstep_snap = bool(self.cfg.get('perstep_snap',
                             _solver_toml.get('perstep_snap', False)))
        # perstep_warm: True → per-ENZ-row WITHOUT cold restarts.  One scipy BDF
        #               instance is driven row by row (t_bound moved to the next
        #               ENZ boundary); the step-size/order history is kept across
        #               boundaries exactly like GEAR in FACSIMILE.  Settings are
        #               still frozen per row via _make_snap; the Jacobian is
        #               refreshed once at every boundary.  Avoids the ~13 LU
        #               refactorisations per row caused by ramping h from
        #               first_step on every cold restart.
        def _tobool(v, default=False):
            if isinstance(v, bool): return v
            if v is None: return default
            return str(v).strip().lower() in ('1', 'true', 'yes', 'on')
        _perstep_warm = _tobool(self.cfg.get('perstep_warm',
                                _solver_toml.get('perstep_warm', False)))
        _mode_lbl = ('per-ENZ-row (warm)' if (_perstep_snap and _perstep_warm)
                     else 'per-ENZ-row' if _perstep_snap else 'continuous')
        print(f"  [SOLVER] mode={_mode_lbl}  "
              f"solver={solver_method}  jac={_jac_mode}")

        # ── Per-row rate-constant cache (per-ENZ-row modes only) ─────────────────
        # Inside a row every setting is frozen by the snapshot, so k depends on y
        # only if a rate expression references an ODE species by name (or a PSS
        # species such as O1D).  If none does, k is evaluated ONCE per row in
        # _make_snap and copied into _k_buf on every RHS / Jacobian call.
        _k_row_dep = (set(_rate_co_names) & _species_set_upper) |                      (set(_rate_co_names) & {s.upper() for s in getattr(self, '_pss_species', set())})
        _k_row_cacheable = (bool(_perstep_snap) and not _k_row_dep
                            and _tobool(self.cfg.get('k_row_cache', True), True))
        _k_row_buf   = np.zeros(len(self.mechanism.reactions))
        _k_row_valid = [False]
        if _perstep_snap:
            print(f"  [SOLVER] per-row k cache: {'ON' if _k_row_cacheable else 'OFF'}"
                  + (f" (rate expressions depend on {sorted(_k_row_dep)[:6]})" if _k_row_dep else ""))

        _atol = float(_solver_toml.get('atol', 1e4))
        # Newton iteration control for the per-row modes.  scipy's defaults
        # (newton_tol = max(10*eps/rtol, min(0.03, sqrt(rtol))), 4 iterations)
        # judge convergence by an RMS norm over all species, which lets trace
        # species near atol through under-converged once the Jacobian is exact
        # (verified on scene_05: HCOCO3H 200x off, OH -0.2 %).  1e-6 / 10
        # reproduces a fully converged reference to <1e-4 at no runtime cost.
        _newton_tol   = float(_solver_toml.get('newton_tol', 1e-6))
        _newton_maxit = int(float(_solver_toml.get('newton_maxiter', 10)))
        import scipy.integrate._ivp.bdf as _bdf_mod
        from scipy.integrate import BDF as _BDFBase
        _bdf_mod.NEWTON_MAXITER = _newton_maxit      # module constant in scipy

        class _BDFNewton(_BDFBase):
            """scipy BDF with a configurable Newton tolerance."""
            def __init__(self, *a, **k):
                super().__init__(*a, **k)
                self.newton_tol = _newton_tol
        if _perstep_snap:
            print(f"  [SOLVER] Newton: tol={_newton_tol:g}  max_iter={_newton_maxit}")
        # Absolute tolerance floor for species-specific atol (option 2).
        # Per step: atol_i = max(_atol_floor, rtol * y_i).
        # This gives O(rtol) relative accuracy for every species regardless of
        # concentration, while the floor prevents pathological step-size reduction
        # for near-zero species.  Default floor = atol/100 (100 molec/cm³ for
        # atol=1e4), configurable via [solver] atol_floor in FR_GUI.toml.
        _atol_floor = float(_solver_toml.get('atol_floor', _atol * 1e-2))
        print(f"  [SOLVER] max_step={step_limit}s  rtol={_rtol}  "
              f"atol_floor={_atol_floor:.0e}  (species-specific: max(floor, rtol·y))  "
              f"({'FACS[HMAX]' if _hmax_easy else 'TOML default'})")
        N = self.solver.num_species
        _held_idx = np.array(sorted(self._held), dtype=np.int64)
        if _held_idx.size and not _perstep_snap:
            print("  [HELD] WARNING: held species need the per-row settings snapshot (perstep_snap); "
                  "they are only frozen at segment starts in this mode")

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
        if _held_idx.size:                     # d(held)/dt = 0: no Jacobian row for a held species
            _keep = ~np.isin(_jac_rows, _held_idx)
            _jac_rows, _jac_cols = _jac_rows[_keep], _jac_cols[_keep]
            _jac_factors, _jac_rxn_idxs = _jac_factors[_keep], _jac_rxn_idxs[_keep]
        _dil_idx_all = (np.setdiff1d(np.arange(N, dtype=np.int32), _held_idx.astype(np.int32))
                        if _held_idx.size else np.arange(N, dtype=np.int32))
        print(f">>> [JAC] Analytical Jacobian: {N}×{N}, {len(_jac_rows)} terms "
              f"({'passed to solver' if _jac_mode == 'analytic' else 'built but not used — FD active'})")

        def analytical_jac(t, y_raw):
            """Analytical Jacobian using the pre-computed settings snapshot."""
            y = np.maximum(y_raw, 0.0)
            t_abs = t + _t_offset
            ctx_j = SafeDict()
            # Same fully-frozen approach as system_ode: all settings from snap
            ctx_j.update(_s_snap)
            ctx_j.update(zip(_species_names, y.tolist()))
            # ENZ-constrained ODE species only
            for k in _snap_active:
                if k.upper() not in _species_set_upper:
                    continue
                try: ctx_j[k] = float(self.ctx_data[k.lower() if k.lower() in self.ctx_data else k](t_abs))
                except: pass
            if 'H2O' in _snap_active and 'H2O' in _species_set_upper:
                h2o_vmr = ctx_j['H2O']
                M_now = float(_s_snap.get('M', self.get_M(t_abs)))
                ctx_j['H2O'] = h2o_vmr * M_now
            # O1D PSS — dynamic from mechanism reactions
            if 'O1D' in _pss_species_set:
                ctx_j['O1D'] = float(_pss_o1d(ctx_j))
            ctx_j['AMAX'] = np.maximum
            ctx_j['AMIN'] = np.minimum
            # 8. Pass 2: J-value var_assignments (INPUT variables protected)
            for var_name, comp_expr in self.var_assignments:
                if var_name in self._ctx_data_protected:
                    continue
                try:
                    ctx_j[var_name] = eval(comp_expr, globals(), ctx_j)
                except Exception:
                    ctx_j[var_name] = 0.0
            # 9b. Recompute O1D PSS (Jacobian must match system_ode)
            if 'O1D' in _pss_species_set:
                ctx_j['O1D'] = float(_pss_o1d(ctx_j))

            # 9c. O atom is kept in the ODE state (not PSS) — see note in _build_mechanism.

            # 10. Rate constants: per-row cache when available, else batch function
            if _k_row_cacheable and _k_row_valid[0]:
                np.copyto(_k_buf, _k_row_buf)
            else:
                try:
                    _compute_rates(ctx_j, _k_buf)
                except Exception:
                    for i, rxn in enumerate(self.mechanism.reactions):
                        try: _k_buf[i] = eval(rxn.safe_expr, globals(), ctx_j)
                        except: _k_buf[i] = 0.0
                np.nan_to_num(_k_buf, nan=0.0, posinf=0.0, neginf=0.0, copy=False)
            # 11. Reaction velocities: v[i] = k[i] * prod(y_ode^c) * prod(env^c)
            if _use_numba:
                # Same kernel as the RHS; dydt scratch output is discarded.
                for _ni, _nn in enumerate(_env_r_unique_names):
                    _env_vals_buf[_ni] = ctx_j[_nn]
                _dydt_stoich_buf[:] = 0.0
                _stoich_eval_numba(
                    _k_buf, y, _env_vals_buf,
                    _ode_r_rxn, _ode_r_sp, _ode_r_coeff,
                    _env_r_rxn, _env_r_name_idx, _env_r_coeff,
                    _stoich_rxn, _stoich_sp, _stoich_sc,
                    _v_buf, _dydt_stoich_buf,
                )
            else:
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
                dil_idx  = _dil_idx_all
                dil_data = np.full(len(dil_idx), diag_rate)
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

        # --- FILTER step_change_times WITH EFFECTIVE VALUES ---
        # load_saphir_environment computed step_change_times from raw ENZ values.
        # After CONST overrides (e.g. HV=CONST(1)), some of those "changes" are
        # no longer real — the effective value stays constant.  Filter them out now
        # using ctx_data, which already incorporates the CONST overrides.
        # Only keep a time t if at least one ENZ-sourced control variable actually
        # changes its effective value at t (checked via ctx_data interpolators).
        # Which ctx_data keys to check for effective value changes:
        # T/P (k-values + thermal expansion) and all step-pattern columns (HV, Q*).
        import fnmatch as _fnm
        _step_pats = _load_step_col_patterns()
        # T/P only need BDF restarts when they are FIXED CONSTANT overrides (an override
        # creates a genuine step change at ENZ boundaries relative to a measured ENZ T/P).
        # When T or P comes from the ENZ (SAPHIR File), it varies smoothly row-to-row;
        # the BDF handles that without restarts.  Including ENZ-sourced T here would
        # generate 1000+ spurious restarts (one per 60 s ENZ row), turning scene_04 into
        # 1078 short segments where the HV-on phase immediately stalls.
        _em = self.cfg.get('env_map', {})
        _t_type = _em.get('TEMP', _em.get('T', {})).get('type', '')
        _p_type = _em.get('PRESS', _em.get('P', {})).get('type', '')
        # T/P trigger restarts only when NOT from ENZ and NOT Fixed Constant.
        # Fixed Constant: k-values don't change at ENZ boundaries → no restart needed.
        # SAPHIR File: T varies smoothly row-to-row → BDF handles without restart.
        _t_no_restart = _t_type in ('SAPHIR File', 'Fixed Constant')
        _p_no_restart = _p_type in ('SAPHIR File', 'Fixed Constant')
        _RESTART_KEYS = {
            k for k in self.ctx_data
            if (k.upper() in {'T', 'TEMP'} and not _t_no_restart)
            or (k.upper() in {'P', 'PRESS'} and not _p_no_restart)
            or any(_fnm.fnmatch(k.upper(), p.upper()) for p in _step_pats)
        }
        # ENZ breakpoints are shifted +1e-9 s in load_saphir_environment so that
        # exactly at the boundary time t, the interpolator returns the OLD value.
        # IMPORTANT: ENZ timestamps are absolute (seconds since 1.1.2000, ~7.74e8).
        # At that scale, float64 ULP ≈ 119 ns.  1 ns offsets round to zero and both
        # queries land on the same float64 value → filter incorrectly sees no change.
        # Use 1e-4 (100 µs) >> 119 ns ULP, and << 60 s ENZ row spacing.
        _eff_eps   = 1e-4   # left-hand query offset (must exceed float64 ULP at ENZ scale)
        _rh_offset = 1e-4   # right-hand query offset (same reasoning)
        raw_sct = getattr(self, '_step_change_times', [])
        filtered_sct = []
        _sct_triggers = []   # which key triggered each boundary (for diagnostic print)
        for _t in raw_sct:
            for _k in _RESTART_KEYS:
                _vfn = self.ctx_data.get(_k)
                if _vfn is None:
                    continue
                try:
                    _v_before = float(_vfn(_t - _eff_eps))
                    _v_after  = float(_vfn(_t + _rh_offset))
                    _scale    = max(abs(_v_before), abs(_v_after), 1e-30)
                    if abs(_v_after - _v_before) / _scale > 1e-8:
                        filtered_sct.append(_t)
                        _sct_triggers.append(_k)
                        break
                except Exception:
                    pass
        if len(filtered_sct) != len(raw_sct) and not _perstep_snap:
            print(f"  [SOLVER] step_change_times: {len(raw_sct)} raw "
                  f"-> {len(filtered_sct)} effective "
                  f"({len(raw_sct)-len(filtered_sct)} suppressed by CONST overrides)")
        if filtered_sct and not _perstep_snap:
            import datetime as _dt
            _epoch = _dt.datetime(2000, 1, 1)
            _show = min(len(filtered_sct), 50)
            print(f"  [SOLVER] Restart boundaries ({len(filtered_sct)}"
                  f"{f', showing first {_show}' if len(filtered_sct) > _show else ''}):")
            for _t, _k in zip(filtered_sct[:_show], _sct_triggers[:_show]):
                _ts = (_epoch + _dt.timedelta(seconds=_t)).strftime("%Y-%m-%d %H:%M")
                print(f"           {_ts}  trigger={_k}")
        self._step_change_times = filtered_sct

        _eps = 1e-6
        event_times = [float(t) - _eps for t in self._step_change_times
                       if self.t_start < t < self.t_end]

        if not _perstep_snap:
            if event_times:
                print(f">>> [SOLVER] {len(event_times)} segment boundary/ies "
                      f"detected — integration split at step changes.")
            else:
                print(f">>> [SOLVER] No step changes detected — single continuous segment.")

        # Build list of (t0, t1) segment boundaries
        boundaries = sorted(set([self.t_start] + event_times + [self.t_end]))
        # Remove duplicate / too-close boundaries
        boundaries = [boundaries[0]] + [b for i, b in enumerate(boundaries[1:], 1)
                                         if b - boundaries[i - 1] > 1.0]

        _n_segs = len(self.t_span)-1 if _perstep_snap else len(boundaries)-1
        print(f">>> [SOLVER] {len(self.mechanism.species)} species  "
              f"solver={solver_method}  segments={_n_segs}  max_step={step_limit}s")

        # ── Integration loop ──────────────────────────────────────────────────────
        seg_t_list = [np.array([self.t_span[0]])]
        seg_y_list = [y0[:, np.newaxis]]
        y_seg = y0.copy()

        # first_step: JSON (GUI combo) takes precedence over TOML default.
        _fs_json = self.cfg.get('first_step', None)
        if _fs_json is not None:
            _fs_val = float(_fs_json)
        else:
            _fs_val = float(_solver_toml.get('first_step', 0.0))
        _fs_kw  = {'first_step': _fs_val} if _fs_val > 0 else {}
        if _fs_kw:
            if _perstep_snap:
                print(f"  [SOLVER] first_step={_fs_val:.2e}s (from TOML, applied per restart)")
            else:
                print(f"  [SOLVER] first_step={_fs_val:.2e}s (from TOML, ignored in continuous mode)")

        if _perstep_snap and _perstep_warm:
            # ── Per-ENZ-row mode, warm continuation (GEAR-like) ─────────────────
            # One BDF instance for the whole run.  At every ENZ row boundary:
            #   1. output y(t1)                     (pre-expansion, FACSIMILE convention)
            #   2. thermal expansion y *= M_new/M_old, also applied to the BDF
            #      difference array D (pure rescaling keeps the history consistent)
            #   3. new settings snapshot (_make_snap) at t1
            #   4. Jacobian refreshed, LU invalidated, t_bound moved to the next row
            # Step size and BDF order are NOT reset → no order-1 ramp, no LU churn.
            _BDFSolver = _BDFNewton
            n_out_steps = len(self.t_span) - 1
            _wsolver = None
            _n_warm_fallback = 0
            _n_soft_restart  = 0
            # Real discontinuities (HV / Q* / overridden T,P changes, as filtered
            # above): at these rows the RHS jumps and the BDF history is invalid,
            # so a cold restart (first_step, order 1) is done IN PLACE — exactly
            # what the cold mode does at every row, but only here (~7 times).
            _sct_arr = np.array(sorted(self._step_change_times), dtype=float)
            def _is_discontinuity(t_abs):
                if _sct_arr.size == 0: return False
                _i = np.searchsorted(_sct_arr, t_abs)
                for _j in (_i - 1, _i):
                    if 0 <= _j < _sct_arr.size and abs(_sct_arr[_j] - t_abs) < 0.5:
                        return True
                return False
            for step_idx in range(n_out_steps):
                t0 = self.t_span[step_idx]
                t1 = self.t_span[step_idx + 1]
                _make_snap(y_seg, t0)
                if _held_idx.size and _wsolver is not None:
                    _wsolver.y[_held_idx] = y_seg[_held_idx]
                    _wsolver.D[0, _held_idx] = y_seg[_held_idx]
                    _wsolver.D[1:, _held_idx] = 0.0
                _atol_vec = np.maximum(_atol_floor, _rtol * np.maximum(y_seg, 0.0))
                if _wsolver is not None and _is_discontinuity(t0):
                    # Soft (in-place) cold restart: same state as a fresh BDF().
                    _h0 = _fs_val if _fs_val > 0 else min(_wsolver.h_abs, 1e-3)
                    _f0 = system_ode(_wsolver.t, _wsolver.y)
                    _wsolver.D[:] = 0.0
                    _wsolver.D[0] = _wsolver.y
                    _wsolver.D[1] = _f0 * _h0 * _wsolver.direction
                    _wsolver.h_abs = _h0
                    _wsolver.order = 1
                    _wsolver.n_equal_steps = 0
                    _n_soft_restart += 1
                if _wsolver is None:
                    _wsolver = _BDFSolver(
                        system_ode, t0 - _t_offset, y_seg,
                        t_bound=t1 - _t_offset, max_step=step_limit,
                        rtol=_rtol, atol=_atol_vec,
                        jac=analytical_jac if _jac_mode == 'analytic' else None,
                        **_fs_kw)
                else:
                    _wsolver.t_bound = t1 - _t_offset
                    _wsolver.status  = 'running'
                    _wsolver.atol    = _atol_vec
                    # Settings changed at the boundary → Jacobian/LU are stale.
                    _wsolver.J  = _wsolver.jac(_wsolver.t, _wsolver.y)
                    _wsolver.LU = None
                while _wsolver.status == 'running':
                    _msg = _wsolver.step()
                if _wsolver.status == 'failed':
                    # Fallback: finish this row with a cold solve_ivp, then
                    # re-create the warm solver at t1 on the next iteration.
                    _n_warm_fallback += 1
                    print(f"\n[WARNING] Warm step {step_idx+1}/{n_out_steps} "
                          f"({t0:.0f}-{t1:.0f}s) failed: {_msg} — cold fallback")
                    _fb = solve_ivp(
                        system_ode, [_wsolver.t, t1 - _t_offset], _wsolver.y,
                        method=solver_method, t_eval=[t1 - _t_offset],
                        max_step=step_limit, rtol=_rtol, atol=_atol_vec,
                        jac=analytical_jac if _jac_mode == 'analytic' else None,
                        **_fs_kw)
                    y_seg = _fb.y[:, -1] if _fb.y.shape[1] > 0 else _wsolver.y.copy()
                    _wsolver = None
                else:
                    y_seg = _wsolver.y.copy()
                if not np.all(np.isfinite(y_seg)):
                    y_seg = np.where(np.isfinite(y_seg), y_seg, 0.0)
                    if _wsolver is not None:
                        _wsolver.y = y_seg.copy()
                        _wsolver.D[0] = y_seg
                seg_t_list.append(np.array([t1]))
                seg_y_list.append(y_seg[:, np.newaxis])
                M_snap_step = float(_s_snap.get('M', self.get_M(t0)))
                M_end = self.get_M(t1)
                if M_snap_step > 0 and abs(M_end / M_snap_step - 1.0) > 1e-10:
                    _xf = M_end / M_snap_step
                    y_seg = y_seg * _xf
                    if _wsolver is not None:
                        _wsolver.y = y_seg.copy()
                        _wsolver.D *= _xf          # D[0] == y, higher rows scale alike
                if step_idx % 200 == 0 or step_idx == n_out_steps - 1:
                    pct = (step_idx + 1) / n_out_steps * 100
                    elapsed = time.time() - self.real_start_time
                    print(f"  ... {pct:.0f}% ({step_idx+1}/{n_out_steps} steps, "
                          f"{elapsed:.1f}s elapsed)", end='\r')
            if _wsolver is not None:
                print(f"\n  [SOLVER] warm mode: {_wsolver.nfev} RHS evals, "
                      f"{_wsolver.njev} Jacobians, {_wsolver.nlu} LU factorisations, "
                      f"{_n_soft_restart} soft restarts at discontinuities, "
                      f"{_n_warm_fallback} cold fallbacks")
        elif _perstep_snap:
            # ── Per-ENZ-row mode: one solve_ivp per output step ──────────────────
            # Settings frozen from snap at each step start.
            n_out_steps = len(self.t_span) - 1
            for step_idx in range(n_out_steps):
                t0 = self.t_span[step_idx]
                t1 = self.t_span[step_idx + 1]
                _make_snap(y_seg, t0)
                _atol_vec = np.maximum(_atol_floor, _rtol * np.maximum(y_seg, 0.0))
                seg_sol = solve_ivp(
                    system_ode, [t0 - _t_offset, t1 - _t_offset], y_seg,
                    method=_BDFNewton if solver_method == 'BDF' else solver_method,
                    t_eval=[t0-_t_offset, t1-_t_offset],
                    max_step=step_limit, rtol=_rtol, atol=_atol_vec,
                    jac=analytical_jac if _jac_mode == 'analytic' else None,
                    **_fs_kw)
                if not seg_sol.success:
                    print(f"\n[WARNING] Step {step_idx+1}/{n_out_steps} "
                          f"({t0:.0f}-{t1:.0f}s): {seg_sol.message}")
                _sol_y = seg_sol.y
                if isinstance(_sol_y, list):
                    _sol_y = np.column_stack(_sol_y) if _sol_y else np.empty((len(y_seg), 1))
                y_seg = _sol_y[:, -1] if _sol_y.shape[1] > 0 else y_seg
                if not np.all(np.isfinite(y_seg)):
                    y_seg = np.where(np.isfinite(y_seg), y_seg, 0.0)
                seg_t_list.append(np.array([t1]))
                seg_y_list.append(y_seg[:, np.newaxis])
                M_snap_step = float(_s_snap.get('M', self.get_M(t0)))
                M_end = self.get_M(t1)
                if M_snap_step > 0 and abs(M_end / M_snap_step - 1.0) > 1e-10:
                    y_seg = y_seg * (M_end / M_snap_step)
                if step_idx % 200 == 0 or step_idx == n_out_steps - 1:
                    pct = (step_idx + 1) / n_out_steps * 100
                    elapsed = time.time() - self.real_start_time
                    print(f"  ... {pct:.0f}% ({step_idx+1}/{n_out_steps} steps, "
                          f"{elapsed:.1f}s elapsed)", end='\r')
        else:
            # ── Hybrid mode: one solve_ivp per ENZ segment (~7 restarts) ────────
            # Settings re-evaluated live within each segment; snap used only for
            # thermal expansion at segment boundaries and for the -1ns export shift.
            n_segs = len(boundaries) - 1
            for seg_idx in range(n_segs):
                t_seg0 = boundaries[seg_idx]
                t_seg1 = boundaries[seg_idx + 1]
                # Snapshot at segment start so analytical_jac has valid k-values.
                # The ODE itself evaluates settings live; the snap is only used by
                # analytical_jac as an approximation (correct at t_seg0, slightly
                # stale within the segment — still far better than an empty snap).
                if _jac_mode == 'analytic':
                    _make_snap(y_seg, t_seg0)
                _M_seg_start = self.get_M(t_seg0)
                t_eval_seg = self.t_span[(self.t_span >= t_seg0) & (self.t_span <= t_seg1)]
                if len(t_eval_seg) < 2 or t_eval_seg[-1] < t_seg1 - 0.5:
                    t_eval_seg = np.union1d(t_eval_seg, [t_seg1])
                _atol_vec = np.maximum(_atol_floor, _rtol * np.maximum(y_seg, 0.0))
                seg_sol = solve_ivp(
                    system_ode,
                    [t_seg0 - _t_offset, t_seg1 - _t_offset], y_seg,
                    method=solver_method_cont,
                    t_eval=t_eval_seg - _t_offset,
                    max_step=step_limit, rtol=_rtol, atol=_atol_vec,
                    jac=analytical_jac if _jac_mode == 'analytic' else None)
                if not seg_sol.success:
                    print(f"\n[WARNING] Segment {seg_idx+1}/{n_segs} "
                          f"({t_seg0:.0f}-{t_seg1:.0f}s): {seg_sol.message}")
                # Skip first point of every segment: the initialization already
                # provides t=t_span[0] for seg_idx=0, and the previous segment's
                # last point is the boundary for seg_idx>0.  Without this skip,
                # t=t_span[0] appears twice (once from init, once from sol.t[0]),
                # the second occurrence with k=0 in the kv export — a spurious artefact.
                # Defensive: if solver returns y as a plain list (can happen when atol
                # contains NaN from a prior stalled segment), convert to numpy array.
                _sol_y = seg_sol.y
                if isinstance(_sol_y, list):
                    _sol_y = np.column_stack(_sol_y) if _sol_y else np.empty((len(y_seg), 0))
                seg_t_list.append(seg_sol.t[1:] + _t_offset)
                seg_y_list.append(_sol_y[:, 1:])
                y_seg = _sol_y[:, -1] if _sol_y.shape[1] > 0 else y_seg
                # Sanitize: replace NaN/Inf from a stalled solver with 0 to prevent
                # cascading failures in subsequent segments.
                if not np.all(np.isfinite(y_seg)):
                    print(f"\n[WARNING] Segment {seg_idx+1}: non-finite concentrations "
                          f"after solver failure — clamping to 0 for next segment.")
                    y_seg = np.where(np.isfinite(y_seg), y_seg, 0.0)

                # Thermal expansion at segment boundary.
                # When T or P jumps discontinuously (piecewise ENZ), all concentrations
                # must be rescaled to preserve mixing ratios immediately after the jump:
                #   y_new = y_old × M(T_new) / M(T_old)
                # The ODE works in molec/cm³ and cannot apply this automatically.
                # Correct check: compare M just BEFORE the boundary (= at t_seg1) with
                # M just AFTER (= at t_seg1 + epsilon). The boundary is at t_seg1 so the
                # NEXT segment starts at T_new. Record output BEFORE expanding (matches
                # FACSIMILE convention: output uses pre-transition concentrations).
                _M_at_end   = self.get_M(t_seg1)
                # Use t_seg1 + 1e-3 (not 1e-6) to avoid floating-point cancellation:
                # t_seg1 = t_boundary - 1e-6, so (t_seg1 + 1e-6) == t_boundary exactly,
                # which is still BEFORE the +1ns-shifted ENZ breakpoint.
                _M_after    = self.get_M(t_seg1 + 1e-3)
                if _M_at_end > 0 and abs(_M_after / _M_at_end - 1.0) > 1e-10:
                    y_seg = y_seg * (_M_after / _M_at_end)

                elapsed = time.time() - self.real_start_time
                print(f"  ... segment {seg_idx+1}/{n_segs} done  "
                      f"({t_seg0:.0f}–{t_seg1:.0f}s, {elapsed:.1f}s elapsed)", end='\r')
        print()

        # ── Profiling summary ─────────────────────────────────────────────────
        _n = max(_prof['calls'], 1)
        _total = (_prof['snap_read'] + _prof['ctx_data'] + _prof['elog']
                  + _prof['rates'] + _prof['stoich'] + _prof['other'])
        print(f"  [PROF] system_ode called {_n}x  total={_total:.2f}s  per_call={_total/_n*1000:.3f}ms")
        for _k, _lbl in [('snap_read','snap/ctx_data'), ('ctx_data','ctx_data(live)'),
                          ('elog','ELOG-evaluator'), ('rates','_compute_rates'),
                          ('stoich','stoichiometry'), ('other','other')]:
            _s = _prof[_k]
            if _s > 0.001:
                print(f"         {_lbl:20s}  {_s:7.3f}s  {_s/_total*100:5.1f}%  {_s/_n*1000:.3f}ms/call")

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
        df.index.name = "TIME_s"

        # Settings reporting convention (matches FACSIMILE):
        # At output time t_i the c-values reflect the step [t_{i-1}, t_i] which was
        # computed using settings from ENZ row i-1 (read at t_{i-1}).  The s-values
        # in the output row should therefore also come from t_{i-1}, so they describe
        # exactly the settings that produced the accompanying c-values.
        # Implementation: evaluate ctx_data at sol.t[k-1] for output row k (k≥1).
        # This is robust at any timestamp scale (no float64 precision issue).
        _s_eval_t = np.concatenate([[self.sol.t[0]], self.sol.t[:-1]])

        # Calculate Air Density (M) at the (shifted) settings time
        df['M'] = [self.get_M(t) for t in _s_eval_t]

        # Build evaluation context: ODE species + M + ENZ inputs under EASY names.
        # Output mirrors the FACSIMILE .cs structure: all c (concentrations) and
        # s (settings) under their plain EASY-script names — no SAPHIR_ prefix.
        _eval_ctx = {col.upper(): df[col].values.astype(float) for col in df.columns}
        _eval_ctx.update({'np': np, 'EXP': np.exp, 'LOG': np.log,
                          'LOG10': np.log10, 'AMAX': np.maximum, 'AMIN': np.minimum,
                          'ABS': np.abs, 'SQRT': np.sqrt,
                          'COS': np.cos, 'SIN': np.sin, 'TAN': np.tan})

        # Only output columns that are defined in the EASY script (c or s entries).
        # Extra ENZ columns (TEMP, PRESS, Q_* switches) are NOT in the EASY script
        # and should not appear in the cs output — keep the output clean and
        # consistent with what FACSIMILE outputs.
        # _mechanism_s_names: s_entry names collected from mechanism_files parsing
        # (in elog mode these are already in self.mechanism.s_entries).
        _easy_names_upper = ({e.name.upper() for e in self.mechanism.s_entries
                               if not getattr(e, 'is_synthetic', False)}
                             | {s.upper() for s in self.solver.species}
                             | getattr(self, '_mechanism_s_names', set()))

        # Add ctx_data (ENZ inputs) under their EASY names (T, HV, JNO2, etc.)
        # Evaluated at _s_eval_t (shifted) so reported values match the step-start snap.
        for k, v_func in self.ctx_data.items():
            k_up = k.upper()
            if k_up == 'TIME':
                continue
            if callable(v_func):
                try:
                    vals = v_func(_s_eval_t)
                    if isinstance(vals, (float, int)):
                        vals = np.full(len(df), float(vals))
                    if k_up not in _eval_ctx:    # don't overwrite ODE species
                        _eval_ctx[k_up] = vals
                    if k_up not in df.columns and k_up in _easy_names_upper:
                        df[k_up] = vals          # only if in EASY script
                except Exception:
                    pass

        # Build export active_constraints: same alias logic as system_ode.
        # ENZ-sourced env_map variables (and their T/TEMP P/PRESS aliases) must not
        # be overwritten by CONST entries from idl_parameters in var_assignments.
        # _ctx_data_protected raw aliases (e.g. H2O_IN) are also protected.
        # Non-ELOG mode: protect env_map SAPHIR File variables from being overwritten.
        # In ELOG mode these aren't in ctx_data (no renaming), so nothing to protect.
        _export_env_guard = set()
        if not self.cfg.get('elog_path'):
            for _ek, _em in self.cfg.get('env_map', {}).items():
                if _em.get('type') == 'SAPHIR File' and _ek.upper() in _eval_ctx:
                    _export_env_guard.add(_ek.upper())
            for _a, _b in (('TEMP', 'T'), ('T', 'TEMP'), ('PRESS', 'P'), ('P', 'PRESS')):
                if _a in _export_env_guard and _b.upper() in _eval_ctx:
                    _export_env_guard.add(_b.upper())

        _ode_species_upper = {s.upper() for s in self.solver.species}

        if self._s_vec:
            # ELOG mode: evaluate s_vec in topological order using pre-compiled code.
            # Build vectorized ENZ arrays over _s_eval_t.
            for _e in self._enz_vec:
                _fn = self._enz_providers.get(_e.name)
                if _fn:
                    try:
                        _arr = _fn(_s_eval_t)
                        _e.value = _arr if isinstance(_arr, np.ndarray) else np.full(len(df), float(_arr))
                    except Exception:
                        _e.value = np.zeros(len(df))
                if _e.name.upper() not in _eval_ctx:   # don't overwrite ODE species
                    _eval_ctx[_e.name.upper()] = _e.value

            # Build a name→value map from the ENZ vec for O(1) INPUT() lookup.
            # This mirrors update_s_vec() in the ODE: INPUT entries always read
            # from the raw ENZ column, never from _eval_ctx (which can be
            # overwritten by a CONST() definition for the same name).
            _enz_export_map = {_e.name.upper(): _e.value for _e in self._enz_vec}

            for _entry in self._s_vec_ordered:
                _input_col = getattr(_entry, '_input_col', None)
                _code      = getattr(_entry, '_code', None)
                try:
                    if _input_col is not None:
                        # Always resolve INPUT from ENZ vec (consistent with ODE)
                        _raw = _enz_export_map.get(_input_col)
                        _result = (_raw if isinstance(_raw, np.ndarray) and _raw.shape == (len(df),)
                                   else np.zeros(len(df)))
                    elif _code is not None:
                        _result = eval(_code, {"__builtins__": {}}, _eval_ctx)
                    else:
                        continue
                    if isinstance(_result, (int, float)):
                        _result = np.full(len(df), float(_result))
                    if isinstance(_result, np.ndarray) and _result.shape == (len(df),):
                        _eval_ctx[_entry.name] = _result
                        if not getattr(_entry, 'is_synthetic', False):
                            df[_entry.name] = _result   # synthetic entries: ODE only, not in cs
                except Exception:
                    pass
        else:
            # Non-ELOG mode: old 3-pass var_assignments loop.
            for _pass in range(3):
                for _prot in self._ctx_data_protected:
                    if _prot in self.ctx_data:
                        try:
                            _eval_ctx[_prot] = self.ctx_data[_prot](_s_eval_t)
                        except Exception:
                            pass
                for var_name, compiled_expr in self.var_assignments:
                    if var_name in self._ctx_data_protected:
                        continue
                    if var_name in _export_env_guard and var_name not in self._idl_params_upper:
                        continue
                    if var_name in _ode_species_upper:
                        continue
                    try:
                        result = eval(compiled_expr, globals(), _eval_ctx)
                        if isinstance(result, (int, float)):
                            result = np.full(len(df), float(result))
                        if isinstance(result, np.ndarray) and result.shape == (len(df),):
                            _eval_ctx[var_name] = result
                            df[var_name] = result
                    except Exception:
                        pass

        # --- O1D pseudo-steady-state correction ------------------------------------
        # O1D has a lifetime of ~1 ns (quenched mainly by N2 and O2).  The ODE
        # solver (atol=1e4) cannot resolve it accurately — atol is ~10⁶× the
        # actual PSS value (~0.01 molec/cm³), so the integrated column carries
        # numerical noise / an initial-condition spike.  Replace with the
        # analytical PSS:  [O1D]_pss = J1·[O3] / (k_N2·[N2] + k_O2·[O2] + k_H2O·[H2O])
        # O1D: always export as PSS — either replacing an ODE column (legacy path where
        # O1D was kept in the ODE) or adding it fresh (current path where O1D is PSS-only).
        if 'O1D' in df.columns or 'O1D' in getattr(self, '_pss_species', set()):
            # Dynamic O1D PSS from mechanism reactions (vectorised over all time steps).
            _o1d_pss = self._pss_o1d(_eval_ctx)
            if isinstance(_o1d_pss, np.ndarray) and _o1d_pss.shape == (len(df),):
                pass
            else:
                _o1d_pss = np.full(len(df), float(_o1d_pss))
            df['O1D'] = _o1d_pss
            _eval_ctx['O1D'] = _o1d_pss
            print("  [O1D] PSS concentration written to output (dynamic from mechanism reactions)")

        # --- O atom export: disabled (was post-process PSS overwrite) ---------------
        # FACSIMILE keeps O in the ODE state and reports values O(~10²) molec/cm³,
        # which is ~100× smaller than the analytical PSS J4·NO2/k_O+O2+M (~10⁴).
        # The PSS overwrite produced spurious +10880% bias in the diagnostic.
        # Per commit 5da582f, JUICE keeps O in the ODE — leave the raw ODE state in
        # the export so comparisons reflect what the solver actually computes.
        if False and 'O' in df.columns:
            _T   = (df['TEMP'].values if 'TEMP' in df.columns else
                    df['T'].values    if 'T'    in df.columns else
                    np.full(len(df), 298.15))
            _M   = df['M'].values
            _HV  = df['HV'].values  if 'HV'  in df.columns else np.ones(len(df))
            _no2 = df['NO2'].values if 'NO2' in df.columns else np.zeros(len(df))
            _o3  = df['O3'].values  if 'O3'  in df.columns else np.zeros(len(df))
            _no3 = df['NO3'].values if 'NO3' in df.columns else np.zeros(len(df))
            _j4  = df['J4'].values  if 'J4'  in df.columns else np.zeros(len(df))
            _j2  = df['J2'].values if 'J2' in df.columns else np.zeros(len(df))   # O3 + hv → O
            _j6  = df['J6'].values if 'J6' in df.columns else np.zeros(len(df))   # NO3 + hv → NO2 + O
            # J4/J2/J6 are already scaled by HV from the ODE evaluation.
            # No threshold: HV=0.5 must give J*0.5, not J*0.
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

        out_path = os.path.join(self.root, f"{self._stem}.cs")
        _write_csv(out_path, df, _juice_header(self._stem, getattr(self, '_ode_backend', 'python')), file_type='cs')
        print(f"  ... Saved to {out_path}")

        # Export kv/pd only when requested ("export_kv_pd" in the scenario JSON: the Run step's
        # "Also write the reaction rates and production/loss budgets" box, run_benchmark.py --kv-pd)
        if self.cfg.get('export_kv_pd', False):
            self._export_kinetics(df)
        else:
            print("  (kv/pd export skipped — tick 'Also write the reaction rates and production/loss "
                  "budgets' on the Run step, or set \"export_kv_pd\": true in the scenario JSON)")

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

        # ── PRE-EXTRACT: species concentrations as contiguous numpy matrix ─────────
        # df[col].iloc[t_i] in a loop is ~100× slower than numpy indexing.
        # Only keep species that are actually ODE columns in df.
        _sp_cols    = [s for s in _species_names if s in df.columns]
        _sp_missing = [s for s in _species_names if s not in df.columns]
        y_matrix = np.maximum(df[_sp_cols].values.astype(np.float64), 0.0)  # (n_t, len(_sp_cols))
        _sp_idx_map = {s: i for i, s in enumerate(_sp_cols)}  # species → column in y_matrix

        # ── PRE-EXTRACT: env time series by calling interpolators vectorized ───────
        # Same 1-ns shift as _export_results: s-values reflect settings from the
        # Same FACSIMILE convention as _export_results: k-values at t_i reflect
        # the ENZ settings from t_{i-1} (the settings that produced the c-values at t_i).
        _kv_eval_t = np.concatenate([[self.sol.t[0]], self.sol.t[:-1]])
        # Non-ELOG mode: env_map SAPHIR File variables must not be overwritten.
        # In ELOG mode _ctx_data_protected already covers all ENZ aliases.
        _kv_env_guard = set()
        if not self.cfg.get('elog_path'):
            for _ke, _em in self.cfg.get('env_map', {}).items():
                if _em.get('type') == 'SAPHIR File':
                    _kv_env_guard.add(_ke.upper())

        env_arrs = {}   # key.upper() → numpy array shape (n_t,)
        for k_env, v_func in self.ctx_data.items():
            try:
                arr = v_func(_kv_eval_t)
                if not isinstance(arr, np.ndarray):
                    arr = np.full(n_t, float(arr))
                env_arrs[k_env.upper()] = arr.astype(np.float64)
            except Exception:
                pass

        # Pre-extract M and T arrays (used every step)
        M_arr = df['M'].values.astype(np.float64) if 'M' in df.columns else \
                np.array([self.get_M(t) for t in _kv_eval_t])
        T_arr = env_arrs.get('TEMP', np.full(n_t, 298.15))

        # O1D PSS array — populated during the k-loop below and injected into
        # env_arrs so the vectorised turnover loop can use it as a reactant array.
        _o1d_pss_arr = np.zeros(n_t) if 'O1D' in getattr(self, '_pss_species', set()) else None

        # Effective env-reactant arrays — filled during the k-loop from the FULL
        # evaluated ctx (after var_assignments / CONST overrides).  This guarantees
        # v = k × [A] × [B] always uses the same context as k, regardless of whether
        # a reactant is driven by ENZ, env_map, CONST override, or PSS formula.
        _all_env_reactant_names = {
            r['name']
            for rxn in self.mechanism.reactions
            for r in rxn.reactants_env
        }
        _eff_env_arrs = {name: np.zeros(n_t) for name in _all_env_reactant_names}

        # H2O array — no automatic unit conversion; EASY script defines units.
        if 'H2O' in df.columns:
            h2o_arr = df['H2O'].values.astype(np.float64)
        else:
            raise KeyError(
                "H2O is not defined in the EASY script and not available from ENZ.")

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
                        "np": np, "exp": np.exp, "EXP": np.exp,
                        "LOG10": np.log10, "log10": np.log10,
                        "LOG": np.log, "log": np.log, "SQRT": np.sqrt,
                        "COS": np.cos, "SIN": np.sin, "TAN": np.tan,
                        "AMAX": np.maximum, "AMIN": np.minimum})
            # Species from numpy (no pandas .iloc)
            for j, s in enumerate(_sp_cols):
                ctx[s] = float(y_matrix[t_i, j])
            # Evaluate settings via new vector system (ELOG) or old 3-pass (non-ELOG)
            if self._s_vec:
                from evaluator import update_s_vec as _us, build_ctx as _bctx
                from mechanism import CEntry
                _uenz_kv = getattr(self, '_enz_vec', [])
                for _e in _uenz_kv:
                    _e.value = float(env_arrs.get(_e.name.upper(), [0]*n_t)[t_i]
                                     if _e.name.upper() in env_arrs
                                     else env_arrs.get(_e.name, [0]*n_t)[t_i]
                                     if _e.name in env_arrs else 0.0)
                _c_kv = [CEntry(name=s, value=float(y_matrix[t_i, j]))
                         for j, s in enumerate(_sp_cols)]
                _us(self._s_vec, self._function_entries, _c_kv, _uenz_kv,
                    _ordered=self._s_vec_ordered)
                ctx.update(_bctx(self._s_vec, _c_kv, _uenz_kv))
            else:
                for _kv_pass in range(3):
                    for _prot in self._ctx_data_protected:
                        if _prot in env_arrs:
                            ctx[_prot] = float(env_arrs[_prot][t_i])
                    for var_name, comp_expr in self.var_assignments:
                        if var_name in self._ctx_data_protected:
                            continue
                        if var_name in _kv_env_guard and var_name not in self._idl_params_upper:
                            continue
                        try: ctx[var_name] = eval(comp_expr, globals(), ctx)
                        except:
                            if _kv_pass == 0: ctx[var_name] = 0.0

            # O1D PSS — dynamic from mechanism reactions
            if 'O1D' in getattr(self, '_pss_species', set()):
                ctx['O1D'] = float(self._pss_o1d(ctx))
                if _o1d_pss_arr is not None:
                    _o1d_pss_arr[t_i] = ctx['O1D']

            try:
                _compute_rates(ctx, _k_loc)
            except Exception:
                for i, rxn in enumerate(self.mechanism.reactions):
                    try: _k_loc[i] = eval(rxn.safe_expr, globals(), ctx)
                    except: _k_loc[i] = 0.0
            np.nan_to_num(_k_loc, nan=0.0, posinf=0.0, neginf=0.0, copy=False)
            k_arr[t_i] = _k_loc

            # Collect effective env-reactant values from the fully-evaluated ctx.
            # This is the same context used for k, so v = k × reactants is always
            # consistent — CONST overrides, PSS, and ENZ-driven values all correct.
            for _ename in _all_env_reactant_names:
                _eff_env_arrs[_ename][t_i] = max(float(ctx.get(_ename, 0.0)), 0.0)

            if (t_i + 1) % 400 == 0 or t_i == n_t - 1:
                print(f"  ... k-values {t_i + 1}/{n_t}")

        # ── VECTORISED TURNOVER COMPUTATION ───────────────────────────────────────
        # v[t,i] = k[t,i] * prod_j( y[t,j]^c_j ) * prod_e( env[t,e]^c_e )
        # All env reactants (HV, DIL, Q*, O1D, …) come from _eff_env_arrs — filled
        # during the k-loop from the fully-evaluated ctx.  No special-casing needed.

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
            # Env reactants: always use _eff_env_arrs — collected from the fully
            # evaluated ctx during the k-loop (same var_assignments / CONST overrides).
            # v = k × [A] × [B] is therefore always consistent with how k was computed.
            for r in rxn.reactants_env:
                r_name = r['name']
                env_col = _eff_env_arrs.get(r_name)
                if env_col is None:
                    if r_name == 'M':
                        env_col = M_arr
                    elif r_name == 'H2O':
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

        # ── SAVE NPZ (internal format) ────────────────────────────────────────────
        npz_path = os.path.join(self.root, f"{self._stem}_kinetics.npz")
        np.savez_compressed(npz_path,
                            time_s=self.sol.t.astype(np.float64),
                            k_values=k_arr.astype(np.float32),
                            turnovers=v_arr.astype(np.float32),
                            reaction_labels=np.array(rxn_labels))
        print(f"  ... Kinetics saved to {npz_path}  "
              f"({n_rxn} reactions × {n_t} timesteps, float32)")

        man_path = os.path.join(self.root, f"{self._stem}_reactions.csv")
        man_df = pd.DataFrame(rows_manifest)
        with open(man_path, 'w', newline='', encoding='utf-8') as _mfh:
            _mfh.write(f"# {_juice_header(self._stem, getattr(self, '_ode_backend', 'python'))}\n")
            man_df.to_csv(_mfh, index=False)
        print(f"  ... Reaction manifest saved to {man_path}")

        # ── KV: rate constants and turnover rates per reaction ────────────────────
        print("  ... writing kv (k-values + turnovers per reaction)...")
        kv_cols = {'TIME_s': self.sol.t}
        for i, lbl in enumerate(rxn_labels):
            r_str = rows_manifest[i]['reactants']
            p_str = rows_manifest[i]['products']
            hdr   = f"{lbl}:{r_str}-->{p_str}"[:60]
            kv_cols[f"k_{hdr}"] = k_arr[:, i]
            kv_cols[f"v_{hdr}"] = v_arr[:, i]
        kv_path = os.path.join(self.root, f"{self._stem}.kv")
        _write_csv(kv_path, pd.DataFrame(kv_cols).set_index('TIME_s'),
                   _juice_header(self._stem, getattr(self, '_ode_backend', 'python')), file_type='kv')
        print(f"  ... kv saved to {kv_path}")

        # ── PD: production and destruction rates per species ──────────────────────
        print("  ... computing pd (production + destruction per species)...")
        _species_list = list(self.solver.species)
        n_sp = len(_species_list)
        P_arr = np.zeros((n_t, n_sp), dtype=np.float64)
        D_arr = np.zeros((n_t, n_sp), dtype=np.float64)

        # kOH = Σ k_i × [non-OH ODE reactants] × [env reactants]  (s⁻¹)
        # Computed without [OH] in the product → valid even when [OH] = 0 (dark phase).
        _oh_idx = _species_list.index('OH') if 'OH' in _species_list else -1
        _kOH_arr = np.zeros(n_t, dtype=np.float64) if _oh_idx >= 0 else None

        for i, rxn in enumerate(self.mechanism.reactions):
            v_i = v_arr[:, i]
            for prod in getattr(rxn, 'products_parsed', []):
                P_arr[:, prod['idx']] += v_i * prod['coeff']
            for react in getattr(rxn, 'reactants_parsed', []):
                D_arr[:, react['idx']] += v_i * react['coeff']

            # kOH: only for reactions where OH is an ODE reactant
            if _kOH_arr is not None:
                _oh_reacts = [r for r in getattr(rxn, 'reactants_parsed', [])
                              if r['idx'] == _oh_idx]
                if _oh_reacts:
                    _ki = k_arr[:, i].copy()
                    _valid = True
                    for r in getattr(rxn, 'reactants_parsed', []):
                        if r['idx'] == _oh_idx:
                            continue  # skip OH itself
                        _ci = _sp_idx_map.get(_species_list[r['idx']])
                        if _ci is None:
                            _valid = False; break
                        _yc = y_matrix[:, _ci]
                        _ki *= _yc if r['coeff'] == 1.0 else _yc ** r['coeff']
                    if _valid:
                        for r in getattr(rxn, 'reactants_env', []):
                            _ec = _eff_env_arrs.get(r['name'])
                            if _ec is None:
                                if r['name'] == 'M':
                                    _ec = M_arr
                                else:
                                    _valid = False; break
                            _ecv = np.maximum(_ec, 0.0)
                            _ki *= _ecv if r['coeff'] == 1.0 else _ecv ** r['coeff']
                    if _valid:
                        _kOH_arr += _ki

        pd_cols = {'TIME_s': self.sol.t}
        for j, sp in enumerate(_species_list):
            pd_cols[f"P_{sp}"] = P_arr[:, j]
            pd_cols[f"D_{sp}"] = D_arr[:, j]
        if _kOH_arr is not None:
            pd_cols['kOH'] = _kOH_arr  # Σ k_i × [X_i] s⁻¹, ohne [OH]
        pd_path = os.path.join(self.root, f"{self._stem}.pd")
        _write_csv(pd_path, pd.DataFrame(pd_cols).set_index('TIME_s'),
                   _juice_header(self._stem, getattr(self, '_ode_backend', 'python')), file_type='pd')
        print(f"  ... pd saved to {pd_path}")

    def _show_results(self, out_path):
        # Launch viewer only when explicitly requested (--gui flag from FR_GUI).
        # Headless / command-line runs never open a GUI.
        if not getattr(self, '_show_viewer', False):
            return
        import subprocess
        viewer = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                              "results_viewer.py")
        _toml_args = ["--toml", self._toml_path] if getattr(self, '_toml_path', '') else []
        subprocess.Popen([sys.executable, viewer, out_path] + _toml_args)

if __name__ == "__main__":
    import argparse as _ap
    _parser = _ap.ArgumentParser(description="Run a JUICE scenario simulation.")
    _parser.add_argument("config", help="Path to scenario_config.json")
    _parser.add_argument("--gui", action="store_true",
                         help="Show the results viewer after the simulation")
    _parser.add_argument("--toml", default="",
                         help="Path to active TOML config (threads to viewer and solver)")
    _parser.add_argument("--solver", choices=["fast", "precise"], default=None,
                         help="solver profile, overriding the one written into the JSON: fast = rtol 1e-4 with "
                              "warm continuation (~2x faster); precise = rtol 1e-6 with a cold restart per ENZ row "
                              "(tightest agreement with FACSIMILE). See FR_GUI.toml [solver]")
    _args = _parser.parse_args()

    runner = ScenarioRunner(_args.config)
    if _args.solver:
        from scenario_prep import resolve_profile
        _sec = {}
        try:
            import tomllib as _tl
            with open(_args.toml or os.path.join(os.path.dirname(os.path.abspath(__file__)), "FR_GUI.toml"), "rb") as _fh:
                _sec = _tl.load(_fh).get("solver", {})
        except Exception:
            pass
        runner.cfg.update(resolve_profile(_sec, _args.solver))
    runner._show_viewer = _args.gui
    runner._toml_path   = _args.toml
    runner.run()