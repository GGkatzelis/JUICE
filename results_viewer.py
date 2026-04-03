import sys
import os
import re
import json
import warnings
import pandas as pd
import numpy as np
from datetime import datetime

# Set Windows AppUserModelID as early as possible so the taskbar uses our icon
if sys.platform == "win32":
    try:
        import ctypes
        ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID("FZJ.JUICE")
    except Exception:
        pass

# Force Python to look in the current script's folder
current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)

try:
    import netCDF4 as nc4
    HAS_NC4 = True
except ImportError:
    HAS_NC4 = False

try:
    from PyQt6.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
                                 QTreeWidget, QTreeWidgetItem, QListWidget, QListWidgetItem,
                                 QPushButton, QLabel, QComboBox, QSplitter, QLineEdit,
                                 QTabWidget, QFileDialog, QDialog, QDialogButtonBox,
                                 QTableWidget, QTableWidgetItem, QCheckBox, QHeaderView,
                                 QInputDialog, QMessageBox, QAbstractItemView)
    from PyQt6.QtCore import Qt
    from PyQt6.QtGui import QColor, QIcon, QPixmap
    from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
    from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
    from matplotlib.figure import Figure
    from matplotlib.ticker import FuncFormatter
    import matplotlib.patches as mpatches
    import matplotlib.cm as _mcm
    from matplotlib.colors import TwoSlopeNorm
    HAS_PYQT = True
except ImportError:
    HAS_PYQT = False

try:
    from easy_parser import EasyParser
    from mechanism import Mechanism, Reaction
    import mcm_kinetics
    HAS_BUDGET = True
except ImportError:
    HAS_BUDGET = False

# Seconds from Unix epoch (1970-01-01) to SAPHIR NC epoch (2000-01-01)
_NC_EPOCH_OFFSET = 946684800.0


def _fmt_hhmm(h, _pos):
    """Matplotlib tick formatter: decimal hours → HH:MM."""
    total_min = int(round(h * 60))
    hh = (total_min // 60) % 24
    mm = total_min % 60
    return f"{hh:02d}:{mm:02d}"

# ── colour palettes ───────────────────────────────────────────────────────────
_PROD_COLORS = ['#2ca02c', '#98df8a', '#17becf', '#9edae5', '#1f77b4',
                '#aec7e8', '#bcbd22', '#dbdb8d', '#c5b0d5', '#f7b6d2']
_DEST_COLORS = ['#d62728', '#ff9896', '#ff7f0e', '#ffbb78', '#9467bd',
                '#c49c94', '#8c564b', '#c7c7c7', '#e377c2', '#7f7f7f']

# ── Nitrogen budget helpers ───────────────────────────────────────────────────
_NB_CLASS_ORDER = [
    'NOx', 'HNO3', 'HONO', 'N2O5', 'PANs', 'Org. Nitrates', 'Other N',
]
_NB_CLASS_COLORS = {
    'NOx':          '#d62728',
    'HNO3':         '#ff7f0e',
    'HONO':         '#bcbd22',
    'N2O5':         '#2ca02c',
    'PANs':         '#1f77b4',
    'Org. Nitrates':'#9467bd',
    'Other N':      '#c7c7c7',
}

_MCM_NITROGEN_LOOKUP = {
    'N2O5': 2, 'HONO': 1, 'HNO3': 1, 'HNO4': 1,
    'NO': 1, 'NO2': 1, 'NO3': 1, 'N2O': 2,
    'PAN': 1, 'MPAN': 1, 'APAN': 1, 'GPAN': 1, 'PBZN': 1,
}

def _nitrogen_count(name):
    """Return the number of nitrogen atoms in a MCM species name."""
    upper = name.upper()
    if upper in _MCM_NITROGEN_LOOKUP:
        return _MCM_NITROGEN_LOOKUP[upper]
    total = 0
    for m in re.finditer(r'N(?!a|i|e)(\d*)', name):
        n = m.group(1)
        total += int(n) if n else 1
    return total

def _nitrogen_class(name):
    """Classify an MCM species into a chemical class for nitrogen budget grouping."""
    n = name.upper()
    if n in ('NO', 'NO2', 'NO3'):
        return 'NOx'
    if n in ('HNO3', 'HNO4'):
        return 'HNO3'
    if n == 'HONO':
        return 'HONO'
    if n == 'N2O5':
        return 'N2O5'
    if n.endswith('CO3NO2') or n in ('PAN', 'MPAN', 'APAN', 'GPAN', 'PBZN'):
        return 'PANs'
    if ('NO3' in n and n != 'NO3') or 'ONO2' in n or n.endswith('NO2'):
        return 'Org. Nitrates'
    return 'Other N'

# ── Carbon budget helpers ─────────────────────────────────────────────────────
# Lookup table for MCM species whose names are NOT condensed molecular formulas.
# Maps uppercase species name → integer carbon count.
_MCM_CARBON_LOOKUP = {
    # Isoprene and its common MCM fragments
    'ISOP': 5, 'ISOPD': 5, 'ISOPB': 5,
    # Monoterpenes (C10)
    'APINENE': 10, 'BPINENE': 10, 'LIMONENE': 10, 'MYRCENE': 10,
    'SABINENE': 10, 'D3CARENE': 10, 'TERPINOLENE': 10, 'CAMPHENE': 10,
    # Aromatics (no 'C' prefix in name → regex would return 0)
    'BENZENE': 6, 'TOLUENE': 7,
    'OXYLENE': 8, 'MXYLENE': 8, 'PXYLENE': 8, 'EBENZ': 8,
    # Common 4-carbon abbreviations
    'MEK': 4, 'MVK': 4, 'MACR': 4,
}

# Ordered list of classes for the "Grouped" stackplot (bottom → top).
# Dominant, slow-decaying groups at the bottom for visual stability.
_CB_CLASS_ORDER = [
    'Parent VOC', 'Aldehydes', 'Ketones', 'Carboxylic Acids',
    'Hydroperoxides', 'Alcohols', 'Peroxy Radicals', 'Acyl Peroxy',
    'Org. Nitrates', 'PANs', 'CO', 'Other',
]
_CB_CLASS_COLORS = {
    'Parent VOC':       '#1f77b4',
    'Aldehydes':        '#ff7f0e',
    'Ketones':          '#bcbd22',
    'Carboxylic Acids': '#2ca02c',
    'Hydroperoxides':   '#17becf',
    'Alcohols':         '#aec7e8',
    'Peroxy Radicals':  '#9467bd',
    'Acyl Peroxy':      '#c5b0d5',
    'Org. Nitrates':    '#ff9896',
    'PANs':             '#d62728',
    'CO':               '#8c564b',
    'Other':            '#c7c7c7',
}


def _carbon_count(name):
    """Return the number of carbon atoms in a MCM species name.

    Strategy:
    1. Check _MCM_CARBON_LOOKUP for abbreviated names (TOLUENE, ISOP, …).
    2. Apply a regex on the condensed-formula name: each occurrence of
       'C' (not followed by 'l' for chlorine) contributes 1 carbon;
       'Cn' (C followed by a digit) contributes n carbons.
    Returns 0 for inorganic/unrecognised species.
    """
    upper = name.upper()
    if upper in _MCM_CARBON_LOOKUP:
        return _MCM_CARBON_LOOKUP[upper]
    total = 0
    for m in re.finditer(r'C(?!l)(\d*)', name):
        n = m.group(1)
        total += int(n) if n else 1
    return total


def _carbon_class(name, parent_voc=None):
    """Classify an MCM species into a chemical class for carbon budget grouping.

    Priority order avoids mis-classification of complex species (e.g. acyl
    peroxy before generic peroxy, peracids before hydroperoxides).
    """
    n = name.upper()
    # Parent VOC (injected organic precursors from config)
    if parent_voc and n in parent_voc:
        return 'Parent VOC'
    # Simple C1 inorganic
    if n == 'CO':
        return 'CO'
    if n == 'CO2':
        return 'Other'          # rarely appears in MCM CSV
    # PAN family: peroxyacyl nitrates (e.g. CH3CO3NO2, MPAN)
    if n.endswith('CO3NO2') or n in ('PAN', 'MPAN', 'APAN', 'GPAN', 'PBZN'):
        return 'PANs'
    # Organic nitrates / organic peroxy nitrates
    if ('NO3' in n and n != 'NO3') or 'ONO2' in n:
        return 'Org. Nitrates'
    # Peracids (–C(=O)OOH written as CO3H in MCM)
    if n.endswith('CO3H'):
        return 'Hydroperoxides'
    # Acyl peroxy radicals (RC(=O)OO• → ends in CO3)
    if n.endswith('CO3'):
        return 'Acyl Peroxy'
    # Alkyl/other peroxy radicals (RO2)
    if (n.endswith('O2') and n != 'CO2') or n.endswith('OO'):
        return 'Peroxy Radicals'
    # Hydroperoxides (ROOH)
    if n.endswith('OOH'):
        return 'Hydroperoxides'
    # Carboxylic acids (RCOOH)
    if n.endswith('CO2H') or n.endswith('COOH') or n == 'HCOOH':
        return 'Carboxylic Acids'
    # Aldehydes (RCHO)
    if n.endswith('CHO') or n == 'HCHO':
        return 'Aldehydes'
    # Ketones / other C=O species (CO still in name after above filters)
    if 'CO' in n:
        return 'Ketones'
    # Alcohols (ROH, after OOH/CO2H/COOH already caught above)
    if n.endswith('OH'):
        return 'Alcohols'
    return 'Other'


# ─────────────────────────────────────────────────────────────────────────────
# BudgetEngine
# ─────────────────────────────────────────────────────────────────────────────
class BudgetEngine:
    """Vectorised reaction-rate budget engine for process analysis."""

    def __init__(self):
        self.cfg = None
        self.mechanism = None
        self.var_assignments = []
        self._ready = False

    # ── public ────────────────────────────────────────────────────────────────
    def load(self, config_path):
        if not HAS_BUDGET:
            raise ImportError("easy_parser / mechanism / mcm_kinetics not importable")

        with open(config_path) as f:
            self.cfg = json.load(f)

        chem_path = os.path.join(self.cfg['folder'], 'chemistry')
        parser = EasyParser(base_path=chem_path)
        self.mechanism = Mechanism()
        self.var_assignments = []

        mech_files = sorted(self.cfg['mechanism_files'], key=lambda x: '.rate' not in x)
        for f_name in mech_files:
            f_path = os.path.join(chem_path, f_name)
            if not os.path.exists(f_path):
                continue
            m_part = parser.parse(f_path)
            self.mechanism.species.extend(m_part.species)
            self.mechanism.reactions.extend(m_part.reactions)
            for p_name, p_expr in m_part.parameters.items():
                p_clean = p_expr.replace('@', '**')
                p_clean = re.sub(r'(?<=\d)[dD](?=[+-]?\d)', 'E', p_clean)
                p_clean = re.sub(r'J\((\d+)\)', r'J\1', p_clean)
                try:
                    self.var_assignments.append(
                        (p_name.upper(), compile(p_clean.strip(), '<string>', 'eval'))
                    )
                except Exception:
                    pass

        # Blocked parameters
        blocked = {p.upper() for p in self.cfg.get('blocked_parameters', [])}
        if blocked:
            self.var_assignments = [(n, c) for n, c in self.var_assignments if n not in blocked]

        # j_scalings
        for j_target, j_expr in self.cfg.get('j_scalings', {}).items():
            j_clean = j_expr.upper().replace('@', '**')
            try:
                self.var_assignments.append((j_target.upper(), compile(j_clean, '<string>', 'eval')))
            except Exception:
                pass

        # Artifact reactions (same as ScenarioRunner._add_artifacts)
        arts = self.cfg.get('artifacts', {})
        if arts.get('active'):
            scale = eval(str(arts.get('scale', '1.0')))
            if arts.get('hono_active'):
                self.mechanism.reactions.append(
                    Reaction([], ['HONO'], f"({arts['hono_eqn']}) * {scale}"))
            if arts.get('hcho_active'):
                self.mechanism.reactions.append(
                    Reaction([], ['HCHO'], arts['hcho_eqn']))
            if arts.get('y_active'):
                self.mechanism.reactions.append(
                    Reaction(['OH'], ['HO2'], f"({arts['y_eqn']}) * 1e-9 * M * 2.4e-13"))
            if arts.get('wall_active'):
                self.mechanism.reactions.append(
                    Reaction([], ['HONO'], str(arts['dark_hono'])))

        # Custom reactions
        for r in self.cfg.get('custom_reactions', []):
            if not r.get('reactants') or not r.get('products') or not r.get('rate'):
                continue
            reacs = [x.strip().upper() for x in r['reactants'].split('+') if x.strip()]
            prods = [x.strip().upper() for x in r['products'].split('+') if x.strip()]
            self.mechanism.reactions.append(Reaction(reacs, prods, r['rate']))

        # Excluded reactions
        excl = self.cfg.get('excluded_reactions', [])
        if excl:
            def _canon(lst):
                s = set()
                for x in lst:
                    name = x.split('*')[-1].strip().upper()
                    if name and not name.replace('.', '').isdigit():
                        s.add(name)
                return frozenset(s)
            excl_keys = set()
            for e in excl:
                if e.get('action', 'Remove') != 'Modify Rate':
                    excl_keys.add((
                        _canon([x.strip() for x in e['reactants'].split('+') if x.strip()]),
                        _canon([x.strip() for x in e['products'].split('+') if x.strip()])
                    ))
            self.mechanism.reactions = [
                rxn for rxn in self.mechanism.reactions
                if (_canon(rxn.reactants), _canon(rxn.products)) not in excl_keys
            ]

        # Compile safe_expr for every reaction
        for rxn in self.mechanism.reactions:
            expr = getattr(rxn, 'rate_expression', '0.0').upper()
            expr = expr.replace('@', '**')
            expr = re.sub(r'(?<=\d)[dD](?=[+-]?\d)', 'E', expr)
            expr = re.sub(r'J\((\d+)\)', r'J\1', expr)
            try:
                rxn.safe_expr = compile(expr, '<string>', 'eval')
            except Exception:
                rxn.safe_expr = compile('0.0', '<string>', 'eval')

        self._ready = True
        print(f"[BudgetEngine] Loaded {len(self.mechanism.reactions)} reactions, "
              f"{len(self.var_assignments)} parameters from {os.path.basename(config_path)}")

    def _prepare_ctx(self, df):
        """Build fully-evaluated vectorised context from a results DataFrame."""
        n = len(df)
        ctx = self._build_vec_ctx(df)

        protected = set()
        for env_key in self.cfg.get('env_map', {}):
            if f'SAPHIR_{env_key}'.upper() in ctx:
                protected.add(env_key.upper())

        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            for _ in range(2):
                for var_name, comp_expr in self.var_assignments:
                    if var_name in protected:
                        continue
                    try:
                        result = eval(comp_expr, globals(), ctx)
                        ctx[var_name] = (np.full(n, float(result))
                                         if np.isscalar(result)
                                         else np.asarray(result, dtype=float))
                    except Exception:
                        ctx.setdefault(var_name, np.zeros(n))

        if 'SAPHIR_JO1D' in ctx:
            ctx['J1'] = ctx['SAPHIR_JO1D']
        # J4, J7, J3 are always overridden by measured values (JSCALE approx in anorg phot.rate
        # sets them non-zero, so the ==0 guard would silently ignore measured data).
        for j_name, j_col in [('J4', 'SAPHIR_JNO2'), ('J7', 'SAPHIR_JHONO'), ('J3', 'SAPHIR_JH2O2')]:
            if j_col in ctx:
                ctx[j_name] = ctx[j_col]
        for j_name, j_col in [('J11', 'SAPHIR_JHCHO_R'), ('J12', 'SAPHIR_JHCHO_M')]:
            if j_col in ctx and np.all(ctx.get(j_name, np.zeros(1)) == 0):
                ctx[j_name] = ctx[j_col]

        T_arr   = ctx.get('TEMP', np.full(n, 298.15))
        k_N2_q  = 2.15e-11 * np.exp(110.0 / T_arr)
        k_O2_q  = 3.30e-11 * np.exp( 55.0 / T_arr)
        k_H2O_q = np.full(n, 2.20e-10)
        k_quench = (k_N2_q  * ctx.get('N2',  0.7809 * ctx['M']) +
                    k_O2_q  * ctx.get('O2',  0.2095 * ctx['M']) +
                    k_H2O_q * ctx.get('H2O', np.zeros(n)))
        hv_arr = ctx.get('HV', np.ones(n))
        j1 = ctx.get('J1', np.zeros(n)) * np.where(hv_arr > 0.5, 1.0, 0.0)
        o3 = ctx.get('O3', np.zeros(n))
        ctx['O1D'] = np.where(k_quench > 0, j1 * o3 / np.maximum(k_quench, 1e-30), 0.0)
        return ctx

    def compute(self, species, df):
        """Return (prod_dict, dest_dict) {label: np.array} fluxes in ppb hr⁻¹."""
        if not self._ready:
            return {}, {}

        spec_up = species.upper()
        n = len(df)
        ctx = self._prepare_ctx(df)
        M_arr = ctx['M']

        prod, dest = {}, {}

        for rxn in self.mechanism.reactions:
            reac_parsed = self._parse_stoich(rxn.reactants)
            prod_parsed = self._parse_stoich(rxn.products)

            reac_names = [name for _, name in reac_parsed]
            prod_names = [name for _, name in prod_parsed]
            is_prod = spec_up in prod_names
            is_dest = spec_up in reac_names
            if not is_prod and not is_dest:
                continue

            try:
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore')
                    k = eval(rxn.safe_expr, globals(), ctx)
            except Exception:
                continue

            if callable(k):
                continue
            k = np.full(n, float(k)) if np.isscalar(k) else np.asarray(k, dtype=float)

            v = np.maximum(k, 0.0).copy()
            for r_coeff, r_name in reac_parsed:
                conc = ctx.get(r_name)
                if conc is None:
                    v = np.zeros(n)
                    break
                c = np.full(n, float(conc)) if np.isscalar(conc) else np.asarray(conc, dtype=float)
                v *= np.maximum(c, 0.0) ** r_coeff

            v_ppb_hr = np.maximum(v / M_arr, 0.0) * 1e9 * 3600.0
            label = self._make_label(rxn)

            if is_prod:
                for coeff, name in prod_parsed:
                    if name == spec_up:
                        prod[label] = prod.get(label, np.zeros(n)) + coeff * v_ppb_hr

            if is_dest:
                for coeff, name in reac_parsed:
                    if name == spec_up:
                        dest[label] = dest.get(label, np.zeros(n)) + coeff * v_ppb_hr

        return prod, dest

    def compute_reactivity(self, radical, df):
        """Return {label: np.array} pseudo-first-order loss rate coefficients (s⁻¹)."""
        if not self._ready:
            return {}

        rad_up = radical.upper()
        n = len(df)
        ctx = self._prepare_ctx(df)

        contributions = {}
        for rxn in self.mechanism.reactions:
            reac_parsed = self._parse_stoich(rxn.reactants)
            if rad_up not in [name for _, name in reac_parsed]:
                continue

            try:
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore')
                    k = eval(rxn.safe_expr, globals(), ctx)
            except Exception:
                continue
            if callable(k):
                continue
            k = np.full(n, float(k)) if np.isscalar(k) else np.asarray(k, dtype=float)

            v = np.maximum(k, 0.0).copy()
            radical_skipped = False
            for r_coeff, r_name in reac_parsed:
                if r_name == rad_up and not radical_skipped:
                    radical_skipped = True
                    if r_coeff > 1.0:
                        conc = ctx.get(r_name)
                        if conc is None:
                            v = np.zeros(n); break
                        c = (np.full(n, float(conc)) if np.isscalar(conc)
                             else np.asarray(conc, dtype=float))
                        v *= np.maximum(c, 0.0) ** (r_coeff - 1.0)
                else:
                    conc = ctx.get(r_name)
                    if conc is None:
                        v = np.zeros(n); break
                    c = (np.full(n, float(conc)) if np.isscalar(conc)
                         else np.asarray(conc, dtype=float))
                    v *= np.maximum(c, 0.0) ** r_coeff

            label = self._make_label(rxn)
            contributions[label] = contributions.get(label, np.zeros(n)) + v

        return contributions

    # ── private helpers ───────────────────────────────────────────────────────
    @staticmethod
    def _parse_stoich(species_list):
        out = []
        for s in species_list:
            parts = s.strip().upper().split('*', 1)
            if len(parts) == 2:
                try:
                    out.append((float(parts[0]), parts[1]))
                    continue
                except ValueError:
                    pass
            out.append((1.0, parts[0]))
        return out

    @staticmethod
    def _make_label(rxn):
        reacs = [r.split('*')[-1].strip() for r in rxn.reactants] or ['hv/wall']
        prods = [p.split('*')[-1].strip() for p in rxn.products]
        lhs = ' + '.join(reacs[:2]) + (' + ...' if len(reacs) > 2 else '')
        rhs = ' + '.join(prods[:3]) + (' + ...' if len(prods) > 3 else '')
        return f"{lhs} → {rhs}"

    def _build_vec_ctx(self, df):
        n = len(df)
        ctx = {}
        for col in df.columns:
            try:
                ctx[col.upper()] = df[col].values.astype(float)
            except Exception:
                pass

        M_arr = ctx.get('M', np.full(n, 2.5e19))

        T_arr = ctx.get('SAPHIR_TEMP', ctx.get('SAPHIR_T', np.full(n, 298.15)))
        HV_arr = ctx.get('SAPHIR_HV', np.ones(n))
        JNO2_arr = ctx.get('SAPHIR_JNO2', np.zeros(n))

        COSX_arr = np.full(n, 0.5)
        if 'SAPHIR_SZ' in ctx:
            COSX_arr = np.maximum(np.cos(ctx['SAPHIR_SZ'] * np.pi / 180.0), 1e-4)
        safe_cosx = np.where(COSX_arr > 0, COSX_arr, 1e-4)

        ctx.update({
            'T': T_arr, 'TEMP': T_arr,
            'M': M_arr,
            'O2': 0.2095 * M_arr,
            'N2': 0.7809 * M_arr,
            'HV': HV_arr,
            'COSX': COSX_arr,
            'JNO2': JNO2_arr,
            'np': np, 'exp': np.exp, 'EXP': np.exp,
            'LOG10': np.log10, 'log10': np.log10,
            'LOG': np.log, 'log': np.log,
            'COS': np.cos, 'SIN': np.sin, 'SQRT': np.sqrt,
            'AMAX': np.maximum, 'AMIN': np.minimum,
        })

        for func_name in dir(mcm_kinetics):
            if func_name.startswith('_'):
                continue
            func = getattr(mcm_kinetics, func_name)
            if func_name.startswith('KMT'):
                ctx[func_name.upper()] = lambda *args, f=func, T=T_arr, M=M_arr: f(*args, T, M)
            else:
                ctx[func_name.upper()] = func

        def _phtfkt(a, b, c, _c=safe_cosx, _j=JNO2_arr, _hv=HV_arr):
            cs_j = a * _c ** b * np.exp(np.clip(c / _c, -700, 700))
            jno2_cs = 1.165e-2 * _c ** 0.244 * np.exp(np.clip(-0.267 / _c, -700, 700))
            scale = np.where(
                (_j > 0) & (jno2_cs > 0),
                np.minimum(_j / np.where(jno2_cs > 0, jno2_cs, 1e-30), 1.0),
                1.0
            )
            return np.where(_hv > 0.5, cs_j * scale, 0.0)

        ctx['PHTFKT'] = _phtfkt
        ctx['PHTFKT_CS'] = lambda a, b, c, _c=safe_cosx: (
            a * _c ** b * np.exp(np.clip(c / _c, -700, 700))
        )

        h2o = ctx.get('SAPHIR_H2O', 0.01 * M_arr)
        h2o_max = np.nanmax(h2o)
        if h2o_max < 1e14:
            if h2o_max > 1.0:
                h2o = (h2o / 100.0) * M_arr
            else:
                h2o = h2o * M_arr
        ctx['H2O'] = h2o

        p_sat = np.exp(21.36469 - 5339.66 / T_arr)
        ctx['RH'] = np.clip(
            (ctx['H2O'] / M_arr) * 100.0 / (p_sat / (1.38e-16 * T_arr * M_arr)), 0, 100
        )

        for env_key in self.cfg.get('env_map', {}):
            if env_key.upper() == 'H2O':
                continue
            saphir_col = f'SAPHIR_{env_key}'.upper()
            if saphir_col in ctx:
                ctx[env_key.upper()] = ctx[saphir_col]

        return ctx


# ─────────────────────────────────────────────────────────────────────────────
# _NcFileSelector  —  checkable list of NC files
# ─────────────────────────────────────────────────────────────────────────────
class _NcFileSelector(QDialog):
    """Dialog showing NC files in a directory with checkboxes for multi-selection."""

    def __init__(self, nc_dir, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Select NC Files to Load")
        self.resize(750, 520)
        layout = QVBoxLayout(self)

        layout.addWidget(QLabel(f"<b>Folder:</b> {nc_dir}"))

        self._list = QListWidget()
        layout.addWidget(self._list)

        try:
            nc_files = sorted(f for f in os.listdir(nc_dir) if f.lower().endswith('.nc'))
        except OSError:
            nc_files = []

        self._full_paths = []
        for fname in nc_files:
            item = QListWidgetItem(fname)
            item.setFlags(item.flags() | Qt.ItemFlag.ItemIsUserCheckable)
            item.setCheckState(Qt.CheckState.Unchecked)
            self._list.addItem(item)
            self._full_paths.append(os.path.join(nc_dir, fname))

        btn_row = QHBoxLayout()
        btn_all = QPushButton("Check All")
        btn_all.clicked.connect(self._check_all)
        btn_none = QPushButton("Uncheck All")
        btn_none.clicked.connect(self._uncheck_all)
        btn_row.addWidget(btn_all)
        btn_row.addWidget(btn_none)
        btn_row.addStretch()
        layout.addLayout(btn_row)

        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def _check_all(self):
        for i in range(self._list.count()):
            self._list.item(i).setCheckState(Qt.CheckState.Checked)

    def _uncheck_all(self):
        for i in range(self._list.count()):
            self._list.item(i).setCheckState(Qt.CheckState.Unchecked)

    def selected_paths(self):
        result = []
        for i in range(self._list.count()):
            if self._list.item(i).checkState() == Qt.CheckState.Checked:
                result.append(self._full_paths[i])
        return result


# ─────────────────────────────────────────────────────────────────────────────
# _MatchDialog  —  verify NC variable → MCM mapping
# ─────────────────────────────────────────────────────────────────────────────
class _MatchDialog(QDialog):
    """Shows loaded NC variables with suggested MCM matches for user verification."""

    def __init__(self, var_rows, mcm_columns, parent=None):
        """
        var_rows   : list of dicts {nc_file, nc_var, units, suggested_mcm, t, v}
        mcm_columns: list of MCM column names for the dropdown
        """
        super().__init__(parent)
        self.setWindowTitle("Verify NC Variable Mapping")
        self.resize(960, 540)
        layout = QVBoxLayout(self)

        info = QLabel(
            "Review the NC variables found. The <b>MCM Match</b> column is informational — "
            "it links an instrument variable to the closest model species for reference. "
            "Check <b>Load</b> for each variable you want to add to the browser tree."
        )
        info.setWordWrap(True)
        layout.addWidget(info)

        self._table = QTableWidget(len(var_rows), 5)
        self._table.setHorizontalHeaderLabels(
            ["NC File", "NC Variable", "Units", "MCM Match", "Load"])
        hh = self._table.horizontalHeader()
        hh.setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch)
        hh.setSectionResizeMode(1, QHeaderView.ResizeMode.ResizeToContents)
        hh.setSectionResizeMode(2, QHeaderView.ResizeMode.ResizeToContents)
        hh.setSectionResizeMode(3, QHeaderView.ResizeMode.Stretch)
        hh.setSectionResizeMode(4, QHeaderView.ResizeMode.ResizeToContents)

        self._var_rows = var_rows
        self._combos = []
        self._checks = []

        mcm_opts = [""] + sorted(mcm_columns)

        for row_idx, row in enumerate(var_rows):
            def _ro(text):
                it = QTableWidgetItem(text)
                it.setFlags(it.flags() & ~Qt.ItemFlag.ItemIsEditable)
                return it

            self._table.setItem(row_idx, 0, _ro(row['nc_file']))
            self._table.setItem(row_idx, 1, _ro(row['nc_var']))
            self._table.setItem(row_idx, 2, _ro(row.get('units', '')))

            combo = QComboBox()
            combo.addItems(mcm_opts)
            if row['suggested_mcm'] in mcm_opts:
                combo.setCurrentText(row['suggested_mcm'])
            self._table.setCellWidget(row_idx, 3, combo)
            self._combos.append(combo)

            chk = QCheckBox()
            chk.setChecked(True)
            wrapper = QWidget()
            wl = QHBoxLayout(wrapper)
            wl.addWidget(chk)
            wl.setAlignment(Qt.AlignmentFlag.AlignCenter)
            wl.setContentsMargins(0, 0, 0, 0)
            self._table.setCellWidget(row_idx, 4, wrapper)
            self._checks.append(chk)

        layout.addWidget(self._table)

        btn_row = QHBoxLayout()
        btn_all = QPushButton("Check All")
        btn_all.clicked.connect(lambda: [c.setChecked(True) for c in self._checks])
        btn_none = QPushButton("Uncheck All")
        btn_none.clicked.connect(lambda: [c.setChecked(False) for c in self._checks])
        btn_row.addWidget(btn_all)
        btn_row.addWidget(btn_none)
        btn_row.addStretch()
        layout.addLayout(btn_row)

        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def confirmed_rows(self):
        """Return list of (var_row_dict, mcm_match_str) for all checked rows."""
        return [
            (row, combo.currentText())
            for row, combo, chk in zip(self._var_rows, self._combos, self._checks)
            if chk.isChecked()
        ]


# ─────────────────────────────────────────────────────────────────────────────
# ResultsViewer
# ─────────────────────────────────────────────────────────────────────────────
class ResultsViewer(QMainWindow):
    _PALETTE = [
        "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
        "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
        "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
        "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5",
    ]

    def __init__(self, df, csv_path=None):
        super().__init__()
        self.df = df.copy()
        self._color_index = 0
        self.active_lines = {}
        self.csv_path = csv_path

        # The CSV index is in SAPHIR epoch (seconds since 2000-01-01), matching the ENZ
        # file header "TIME / SEC SINCE 1.1.2000 00:00:00".  NC archive files use the
        # same epoch.  Convert to true Unix epoch (seconds since 1970-01-01) once here
        # so that _t0_unix / _midnight_unix are consistent with datetime.utcfromtimestamp
        # and the NC time formula  t_hours = (t_nc + _NC_EPOCH_OFFSET - _midnight_unix) / 3600.
        self._t0_unix = float(self.df.index[0]) + _NC_EPOCH_OFFSET  # SAPHIR → Unix
        self._midnight_unix = self._t0_unix - (self._t0_unix % 86400)
        self._last_nc_dir = None

        # Instrument data loaded from NC files: display_key -> {t, v, nc_file, nc_var, mcm_match}
        self._inst_data = {}

        self.setWindowTitle("SAPHIR/MCM Results Viewer")
        self.setGeometry(50, 50, 1600, 900)
        _icon_pixmap = QPixmap(os.path.join(os.path.dirname(os.path.abspath(__file__)), "Juelich Logo.png"))
        if not _icon_pixmap.isNull():
            self.setWindowIcon(QIcon(_icon_pixmap))

        # Convert MCM index from SAPHIR epoch (sec since 2000-01-01) to hours since
        # UTC midnight of the experiment day, same reference as NC instrument data.
        self.df.index = ((self.df.index.astype(float) + _NC_EPOCH_OFFSET) - self._midnight_unix) / 3600.0

        self.budget_engine = BudgetEngine() if HAS_BUDGET else None
        self._config_label_text = "No config loaded"

        # ── Central layout ──────────────────────────────────────────────────
        central = QWidget()
        self.setCentralWidget(central)
        main_layout = QHBoxLayout(central)
        splitter = QSplitter(Qt.Orientation.Horizontal)
        main_layout.addWidget(splitter)

        # ── LEFT PANEL ──────────────────────────────────────────────────────
        left_panel = QWidget()
        left_layout = QVBoxLayout(left_panel)
        left_layout.setContentsMargins(4, 4, 4, 4)

        self.search_box = QLineEdit()
        self.search_box.setPlaceholderText("Search variable...")
        self.search_box.textChanged.connect(self.filter_list)
        left_layout.addWidget(self.search_box)

        # Tree widget replacing the old flat list
        self.var_tree = QTreeWidget()
        self.var_tree.setHeaderHidden(True)
        self.var_tree.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        self.var_tree.itemDoubleClicked.connect(self._on_tree_double_click)
        left_layout.addWidget(self.var_tree)

        # Top-level folder nodes
        self._mcm_folder    = QTreeWidgetItem(self.var_tree, ["MCM Output"])
        self._saphir_folder = QTreeWidgetItem(self.var_tree, ["SAPHIR Environment"])
        self._inst_folder   = QTreeWidgetItem(self.var_tree, ["SAPHIR Instruments"])

        self._mcm_folder.setExpanded(True)
        self._saphir_folder.setExpanded(True)
        self._inst_folder.setExpanded(True)

        self._populate_tree()

        # NC data loader button
        btn_nc = QPushButton("Load SAPHIR NC Data...")
        btn_nc.clicked.connect(self._load_nc_data)
        left_layout.addWidget(btn_nc)

        # Config loader
        cfg_row = QHBoxLayout()
        btn_cfg = QPushButton("Load Scenario Config (.json)")
        btn_cfg.clicked.connect(self._browse_config)
        cfg_row.addWidget(btn_cfg)
        left_layout.addLayout(cfg_row)

        self.cfg_status = QLabel(self._config_label_text)
        self.cfg_status.setWordWrap(True)
        self.cfg_status.setStyleSheet("color: #555; font-size: 10px;")
        left_layout.addWidget(self.cfg_status)

        # Timeseries controls
        control_layout = QVBoxLayout()
        control_layout.addWidget(QLabel("Target Axis:"))
        self.axis_combo = QComboBox()
        self.axis_combo.addItems(["Left Axis (Main)", "Right Axis (Main)",
                                  "Left Axis (Offset)", "Right Axis (Offset)"])
        control_layout.addWidget(self.axis_combo)

        control_layout.addWidget(QLabel("Conversion (MCM Data):"))
        self.unit_combo = QComboBox()
        self.unit_combo.addItems(["Convert to ppb", "Native (molec/cm³ or SAPHIR unit)"])
        control_layout.addWidget(self.unit_combo)

        btn_add = QPushButton("Plot Selected")
        btn_add.clicked.connect(self.add_selected)
        control_layout.addWidget(btn_add)

        control_layout.addWidget(QLabel("Active Plots (Double-Click to Remove):"))
        self.active_list = QListWidget()
        self.active_list.setSelectionMode(QListWidget.SelectionMode.ExtendedSelection)
        self.active_list.itemDoubleClicked.connect(self.remove_selected)
        control_layout.addWidget(self.active_list)

        btn_remove = QPushButton("Remove Selected from Plot")
        btn_remove.clicked.connect(lambda: self.remove_selected(None))
        control_layout.addWidget(btn_remove)

        btn_clear = QPushButton("Clear Canvas")
        btn_clear.clicked.connect(self.clear_plot)
        control_layout.addWidget(btn_clear)

        left_layout.addLayout(control_layout)
        splitter.addWidget(left_panel)

        # ── RIGHT PANEL: tabbed ─────────────────────────────────────────────
        right_panel = QWidget()
        right_layout = QVBoxLayout(right_panel)
        right_layout.setContentsMargins(0, 0, 0, 0)

        self.right_tabs = QTabWidget()
        right_layout.addWidget(self.right_tabs)

        # Tab 0 – Timeseries
        ts_widget = QWidget()
        ts_layout = QVBoxLayout(ts_widget)
        ts_layout.setContentsMargins(0, 0, 0, 0)
        self.fig = Figure(figsize=(10, 6))
        self.canvas = FigureCanvas(self.fig)
        self.toolbar = NavigationToolbar(self.canvas, self)
        ts_layout.addWidget(self.toolbar)
        ts_layout.addWidget(self.canvas)
        self.right_tabs.addTab(ts_widget, "Timeseries")

        # Tab 1 – Process Analysis
        pa_widget = QWidget()
        pa_layout = QVBoxLayout(pa_widget)
        pa_layout.setContentsMargins(0, 0, 0, 0)
        self.pa_fig = Figure(figsize=(12, 6))
        self.pa_canvas = FigureCanvas(self.pa_fig)
        self.pa_toolbar = NavigationToolbar(self.pa_canvas, self)
        pa_layout.addWidget(self.pa_toolbar)
        pa_layout.addWidget(self.pa_canvas)
        self.right_tabs.addTab(pa_widget, "Process Analysis")

        # Tab 2 – Reactivity
        rx_widget = QWidget()
        rx_layout = QVBoxLayout(rx_widget)
        rx_layout.setContentsMargins(2, 2, 2, 0)
        rx_layout.setSpacing(2)
        rx_ctrl = QHBoxLayout()
        rx_ctrl.setContentsMargins(0, 0, 0, 0)
        rx_ctrl.addWidget(QLabel("Radical:"))
        self.rx_combo = QComboBox()
        self.rx_combo.addItems(["OH", "O3", "NO3"])
        rx_ctrl.addWidget(self.rx_combo)
        rx_btn = QPushButton("Compute Reactivity")
        rx_btn.clicked.connect(self._update_reactivity)
        rx_ctrl.addWidget(rx_btn)
        rx_meas_btn = QPushButton("Load k_OH (NC file)")
        rx_meas_btn.clicked.connect(self._load_koh_nc)
        rx_ctrl.addWidget(rx_meas_btn)
        self.rx_meas_label = QLabel("No measured k_OH loaded")
        self.rx_meas_label.setStyleSheet("color: #555; font-size: 10px;")
        rx_ctrl.addWidget(self.rx_meas_label)
        rx_ctrl.addStretch()
        rx_layout.addLayout(rx_ctrl, 0)
        self._measured_koh = None
        self.rx_fig = Figure(figsize=(12, 6))
        self.rx_canvas = FigureCanvas(self.rx_fig)
        self.rx_toolbar = NavigationToolbar(self.rx_canvas, self)
        rx_layout.addWidget(self.rx_toolbar, 0)
        rx_layout.addWidget(self.rx_canvas, 1)
        self.right_tabs.addTab(rx_widget, "Reactivity")

        # Tab 3 – Carbon Budget
        self._build_carbon_tab()

        # Tab 4 – Comparison (last)
        self._build_comparison_tab()

        splitter.addWidget(right_panel)
        splitter.setSizes([300, 1300])

        self.init_axes()
        self._plot_defaults()

        # Auto-discover scenario_config.json — must come after UI is fully built
        self._auto_load_config(csv_path)

    # ── Tree building ─────────────────────────────────────────────────────────
    def _populate_tree(self):
        """Fill MCM Output and SAPHIR Environment folders from self.df columns."""
        self._mcm_folder.takeChildren()
        self._saphir_folder.takeChildren()

        mcm_active, mcm_empty = [], []
        for col in sorted(self.df.columns):
            if col == 'M':
                continue
            if col.upper().startswith('SAPHIR_'):
                item = QTreeWidgetItem([col])
                self._saphir_folder.addChild(item)
            else:
                vals = self.df[col].values
                has_data = np.nanmax(np.abs(vals)) > 0 if len(vals) > 0 else False
                if has_data:
                    mcm_active.append(col)
                else:
                    mcm_empty.append(col)

        for col in mcm_active:
            self._mcm_folder.addChild(QTreeWidgetItem([col]))

        for col in mcm_empty:
            item = QTreeWidgetItem([col])
            item.setForeground(0, QColor('#aaaaaa'))
            item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEnabled)
            self._mcm_folder.addChild(item)

    def _get_or_create_inst_subfolder(self, nc_basename):
        """Find or create a sub-folder under SAPHIR Instruments for nc_basename."""
        for i in range(self._inst_folder.childCount()):
            child = self._inst_folder.child(i)
            if child.text(0) == nc_basename:
                return child
        subfolder = QTreeWidgetItem(self._inst_folder, [nc_basename])
        subfolder.setExpanded(True)
        return subfolder

    # ── Tree interaction ──────────────────────────────────────────────────────
    def _on_tree_single_click(self, item, _column):
        """Single click on a leaf: run budget for MCM/SAPHIR, show info for instruments."""
        if item.childCount() != 0:
            return  # folder node — ignore

        display_key = item.data(0, Qt.ItemDataRole.UserRole)
        self.right_tabs.setCurrentIndex(1)

        if display_key and display_key in self._inst_data:
            data = self._inst_data[display_key]
            self.pa_fig.clear()
            ax = self.pa_fig.add_subplot(111)
            ax.text(0.5, 0.55,
                    f"Instrument variable\n\n"
                    f"File:  {data['nc_file']}\n"
                    f"Variable:  {data['nc_var']}\n"
                    f"MCM match:  {data['mcm_match'] or '—'}\n\n"
                    f"No mechanism budget for instrument data.",
                    ha='center', va='center', transform=ax.transAxes,
                    fontsize=11, family='monospace')
            ax.axis('off')
            self.pa_canvas.draw()
        else:
            self._update_budget(item.text(0))

    def _on_tree_double_click(self, item, _column):
        """Double-click: context-aware action based on the active right tab.
        - Timeseries tab  → plot the species on the timeseries canvas.
        - Process Analysis tab → run the budget for this species.
        - Other tabs → no-op.
        """
        if item.childCount() != 0:
            return  # folder node
        active = self.right_tabs.currentIndex()
        if active == 0:          # Timeseries
            self.add_selected()
        elif active == 1:        # Process Analysis
            self._update_budget(item.text(0))

    # ── Search / filter ───────────────────────────────────────────────────────
    def filter_list(self, text):
        text_low = text.lower()

        for folder in [self._mcm_folder, self._saphir_folder]:
            any_visible = False
            for i in range(folder.childCount()):
                child = folder.child(i)
                enabled = bool(child.flags() & Qt.ItemFlag.ItemIsEnabled)
                visible = enabled and ((not text_low) or (text_low in child.text(0).lower()))
                child.setHidden(not visible)
                if visible:
                    any_visible = True
            folder.setHidden(bool(text_low) and not any_visible)

        inst_any = False
        for i in range(self._inst_folder.childCount()):
            subfolder = self._inst_folder.child(i)
            sf_any = False
            for j in range(subfolder.childCount()):
                child = subfolder.child(j)
                visible = (not text_low) or (text_low in child.text(0).lower())
                child.setHidden(not visible)
                if visible:
                    sf_any = True
            subfolder.setHidden(bool(text_low) and not sf_any)
            if sf_any:
                inst_any = True
        self._inst_folder.setHidden(bool(text_low) and not inst_any)

    # ── Add / remove / clear ─────────────────────────────────────────────────
    def add_selected(self, *_):
        # Collect all selected leaves (items with a parent and no children)
        selected = [
            item for item in self.var_tree.selectedItems()
            if item.parent() is not None and item.childCount() == 0
        ]

        ax_name = self.axis_combo.currentText()
        unit = self.unit_combo.currentText()

        for item in selected:
            display_key = item.data(0, Qt.ItemDataRole.UserRole)
            if display_key and display_key in self._inst_data:
                self._plot_instrument(display_key, ax_name)
            else:
                self.plot_variable(item.text(0), ax_name, unit)

        self._update_legends()
        self._update_ylabels()
        self.canvas.draw()
        self.var_tree.clearSelection()
        self.right_tabs.setCurrentIndex(0)

    def remove_selected(self, item=None):
        items = [item] if item is not None else self.active_list.selectedItems()
        for i in items:
            list_label = i.text()
            if list_label in self.active_lines:
                line, _ = self.active_lines[list_label]
                try:
                    line.remove()
                except ValueError:
                    pass
                del self.active_lines[list_label]
            self.active_list.takeItem(self.active_list.row(i))

        self._rescale_axes()
        self._update_legends()
        self._update_ylabels()
        self.canvas.draw()

    def _rescale_axes(self):
        for ax in self.axes_map.values():
            ax.relim()
            ax.autoscale_view()

    def clear_plot(self):
        self.active_list.clear()
        self.active_lines.clear()
        self.init_axes()
        self.canvas.draw()

    # ── MCM variable plotting ─────────────────────────────────────────────────
    def plot_variable(self, var_name, axis_name, unit):
        ax = self.axes_map[axis_name]
        data = self.df[var_name].copy()

        if "ppb" in unit and not var_name.upper().startswith("SAPHIR_") and 'M' in self.df.columns:
            data = (data / self.df['M']) * 1e9
            label = f"{var_name} [ppb]"
        else:
            if var_name.upper().startswith("SAPHIR_"):
                label = f"{var_name}"
            else:
                label = f"{var_name} [molec/cm³]"

        list_label = f"{label} [{axis_name}]"

        for i in range(self.active_list.count()):
            if self.active_list.item(i).text() == list_label:
                return

        color = self._PALETTE[self._color_index % len(self._PALETTE)]
        self._color_index += 1
        line, = ax.plot(self.df.index, data, label=label, linewidth=2.0, color=color)

        self.active_lines[list_label] = (line, ax)
        item = QListWidgetItem(list_label)
        item.setForeground(QColor(color))
        self.active_list.addItem(item)

        # Auto-overlay matching SAPHIR observation as bullets (same colour).
        # _plot_instrument handles active_list insertion internally.
        for display_key, inst in self._inst_data.items():
            if inst.get('mcm_match', '').upper() == var_name.upper():
                self._plot_instrument(display_key, axis_name, color=color)
                break

    # ── Instrument data plotting ──────────────────────────────────────────────
    def _plot_instrument(self, display_key, axis_name, color=None):
        """Plot a variable from self._inst_data as filled circles (observations).

        When *color* is supplied the instrument shares the model's colour so
        the pair is visually linked.  A new palette colour is allocated only
        when plotting an instrument that has no model counterpart on the canvas.
        """
        data = self._inst_data[display_key]
        nc_var  = data['nc_var']
        nc_file = data['nc_file']

        label = f"{nc_var} (obs)"
        list_label = f"{label} [{axis_name}]"

        for i in range(self.active_list.count()):
            if self.active_list.item(i).text() == list_label:
                return

        ax = self.axes_map[axis_name]
        if color is None:
            color = self._PALETTE[self._color_index % len(self._PALETTE)]
            self._color_index += 1

        t = data['t']
        v = data['v']

        mask = np.isfinite(v)
        line, = ax.plot(t[mask], v[mask], label=label, color=color,
                        marker='o', markersize=4, alpha=0.65,
                        linestyle='none', linewidth=0)

        self.active_lines[list_label] = (line, ax)
        item = QListWidgetItem(list_label)
        item.setForeground(QColor(color))
        self.active_list.addItem(item)

    # ── NC data loading workflow ──────────────────────────────────────────────
    def _load_nc_data(self):
        """Full workflow: infer date → confirm path → select files → match dialog → add to tree."""
        if not HAS_NC4:
            QMessageBox.warning(self, "netCDF4 Unavailable",
                                "The netCDF4 package is not installed.\n"
                                "Install it with:  pip install netCDF4")
            return

        # Infer date from t0_unix
        try:
            t0_dt = datetime.utcfromtimestamp(self._t0_unix)
            inferred_path = (
                f"Z:/IEK8-SAPHIR/data/"
                f"{t0_dt.year:04d}/{t0_dt.month:02d}/{t0_dt.day:02d}"
            ).replace('/', os.sep)
        except Exception:
            inferred_path = r"Z:\IEK8-SAPHIR\data"

        # Ask user to confirm / edit the path
        nc_dir, ok = QInputDialog.getText(
            self, "SAPHIR Archive Path",
            "Path to NC files (edit if the inferred date is wrong):",
            text=inferred_path
        )
        if not ok or not nc_dir.strip():
            return

        nc_dir = nc_dir.strip()
        if not os.path.isdir(nc_dir):
            QMessageBox.warning(self, "Path Not Found",
                                f"Directory not found:\n{nc_dir}")
            return

        self._last_nc_dir = nc_dir

        # Let user pick files from the directory
        selector = _NcFileSelector(nc_dir, parent=self)
        if selector.exec() != QDialog.DialogCode.Accepted:
            return

        paths = selector.selected_paths()
        if not paths:
            QMessageBox.information(self, "No Files Selected",
                                    "No NC files were checked. Nothing to load.")
            return

        # Load variables from selected files
        var_rows = self._load_nc_vars(paths)
        if not var_rows:
            QMessageBox.warning(self, "No Variables Found",
                                "Could not extract any variables from the selected files.\n"
                                "Check the console for details.")
            return

        # Show match/verification dialog
        mcm_cols = [c for c in self.df.columns if c != 'M']
        dlg = _MatchDialog(var_rows, mcm_cols, parent=self)
        if dlg.exec() != QDialog.DialogCode.Accepted:
            return

        confirmed = dlg.confirmed_rows()
        if not confirmed:
            return

        self._add_confirmed_rows(confirmed)
        print(f"[NC Loader] Added {len(confirmed)} instrument variable(s) to the browser.")

    def _load_nc_vars(self, paths):
        """Open each NC file and return a list of variable descriptor dicts."""
        var_rows = []
        for path in paths:
            fname = os.path.basename(path)
            nc_basename = fname[:-3] if fname.lower().endswith('.nc') else fname
            try:
                ds = nc4.Dataset(path)
            except Exception as e:
                print(f"[NC Loader] Cannot open {fname}: {e}")
                continue

            if 'time' not in ds.variables:
                print(f"[NC Loader] No 'time' variable in {fname} — skipping")
                ds.close()
                continue

            try:
                t_nc = np.ma.filled(ds.variables['time'][:], np.nan).astype(float)
                t_hours = (t_nc + _NC_EPOCH_OFFSET - self._midnight_unix) / 3600.0
            except Exception as e:
                print(f"[NC Loader] Cannot read time from {fname}: {e}")
                ds.close()
                continue

            for var in ds.variables:
                if var.lower().startswith('time') or '@' in var:
                    continue
                try:
                    raw = np.ma.filled(ds.variables[var][:], np.nan)
                    if raw.ndim > 1:
                        raw = raw.squeeze()
                    if raw.ndim != 1:
                        print(f"[NC Loader] Skipping {var} in {fname}: not 1-D after squeeze")
                        continue
                    if len(raw) != len(t_hours):
                        print(f"[NC Loader] Length mismatch for {var} in {fname} "
                              f"({len(raw)} vs {len(t_hours)}) — skipping")
                        continue
                    values = raw.astype(float)
                    units = getattr(ds.variables[var], 'units', '')
                    var_rows.append({
                        'nc_file':       nc_basename,
                        'nc_var':        var,
                        'units':         units,
                        'suggested_mcm': self._find_mcm_match(var),
                        't':             t_hours,
                        'v':             values,
                    })
                except Exception as e:
                    print(f"[NC Loader] Cannot read {var} from {fname}: {e}")

            ds.close()

        return var_rows

    def _find_mcm_match(self, nc_var):
        """Case-insensitive match: exact, then SAPHIR_-prefixed."""
        nc_up = nc_var.upper()
        mcm_cols = [c for c in self.df.columns if c != 'M']
        for c in mcm_cols:
            if c.upper() == nc_up:
                return c
        for c in mcm_cols:
            if c.upper() == f'SAPHIR_{nc_up}':
                return c
        return ""

    def _add_confirmed_rows(self, confirmed):
        """Store confirmed variables in self._inst_data and add them to the tree."""
        for row, mcm_match in confirmed:
            nc_basename = row['nc_file']
            nc_var      = row['nc_var']
            display_key = f"{nc_basename}/{nc_var}"

            if display_key in self._inst_data:
                continue  # already loaded

            self._inst_data[display_key] = {
                't':         row['t'],
                'v':         row['v'],
                'nc_file':   nc_basename,
                'nc_var':    nc_var,
                'mcm_match': mcm_match,
            }

            subfolder = self._get_or_create_inst_subfolder(nc_basename)
            child = QTreeWidgetItem(subfolder, [nc_var])
            child.setData(0, Qt.ItemDataRole.UserRole, display_key)
            if mcm_match:
                child.setToolTip(0, f"MCM match: {mcm_match}")

        self._inst_folder.setExpanded(True)

    # ── config helpers ────────────────────────────────────────────────────────
    def _auto_load_config(self, csv_path):
        if not csv_path or not HAS_BUDGET:
            return
        candidates = []
        d = os.path.dirname(os.path.abspath(csv_path))
        for root in [d, os.path.dirname(d)]:
            for fname in os.listdir(root):
                if fname == 'scenario_config.json':
                    candidates.append(os.path.join(root, fname))
        if candidates:
            self._load_config(candidates[0], silent=True)

    def _browse_config(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Open Scenario Config", "", "JSON Files (*.json)"
        )
        if path:
            self._load_config(path)

    def _load_config(self, path, silent=False):
        if not HAS_BUDGET:
            self.cfg_status.setText("Budget engine unavailable (missing easy_parser)")
            return
        try:
            self.budget_engine.load(path)
            name = os.path.basename(path)
            self._config_label_text = f"✓ {name}"
            self.cfg_status.setText(self._config_label_text)
            if not silent:
                print(f"[BudgetEngine] Config loaded: {path}")
        except Exception as e:
            self.cfg_status.setText(f"Error: {e}")
            print(f"[BudgetEngine] Load failed: {e}")

    # ── Process Analysis tab ──────────────────────────────────────────────────
    def _update_budget(self, species):
        self.pa_fig.clear()
        ax_msg = self.pa_fig.add_subplot(111)

        if not HAS_BUDGET:
            ax_msg.text(0.5, 0.5, "Budget engine unavailable\n(install easy_parser)",
                        ha='center', va='center', transform=ax_msg.transAxes, fontsize=12)
            self.pa_canvas.draw()
            return

        if not (self.budget_engine and self.budget_engine._ready):
            ax_msg.text(0.5, 0.5,
                        "No scenario config loaded.\n"
                        "Use 'Load Scenario Config (.json)' in the left panel.",
                        ha='center', va='center', transform=ax_msg.transAxes, fontsize=12)
            self.pa_canvas.draw()
            return

        try:
            prod, dest = self.budget_engine.compute(species, self.df)
        except Exception as e:
            ax_msg.text(0.5, 0.5, f"Budget compute failed:\n{e}",
                        ha='center', va='center', transform=ax_msg.transAxes, fontsize=10)
            self.pa_canvas.draw()
            return

        if not prod and not dest:
            ax_msg.text(0.5, 0.5,
                        f"'{species}' not found in mechanism reactions.\n"
                        f"(It may be an env/SAPHIR variable, not a modelled species.)",
                        ha='center', va='center', transform=ax_msg.transAxes, fontsize=11)
            self.pa_canvas.draw()
            return

        self._draw_budget(species, prod, dest)

    def _draw_budget(self, species, prod, dest):
        """Draw stacked production/destruction timeseries + mean bar chart."""
        self.pa_fig.clear()
        t = self.df.index.values

        THRESH = 0.02

        def _group(pathways):
            if not pathways:
                return {}, np.zeros(len(t))
            total = max(sum(np.nanmean(v) for v in pathways.values()), 1e-30)
            major, other = {}, None
            for lbl, arr in sorted(pathways.items(), key=lambda x: -np.nanmean(x[1])):
                if np.nanmean(arr) / total >= THRESH:
                    major[lbl] = arr
                else:
                    other = arr if other is None else other + arr
            return major, (other if other is not None else np.zeros(len(t)))

        major_prod, other_prod = _group(prod)
        major_dest, other_dest = _group(dest)

        if np.nanmean(other_prod) > 0:
            major_prod['Other production'] = other_prod
        if np.nanmean(other_dest) > 0:
            major_dest['Other destruction'] = other_dest

        gs = self.pa_fig.add_gridspec(1, 2, width_ratios=[3, 1], wspace=0.35,
                                      left=0.07, right=0.93, top=0.88, bottom=0.08)
        ax_ts  = self.pa_fig.add_subplot(gs[0, 0])
        ax_bar = self.pa_fig.add_subplot(gs[0, 1])

        # HV shading
        hv_col = next((c for c in self.df.columns if c.upper() in ('SAPHIR_HV', 'HV')), None)
        if hv_col is not None:
            hv = self.df[hv_col].values
            in_open, t0 = False, None
            for ti, vi in zip(t, hv):
                if vi > 0.5 and not in_open:
                    in_open, t0 = True, ti
                elif vi <= 0.5 and in_open:
                    ax_ts.axvspan(t0, ti, color='#fffacd', alpha=0.4, zorder=0, linewidth=0)
                    in_open = False
            if in_open:
                ax_ts.axvspan(t0, t[-1], color='#fffacd', alpha=0.4, zorder=0, linewidth=0)

        prod_labels = list(major_prod.keys())
        prod_arrays = [major_prod[k] for k in prod_labels]
        prod_colors = _PROD_COLORS[:len(prod_arrays)]
        if prod_arrays:
            ax_ts.stackplot(t, prod_arrays, labels=prod_labels, colors=prod_colors, alpha=0.80)

        dest_labels = list(major_dest.keys())
        dest_arrays = [-major_dest[k] for k in dest_labels]
        dest_colors = _DEST_COLORS[:len(dest_arrays)]
        if dest_arrays:
            ax_ts.stackplot(t, dest_arrays, labels=dest_labels, colors=dest_colors, alpha=0.80)

        net = (sum(prod.values()) if prod else np.zeros(len(t))) \
            - (sum(dest.values()) if dest else np.zeros(len(t)))
        ax_ts.plot(t, net, color='black', linewidth=1.2, linestyle='--',
                   label='Net dC/dt', zorder=5)
        ax_ts.axhline(0, color='black', linewidth=0.6, alpha=0.4)

        # Shade the destruction region (below zero) with a light red tint.
        # Save and restore y-limits so the span's infinite extent doesn't
        # override the data-driven autoscale.
        _ylim = ax_ts.get_ylim()
        ax_ts.axhspan(_ylim[0], 0, color='#ffdddd', alpha=0.25, zorder=0, linewidth=0)
        ax_ts.set_ylim(_ylim)

        ax_ts.set_xlabel("Time of Day (UTC)", fontweight='bold')
        ax_ts.xaxis.set_major_formatter(FuncFormatter(_fmt_hhmm))
        ax_ts.set_ylabel("Flux (ppb hr⁻¹)", fontweight='bold')
        ax_ts.set_title(f"Process Analysis: {species}", fontweight='bold')
        ax_ts.grid(True, linestyle='--', alpha=0.35)

        handles, labels = ax_ts.get_legend_handles_labels()
        if handles:
            ax_ts.legend(handles, labels, loc='upper left',
                         fontsize=7, framealpha=0.75, ncol=1,
                         handlelength=1.2, handleheight=0.8,
                         borderpad=0.5, labelspacing=0.3)

        # ── Bar chart: all positive, destruction at bottom, production at top ──
        # Destruction: sorted ascending → smallest at bottom, largest just below middle
        dest_sorted = sorted(major_dest.items(), key=lambda x: np.nanmean(x[1]))
        # Production: sorted descending → largest just above middle, smallest at top
        prod_sorted = sorted(major_prod.items(), key=lambda x: -np.nanmean(x[1]))

        bar_labels, bar_vals, bar_colors = [], [], []
        n_dest = len(dest_sorted)

        for lbl, arr in dest_sorted:
            bar_labels.append(lbl)
            bar_vals.append(np.nanmean(arr))          # positive (absolute)
            idx = next((i for i, k in enumerate(dest_labels) if k == lbl), 0)
            bar_colors.append(_DEST_COLORS[idx % len(_DEST_COLORS)])

        for lbl, arr in prod_sorted:
            bar_labels.append(lbl)
            bar_vals.append(np.nanmean(arr))          # positive
            idx = next((i for i, k in enumerate(prod_labels) if k == lbl), 0)
            bar_colors.append(_PROD_COLORS[idx % len(_PROD_COLORS)])

        y_pos = np.arange(len(bar_labels))
        if len(bar_labels) > 0:
            # Shade destruction rows (bottom block) with a light red background
            if n_dest > 0:
                ax_bar.axhspan(-0.5, n_dest - 0.5, color='#ffe8e8', alpha=0.6,
                               zorder=0, linewidth=0)
            # Shade production rows (top block) with a light green background
            if len(prod_sorted) > 0:
                ax_bar.axhspan(n_dest - 0.5, len(bar_labels) - 0.5,
                               color='#e8f5e8', alpha=0.6, zorder=0, linewidth=0)

            ax_bar.barh(y_pos, bar_vals, color=bar_colors, alpha=0.85,
                        edgecolor='white', linewidth=0.5)
            ax_bar.set_yticks(y_pos)
            ax_bar.set_yticklabels(bar_labels, fontsize=7)

            # Divider line between destruction and production blocks
            if n_dest > 0 and len(prod_sorted) > 0:
                ax_bar.axhline(n_dest - 0.5, color='#888888', linewidth=1.2,
                               linestyle='-', zorder=3)

            ax_bar.set_xlabel("|Mean flux|\n(ppb hr⁻¹)", fontsize=8)
            ax_bar.set_title("Time-average", fontsize=9)
            ax_bar.grid(True, axis='x', linestyle='--', alpha=0.35)

            # Section labels outside the axes box, centred on each block.
            # y in axes-fraction: each bar occupies 1 data unit, axis spans
            # from -0.5 to len(bar_labels)-0.5, so fraction = (y+0.5)/n_total
            n_total = len(bar_labels)
            n_prod  = len(prod_sorted)
            if n_dest > 0:
                dest_frac = n_dest / (2.0 * n_total)
                ax_bar.annotate("DEST", xy=(1.04, dest_frac),
                                xycoords='axes fraction',
                                fontsize=7.5, color='#c0392b', fontweight='bold',
                                ha='center', va='center', rotation=90,
                                clip_on=False)
            if n_prod > 0:
                prod_frac = (n_dest + n_prod / 2.0) / n_total
                ax_bar.annotate("PROD", xy=(1.04, prod_frac),
                                xycoords='axes fraction',
                                fontsize=7.5, color='#27ae60', fontweight='bold',
                                ha='center', va='center', rotation=90,
                                clip_on=False)

        self.pa_canvas.draw()

    # ── Timeseries axis helpers ───────────────────────────────────────────────
    def init_axes(self):
        self._color_index = 0
        self.fig.clear()
        self.ax1 = self.fig.add_subplot(111)
        self.ax1.set_xlabel("Time of Day (UTC)", fontweight='bold')
        self.ax1.xaxis.set_major_formatter(FuncFormatter(_fmt_hhmm))
        self.ax1.grid(True, linestyle='--', alpha=0.4)
        self.ax2 = self.ax1.twinx()
        self.ax3 = self.ax1.twinx()
        self.ax3.spines['right'].set_position(('outward', 80))
        self.ax4 = self.ax1.twinx()
        self.ax4.spines['left'].set_position(('outward', 80))
        self.ax4.yaxis.set_ticks_position('left')
        self.ax4.yaxis.set_label_position('left')

        self.axes_map = {
            "Left Axis (Main)":   self.ax1,
            "Right Axis (Main)":  self.ax2,
            "Right Axis (Offset)": self.ax3,
            "Left Axis (Offset)": self.ax4,
        }
        self.fig.subplots_adjust(left=0.22, right=0.78, top=0.85)
        self._update_ylabels()

    def _update_ylabels(self):
        spine_map = {
            "Left Axis (Main)":    ('left',  0),
            "Right Axis (Main)":   ('right', 0),
            "Right Axis (Offset)": ('right', 80),
            "Left Axis (Offset)":  ('left',  80),
        }
        for ax_name, ax in self.axes_map.items():
            lines = ax.get_lines()
            spine_side, _ = spine_map[ax_name]

            if hasattr(ax, 'color_dummies'):
                for dummy in ax.color_dummies:
                    try:
                        dummy.remove()
                    except Exception:
                        pass
                delattr(ax, 'color_dummies')

            if not lines:
                ax.set_ylabel("")
                ax.tick_params(axis='y', colors='black')
                ax.spines[spine_side].set_color('black')
                if ax != self.ax1:
                    ax.yaxis.set_visible(False)
                    ax.spines[spine_side].set_visible(False)
                continue

            ax.yaxis.set_visible(True)
            ax.spines[spine_side].set_visible(True)

            labels = [line.get_label() for line in lines]
            colors = [line.get_color() for line in lines]

            if len(lines) == 1:
                ax.set_ylabel(labels[0], color=colors[0], fontweight='bold')
                ax.tick_params(axis='y', colors=colors[0])
                ax.spines[spine_side].set_color(colors[0])
            else:
                ylabel_text = ", ".join(labels)
                ax.set_ylabel(ylabel_text, color='black', fontweight='bold')
                ax.tick_params(axis='y', colors='black')
                ax.spines[spine_side].set_color('black')

    def _update_legends(self):
        lines, labels = [], []
        for ax in self.axes_map.values():
            h, l = ax.get_legend_handles_labels()
            lines.extend(h)
            labels.extend(l)
        if self.fig.legends:
            self.fig.legends[0].remove()
        if lines:
            leg = self.fig.legend(lines, labels, loc='upper center',
                                  bbox_to_anchor=(0.5, 0.98), ncol=6)
            for text, line in zip(leg.get_texts(), lines):
                text.set_color(line.get_color())
                text.set_fontweight('bold')

    # ── Reactivity tab ────────────────────────────────────────────────────────
    def _load_koh_nc(self):
        """Load measured k_OH from a SAPHIR NC file."""
        if not HAS_NC4:
            self.rx_meas_label.setText("netCDF4 not available")
            self.rx_meas_label.setStyleSheet("color: #cc0000; font-size: 10px;")
            return

        start_dir = self._last_nc_dir or r"Z:\IEK8-SAPHIR\data"
        path, _ = QFileDialog.getOpenFileName(
            self, "Load k_OH from NC File", start_dir,
            "NetCDF Files (*.nc);;All Files (*)"
        )
        if not path:
            return

        try:
            ds = nc4.Dataset(path)
            if 'time' not in ds.variables:
                raise ValueError("No 'time' variable found")

            t_nc = np.ma.filled(ds.variables['time'][:], np.nan).astype(float)
            t_hours = (t_nc + _NC_EPOCH_OFFSET - self._midnight_unix) / 3600.0

            # Prefer a variable whose name contains kOH / k_OH; fall back to first primary var
            koh_var = None
            for v in ds.variables:
                if v.lower().startswith('time') or '@' in v:
                    continue
                if re.search(r'(?i)k_?oh', v):
                    koh_var = v
                    break
            if koh_var is None:
                for v in ds.variables:
                    if not v.lower().startswith('time') and '@' not in v:
                        koh_var = v
                        break
            if koh_var is None:
                raise ValueError("No suitable variable found")

            k_vals = np.ma.filled(ds.variables[koh_var][:], np.nan).squeeze().astype(float)
            ds.close()

            self._measured_koh = {'t': t_hours, 'k': k_vals}
            fname = os.path.basename(path)
            self.rx_meas_label.setText(f"k_OH: {fname} / {koh_var}")
            self.rx_meas_label.setStyleSheet(
                "color: #1a7a1a; font-size: 10px; font-weight: bold;")
            print(f"[Reactivity] k_OH loaded from NC: {len(t_hours)} points, "
                  f"range {np.nanmin(k_vals):.2f}–{np.nanmax(k_vals):.2f} s⁻¹")
        except Exception as e:
            self.rx_meas_label.setText(f"Load failed: {e}")
            self.rx_meas_label.setStyleSheet("color: #cc0000; font-size: 10px;")

    def _auto_measured_koh_from_df(self):
        """Return (t_hours, k_s⁻¹) if k_OH is available in self.df or self._inst_data."""
        # 1. Check self.df (SAPHIR_ env columns from CSV)
        candidates = [c for c in self.df.columns
                      if re.search(r'(?i)^(saphir_)?k_?oh$', c)]
        if candidates:
            col = candidates[0]
            k = self.df[col].values.astype(float)
            t = self.df.index.values.astype(float)
            mask = np.isfinite(k) & (k > 0)
            if mask.sum() >= 5:
                return t, k

        # 2. Check instrument data loaded from NC files
        for key, data in self._inst_data.items():
            if re.search(r'(?i)k_?oh', key.split('/')[-1]):
                t = data['t']
                k = data['v']
                mask = np.isfinite(k) & (k > 0)
                if mask.sum() >= 5:
                    return t, k

        return None

    def _update_reactivity(self):
        radical = self.rx_combo.currentText()
        self.rx_fig.clear()
        ax_msg = self.rx_fig.add_subplot(111)

        if not (self.budget_engine and self.budget_engine._ready):
            ax_msg.text(0.5, 0.5,
                        "No scenario config loaded.\nUse 'Load Scenario Config (.json)'.",
                        ha='center', va='center', transform=ax_msg.transAxes, fontsize=12)
            self.rx_canvas.draw()
            return

        try:
            contribs = self.budget_engine.compute_reactivity(radical, self.df)
        except Exception as e:
            ax_msg.text(0.5, 0.5, f"Reactivity compute failed:\n{e}",
                        ha='center', va='center', transform=ax_msg.transAxes, fontsize=10)
            self.rx_canvas.draw()
            return

        if not contribs:
            ax_msg.text(0.5, 0.5,
                        f"No reactions found consuming {radical}.\n"
                        f"Check that {radical} is present in the loaded mechanism.",
                        ha='center', va='center', transform=ax_msg.transAxes, fontsize=11)
            self.rx_canvas.draw()
            return

        self._draw_reactivity(radical, contribs)

    def _draw_reactivity(self, radical, contribs):
        """Draw stacked reactivity contributions (s⁻¹) + total + mean bar chart."""
        self.rx_fig.clear()
        t = self.df.index.values
        n = len(t)
        THRESH = 0.02

        total_mean = max(sum(np.nanmean(v) for v in contribs.values()), 1e-30)
        major, other = {}, None
        for lbl, arr in sorted(contribs.items(), key=lambda x: -np.nanmean(x[1])):
            if np.nanmean(arr) / total_mean >= THRESH:
                major[lbl] = arr
            else:
                other = arr if other is None else other + arr
        if other is not None and np.nanmean(other) > 0:
            major['Other'] = other

        gs = self.rx_fig.add_gridspec(1, 2, width_ratios=[3, 1], wspace=0.35,
                                      left=0.07, right=0.97, top=0.88, bottom=0.12)
        ax_ts  = self.rx_fig.add_subplot(gs[0, 0])
        ax_bar = self.rx_fig.add_subplot(gs[0, 1])

        hv_col = next((c for c in self.df.columns if c.upper() in ('SAPHIR_HV', 'HV')), None)
        if hv_col is not None:
            hv = self.df[hv_col].values
            in_open, t0 = False, None
            for ti, vi in zip(t, hv):
                if vi > 0.5 and not in_open:
                    in_open, t0 = True, ti
                elif vi <= 0.5 and in_open:
                    ax_ts.axvspan(t0, ti, color='#fffacd', alpha=0.4, zorder=0, linewidth=0)
                    in_open = False
            if in_open:
                ax_ts.axvspan(t0, t[-1], color='#fffacd', alpha=0.4, zorder=0, linewidth=0)

        labels  = list(major.keys())
        arrays  = [major[k] for k in labels]
        colors  = _PROD_COLORS[:len(arrays)]
        if arrays:
            ax_ts.stackplot(t, arrays, labels=labels, colors=colors, alpha=0.80)

        k_total = sum(contribs.values()) if contribs else np.zeros(n)
        ax_ts.plot(t, k_total, color='black', lw=1.5, ls='--',
                   label=f'Total k_{radical} (model)', zorder=5)

        # Measured k_OH overlay
        if radical.upper() == 'OH':
            meas = self._measured_koh
            if meas is None:
                meas_from_df = self._auto_measured_koh_from_df()
                if meas_from_df is not None:
                    meas = {'t': meas_from_df[0], 'k': meas_from_df[1]}
            if meas is not None:
                t_m = meas['t']
                k_m = meas['k']
                mask = np.isfinite(k_m)
                ax_ts.plot(t_m[mask], k_m[mask], color='#e31a1c', lw=2.0,
                           ls='-', marker='o', markersize=2, alpha=0.85,
                           label='Measured k_OH (SAPHIR)', zorder=6)

        ax_ts.set_xlabel("Time of Day (UTC)", fontweight='bold')
        ax_ts.xaxis.set_major_formatter(FuncFormatter(_fmt_hhmm))
        ax_ts.set_ylabel(f"k_{radical} (s⁻¹)", fontweight='bold')
        ax_ts.set_title(f"Total {radical} Reactivity", fontweight='bold')
        ax_ts.grid(True, linestyle='--', alpha=0.35)
        handles, lbls = ax_ts.get_legend_handles_labels()
        if handles:
            ax_ts.legend(handles, lbls, loc='upper center',
                         bbox_to_anchor=(0.5, -0.14), ncol=3, fontsize=7, framealpha=0.7)

        bar_labels = list(major.keys())
        bar_vals   = [np.nanmean(major[k]) for k in bar_labels]
        bar_colors = colors[:len(bar_labels)]
        y_pos = np.arange(len(bar_labels))
        if bar_labels:
            ax_bar.barh(y_pos, bar_vals, color=bar_colors, alpha=0.85,
                        edgecolor='white', linewidth=0.5)
            ax_bar.set_yticks(y_pos)
            ax_bar.set_yticklabels(bar_labels, fontsize=7)
            ax_bar.axvline(0, color='black', linewidth=0.8)
            ax_bar.set_xlabel("Mean k\n(s⁻¹)", fontsize=8)
            ax_bar.set_title("Time-average", fontsize=9)
            ax_bar.grid(True, axis='x', linestyle='--', alpha=0.35)

        self.rx_canvas.draw()

    # ─────────────────────────────────────────────────────────────────────────
    # Comparison tab  (Tab 3)
    # ─────────────────────────────────────────────────────────────────────────
    _NOISE_EXCEPTIONS = {'OH', 'HO2', 'NO3', 'O1D', 'CH3O', 'C2H5O', 'HCO', 'CHO'}

    def _build_comparison_tab(self):
        """Build the Comparison tab: reference loader + slope matrix + two timeseries panels."""
        # State
        self._cmp_df_ref        = None
        self._cmp_ref_label     = ""
        self._cmp_slope_results = {}
        self._cmp_slope_box_map = []
        self._cmp_selected      = [None, None]   # species selected per slot
        self._cmp_ts_data       = [None, None]   # (py, ref, common_t, sp) per slot
        self._cmp_crosshair     = [
            {'vline': None, 'text': None},
            {'vline': None, 'text': None},
        ]
        self._cmp_mirror_vlines = [None, None]   # vline drawn in canvas s by the OTHER slot's hover
        self._cmp_mirror_texts  = [None, None]   # text drawn in canvas s by the OTHER slot's hover

        outer = QWidget()
        outer_layout = QVBoxLayout(outer)
        outer_layout.setContentsMargins(4, 4, 4, 4)

        # ── Reference source bar ─────────────────────────────────────────────
        src_row = QHBoxLayout()
        src_row.addWidget(QLabel("Reference:"))
        self._cmp_src_combo = QComboBox()
        self._cmp_src_combo.addItems([
            "FACSIMILE output (.cs)",
            "Previous simulation (.csv)",
            "Loaded SAPHIR instruments (NC)",
        ])
        src_row.addWidget(self._cmp_src_combo)
        load_btn = QPushButton("Load Reference…")
        load_btn.clicked.connect(self._cmp_load_reference)
        src_row.addWidget(load_btn)

        src_row.addSpacing(12)
        src_row.addWidget(QLabel("Noise threshold:"))
        self._cmp_noise_input = QLineEdit("1e-6")
        self._cmp_noise_input.setMaximumWidth(70)
        src_row.addWidget(self._cmp_noise_input)
        btn_slopes = QPushButton("Compute Slopes")
        btn_slopes.clicked.connect(self._cmp_compute_slopes)
        src_row.addWidget(btn_slopes)

        self._cmp_ref_lbl = QLabel("No reference loaded")
        self._cmp_ref_lbl.setStyleSheet("color: #666; font-style: italic;")
        src_row.addWidget(self._cmp_ref_lbl, stretch=1)
        outer_layout.addLayout(src_row)

        # ── Splitter: 1/3 matrix | 2/3 detail ────────────────────────────────
        splitter = QSplitter(Qt.Orientation.Horizontal)

        left = QWidget()
        ll = QVBoxLayout(left)
        ll.setContentsMargins(0, 0, 0, 0)
        self._cmp_slope_fig = Figure()
        self._cmp_slope_canvas = FigureCanvas(self._cmp_slope_fig)
        self._cmp_slope_canvas.mpl_connect('button_press_event', self._cmp_on_slope_click)
        ll.addWidget(NavigationToolbar(self._cmp_slope_canvas, self))
        ll.addWidget(self._cmp_slope_canvas)
        splitter.addWidget(left)

        # ── Right panel: two stacked timeseries slots ─────────────────────────
        right = QWidget()
        rl = QVBoxLayout(right)
        rl.setContentsMargins(0, 0, 0, 0)

        right_vsplit = QSplitter(Qt.Orientation.Vertical)

        self._cmp_ts_figs     = []
        self._cmp_ts_canvases = []
        self._cmp_ts_combos   = []

        for slot in range(2):
            slot_widget = QWidget()
            slot_layout = QVBoxLayout(slot_widget)
            slot_layout.setContentsMargins(2, 2, 2, 2)
            slot_layout.setSpacing(2)

            ctrl = QHBoxLayout()
            lbl_txt = "Left-click a box to select:" if slot == 0 else "Right-click a box to select:"
            ctrl.addWidget(QLabel(lbl_txt))
            combo = QComboBox()
            combo.setMinimumWidth(120)
            s = slot
            combo.currentTextChanged.connect(
                lambda sp, s=s: self._cmp_draw_timeseries(s, sp) if sp else None)
            self._cmp_ts_combos.append(combo)
            ctrl.addWidget(combo, stretch=1)
            slot_layout.addLayout(ctrl)

            fig = Figure()
            canvas = FigureCanvas(fig)
            self._cmp_ts_figs.append(fig)
            self._cmp_ts_canvases.append(canvas)
            slot_layout.addWidget(NavigationToolbar(canvas, self))
            slot_layout.addWidget(canvas, stretch=1)

            canvas.mpl_connect('motion_notify_event',
                               lambda ev, s=s: self._cmp_on_mouse_move(ev, s))
            canvas.mpl_connect('axes_leave_event',
                               lambda ev, s=s: self._cmp_clear_crosshair(s))

            right_vsplit.addWidget(slot_widget)

        rl.addWidget(right_vsplit, stretch=1)
        splitter.addWidget(right)

        splitter.setSizes([500, 1000])
        outer_layout.addWidget(splitter, 1)

        self.right_tabs.addTab(outer, "Comparison")

    # ── Reference loading ─────────────────────────────────────────────────────
    def _cmp_load_reference(self):
        src = self._cmp_src_combo.currentText()

        if src.startswith("FACSIMILE"):
            path, _ = QFileDialog.getOpenFileName(
                self, "Load FACSIMILE Output", "", "CS Files (*.cs);;All Files (*)")
            if not path:
                return
            try:
                df = self._cmp_parse_facsimile(path)
                self._cmp_df_ref    = df
                self._cmp_ref_label = os.path.basename(path)
            except Exception as e:
                QMessageBox.warning(self, "Load Error", str(e))
                return

        elif src.startswith("Previous simulation"):
            path, _ = QFileDialog.getOpenFileName(
                self, "Load Previous Simulation CSV", "", "CSV Files (*.csv)")
            if not path:
                return
            try:
                df = pd.read_csv(path, index_col=0)
                # Same SAPHIR-epoch → hours-since-midnight conversion as main CSV
                t_unix = df.index.astype(float) + _NC_EPOCH_OFFSET
                df.index = (t_unix - self._midnight_unix) / 3600.0
                self._cmp_df_ref    = df
                self._cmp_ref_label = os.path.basename(path)
            except Exception as e:
                QMessageBox.warning(self, "Load Error", str(e))
                return

        else:  # SAPHIR Instruments (already loaded NC data)
            if not self._inst_data:
                QMessageBox.information(self, "No Instruments Loaded",
                                        "Load SAPHIR NC data first (left panel button).")
                return
            # Build a DataFrame from _inst_data on a shared time grid
            t_min = min(d['t'].min() for d in self._inst_data.values())
            t_max = max(d['t'].max() for d in self._inst_data.values())
            common_t = np.linspace(t_min, t_max, 500)
            rows = {'__t__': common_t}
            for key, data in self._inst_data.items():
                col = data.get('nc_var', key)
                rows[col] = np.interp(common_t, data['t'], data['v'],
                                      left=np.nan, right=np.nan)
            df = pd.DataFrame(rows).set_index('__t__')
            df.index.name = 'Time_hours'
            self._cmp_df_ref    = df
            self._cmp_ref_label = "Loaded SAPHIR instruments"

        self._cmp_ref_lbl.setText(f"Ref: {self._cmp_ref_label}")
        self._cmp_ref_lbl.setStyleSheet("color: black; font-style: normal;")

    def _cmp_parse_facsimile(self, path):
        """Parse a FACSIMILE .cs file; align to hours-since-midnight of this experiment."""
        col_names = {}
        data_lines = []
        with open(path, 'r') as f:
            for line in f:
                s = line.strip()
                if s.startswith('COLUMN'):
                    m = re.match(r'COLUMN\s+(\d+)=(.+)', s)
                    if m:
                        idx  = int(m.group(1))
                        name = m.group(2).split('<')[0].strip()
                        col_names[idx] = 'TIME' if idx == 1 else name
                elif s and s[0].isdigit():
                    # Fix missing-E scientific notation: 3.84-116 → 3.84E-116
                    fixed = re.sub(r'(\d)([-+])(\d)', r'\1E\2\3', s)
                    data_lines.append(fixed)

        if not col_names or not data_lines:
            raise ValueError("No column definitions or data rows found in FACSIMILE file.")

        n_cols  = max(col_names.keys())
        headers = [col_names.get(i, f'COL{i}') for i in range(1, n_cols + 1)]
        rows    = [list(map(float, ln.split())) for ln in data_lines]
        df      = pd.DataFrame(rows, columns=headers[:len(rows[0])])

        # TIME column is SAPHIR epoch (sec since 2000-01-01) — same as ENZ/CSV
        midnight_saphir = self._midnight_unix - _NC_EPOCH_OFFSET
        df.index = (df['TIME'] - midnight_saphir) / 3600.0
        df.index.name = 'Time_hours'
        df = df.drop(columns=['TIME'])
        return df

    # ── Slope analysis ────────────────────────────────────────────────────────
    def _cmp_compute_slopes(self):
        if self._cmp_df_ref is None:
            QMessageBox.information(self, "No Reference", "Load a reference dataset first.")
            return
        try:
            threshold = float(self._cmp_noise_input.text())
        except ValueError:
            threshold = 1e-6

        mcm_cols = set(self.df.columns) - {'M'}
        ref_cols = set(self._cmp_df_ref.columns) - {'M'}
        common   = sorted(mcm_cols & ref_cols)
        if not common:
            QMessageBox.warning(self, "No Common Species",
                                "No overlapping species found between MCM and reference.")
            return

        t_min = max(self.df.index.min(), self._cmp_df_ref.index.min())
        t_max = min(self.df.index.max(), self._cmp_df_ref.index.max())
        n_pts = min(500, len(self.df))
        common_t = np.linspace(t_min, t_max, n_pts)

        def interp(df, col):
            s = df[col].dropna()
            return np.interp(common_t, s.index.values.astype(float), s.values)

        # Global noise scale: 99th percentile across all common species
        py_scale = float(np.nanpercentile(
            np.abs(self.df[common].select_dtypes(include='number').values), 99))
        ref_scale = float(np.nanpercentile(
            np.abs(self._cmp_df_ref[common].select_dtypes(include='number').values), 99))
        py_scale  = py_scale  if py_scale  > 0 else 1.0
        ref_scale = ref_scale if ref_scale > 0 else 1.0

        self._cmp_slope_results = {}
        for sp in common:
            py_s  = interp(self.df, sp)
            ref_s = interp(self._cmp_df_ref, sp)

            if sp not in self._NOISE_EXCEPTIONS:
                if (np.nanmax(np.abs(py_s))  < threshold * py_scale and
                        np.nanmax(np.abs(ref_s)) < threshold * ref_scale):
                    continue

            mask = np.isfinite(py_s) & np.isfinite(ref_s)
            if mask.sum() < 5:
                continue
            px, rx = py_s[mask], ref_s[mask]
            if np.ptp(rx) == 0 and np.ptp(px) == 0:
                continue

            # Robust slope: median of point-wise ratios where ref exceeds noise floor.
            # This avoids the polyfit being dominated by a single outlier.
            noise_floor = threshold * ref_scale
            ratio_mask = (rx > noise_floor) & (px > 0)
            if ratio_mask.sum() >= 3:
                slope = float(np.median(px[ratio_mask] / rx[ratio_mask]))
            elif rx.mean() > 0:
                slope = float(px.mean() / rx.mean())
            else:
                slope = np.nan

            if not np.isfinite(slope):
                continue

            # theta: angle of the mean ratio line from the x-axis (degrees).
            # 45° = perfect 1:1, >45° MCM overestimates, <45° underestimates.
            theta = float(np.degrees(np.arctan(slope)))

            self._cmp_slope_results[sp] = {
                'slope': slope, 'theta': theta,
                'py': px, 'ref': rx, 'common_t': common_t[mask],
            }

        print(f"[Comparison] {len(self._cmp_slope_results)} species passed noise filter.")

        # Populate species combos with the result set (block signals to avoid spurious draws)
        species_list = sorted(self._cmp_slope_results.keys())
        for slot, combo in enumerate(self._cmp_ts_combos):
            prev = combo.currentText()
            combo.blockSignals(True)
            combo.clear()
            combo.addItems(species_list)
            if prev in species_list:
                combo.setCurrentText(prev)
            combo.blockSignals(False)
            # Redraw whatever is currently selected
            sp = combo.currentText()
            if sp:
                self._cmp_draw_timeseries(slot, sp)

        self._cmp_draw_slope_matrix()

    def _cmp_draw_slope_matrix(self):
        self._cmp_slope_fig.clear()
        ax = self._cmp_slope_fig.add_subplot(111)
        ax.axis('off')

        species_list = sorted(self._cmp_slope_results.keys())
        n = len(species_list)
        if n == 0:
            ax.text(0.5, 0.5, "No significant species found.\nTry lowering the noise threshold.",
                    ha='center', va='center', transform=ax.transAxes, fontsize=11)
            self._cmp_slope_canvas.draw()
            return

        # Choose a column count that gives a compact grid close to the canvas aspect.
        # Prefer slightly more columns than rows so the grid fills a landscape panel.
        ncols = max(1, min(14, int(np.ceil(np.sqrt(n * 1.6)))))
        nrows = int(np.ceil(n / ncols))

        # Fixed data-unit box size — the axes auto-stretches to fill the canvas
        # (no set_aspect('equal')), so the visual box size is determined by
        # ncols/nrows vs. canvas size, not by these data-unit values.
        box_w, box_h, pad = 1.0, 0.75, 0.08

        # Font size based on grid density: more boxes → smaller text
        density = max(ncols, nrows)
        fsize_name = max(4.0, min(9.0, 72.0 / density))
        fsize_val  = max(3.5, min(7.5, 60.0 / density))

        slopes = np.array([self._cmp_slope_results[s]['slope'] for s in species_list])
        s5, s95 = np.percentile(slopes, 5), np.percentile(slopes, 95)
        vmin = min(s5, 0.75); vmax = max(s95, 1.25)
        if vmin >= 1.0: vmin = 0.5
        if vmax <= 1.0: vmax = 1.5
        norm = TwoSlopeNorm(vcenter=1.0, vmin=vmin, vmax=vmax)
        cmap = _mcm.RdBu_r

        self._cmp_slope_box_map = []
        for i, sp in enumerate(species_list):
            row = i // ncols; col = i % ncols
            x = col * (box_w + pad)
            y = (nrows - 1 - row) * (box_h + pad)
            slope = self._cmp_slope_results[sp]['slope']
            theta = self._cmp_slope_results[sp]['theta']
            rgba  = cmap(norm(float(np.clip(slope, vmin, vmax))))
            # Highlight border if this species is selected in either slot
            edge_color = '#333'
            edge_lw    = 0.5
            if sp == self._cmp_selected[0]:
                edge_color, edge_lw = 'steelblue', 2.0
            elif sp == self._cmp_selected[1]:
                edge_color, edge_lw = 'darkorange', 2.0
            rect  = mpatches.FancyBboxPatch((x, y), box_w, box_h,
                        boxstyle="round,pad=0.02",
                        facecolor=rgba, edgecolor=edge_color, linewidth=edge_lw)
            ax.add_patch(rect)
            lum = 0.299*rgba[0] + 0.587*rgba[1] + 0.114*rgba[2]
            tc  = 'black' if lum > 0.45 else 'white'
            ax.text(x + box_w/2, y + box_h*0.67, sp,
                    ha='center', va='center', fontsize=fsize_name, fontweight='bold', color=tc)
            theta_str = f"θ={theta:.0f}°" if np.isfinite(theta) else ""
            ax.text(x + box_w/2, y + box_h*0.30, f"m={slope:.2f}  {theta_str}",
                    ha='center', va='center', fontsize=fsize_val, color=tc)
            self._cmp_slope_box_map.append((x, y, x+box_w, y+box_h, sp))

        # Tight limits — no extra whitespace around the grid
        ax.set_xlim(-pad/2, ncols*(box_w+pad) - pad/2)
        ax.set_ylim(-pad/2, nrows*(box_h+pad) - pad/2)

        sm = _mcm.ScalarMappable(cmap=cmap, norm=norm); sm.set_array([])
        self._cmp_slope_fig.colorbar(sm, ax=ax, orientation='horizontal',
                                     fraction=0.03, pad=0.01,
                                     label=f'Slope  (MCM / {self._cmp_ref_label})')
        self._cmp_slope_fig.tight_layout(pad=0.3)
        self._cmp_slope_canvas.draw()

    def _cmp_on_slope_click(self, event):
        if event.inaxes is None or event.xdata is None:
            return
        x, y = event.xdata, event.ydata
        # Left click (button 1) → slot 0 (top panel); right click (button 3) → slot 1 (bottom)
        slot = 0 if event.button == 1 else (1 if event.button == 3 else None)
        if slot is None:
            return
        for (x1, y1, x2, y2, sp) in self._cmp_slope_box_map:
            if x1 <= x <= x2 and y1 <= y <= y2:
                self._cmp_selected[slot] = sp
                combo = self._cmp_ts_combos[slot]
                combo.blockSignals(True)
                combo.setCurrentText(sp)
                combo.blockSignals(False)
                self._cmp_draw_timeseries(slot, sp)
                self._cmp_draw_slope_matrix()
                break

    def _cmp_draw_timeseries(self, slot, sp):
        """Draw a single timeseries panel for the given slot (0=top, 1=bottom)."""
        data = self._cmp_slope_results.get(sp)
        if data is None:
            return
        slope, theta = data['slope'], data['theta']
        py_vals, ref_vals, common_t = data['py'], data['ref'], data['common_t']

        # Store for crosshair lookups
        self._cmp_ts_data[slot] = (py_vals, ref_vals, common_t, sp)
        # Clear any stale crosshair artists before redrawing
        self._cmp_crosshair[slot] = {'vline': None, 'text': None}

        fig    = self._cmp_ts_figs[slot]
        canvas = self._cmp_ts_canvases[slot]
        fig.clear()
        ax = fig.add_subplot(111)

        hv_col = next((c for c in self.df.columns if c.upper() in ('SAPHIR_HV', 'HV')), None)
        if hv_col is not None:
            hv_interp = np.interp(common_t, self.df.index.values.astype(float),
                                  self.df[hv_col].values.astype(float),
                                  left=0.0, right=0.0)
            open_mask = hv_interp > 0.5
            transitions = np.diff(open_mask.astype(int))
            starts = list(np.where(transitions == 1)[0] + 1)
            ends   = list(np.where(transitions == -1)[0] + 1)
            if open_mask[0]:  starts = [0] + starts
            if open_mask[-1]: ends   = ends + [len(common_t) - 1]
            first = True
            for si, ei in zip(starts, ends):
                ax.axvspan(common_t[si], common_t[ei], alpha=0.18, color='gold',
                           label='Roof open' if first else '', zorder=0)
                first = False

        theta_str = f"θ={theta:.0f}°" if np.isfinite(theta) else ""
        ax.plot(common_t, py_vals,  label=f"{sp} MCM",
                lw=1.5, color='steelblue', zorder=2)
        ax.plot(common_t, ref_vals, label=f"{sp} {self._cmp_ref_label}",
                lw=1.5, color='darkorange', ls='--', zorder=2)
        ax.xaxis.set_major_formatter(FuncFormatter(_fmt_hhmm))
        ax.set_xlabel("Time of Day (UTC)", fontsize=9)
        ax.set_ylabel(sp, fontsize=9)
        ax.set_title(f"{sp}   m={slope:.3f}   {theta_str}", fontweight='bold', fontsize=10)
        ax.legend(loc='upper right', fontsize=7)
        ax.grid(True, ls='--', alpha=0.35)

        fig.tight_layout(pad=0.8)
        canvas.draw()

    # ── Crosshair cursor ─────────────────────────────────────────────────────

    def _cmp_on_mouse_move(self, event, slot):
        if event.inaxes is None:
            self._cmp_clear_crosshair(slot)
            return
        other  = 1 - slot
        ax     = event.inaxes
        canvas = self._cmp_ts_canvases[slot]
        ch     = self._cmp_crosshair[slot]
        x      = event.xdata

        # Remove previous own artists
        if ch['vline'] is not None:
            try: ch['vline'].remove()
            except Exception: pass
        if ch['text'] is not None:
            try: ch['text'].remove()
            except Exception: pass

        # Remove previous mirror artists in the peer canvas
        if self._cmp_mirror_vlines[other] is not None:
            try: self._cmp_mirror_vlines[other].remove()
            except Exception: pass
            self._cmp_mirror_vlines[other] = None
        if self._cmp_mirror_texts[other] is not None:
            try: self._cmp_mirror_texts[other].remove()
            except Exception: pass
            self._cmp_mirror_texts[other] = None

        ch['vline'] = ax.axvline(x, color='gray', lw=0.8, ls='--', zorder=10)

        data = self._cmp_ts_data[slot]
        if data is not None:
            py_vals, ref_vals, common_t, sp = data
            py_at  = float(np.interp(x, common_t, py_vals))
            ref_at = float(np.interp(x, common_t, ref_vals))
            h, rem = divmod(max(0, int(x * 3600)), 3600)
            m      = rem // 60
            pct_txt = f"   Agr={py_at/ref_at*100:.1f}%" if ref_at != 0 else ""
            txt    = f"{h:02d}:{m:02d}   MCM={py_at:.3g}   Ref={ref_at:.3g}{pct_txt}"
            ch['text'] = ax.text(
                0.02, 0.98, txt, transform=ax.transAxes,
                va='top', ha='left', fontsize=8,
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.85),
                zorder=11)

        # Draw mirror vline + text in the peer canvas
        other_axes = self._cmp_ts_figs[other].axes
        if other_axes:
            other_ax = other_axes[0]
            self._cmp_mirror_vlines[other] = other_ax.axvline(
                x, color='gray', lw=0.8, ls='--', zorder=10)
            other_data = self._cmp_ts_data[other]
            if other_data is not None:
                py_vals2, ref_vals2, common_t2, sp2 = other_data
                py_at2  = float(np.interp(x, common_t2, py_vals2))
                ref_at2 = float(np.interp(x, common_t2, ref_vals2))
                h2, rem2 = divmod(max(0, int(x * 3600)), 3600)
                m2       = rem2 // 60
                pct_txt2 = f"   Agr={py_at2/ref_at2*100:.1f}%" if ref_at2 != 0 else ""
                txt2     = f"{h2:02d}:{m2:02d}   MCM={py_at2:.3g}   Ref={ref_at2:.3g}{pct_txt2}"
                self._cmp_mirror_texts[other] = other_ax.text(
                    0.02, 0.98, txt2, transform=other_ax.transAxes,
                    va='top', ha='left', fontsize=8,
                    bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.85),
                    zorder=11)
            self._cmp_ts_canvases[other].draw_idle()

        canvas.draw_idle()

    def _cmp_clear_crosshair(self, slot):
        other  = 1 - slot
        ch     = self._cmp_crosshair[slot]
        canvas = self._cmp_ts_canvases[slot]
        if ch['vline'] is not None:
            try: ch['vline'].remove()
            except Exception: pass
            ch['vline'] = None
        if ch['text'] is not None:
            try: ch['text'].remove()
            except Exception: pass
            ch['text'] = None
        # Also remove the mirror vline + text we drew in the peer canvas
        peer_dirty = False
        if self._cmp_mirror_vlines[other] is not None:
            try: self._cmp_mirror_vlines[other].remove()
            except Exception: pass
            self._cmp_mirror_vlines[other] = None
            peer_dirty = True
        if self._cmp_mirror_texts[other] is not None:
            try: self._cmp_mirror_texts[other].remove()
            except Exception: pass
            self._cmp_mirror_texts[other] = None
            peer_dirty = True
        if peer_dirty:
            self._cmp_ts_canvases[other].draw_idle()
        canvas.draw_idle()

    # ─────────────────────────────────────────────────────────────────────────
    # Budgets tab  — Carbon (top) + Nitrogen (bottom)
    # ─────────────────────────────────────────────────────────────────────────
    def _build_carbon_tab(self):
        """Build the unified Budgets tab: Carbon on top, Nitrogen on bottom."""
        outer = QWidget()
        outer_layout = QVBoxLayout(outer)
        outer_layout.setContentsMargins(2, 2, 2, 2)
        outer_layout.setSpacing(2)

        # ── Shared controls row ───────────────────────────────────────────
        ctrl_layout = QHBoxLayout()
        ctrl_layout.setContentsMargins(0, 0, 0, 0)
        ctrl_layout.setSpacing(6)

        ctrl_layout.addWidget(QLabel("Display:"))
        self._cb_mode_combo = QComboBox()
        self._cb_mode_combo.addItems(["All Species", "Grouped by Class"])
        ctrl_layout.addWidget(self._cb_mode_combo)

        self._cb_total_chk = QCheckBox("Total C")
        self._cb_total_chk.setChecked(True)
        ctrl_layout.addWidget(self._cb_total_chk)

        self._nb_total_chk = QCheckBox("Total N")
        self._nb_total_chk.setChecked(True)
        ctrl_layout.addWidget(self._nb_total_chk)

        btn = QPushButton("Plot")
        btn.clicked.connect(self._update_both_budgets)
        ctrl_layout.addWidget(btn)

        self._cb_status_lbl = QLabel("")
        self._cb_status_lbl.setStyleSheet("color: #777; font-size: 10px;")
        ctrl_layout.addWidget(self._cb_status_lbl)

        ctrl_layout.addStretch()
        outer_layout.addLayout(ctrl_layout, 0)

        # ── Vertical splitter: Carbon top / Nitrogen bottom ───────────────
        vsplit = QSplitter(Qt.Orientation.Vertical)

        # Carbon panel
        carbon_widget = QWidget()
        carbon_layout = QVBoxLayout(carbon_widget)
        carbon_layout.setContentsMargins(0, 0, 0, 0)
        carbon_layout.setSpacing(0)
        self._cb_fig = Figure(figsize=(12, 4))
        self._cb_canvas = FigureCanvas(self._cb_fig)
        self._cb_toolbar = NavigationToolbar(self._cb_canvas, carbon_widget)
        carbon_layout.addWidget(self._cb_toolbar)
        carbon_layout.addWidget(self._cb_canvas, 1)
        vsplit.addWidget(carbon_widget)

        # Nitrogen panel
        nitrogen_widget = QWidget()
        nitrogen_layout = QVBoxLayout(nitrogen_widget)
        nitrogen_layout.setContentsMargins(0, 0, 0, 0)
        nitrogen_layout.setSpacing(0)
        self._nb_fig = Figure(figsize=(12, 4))
        self._nb_canvas = FigureCanvas(self._nb_fig)
        self._nb_toolbar = NavigationToolbar(self._nb_canvas, nitrogen_widget)
        nitrogen_layout.addWidget(self._nb_toolbar)
        nitrogen_layout.addWidget(self._nb_canvas, 1)
        vsplit.addWidget(nitrogen_widget)

        vsplit.setSizes([500, 500])
        outer_layout.addWidget(vsplit, 1)

        # Auto-replot when controls change
        self._cb_mode_combo.currentIndexChanged.connect(self._update_both_budgets)
        self._cb_total_chk.stateChanged.connect(self._update_carbon_budget)
        self._nb_total_chk.stateChanged.connect(self._update_nitrogen_budget)

        self.right_tabs.addTab(outer, "Budgets")

    def _update_both_budgets(self):
        self._update_carbon_budget()
        self._update_nitrogen_budget()

    def _update_carbon_budget(self):
        """Compute carbon-weighted concentrations and draw the stacked timeseries."""
        self._cb_fig.clear()
        ax = self._cb_fig.add_subplot(111)

        t = self.df.index.values          # hours since midnight
        M = self.df['M'].values if 'M' in self.df.columns else None

        # ── Parent VOC set from config (burst injections only) ────────────
        parent_voc = set()
        if self.budget_engine and self.budget_engine.cfg:
            for inj in self.budget_engine.cfg.get('burst_injections', []):
                sp = inj.get('species', '')
                if sp and _carbon_count(sp) > 0:
                    parent_voc.add(sp.upper())

        # ── Collect carbon-containing MCM species ─────────────────────────
        species_data = {}    # name → (n_carbons, ppb_C array)
        unknown_species = [] # non-zero species with unresolvable C count
        co2_present    = False

        for col in self.df.columns:
            if col == 'M' or col.upper().startswith('SAPHIR_'):
                continue
            # Always exclude CO2 — MCM treats it as a sink, not a tracked tracer
            if col.upper() == 'CO2':
                vals = self.df[col].values
                ppb  = (vals / M) * 1e9 if M is not None else vals.copy()
                if np.nanmax(np.abs(ppb)) > 1e-3:
                    co2_present = True
                continue
            n_c = _carbon_count(col)
            vals = self.df[col].values
            ppb = (vals / M) * 1e9 if M is not None else vals.copy()
            if n_c == 0:
                if np.nanmax(np.abs(ppb)) > 1e-3:
                    unknown_species.append(col)
                continue
            ppb_c = np.clip(ppb, 0.0, None) * n_c
            if np.nanmax(ppb_c) < 1e-9:
                continue
            species_data[col] = (n_c, ppb_c)

        # ── Status label ──────────────────────────────────────────────────
        if unknown_species:
            shown = ', '.join(unknown_species[:6])
            tail = f' +{len(unknown_species)-6} more' if len(unknown_species) > 6 else ''
            self._cb_status_lbl.setText(f"⚠ Unknown C count: {shown}{tail}")
            self._cb_status_lbl.setToolTip(
                "Non-zero species with unrecognised carbon count "
                "(add to _MCM_CARBON_LOOKUP if needed):\n" +
                '\n'.join(unknown_species))
        else:
            self._cb_status_lbl.setText("✓ All non-zero species accounted for")
            self._cb_status_lbl.setToolTip("")

        if not species_data:
            ax.text(0.5, 0.5, "No carbon-containing species found in simulation output.",
                    ha='center', va='center', transform=ax.transAxes, fontsize=11)
            self._cb_canvas.draw()
            return

        # ── HV open-roof shading ──────────────────────────────────────────
        hv_col = next((c for c in self.df.columns if c.upper() in ('SAPHIR_HV', 'HV')), None)
        if hv_col is not None:
            hv = self.df[hv_col].values
            in_open, open_start = False, None
            for ti, vi in zip(t, hv):
                if vi > 0.5 and not in_open:
                    in_open, open_start = True, ti
                elif vi <= 0.5 and in_open:
                    ax.axvspan(open_start, ti, color='#fffacd', alpha=0.35,
                               zorder=0, linewidth=0)
                    in_open = False
            if in_open:
                ax.axvspan(open_start, t[-1], color='#fffacd', alpha=0.35,
                           zorder=0, linewidth=0)

        # ── Draw stacked areas ────────────────────────────────────────────
        mode = self._cb_mode_combo.currentText()
        if mode == "All Species":
            self._cb_draw_all_species(ax, t, species_data)
        else:
            self._cb_draw_grouped(ax, t, species_data, parent_voc)

        # ── Total C overlay ───────────────────────────────────────────────
        if self._cb_total_chk.isChecked():
            total = sum(np.nan_to_num(ppb_c) for _, ppb_c in species_data.values())
            ax.plot(t, total, color='black', lw=1.5, ls='--',
                    label='Total tracked C', zorder=10)

        # ── Exclusion annotation (bottom-left corner) ─────────────────────
        excl_lines = []
        if co2_present:
            excl_lines.append('• CO₂  (MCM sink — not a tracked tracer)')
        for sp in unknown_species:
            excl_lines.append(f'• {sp}  (carbon count unknown)')
        if excl_lines:
            note = 'Excluded from budget:\n' + '\n'.join(excl_lines)
            ax.text(0.02, 0.02, note, transform=ax.transAxes, fontsize=7.5,
                    va='bottom', ha='left', family='monospace',
                    bbox=dict(boxstyle='round,pad=0.4', facecolor='#fff8dc',
                              alpha=0.88, edgecolor='#ccaa00'))

        # ── Axes formatting ───────────────────────────────────────────────
        ax.set_xlabel('Time of Day (UTC)')
        ax.set_ylabel('Carbon content (ppb-C)')
        ax.xaxis.set_major_formatter(FuncFormatter(_fmt_hhmm))
        ax.set_xlim(t[0], t[-1])
        ax.set_ylim(bottom=0)
        ax.grid(True, linestyle='--', alpha=0.4)

        handles, labels = ax.get_legend_handles_labels()
        if handles:
            if mode == "All Species" and len(handles) > 15:
                ncol = max(1, len(handles) // 20 + 1)
                ax.legend(handles, labels, loc='upper right',
                          fontsize=6, ncol=ncol, framealpha=0.85)
            else:
                ax.legend(handles, labels, loc='upper right',
                          fontsize=8, framealpha=0.85)

        self._cb_fig.tight_layout()
        self._cb_canvas.draw()

    def _cb_draw_all_species(self, ax, t, species_data):
        """Stacked area with one band per species, sorted by mean contribution."""
        sorted_items = sorted(species_data.items(),
                              key=lambda x: np.nanmean(x[1][1]), reverse=True)
        names  = [item[0] for item in sorted_items]
        arrays = [np.nan_to_num(item[1][1]) for item in sorted_items]

        cmap   = _mcm.get_cmap('tab20')
        colors = [cmap((i % 20) / 20.0) for i in range(len(names))]

        ax.stackplot(t, *arrays, labels=names, colors=colors, alpha=0.85)

    def _cb_draw_grouped(self, ax, t, species_data, parent_voc):
        """Stacked area with one band per chemical class using fixed class colours."""
        class_totals = {}
        for sp, (_, ppb_c) in species_data.items():
            cls = _carbon_class(sp, parent_voc)
            if cls not in class_totals:
                class_totals[cls] = np.zeros(len(t))
            class_totals[cls] += np.nan_to_num(ppb_c)

        names, arrays, colors = [], [], []
        # Emit in defined order first
        for cls in _CB_CLASS_ORDER:
            if cls in class_totals and np.nanmax(class_totals[cls]) > 1e-9:
                names.append(cls)
                arrays.append(class_totals[cls])
                colors.append(_CB_CLASS_COLORS.get(cls, '#aaaaaa'))
        # Append any classes not in the predefined order
        for cls, arr in class_totals.items():
            if cls not in _CB_CLASS_ORDER and np.nanmax(arr) > 1e-9:
                names.append(cls)
                arrays.append(arr)
                colors.append('#aaaaaa')

        if arrays:
            ax.stackplot(t, *arrays, labels=names, colors=colors, alpha=0.85)

    # ── Nitrogen budget ───────────────────────────────────────────────────────
    def _update_nitrogen_budget(self):
        """Compute nitrogen-weighted concentrations and draw the stacked timeseries."""
        self._nb_fig.clear()
        ax = self._nb_fig.add_subplot(111)

        t = self.df.index.values
        M = self.df['M'].values if 'M' in self.df.columns else None

        species_data = {}   # name → (n_N, ppb_N array)
        unknown_species = []

        for col in self.df.columns:
            if col == 'M' or col.upper().startswith('SAPHIR_'):
                continue
            n_n = _nitrogen_count(col)
            if n_n == 0:
                continue
            vals = self.df[col].values
            ppb  = (vals / M) * 1e9 if M is not None else vals.copy()
            ppb_n = np.clip(ppb, 0.0, None) * n_n
            if np.nanmax(ppb_n) < 1e-9:
                continue
            species_data[col] = (n_n, ppb_n)

        # HV shading
        hv_col = next((c for c in self.df.columns if c.upper() in ('SAPHIR_HV', 'HV')), None)
        if hv_col is not None:
            hv = self.df[hv_col].values
            in_open, open_start = False, None
            for ti, vi in zip(t, hv):
                if vi > 0.5 and not in_open:
                    in_open, open_start = True, ti
                elif vi <= 0.5 and in_open:
                    ax.axvspan(open_start, ti, color='#fffacd', alpha=0.35,
                               zorder=0, linewidth=0)
                    in_open = False
            if in_open:
                ax.axvspan(open_start, t[-1], color='#fffacd', alpha=0.35,
                           zorder=0, linewidth=0)

        if not species_data:
            ax.text(0.5, 0.5, "No nitrogen-containing species found in simulation output.",
                    ha='center', va='center', transform=ax.transAxes, fontsize=11)
            self._nb_canvas.draw()
            return

        mode = self._cb_mode_combo.currentText()
        if mode == "All Species":
            sorted_items = sorted(species_data.items(),
                                  key=lambda x: np.nanmean(x[1][1]), reverse=True)
            names  = [item[0] for item in sorted_items]
            arrays = [np.nan_to_num(item[1][1]) for item in sorted_items]
            cmap   = _mcm.get_cmap('tab20')
            colors = [cmap((i % 20) / 20.0) for i in range(len(names))]
            ax.stackplot(t, *arrays, labels=names, colors=colors, alpha=0.85)
        else:
            class_totals = {}
            for sp, (_, ppb_n) in species_data.items():
                cls = _nitrogen_class(sp)
                if cls not in class_totals:
                    class_totals[cls] = np.zeros(len(t))
                class_totals[cls] += np.nan_to_num(ppb_n)
            names, arrays, colors = [], [], []
            for cls in _NB_CLASS_ORDER:
                if cls in class_totals and np.nanmax(class_totals[cls]) > 1e-9:
                    names.append(cls)
                    arrays.append(class_totals[cls])
                    colors.append(_NB_CLASS_COLORS.get(cls, '#aaaaaa'))
            for cls, arr in class_totals.items():
                if cls not in _NB_CLASS_ORDER and np.nanmax(arr) > 1e-9:
                    names.append(cls)
                    arrays.append(arr)
                    colors.append('#aaaaaa')
            if arrays:
                ax.stackplot(t, *arrays, labels=names, colors=colors, alpha=0.85)

        if self._nb_total_chk.isChecked():
            total = sum(np.nan_to_num(ppb_n) for _, ppb_n in species_data.values())
            ax.plot(t, total, color='black', lw=1.5, ls='--',
                    label='Total tracked N', zorder=10)

        ax.set_xlabel('Time of Day (UTC)')
        ax.set_ylabel('Nitrogen content (ppb-N)')
        ax.xaxis.set_major_formatter(FuncFormatter(_fmt_hhmm))
        ax.set_xlim(t[0], t[-1])
        ax.set_ylim(bottom=0)
        ax.grid(True, linestyle='--', alpha=0.4)

        handles, labels = ax.get_legend_handles_labels()
        if handles:
            if mode == "All Species" and len(handles) > 15:
                ax.legend(handles, labels, loc='upper right',
                          fontsize=6, ncol=max(1, len(handles) // 20 + 1),
                          framealpha=0.85)
            else:
                ax.legend(handles, labels, loc='upper right',
                          fontsize=8, framealpha=0.85)

        self._nb_fig.tight_layout()
        self._nb_canvas.draw()

    # ── Default plots on startup ──────────────────────────────────────────────
    def _plot_defaults(self):
        _voc_candidates = ['C2H5OH', 'TOLUENE', 'ISOPRENE', 'ALPHA_PINENE', 'BENZENE',
                           'PROPANE', 'BUTANE', 'PENTANE', 'HEXANE']
        for _voc in _voc_candidates:
            if _voc in self.df.columns:
                self.plot_variable(_voc, "Left Axis (Offset)", "Convert to ppb")
                break

        for _sp in ['NO', 'NO2', 'O3']:
            if _sp in self.df.columns:
                self.plot_variable(_sp, "Left Axis (Main)", "Convert to ppb")
        for _sp in ['SAPHIR_NO', 'SAPHIR_NO2', 'SAPHIR_O3']:
            if _sp in self.df.columns:
                self.plot_variable(_sp, "Left Axis (Main)", "Native (molec/cm³ or SAPHIR unit)")

        for _sp in ['OH', 'HO2']:
            if _sp in self.df.columns:
                self.plot_variable(_sp, "Right Axis (Main)", "Native (molec/cm³ or SAPHIR unit)")

        if any(sp in self.df.columns for sp in ['OH', 'HO2']):
            self.ax2.set_yscale('log')
            self.ax2.set_ylim(bottom=1e4)

        _hv_col = next((c for c in self.df.columns if c.upper() in ('SAPHIR_HV', 'HV')), None)
        if _hv_col is not None:
            _hv = self.df[_hv_col].values
            _t  = self.df.index.values
            _in_open, _open_start = False, None
            for t, v in zip(_t, _hv):
                if v > 0.5 and not _in_open:
                    _in_open, _open_start = True, t
                elif v <= 0.5 and _in_open:
                    self.ax1.axvspan(_open_start, t, color='#fffacd', alpha=0.35,
                                     zorder=0, linewidth=0)
                    _in_open = False
            if _in_open:
                self.ax1.axvspan(_open_start, _t[-1], color='#fffacd', alpha=0.35,
                                 zorder=0, linewidth=0)

        self._update_legends()
        self._update_ylabels()
        self.canvas.draw()


# ─────────────────────────────────────────────────────────────────────────────
# Entry points
# ─────────────────────────────────────────────────────────────────────────────
def show_results(out_path):
    if HAS_PYQT:
        print(">>> [VISUALIZATION] Launching Interactive PyQt Viewer...")
        df = pd.read_csv(out_path, index_col=0)
        app = QApplication.instance()
        if app is None:
            app = QApplication(sys.argv)
        viewer = ResultsViewer(df, csv_path=out_path)
        viewer.show()
        app.exec()
        return

    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("[WARNING] Missing matplotlib. Visualization skipped.")
        return

    print("[WARNING] PyQt6 is not installed. Falling back to simple Matplotlib visualization.")
    df = pd.read_csv(out_path, index_col=0)
    midnight_unix = df.index[0] - (df.index[0] % 86400)
    t_hours = (df.index - midnight_unix) / 3600.0

    fig, ax1 = plt.subplots(figsize=(12, 6))
    ax1.set_xlabel('Time of Day (UTC)')
    ax1.xaxis.set_major_formatter(plt.FuncFormatter(_fmt_hhmm))
    ax1.set_ylabel('Mixing Ratio (ppb)', color='black')

    for spec in ['O3', 'NO', 'NO2', 'HONO', 'C2H5OH', 'HCHO']:
        if spec in df.columns:
            ax1.plot(t_hours, (df[spec] / df['M']) * 1e9, label=f"{spec} (ppb)")

    ax1.grid(True, linestyle='--', alpha=0.6)
    if any(r in df.columns for r in ['OH', 'HO2', 'RO2']):
        ax2 = ax1.twinx()
        ax2.set_ylabel('Radicals (molec/cm³)', color='red')
        for r in ['OH', 'HO2', 'RO2']:
            if r in df.columns:
                ax2.plot(t_hours, df[r], label=f"{r} (molec/cm³)", linestyle='--')

    fig.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05), ncol=5)
    plt.title('SAPHIR Scenario Simulation')
    plt.tight_layout()
    plt.show()


def main():
    """Standalone viewer for any simulation_results.csv file."""
    import argparse
    parser = argparse.ArgumentParser(
        description="Interactive viewer for SAPHIR-MCM simulation results.")
    parser.add_argument("results_csv_path", nargs='?', default=None,
                        help="Path to the simulation_results.csv file to view.")
    args = parser.parse_args()

    if args.results_csv_path and os.path.exists(args.results_csv_path):
        print(f">>> [VIEWER] Loading results from: {args.results_csv_path}")
        show_results(args.results_csv_path)
    else:
        if args.results_csv_path:
            print(f"[ERROR] File not found: {args.results_csv_path}")
        print("\nUsage: python results_viewer.py <path_to_simulation_results.csv>")
        print(">>> [VIEWER] No path provided or path is invalid. Exiting.")


if __name__ == "__main__":
    if HAS_PYQT:
        _app = QApplication(sys.argv)
        _pm = QPixmap(os.path.join(os.path.dirname(os.path.abspath(__file__)), "Juelich Logo.png"))
        if not _pm.isNull():
            _app.setWindowIcon(QIcon(_pm))
    main()
