import sys
import os
import re
import json
import warnings

# Started as a program: say what is happening before pandas / matplotlib / the viewer load
_LOADER = None
if __name__ == "__main__":
    try:
        from loading_box import LoadingBox
        _LOADER = LoadingBox.start("opening results")
        _LOADER.step("Loading the plotting and data libraries", None, "pandas, matplotlib and Qt are imported once per window; this takes a few seconds")
    except Exception:
        _LOADER = None

import pandas as pd
import numpy as np
from datetime import datetime
from column_file_reader import load_column_file


def _interact(canvas, hover: bool = True):
    """Overview-style mouse interaction (wheel/Ctrl+wheel/right-drag zoom, drag pan, double-click reset,
    hover readout) on any canvas of the viewer; silently skipped if the helper is unavailable."""
    try:
        from plot_interaction import attach
        attach(canvas, hover=hover)
    except Exception as e:
        print(f"[viewer] interaction not attached: {e}")


def _progress(message: str, percent=None, why: str = "") -> None:
    """Report a start-up stage to the loading window when there is one (no-op otherwise)."""
    if _LOADER is not None:
        try:
            _LOADER.step(message, percent, why)
        except Exception:
            pass


_PARSE_CACHE: dict = {}


def _parse_easy_cached(parser, path: str):
    """EasyParser.parse(path) once per file per process: the viewer needs the same mechanism for
    the species classes and for the budget engine (3–4 s each for a full elog)."""
    key = (os.path.abspath(path),)                  # same file ⇒ same mechanism (elogs are flat; include files resolve by their own folder)
    try:
        st = os.stat(path); key += (st.st_mtime_ns, st.st_size)
    except OSError:
        pass
    m = _PARSE_CACHE.get(key)
    if m is None:
        m = parser.parse(path); _PARSE_CACHE[key] = m
    return m


def _is_juice_output(path: str) -> bool:
    """Return True if path is a JUICE output file (first line starts with '# Calculated by JUICE')."""
    try:
        with open(path, encoding='latin-1') as _f:
            return _f.readline().strip().startswith('# Calculated by JUICE')
    except Exception:
        return False

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
    # single source of truth for formula/MW/RO2/SMILES: species_db (SpeciesDB is a drop-in for the old
    # MCM_DB -- .get(name) returns a Species with .formula/.n_C/.n_H/.n_O/.n_N/.mass/.is_peroxy, and
    # `name in MCM_DB` works). KNOWN_ZERO_CARBON stays as the inorganic-constant set.
    from species_db import load_species_db
    from mcm_species_db import KNOWN_ZERO_CARBON
    MCM_DB = load_species_db()
    if len(getattr(MCM_DB, "_t", {})) == 0:          # no species table found at all → the MCM export
        import mcm_species_db as _msd
        MCM_DB = _msd.MCM_DB
        print(f"[viewer] species_db empty — using mcm_species_db ({len(MCM_DB)} species)")
except ImportError:
    MCM_DB = {}
    KNOWN_ZERO_CARBON = frozenset()

try:
    from PyQt6.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
                                 QTreeWidget, QTreeWidgetItem, QListWidget, QListWidgetItem,
                                 QPushButton, QLabel, QComboBox, QSplitter, QLineEdit,
                                 QTabWidget, QFileDialog, QDialog, QDialogButtonBox,
                                 QTableWidget, QTableWidgetItem, QCheckBox, QHeaderView,
                                 QInputDialog, QMessageBox, QAbstractItemView, QSpinBox, QDoubleSpinBox, QToolButton, QSizePolicy)
    from PyQt6.QtCore import Qt, QTimer
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


def _distinct_colors(labels):
    """One clearly different colour per label (up to 60 without repeating: tab10, then the
    tab20b/tab20c dark members, then their light members); 'Other' is always grey."""
    import matplotlib.cm as _cm
    base = [_cm.tab10(i) for i in range(10)]
    for name in ("tab20b", "tab20c"):
        base += [_cm.get_cmap(name)(i) for i in range(0, 20, 4)]      # dark members
    for name in ("tab20b", "tab20c"):
        base += [_cm.get_cmap(name)(i) for i in range(2, 20, 4)]      # light members
    out, i = [], 0
    for lbl in labels:
        if str(lbl).lower().startswith("other"):
            out.append("#9a9a9a")
        else:
            out.append(base[i % len(base)]); i += 1
    return out

_HATCHES = ['']            # no patterns in the budgets (Georgios 2026-10-22): shades only


def _shade(color, k: int, n: int):
    """Member k of n of one colour family as a ramp from very light (near white) through the base
    colour to very dark (near black) — no patterns, so many injected precursors stay distinguishable."""
    import matplotlib.colors as _mcolors
    r, g, b = _mcolors.to_rgb(color)
    w = min(0.5, 0.12 * (n - 1))                               # few members: stay near the base colour
    p = 0.5 if n <= 1 else 0.5 - w + 2 * w * k / (n - 1)       # 0 = lightest … 1 = darkest
    if p <= 0.5:
        f = 0.88 * (1 - p / 0.5)                               # mix with white
        return (r + (1 - r) * f, g + (1 - g) * f, b + (1 - b) * f)
    f = 0.8 * (p - 0.5) / 0.5                                  # mix with black
    return (r * (1 - f), g * (1 - f), b * (1 - f))


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

def _nitrogen_count(name: str) -> int:
    """Return the number of nitrogen atoms for an MCM species.

    Primary source: MCM_DB formula (exact, from mcm_export_species.tsv).
    Fallback: sequential functional-group scan on the MCM name (longest
    pattern first, to avoid matching N in ETHANOL/METHANOL etc.).
    """
    upper = name.upper()
    entry = MCM_DB.get(upper)
    if entry is not None:
        return entry.n_N

    # Fallback — functional-group scan on name for species not in DB
    _N_GROUPS = [
        ('N2O5', 2), ('CO3NO2', 1), ('ONO2', 1),
        ('HNO3', 1), ('HNO4', 1), ('HONO', 1),
        ('NO3',  1), ('NO2',  1), ('NO',   1),
    ]
    tmp, total = upper, 0
    for pattern, n_atoms in _N_GROUPS:
        while pattern in tmp:
            total += n_atoms
            tmp = tmp.replace(pattern, '\x00' * len(pattern), 1)
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
    if re.search(r'PAN\d*$', n) or n.endswith('CO3NO2') or n in ('PPN', 'PBZN', 'PHAN', 'APAN', 'GPAN'):   # MCM peroxyacyl nitrates (…PAN, …PAN5, PPN, PBZN, PHAN); A2PANO/…PANOO are radicals
        return 'PANs'
    if ('NO3' in n and n != 'NO3') or 'ONO2' in n or n.endswith('NO2'):
        return 'Org. Nitrates'
    return 'Other N'

# ── Carbon budget helpers ─────────────────────────────────────────────────────

_CB_OTHER_SPLIT = 0.05   # a species classed "Other" with a peak ≥ 5 % of the total carbon peak is shown on its own
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


def _carbon_count(name: str) -> int:
    """Return the number of carbon atoms for an MCM species.

    Primary source: MCM_DB formula (exact, from mcm_export_species.tsv).
    Fallback: regex on the MCM name — C(n) pattern, skipping Cl.
    Returns 0 for inorganic/unrecognised species.
    """
    upper = name.upper()
    entry = MCM_DB.get(upper)
    if entry is not None:
        return entry.n_C

    # Fallback — regex on name for species not in DB
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
    # inorganic nitrogen is never an organic nitrate, whatever the name suggests
    if n in ('HNO3', 'HO2NO2', 'N2O5', 'NO3', 'NO2', 'NO', 'HONO', 'NA', 'SA', 'SO2', 'SO3', 'HSO3', 'H2O2', 'O3', 'OH', 'HO2', 'H2', 'H2O', 'CL', 'CLNO2', 'CLONO2'):
        return 'Other'
    # Parent VOC (injected organic precursors from config)
    if parent_voc and n in parent_voc:
        return 'Parent VOC'
    # Simple C1 inorganic
    if n == 'CO':
        return 'CO'
    if n == 'CO2':
        return 'Other'          # rarely appears in MCM CSV
    # PAN family: peroxyacyl nitrates (e.g. CH3CO3NO2, MPAN)
    if re.search(r'PAN\d*$', n) or n.endswith('CO3NO2') or n in ('PPN', 'PBZN', 'PHAN', 'APAN', 'GPAN'):   # MCM peroxyacyl nitrates (…PAN, …PAN5, PPN, PBZN, PHAN); A2PANO/…PANOO are radicals
        return 'PANs'
    # Organic nitrates (RONO2, names carry NO3 or ONO2) and peroxy nitrates (RO2NO2)
    if ('NO3' in n and n != 'NO3') or 'ONO2' in n or n.endswith('O2NO2'):
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
    # Alcohols (ROH, after OOH/CO2H/COOH already caught above); MCM also names them …OL / …GLY
    if n.endswith('OH') or n.endswith('OL') or n.endswith('GLY'):
        return 'Alcohols'
    # MCM aldehyde names ending in …AL (PINAL, MALDIAL, C4MDIAL), ketones …ONE / …CARB (dicarbonyls)
    if n.endswith('AL') and not n.endswith('SAL'):
        return 'Aldehydes'
    if n.endswith('ONE') or n.endswith('CARB'):
        return 'Ketones'
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
    def load(self, config_path, include_dirs=None):
        if not HAS_BUDGET:
            raise ImportError("easy_parser / mechanism / mcm_kinetics not importable")
        self._include_dirs = list(include_dirs or [])

        with open(config_path) as f:
            self.cfg = json.load(f)

        self.root = self.cfg.get('folder') or os.path.dirname(os.path.abspath(config_path))   # IDL-imported JSONs carry folder=''
        chem_path = os.path.join(self.root, 'chemistry')
        parser = EasyParser(base_path=chem_path)
        self.mechanism = Mechanism()
        self.var_assignments = []

        mech_files = sorted(self.cfg['mechanism_files'], key=lambda x: '.rate' not in x)
        for f_name in mech_files:
            f_path = os.path.join(chem_path, f_name)
            if not os.path.exists(f_path):
                continue
            m_part = _parse_easy_cached(parser, f_path)
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

        # Elog-mode JSONs: parse the canonical elog (a single flat-expanded file).
        # Mirrors scenario_runner's elog_path handling; reactions/params merge in.
        elog_path = self.cfg.get('elog_path', '')
        if elog_path:
            if not os.path.isabs(elog_path):
                elog_path = os.path.join(self.root, elog_path)
            if not os.path.exists(elog_path):
                elog_path = self._regenerate_elog(config_path, elog_path) or elog_path
            if os.path.exists(elog_path):
                m_part = _parse_easy_cached(parser, elog_path)
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

        # ELOG mode (every run written by the step window / IDL import): the canonical elog IS the
        # whole mechanism.  The JSON's custom_reactions, chamber_source_reactions, injection_reactions
        # and artifacts are the parser's informational copies of the same EASY statements — the runner
        # skips them in this mode (scenario_runner: `elif not _using_elog`, `_add_artifacts`), so the
        # budget engine must too, or every reaction is counted twice (kOH stack 2× the .pd kOH).
        _using_elog = bool(elog_path) and os.path.exists(elog_path)

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
        arts = self.cfg.get('artifacts', {}) if not _using_elog else {}
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
        for r in (self.cfg.get('custom_reactions', []) if not _using_elog else []):
            if not r.get('reactants') or not r.get('products') or not r.get('rate'):
                continue
            reacs = [x.strip().upper() for x in r['reactants'].split('+') if x.strip()]
            prods = [x.strip().upper() for x in r['products'].split('+') if x.strip()]
            self.mechanism.reactions.append(Reaction(reacs, prods, r['rate']))

        # IDL chamber_source_reactions (SHONO/SHCHO/H2O)
        for r in (self.cfg.get('chamber_source_reactions', []) if not _using_elog else []):
            if not r.get('products') or not r.get('rate'):
                continue
            reacts_str = r.get('reactants', '')
            reacs = [x.strip().upper() for x in reacts_str.split('+') if x.strip()] if reacts_str.strip() else []
            prods = [x.strip().upper() for x in r['products'].split('+') if x.strip()]
            self.mechanism.reactions.append(Reaction(reacs, prods, r['rate']))

        # IDL injection_reactions
        enz_cols_upper = set()
        enz_path_rv = next((os.path.join(self.root, f) for f in os.listdir(self.root)
                            if f.lower().endswith('.enz') and 'kammer' in f.lower()), None)
        if enz_path_rv and os.path.exists(enz_path_rv):
            try:
                with open(enz_path_rv, 'r', encoding='latin-1') as _fh:
                    for _line in _fh:
                        if _line.startswith('COLUMN'):
                            _parts = _line.split('=', 1)
                            _name = _parts[1].split('/')[0].split('(')[0].strip().upper()
                            enz_cols_upper.add(_name)
                        if '&&' in _line:
                            break
            except Exception:
                pass
        for r in (self.cfg.get('injection_reactions', []) if not _using_elog else []):
            if not r.get('products') or not r.get('rate'):
                continue
            reacts_str = r.get('reactants', '')
            reacs_raw = [x.strip().upper() for x in reacts_str.split('+') if x.strip()] if reacts_str.strip() else []
            # Skip ENZ switch columns as reactants — they are constrained inputs, not ODE species
            reacs = [s for s in reacs_raw
                     if re.sub(r'^[\d\.]+\s*\*\s*', '', s) not in enz_cols_upper]
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
            excl_keys = set(); overrides = []
            for e in excl:
                key = (_canon([x.strip() for x in e['reactants'].split('+') if x.strip()]),
                       _canon([x.strip() for x in e['products'].split('+') if x.strip()]))
                if e.get('action', 'Remove') == 'Modify Rate' and e.get('rate', '').strip():
                    overrides.append((key, e['rate'].strip()))      # as the runner: patch the rate
                else:
                    excl_keys.add(key)
            self.mechanism.reactions = [
                rxn for rxn in self.mechanism.reactions
                if (_canon(rxn.reactants), _canon(rxn.products)) not in excl_keys
            ]
            for rxn in self.mechanism.reactions:
                rk = (_canon(rxn.reactants), _canon(rxn.products))
                for key, new_rate in overrides:
                    if rk == key:
                        rxn.rate_expression = new_rate

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

        if 'JO1D' in ctx and 'J1' not in ctx:
            ctx['J1'] = ctx['JO1D']
        elif 'SAPHIR_JO1D' in ctx and 'J1' not in ctx:
            ctx['J1'] = ctx['SAPHIR_JO1D']
        # J4/J7/J3/J11/J12: use plain names (from EASY script), fallback to SAPHIR_ names
        for j_name, plain, saphir in [('J4','JNO2','SAPHIR_JNO2'), ('J7','JHONO','SAPHIR_JHONO'),
                                       ('J3','JH2O2','SAPHIR_JH2O2'), ('J11','JHCHO_R','SAPHIR_JHCHO_R'),
                                       ('J12','JHCHO_M','SAPHIR_JHCHO_M')]:
            src = plain if plain in ctx else (saphir if saphir in ctx else None)
            if src and j_name not in ctx:
                ctx[j_name] = ctx[src]

        T_arr   = ctx.get('TEMP', np.full(n, 298.15))
        k_N2_q  = 2.15e-11 * np.exp(110.0 / T_arr)
        k_O2_q  = 3.30e-11 * np.exp( 55.0 / T_arr)
        k_H2O_q = np.full(n, 2.20e-10)
        k_quench = (k_N2_q  * ctx.get('N2',  0.7809 * ctx['M']) +
                    k_O2_q  * ctx.get('O2',  0.2095 * ctx['M']) +
                    k_H2O_q * ctx.get('H2O', np.zeros(n)))
        hv_arr = ctx.get('HV', np.ones(n))
        j1 = ctx.get('J1', np.zeros(n)) * hv_arr
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

    def _regenerate_elog(self, config_path: str, wanted: str) -> str:
        """The JSON names a JUICE elog that is not on disk (FACSIMILE reference folders: the elog is a
        derived file the GUI writes at run time).  Rebuild it from the .easy next to the JSON into the
        user's temp folder — never into the scenario folder — and return its path ('' if impossible)."""
        try:
            import tempfile
            from elog_generator import generate_both
            stem = os.path.splitext(os.path.basename(config_path))[0]
            easy = os.path.join(self.root, stem + '.easy')
            if not os.path.isfile(easy):
                return ''
            out = os.path.join(tempfile.gettempdir(), 'JUICE_elog', stem)
            os.makedirs(out, exist_ok=True)
            dirs = [self.root, os.path.join(self.root, 'chemistry')] + [d for d in self._include_dirs if d and os.path.isdir(d)]
            _raw, canon = generate_both(easy, search_dirs=dirs, output_dir=out)
            print(f"[BudgetEngine] {os.path.basename(wanted)} not found — rebuilt from {os.path.basename(easy)} in {out}")
            return canon
        except Exception as e:
            print(f"[BudgetEngine] could not rebuild the elog: {e}")
            return ''

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

        T_arr = ctx.get('T', ctx.get('TEMP', ctx.get('SAPHIR_T', ctx.get('SAPHIR_TEMP', np.full(n, 298.15)))))
        HV_arr = ctx.get('HV', ctx.get('SAPHIR_HV', np.ones(n)))
        JNO2_arr = ctx.get('JNO2', ctx.get('SAPHIR_JNO2', np.zeros(n)))

        COSX_arr = np.full(n, 0.5)
        _sz = ctx.get('SZ', ctx.get('SAPHIR_SZ', None))
        if _sz is not None:
            COSX_arr = np.maximum(np.cos(_sz * np.pi / 180.0), 1e-4)
        safe_cosx = np.where(COSX_arr > 0, COSX_arr, 1e-4)

        # O2/N2: respect env_map/config overrides. Atmospheric fractions are only
        # applied if not already provided by the user (mirrors scenario_runner fix).
        ctx.setdefault('O2', 0.2095 * M_arr)
        ctx.setdefault('N2', 0.7809 * M_arr)
        ctx.update({
            'T': T_arr, 'TEMP': T_arr,
            'M': M_arr,
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
            return cs_j * scale * _hv

        ctx['PHTFKT'] = _phtfkt
        ctx['PHTFKT_CS'] = lambda a, b, c, _c=safe_cosx: (
            a * _c ** b * np.exp(np.clip(c / _c, -700, 700))
        )

        h2o = ctx.get('H2O', ctx.get('SAPHIR_H2O', 0.01 * M_arr))
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

    def __init__(self, df, csv_path=None, toml_path=None):
        super().__init__()
        self.df = df.copy()
        self._toml_path = toml_path
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

        # Load s_entries from the mechanism files to classify columns structurally.
        # s_entries = anything defined with a function call in the EASY script.
        self._easy_s_names = set()   # lowercase names of s_entries
        self._easy_c_names = set()   # lowercase names of ODE species
        try:
            _jui_stem = os.path.splitext(csv_path)[0].replace('_cs', '')
            _jui_json = _jui_stem + '.json'
            if os.path.isfile(_jui_json):
                import json as _json
                with open(_jui_json) as _fh:
                    _jcfg = _json.load(_fh)
                from easy_parser import EasyParser as _EP
                _ep = _EP()
                _s_entries, _species = [], []
                _elog = _jcfg.get('elog_path', '')
                _progress("Reading the mechanism", 25, "species are classified from the chemistry the run used (one parse, reused by the budget engine)")
                if _elog and os.path.isfile(_elog):
                    _m = _parse_easy_cached(_ep, _elog)
                    _s_entries.extend(_m.s_entries)
                    _species.extend(_m.species)
                else:
                    for _f in _jcfg.get('mechanism_files', []):
                        if os.path.isfile(_f):
                            _m = _parse_easy_cached(_ep, _f)
                            _s_entries.extend(_m.s_entries)
                            _species.extend(_m.species)
                self._easy_s_names = {e.name.lower() for e in _s_entries
                                      if not getattr(e, 'is_synthetic', False)}
                self._easy_c_names = {sp.lower() for sp in _species}
        except Exception:
            pass

        # Instrument data loaded from NC files: display_key -> {t, v, nc_file, nc_var, mcm_match}
        self._inst_data = {}

        self.setWindowTitle("JUICE — Results" + (f"  ·  {os.path.basename(csv_path)}" if csv_path else ""))
        self.setGeometry(50, 50, 1600, 900)
        # About as standalone top-left menu bar action (no submenu)
        self.menuBar().addAction("About JUICE…", self._show_about)
        try:
            from juice_theme import app_icon
            self.setWindowIcon(app_icon())
        except Exception:
            _icon_pixmap = QPixmap(os.path.join(os.path.dirname(os.path.abspath(__file__)), "JUICE Logo.png"))
            if not _icon_pixmap.isNull():
                self.setWindowIcon(QIcon(_icon_pixmap))

        # Convert MCM index from SAPHIR epoch (sec since 2000-01-01) to hours since
        # UTC midnight of the experiment day, same reference as NC instrument data.
        self.df.index = ((self.df.index.astype(float) + _NC_EPOCH_OFFSET) - self._midnight_unix) / 3600.0

        # Add virtual derived species (in-memory only, not written to CSV)
        self.df = self._add_virtual_species(self.df)

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
        self._mcm_folder    = QTreeWidgetItem(self.var_tree, ["Species (c)"])
        self._saphir_folder = QTreeWidgetItem(self.var_tree, ["Settings (s)"])
        self._inst_folder   = QTreeWidgetItem(self.var_tree, ["SAPHIR Instruments"])

        self._mcm_folder.setExpanded(True)
        self._saphir_folder.setExpanded(True)
        self._inst_folder.setExpanded(True)

        self._populate_tree()

        # Scenario (found automatically next to the result; '…' only to pick another one)
        cfg_row = QHBoxLayout(); cfg_row.setSpacing(4)
        self.cfg_status = QLabel(self._config_label_text)
        self.cfg_status.setWordWrap(True)
        self.cfg_status.setStyleSheet("color: #555; font-size: 10px;")
        self.cfg_status.setToolTip("scenario JSON used for the process analysis, reactivity and budgets — loaded automatically next to the result")
        cfg_row.addWidget(self.cfg_status, stretch=1)
        btn_cfg = QToolButton(); btn_cfg.setText("…"); btn_cfg.setAutoRaise(True); btn_cfg.setToolTip("choose another scenario JSON")
        btn_cfg.clicked.connect(self._browse_config)
        cfg_row.addWidget(btn_cfg)
        left_layout.addLayout(cfg_row)

        btn_exit = QPushButton("✕ Exit")
        btn_exit.setObjectName("primary")
        btn_exit.clicked.connect(self.close)
        left_layout.addWidget(btn_exit)

        # Legacy single-canvas controls: kept alive (plot_variable & co. still work) but no longer shown —
        # the Timeseries tab is a panel grid now (timeseries_tab.TimeseriesTab).
        self._legacy_controls = QWidget(); self._legacy_controls.hide()
        control_layout = QVBoxLayout(self._legacy_controls)
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

        splitter.addWidget(left_panel)

        # ── RIGHT PANEL: tabbed ─────────────────────────────────────────────
        right_panel = QWidget()
        right_layout = QVBoxLayout(right_panel)
        right_layout.setContentsMargins(0, 0, 0, 0)

        self.right_tabs = QTabWidget()
        right_layout.addWidget(self.right_tabs)
        self._init_figure_list()                      # "Add to figure list" / "Figure list (n)" in the tab-bar corner

        # Tab 0 – Timeseries: a grid of linked panels (time series and correlations)
        self.fig = Figure(figsize=(10, 6))                 # legacy canvas, hidden (plot_variable & co.)
        self.canvas = FigureCanvas(self.fig)
        self.toolbar = NavigationToolbar(self.canvas, self._legacy_controls)
        try:
            from timeseries_tab import TimeseriesTab
            self.ts_tab = TimeseriesTab(self)
            self.right_tabs.addTab(self.ts_tab, "Timeseries")
        except Exception as _e:
            print(f"[viewer] Timeseries panels not available ({_e}) — legacy canvas shown")
            self.ts_tab = None
            ts_widget = QWidget(); ts_layout = QVBoxLayout(ts_widget); ts_layout.setContentsMargins(0, 0, 0, 0)
            ts_layout.addWidget(self.toolbar); ts_layout.addWidget(self.canvas)
            self.right_tabs.addTab(ts_widget, "Timeseries")

        _progress("Building the tabs: Timeseries, Process analysis, Reactivity", 55)
        # Tab 1 – Process Analysis
        pa_widget = QWidget()
        pa_layout = QVBoxLayout(pa_widget)
        pa_layout.setContentsMargins(0, 0, 0, 0)
        self.pa_fig = Figure(figsize=(12, 6))
        self.pa_canvas = FigureCanvas(self.pa_fig); _interact(self.pa_canvas)
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
        rx_ctrl.addWidget(QLabel("group into Other below"))
        self.rx_thresh = QDoubleSpinBox(); self.rx_thresh.setRange(0.0, 50.0); self.rx_thresh.setDecimals(1)
        self.rx_thresh.setSingleStep(0.5); self.rx_thresh.setValue(2.0); self.rx_thresh.setSuffix(" %")
        self.rx_thresh.setToolTip("contributors whose time-average share of the total reactivity is below this are summed into 'Other'")
        self.rx_thresh.valueChanged.connect(lambda _v: self._update_reactivity())
        rx_ctrl.addWidget(self.rx_thresh)
        self.rx_meas_label = QLabel("No measured k_OH loaded")
        self.rx_meas_label.setStyleSheet("color: #555; font-size: 10px;")
        self.rx_meas_label.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Preferred)
        rx_ctrl.addWidget(self.rx_meas_label)
        rx_ctrl.addStretch()
        rx_layout.addLayout(rx_ctrl, 0)
        self._measured_koh = None
        self.rx_fig = Figure(figsize=(12, 6))
        self.rx_canvas = FigureCanvas(self.rx_fig); _interact(self.rx_canvas)
        self.rx_toolbar = NavigationToolbar(self.rx_canvas, self)
        rx_layout.addWidget(self.rx_toolbar, 0)
        rx_layout.addWidget(self.rx_canvas, 1)
        self.right_tabs.addTab(rx_widget, "Reactivity")

        # Tab 3 – Carbon Budget
        _progress("Building the tabs: Carbon budget, Chemical space", 65)
        self._build_carbon_tab()

        # Tab 4 – Chemical Space (VBS + molecular distribution)
        self._build_chemspace_tab()

        # Tab 5 – Comparison (last)
        _progress("Building the tabs: Comparison, Overview", 75, "the Overview reads the day's observation file next to the result, or the archive the first time")
        self._build_comparison_tab()

        splitter.addWidget(right_panel)
        splitter.setSizes([300, 1300])

        _progress("Drawing the default plot", 85)
        self.init_axes()
        self._plot_defaults()

        # Auto-discover scenario_config.json — must come after UI is fully built
        _progress("Loading the scenario for the process analysis", 92, "rate constants and reactions of the run are prepared for the budgets")
        self._auto_load_config(csv_path)
        _progress("Opening the window", 99)
        self.showMaximized()

    # ── Virtual / derived species ────────────────────────────────────────────
    def _add_virtual_species(self, df):
        """Add derived columns that are not in the data but can be computed from it.
        These exist only in memory — never written to the output CSV.

        In ELOG mode (when _easy_s_names is populated), virtual species are only
        added if they are NOT defined in the EASY script — otherwise they would
        already be in the data as s_entries or c_entries.
        If they ARE in the EASY script, they come from the cs.csv directly.
        """
        df = df.copy()
        col_upper = {c.upper(): c for c in df.columns}   # upper→actual name map

        def _get(name):
            """Return column array for name (case-insensitive) or None."""
            actual = col_upper.get(name.upper())
            return df[actual].values.astype(float) if actual is not None else None

        # Only add a virtual species if it is NOT defined in the EASY script.
        # If it IS defined in the script it already appears in the cs.csv as an s-entry.
        # When _easy_s_names is empty (no elog loaded) we add all virtual species.
        easy_s = getattr(self, '_easy_s_names', set())

        def _virtual_sum(label, constituents, df_out, col_upper_out):
            """Sum constituents into label if ≥2 found, label not already present,
            and label not defined in the EASY script (those come from cs directly)."""
            if label.upper() in col_upper_out:
                return df_out, col_upper_out          # already in data
            if easy_s and label.lower() in easy_s:
                return df_out, col_upper_out          # defined in EASY script
            total = None
            found = []
            for sp, factor in constituents:
                v = _get(sp)
                if v is not None:
                    total = v * factor if total is None else total + v * factor
                    found.append(f"{factor}×{sp}" if factor != 1 else sp)
            if total is not None and len(found) >= 2:
                df_out = df_out.copy() if id(df_out) == id(df) else df_out
                df_out[label] = total
                col_upper_out = col_upper_out | {label.upper(): label}
            return df_out, col_upper_out

        # NOx = NO + NO2
        df, col_upper = _virtual_sum('NOx',
            [('NO', 1), ('NO2', 1)], df, col_upper)

        # Ox = O3 + NO2
        df, col_upper = _virtual_sum('Ox',
            [('O3', 1), ('NO2', 1)], df, col_upper)

        # NOy = NO + NO2 + NO3 + 2×N2O5 + HNO3 + HONO + HNO4 + PAN
        df, col_upper = _virtual_sum('NOy',
            [('NO', 1), ('NO2', 1), ('NO3', 1), ('N2O5', 2),
             ('HNO3', 1), ('HONO', 1), ('HNO4', 1), ('PAN', 1)],
            df, col_upper)

        return df

    # ── Tree building ─────────────────────────────────────────────────────────

    def _is_setting_col(self, col):
        """True if col is a setting (s-value) from the EASY script.
        Classification is based exclusively on s_entries parsed from the mechanism
        files at load time (_easy_s_names). A JSON file is always present.
        """
        return col.lower() in self._easy_s_names

    def _populate_tree(self):
        """Fill Species (c) and Settings (s) folders from self.df columns."""
        self._mcm_folder.takeChildren()
        self._saphir_folder.takeChildren()
        self._saphir_folder.setText(0, "Settings (s)")

        mcm_active, mcm_empty = [], []
        for col in sorted(self.df.columns):
            if self._is_setting_col(col):
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
        self.right_tabs.setCurrentIndex(self._tab("Process Analysis"))

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

    def _tab(self, name: str) -> int:
        """Index of a right-hand tab by its title (tabs may be inserted in front, e.g. Overview)."""
        for i in range(self.right_tabs.count()):
            if self.right_tabs.tabText(i) == name:
                return i
        return 0

    def _on_tree_double_click(self, item, _column):
        """Double-click: context-aware action based on the active right tab.
        - Timeseries tab  → plot the species on the timeseries canvas.
        - Process Analysis tab → run the budget for this species.
        - Other tabs → no-op.
        """
        if item.childCount() != 0:
            return  # folder node
        active = self.right_tabs.tabText(self.right_tabs.currentIndex())
        if active in ("Timeseries", "Overview"):      # plot it (the Overview is only a landing page)
            self.right_tabs.setCurrentIndex(self._tab("Timeseries"))
            self.add_selected()
        elif active == "Process Analysis":
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
        if getattr(self, "ts_tab", None) is not None:
            self.ts_tab.add_selected()
            self.right_tabs.setCurrentIndex(self._tab("Timeseries"))
            return
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
        self.right_tabs.setCurrentIndex(self._tab("Timeseries"))

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

    def _show_about(self):
        import os
        from PyQt6.QtWidgets import QDialog, QVBoxLayout, QLabel, QDialogButtonBox
        from PyQt6.QtCore import Qt
        from PyQt6.QtGui import QPixmap
        try:
            from version import __version__, __date__, __authors__
        except ImportError:
            __version__, __date__, __authors__ = "unknown", "", []
        dlg = QDialog(self)
        dlg.setWindowTitle("About JUICE")
        dlg.setFixedWidth(380)
        lay = QVBoxLayout(dlg)
        lay.setSpacing(10)
        lay.setContentsMargins(24, 20, 24, 16)
        logo_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "JUICE Logo.png")
        logo_lbl = QLabel()
        px = QPixmap(logo_path)
        if not px.isNull():
            logo_lbl.setPixmap(px.scaledToWidth(120, Qt.TransformationMode.SmoothTransformation))
        logo_lbl.setAlignment(Qt.AlignmentFlag.AlignHCenter)
        lay.addWidget(logo_lbl)
        title = QLabel(f"<h2>JUICE {__version__}</h2>")
        title.setAlignment(Qt.AlignmentFlag.AlignHCenter)
        lay.addWidget(title)
        subtitle = QLabel("Jülich Utility for Integrating Chemical Equations")
        subtitle.setAlignment(Qt.AlignmentFlag.AlignHCenter)
        subtitle.setWordWrap(True)
        lay.addWidget(subtitle)
        lay.addSpacing(6)
        info = QLabel(
            f"<b>Version:</b> {__version__}<br>"
            f"<b>Date:</b> {__date__}<br>"
            f"<b>Authors:</b><br>" +
            "<br>".join(f"&nbsp;&nbsp;{a}" for a in __authors__)
        )
        info.setAlignment(Qt.AlignmentFlag.AlignLeft)
        lay.addWidget(info)
        lay.addSpacing(4)
        bb = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok)
        bb.accepted.connect(dlg.accept)
        lay.addWidget(bb)
        dlg.exec()

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
        csv_path = os.path.abspath(csv_path)
        csv_dir  = os.path.dirname(csv_path)
        csv_stem = os.path.splitext(os.path.basename(csv_path))[0]
        # Strip "_juice" suffix so scene_05_juice.cs matches scene_05.json
        if csv_stem.endswith('_juice'):
            csv_stem = csv_stem[:-len('_juice')]
        candidates = []
        for root in [csv_dir, os.path.dirname(csv_dir)]:
            if not root or not os.path.isdir(root):
                continue
            stem_match = os.path.join(root, f"{csv_stem}.json")
            if os.path.isfile(stem_match):
                candidates.append(stem_match)
            legacy = os.path.join(root, 'scenario_config.json')
            if os.path.isfile(legacy):
                candidates.append(legacy)
        if candidates:
            self._load_config(candidates[0], silent=True)

    # ── figure list (report_basket / report_dock) ─────────────────────────────
    def _init_figure_list(self):
        """Corner buttons of the tab bar + the dock; the list lives next to the result file."""
        self.figure_list = None; self.figure_dock = None
        try:
            from report_basket import FigureList
            from report_dock import FigureListDock
            if getattr(self, "csv_path", ""):
                self.figure_list = FigureList.for_result(self.csv_path)
            else:
                import tempfile
                self.figure_list = FigureList(os.path.join(tempfile.gettempdir(), "juice_figures"), "results")
            corner = QWidget(); cl = QHBoxLayout(corner); cl.setContentsMargins(0, 0, 4, 0); cl.setSpacing(4)
            add = QToolButton(); add.setText("➕ Add to figure list"); add.setAutoRaise(True)
            add.setToolTip("snapshot the figures of this tab (SVG + PNG) into the figure list for the document")
            add.clicked.connect(self.add_current_to_figure_list); cl.addWidget(add)
            self._figlist_btn = QToolButton(); self._figlist_btn.setText("Figure list (0)"); self._figlist_btn.setAutoRaise(True)
            self._figlist_btn.setToolTip("show the collected figures; reorder, caption, create the document")
            self._figlist_btn.clicked.connect(self._toggle_figure_dock); cl.addWidget(self._figlist_btn)
            self.right_tabs.setCornerWidget(corner, Qt.Corner.TopRightCorner)
            self.figure_dock = FigureListDock(self, self.figure_list, self)
            self.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, self.figure_dock)
            self.figure_dock.hide()
        except Exception as e:
            print(f"[viewer] figure list not available: {e}")

    def _toggle_figure_dock(self):
        if self.figure_dock is not None:
            self.figure_dock.setVisible(not self.figure_dock.isVisible())

    def _report_figures(self):
        """[(title, Figure)] of the tab currently shown — what 'Add to figure list' snapshots."""
        name = self.right_tabs.tabText(self.right_tabs.currentIndex()); tab = self.right_tabs.currentWidget()
        figs = []; self._report_grid = None; self._report_panel = (6.4, 3.8)
        grid = getattr(tab, "grid", None)
        if grid is not None:                                           # Overview / Timeseries panel grids
            for i in range(grid.count()):
                w = grid.itemAt(i).widget()
                if w is not None and getattr(w, "fig", None) is not None:
                    title = getattr(w, "title", "") or f"{name} panel {i + 1}"
                    figs.append((title, w.fig))
            try:
                self._report_grid = (grid.rowCount(), grid.columnCount())      # same arrangement as on screen
                if grid.columnCount() == 1:
                    self._report_panel = (10.0, 3.6)                            # stacked panels: wide
            except Exception:
                pass
        elif name == "Process Analysis":
            figs.append((self.pa_fig.axes[0].get_title() if self.pa_fig.axes else "Process analysis", self.pa_fig))
        elif name == "Reactivity":
            figs.append((f"{self.rx_combo.currentText()} reactivity", self.rx_fig))
        elif name == "Budgets":
            self._report_grid = (2, 1); self._report_panel = (11.0, 4.0)
            figs.append(("Carbon budget", self._cb_fig))
            figs.append(((self._nb_view.currentText() if hasattr(self, "_nb_view") else "Nitrogen") + " budget", self._nb_fig))
        elif name == "Chemical Space":
            self._report_grid = (2, 2); self._report_panel = (6.4, 4.4)
            for title, fig in (("Injections and roof status", self._cs_ts_fig), ("Volatility basis set", self._cs_vbs_fig),
                               ("Composition by carbon number", self._cs_co_fig), ("Mass defect", self._cs_fam_fig)):
                figs.append((title, fig))
        elif name == "Comparison":
            figs.append(("Comparison slopes", self._cmp_slope_fig))
            for i, f in enumerate(getattr(self, "_cmp_ts_figs", [])):
                figs.append((f"Comparison timeseries {i + 1}", f))
        return [(t, f) for t, f in figs if f is not None and f.axes]

    def add_current_to_figure_list(self):
        if self.figure_dock is None:
            return
        name = self.right_tabs.tabText(self.right_tabs.currentIndex())
        figs = self._report_figures()
        if not figs:
            self.statusBar().showMessage("nothing drawn on this tab yet", 4000); return
        if len(figs) > 1:                       # several figures on the tab → one figure with sub-panels, as on screen
            try:
                cap = f"{self.figure_list.run_name}: {name} — " + "; ".join(t for t, _f in figs)
                self.figure_list.add_composite([f for _t, f in figs], name, source=name, caption=cap, grid=getattr(self, "_report_grid", None),
                                               panel_size=getattr(self, "_report_panel", (6.4, 3.8)))
                self.figure_dock.refresh(); n = 1
            except Exception as e:
                print(f"[figure list] composite failed ({e}); adding the panels separately")
                n = self.figure_dock.add_figures(figs, source=name)
        else:
            n = self.figure_dock.add_figures(figs, source=name)
        self.statusBar().showMessage(f"{n} figure(s) from {name} added to the list ({len(self.figure_list.items)} in total)", 5000)
        if not self.figure_dock.isVisible():
            self.figure_dock.show()

    def _toml_include_dirs(self):
        """[paths] include_dir of the active TOML (where the EASY #include files live), if any."""
        try:
            import tomllib as _tl
            p = getattr(self, '_toml_path', None)
            if not (p and os.path.isfile(p)):                      # started without --toml: the project TOMLs
                here = os.path.dirname(os.path.abspath(__file__))
                p = next((q for q in (os.path.join(here, 'GG_GUI.toml'), os.path.join(here, 'FR_GUI.toml')) if os.path.isfile(q)), None)
            if p and os.path.isfile(p):
                with open(p, 'rb') as fh:
                    t = _tl.load(fh)
                d = (t.get('paths', {}) or {}).get('include_dir') or t.get('include_dir')
                return [d] if d else []
        except Exception:
            pass
        return []

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
            self.budget_engine.load(path, include_dirs=self._toml_include_dirs())
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
                if vi > 0 and not in_open:
                    in_open, t0 = True, ti
                elif vi <= 0 and in_open:
                    ax_ts.axvspan(t0, ti, color='#fffacd', alpha=0.25, zorder=0, linewidth=0)
                    in_open = False
            if in_open:
                ax_ts.axvspan(t0, t[-1], color='#fffacd', alpha=0.25, zorder=0, linewidth=0)

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
                         borderpad=0.5, labelspacing=0.3).set_draggable(True)

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
            leg.set_draggable(True)

    # ── Reactivity tab ────────────────────────────────────────────────────────
    def _auto_load_observations(self):
        """The day's SAPHIR observations (observation file next to the result, else the archive) as
        'SAPHIR Instruments' entries of the tree — what the upload button used to do by hand."""
        try:
            import experiment_plots as _ep, datetime as _dt
            spec = _ep.spec_for_run(self.csv_path) if getattr(self, "csv_path", "") else None
            root = spec.conditions.archive_root if spec else r"Z:\IEK8-SAPHIR\data"
            date = (spec.conditions.date if spec and spec.conditions.date
                    else _dt.datetime.fromtimestamp(self._midnight_unix, _dt.timezone.utc).strftime("%Y-%m-%d"))
            obs = _ep.load_observations(root, date, None, self._midnight_unix,
                                        cache_dir=os.path.dirname(os.path.abspath(self.csv_path)) if getattr(self, "csv_path", "") else "")
            rows = []
            for key, (h, v, label, kind) in obs.items():
                suffix = _ep.OBS_MAP[key][0] if key in _ep.OBS_MAP else "observations"
                rows.append(({'nc_file': f"{date}.{suffix}", 'nc_var': key, 't': np.asarray(h, float), 'v': np.asarray(v, float)},
                             key if key in self.df.columns else ''))
            if rows:
                self._add_confirmed_rows(rows)
        except Exception as e:
            print(f"[viewer] observations not added to the tree: {e}")

    def _auto_load_koh(self):
        """At start-up: the measured k_OH of the experiment day straight from the archive, if it exists."""
        try:
            import datetime as _dt
            import experiment_plots as _ep
            spec = _ep.spec_for_run(self.csv_path) if getattr(self, "csv_path", "") else None
            root = spec.conditions.archive_root if spec else r"Z:\IEK8-SAPHIR\data"
            date = (spec.conditions.date if spec and spec.conditions.date
                    else _dt.datetime.fromtimestamp(self._midnight_unix, _dt.timezone.utc).strftime("%Y-%m-%d"))
            note: list = []
            obs = _ep.load_observations(root, date, ["kOH"], self._midnight_unix,
                                        cache_dir=os.path.dirname(os.path.abspath(self.csv_path)) if getattr(self, "csv_path", "") else "",
                                        status=note)
            if "kOH" in obs:
                h, k, label, _kind = obs["kOH"]
                self._measured_koh = {'t': np.asarray(h, dtype=float), 'k': np.asarray(k, dtype=float)}
                self.rx_meas_label.setText(f"k_OH: {label}, {note[0] if note else 'archive'}")
                self.rx_meas_label.setStyleSheet("color: #1a7a1a; font-size: 10px; font-weight: bold;")
                for name in ("_update_reactivity", "_plot_reactivity", "_refresh_reactivity", "_rx_update"):
                    fn = getattr(self, name, None)
                    if callable(fn):
                        try:
                            fn()
                        except Exception:
                            pass
                        break
        except Exception as e:
            print(f"[Reactivity] auto-load skipped: {e}")

    def _load_koh_nc(self, path: str = ""):
        """Load measured k_OH from a SAPHIR NC file (dialog when no path is given)."""
        if not HAS_NC4:
            self.rx_meas_label.setText("netCDF4 not available")
            self.rx_meas_label.setStyleSheet("color: #cc0000; font-size: 10px;")
            return

        if not path:
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
        THRESH = float(self.rx_thresh.value()) / 100.0 if hasattr(self, "rx_thresh") else 0.02

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
                if vi > 0 and not in_open:
                    in_open, t0 = True, ti
                elif vi <= 0 and in_open:
                    ax_ts.axvspan(t0, ti, color='#fffacd', alpha=0.25, zorder=0, linewidth=0)
                    in_open = False
            if in_open:
                ax_ts.axvspan(t0, t[-1], color='#fffacd', alpha=0.25, zorder=0, linewidth=0)

        labels  = list(major.keys())
        arrays  = [major[k] for k in labels]
        colors  = _distinct_colors(labels)
        if arrays:
            ax_ts.stackplot(t, arrays, labels=labels, colors=colors, alpha=0.85)

        k_total = sum(contribs.values()) if contribs else np.zeros(n)
        ax_ts.plot(t, k_total, color='black', lw=1.5, ls='--',
                   label=f'Total k_{radical} (model)', zorder=5)

        # kOH from pd file: direct kOH column (Σ k_i×[X_i], no division by [OH])
        if radical.upper() == 'OH' and self.csv_path:
            import os as _os
            _pd_path = (self.csv_path.replace('_cs.csv', '_pd.csv')
                        if self.csv_path.endswith('_cs.csv')
                        else os.path.splitext(self.csv_path)[0] + '.pd')
            if _os.path.isfile(_pd_path):
                try:
                    _pd_df = load_column_file(_pd_path)
                    if 'kOH' in _pd_df.columns:
                        _koh_pd = np.interp(t, _pd_df.index.values.astype(float),
                                            _pd_df['kOH'].values.astype(float))
                    elif 'D_OH' in _pd_df.columns and 'OH' in self.df.columns:
                        _d_oh = np.interp(t, _pd_df.index.values.astype(float),
                                          _pd_df['D_OH'].values.astype(float))
                        _oh = self.df['OH'].values.astype(float)
                        with np.errstate(divide='ignore', invalid='ignore'):
                            _koh_pd = np.where(_oh > 1e2, _d_oh / _oh, np.nan)
                    else:
                        _koh_pd = None
                    if _koh_pd is not None:
                        ax_ts.plot(t, _koh_pd, color='#2ca02c', lw=1.8, ls='-.',
                                   label='kOH JUICE (Σk·[X])', zorder=5)
                except Exception:
                    pass

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
                         bbox_to_anchor=(0.5, -0.14), ncol=min(6, max(3, (len(handles) + 3) // 4)), fontsize=7, framealpha=0.7).set_draggable(True)

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
        self._cmp_load_btn = QPushButton("Load Reference…")
        self._cmp_load_btn.clicked.connect(self._cmp_load_reference)
        src_row.addWidget(self._cmp_load_btn)

        src_row.addSpacing(12)
        src_row.addWidget(QLabel("Noise threshold:"))
        self._cmp_noise_input = QLineEdit("1e-6")
        self._cmp_noise_input.setMaximumWidth(70)
        src_row.addWidget(self._cmp_noise_input)
        btn_slopes = QPushButton("Compute Slopes")
        btn_slopes.clicked.connect(self._cmp_compute_slopes)
        src_row.addWidget(btn_slopes)
        src_row.addStretch()
        outer_layout.addLayout(src_row)              # row 1: reference + slopes
        src_row = QHBoxLayout()                      # row 2: what to compare + status (one row was 2000 px wide)

        # 3-state species selector: Favoriten | c-Spezies | Alle
        from PyQt6.QtWidgets import QButtonGroup
        _btn_style = ("QPushButton { border: 1px solid #888; padding: 3px 8px; }"
                      "QPushButton:checked { background-color: #2E5A88; color: white;"
                      "  font-weight: bold; border: 1px solid #1a3f66; }")
        self._cmp_mode_group = QButtonGroup(self)
        self._cmp_mode_group.setExclusive(True)

        # Load favorites from TOML (optional section)
        _toml_favs = []
        try:
            import tomllib as _tl, os as _os
            _tp = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "FR_GUI.toml")
            with open(_tp, 'rb') as _fh:
                _tcfg = _tl.load(_fh)
            _toml_favs = _tcfg.get('comparison', {}).get('favorites', [])
            # Handle legacy string format e.g. "['OH','HO2',...]"
            if isinstance(_toml_favs, str):
                import ast as _ast
                try: _toml_favs = _ast.literal_eval(_toml_favs)
                except Exception: _toml_favs = []
        except Exception:
            pass
        self._cmp_favorites = [str(s).strip() for s in _toml_favs
                               if isinstance(s, str) and len(s) > 1]

        # 4 exclusive mode buttons (All removed)
        _defs = [
            (0, "★ Favorites", "Parameters listed in TOML [comparison] favorites",      bool(self._cmp_favorites)),
            (1, "c-Species",   "Concentrations: ODE species + derived (NOx, NOy, Ox, kOH)", not bool(self._cmp_favorites)),
            (2, "s-Species",   "Settings: T, P, HV, SZ, RH, DILUTE, …",                False),
            (3, "k-Values",    "Reaction turnover rates from kv files",                  False),
        ]
        for _id, _label, _tip, _checked in _defs:
            _b = QPushButton(_label)
            _b.setCheckable(True)
            _b.setChecked(_checked)
            _b.setToolTip(_tip)
            _b.setStyleSheet(_btn_style)
            if _id == 0:
                _b.setEnabled(bool(self._cmp_favorites))
                self._cmp_btn_fav = _b
            elif _id == 1:
                self._cmp_btn_c = _b
            elif _id == 2:
                self._cmp_btn_s = _b
            else:
                self._cmp_btn_k = _b
            self._cmp_mode_group.addButton(_b, _id)
            src_row.addWidget(_b)

        # Keep legacy attribute so existing code referencing _cmp_react_only_btn still works
        self._cmp_react_only_btn = self._cmp_btn_c

        # Re-run slope computation whenever the mode changes (if reference is loaded)
        self._cmp_mode_group.idClicked.connect(
            lambda _id: self._cmp_compute_slopes()
            if getattr(self, '_cmp_df_ref', None) is not None else None)

        self._cmp_ref_lbl = QLabel("No reference loaded")
        self._cmp_ref_lbl.setStyleSheet("color: #666; font-style: italic;")
        self._cmp_ref_lbl.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Preferred)
        src_row.addWidget(self._cmp_ref_lbl, stretch=1)
        outer_layout.addLayout(src_row)

        # ── Splitter: 1/3 matrix | 2/3 detail ────────────────────────────────
        splitter = QSplitter(Qt.Orientation.Horizontal)

        left = QWidget()
        ll = QVBoxLayout(left)
        ll.setContentsMargins(0, 0, 0, 0)
        ll.setSpacing(2)

        # ── Filter bar (above toolbar, above canvas) ──────────────────────────
        filter_row = QHBoxLayout()
        filter_row.setSpacing(4)

        filter_row.addWidget(QLabel("Filter:"))
        self._cmp_filter_combo = QComboBox()
        self._cmp_filter_combo.addItems([
            "Show all",
            "Best 50  (slope closest to 1)",
            "Best 25  (slope closest to 1)",
            "Worst 50  (slope furthest from 1)",
            "Worst 25  (slope furthest from 1)",
        ])
        self._cmp_filter_combo.setMinimumWidth(190)
        self._cmp_filter_combo.currentIndexChanged.connect(self._cmp_apply_filter)
        filter_row.addWidget(self._cmp_filter_combo)

        filter_row.addSpacing(8)
        filter_row.addWidget(QLabel("Sort:"))
        self._cmp_sort_combo = QComboBox()
        self._cmp_sort_combo.addItems([
            "Alphabetical",
            "Time of peak  (earliest first)",
            "Time of peak  (latest first)",
            "Carbon count  (low → high)",
            "Carbon count  (high → low)",
            "Nitrogen count  (N≥1, low → high)",
            "Radicals first",
        ])
        self._cmp_sort_combo.setMinimumWidth(200)
        self._cmp_sort_combo.currentIndexChanged.connect(self._cmp_apply_filter)
        filter_row.addWidget(self._cmp_sort_combo)

        filter_row.addSpacing(8)
        filter_row.addWidget(QLabel("Search:"))
        self._cmp_search_box = QLineEdit()
        self._cmp_search_box.setPlaceholderText("type species name…")
        self._cmp_search_box.setMaximumWidth(160)
        self._cmp_search_box.textChanged.connect(self._cmp_apply_filter)
        filter_row.addWidget(self._cmp_search_box)

        filter_row.addSpacing(8)
        filter_row.addWidget(QLabel("Show:"))
        self._cmp_limit_spin = QSpinBox()
        self._cmp_limit_spin.setRange(1, 9999)
        self._cmp_limit_spin.setValue(9999)
        self._cmp_limit_spin.setSpecialValueText("All")   # shows "All" when at minimum... use max trick below
        self._cmp_limit_spin.setSuffix("  boxes")
        self._cmp_limit_spin.setMaximumWidth(105)
        self._cmp_limit_spin.setToolTip("Maximum number of boxes to display after filter/sort/search.\nSet to 9999 to show all.")
        self._cmp_limit_spin.valueChanged.connect(self._cmp_apply_filter)
        filter_row.addWidget(self._cmp_limit_spin)
        filter_row.addStretch(1)

        ll.addLayout(filter_row)

        self._cmp_slope_fig = Figure()
        self._cmp_slope_canvas = FigureCanvas(self._cmp_slope_fig)
        self._cmp_slope_canvas.mpl_connect('button_press_event', self._cmp_on_slope_click)
        ll.addWidget(NavigationToolbar(self._cmp_slope_canvas, self))
        ll.addWidget(self._cmp_slope_canvas, stretch=1)
        splitter.addWidget(left)

        # ── Right panel: two stacked timeseries slots ─────────────────────────
        right = QWidget()
        rl = QVBoxLayout(right)
        rl.setContentsMargins(0, 0, 0, 0)

        right_vsplit = QSplitter(Qt.Orientation.Vertical)

        self._cmp_ts_figs      = []
        self._cmp_ts_canvases  = []
        self._cmp_ts_combos    = []
        self._cmp_ts_ppb_btns  = []
        self._cmp_ts_diff_btns = []

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

            ppb_btn = QPushButton("ppb")
            ppb_btn.setCheckable(True)
            ppb_btn.setChecked(True)
            ppb_btn.setFixedWidth(44)
            ppb_btn.setToolTip("Toggle between molec/cm³ and ppb (÷M×10⁹)")
            ppb_btn.toggled.connect(
                lambda checked, s=s: self._cmp_draw_timeseries(
                    s, self._cmp_ts_combos[s].currentText()))
            self._cmp_ts_ppb_btns.append(ppb_btn)
            ctrl.addWidget(ppb_btn)

            diff_btn = QPushButton("Δ")
            diff_btn.setCheckable(True)
            diff_btn.setFixedWidth(36)
            diff_btn.setToolTip(
                "Show difference JUICE − Reference on right axis\n"
                "(left axis: both timeseries, right axis: Δ)")
            diff_btn.toggled.connect(
                lambda checked, s=s: self._cmp_draw_timeseries(
                    s, self._cmp_ts_combos[s].currentText()))
            self._cmp_ts_diff_btns.append(diff_btn)
            ctrl.addWidget(diff_btn)
            slot_layout.addLayout(ctrl)

            fig = Figure()
            canvas = FigureCanvas(fig); _interact(canvas, hover=False)      # has its own crosshair
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
        # Landing page: built from the experiment description + archive observations
        try:
            from overview_tab import OverviewTab
            self.right_tabs.insertTab(0, OverviewTab(self), "Overview")
            self.right_tabs.setCurrentIndex(0)
        except Exception as _e:
            print(f"[viewer] Overview tab not available: {_e}")
        self._auto_load_koh()
        self._auto_load_observations()

    # ── Reference loading ─────────────────────────────────────────────────────
    def _cmp_easy_scenario_dir(self) -> str:
        """Return the EASY scenario subfolder for the current simulation.

        Derives the scenario name from the loaded CSV path (e.g. scene_01 from
        …/JUICE/scene_01/simulation_results.csv), reads easy_root from
        FR_GUI.toml, and returns easy_root/scene_01 if it exists and contains
        at least one .cs file.  Falls back to "" so the file dialog opens in
        its default location.
        """
        try:
            import tomllib
            toml_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                     "FR_GUI.toml")
            if os.path.isfile(toml_path):
                with open(toml_path, "rb") as fh:
                    toml = tomllib.load(fh)
                easy_root = toml.get("paths", {}).get("easy_root", "")
            else:
                easy_root = ""
        except Exception:
            easy_root = ""

        scenario_name = ""
        if hasattr(self, "csv_path") and self.csv_path:
            scenario_name = os.path.basename(os.path.dirname(
                os.path.abspath(self.csv_path)))

        if easy_root and scenario_name:
            candidate = os.path.join(easy_root, scenario_name)
            if os.path.isdir(candidate) and any(
                    f.endswith(".cs") for f in os.listdir(candidate)):
                return candidate

        return easy_root if easy_root else ""

    def _cmp_load_reference(self):
        src = self._cmp_src_combo.currentText()

        if src.startswith("FACSIMILE"):
            # Auto-load if exactly one .cs file exists in the EASY scenario folder.
            start_dir = self._cmp_easy_scenario_dir()
            auto_path = None
            if start_dir and os.path.isdir(start_dir):
                cs_files = [os.path.join(start_dir, f)
                            for f in os.listdir(start_dir) if f.endswith(".cs")]
                if len(cs_files) == 1:
                    auto_path = cs_files[0]
            if auto_path:
                path = auto_path
            else:
                path, _ = QFileDialog.getOpenFileName(
                    self, "Load FACSIMILE Output", start_dir,
                    "CS Files (*.cs);;All Files (*)")
            if not path:
                return
            try:
                df = self._cmp_parse_facsimile(path)
                self._cmp_df_ref    = self._add_virtual_species(df)
                self._cmp_ref_label = os.path.basename(path)
                self._cmp_kv_ref        = None
                self._cmp_kv_jui        = None   # reset so JUICE kv reloads fresh
                self._cmp_jui_aux_loaded = False  # allow lazy re-load of JUICE aux
                self._cmp_load_aux(path, ref=True)   # load kv/pd alongside
            except Exception as e:
                QMessageBox.warning(self, "Load Error", str(e))
                return

        elif src.startswith("Previous simulation"):
            path, _ = QFileDialog.getOpenFileName(
                self, "Load Previous Simulation CSV", "", "CSV Files (*.csv)")
            if not path:
                return
            try:
                df = load_column_file(path)
                # Same SAPHIR-epoch → hours-since-midnight conversion as main CSV
                t_unix = df.index.astype(float) + _NC_EPOCH_OFFSET
                df.index = (t_unix - self._midnight_unix) / 3600.0
                self._cmp_df_ref    = self._add_virtual_species(df)
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
            self._cmp_df_ref    = self._add_virtual_species(df)
            self._cmp_ref_label = "Loaded SAPHIR instruments"

        self._cmp_ref_lbl.setText(f"Ref: {self._cmp_ref_label}")
        self._cmp_ref_lbl.setStyleSheet("color: black; font-style: normal;")
        self._cmp_load_btn.setStyleSheet(
            "font-weight: bold; background-color: #1e7e34; color: white;")

    # ── Derived quantities and auxiliary data ─────────────────────────────────

    def _cmp_add_derived(self, df, pd_df=None):
        """Add NOx, NOy, Ox, kOH to a cs DataFrame in-place."""
        c = df.columns
        if 'NO' in c and 'NO2' in c:
            df['NOx'] = df['NO'] + df['NO2']
        if 'O3' in c and 'NO2' in c:
            df['Ox'] = df['O3'] + df['NO2']
        if 'NO' in c and 'NO2' in c:
            noy = df['NO'] + df['NO2']
            for sp, coeff in [('HNO3',1),('HNO4',1),('HONO',1),('NO3',1),
                               ('N2O5',2),('PAN',1),('PHAN',1),('CH3NO3',1)]:
                if sp in c:
                    noy = noy + coeff * df[sp]
            df['NOy'] = noy
        # kOH from pd file
        if pd_df is not None and 'kOH' in pd_df.columns:
            df['kOH'] = np.interp(df.index.values.astype(float),
                                  pd_df.index.values.astype(float),
                                  pd_df['kOH'].values.astype(float))
        elif pd_df is not None and 'D_OH' in pd_df.columns and 'OH' in c:
            d_oh = np.interp(df.index.values.astype(float),
                             pd_df.index.values.astype(float),
                             pd_df['D_OH'].values.astype(float))
            df['kOH'] = d_oh / (df['OH'].values + 1e-20)
        return df

    def _cmp_load_aux(self, cs_path, ref=False):
        """Load kv and pd files alongside a cs file; augment the stored cs DataFrame."""
        import sys as _sys
        _sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        try:
            from compare_juice_facsimile import (load_facsimile_kv,
                                                 load_facsimile_pd,
                                                 _trim_boundary_rows)
        except ImportError:
            return

        base = os.path.splitext(cs_path)[0]
        # FACSIMILE .cs: ends in .cs but NOT a JUICE output file
        is_fac = cs_path.endswith('.cs') and not _is_juice_output(cs_path)

        # ── pd ────────────────────────────────────────────────────────────────
        if cs_path.endswith('_cs.csv'):         # old JUICE format (backward compat)
            pd_path = base.replace('_cs', '_pd') + '.csv'
        else:
            pd_path = base + '.pd'
        pd_df = None
        if os.path.isfile(pd_path):
            try:
                if is_fac:
                    pd_df = load_facsimile_pd(pd_path)
                    # align to hours-since-midnight like the cs data
                    pd_df.index = (pd_df.index.values.astype(float)
                                   + _NC_EPOCH_OFFSET - self._midnight_unix) / 3600.0
                else:
                    pd_df = load_column_file(pd_path)
                    pd_df.index = (pd_df.index.values.astype(float)
                                   + _NC_EPOCH_OFFSET - self._midnight_unix) / 3600.0
            except Exception:
                pd_df = None

        # ── augment cs DataFrame ──────────────────────────────────────────────
        if ref:
            if self._cmp_df_ref is not None:
                self._cmp_add_derived(self._cmp_df_ref, pd_df)
        else:
            self._cmp_add_derived(self.df, pd_df)

        # ── kv ────────────────────────────────────────────────────────────────
        if cs_path.endswith('_cs.csv'):          # old JUICE format (backward compat)
            kv_path = base.replace('_cs', '_kv') + '.csv'
        else:
            kv_path = base + '.kv'
        if not os.path.isfile(kv_path):
            return
        try:
            if is_fac:
                k_df, _v = load_facsimile_kv(kv_path)
                # use k (rate constants); align time
                k_df.index = (k_df.index.values.astype(float)
                              + _NC_EPOCH_OFFSET - self._midnight_unix) / 3600.0
                kv_out = k_df
            else:
                raw = load_column_file(kv_path)
                raw.index = (raw.index.values.astype(float)
                             + _NC_EPOCH_OFFSET - self._midnight_unix) / 3600.0
                # keep only k_ columns (rate constants); strip k_ prefix and R#####: tag
                from compare_juice_facsimile import _canon_kv_label as _ckl
                k_map = {}
                for c in raw.columns:
                    if not c.startswith('k_'):
                        continue
                    k_map[c] = _ckl(c)   # canonical label matching FACSIMILE
                kv_out = raw[list(k_map)].rename(columns=k_map)
            if ref:
                self._cmp_kv_ref = kv_out
            else:
                self._cmp_kv_jui = kv_out
        except Exception as _e:
            import traceback
            print(f"[kv load ERROR] {kv_path}:\n{traceback.format_exc()}")

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
        try:
            from compare_juice_facsimile import _trim_boundary_rows
            df = _trim_boundary_rows(df)
        except ImportError:
            pass
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

        # Determine mode (0=Favorites, 1=c-Species, 2=s-Species, 3=k-Values)
        _mode = 1
        if hasattr(self, '_cmp_mode_group'):
            _mode = self._cmp_mode_group.checkedId()
            if _mode < 0:
                _mode = 1

        # Lazily load JUICE aux data (kv/pd/derived) on first slope compute
        if not getattr(self, '_cmp_jui_aux_loaded', False) and self.csv_path:
            self._cmp_kv_jui = None
            self._cmp_load_aux(self.csv_path, ref=False)
            self._cmp_jui_aux_loaded = True

        # For k-Values mode: use kv turnover-rate DataFrames
        if _mode == 3:
            jui_kv = getattr(self, '_cmp_kv_jui', None)
            ref_kv = getattr(self, '_cmp_kv_ref', None)
            if jui_kv is None or ref_kv is None:
                QMessageBox.information(self, "k-Values",
                    "No kv files found. Load a FACSIMILE reference whose .kv file\n"
                    "sits in the same folder as the .cs file.")
                return
            mcm_cols = set(jui_kv.columns)
            ref_cols = set(ref_kv.columns)
        else:
            mcm_cols = set(self.df.columns)
            ref_cols = set(self._cmp_df_ref.columns)

        if _mode == 0:
            # Favoriten: parameters listed in TOML, intersected with available columns
            favs = getattr(self, '_cmp_favorites', [])
            fav_upper = {f.upper(): f for f in favs}
            # match case-insensitively against available columns
            common = sorted(
                c for c in mcm_cols & ref_cols
                if c.upper() in fav_upper
            )
            # preserve TOML order where possible
            _order = {v.upper(): i for i, v in enumerate(favs)}
            common = sorted(common, key=lambda c: _order.get(c.upper(), 9999))
            self._cmp_only_juice = []
            self._cmp_only_ref   = []
        elif _mode == 1:
            # c-Species: ODE concentrations — everything that is not a setting
            common = sorted(c for c in mcm_cols & ref_cols if not self._is_setting_col(c))
            self._cmp_only_juice = []
            self._cmp_only_ref   = []
        elif _mode == 2:
            # s-Species: use elog-aware _is_setting_col (falls back to heuristic when no elog)
            common = sorted(c for c in mcm_cols & ref_cols if self._is_setting_col(c))
            self._cmp_only_juice = []
            self._cmp_only_ref   = []
        elif _mode == 3:
            # k-Values: jui_kv/ref_kv already contain only rate constants (k_ columns),
            # canonically labelled as reaction strings — just intersect, no extra filter.
            common = sorted(mcm_cols & ref_cols)
            self._cmp_only_juice = []
            self._cmp_only_ref   = []

        if not common and not getattr(self, '_cmp_only_juice', []) and not getattr(self, '_cmp_only_ref', []):
            QMessageBox.warning(self, "No Species", "No matching parameters found.")
            return

        # Choose data source: kv DataFrames for mode 3, cs DataFrames otherwise
        if _mode == 3:
            _df_jui = self._cmp_kv_jui
            _df_ref = self._cmp_kv_ref
        else:
            _df_jui = self.df
            _df_ref = self._cmp_df_ref

        t_min = max(_df_jui.index.min(), _df_ref.index.min())
        t_max = min(_df_jui.index.max(), _df_ref.index.max())
        common_t = _df_jui.index.values.astype(float)
        mask = (common_t >= t_min) & (common_t <= t_max)
        common_t = common_t[mask]

        def interp(df, col):
            s = df[col].dropna()
            return np.interp(common_t, s.index.values.astype(float), s.values)

        # Global noise scale: 99th percentile across chemical-species common columns.
        # Physical/environmental parameters (M, T, P, DIL, DILUTE, SAPHIR_*) are
        # excluded: M (~2.5e19) or varying T/P would inflate the scale and cause
        # real chemical species to be filtered out as noise.
        _NON_SCALE = {'M', 'T', 'TEMP', 'P', 'PRESS', 'DIL', 'DILUTE', 'MR_X', 'MR_Y', 'Q', 'HV', 'SZ'}
        scale_cols = [c for c in common
                      if c.upper() not in _NON_SCALE
                      and not c.upper().startswith('SAPHIR_')
                      and not c.upper().startswith('MR_')]
        if not scale_cols:
            scale_cols = common
        _sc_jui = [c for c in scale_cols if c in _df_jui.columns]
        _sc_ref = [c for c in scale_cols if c in _df_ref.columns]
        py_scale = float(np.nanpercentile(
            np.abs(_df_jui[_sc_jui].select_dtypes(include='number').values), 99)) if _sc_jui else 1.0
        ref_scale = float(np.nanpercentile(
            np.abs(_df_ref[_sc_ref].select_dtypes(include='number').values), 99)) if _sc_ref else 1.0
        py_scale  = py_scale  if py_scale  > 0 else 1.0
        ref_scale = ref_scale if ref_scale > 0 else 1.0

        self._cmp_slope_results = {}
        for sp in common:
            if sp not in _df_jui.columns or sp not in _df_ref.columns:
                continue
            py_s  = interp(_df_jui, sp)
            ref_s = interp(_df_ref, sp)

            mask = np.isfinite(py_s) & np.isfinite(ref_s)
            if mask.sum() < 5:
                continue
            px, rx = py_s[mask], ref_s[mask]

            # Skip only if both are identically zero throughout
            if np.nanmax(np.abs(px)) == 0 and np.nanmax(np.abs(rx)) == 0:
                continue

            # Slope filter depends on the quantity type:
            #   c-Species (mode 1): filter by JUICE >= atol=1e4.
            #     The ODE cannot reliably resolve concentrations below its
            #     absolute tolerance (e.g. OH at night → JUICE=0, FAC=1e-13
            #     GEAR noise → ratio=0 biases the median).
            #   s-Settings (mode 2) and k-Values (mode 3): these have units
            #     unrelated to molec/cm³ (K, s⁻¹, cm³ molec⁻¹ s⁻¹, …).
            #     atol=1e4 is meaningless here; use rx > 0 (physically, settings
            #     and rate constants are always positive when chemically active).
            _atol = 1e4  # ODE absolute tolerance — only valid for c-species
            if _mode == 1:  # c-Species: concentrations in molec/cm³
                ratio_mask = (px >= _atol) & (rx > 0) & np.isfinite(px / np.where(rx > 0, rx, np.nan))
            else:           # s-Settings, k-Values: use reference > 0
                ratio_mask = (rx > 0) & np.isfinite(px / np.where(rx > 0, rx, np.nan))
            if ratio_mask.sum() >= 3:
                slope = float(np.median(px[ratio_mask] / rx[ratio_mask]))
            elif np.nanmean(np.abs(rx)) > 0:
                slope = float(np.nanmean(px) / np.nanmean(rx))
            else:
                slope = np.nan

            if not np.isfinite(slope):
                continue

            # σ = √(median((F-J)²)) / mean([F,J])  — median-based relative deviation.
            # Using the median of squared differences (instead of the mean) makes
            # this robust against outlier time steps (e.g. phase transitions where
            # one model drops faster). Normalized by the mean concentration over
            # the same filtered points so the result is dimensionless (%).
            _diff = px[ratio_mask] - rx[ratio_mask]
            _rmeds = float(np.sqrt(np.median(_diff ** 2)))
            _mean_fj = float(np.mean(np.concatenate([px[ratio_mask], rx[ratio_mask]])))
            sigma = _rmeds / _mean_fj if _mean_fj > 0 else np.nan

            # theta: angle of the mean ratio line from the x-axis (degrees).
            # 45° = perfect 1:1, >45° MCM overestimates, <45° underestimates.
            theta = float(np.degrees(np.arctan(slope)))

            self._cmp_slope_results[sp] = {
                'slope': slope, 'theta': theta, 'sigma': sigma,
                'py': px, 'ref': rx, 'common_t': common_t[mask],
            }

        only_juice = getattr(self, '_cmp_only_juice', [])
        only_ref   = getattr(self, '_cmp_only_ref',   [])

        # Populate species combos — filtered to the current mode:
        #   mode 0 (Favorites): all matched columns
        #   mode 1 (c-Species): only c-columns (not settings)
        #   mode 2 (s-Species): only s-columns (settings)
        #   mode 3 (k-Values): reaction strings from kv files only
        species_list = sorted(self._cmp_slope_results.keys())
        juice_only_items = [f"{c} (nur JUICE)" for c in only_juice]
        ref_only_items   = [f"{c} (nur Ref)"   for c in only_ref]
        if _mode == 3:
            full_list = species_list + juice_only_items + ref_only_items
        elif _mode == 2:
            # s-Species: only settings columns
            full_list = [c for c in species_list if self._is_setting_col(c)]
            full_list += juice_only_items + ref_only_items
        elif _mode == 1:
            # c-Species: only non-settings columns
            full_list = [c for c in species_list if not self._is_setting_col(c)]
            full_list += juice_only_items + ref_only_items
        else:
            full_list = species_list + juice_only_items + ref_only_items
        for slot, combo in enumerate(self._cmp_ts_combos):
            prev = combo.currentText()
            combo.blockSignals(True)
            combo.clear()
            combo.addItems(full_list)
            if prev in full_list:
                combo.setCurrentText(prev)
            combo.blockSignals(False)
            # Redraw whatever is currently selected
            sp = combo.currentText()
            if sp:
                self._cmp_draw_timeseries(slot, sp)

        self._cmp_draw_slope_matrix()

    # ── Sorting helpers ───────────────────────────────────────────────────────
    @staticmethod
    def _atom_count(sp, atom):
        """Count occurrences of a given atom letter in a species name.
        Handles C2, C12, C (=1), etc.  atom must be 'C' or 'N'."""
        import re
        # Find atom followed by a digit (e.g. C2, N12) or atom at end/before
        # uppercase/non-alpha (e.g. lone C in CH3, N in NO2).
        total = 0
        for m in re.finditer(rf'{atom}(\d+)', sp):
            total += int(m.group(1))
        # Count bare atom letters not already captured by the digit pattern
        bare = len(re.findall(rf'{atom}(?!\d)', sp))
        total += bare
        return total

    # Radical patterns: peroxy (ends RO2 / contains O2 after R-group),
    # alkoxy (RO•), and named small radicals.
    _RADICAL_NAMES = {'OH', 'HO2', 'NO3', 'O1D', 'O3', 'HCO', 'CHO',
                      'CH3O', 'C2H5O', 'CH3O2', 'C2H5O2'}
    _RADICAL_SUFFIXES = ('O2', 'OO', 'O•', 'CHO')
    _RADICAL_PREFIXES = ('HO', 'RO')

    @classmethod
    def _is_radical(cls, sp):
        if sp in cls._RADICAL_NAMES:
            return True
        # Peroxy / alkoxy: most MCM names end in O2 (e.g. CH3O2, HOCH2CH2O2)
        if sp.endswith('O2') or sp.endswith('OO'):
            return True
        # Alkoxy radicals end in a single O after a carbon chain
        if len(sp) > 2 and sp.endswith('O') and sp[-2].isdigit():
            return True
        return False

    def _cmp_apply_filter(self):
        """Called whenever filter combo, sort combo, or search box changes."""
        self._cmp_draw_slope_matrix()

    def _cmp_get_filtered_species(self):
        """Return the species list to display, applying mode → filter → sort → search."""
        all_sp = sorted(self._cmp_slope_results.keys())
        if not all_sp:
            return all_sp

        # ── 0. c/s mode filter (same as the dropdown) ────────────────────────
        _mode = getattr(self._cmp_mode_group, 'checkedId', lambda: 1)()
        if _mode == 2:   # s-Species: settings only
            all_sp = [s for s in all_sp if self._is_setting_col(s)]
        elif _mode == 1:  # c-Species: concentrations only
            all_sp = [s for s in all_sp if not self._is_setting_col(s)]
        # mode 0 (all) or 3 (k-Values): no c/s filter

        # ── 1. Filter by agreement quality ───────────────────────────────────
        choice = self._cmp_filter_combo.currentText() if hasattr(self, '_cmp_filter_combo') else "Show all"
        if "Best 50" in choice:
            all_sp = sorted(all_sp, key=lambda s: abs(self._cmp_slope_results[s]['slope'] - 1.0))[:50]
        elif "Best 25" in choice:
            all_sp = sorted(all_sp, key=lambda s: abs(self._cmp_slope_results[s]['slope'] - 1.0))[:25]
        elif "Worst 50" in choice:
            all_sp = sorted(all_sp, key=lambda s: abs(self._cmp_slope_results[s]['slope'] - 1.0), reverse=True)[:50]
        elif "Worst 25" in choice:
            all_sp = sorted(all_sp, key=lambda s: abs(self._cmp_slope_results[s]['slope'] - 1.0), reverse=True)[:25]

        # ── 2. Search box: substring filter (case-insensitive) ────────────────
        query = self._cmp_search_box.text().strip().upper() if hasattr(self, '_cmp_search_box') else ""
        if query:
            all_sp = [s for s in all_sp if query in s.upper()]

        # ── 3. Sort ───────────────────────────────────────────────────────────
        sort_choice = self._cmp_sort_combo.currentText() if hasattr(self, '_cmp_sort_combo') else "Alphabetical"

        if "Alphabetical" in sort_choice:
            all_sp = sorted(all_sp)

        elif "Time of peak" in sort_choice:
            def _t_peak(sp):
                d = self._cmp_slope_results[sp]
                t_arr, py = d['common_t'], d['py']
                if len(py) == 0:
                    return np.inf
                return float(t_arr[np.argmax(np.abs(py))])
            reverse = "latest" in sort_choice
            all_sp = sorted(all_sp, key=_t_peak, reverse=reverse)

        elif "Carbon count" in sort_choice:
            reverse = "high → low" in sort_choice
            all_sp = sorted(all_sp,
                            key=lambda s: (self._atom_count(s, 'C'), s),
                            reverse=reverse)

        elif "Nitrogen count" in sort_choice:
            # Show N≥1 first (sorted low→high), then N=0 alphabetically at the end
            def _n_key(sp):
                n = self._atom_count(sp, 'N')
                return (0 if n >= 1 else 1, n, sp)
            all_sp = sorted(all_sp, key=_n_key)

        elif "Radicals first" in sort_choice:
            # Radicals (+ O3) first alphabetically, then non-radicals alphabetically
            all_sp = sorted(all_sp, key=lambda s: (0 if self._is_radical(s) else 1, s))

        # ── 4. Hard cap: user-chosen number of boxes ──────────────────────────
        if hasattr(self, '_cmp_limit_spin'):
            limit = self._cmp_limit_spin.value()
            if limit < len(all_sp):
                all_sp = all_sp[:limit]

        return all_sp

    def _cmp_draw_slope_matrix(self):
        self._cmp_slope_fig.clear()
        ax = self._cmp_slope_fig.add_subplot(111)
        ax.axis('off')

        species_list = self._cmp_get_filtered_species()
        total_n = len(self._cmp_slope_results)

        # k-Values mode: show all reactions with |slope-1| > 0.2, plus always
        # the 20 worst matches (so the grid is never empty even if all agree).
        _mode = getattr(self._cmp_mode_group, 'checkedId', lambda: 1)()
        if _mode == 3 and self._cmp_slope_results:
            _in_results = [sp for sp in species_list if sp in self._cmp_slope_results]
            _exceeding  = {sp for sp in _in_results
                           if not (0.8 <= self._cmp_slope_results[sp]['slope'] <= 1.2)}
            _worst20    = set(sorted(_in_results,
                                     key=lambda s: abs(self._cmp_slope_results[s]['slope'] - 1.0),
                                     reverse=True)[:20])
            _show = _exceeding | _worst20
            species_list = [sp for sp in species_list if sp in _show]

        n = len(species_list)
        if n == 0:
            msg = ("No species match the current filter / search.\n"
                   "Try 'Show all' or clear the search box."
                   if total_n > 0 else
                   "No significant species found.\nTry lowering the noise threshold.")
            ax.text(0.5, 0.5, msg,
                    ha='center', va='center', transform=ax.transAxes, fontsize=11)
            self._cmp_slope_canvas.draw()
            return

        # Title line: count + active sort + search query + box cap
        sort_label = (self._cmp_sort_combo.currentText().split("(")[0].strip()
                      if hasattr(self, '_cmp_sort_combo') else "")
        search_q   = (self._cmp_search_box.text().strip()
                      if hasattr(self, '_cmp_search_box') else "")
        limit_val  = (self._cmp_limit_spin.value()
                      if hasattr(self, '_cmp_limit_spin') else 9999)
        count_str  = (f"Showing {n} of {total_n} reactions (worst 20 + |slope−1|>0.2)"
                      if _mode == 3 else f"Showing {n} of {total_n} species")
        if sort_label and sort_label != "Alphabetical":
            count_str += f"  ·  {sort_label}"
        if search_q:
            count_str += f"  ·  search: '{search_q}'"
        if limit_val < total_n and n == limit_val:
            count_str += f"  ·  capped at {limit_val}"
        ax.set_title(count_str, fontsize=8, color='#444', pad=2)

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
            sigma = self._cmp_slope_results[sp].get('sigma', np.nan)
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
            # k-Values mode: split "A+B-->C+D" into three lines with --> on line 2
            if _mode == 3 and '-->' in sp:
                parts = sp.split('-->', 1)
                display_name = parts[0] + '\n-->\n' + parts[1]
                name_y  = y + box_h * 0.62
                slope_y = y + box_h * 0.20
                range_y = y + box_h * 0.05
            else:
                display_name = sp
                name_y  = y + box_h * 0.70
                slope_y = y + box_h * 0.42
                range_y = y + box_h * 0.15
            ax.text(x + box_w/2, name_y, display_name,
                    ha='center', va='center', fontsize=fsize_name, fontweight='bold',
                    color=tc, linespacing=1.1)
            ax.text(x + box_w/2, slope_y, f"m={slope:.2f}",
                    ha='center', va='center', fontsize=fsize_val, color=tc)
            if np.isfinite(sigma):
                sigma_str = f"σ={sigma*100:.3g}%"
                ax.text(x + box_w/2, range_y, sigma_str,
                        ha='center', va='center', fontsize=fsize_val * 0.85, color=tc)
            self._cmp_slope_box_map.append((x, y, x+box_w, y+box_h, sp))

        # Tight limits — no extra whitespace around the grid
        ax.set_xlim(-pad/2, ncols*(box_w+pad) - pad/2)
        ax.set_ylim(-pad/2, nrows*(box_h+pad) - pad/2)

        sm = _mcm.ScalarMappable(cmap=cmap, norm=norm); sm.set_array([])
        self._cmp_slope_fig.colorbar(sm, ax=ax, orientation='horizontal',
                                     fraction=0.03, pad=0.01,
                                     label=f'Slope  (JUICE / {self._cmp_ref_label})')
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

        fig    = self._cmp_ts_figs[slot]
        canvas = self._cmp_ts_canvases[slot]
        fig.clear()
        ax = fig.add_subplot(111)

        # ppb toggle state and conversion helper
        use_ppb = (len(self._cmp_ts_ppb_btns) > slot
                   and self._cmp_ts_ppb_btns[slot].isChecked())
        # Quantities that are not mixing ratios — ppb conversion is meaningless for them.
        _ABSOLUTE_QUANTITIES = {
            'M', 'T', 'TEMP', 'P', 'PRESS', 'HV', 'SZ', 'SZA', 'COSX',
            'DIL', 'DILUTE', 'FLOW', 'RH', 'N2', 'O2', 'JSCALE',
        }
        _sp_is_absolute = sp.upper() in _ABSOLUTE_QUANTITIES or sp.upper().startswith('KMT')
        # Rate coefficients (s⁻¹): kOH and named reaction rates k(...)
        _sp_is_rate = sp.upper() == 'KOH' or sp.upper().startswith('K(')

        unit_lbl = ("s⁻¹" if _sp_is_rate else
                    "ppb" if use_ppb and not _sp_is_absolute else "molec/cm³")

        def _to_unit(arr, t_arr_ref):
            if not use_ppb or _sp_is_absolute or _sp_is_rate or 'M' not in self.df.columns:
                return arr
            M_interp = np.interp(t_arr_ref,
                                 self.df.index.values.astype(float),
                                 self.df['M'].values.astype(float))
            M_interp = np.where(M_interp > 0, M_interp, 1.0)
            return arr / M_interp * 1e9

        # Variable not in slope results — plot what we have
        if data is None:
            if sp not in self.df.columns:
                canvas.draw()
                return
            t_arr = self.df.index.values.astype(float)
            vals  = _to_unit(self.df[sp].values.astype(float), t_arr)
            ax.plot(t_arr, vals, lw=1.5, color='steelblue', label=f"{sp} JUICE")
            # Also plot reference if available under the same name
            ref_vals_raw = None
            if self._cmp_df_ref is not None and sp in self._cmp_df_ref.columns:
                t_ref = self._cmp_df_ref.index.values.astype(float)
                rv = _to_unit(self._cmp_df_ref[sp].values.astype(float), t_ref)
                ax.plot(t_ref, rv, lw=1.5, color='darkorange', ls='--',
                        label=f"{sp} {self._cmp_ref_label}")
                ref_vals_raw = np.interp(t_arr, t_ref, rv)
            ax.set_ylabel(f"[{unit_lbl}]", fontsize=9)
            ax.set_title(sp, fontweight='bold', fontsize=10, loc='left')
            ax.xaxis.set_major_formatter(FuncFormatter(_fmt_hhmm))
            ax.set_xlabel("Time of Day (UTC)")
            ax.yaxis.get_major_formatter().set_useOffset(False)
            ax.legend(loc='lower right', bbox_to_anchor=(1, 1.01), ncol=1,
                      fontsize=7, borderaxespad=0, framealpha=0.85).set_draggable(True)
            ax.grid(True, ls='--', alpha=0.35)
            fig.tight_layout(pad=0.8)
            fig.subplots_adjust(left=0.18)
            canvas.draw()
            self._cmp_ts_data[slot] = (vals, ref_vals_raw, t_arr, sp)
            return

        slope, theta = data['slope'], data['theta']
        sigma = data.get('sigma', np.nan)
        py_vals  = _to_unit(data['py'],  data['common_t'])
        ref_vals = _to_unit(data['ref'], data['common_t'])
        common_t = data['common_t']

        # Clear any stale crosshair artists before redrawing
        self._cmp_crosshair[slot] = {'vline': None, 'text': None}

        hv_col = next((c for c in self.df.columns if c.upper() in ('SAPHIR_HV', 'HV')), None)
        if hv_col is not None:
            hv_interp = np.interp(common_t, self.df.index.values.astype(float),
                                  self.df[hv_col].values.astype(float),
                                  left=0.0, right=0.0)
            open_mask = hv_interp > 0
            transitions = np.diff(open_mask.astype(int))
            starts = list(np.where(transitions == 1)[0] + 1)
            ends   = list(np.where(transitions == -1)[0] + 1)
            if open_mask[0]:  starts = [0] + starts
            if open_mask[-1]: ends   = ends + [len(common_t) - 1]
            first = True
            for si, ei in zip(starts, ends):
                ax.axvspan(common_t[si], common_t[ei], alpha=0.12, color='gold',
                           label='Roof open' if first else '', zorder=0)
                first = False

        def _scale_for_label(arr):
            """Return (scaled_array, scale_str) so tick values are in 0.1–1000 range."""
            max_abs = np.nanmax(np.abs(arr)) if len(arr) else 0.0
            if max_abs == 0:
                return arr, ""
            exp = int(np.floor(np.log10(max_abs)))
            if -2 <= exp <= 3:          # already readable — no scaling
                return arr, ""
            # Build superscript string using LaTeX mathtext
            scale = 10.0 ** exp
            sup = f"$\\times 10^{{{exp}}}$"
            return arr / scale, sup

        show_diff = (len(self._cmp_ts_diff_btns) > slot
                     and self._cmp_ts_diff_btns[slot].isChecked())

        # Scale left-axis data: JUICE and FAC must use the SAME scale factor so
        # both curves are comparable on the same axis.  Applying _scale_for_label
        # independently gives each its own exponent — a 10× difference in values
        # then looks like a match because FAC at 4 (×10⁻⁶) plots at the same
        # height as JUICE at 3.5 (×10⁻⁵) while the axis only shows ×10⁻⁵.
        _combined = np.concatenate([
            py_vals,
            ref_vals if ref_vals is not None else np.array([])
        ])
        _combined_finite = _combined[np.isfinite(_combined)]

        # Decide whether to subtract a mean offset before scaling.
        # When data clusters around a large non-zero mean (like NO2 ≈ 10 ppb
        # with ±1e-7 ppb variation), the variation/mean ratio is tiny.
        # matplotlib's useOffset=True then injects "+9.999999" into the tick
        # exponent area, creating confusing "1e-7+9.999999 ppb" labels.
        # Fix: subtract the mean explicitly; show offset in the axis label;
        # the unit bracket then contains ONLY the pure scale factor (×10^n).
        _lax_offset = 0.0
        if len(_combined_finite) > 0:
            _cmean    = float(np.nanmean(_combined_finite))
            _cvariation = float(np.nanmax(np.abs(_combined_finite - _cmean)))
            _cmax_abs   = float(np.nanmax(np.abs(_combined_finite)))
            if _cmax_abs > 0 and _cvariation / _cmax_abs < 0.01:
                # Clustered: subtract mean, scale the variation
                _lax_offset = _cmean
                _work = _combined_finite - _cmean
            else:
                _work = _combined_finite
        else:
            _work = _combined_finite

        _max_abs = float(np.nanmax(np.abs(_work))) if len(_work) > 0 else 0.0
        if _max_abs > 0:
            _exp = int(np.floor(np.log10(_max_abs)))
            if -2 <= _exp <= 3:
                _scale, lax_sup = 1.0, ""
            else:
                _scale = 10.0 ** _exp
                lax_sup = f"$\\times 10^{{{_exp}}}$"
        else:
            _scale, lax_sup = 1.0, ""

        py_plot  = (py_vals  - _lax_offset) / _scale
        ref_plot = (ref_vals - _lax_offset) / _scale if ref_vals is not None else None

        # Axis label: scale in unit bracket; offset shown as separate annotation
        if lax_sup:
            lax_lbl = f"{sp}  [{lax_sup} {unit_lbl}]"
        else:
            lax_lbl = f"{sp}  [{unit_lbl}]"
        if abs(_lax_offset) > 0:
            # Format offset compactly — e.g. "+10.0 ppb" or "+2500 ppb"
            _off_str = f"{_lax_offset:.5g} {unit_lbl}"
            lax_lbl += f"  (+{_off_str})"

        ax.plot(common_t, py_plot,  label=f"{sp} JUICE",
                lw=1.5, color='steelblue', zorder=2)
        if ref_plot is not None:
            ax.plot(common_t, ref_plot, label=f"{sp} {self._cmp_ref_label}",
                    lw=1.5, color='darkorange', ls='--', zorder=2)
        ax.xaxis.set_major_formatter(FuncFormatter(_fmt_hhmm))
        ax.set_xlabel("Time of Day (UTC)", fontsize=9)
        ax.set_ylabel(lax_lbl, fontsize=9)
        _sig_str = f"   σ={sigma*100:.3g}%" if np.isfinite(sigma) else ""
        ax.set_title(f"{sp}   m={slope:.3f}{_sig_str}",
                     fontweight='bold', fontsize=10, loc='left')
        # Prevent matplotlib from adding its own "+offset" to the tick exponent area.
        # Our lax_lbl already carries the only scale/offset information shown.
        # Use a fresh ScalarFormatter with useOffset=False — more reliable than
        # modifying the existing formatter (which may be a different type).
        import matplotlib.ticker as _mticker
        ax.yaxis.set_major_formatter(_mticker.ScalarFormatter(useOffset=False))

        if show_diff and ref_vals is not None:
            diff = py_vals - ref_vals
            diff_plot, diff_sup = _scale_for_label(diff)
            diff_lbl = f"Δ  [{diff_sup} {unit_lbl}]" if diff_sup else f"Δ  [{unit_lbl}]"
            ax2 = ax.twinx()
            ax2.plot(common_t, diff_plot, lw=1.2, color='green', ls=':', zorder=3,
                     label="Δ (JUICE − Ref)")
            ax2.axhline(0, color='green', lw=0.5, ls='--', alpha=0.5)
            ax2.set_ylabel(diff_lbl, fontsize=9, color='green')
            ax2.tick_params(axis='y', labelcolor='green')
            ax2.xaxis.set_major_formatter(FuncFormatter(_fmt_hhmm))
            h1, l1 = ax.get_legend_handles_labels()
            h2, l2 = ax2.get_legend_handles_labels()
            ax.legend(h1 + h2, l1 + l2, loc='lower right',
                      bbox_to_anchor=(1, 1.01), ncol=1,
                      fontsize=7, borderaxespad=0, framealpha=0.85).set_draggable(True)
        else:
            ax.legend(loc='lower right', bbox_to_anchor=(1, 1.01), ncol=1,
                      fontsize=7, borderaxespad=0, framealpha=0.85).set_draggable(True)

        self._cmp_ts_data[slot] = (py_vals, ref_vals, common_t, sp)
        ax.grid(True, ls='--', alpha=0.35)
        fig.tight_layout(pad=0.8)
        fig.subplots_adjust(left=0.18)   # prevent y-label clipping
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

        def _crosshair_text(ax_ref, t_arr, py_v, ref_v, x_val):
            py_at = float(np.interp(x_val, t_arr, py_v))
            h, rem = divmod(max(0, int(x_val * 3600)), 3600)
            m = rem // 60
            if ref_v is not None:
                ref_at = float(np.interp(x_val, t_arr, ref_v))
                pct_txt = f"   Agr={py_at/ref_at*100:.1f}%" if ref_at != 0 else ""
                txt = f"{h:02d}:{m:02d}   MCM={py_at:.3g}   Ref={ref_at:.3g}{pct_txt}"
            else:
                txt = f"{h:02d}:{m:02d}   {py_at:.3g}"
            return ax_ref.text(0.0, 1.01, txt, transform=ax_ref.transAxes,
                               va='bottom', ha='left', fontsize=8,
                               bbox=dict(boxstyle='round,pad=0.2',
                                         facecolor='white', alpha=0.9),
                               zorder=11, clip_on=False)

        data = self._cmp_ts_data[slot]
        if data is not None:
            py_vals, ref_vals, common_t, sp = data
            ch['text'] = _crosshair_text(ax, common_t, py_vals, ref_vals, x)

        # Draw mirror vline + text in the peer canvas
        other_axes = self._cmp_ts_figs[other].axes
        if other_axes:
            other_ax = other_axes[0]
            self._cmp_mirror_vlines[other] = other_ax.axvline(
                x, color='gray', lw=0.8, ls='--', zorder=10)
            other_data = self._cmp_ts_data[other]
            if other_data is not None:
                py_vals2, ref_vals2, common_t2, sp2 = other_data
                self._cmp_mirror_texts[other] = _crosshair_text(
                    other_ax, common_t2, py_vals2, ref_vals2, x)
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

        ctrl_layout.addWidget(QLabel("N view:"))
        self._nb_view = QComboBox(); self._nb_view.addItems(["All nitrogen", "PANs", "Org. Nitrates"])
        self._nb_view.setToolTip("nitrogen figure: every class, or the species of one class (PANs / organic nitrates) individually")
        self._nb_view.currentIndexChanged.connect(self._update_nitrogen_budget)
        ctrl_layout.addWidget(self._nb_view)
        ctrl_layout.addWidget(QLabel("Exclude:"))
        self._cb_exclude = QLineEdit("CO2"); self._cb_exclude.setFixedWidth(220)
        self._cb_exclude.setPlaceholderText("names or formulas, e.g. CO, CH4, C2H6O")
        self._cb_exclude.setToolTip("Species left out of both budgets: MCM names or chemical formulas, comma-separated. "
                                    "A formula removes every species with that formula (isomers included).")
        self._cb_exclude.editingFinished.connect(self._update_both_budgets)
        ctrl_layout.addWidget(self._cb_exclude)

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
        self._cb_canvas = FigureCanvas(self._cb_fig); _interact(self._cb_canvas)
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
        self._nb_canvas = FigureCanvas(self._nb_fig); _interact(self._nb_canvas)
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

        parent_voc = [sp for sp in self._parent_species(_carbon_count) if sp not in ('CO', 'CO2', 'CH4')]   # tracers are not VOCs: CO keeps its own class

        # ── Collect carbon-containing MCM species ─────────────────────────
        species_data = {}    # name → (n_carbons, ppb_C array)
        unknown_species = [] # non-zero species with truly unresolvable C count
        co2_present    = False

        for col in self.df.columns:
            if col == 'M' or col.upper().startswith('SAPHIR_'):
                continue
            if not self._is_model_species(col) or self._budget_excluded(col):
                continue
            n_c = _carbon_count(col)
            vals = self.df[col].values
            ppb = (vals / M) * 1e9 if M is not None else vals.copy()
            if n_c == 0:
                # Distinguish known zero-carbon species (inorganic, in DB or explicit
                # inorganic list) from truly unidentified species not resolvable at all.
                upper = col.upper()
                known_zero = upper in MCM_DB or upper in KNOWN_ZERO_CARBON
                if not known_zero and np.nanmax(np.abs(ppb)) > 1e-3:
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
            self._cb_status_lbl.setText(f"⚠ C count unknown: {shown}{tail}")
            self._cb_status_lbl.setToolTip(
                "Non-zero species not in MCM DB — carbon count unresolvable:\n" +
                '\n'.join(unknown_species))
        else:
            self._cb_status_lbl.setText("✓ All species accounted for")
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
            names, colors = self._cb_draw_all_species(ax, t, species_data)
        else:
            names, colors = self._cb_draw_grouped(ax, t, species_data, parent_voc)

        # ── Total C overlay ───────────────────────────────────────────────
        if self._cb_total_chk.isChecked():
            total = sum(np.nan_to_num(ppb_c) for _, ppb_c in species_data.values())
            ax.plot(t, total, color='black', lw=1.5, ls='--',
                    label='Total tracked C', zorder=10)

        # Species with no resolvable carbon count cannot be budgeted: they go into the Exclude field
        # (visible and editable there) instead of a note in the figure.
        if unknown_species and hasattr(self, '_cb_exclude'):
            have = [x.strip() for x in self._cb_exclude.text().split(',') if x.strip()]
            new = [sp for sp in unknown_species if sp.upper() not in {h.upper() for h in have}]
            if new:
                self._cb_exclude.setText(", ".join(have + new))

        # ── Axes formatting ───────────────────────────────────────────────
        ax.set_xlabel('Time of Day (UTC)')
        ax.set_ylabel('Carbon content (ppb-C)')
        ax.xaxis.set_major_formatter(FuncFormatter(_fmt_hhmm))
        ax.set_xlim(t[0], t[-1])
        ax.set_ylim(bottom=0)
        ax.grid(True, linestyle='--', alpha=0.4)

        handles, labels = ax.get_legend_handles_labels()
        if handles:
            ncol = max(1, len(handles) // 20 + 1) if len(handles) > 15 else 1
            ax.legend(handles, labels, loc='upper right',
                      fontsize=6 if len(handles) > 15 else 8,
                      ncol=ncol, framealpha=0.85).set_draggable(True)

        self._cb_fig.tight_layout()
        self._cb_canvas.draw()

    def _budget_excluded(self, col: str) -> bool:
        """True when the Budgets 'Exclude:' field names this species or its formula."""
        txt = self._cb_exclude.text() if hasattr(self, '_cb_exclude') else ''
        if not txt.strip():
            return False
        excl = {x.strip().upper() for x in txt.replace(';', ',').split(',') if x.strip()}
        u = col.upper()
        if u in excl:
            return True
        e = MCM_DB.get(u)
        return bool(e is not None and (getattr(e, 'formula', '') or '').upper() in excl)

    def _is_model_species(self, col: str) -> bool:
        """A result column that is a chemical species: in the species table, or an ODE species of
        the mechanism the run used.  Everything else in a .cs (ENZ switch columns, J values, rate
        constants, diagnostic sums such as XNOY/NOx/Ox) must not enter an element budget."""
        u = col.upper()
        if u in ('XNOY', 'NOX', 'NOY', 'OX'):          # bookkeeping tracers of the EASY script (injected NOy, sums), not molecules
            return False
        if MCM_DB.get(u) is not None:
            return True
        c_names = getattr(self, '_easy_c_names', None)
        if c_names:
            return u.lower() in c_names
        return not (u.startswith('Q') or u.startswith('J') or u.startswith('K') or u in ('XNOY', 'NOX', 'NOY', 'OX', 'Y', 'QY'))

    def _parent_species(self, count_fn):
        """Injected precursors that carry the budget's element (count_fn > 0): from the scenario's
        burst injections and from the experiment description next to the result (injected products
        and chosen compounds).  Ordered as injected."""
        out = []
        if self.budget_engine and self.budget_engine.cfg:
            for inj in self.budget_engine.cfg.get('burst_injections', []):
                sp = (inj.get('species', '') or '').upper()
                if sp and count_fn(sp) > 0 and sp not in out:
                    out.append(sp)
        try:
            import experiment_plots as _ep
            _spec = _ep.spec_for_run(self.csv_path) if getattr(self, "csv_path", "") else None
            if _spec is not None:
                for sp in list(_ep.injected_species(_spec)) + [p.upper() for p in _spec.precursors]:
                    if count_fn(sp) > 0 and sp not in out:
                        out.append(sp)
        except Exception:
            pass
        return out

    @staticmethod
    def _parent_bands(species_data, parents, base_color, prefix):
        """One band per injected precursor present in the data, in the family of `base_color`:
        shades + hatches so the precursors can be told apart.  Precursors whose peak is below
        _CB_OTHER_SPLIT of the summed parent peak share one band '<prefix> (minor)' — a mixture of
        30 compounds would otherwise drown the plot.  Returns (names, arrays, colors, hatches)."""
        present = [(sp, np.nan_to_num(species_data[sp][1])) for sp in parents
                   if sp in species_data and np.nanmax(np.nan_to_num(species_data[sp][1])) > 1e-9]
        names, arrays, colors, hatches = [], [], [], []
        if not present:
            return names, arrays, colors, hatches
        peak = float(np.nanmax(sum(a for _, a in present)))
        thr = 0.01 * peak                                          # own band down to 1 % of the parent peak
        big = [(sp, a) for sp, a in present if np.nanmax(a) >= thr]
        minor = [a for sp, a in present if np.nanmax(a) < thr]
        big.sort(key=lambda kv: -float(np.nanmean(kv[1])))         # largest first = lightest … darkest
        n = len(big) + (1 if minor else 0)
        for i, (sp, a) in enumerate(big):
            names.append(f"{prefix}: {sp}"); arrays.append(a)
            colors.append(_shade(base_color, i, n)); hatches.append('')
        if minor:
            names.append(f"{prefix} (minor, {len(minor)} compounds)"); arrays.append(sum(minor))
            colors.append(_shade(base_color, n - 1, n)); hatches.append('')
        return names, arrays, colors, hatches

    def _cb_draw_all_species(self, ax, t, species_data):
        """Stacked area with one band per species, sorted by mean contribution.
        Returns (names, colors) for legend building."""
        sorted_items = sorted(species_data.items(),
                              key=lambda x: np.nanmean(x[1][1]), reverse=True)
        names  = [item[0] for item in sorted_items]
        arrays = [np.nan_to_num(item[1][1]) for item in sorted_items]
        cmap   = _mcm.get_cmap('tab20')
        colors = [cmap((i % 20) / 20.0) for i in range(len(names))]
        ax.stackplot(t, *arrays, labels=names, colors=colors, alpha=0.85)
        return names, colors

    def _cb_draw_grouped(self, ax, t, species_data, parent_voc):
        """Stacked area with one band per chemical class using fixed class colours.
        Returns (names, colors) for legend building."""
        parent_voc = list(parent_voc)
        p_names, p_arrays, p_colors, p_hatches = self._parent_bands(species_data, parent_voc, _CB_CLASS_COLORS['Parent VOC'], "Parent VOC")
        class_totals = {}
        for sp, (_, ppb_c) in species_data.items():
            cls = _carbon_class(sp, parent_voc)
            if cls == 'Parent VOC':
                continue                                   # drawn as its own band below
            if cls not in class_totals:
                class_totals[cls] = np.zeros(len(t))
            class_totals[cls] += np.nan_to_num(ppb_c)

        # "Other" must not hide a major contributor: any species classed Other whose peak
        # exceeds _CB_OTHER_SPLIT of the total carbon peak gets its own band (named after it)
        total = sum(np.nan_to_num(v[1]) for v in species_data.values()) if species_data else np.zeros(len(t))
        split, rest = [], np.zeros(len(t))
        if 'Other' in class_totals and np.nanmax(total) > 0:
            thr = _CB_OTHER_SPLIT * float(np.nanmax(total))
            for sp, (_, ppb_c) in species_data.items():
                if _carbon_class(sp, parent_voc) == 'Other':
                    arr = np.nan_to_num(ppb_c)
                    (split.append((sp, arr)) if np.nanmax(arr) >= thr else rest.__iadd__(arr))
            split.sort(key=lambda kv: -float(np.nanmean(kv[1])))
            class_totals['Other'] = rest
        names, arrays, colors, hatches = [], [], [], []
        _split_cmap = _mcm.get_cmap('Dark2')
        for cls in _CB_CLASS_ORDER:
            if cls == 'Parent VOC':
                names += p_names; arrays += p_arrays; colors += p_colors; hatches += p_hatches
                continue
            if cls == 'Other':
                for i, (sp, arr) in enumerate(split):          # just below the remaining "Other"
                    names.append(f"{sp} (from Other)"); arrays.append(arr); colors.append(_split_cmap(i % 8)); hatches.append('')
            if cls in class_totals and np.nanmax(class_totals[cls]) > 1e-9:
                names.append(cls)
                arrays.append(class_totals[cls])
                colors.append(_CB_CLASS_COLORS.get(cls, '#aaaaaa')); hatches.append('')
        for cls, arr in class_totals.items():
            if cls not in _CB_CLASS_ORDER and np.nanmax(arr) > 1e-9:
                names.append(cls)
                arrays.append(arr)
                colors.append('#aaaaaa'); hatches.append('')

        if arrays:
            ax.stackplot(t, *arrays, labels=names, colors=colors, hatch=hatches, alpha=0.85)
        return names, colors

    # ── Nitrogen budget ───────────────────────────────────────────────────────
    def _update_nitrogen_budget(self):
        """Compute nitrogen-weighted concentrations and draw the stacked timeseries."""
        self._nb_fig.clear()
        ax = self._nb_fig.add_subplot(111)

        t = self.df.index.values
        M = self.df['M'].values if 'M' in self.df.columns else None

        species_data = {}   # name → (n_N, ppb_N array)

        for col in self.df.columns:
            if col == 'M' or col.upper().startswith('SAPHIR_'):
                continue
            if not self._is_model_species(col) or self._budget_excluded(col):
                continue  # ENZ switches (QNO), J values (JNO2), rate constants (KRO2NO), sums (XNOY) are not species
            n_n = _nitrogen_count(col)
            if n_n == 0:
                continue  # known zero-N (CO, C2H5OH, etc.) — silently skip
            vals = self.df[col].values
            ppb  = (vals / M) * 1e9 if M is not None else vals.copy()
            ppb_n = np.clip(np.nan_to_num(ppb, nan=0.0, posinf=0.0, neginf=0.0), 0.0, None) * n_n
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

        names, colors = [], []
        mode = self._cb_mode_combo.currentText()
        view = self._nb_view.currentText() if hasattr(self, '_nb_view') else "All nitrogen"
        if view in ("PANs", "Org. Nitrates"):
            # one class, species by species: the big ones named, the rest summed into Other
            members = {sp: d for sp, d in species_data.items() if _nitrogen_class(sp) == view}
            if not members:
                ax.text(0.5, 0.5, f"No {view} in this run.", ha='center', va='center', transform=ax.transAxes, fontsize=11)
                self._nb_canvas.draw(); return
            species_data = members
            total = sum(np.nan_to_num(d[1]) for d in members.values()); thr = _CB_OTHER_SPLIT * float(np.nanmax(total)) if np.nanmax(total) > 0 else 0.0
            big = sorted(((sp, np.nan_to_num(d[1])) for sp, d in members.items() if np.nanmax(np.nan_to_num(d[1])) >= thr), key=lambda kv: -float(np.nanmean(kv[1])))
            rest = [np.nan_to_num(d[1]) for sp, d in members.items() if np.nanmax(np.nan_to_num(d[1])) < thr]
            names = [sp for sp, _a in big] + ([f"Other {view} ({len(rest)})"] if rest else [])
            arrays = [a for _sp, a in big] + ([sum(rest)] if rest else [])
            colors = _distinct_colors(names)
            ax.stackplot(t, *arrays, labels=names, colors=colors, alpha=0.85)
            ax.set_title(f"{view}: species contributions (ppb-N)", fontsize=10)
        elif mode == "All Species":
            sorted_items = sorted(species_data.items(),
                                  key=lambda x: np.nanmean(x[1][1]), reverse=True)
            names  = [item[0] for item in sorted_items]
            arrays = [np.nan_to_num(item[1][1]) for item in sorted_items]
            cmap   = _mcm.get_cmap('tab20')
            colors = [cmap((i % 20) / 20.0) for i in range(len(names))]
            ax.stackplot(t, *arrays, labels=names, colors=colors, alpha=0.85)
        else:
            parents = self._parent_species(_nitrogen_count)          # injected NO, NO2, … as own bands
            class_totals = {}
            for sp, (_, ppb_n) in species_data.items():
                if sp in parents:
                    continue
                cls = _nitrogen_class(sp)
                if cls not in class_totals:
                    class_totals[cls] = np.zeros(len(t))
                class_totals[cls] += np.nan_to_num(ppb_n)
            names, arrays, colors, hatches = [], [], [], []
            for cls in _NB_CLASS_ORDER:
                p_in_cls = [sp for sp in parents if _nitrogen_class(sp) == cls]
                if p_in_cls:                                            # injected members first, in the class colour family
                    pn, pa, pc, ph = self._parent_bands(species_data, p_in_cls, _NB_CLASS_COLORS.get(cls, '#aaaaaa'), "Injected")
                    names += pn; arrays += pa; colors += pc; hatches += ph
                if cls in class_totals and np.nanmax(class_totals[cls]) > 1e-9:
                    names.append(cls + (" (other)" if p_in_cls else ""))
                    arrays.append(class_totals[cls])
                    colors.append(_NB_CLASS_COLORS.get(cls, '#aaaaaa')); hatches.append('')
            for cls, arr in class_totals.items():
                if cls not in _NB_CLASS_ORDER and np.nanmax(arr) > 1e-9:
                    names.append(cls)
                    arrays.append(arr)
                    colors.append('#aaaaaa'); hatches.append('')
            if arrays:
                ax.stackplot(t, *arrays, labels=names, colors=colors, hatch=hatches, alpha=0.85)

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
            ncol = max(1, len(handles) // 20 + 1) if len(handles) > 15 else 1
            ax.legend(handles, labels, loc='upper right',
                      fontsize=6 if len(handles) > 15 else 8,
                      ncol=ncol, framealpha=0.85).set_draggable(True)

        self._nb_fig.tight_layout()
        self._nb_canvas.draw()

    # ═════════════════════════════════════════════════════════════════════════
    # Chemical Space tab (VBS + molecular distribution)
    # ═════════════════════════════════════════════════════════════════════════
    # Formulas:
    #   • Volatility  LogC*  — Li et al., ACP 16, 3327 (2016), doi:10.5194/acp-16-3327-2016
    #                          "Molecular corridors and parameterizations of volatility..."
    #   • OSc         — Kroll et al., Nat. Chem. 3, 133 (2011), OSc ≈ 2·nO/nC − nH/nC
    #   • Only organic species (nC ≥ 1); radicals (RO2, RO, OH, HO2) excluded by default
    #     because (a) they are not realistic condensable vapours and (b) their
    #     lifetimes are too short to define a meaningful mixing ratio.  Average radical
    #     concentrations (OH, HO2, O3, RO2) are shown in a side panel instead.
    # ─────────────────────────────────────────────────────────────────────────

    def _cs_calc_logc_star(self, nC, nO, nN, nS=0):
        """Li et al. 2016 LogC* parameterization (µg/m³ at 300 K).
        Choses coefficient set based on which heteroatoms are present."""
        if nC <= 0:
            return np.nan
        # CH only
        if nO == 0 and nN == 0 and nS == 0:
            n0C, bC, bO, bCO, bN, bS = 23.80, 0.4861, 0.0, 0.0, 0.0, 0.0
        elif nO > 0 and nN == 0 and nS == 0:      # CHO
            n0C, bC, bO, bCO, bN, bS = 22.66, 0.4481, 1.656, -0.7790, 0.0, 0.0
        elif nN > 0 and nO == 0 and nS == 0:      # CHN
            n0C, bC, bO, bCO, bN, bS = 24.59, 0.4066, 0.0, 0.0, 0.9619, 0.0
        elif nO > 0 and nN > 0 and nS == 0:       # CHON
            n0C, bC, bO, bCO, bN, bS = 24.13, 0.3667, 0.7732, -0.07790, 1.114, 0.0
        elif nO > 0 and nS > 0 and nN == 0:       # CHOS
            n0C, bC, bO, bCO, bN, bS = 24.06, 0.3637, 1.327, -0.3988, 0.0, 0.7579
        elif nO > 0 and nN > 0 and nS > 0:        # CHONS
            n0C, bC, bO, bCO, bN, bS = 28.50, 0.3848, 1.011, 0.2921, 1.053, 1.316
        else:
            return np.nan
        return (n0C - nC)*bC - nO*bO - 2.0*(nC*nO/max(nC+nO, 1))*bCO - nN*bN - nS*bS

    def _cs_calc_osc(self, nC, nH, nO):
        if nC <= 0:
            return np.nan
        return (2.0*nO - nH) / nC

    def _cs_classify(self, nC, nH, nO, nN):
        """CH / CHO / CHON1 / CHON2 / CHON3+ family classification."""
        if nC == 0:
            return "Other"
        if nO == 0 and nN == 0:
            return "CH"
        if nO > 0 and nN == 0:
            return "CHO"
        if nN == 1:
            return "CHON1"
        if nN == 2:
            return "CHON2"
        if nN >= 3:
            return "CHON3+"
        return "Other"

    _CS_CLASS_ORDER = ["CH", "CHO", "CHON1", "CHON2", "CHON3+", "Other"]
    _CS_CLASS_COLORS = {
        "CH":     "#1f77b4",
        "CHO":    "#2ca02c",
        "CHON1":  "#ff7f0e",
        "CHON2":  "#d62728",
        "CHON3+": "#9467bd",
        "Other":  "#7f7f7f",
    }

    def _cs_build_species_table(self):
        """Build a table of organic species with their atom counts and derived properties.
        Returns a pandas DataFrame with columns: name, nC, nH, nO, nN, mass, class, logc, osc, is_radical.
        Cached on self._cs_species_table after first call."""
        if hasattr(self, '_cs_species_table'):
            return self._cs_species_table
        rows = []
        for col in self.df.columns:
            if col.startswith('SAPHIR_') or col in ('M', 'TIME_s'):
                continue
            entry = MCM_DB.get(col.upper())
            if entry is None or entry.n_C < 1:
                continue
            is_radical = bool(entry.is_peroxy) or col.upper().endswith('O') and col.upper() not in {'HCHO', 'CH3CHO'}
            # Simpler radical rule: ends with O2 / has is_peroxy flag / is OH/HO2
            is_radical = bool(entry.is_peroxy) or col.upper() in {'OH', 'HO2'} or col.upper().endswith('O2')
            rows.append({
                'name':    col,
                'formula': entry.formula or '',
                'nC':      entry.n_C, 'nH': entry.n_H, 'nO': entry.n_O, 'nN': entry.n_N,
                'mass':    entry.mass,
                'class':   self._cs_classify(entry.n_C, entry.n_H, entry.n_O, entry.n_N),
                'logc':    self._cs_calc_logc_star(entry.n_C, entry.n_O, entry.n_N),
                'osc':     self._cs_calc_osc(entry.n_C, entry.n_H, entry.n_O),
                'is_radical': is_radical,
            })
        self._cs_species_table = pd.DataFrame(rows)
        return self._cs_species_table

    def _cs_current_window(self):
        """Return (t_lo, t_hi) in hours since midnight for the current time selector.
        Snapshot mode → zero-width window at _cs_time_point.
        Window mode    → the draggable _cs_time_range."""
        mode_window = (hasattr(self, '_cs_mode_combo')
                       and self._cs_mode_combo.currentIndex() == 1)
        if mode_window:
            if not hasattr(self, '_cs_time_range'):
                t = self.df.index.values.astype(float)
                return (float(t[0]), float(t[-1]))
            return self._cs_time_range
        # Snapshot
        if not hasattr(self, '_cs_time_point'):
            t = self.df.index.values.astype(float)
            return (float(t[len(t)//2]), float(t[len(t)//2]))
        return (self._cs_time_point, self._cs_time_point)

    def _cs_avg_conc_ppb(self, species):
        """Concentration (ppb) of a species, averaged over the current window
        or nearest-timestep lookup in snapshot mode."""
        if species not in self.df.columns or 'M' not in self.df.columns:
            return np.nan
        t = self.df.index.values.astype(float)
        t_lo, t_hi = self._cs_current_window()
        if t_lo == t_hi:
            # Snapshot: nearest row
            idx = int(np.argmin(np.abs(t - t_lo)))
            v = float(self.df[species].values[idx])
            m = float(self.df['M'].values[idx])
            return v / m * 1e9 if m > 0 else np.nan
        mask = (t >= t_lo) & (t <= t_hi)
        if mask.sum() == 0:
            return np.nan
        vals = self.df[species].values[mask]   # molec/cm³
        M    = self.df['M'].values[mask]
        with np.errstate(divide='ignore', invalid='ignore'):
            ppb = np.where(M > 0, vals / M * 1e9, np.nan)
        return float(np.nanmean(ppb))

    def _cs_avg_conc_molec(self, species):
        """Concentration (molec/cm³) of a species, averaged over the current window
        or nearest-timestep lookup in snapshot mode."""
        if species not in self.df.columns:
            return np.nan
        t = self.df.index.values.astype(float)
        t_lo, t_hi = self._cs_current_window()
        if t_lo == t_hi:
            idx = int(np.argmin(np.abs(t - t_lo)))
            return float(self.df[species].values[idx])
        mask = (t >= t_lo) & (t <= t_hi)
        if mask.sum() == 0:
            return np.nan
        return float(np.nanmean(self.df[species].values[mask]))

    def _build_chemspace_tab(self):
        """Chemical Space tab: injection timeseries (top-left), VBS (top-right),
        C/O distribution (bottom-left), CHO family groups (bottom-right)."""
        outer = QWidget()
        outer_layout = QVBoxLayout(outer)
        outer_layout.setContentsMargins(2, 2, 2, 2)
        outer_layout.setSpacing(2)

        # ── Controls row ──────────────────────────────────────────────────
        ctrl = QHBoxLayout()
        ctrl.setSpacing(5)
        ctrl.addWidget(QLabel("Mode:"))
        self._cs_mode_combo = QComboBox()
        self._cs_mode_combo.addItems(["Snapshot (single time)", "Window (time range)"])
        self._cs_mode_combo.setCurrentIndex(0)  # Snapshot is the default
        self._cs_mode_combo.setToolTip(
            "Snapshot: all plots reflect chamber state at one instant.\n"
            "Window: all plots show averages over a draggable time range."
        )
        ctrl.addWidget(self._cs_mode_combo)

        ctrl.addWidget(QLabel("Units:"))
        self._cs_units_combo = QComboBox()
        self._cs_units_combo.addItems(["ppb", "µg/m³"])
        ctrl.addWidget(self._cs_units_combo)

        ctrl.addWidget(QLabel("Exclude:"))
        self._cs_exclude = QLineEdit()
        self._cs_exclude.setPlaceholderText("names or formulas, e.g. CH4, CO, C2H6O")
        self._cs_exclude.setToolTip("Species to leave out of every chemical-space panel: MCM names or chemical formulas, comma-separated. "
                                    "A formula removes every species with that formula (isomers included).")
        self._cs_exclude.setFixedWidth(220)
        self._cs_exclude.editingFinished.connect(lambda: self._cs_redraw_all() if hasattr(self, '_cs_redraw_all') else None)
        ctrl.addWidget(self._cs_exclude)

        ctrl.addWidget(QLabel("VBS size scale:"))
        self._cs_size_combo = QComboBox()
        self._cs_size_combo.addItems(["linear", "log10"])
        self._cs_size_combo.setCurrentIndex(1)
        ctrl.addWidget(self._cs_size_combo)

        # Video export controls (no live playback — it crashed the GUI).
        # The cursor remains fully interactive manually; clicking or dragging
        # on the timeseries in Snapshot mode updates all panels instantly.
        ctrl.addStretch(1)                               # keep every control right next to its label
        outer_layout.addLayout(ctrl, 0)                 # row 1: what is shown
        ctrl = QHBoxLayout(); ctrl.setSpacing(5)         # row 2: video export + readout (one row pushed the tab off-screen)
        ctrl.addWidget(QLabel("Video from"))
        self._cs_vid_start = QLineEdit(); self._cs_vid_start.setFixedWidth(52); self._cs_vid_start.setPlaceholderText("HH:MM")
        self._cs_vid_start.setToolTip("first frame of the video (HH:MM UTC); empty = start of the run")
        ctrl.addWidget(self._cs_vid_start)
        _b = QToolButton(); _b.setText("⌖"); _b.setAutoRaise(True); _b.setToolTip("take the start from the cursor (snapshot time, or the left edge of the window)")
        _b.clicked.connect(lambda: self._cs_vid_from_cursor("start")); ctrl.addWidget(_b)
        ctrl.addWidget(QLabel("to"))
        self._cs_vid_end = QLineEdit(); self._cs_vid_end.setFixedWidth(52); self._cs_vid_end.setPlaceholderText("HH:MM")
        self._cs_vid_end.setToolTip("last frame of the video (HH:MM UTC); empty = end of the run")
        ctrl.addWidget(self._cs_vid_end)
        _b = QToolButton(); _b.setText("⌖"); _b.setAutoRaise(True); _b.setToolTip("take the end from the cursor (snapshot time, or the right edge of the window)")
        _b.clicked.connect(lambda: self._cs_vid_from_cursor("end")); ctrl.addWidget(_b)
        ctrl.addWidget(QLabel("duration:"))
        self._cs_duration_combo = QComboBox()
        self._cs_duration_combo.addItems(["5 s", "10 s", "15 s", "30 s", "45 s", "60 s"])
        self._cs_duration_combo.setCurrentIndex(2)  # 15s default
        self._cs_duration_combo.setToolTip(
            "Wall-clock duration of the exported MP4/GIF.\n"
            "The full experiment is mapped to this duration."
        )
        ctrl.addWidget(self._cs_duration_combo)

        ctrl.addWidget(QLabel("fps:"))
        self._cs_fps_combo = QComboBox()
        self._cs_fps_combo.addItems(["10", "15", "20", "24", "30"])
        self._cs_fps_combo.setCurrentIndex(1)  # 15 fps default
        self._cs_fps_combo.setToolTip(
            "Frames per second in the output file. Lower values render faster\n"
            "(fewer frames) but look choppier. 15 fps is fine for data animation."
        )
        ctrl.addWidget(self._cs_fps_combo)

        self._cs_save_video_btn = QPushButton("🎥 Save video…")
        self._cs_save_video_btn.setToolTip(
            "Render an MP4 (or GIF) of the full experiment evolution.\n"
            "Runs frame-by-frame off-screen; the GUI itself does not animate."
        )
        self._cs_save_video_btn.clicked.connect(self._cs_save_video)
        ctrl.addWidget(self._cs_save_video_btn)

        ctrl.addSpacing(10)
        self._cs_radical_lbl = QLabel("")
        self._cs_radical_lbl.setStyleSheet("color: #333; font-family: monospace; font-size: 10px;")
        # a QLabel never shrinks below its text: with the long radical readout the control row pushed
        # the whole tab beyond the screen.  Let the layout clip it instead.
        self._cs_radical_lbl.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Preferred)
        self._cs_radical_lbl.setMinimumWidth(0)
        ctrl.addWidget(self._cs_radical_lbl, 1)
        outer_layout.addLayout(ctrl, 0)

        # Tab-level widget for video frame grabbing (grabs outer + toolbar)
        self._cs_tab_widget = outer

        # ── 2×2 grid of plots ─────────────────────────────────────────────
        grid_widget = QWidget()
        from PyQt6.QtWidgets import QGridLayout
        grid = QGridLayout(grid_widget)
        grid.setContentsMargins(0, 0, 0, 0)
        grid.setSpacing(2)

        # Top-left: injection timeseries with slider
        self._cs_ts_fig = Figure(figsize=(8, 3))
        self._cs_ts_canvas = FigureCanvas(self._cs_ts_fig)
        grid.addWidget(self._cs_ts_canvas, 0, 0)

        # Top-right: VBS
        self._cs_vbs_fig = Figure(figsize=(8, 6))
        self._cs_vbs_canvas = FigureCanvas(self._cs_vbs_fig); _interact(self._cs_vbs_canvas, hover=False)
        grid.addWidget(self._cs_vbs_canvas, 0, 1)

        # Bottom-left: C/O distribution
        self._cs_co_fig = Figure(figsize=(8, 4))
        self._cs_co_canvas = FigureCanvas(self._cs_co_fig); _interact(self._cs_co_canvas, hover=False)
        grid.addWidget(self._cs_co_canvas, 1, 0)

        # Bottom-right: family pie
        self._cs_fam_fig = Figure(figsize=(8, 4))
        self._cs_fam_canvas = FigureCanvas(self._cs_fam_fig); _interact(self._cs_fam_canvas, hover=False)
        grid.addWidget(self._cs_fam_canvas, 1, 1)

        # Stretch: top-left is 1/4 of space (small), top-right is dominant
        grid.setColumnStretch(0, 1)
        grid.setColumnStretch(1, 1)
        grid.setRowStretch(0, 1)
        grid.setRowStretch(1, 1)
        outer_layout.addWidget(grid_widget, 1)

        # ── Initial draw ──────────────────────────────────────────────────
        # Default time window: full experiment range
        t = self.df.index.values.astype(float)
        self._cs_time_range = (float(t[0]), float(t[-1]))
        self._cs_time_point = float(t[len(t)//2])

        # Interactive state
        self._cs_drag = None   # 'left', 'right', 'point', 'pan'
        self._cs_pan_anchor = None
        self._cs_vbs_annot = None
        self._cs_md_annot = None

        # Fallback video fps (user can pick in the dropdown)
        self._cs_video_fps = 15

        # Wire up controls
        self._cs_mode_combo.currentIndexChanged.connect(self._cs_redraw_all)
        self._cs_units_combo.currentIndexChanged.connect(self._cs_redraw_all)
        self._cs_size_combo.currentIndexChanged.connect(self._cs_redraw_all)

        # Mouse interactions on the timeseries panel
        self._cs_ts_canvas.mpl_connect('button_press_event',   self._cs_on_press)
        self._cs_ts_canvas.mpl_connect('motion_notify_event',  self._cs_on_motion)
        self._cs_ts_canvas.mpl_connect('button_release_event', self._cs_on_release)

        # Hover tooltips on VBS and mass defect plots
        self._cs_vbs_canvas.mpl_connect(
            'motion_notify_event', lambda e: self._cs_on_scatter_hover(e, 'vbs'))
        self._cs_fam_canvas.mpl_connect(
            'motion_notify_event', lambda e: self._cs_on_scatter_hover(e, 'md'))

        self.right_tabs.addTab(outer, "Chemical Space")
        self._cs_redraw_all()

    # ── Drawing routines ──────────────────────────────────────────────────

    def _cs_draw_injection_ts(self):
        """Top-left panel: injected species timeseries + roof-open shading + time cursor/range."""
        self._cs_ts_fig.clear()
        ax = self._cs_ts_fig.add_subplot(111)

        cfg = self.budget_engine.cfg if (self.budget_engine and self.budget_engine.cfg) else {}
        bursts = cfg.get('burst_injections', [])
        # Collect injected target species present in df
        inj_species = []
        for b in bursts:
            tgt = (b.get('target') or '').upper()
            if tgt and tgt in self.df.columns and tgt not in inj_species:
                inj_species.append(tgt)
        # Limit to the top few by peak concentration to avoid clutter
        if len(inj_species) > 6:
            peaks = [(sp, np.nanmax(self.df[sp].values / self.df['M'].values * 1e9)) for sp in inj_species]
            peaks.sort(key=lambda x: -x[1])
            inj_species = [sp for sp, _ in peaks[:6]]

        t = self.df.index.values.astype(float)

        # Roof-open shading (gold)
        hv_col = next((c for c in self.df.columns if c.upper() in ('SAPHIR_HV', 'HV')), None)
        if hv_col is not None:
            hv = self.df[hv_col].values.astype(float)
            in_open, start = False, None
            first_span = True
            for ti, vi in zip(t, hv):
                if vi > 0 and not in_open:
                    in_open, start = True, ti
                elif vi <= 0 and in_open:
                    ax.axvspan(start, ti, color='gold', alpha=0.15, zorder=0,
                               label='Roof open' if first_span else None)
                    first_span = False
                    in_open = False
            if in_open:
                ax.axvspan(start, t[-1], color='gold', alpha=0.15, zorder=0,
                           label='Roof open' if first_span else None)

        # Injection start markers
        first_inj = True
        for b in bursts:
            start_utc = b.get('start_utc')
            if not start_utc:
                continue
            try:
                ts = pd.to_datetime(start_utc, dayfirst=True)
                t_inj = (ts.timestamp() - self._midnight_unix) / 3600.0
                ax.axvline(t_inj, ls=':', color='firebrick', alpha=0.8, lw=1.2,
                           label='Injection' if first_inj else None, zorder=1)
                first_inj = False
            except Exception:
                pass

        # Plot the injected species
        if inj_species:
            if 'M' in self.df.columns:
                M = self.df['M'].values
                for sp in inj_species:
                    ppb = self.df[sp].values / M * 1e9
                    ax.plot(t, ppb, lw=1.3, label=sp, zorder=2)
            ax.set_ylabel("ppb", fontsize=9)
        else:
            ax.text(0.5, 0.5, "No burst_injections in config.\n"
                    "Drag on the x-axis to select a time window.",
                    transform=ax.transAxes, ha='center', va='center', fontsize=9, color='grey')

        # Time selector overlay (current window or point)
        t_lo, t_hi = self._cs_time_range
        mode_window = self._cs_mode_combo.currentIndex() == 1
        if mode_window:
            ax.axvspan(t_lo, t_hi, color='steelblue', alpha=0.18, zorder=3,
                       label='Selected window')
            ax.axvline(t_lo, color='steelblue', lw=1.8, zorder=4)
            ax.axvline(t_hi, color='steelblue', lw=1.8, zorder=4)
        else:
            ax.axvline(self._cs_time_point, color='black', lw=1.8, ls='--', zorder=4,
                       label='Selected time')

        ax.xaxis.set_major_formatter(FuncFormatter(_fmt_hhmm))
        ax.set_xlabel("Time of day (UTC)", fontsize=9)
        ax.set_title("Injections & roof status   (drag to select time)", fontsize=10, fontweight='bold')
        ax.grid(True, ls='--', alpha=0.35)
        if ax.get_legend_handles_labels()[0]:
            ax.legend(loc='upper right', fontsize=7, framealpha=0.85).set_draggable(True)
        # Fixed margins (faster than tight_layout, which re-computes every call)
        self._cs_ts_fig.subplots_adjust(left=0.08, right=0.98, top=0.88, bottom=0.18)
        self._cs_ts_canvas.draw()

    def _cs_global_conc_range(self, unit_key):
        """Return (v_min_nonzero, v_max) of the selected unit across ALL
        organic species and ALL time points. Used so that marker sizes in
        VBS and mass-defect plots reflect *absolute* concentration rather
        than normalizing to the current snapshot's maximum (which would
        pin the dominant species at maximum size regardless of decay).
        Cached on self per unit key; invalidate if the dataframe changes."""
        cache_attr = f'_cs_global_range_{unit_key}'
        if hasattr(self, cache_attr):
            return getattr(self, cache_attr)
        spc_all = self._cs_build_species_table()
        spc_all = spc_all[~spc_all['is_radical']]
        if spc_all.empty or 'M' not in self.df.columns:
            rv = (1e-6, 1.0)
            setattr(self, cache_attr, rv)
            return rv
        M = self.df['M'].values.astype(float)
        vmin, vmax = np.inf, -np.inf
        for _, row in spc_all.iterrows():
            name = row['name']
            if name not in self.df.columns:
                continue
            with np.errstate(divide='ignore', invalid='ignore'):
                ppb = np.where(M > 0, self.df[name].values / M * 1e9, np.nan)
            vals = ppb if unit_key == 'ppb' else ppb * float(row['mass']) / 24.45
            valid = vals[np.isfinite(vals) & (vals > 0)]
            if valid.size:
                vmin = min(vmin, float(valid.min()))
                vmax = max(vmax, float(valid.max()))
        if not np.isfinite(vmin) or not np.isfinite(vmax) or vmax <= 0:
            rv = (1e-6, 1.0)
        else:
            rv = (vmin, vmax)
        setattr(self, cache_attr, rv)
        return rv

    def _cs_compute_sizes(self, vals):
        """Map concentrations (in currently selected unit) to scatter marker
        areas [10, 500] pt².

        Uses a global dynamic range anchored at the experiment max, and
        clamps the floor to four orders of magnitude below that. This
        choice matches typical chamber physics (real concentration changes
        span 3–4 orders; anything smaller is numerical noise), and keeps
        a 3–10× concentration change visibly different.

        Log mode: size ∝ clipped log10(conc).
        Linear mode: size ∝ sqrt(conc) so radius is proportional to
        concentration and area is proportional to conc (nicer perception
        of small changes than the raw conc mapping)."""
        unit_key = 'ppb' if self._cs_units_combo.currentIndex() == 0 else 'ugm3'
        _, v_max = self._cs_global_conc_range(unit_key)
        vals = np.asarray(vals, dtype=float)
        # 4-decade dynamic range below the experiment max
        DYN_RANGE = 4.0
        v_floor = max(v_max * (10.0 ** (-DYN_RANGE)), 1e-30)
        if self._cs_size_combo.currentText() == 'log10':
            s_raw = np.log10(np.maximum(vals, v_floor))
            s_lo  = np.log10(v_floor)
            s_hi  = np.log10(max(v_max, v_floor * 10))
            s_norm = (s_raw - s_lo) / max(s_hi - s_lo, 1e-9)
        else:
            # Linear mode with sqrt scaling → area ∝ conc but the scale is
            # compressed enough that mid-range species are still readable.
            ratio = np.clip(vals / max(v_max, 1e-30), 0.0, 1.0)
            s_norm = np.sqrt(ratio)
        s_norm = np.clip(s_norm, 0.0, 1.0)
        return 10.0 + 490.0 * s_norm

    def _cs_species_with_conc(self):
        """Species table + per-species concentration (ppb and µg/m³) at the
        current cursor time (snapshot) or averaged over the window.

        All four draw methods call this every frame, so we cache the result
        keyed on (mode, cursor-time-or-range, units). Without the cache a
        rendered video would recompute it 4× per frame for ~200 species."""
        mode_idx = self._cs_mode_combo.currentIndex() if hasattr(self, '_cs_mode_combo') else 0
        if mode_idx == 0:
            time_key = round(float(self._cs_time_point), 6)
        else:
            lo, hi = self._cs_current_window()
            time_key = (round(float(lo), 6), round(float(hi), 6))
        unit_idx = self._cs_units_combo.currentIndex() if hasattr(self, '_cs_units_combo') else 0
        excl_txt = self._cs_exclude.text() if hasattr(self, '_cs_exclude') else ''
        cache_key = (mode_idx, time_key, unit_idx, excl_txt)
        if getattr(self, '_cs_spc_cache_key', None) == cache_key:
            return self._cs_spc_cache

        spc = self._cs_build_species_table()
        excl = {x.strip().upper() for x in excl_txt.replace(';', ',').split(',') if x.strip()}
        if excl and not spc.empty:
            drop = spc['name'].str.upper().isin(excl) | spc['formula'].str.upper().isin(excl)
            spc = spc[~drop].copy()
        if spc.empty:
            self._cs_spc_cache_key = cache_key
            self._cs_spc_cache = spc
            return spc
        spc = spc[~spc['is_radical']].copy()
        conc_ppb = np.array([self._cs_avg_conc_ppb(n) for n in spc['name']])
        conc_ug  = conc_ppb * spc['mass'].values / 24.45   # ppb → µg/m³ at 1 atm, 298 K
        spc['conc_ppb']  = conc_ppb
        spc['conc_ugm3'] = conc_ug
        unit = 'conc_ppb' if unit_idx == 0 else 'conc_ugm3'
        spc = spc[np.isfinite(spc[unit]) & (spc[unit] > 1e-6)].copy()
        spc['unit_val'] = spc[unit]
        self._cs_spc_cache_key = cache_key
        self._cs_spc_cache = spc
        return spc

    def _cs_draw_vbs(self):
        """Top-right panel: LogC* vs OSc, circles sized/coloured by concentration."""
        self._cs_vbs_fig.clear()
        ax = self._cs_vbs_fig.add_subplot(111)

        spc = self._cs_species_with_conc()
        if spc.empty or not np.any(np.isfinite(spc['logc'])):
            ax.text(0.5, 0.5, "No organic species with valid LogC* in this time window.",
                    ha='center', va='center', transform=ax.transAxes, color='grey')
            ax.set_axis_off()
            self._cs_vbs_canvas.draw()
            return

        # Volatility class background shading
        for (lo, hi), color in [
            ((-3.5, -0.51), 'dimgrey'),
            ((-0.5,  2.49), 'darkgrey'),
            (( 2.5,  6.49), 'lightgrey'),
        ]:
            ax.axvspan(lo, hi, color=color, alpha=0.18, zorder=0)
        for x, label in [(-2, 'LVOC'), (1, 'SVOC'), (4.5, 'IVOC'), (8.5, 'VOC')]:
            ax.text(x, 2.25, label, ha='center', va='bottom', fontsize=9, color='#444')

        # Size scaling — global (across all species, all times) so markers
        # shrink visibly as concentrations decay instead of staying pinned
        # at the per-frame maximum.
        sizes = self._cs_compute_sizes(spc['unit_val'].values)

        # Color by class
        colors = [self._CS_CLASS_COLORS.get(c, '#7f7f7f') for c in spc['class']]
        sc = ax.scatter(spc['logc'], spc['osc'], s=sizes, c=colors,
                        alpha=0.55, edgecolors='k', linewidths=0.4, zorder=2)
        self._cs_vbs_scatter = sc
        self._cs_vbs_species = spc.reset_index(drop=True)

        ax.set_xlim(-4, 11)
        ax.set_ylim(-2.5, 2.5)
        ax.set_xticks(range(-3, 12))
        ax.set_xticklabels([str(x) if x % 2 == 0 else '' for x in range(-3, 12)])
        ax.set_xlabel(r"Volatility  $\log_{10}(C^*)$  [µg m$^{-3}$]", fontsize=9)
        ax.set_ylabel(r"Carbon Oxidation State  $\overline{\mathrm{OS}}_\mathrm{C}$", fontsize=9)
        unit_lbl = 'ppb' if self._cs_units_combo.currentIndex() == 0 else 'µg/m³'
        ax.set_title(f"Volatility basis set  (N={len(spc)} species, sized by {unit_lbl})",
                     fontsize=10, fontweight='bold')

        # Legend for classes
        from matplotlib.lines import Line2D
        legend_handles = [Line2D([0],[0], marker='o', color='w',
                                 markerfacecolor=self._CS_CLASS_COLORS[c],
                                 markersize=8, label=c)
                          for c in self._CS_CLASS_ORDER if c in spc['class'].values]
        if legend_handles:
            ax.legend(handles=legend_handles, loc='lower right', fontsize=8,
                      framealpha=0.85).set_draggable(True)

        ax.grid(True, ls='--', alpha=0.3)
        ax.text(0.01, -0.14, "Li et al., ACP 16, 3327 (2016) · Kroll et al., Nat. Chem. 3, 133 (2011)",
                transform=ax.transAxes, fontsize=6, color='#666', style='italic')

        self._cs_vbs_fig.subplots_adjust(left=0.10, right=0.96, top=0.90, bottom=0.14)
        self._cs_vbs_canvas.draw()

    def _cs_draw_c_o_distribution(self):
        """Bottom-left panel: two stacked bar plots, both with nC on the x-axis.
        Left subplot: bars stacked by oxygen number (nO).
        Right subplot: bars stacked by CHO family class.
        Matches the NNMF project convention for compound composition breakdown."""
        self._cs_co_fig.clear()
        gs = self._cs_co_fig.add_gridspec(1, 2, wspace=0.3)
        ax_no = self._cs_co_fig.add_subplot(gs[0])  # stacked by nO
        ax_cl = self._cs_co_fig.add_subplot(gs[1])  # stacked by class

        spc = self._cs_species_with_conc()
        if spc.empty:
            for ax in (ax_no, ax_cl):
                ax.text(0.5, 0.5, "No data in this time window.", ha='center', va='center',
                        transform=ax.transAxes, color='grey')
                ax.set_axis_off()
            self._cs_co_canvas.draw()
            return

        unit_lbl = 'ppb' if self._cs_units_combo.currentIndex() == 0 else 'µg/m³'
        x_all = sorted(spc['nC'].unique())

        # ── LEFT: nC vs Σconc, stacked by nO ─────────────────────────
        grp_no = (spc.groupby(['nC', 'nO'])['unit_val'].sum().unstack(fill_value=0))
        grp_no = grp_no.reindex(x_all, fill_value=0)
        n_levels = sorted(grp_no.columns)
        cmap = _mcm.get_cmap('viridis')
        bottoms = np.zeros(len(x_all))
        for no in n_levels:
            color = cmap(no / max(max(n_levels), 1))
            vals = grp_no[no].values
            ax_no.bar(x_all, vals, bottom=bottoms, color=color,
                      edgecolor='k', linewidth=0.3, label=f"nO={int(no)}")
            bottoms = bottoms + vals
        ax_no.set_xlabel("Number of carbon atoms (nC)", fontsize=9)
        ax_no.set_ylabel(f"Σ concentration ({unit_lbl})", fontsize=9)
        ax_no.set_title("Stacked by oxygen number", fontsize=9)
        ax_no.set_xticks(x_all)
        ax_no.grid(True, axis='y', ls='--', alpha=0.35)
        if n_levels:
            ax_no.legend(fontsize=6, loc='upper right', ncol=2, framealpha=0.85,
                         title='nO', title_fontsize=7).set_draggable(True)

        # ── RIGHT: nC vs Σconc, stacked by CHO family class ──────────
        grp_cl = (spc.groupby(['nC', 'class'])['unit_val'].sum().unstack(fill_value=0))
        grp_cl = grp_cl.reindex(x_all, fill_value=0)
        bottoms = np.zeros(len(x_all))
        for cls in self._CS_CLASS_ORDER:
            if cls not in grp_cl.columns:
                continue
            vals = grp_cl[cls].values
            ax_cl.bar(x_all, vals, bottom=bottoms,
                      color=self._CS_CLASS_COLORS[cls],
                      edgecolor='k', linewidth=0.3, label=cls)
            bottoms = bottoms + vals
        ax_cl.set_xlabel("Number of carbon atoms (nC)", fontsize=9)
        ax_cl.set_ylabel(f"Σ concentration ({unit_lbl})", fontsize=9)
        ax_cl.set_title("Stacked by CHO family", fontsize=9)
        ax_cl.set_xticks(x_all)
        ax_cl.grid(True, axis='y', ls='--', alpha=0.35)
        if ax_cl.get_legend_handles_labels()[0]:
            ax_cl.legend(fontsize=7, loc='upper right', framealpha=0.85).set_draggable(True)

        self._cs_co_fig.suptitle("Composition by carbon number",
                                  fontsize=10, fontweight='bold', y=0.99)
        self._cs_co_fig.subplots_adjust(left=0.07, right=0.98, top=0.88, bottom=0.16, wspace=0.28)
        self._cs_co_canvas.draw()

    def _cs_draw_family_pie(self):
        """Bottom-right panel: mass defect plot — molecular mass vs mass defect
        (mass − round(mass)). Circles sized and coloured like the VBS panel.
        Classic tool for distinguishing CH/CHO/CHON families by their intrinsic
        hydrogen/oxygen content relative to integer mass."""
        self._cs_fam_fig.clear()
        ax = self._cs_fam_fig.add_subplot(111)

        spc = self._cs_species_with_conc()
        if spc.empty or not np.any(np.isfinite(spc['mass'])):
            ax.text(0.5, 0.5, "No data in this time window.", ha='center', va='center',
                    transform=ax.transAxes, color='grey')
            ax.set_axis_off()
            self._cs_fam_canvas.draw()
            return

        # Mass defect: exact mass − nominal integer mass
        masses = spc['mass'].values
        defect = masses - np.round(masses)

        # Size scaling — same global rule as VBS
        sizes = self._cs_compute_sizes(spc['unit_val'].values)

        # Color by family class (consistent with other panels)
        colors = [self._CS_CLASS_COLORS.get(c, '#7f7f7f') for c in spc['class']]
        sc = ax.scatter(masses, defect, s=sizes, c=colors,
                        alpha=0.55, edgecolors='k', linewidths=0.4, zorder=2)

        # Cache scatter for hover tooltips
        self._cs_md_scatter = sc
        self._cs_md_species = spc.reset_index(drop=True)

        # Reference line at zero defect
        ax.axhline(0.0, color='k', lw=0.6, ls='--', alpha=0.5, zorder=1)

        unit_lbl = 'ppb' if self._cs_units_combo.currentIndex() == 0 else 'µg/m³'
        ax.set_xlabel("Molecular mass  (g/mol)", fontsize=9)
        ax.set_ylabel("Mass defect  (m − round(m))", fontsize=9)
        ax.set_title(f"Mass defect plot  (N={len(spc)} species, sized by {unit_lbl})",
                     fontsize=10, fontweight='bold')
        ax.grid(True, ls='--', alpha=0.3)

        from matplotlib.lines import Line2D
        legend_handles = [Line2D([0],[0], marker='o', color='w',
                                 markerfacecolor=self._CS_CLASS_COLORS[c],
                                 markersize=8, label=c)
                          for c in self._CS_CLASS_ORDER if c in spc['class'].values]
        if legend_handles:
            ax.legend(handles=legend_handles, loc='lower right', fontsize=8,
                      framealpha=0.85).set_draggable(True)

        self._cs_fam_fig.subplots_adjust(left=0.10, right=0.96, top=0.90, bottom=0.14)
        self._cs_fam_canvas.draw()

    def _cs_update_radical_panel(self):
        """Update the radical-averages label in the control row."""
        lines = []
        for sp, lbl, unit in [('OH',   'OH',   'molec/cm³'),
                              ('HO2',  'HO2',  'molec/cm³'),
                              ('O3',   'O3',   'ppb'),
                              ('RO2',  'RO2',  'ppb'),
                              ('NO',   'NO',   'ppb'),
                              ('NO2',  'NO2',  'ppb')]:
            if sp not in self.df.columns:
                continue
            if unit == 'molec/cm³':
                v = self._cs_avg_conc_molec(sp)
                if np.isfinite(v):
                    lines.append(f"{lbl}={v:.2e} {unit}")
            else:
                v = self._cs_avg_conc_ppb(sp)
                if np.isfinite(v):
                    lines.append(f"{lbl}={v:.3g} {unit}")
        # Show time window context
        t_lo, t_hi = self._cs_current_window()
        mode = self._cs_mode_combo.currentText()
        if mode.startswith("Snapshot"):
            header = f"@ t={self._cs_time_point:.2f}h  |  "
        else:
            header = f"window {t_lo:.2f}–{t_hi:.2f}h  |  "
        self._cs_radical_lbl.setText(header + "  ".join(lines))

    def _cs_redraw_all(self):
        """Redraw all four panels from scratch (call after any state change)."""
        try:
            self._cs_draw_injection_ts()
            self._cs_draw_vbs()
            self._cs_draw_c_o_distribution()
            self._cs_draw_family_pie()
            self._cs_update_radical_panel()
        except Exception as e:
            import traceback
            traceback.print_exc()
            QMessageBox.warning(self, "Chemical Space", f"Redraw failed:\n{e}")

    # ── Mouse interaction on timeseries ──────────────────────────────────

    # ── Video export ────────────────────────────────────────────────────
    # Live playback has been removed (it crashed the GUI on Windows and was
    # slow due to the cost of redrawing four matplotlib figures per frame).
    # Instead: render the animation off-screen to an MP4/GIF file and offer
    # to open it in the system default player.

    @staticmethod
    def _cs_parse_duration(s):
        """'30 s' / '1 min' / '2 min' → wall-clock seconds."""
        s = s.strip().lower().replace(' ', '')
        if s.endswith('min'):
            return float(s[:-3]) * 60.0
        if s.endswith('s'):
            return float(s[:-1])
        try:
            return float(s)
        except ValueError:
            return 15.0

    def _cs_video_range(self):
        """(start, end) in hours for the video: the 'Video from … to …' fields, each falling back to the
        run's first/last time; clamped to the run."""
        t = self.df.index.values.astype(float); lo, hi = float(t[0]), float(t[-1])

        def parse(txt, default):
            txt = (txt or "").strip()
            if not txt:
                return default
            try:
                if ":" in txt:
                    h, m = txt.split(":", 1); v = int(h) + int(m) / 60.0
                else:
                    v = float(txt)
            except ValueError:
                return default
            return min(max(v, lo), hi)
        return parse(self._cs_vid_start.text(), lo), parse(self._cs_vid_end.text(), hi)

    def _cs_vid_from_cursor(self, which: str):
        """Fill the start or end field from the cursor: the snapshot time, or the window's edge."""
        if self._cs_mode_combo.currentIndex() == 1:
            v = self._cs_time_range[0 if which == "start" else 1]
        else:
            v = self._cs_time_point
        txt = f"{int(v):02d}:{int(round((v - int(v)) * 60)) % 60:02d}"
        (self._cs_vid_start if which == "start" else self._cs_vid_end).setText(txt)

    def _cs_save_video(self):
        """Render the full experiment evolution to an MP4 (or GIF) file.
        The tab widget is grabbed off-screen frame by frame; the GUI itself
        does not animate. A QProgressDialog provides feedback and allows
        cancellation. After rendering, the user is offered to open the file
        in the system default player."""
        try:
            import imageio.v2 as imageio
        except ImportError:
            QMessageBox.warning(self, "Missing package",
                "Saving video requires the `imageio` package.\n"
                "Install with:  pip install imageio imageio-ffmpeg")
            return

        # Output file
        start_dir = (os.path.dirname(self.csv_path) if self.csv_path else os.getcwd())
        default_name = "chemspace_animation.mp4"
        out_path, _ = QFileDialog.getSaveFileName(
            self, "Save animation", os.path.join(start_dir, default_name),
            "MP4 video (*.mp4);;GIF animation (*.gif)"
        )
        if not out_path:
            return

        # Build frame timeline — user-chosen wall-clock duration & fps
        duration_s = self._cs_parse_duration(self._cs_duration_combo.currentText())
        try:
            fps = int(self._cs_fps_combo.currentText())
        except Exception:
            fps = self._cs_video_fps
        n_frames = max(int(duration_s * fps), 2)
        t0, t1 = self._cs_video_range()
        if t1 <= t0:
            QMessageBox.warning(self, "Video range", "The end time must be after the start time."); return
        frame_times = np.linspace(t0, t1, n_frames)

        # Remember state so we can restore after rendering
        saved_mode = self._cs_mode_combo.currentIndex()
        saved_tp   = self._cs_time_point
        # Snapshot mode for rendering
        if saved_mode != 0:
            self._cs_mode_combo.setCurrentIndex(0)

        # Choose writer
        is_gif = out_path.lower().endswith('.gif')
        try:
            if is_gif:
                writer = imageio.get_writer(out_path, mode='I',
                                            duration=1.0 / fps, loop=0)
            else:
                writer = imageio.get_writer(out_path, fps=fps,
                                            codec='libx264', quality=8,
                                            pixelformat='yuv420p',
                                            macro_block_size=2)
        except Exception as e:
            QMessageBox.critical(self, "Video writer failed",
                f"Could not open {out_path} for writing:\n{e}\n\n"
                "For MP4 output, ensure `imageio-ffmpeg` is installed:\n"
                "  pip install imageio-ffmpeg")
            return

        from PyQt6.QtWidgets import QProgressDialog

        # Bypass Qt's widget.grab() entirely — it was crashing on repeated calls
        # (likely QImage buffer ownership + Python GC racing). Read pixel data
        # straight out of each matplotlib canvas's Agg buffer and composite the
        # four panels into one frame with numpy/PIL. No Qt image objects touched.
        try:
            from PIL import Image as PILImage
        except ImportError:
            QMessageBox.warning(self, "Missing package",
                "Saving video requires Pillow (installed with matplotlib).")
            return

        # Fixed output layout: columns 720 px wide, top row 540 px tall,
        # bottom row 420 px tall → 1440 × 960 frame.
        COL_W, ROW_H_TOP, ROW_H_BOT = 720, 540, 420

        def _fig_rgb(fig):
            """Return the figure's Agg buffer as a contiguous (H,W,3) uint8 array."""
            fig.canvas.draw()
            w, h = fig.canvas.get_width_height()
            # buffer_rgba is a memoryview; wrap it and copy to own memory
            buf = np.asarray(fig.canvas.buffer_rgba()).reshape(h, w, 4)
            return np.ascontiguousarray(buf[:, :, :3])

        def _resize_rgb(arr, w, h):
            if arr.shape[0] == h and arr.shape[1] == w:
                return arr
            img = PILImage.fromarray(arr).resize((w, h), PILImage.BILINEAR)
            return np.asarray(img, dtype=np.uint8)

        def _composite_frame():
            ts  = _resize_rgb(_fig_rgb(self._cs_ts_fig),  COL_W, ROW_H_TOP)
            vbs = _resize_rgb(_fig_rgb(self._cs_vbs_fig), COL_W, ROW_H_TOP)
            co  = _resize_rgb(_fig_rgb(self._cs_co_fig),  COL_W, ROW_H_BOT)
            fam = _resize_rgb(_fig_rgb(self._cs_fam_fig), COL_W, ROW_H_BOT)
            top = np.hstack([ts, vbs])
            bot = np.hstack([co, fam])
            return np.vstack([top, bot])

        # Progress dialog
        progress = QProgressDialog(f"Rendering {n_frames} frames…",
                                   "Cancel", 0, n_frames, self)
        progress.setWindowTitle("Saving video")
        progress.setMinimumDuration(0)
        progress.setAutoClose(True)

        rendered_ok = False
        cancelled = False
        frame_w = COL_W * 2
        frame_h = ROW_H_TOP + ROW_H_BOT
        try:
            for i, tp in enumerate(frame_times):
                if progress.wasCanceled():
                    cancelled = True
                    break
                self._cs_time_point = float(tp)
                self._cs_redraw_all()
                arr = _composite_frame()
                writer.append_data(arr)
                progress.setValue(i + 1)
                # Let the progress dialog paint + stay responsive to Cancel.
                # Keep this OUTSIDE the critical buffer-access section above
                # so nothing re-enters while Agg is reading pixels.
                QApplication.processEvents()
            if not cancelled:
                rendered_ok = True
        except Exception as e:
            import traceback
            traceback.print_exc()
            QMessageBox.critical(self, "Render error", f"Frame export failed:\n{e}")
        finally:
            try: writer.close()
            except Exception: pass
            progress.close()
            # Restore UI state
            self._cs_mode_combo.setCurrentIndex(saved_mode)
            self._cs_time_point = saved_tp
            self._cs_redraw_all()

        if cancelled:
            # Clean up the incomplete file
            try: os.remove(out_path)
            except Exception: pass
            return

        if rendered_ok:
            ans = QMessageBox.question(
                self, "Animation saved",
                f"Wrote {n_frames} frames at {fps} fps ({frame_w}×{frame_h}) to:\n\n"
                f"{out_path}\n\nOpen it now?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                QMessageBox.StandardButton.Yes
            )
            if ans == QMessageBox.StandardButton.Yes:
                try:
                    if sys.platform.startswith('win'):
                        os.startfile(out_path)
                    elif sys.platform == 'darwin':
                        import subprocess; subprocess.Popen(['open', out_path])
                    else:
                        import subprocess; subprocess.Popen(['xdg-open', out_path])
                except Exception as e:
                    QMessageBox.warning(self, "Open failed",
                        f"Could not open the file:\n{e}")

    # ── Mouse interaction on timeseries ──────────────────────────────────

    def _cs_on_press(self, event):
        if event.inaxes is None or event.xdata is None:
            return
        x = float(event.xdata)
        mode_window = self._cs_mode_combo.currentIndex() == 1
        if mode_window:
            t_lo, t_hi = self._cs_time_range
            ax = event.inaxes
            # Pixels-per-hour for click tolerance
            width_px = ax.bbox.width
            xlim = ax.get_xlim()
            hr_per_px = (xlim[1] - xlim[0]) / max(width_px, 1)
            tol = hr_per_px * 5  # 5 pixels
            # Check which edge (or inside)
            if abs(x - t_lo) < tol:
                self._cs_drag = 'left'
            elif abs(x - t_hi) < tol:
                self._cs_drag = 'right'
            elif t_lo < x < t_hi:
                self._cs_drag = 'pan'
                self._cs_pan_anchor = (x, t_lo, t_hi)
            else:
                # Click outside: move the nearest edge there
                if abs(x - t_lo) < abs(x - t_hi):
                    self._cs_time_range = (x, t_hi)
                    self._cs_drag = 'left'
                else:
                    self._cs_time_range = (t_lo, x)
                    self._cs_drag = 'right'
                self._cs_redraw_all()
        else:
            self._cs_time_point = x
            self._cs_drag = 'point'
            self._cs_redraw_all()

    def _cs_on_motion(self, event):
        if self._cs_drag is None or event.inaxes is None or event.xdata is None:
            return
        x = float(event.xdata)
        t_arr = self.df.index.values.astype(float)
        t_min, t_max = float(t_arr[0]), float(t_arr[-1])
        x = max(t_min, min(t_max, x))
        t_lo, t_hi = self._cs_time_range
        if self._cs_drag == 'left':
            t_lo = min(x, t_hi - 0.01)
            self._cs_time_range = (t_lo, t_hi)
        elif self._cs_drag == 'right':
            t_hi = max(x, t_lo + 0.01)
            self._cs_time_range = (t_lo, t_hi)
        elif self._cs_drag == 'pan':
            anchor_x, anchor_lo, anchor_hi = self._cs_pan_anchor
            dx = x - anchor_x
            width = anchor_hi - anchor_lo
            new_lo = max(t_min, anchor_lo + dx)
            new_hi = new_lo + width
            if new_hi > t_max:
                new_hi = t_max
                new_lo = new_hi - width
            self._cs_time_range = (new_lo, new_hi)
        elif self._cs_drag == 'point':
            self._cs_time_point = x
        # Light redraw during drag: just the timeseries overlay
        self._cs_draw_injection_ts()

    def _cs_on_release(self, event):
        if self._cs_drag is None:
            return
        self._cs_drag = None
        self._cs_pan_anchor = None
        # Heavy redraw at end: all four panels now reflect the new window
        self._cs_redraw_all()

    def _cs_on_scatter_hover(self, event, kind):
        """Shared hover handler for VBS ('vbs') and mass-defect ('md') plots.
        Shows a yellow tooltip with the MCM name, chemical formula, key
        coordinates, and concentration in the selected unit."""
        if kind == 'vbs':
            sc_attr, spc_attr, annot_attr, canvas = ('_cs_vbs_scatter',
                '_cs_vbs_species', '_cs_vbs_annot', self._cs_vbs_canvas)
            xy_from = lambda row: (row['logc'], row['osc'])
            coord_line = lambda row: f"LogC*={row['logc']:.2f}  OSc={row['osc']:.2f}"
        else:  # 'md'
            sc_attr, spc_attr, annot_attr, canvas = ('_cs_md_scatter',
                '_cs_md_species', '_cs_md_annot', self._cs_fam_canvas)
            xy_from = lambda row: (row['mass'], row['mass'] - round(row['mass']))
            coord_line = lambda row: (f"mass={row['mass']:.3f}  "
                                      f"defect={row['mass']-round(row['mass']):+.3f}")

        if event.inaxes is None or not hasattr(self, sc_attr):
            return
        sc = getattr(self, sc_attr)
        cont, info = sc.contains(event)
        ax = sc.axes
        # Clear previous annotation
        old = getattr(self, annot_attr, None)
        if old is not None:
            try: old.remove()
            except Exception: pass
            setattr(self, annot_attr, None)
        if cont and info.get('ind') is not None and len(info['ind']):
            idx = int(info['ind'][0])
            spc = getattr(self, spc_attr)
            if idx < len(spc):
                row = spc.iloc[idx]
                unit_lbl = 'ppb' if self._cs_units_combo.currentIndex() == 0 else 'µg/m³'
                formula = row.get('formula') or row['name']
                txt = (f"{row['name']}  ({formula})\n"
                       f"{coord_line(row)}\n"
                       f"conc={row['unit_val']:.3g} {unit_lbl}")
                annot = ax.annotate(
                    txt, xy=xy_from(row),
                    xytext=(10, 10), textcoords='offset points',
                    fontsize=8, bbox=dict(boxstyle='round,pad=0.4',
                                          fc='lightyellow', ec='k', alpha=0.9))
                setattr(self, annot_attr, annot)
        canvas.draw_idle()

    # ── Default plots on startup ──────────────────────────────────────────────
    def _plot_defaults(self):
        _script_dir = os.path.dirname(os.path.abspath(__file__))
        _tp = self._toml_path or os.path.join(_script_dir, "FR_GUI.toml")
        _ql = ['NO', 'NO2', 'O3']
        _qr = ['OH', 'HO2']
        try:
            import tomllib as _tl
            import ast as _ast
            with open(_tp, 'rb') as _fh:
                _vt = _tl.load(_fh).get('viewer', {})
            _ql_raw = _vt.get('quick_plot_left', None)
            _qr_raw = _vt.get('quick_plot_right', None)
            if _ql_raw:
                _ql = _ast.literal_eval(_ql_raw) if isinstance(_ql_raw, str) else list(_ql_raw)
            if _qr_raw:
                _qr = _ast.literal_eval(_qr_raw) if isinstance(_qr_raw, str) else list(_qr_raw)
        except Exception:
            pass

        if getattr(self, "ts_tab", None) is not None:
            self.ts_tab.apply_defaults(list(_ql), list(_qr))
            return
        for _sp in _ql:
            if _sp in self.df.columns:
                self.plot_variable(_sp, "Left Axis (Main)", "Convert to ppb")
        for _sp in self.df.columns:
            if _sp.upper() in {'SAPHIR_NO', 'SAPHIR_NO2', 'SAPHIR_O3'}:
                self.plot_variable(_sp, "Left Axis (Main)", "Native (molec/cm³ or SAPHIR unit)")
        for _sp in _qr:
            if _sp in self.df.columns:
                self.plot_variable(_sp, "Right Axis (Main)", "Native (molec/cm³ or SAPHIR unit)")

        if any(sp in self.df.columns for sp in _qr):
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
def show_results(out_path, toml_path=None):
    if HAS_PYQT:
        print(">>> [VISUALIZATION] Launching Interactive PyQt Viewer...")
        global _LOADER
        if _LOADER is None:
            try:
                from loading_box import LoadingBox
                _LOADER = LoadingBox.start("opening results")
            except Exception:
                _LOADER = None
        app = QApplication.instance()
        if app is None:
            app = QApplication(sys.argv)
        try:
            from juice_theme import apply_theme
            apply_theme(app)                       # same look as the step window
        except Exception as _e:
            print(f"[viewer] theme not applied: {_e}")
        try:
            _mb = os.path.getsize(out_path) / 1e6
        except OSError:
            _mb = 0.0
        _progress(f"Reading {os.path.basename(out_path)} ({_mb:.0f} MB)", 5,
                  "every species at every output time is loaded into memory so that all tabs respond instantly")
        df = load_column_file(out_path)
        _progress("Building the viewer", 20)
        viewer = ResultsViewer(df, csv_path=out_path, toml_path=toml_path)
        viewer.show()
        try:
            from juice_theme import force_taskbar_icon
            force_taskbar_icon(viewer)
        except Exception:
            pass
        if _LOADER is not None:
            _LOADER.done(); _LOADER = None
        app.exec()
        return

    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("[WARNING] Missing matplotlib. Visualization skipped.")
        return

    print("[WARNING] PyQt6 is not installed. Falling back to simple Matplotlib visualization.")
    df = load_column_file(out_path)
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
    parser.add_argument("--toml", default=None,
                        help="Path to active TOML config (sets quick_plot defaults).")
    args = parser.parse_args()

    csv_path = args.results_csv_path
    if csv_path and not os.path.exists(csv_path):
        print(f"[ERROR] File not found: {csv_path}")
        csv_path = None

    # No path (or bad path): open a file picker if PyQt is available
    if not csv_path and HAS_PYQT:
        app = QApplication.instance() or QApplication(sys.argv)
        scenarios_root = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                      "MCM_tested_scenarios")
        start_dir = scenarios_root if os.path.isdir(scenarios_root) else os.getcwd()
        csv_path, _ = QFileDialog.getOpenFileName(
            None, "Select simulation_results CSV", start_dir,
            "Simulation results (simulation_results*.csv);;All CSV (*.csv);;All files (*)"
        )
        if not csv_path:
            print(">>> [VIEWER] No file selected. Exiting.")
            return

    if csv_path and os.path.exists(csv_path):
        print(f">>> [VIEWER] Loading results from: {csv_path}")
        show_results(csv_path, toml_path=args.toml)
    else:
        print("\nUsage: python results_viewer.py <path_to_simulation_results.csv>")
        print(">>> [VIEWER] No path provided. Exiting.")


if __name__ == "__main__":
    if HAS_PYQT:
        _app = QApplication.instance() or QApplication(sys.argv)
        try:
            from juice_theme import apply_theme
            apply_theme(_app)
        except Exception as _e:
            print(f"[viewer] theme not applied: {_e}")
            _pm = QPixmap(os.path.join(os.path.dirname(os.path.abspath(__file__)), "JUICE Logo.png"))
            if not _pm.isNull():
                _app.setWindowIcon(QIcon(_pm))
    main()
