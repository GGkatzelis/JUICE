"""
enz_builder.py — SAPHIR Environment File (.enz) Builder

Generates kammer.enz from:
  • SAPHIR archive NetCDF files (auto-located by date)
  • Computed quantities (SZ — solar zenith angle, pure Python, no extra deps)
  • User-defined step-function protocol (HV roof times, injections)
  • Wiki-parsed events (optional; requires VPN to saphir.fz-juelich.de)
"""

import os
import re
import math
from datetime import datetime, timedelta, timezone

import numpy as np
from PyQt6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QGridLayout, QGroupBox, QSplitter,
    QLabel, QLineEdit, QPushButton, QTableWidget, QTableWidgetItem,
    QRadioButton, QWidget, QFileDialog, QMessageBox,
    QCheckBox, QHeaderView, QAbstractItemView, QSizePolicy,
)
from PyQt6.QtCore import Qt, QThread, pyqtSignal
from PyQt6.QtGui import QColor

# ── SAPHIR fixed constants ─────────────────────────────────────────────────────
SAPHIR_LAT = 50.9086       # degrees N, Jülich Research Centre
SAPHIR_LON = 6.3639        # degrees E
SAPHIR_ARCHIVE_ROOT = r"Z:\IEK8-SAPHIR\data"

# ENZ time epoch: seconds since 2000-01-01 00:00:00 UTC
ENZ_EPOCH = datetime(2000, 1, 1, 0, 0, 0, tzinfo=timezone.utc)

# ── Archive variable mapping ───────────────────────────────────────────────────
# enz_column → (nc_filename_suffix, variable_name_inside_nc, fill_value)
# Filename in archive: Z:\...\YYYY\MM\DD\YYYY-MM-DD.<suffix>
ARCHIVE_MAP = {
    'dilute':    ('SAPHIR.Flux.nc',                               'TOTAL_FLOW',  0.0),
    't':         ('SAPHIR.USA.TEMP.nc',                           'TK_M_SAPHIR', 284.0),
    'H2O':       ('SAPHIR.PICARRO_CFKADS2270.CO_CO2_CH4_H2O.nc', 'H2O',         0.0),
    'NO2':       ('SAPHIR.TR43.NO_NO2_NOx.nc',                    'NO2',         0.0),
    'NO':        ('SAPHIR.TR43.NO_NO2_NOx.nc',                    'NO',          0.0),
    'O3':        ('SAPHIR.TwoBTech_1014.O3.nc',                   'O3',          0.0),
    'jO1D':      ('SAPHIR.SR.jo1d.nc',                            'jO1D',        0.0),
    'jNO2':      ('SAPHIR.SR.jno2.nc',                            'jNO2',        0.0),
    'jHONO':     ('SAPHIR.SR.jhono.nc',                           'jHONO',       0.0),
    'jH2O2':     ('SAPHIR.SR.jh2o2.nc',                           'jH2O2',       0.0),
    'jHCHO_r':   ('SAPHIR.SR.jhcho_r.nc',                         'jHCHO_R',     0.0),
    'jHCHO_m':   ('SAPHIR.SR.jhcho_m.nc',                         'jHCHO_M',     0.0),
}

# Alternate variable name spellings to try if primary fails (case variants etc.)
ARCHIVE_MAP_ALTS = {
    'TOTAL_FLOW': ['total_flow', 'TOTALFLOW', 'Flow', 'flow'],
    'TK_M_SAPHIR': ['TK_SAPHIR', 'temp', 'TEMP', 'temperature', 'TK'],
    'H2O': ['h2o', 'water', 'H2O_mol'],
    'NO2': ['no2', 'No2'],
    'NO':  ['no', 'No'],
    'O3':  ['o3', 'ozone'],
    'jO1D': ['jo1d', 'JO1D', 'j_o1d'],
    'jNO2': ['jno2', 'JNO2', 'j_no2'],
    'jHONO': ['jhono', 'JHONO'],
    'jH2O2': ['jh2o2', 'JH2O2'],
    'jHCHO_R': ['jHCHO_r', 'jhcho_r', 'JHCHO_R'],
    'jHCHO_M': ['jHCHO_m', 'jhcho_m', 'JHCHO_M'],
}

# Wiki action code → default ENZ injection column
WIKI_CODE_TO_COL = {
    'O': 'QO3',
    'M': 'QNO',
    'D': 'QNO2',
    'C': 'QCO2',
    'X': 'QVOC',
    'H': 'QH2O',
    'S': 'Q_seed',
}

# Standard ENZ column order (TIME is first, written separately)
ENZ_COLUMN_ORDER = [
    'sz', 'hv', 'dilute', 't', 'H2O',
    'QH2O', 'qVCP1', 'qVCP2', 'qDIESEL', 'qGASOLINE',
    'qCOOKING', 'qBIOGENIC', 'q_isobutane',
    'QO3', 'QCO2', 'QCl2', 'QNO', 'Q_seed',
    'NO2', 'NO', 'O3',
    'jO1D', 'jNO2', 'jHONO', 'jH2O2', 'jHCHO_r', 'jHCHO_m',
]

# Default injection / manual-protocol columns shown in the injections table
DEFAULT_MANUAL_COLS = [
    ('QH2O',        '00:00,0'),
    ('QO3',         '00:00,0'),
    ('QNO',         '00:00,0'),
    ('QCO2',        '00:00,0'),
    ('QCl2',        '00:00,0'),
    ('Q_seed',      '00:00,0'),
    ('qVCP1',       '00:00,0'),
    ('qVCP2',       '00:00,0'),
    ('qDIESEL',     '00:00,0'),
    ('qGASOLINE',   '00:00,0'),
    ('qCOOKING',    '00:00,0'),
    ('qBIOGENIC',   '00:00,0'),
    ('q_isobutane', '00:00,0'),
]


# ══════════════════════════════════════════════════════════════════════════════
# Pure helper functions
# ══════════════════════════════════════════════════════════════════════════════

def compute_solar_zenith(grid_seconds):
    """
    Compute solar zenith angle (degrees) for SAPHIR (Jülich) at each time step.
    Uses Spencer (1971) / Iqbal (1983) algorithm. No external dependencies.

    Parameters
    ----------
    grid_seconds : array-like
        Seconds since ENZ_EPOCH (2000-01-01 00:00:00 UTC).

    Returns
    -------
    np.ndarray of SZA in degrees (0 = overhead, 90 = horizon).
    """
    lat_rad = math.radians(SAPHIR_LAT)
    results = []
    for s in grid_seconds:
        dt = ENZ_EPOCH + timedelta(seconds=float(s))
        doy = dt.timetuple().tm_yday
        B = 2.0 * math.pi * (doy - 1) / 365.0

        # Solar declination (radians) — Spencer (1971)
        dec = (0.006918
               - 0.399912 * math.cos(B)    + 0.070257 * math.sin(B)
               - 0.006758 * math.cos(2*B)  + 0.000907 * math.sin(2*B)
               - 0.002697 * math.cos(3*B)  + 0.001480 * math.sin(3*B))

        # Equation of time (minutes)
        eot = 229.18 * (0.000075
                        + 0.001868 * math.cos(B)  - 0.032077 * math.sin(B)
                        - 0.014615 * math.cos(2*B) - 0.040890 * math.sin(2*B))

        # Solar time (minutes from midnight)
        utc_min = dt.hour * 60.0 + dt.minute + dt.second / 60.0
        solar_time = utc_min + eot + 4.0 * SAPHIR_LON

        # Hour angle (radians): 0 at solar noon, ±π at midnight
        hour_angle = math.radians((solar_time / 4.0) - 180.0)

        cos_sza = (math.sin(lat_rad) * math.sin(dec)
                   + math.cos(lat_rad) * math.cos(dec) * math.cos(hour_angle))
        cos_sza = max(-1.0, min(1.0, cos_sza))
        results.append(math.degrees(math.acos(cos_sza)))

    return np.array(results)


def fetch_wiki_events(date_str):
    """
    Fetch and parse the SAPHIR action log from the wiki for a given date.

    Parameters
    ----------
    date_str : str  e.g. '2024-07-09'

    Returns
    -------
    list of dicts: {'time': 'HH:MM', 'code': 'R', 'description': '...'}

    Raises
    ------
    ImportError  if `requests` is not installed.
    requests.RequestException  if the page is unreachable.
    """
    try:
        import requests
    except ImportError:
        raise ImportError(
            "'requests' is not installed.\nRun:  pip install requests"
        )

    url = f"https://saphir.fz-juelich.de/saphirwiki/{date_str}"
    resp = requests.get(url, timeout=15, verify=True)
    resp.raise_for_status()

    # Strip HTML tags and decode common entities
    text = re.sub(r'<[^>]+>', ' ', resp.text)
    text = re.sub(r'&nbsp;', ' ', text)
    text = re.sub(r'&amp;', '&', text)
    text = re.sub(r'&lt;', '<', text)
    text = re.sub(r'&gt;', '>', text)

    # Action log lines: HH:MM  <single-char code>  <rest>
    # The code is one of: H h F f Z z P p R r U u V v O D M c C X Y y S s J j $
    pattern = re.compile(
        r'^\s*(\d{2}:\d{2})\s+([A-Za-z$])\b(.*)', re.MULTILINE
    )
    events = []
    for m in pattern.finditer(text):
        time_str = m.group(1)
        code     = m.group(2)
        desc     = m.group(3).strip()
        # Clean up excess whitespace from stripped HTML
        desc = re.sub(r'\s{2,}', ' ', desc)
        events.append({'time': time_str, 'code': code, 'description': desc})

    return events


def _convert_nc_time_to_enz_seconds(time_array, units_str):
    """
    Convert a NetCDF time array (with its 'units' attribute) to seconds
    since ENZ_EPOCH (2000-01-01 00:00:00 UTC).
    """
    m = re.match(
        r'(\w+)\s+since\s+'
        r'(\d{4})[-/](\d{1,2})[-/](\d{1,2})'
        r'(?:\s+(\d{1,2}):(\d{2})(?::(\d{2}))?)?',
        units_str.strip(), re.IGNORECASE
    )
    if not m:
        # Assume already in ENZ seconds
        return np.asarray(time_array, dtype=float)

    unit = m.group(1).lower()
    yr, mo, dy = int(m.group(2)), int(m.group(3)), int(m.group(4))
    hh = int(m.group(5) or 0)
    mm = int(m.group(6) or 0)
    ss = int(m.group(7) or 0)
    ref_dt = datetime(yr, mo, dy, hh, mm, ss, tzinfo=timezone.utc)
    ref_offset_s = (ref_dt - ENZ_EPOCH).total_seconds()

    factor = {'seconds': 1.0, 'minutes': 60.0, 'hours': 3600.0, 'days': 86400.0}.get(unit, 1.0)
    return np.asarray(time_array, dtype=float) * factor + ref_offset_s


def read_nc_variable(nc_path, var_name):
    """
    Read a variable from a NetCDF file.

    Returns
    -------
    (times_s, values) : tuple of np.ndarray
        times_s  — seconds since ENZ_EPOCH
        values   — data array (NaN where masked)

    Raises
    ------
    ImportError   if netCDF4 is not installed.
    FileNotFoundError  if nc_path does not exist.
    KeyError      if var_name is not in the file (message lists available vars).
    """
    try:
        import netCDF4 as nc4
    except ImportError:
        raise ImportError(
            "'netCDF4' is not installed.\nRun:  pip install netCDF4"
        )

    if not os.path.isfile(nc_path):
        raise FileNotFoundError(f"NC file not found: {nc_path}")

    with nc4.Dataset(nc_path, 'r') as ds:
        # Locate the time variable
        time_var = None
        for candidate in ('time', 'TIME', 'Time', 'epoch', 'seconds'):
            if candidate in ds.variables:
                time_var = ds.variables[candidate]
                break
        if time_var is None:
            raise KeyError(
                f"No time variable found in {os.path.basename(nc_path)}. "
                f"Variables: {list(ds.variables.keys())}"
            )

        units_str = getattr(time_var, 'units', 'seconds since 2000-1-1 00:00:00 UTC')
        raw_time  = np.array(time_var[:], dtype=float)
        times_s   = _convert_nc_time_to_enz_seconds(raw_time, units_str)

        # Locate the data variable (try primary name, then alternatives)
        candidates = [var_name] + ARCHIVE_MAP_ALTS.get(var_name, [])
        found_var = None
        for cand in candidates:
            if cand in ds.variables:
                found_var = cand
                break

        if found_var is None:
            avail = [v for v in ds.variables
                     if v not in ('time', 'TIME', 'Time', 'epoch', 'seconds')]
            raise KeyError(
                f"Variable '{var_name}' not found in {os.path.basename(nc_path)}.\n"
                f"Available: {avail}"
            )

        raw_data = ds.variables[found_var][:]
        if hasattr(raw_data, 'filled'):
            data = raw_data.filled(np.nan).astype(float).flatten()
        else:
            data = np.array(raw_data, dtype=float).flatten()

        # Trim to matching lengths
        n = min(len(times_s), len(data))
        return times_s[:n], data[:n]


def interpolate_to_grid(src_times, src_vals, grid_s, fill_value=0.0):
    """
    Linearly interpolate src data onto the ENZ time grid.
    Values outside the source range are extrapolated with fill_value.
    NaN gaps in source are filled by nearest valid neighbour then fill_value.
    """
    # Replace NaN with fill_value for interpolation
    valid = np.isfinite(src_vals)
    if not np.any(valid):
        return np.full(len(grid_s), fill_value)

    src_times_v = src_times[valid]
    src_vals_v  = src_vals[valid]

    result = np.interp(grid_s, src_times_v, src_vals_v,
                       left=fill_value, right=fill_value)
    return result


def step_to_grid(time_val_str, grid_s, date_str):
    """
    Convert a "HH:MM,value;HH:MM,value;..." string to a step-function array.

    The value at each time step holds until the next entry (piecewise-constant).
    Before the first entry the value is 0.

    Parameters
    ----------
    time_val_str : str  e.g. "00:00,0;11:49,0.2;11:54,0"
    grid_s       : array of seconds since ENZ_EPOCH
    date_str     : str  e.g. "2024-07-09"
    """
    date_obj = datetime.strptime(date_str, '%Y-%m-%d').replace(tzinfo=timezone.utc)
    day_start_s = (date_obj - ENZ_EPOCH).total_seconds()

    pairs = []
    for token in time_val_str.split(';'):
        token = token.strip()
        if not token:
            continue
        parts = token.split(',')
        if len(parts) < 2:
            continue
        try:
            h, mn = map(int, parts[0].strip().split(':'))
            val = float(parts[1].strip())
            t_s = day_start_s + h * 3600 + mn * 60
            pairs.append((t_s, val))
        except (ValueError, TypeError):
            continue

    if not pairs:
        return np.zeros(len(grid_s))

    pairs.sort(key=lambda x: x[0])
    result = np.zeros(len(grid_s))
    for i, t_grid in enumerate(grid_s):
        val = 0.0
        for t_step, v_step in pairs:
            if t_grid >= t_step:
                val = v_step
        result[i] = val
    return result


def hv_intervals_to_grid(intervals, grid_s, date_str):
    """
    Convert a list of (open_hhmm, close_hhmm) strings to a 0/1 step array.

    Parameters
    ----------
    intervals : list of (str, str)  e.g. [('09:00', '14:30')]
    grid_s    : array of seconds since ENZ_EPOCH
    date_str  : str  e.g. "2024-07-09"
    """
    date_obj = datetime.strptime(date_str, '%Y-%m-%d').replace(tzinfo=timezone.utc)
    day_start_s = (date_obj - ENZ_EPOCH).total_seconds()

    def hhmm_to_s(hhmm):
        h, mn = map(int, hhmm.strip().split(':'))
        return day_start_s + h * 3600 + mn * 60

    hv = np.zeros(len(grid_s))
    for (open_t, close_t) in intervals:
        try:
            t_open  = hhmm_to_s(open_t)
            t_close = hhmm_to_s(close_t)
            mask = (grid_s >= t_open) & (grid_s < t_close)
            hv[mask] = 1.0
        except (ValueError, TypeError):
            continue
    return hv


def write_enz(output_path, time_grid_s, column_data):
    """
    Write a SAPHIR .enz file.

    Parameters
    ----------
    output_path   : str
    time_grid_s   : np.ndarray  (seconds since ENZ_EPOCH)
    column_data   : dict {col_name: np.ndarray}  ordered; TIME excluded
    """
    cols = list(column_data.keys())
    n_cols = len(cols) + 1  # +1 for TIME

    def fmt(v):
        if not math.isfinite(v):
            return '0.00000000'
        if v == 0.0:
            return '0'
        if 0 < abs(v) < 1e-2:
            return f'{v:.8E}'
        return f'{v:.8f}'

    with open(output_path, 'w') as f:
        f.write('COLUMN 1=TIME / SEC SINCE 1.1.2000 00:00:00\n')
        for i, col in enumerate(cols, 2):
            f.write(f'COLUMN {i}={col}\n')
        f.write(f'NUMBER OF COLUMNS={n_cols}\n')
        f.write('&&&&&&&&&&&&&&&&&&&&&\n')
        for idx, t in enumerate(time_grid_s):
            vals = [fmt(float(t))] + [fmt(float(column_data[c][idx])) for c in cols]
            f.write(' ' + ', '.join(vals) + '\n')


# ══════════════════════════════════════════════════════════════════════════════
# Background thread for wiki fetching
# ══════════════════════════════════════════════════════════════════════════════

class WikiFetchThread(QThread):
    finished = pyqtSignal(list)   # list of event dicts
    error    = pyqtSignal(str)    # error message

    def __init__(self, date_str):
        super().__init__()
        self.date_str = date_str

    def run(self):
        try:
            events = fetch_wiki_events(self.date_str)
            self.finished.emit(events)
        except Exception as exc:
            self.error.emit(str(exc))


# ══════════════════════════════════════════════════════════════════════════════
# Main dialog
# ══════════════════════════════════════════════════════════════════════════════

class EnzBuilderDialog(QDialog):
    """
    Sub-GUI for building a SAPHIR .enz environment file.
    Opens from the main Scenario Builder via the "Prepare Chamber Environment" button (Step 1).
    Emits project_root_selected(str) with the output folder path after
    successful ENZ generation so the main GUI can auto-load the project root.
    """

    # Emitted with the output folder path when ENZ is written successfully
    project_root_selected = pyqtSignal(str)

    def __init__(self, parent=None, scenario_path=""):
        super().__init__(parent)
        self.setWindowTitle("Step 1 — Prepare Chamber Environment (.enz)")
        self.resize(1100, 750)
        self.scenario_path = scenario_path
        self._wiki_thread = None
        self._build_ui()

    # ── UI construction ────────────────────────────────────────────────────────

    def _build_ui(self):
        root = QVBoxLayout(self)
        root.setSpacing(8)

        # ── 1. Data source & date ─────────────────────────────────────────────
        src_group = QGroupBox("Data Source")
        src_grid = QGridLayout()
        src_grid.setColumnStretch(3, 1)

        self.rb_archive = QRadioButton("Archive date:")
        self.rb_custom  = QRadioButton("Custom NC folder:")
        self.rb_manual  = QRadioButton("All manual (no NC files)")
        self.rb_archive.setChecked(True)
        self.rb_archive.toggled.connect(self._on_source_changed)
        self.rb_custom.toggled.connect(self._on_source_changed)
        self.rb_manual.toggled.connect(self._on_source_changed)

        self.date_edit = QLineEdit("2024-07-09")
        self.date_edit.setFixedWidth(110)
        self.date_edit.setPlaceholderText("YYYY-MM-DD")

        self.archive_root_edit = QLineEdit(SAPHIR_ARCHIVE_ROOT)
        arc_browse = QPushButton("Browse"); arc_browse.setFixedWidth(70)
        arc_browse.clicked.connect(lambda: self._browse_dir(self.archive_root_edit))

        self.custom_nc_edit = QLineEdit()
        self.custom_nc_edit.setPlaceholderText(
            "Folder containing  YYYY-MM-DD.SAPHIR.*.nc  files"
        )
        cust_browse = QPushButton("Browse"); cust_browse.setFixedWidth(70)
        cust_browse.clicked.connect(lambda: self._browse_dir(self.custom_nc_edit))

        src_grid.addWidget(self.rb_archive,       0, 0)
        src_grid.addWidget(self.date_edit,        0, 1)
        src_grid.addWidget(QLabel("Archive root:"), 0, 2)
        src_grid.addWidget(self.archive_root_edit, 0, 3)
        src_grid.addWidget(arc_browse,            0, 4)

        src_grid.addWidget(self.rb_custom,        1, 0)
        src_grid.addWidget(self.custom_nc_edit,   1, 1, 1, 3)
        src_grid.addWidget(cust_browse,           1, 4)

        src_grid.addWidget(self.rb_manual,        2, 0, 1, 5)
        src_group.setLayout(src_grid)
        root.addWidget(src_group)

        # ── 2. Time range ─────────────────────────────────────────────────────
        time_group = QGroupBox("Time Range (UTC)")
        time_hbox = QHBoxLayout()
        time_hbox.addWidget(QLabel("Start:"))
        self.start_edit = QLineEdit("04:00"); self.start_edit.setFixedWidth(60)
        time_hbox.addWidget(self.start_edit)
        time_hbox.addWidget(QLabel("  End:"))
        self.end_edit = QLineEdit("22:00"); self.end_edit.setFixedWidth(60)
        time_hbox.addWidget(self.end_edit)
        time_hbox.addWidget(QLabel("  Step (s):"))
        self.step_edit = QLineEdit("60"); self.step_edit.setFixedWidth(55)
        time_hbox.addWidget(self.step_edit)
        time_hbox.addStretch()
        time_group.setLayout(time_hbox)
        root.addWidget(time_group)

        # ── 3. Main content: archive table (left) + manual protocol (right) ───
        splitter = QSplitter(Qt.Orientation.Horizontal)

        # ── Left: archive variables ───────────────────────────────────────────
        arc_group = QGroupBox("Archive Variables (from NC files)")
        arc_vbox = QVBoxLayout()

        self.archive_table = QTableWidget(len(ARCHIVE_MAP), 4)
        self.archive_table.setHorizontalHeaderLabels(
            ["", "ENZ column", "NC file suffix", "Variable in NC"]
        )
        self.archive_table.setColumnWidth(0, 22)
        self.archive_table.setColumnWidth(1, 80)
        self.archive_table.horizontalHeader().setSectionResizeMode(
            2, QHeaderView.ResizeMode.Stretch
        )
        self.archive_table.setColumnWidth(3, 110)
        self.archive_table.setSelectionMode(QAbstractItemView.SelectionMode.NoSelection)
        self.archive_table.verticalHeader().setDefaultSectionSize(22)

        for row, (col, (suffix, var, _fill)) in enumerate(ARCHIVE_MAP.items()):
            dot = QTableWidgetItem("●")
            dot.setForeground(QColor("#888"))
            dot.setFlags(Qt.ItemFlag.ItemIsEnabled)
            self.archive_table.setItem(row, 0, dot)

            col_item = QTableWidgetItem(col)
            col_item.setFlags(Qt.ItemFlag.ItemIsEnabled)
            self.archive_table.setItem(row, 1, col_item)

            suf_item = QTableWidgetItem(suffix)
            suf_item.setFlags(Qt.ItemFlag.ItemIsEnabled)
            self.archive_table.setItem(row, 2, suf_item)

            # NC variable name: editable so the user can override
            self.archive_table.setItem(row, 3, QTableWidgetItem(var))

        probe_btn = QPushButton("Probe NC Files")
        probe_btn.setToolTip(
            "Check which NC files exist for the selected date/folder "
            "and verify variable names."
        )
        probe_btn.clicked.connect(self._probe_nc_files)

        arc_vbox.addWidget(self.archive_table)
        arc_vbox.addWidget(probe_btn)
        arc_group.setLayout(arc_vbox)
        splitter.addWidget(arc_group)

        # ── Right: manual protocol ────────────────────────────────────────────
        manual_widget = QWidget()
        manual_vbox = QVBoxLayout(manual_widget)
        manual_vbox.setContentsMargins(0, 0, 0, 0)
        manual_vbox.setSpacing(6)

        # Roof / HV section
        hv_group = QGroupBox("Roof / UV lamp (HV)")
        hv_vbox = QVBoxLayout()

        wiki_hbox = QHBoxLayout()
        self.wiki_fetch_btn = QPushButton("Fetch from Wiki")
        self.wiki_fetch_btn.setToolTip(
            "Read action log from saphirwiki/YYYY-MM-DD.\n"
            "Requires VPN connection to Forschungszentrum Jülich."
        )
        self.wiki_fetch_btn.clicked.connect(self._fetch_wiki)
        self.wiki_status_lbl = QLabel("Not fetched")
        self.wiki_status_lbl.setStyleSheet("color: #888; font-size: 10px;")
        wiki_hbox.addWidget(self.wiki_fetch_btn)
        wiki_hbox.addWidget(self.wiki_status_lbl)
        wiki_hbox.addStretch()
        hv_vbox.addLayout(wiki_hbox)

        self.hv_table = QTableWidget(0, 3)
        self.hv_table.setHorizontalHeaderLabels(
            ["Open (HH:MM)", "Close (HH:MM)", "Notes from wiki"]
        )
        self.hv_table.horizontalHeader().setSectionResizeMode(
            2, QHeaderView.ResizeMode.Stretch
        )
        self.hv_table.setColumnWidth(0, 100)
        self.hv_table.setColumnWidth(1, 100)
        self.hv_table.setFixedHeight(100)
        self._add_hv_row("09:00", "18:00", "")   # placeholder default

        hv_btn_hbox = QHBoxLayout()
        add_hv_btn = QPushButton("+ Add interval"); add_hv_btn.setFixedWidth(110)
        add_hv_btn.clicked.connect(lambda: self._add_hv_row("", "", ""))
        del_hv_btn = QPushButton("Remove"); del_hv_btn.setFixedWidth(80)
        del_hv_btn.clicked.connect(self._remove_selected_hv)
        hv_btn_hbox.addWidget(add_hv_btn)
        hv_btn_hbox.addWidget(del_hv_btn)
        hv_btn_hbox.addStretch()
        hv_vbox.addWidget(self.hv_table)
        hv_vbox.addLayout(hv_btn_hbox)
        hv_group.setLayout(hv_vbox)
        manual_vbox.addWidget(hv_group)

        # Injections / boundary conditions section
        inj_group = QGroupBox("Injections & Boundary Conditions")
        inj_vbox = QVBoxLayout()

        hint = QLabel(
            "Format: HH:MM,value;HH:MM,value   "
            "(step function — value holds until next entry)"
        )
        hint.setStyleSheet("color: #888; font-size: 10px;")
        inj_vbox.addWidget(hint)

        self.inj_table = QTableWidget(0, 4)
        self.inj_table.setHorizontalHeaderLabels(
            ["☑", "Column", "Time-series  (HH:MM,val; ...)", "Notes from wiki"]
        )
        self.inj_table.setColumnWidth(0, 28)
        self.inj_table.setColumnWidth(1, 115)
        self.inj_table.horizontalHeader().setSectionResizeMode(
            2, QHeaderView.ResizeMode.Stretch
        )
        self.inj_table.horizontalHeader().setSectionResizeMode(
            3, QHeaderView.ResizeMode.ResizeToContents
        )
        self.inj_table.verticalHeader().setDefaultSectionSize(24)

        for col, series in DEFAULT_MANUAL_COLS:
            self._add_inj_row(col, series, "")

        inj_btn_hbox = QHBoxLayout()
        add_inj_btn = QPushButton("+ Add row"); add_inj_btn.setFixedWidth(80)
        add_inj_btn.clicked.connect(lambda: self._add_inj_row("", "00:00,0", ""))
        del_inj_btn = QPushButton("Remove"); del_inj_btn.setFixedWidth(80)
        del_inj_btn.clicked.connect(self._remove_selected_inj)
        inj_btn_hbox.addWidget(add_inj_btn)
        inj_btn_hbox.addWidget(del_inj_btn)
        inj_btn_hbox.addStretch()

        inj_vbox.addWidget(self.inj_table)
        inj_vbox.addLayout(inj_btn_hbox)
        inj_group.setLayout(inj_vbox)
        manual_vbox.addWidget(inj_group)

        splitter.addWidget(manual_widget)
        splitter.setSizes([400, 700])
        root.addWidget(splitter)

        # ── 4. Project Root (output folder) ───────────────────────────────────
        out_group = QGroupBox("2. Project Root  —  kammer.enz will be saved here")
        out_group.setStyleSheet("QGroupBox { font-weight: bold; }")
        out_hbox = QHBoxLayout()
        out_hbox.addWidget(QLabel("Folder:"))
        self.output_edit = QLineEdit(self.scenario_path or "")
        self.output_edit.setPlaceholderText("Browse to your scenario folder …")
        out_hbox.addWidget(self.output_edit, stretch=1)
        out_browse = QPushButton("Browse"); out_browse.setFixedWidth(70)
        out_browse.clicked.connect(self._browse_output)
        out_hbox.addWidget(out_browse)

        gen_btn = QPushButton("  Generate ENZ  ")
        gen_btn.setFixedHeight(42)
        gen_btn.setStyleSheet(
            "font-weight: bold; background-color: #1e7e34; color: white; padding: 0 18px;"
        )
        gen_btn.clicked.connect(self._generate_enz)
        out_hbox.addWidget(gen_btn)
        out_group.setLayout(out_hbox)
        root.addWidget(out_group)

        self._on_source_changed()

    # ── Helper: table row adders ───────────────────────────────────────────────

    def _add_hv_row(self, open_t, close_t, notes):
        r = self.hv_table.rowCount()
        self.hv_table.insertRow(r)
        self.hv_table.setItem(r, 0, QTableWidgetItem(open_t))
        self.hv_table.setItem(r, 1, QTableWidgetItem(close_t))
        notes_item = QTableWidgetItem(notes)
        notes_item.setForeground(QColor("#aaa"))
        self.hv_table.setItem(r, 2, notes_item)

    def _remove_selected_hv(self):
        rows = sorted({i.row() for i in self.hv_table.selectedItems()}, reverse=True)
        for r in rows:
            self.hv_table.removeRow(r)

    def _add_inj_row(self, col_name, series, notes, checked=True):
        r = self.inj_table.rowCount()
        self.inj_table.insertRow(r)

        chk = QTableWidgetItem()
        chk.setFlags(Qt.ItemFlag.ItemIsUserCheckable | Qt.ItemFlag.ItemIsEnabled)
        chk.setCheckState(
            Qt.CheckState.Checked if checked else Qt.CheckState.Unchecked
        )
        self.inj_table.setItem(r, 0, chk)
        self.inj_table.setItem(r, 1, QTableWidgetItem(col_name))
        self.inj_table.setItem(r, 2, QTableWidgetItem(series))
        notes_item = QTableWidgetItem(notes)
        notes_item.setForeground(QColor("#aaa"))
        self.inj_table.setItem(r, 3, notes_item)

    def _remove_selected_inj(self):
        rows = sorted({i.row() for i in self.inj_table.selectedItems()}, reverse=True)
        for r in rows:
            self.inj_table.removeRow(r)

    # ── Helper: browse dialogs ─────────────────────────────────────────────────

    def _browse_dir(self, line_edit):
        d = QFileDialog.getExistingDirectory(self, "Select folder")
        if d:
            line_edit.setText(d)

    def _browse_output(self):
        folder = QFileDialog.getExistingDirectory(
            self, "Select Project Root folder", self.output_edit.text()
        )
        if folder:
            self.output_edit.setText(folder)

    # ── Source toggle ──────────────────────────────────────────────────────────

    def _on_source_changed(self):
        archive_on = self.rb_archive.isChecked()
        custom_on  = self.rb_custom.isChecked()

        self.date_edit.setEnabled(archive_on)
        self.archive_root_edit.setEnabled(archive_on)

        self.custom_nc_edit.setEnabled(custom_on)

        # Dim the archive table when all-manual
        self.archive_table.setEnabled(not self.rb_manual.isChecked())

    # ── NC folder resolution ───────────────────────────────────────────────────

    def _get_nc_folder(self):
        """Return the folder to look for NC files in, or None if all-manual."""
        if self.rb_manual.isChecked():
            return None
        if self.rb_custom.isChecked():
            return self.custom_nc_edit.text().strip() or None

        # Archive date mode
        date_str = self.date_edit.text().strip()
        try:
            dt = datetime.strptime(date_str, '%Y-%m-%d')
        except ValueError:
            return None
        root = self.archive_root_edit.text().strip()
        return os.path.join(root, dt.strftime('%Y'), dt.strftime('%m'), dt.strftime('%d'))

    def _nc_path(self, suffix):
        """Build full path to a NC file: <folder>/<YYYY-MM-DD>.<suffix>"""
        folder = self._get_nc_folder()
        if not folder:
            return None
        date_str = self.date_edit.text().strip()
        return os.path.join(folder, f"{date_str}.{suffix}")

    # ── Probe NC files ─────────────────────────────────────────────────────────

    def _probe_nc_files(self):
        """Check each archive variable's NC file and update status dots."""
        nc4_available = True
        try:
            import netCDF4  # noqa: F401
        except ImportError:
            nc4_available = False

        for row, (col, (suffix, _var, _fill)) in enumerate(ARCHIVE_MAP.items()):
            var_item = self.archive_table.item(row, 3)
            var_name = var_item.text().strip() if var_item else _var
            dot = self.archive_table.item(row, 0)
            if not dot:
                continue

            if self.rb_manual.isChecked():
                dot.setForeground(QColor("#888"))
                dot.setToolTip("All-manual mode — NC files not used")
                continue

            nc_path = self._nc_path(suffix)
            if not nc_path or not os.path.isfile(nc_path):
                dot.setForeground(QColor("#e05252"))  # red = file missing
                dot.setToolTip(f"File not found:\n{nc_path}")
                continue

            if not nc4_available:
                dot.setForeground(QColor("#f0a500"))  # orange = can't read
                dot.setToolTip("File found but netCDF4 is not installed")
                continue

            # Try to open and find variable
            try:
                import netCDF4 as nc4
                with nc4.Dataset(nc_path, 'r') as ds:
                    candidates = [var_name] + ARCHIVE_MAP_ALTS.get(var_name, [])
                    found = any(c in ds.variables for c in candidates)
                    avail = [v for v in ds.variables
                             if v not in ('time', 'TIME', 'Time', 'epoch')]
                if found:
                    dot.setForeground(QColor("#4caf50"))  # green = OK
                    dot.setToolTip(f"OK  —  {os.path.basename(nc_path)}")
                else:
                    dot.setForeground(QColor("#f0a500"))  # orange = wrong var
                    dot.setToolTip(
                        f"File found but '{var_name}' not present.\n"
                        f"Available variables: {avail}\n"
                        f"Double-click the 'Variable in NC' cell to override."
                    )
            except Exception as exc:
                dot.setForeground(QColor("#e05252"))
                dot.setToolTip(f"Error reading file:\n{exc}")

        QMessageBox.information(
            self, "Probe complete",
            "● Green  = file & variable found\n"
            "● Orange = file found but variable name mismatch (edit the NC variable column)\n"
            "● Red    = file not found"
        )

    # ── Wiki fetch ─────────────────────────────────────────────────────────────

    def _fetch_wiki(self):
        date_str = self.date_edit.text().strip()
        if not date_str:
            QMessageBox.warning(self, "No date", "Please enter an experiment date first.")
            return

        self.wiki_fetch_btn.setEnabled(False)
        self.wiki_status_lbl.setText("Fetching…")
        self.wiki_status_lbl.setStyleSheet("color: #f0a500;")

        self._wiki_thread = WikiFetchThread(date_str)
        self._wiki_thread.finished.connect(self._on_wiki_fetched)
        self._wiki_thread.error.connect(self._on_wiki_error)
        self._wiki_thread.start()

    def _on_wiki_fetched(self, events):
        self.wiki_fetch_btn.setEnabled(True)
        if not events:
            self.wiki_status_lbl.setText("No events found")
            self.wiki_status_lbl.setStyleSheet("color: #888;")
            return

        self.wiki_status_lbl.setText(f"Fetched {len(events)} events — review & confirm below")
        self.wiki_status_lbl.setStyleSheet("color: #4caf50;")
        self._populate_from_wiki(events)

    def _on_wiki_error(self, msg):
        self.wiki_fetch_btn.setEnabled(True)
        self.wiki_status_lbl.setText(f"Error: {msg[:80]}")
        self.wiki_status_lbl.setStyleSheet("color: #e05252;")
        QMessageBox.warning(
            self, "Wiki fetch failed",
            f"Could not fetch wiki page.\n\n{msg}\n\n"
            "Make sure you are connected to the Forschungszentrum Jülich VPN."
        )

    def _populate_from_wiki(self, events):
        """
        Parse fetched wiki events and pre-fill HV table and injection rows.
        All auto-filled values remain fully editable by the user.
        """
        hv_opens  = []
        hv_closes = []
        injections = {}   # col_name → list of (time_str, description)

        for ev in events:
            code = ev['code']
            t    = ev['time']
            desc = ev['description']

            if code == 'R':                      # Roof open
                hv_opens.append((t, desc))
            elif code == 'r':                    # Roof close
                hv_closes.append((t, desc))
            elif code in WIKI_CODE_TO_COL:
                col = WIKI_CODE_TO_COL[code]
                injections.setdefault(col, []).append((t, desc))

        # ── Update HV table ────────────────────────────────────────────────────
        if hv_opens or hv_closes:
            self.hv_table.setRowCount(0)
            # Pair opens and closes in order
            n_pairs = max(len(hv_opens), len(hv_closes))
            for i in range(n_pairs):
                open_t  = hv_opens[i][0]  if i < len(hv_opens)  else ""
                close_t = hv_closes[i][0] if i < len(hv_closes) else ""
                notes   = []
                if i < len(hv_opens):
                    notes.append(f"Open: {hv_opens[i][1][:60]}")
                if i < len(hv_closes):
                    notes.append(f"Close: {hv_closes[i][1][:60]}")
                self._add_hv_row(open_t, close_t, " | ".join(notes))

        # ── Update injection rows ──────────────────────────────────────────────
        for col, entries in injections.items():
            # Find existing row for this column (if any)
            existing_row = None
            for r in range(self.inj_table.rowCount()):
                item = self.inj_table.item(r, 1)
                if item and item.text() == col:
                    existing_row = r
                    break

            # Build notes string from wiki descriptions
            notes_str = " | ".join(f"{t}: {d[:50]}" for t, d in entries)

            # Build a time-series hint: just set the times, user fills values
            # Format: first entry time with placeholder value 0 (user edits)
            hint_series = ";".join(f"{t},?" for t, _ in entries)

            if existing_row is not None:
                # Update notes column; mark series as needing user input
                series_item = self.inj_table.item(existing_row, 2)
                if series_item and series_item.text().strip() in ('00:00,0', ''):
                    series_item.setText(hint_series)
                    series_item.setForeground(QColor("#f0a500"))  # orange = needs value
                notes_item = self.inj_table.item(existing_row, 3)
                if notes_item:
                    notes_item.setText(notes_str)
            else:
                # Add new row for this column
                self._add_inj_row(col, hint_series, notes_str)
                # Colour the series orange to indicate values need filling
                r = self.inj_table.rowCount() - 1
                ser_item = self.inj_table.item(r, 2)
                if ser_item:
                    ser_item.setForeground(QColor("#f0a500"))

    # ── ENZ generation ─────────────────────────────────────────────────────────

    def _generate_enz(self):
        date_str = self.date_edit.text().strip()
        try:
            datetime.strptime(date_str, '%Y-%m-%d')
        except ValueError:
            QMessageBox.warning(self, "Bad date", "Enter the date as YYYY-MM-DD.")
            return

        project_root = self.output_edit.text().strip()
        if not project_root:
            QMessageBox.warning(self, "No project root", "Specify the Project Root folder.")
            return
        output_path = os.path.join(project_root, "kammer.enz")

        # Build time grid
        try:
            date_obj = datetime.strptime(date_str, '%Y-%m-%d').replace(tzinfo=timezone.utc)
            day_s = (date_obj - ENZ_EPOCH).total_seconds()
            sh, sm = map(int, self.start_edit.text().strip().split(':'))
            eh, em = map(int, self.end_edit.text().strip().split(':'))
            step_s = int(self.step_edit.text().strip())
            t_start = day_s + sh * 3600 + sm * 60
            t_end   = day_s + eh * 3600 + em * 60
            grid_s  = np.arange(t_start, t_end + 1, step_s)
        except Exception as exc:
            QMessageBox.warning(self, "Time range error", str(exc))
            return

        column_data = {}
        errors = []

        # ── SZ (computed) ──────────────────────────────────────────────────────
        column_data['sz'] = compute_solar_zenith(grid_s)

        # ── HV (from table) ────────────────────────────────────────────────────
        intervals = []
        for r in range(self.hv_table.rowCount()):
            o_item = self.hv_table.item(r, 0)
            c_item = self.hv_table.item(r, 1)
            open_t  = o_item.text().strip() if o_item else ""
            close_t = c_item.text().strip() if c_item else ""
            if open_t and close_t:
                intervals.append((open_t, close_t))
        column_data['hv'] = hv_intervals_to_grid(intervals, grid_s, date_str)

        # ── Archive variables ──────────────────────────────────────────────────
        if not self.rb_manual.isChecked():
            for row, (col, (suffix, _var, fill)) in enumerate(ARCHIVE_MAP.items()):
                var_item = self.archive_table.item(row, 3)
                var_name = var_item.text().strip() if var_item else _var
                nc_path = self._nc_path(suffix)
                try:
                    times_s, vals = read_nc_variable(nc_path, var_name)
                    column_data[col] = interpolate_to_grid(times_s, vals, grid_s, fill)
                except FileNotFoundError:
                    errors.append(f"{col}: file not found ({nc_path})")
                    column_data[col] = np.full(len(grid_s), fill)
                except (KeyError, ImportError, Exception) as exc:
                    errors.append(f"{col}: {exc}")
                    column_data[col] = np.full(len(grid_s), fill)
        else:
            # All-manual: fill archive columns with their fill values
            for col, (_suffix, _var, fill) in ARCHIVE_MAP.items():
                column_data[col] = np.full(len(grid_s), fill)

        # ── Manual injection columns (from table) ──────────────────────────────
        for r in range(self.inj_table.rowCount()):
            chk_item = self.inj_table.item(r, 0)
            if chk_item and chk_item.checkState() != Qt.CheckState.Checked:
                continue
            col_item = self.inj_table.item(r, 1)
            ser_item = self.inj_table.item(r, 2)
            if not col_item:
                continue
            col    = col_item.text().strip()
            series = ser_item.text().strip() if ser_item else "00:00,0"
            if not col:
                continue
            # Skip if the series still contains placeholder '?'
            if '?' in series:
                errors.append(
                    f"{col}: time-series contains '?' placeholders — fill in values before generating."
                )
                column_data[col] = np.zeros(len(grid_s))
                continue
            column_data[col] = step_to_grid(series, grid_s, date_str)

        # ── Assemble in ENZ column order ───────────────────────────────────────
        # Collect any extra columns not in the standard order
        standard = [c for c in ENZ_COLUMN_ORDER if c in column_data]
        extra    = [c for c in column_data if c not in ENZ_COLUMN_ORDER]
        ordered  = {c: column_data[c] for c in standard + extra}

        # ── Write ──────────────────────────────────────────────────────────────
        try:
            os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
            write_enz(output_path, grid_s, ordered)
        except Exception as exc:
            QMessageBox.critical(self, "Write error", str(exc))
            return

        # ── Summary ────────────────────────────────────────────────────────────
        msg = (
            f"ENZ file written successfully.\n\n"
            f"  Path     : {output_path}\n"
            f"  Date     : {date_str}\n"
            f"  Time     : {self.start_edit.text()} – {self.end_edit.text()} UTC\n"
            f"  Timestep : {step_s} s\n"
            f"  Rows     : {len(grid_s)}\n"
            f"  Columns  : {len(ordered) + 1}  (TIME + {len(ordered)})\n"
        )
        if errors:
            msg += (
                f"\n⚠  {len(errors)} column(s) used fill values due to errors:\n"
                + "\n".join(f"    • {e}" for e in errors)
            )
            QMessageBox.warning(self, "ENZ written (with warnings)", msg)
        else:
            QMessageBox.information(self, "ENZ written", msg)

        # Notify main GUI so it can auto-load this folder as Project Root (Step 2)
        self.project_root_selected.emit(project_root)
