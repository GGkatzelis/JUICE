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
    QRadioButton, QFileDialog, QMessageBox, QFrame,
    QHeaderView, QAbstractItemView,
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
    'dilute':    ('SAPHIR.collected.all_param.nc',                 'TOTAL_FLOW',  0.0),
    't':         ('SAPHIR.collected.all_param.nc',                 'TK_M_SAPHIR', 284.0),
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

# Wiki action code → ENZ injection column
# Fast codes: single event, 60 s default duration
WIKI_FAST_CODE_TO_COL = {
    'O': 'O3',
    'M': 'NO',
    'D': 'NO2',
    'C': 'CO2',
    'c': 'CO',
    'X': 'VOC',
}
# Slow codes: paired start (uppercase) / stop (lowercase) events
WIKI_SLOW_START_CODE = {'H': 'H2O', 'S': 'seed'}
WIKI_SLOW_STOP_CODE  = {'h': 'H2O', 's': 'seed'}
# Combined map kept for backwards compat
WIKI_CODE_TO_COL = {**WIKI_FAST_CODE_TO_COL,
                    **WIKI_SLOW_START_CODE,
                    **WIKI_SLOW_STOP_CODE}

# Fast injections default to 60 s; slow ones (H2O, seed) use wiki start/stop
FAST_INJECTION_COLS    = set(WIKI_FAST_CODE_TO_COL.values()) | {'Cl2'}
SLOW_INJECTION_COLS    = {'H2O', 'seed'}
FAST_DEFAULT_DURATION  = 60  # seconds

# Standard ENZ column order (TIME is first, written separately)
ENZ_COLUMN_ORDER = [
    'sz', 'hv', 'dilute', 't', 'H2O',
    'Injection_H2O', 'Injection_VCP1', 'Injection_VCP2', 'Injection_DIESEL', 'Injection_GASOLINE',
    'Injection_COOKING', 'Injection_BIOGENIC', 'Injection_isobutane',
    'Injection_O3', 'Injection_CO2', 'Injection_Cl2', 'Injection_NO', 'Injection_seed',
    'NO2', 'NO', 'O3',
    'jO1D', 'jNO2', 'jHONO', 'jH2O2', 'jHCHO_r', 'jHCHO_m',
]

# Default injection / manual-protocol columns shown in the injections table
# (col_name, is_fast)  — all unchecked by default; activated when wiki-populated
DEFAULT_MANUAL_COLS = [
    ('H2O',   False),
    ('O3',    True),
    ('NO',    True),
    ('NO2',   True),
    ('Cl2',   True),
    ('seed', False),
]

# Known VOC compounds for wiki description parsing
# lower-case key → (enz_column_suffix, formula, molar_mass_g_mol, liquid_density_g_mL or None)
VOC_KNOWN = {
    'limonene':           ('Limonene',    'C10H16', 136.23, 0.841),
    'd-limonene':         ('Limonene',    'C10H16', 136.23, 0.841),
    'alpha-pinene':       ('aPinene',     'C10H16', 136.23, 0.858),
    'a-pinene':           ('aPinene',     'C10H16', 136.23, 0.858),
    '\u03b1-pinene':      ('aPinene',     'C10H16', 136.23, 0.858),
    'beta-pinene':        ('bPinene',     'C10H16', 136.23, 0.872),
    'b-pinene':           ('bPinene',     'C10H16', 136.23, 0.872),
    '\u03b2-pinene':      ('bPinene',     'C10H16', 136.23, 0.872),
    'nopinene':           ('bPinene',     'C10H16', 136.23, 0.872),
    'toluene':            ('Toluene',     'C7H8',    92.14, 0.867),
    'isoprene':           ('Isoprene',    'C5H8',    68.12, 0.681),
    'ethanol':            ('Ethanol',     'C2H6O',   46.07, 0.789),
    'propene':            ('Propene',     'C3H6',    42.08, None),
    'propylene':          ('Propene',     'C3H6',    42.08, None),
    'isobutene':          ('Isobutene',   'C4H8',    56.11, 0.594),
    'isobutylene':        ('Isobutene',   'C4H8',    56.11, 0.594),
    'isobutane':          ('Isobutane',   'C4H10',   58.12, 0.551),
    'hexane':             ('Hexane',      'C6H14',   86.18, 0.659),
    'heptane':            ('Heptane',     'C7H16',  100.20, 0.684),
    'octane':             ('Octane',      'C8H18',  114.23, 0.703),
    'benzene':            ('Benzene',     'C6H6',    78.11, 0.879),
    'o-xylene':           ('oXylene',     'C8H10',  106.17, 0.897),
    'm-xylene':           ('mXylene',     'C8H10',  106.17, 0.860),
    'p-xylene':           ('pXylene',     'C8H10',  106.17, 0.861),
    'xylene':             ('Xylene',      'C8H10',  106.17, 0.861),
    'acetone':            ('Acetone',     'C3H6O',   58.08, 0.791),
    'methanol':           ('Methanol',    'CH4O',    32.04, 0.791),
    'mvk':                ('MVK',         'C4H6O',   70.09, 0.864),
    'methyl vinyl ketone':('MVK',         'C4H6O',   70.09, 0.864),
    'macr':               ('MACR',        'C4H6O',   70.09, 0.847),
    'methacrolein':       ('MACR',        'C4H6O',   70.09, 0.847),
    'mek':                ('MEK',         'C4H8O',   72.11, 0.805),
    'butanone':           ('MEK',         'C4H8O',   72.11, 0.805),
    'acetaldehyde':       ('AcetAldehyde','C2H4O',   44.05, 0.788),
    'formaldehyde':       ('HCHO',        'CH2O',    30.03, None),
    'cyclohexene':        ('Cyclohexene', 'C6H10',   82.14, 0.811),
    'sabinene':           ('Sabinene',    'C10H16', 136.23, 0.844),
    'camphene':           ('Camphene',    'C10H16', 136.23, 0.842),
    'myrcene':            ('Myrcene',     'C10H16', 136.23, 0.794),
    'ocimene':            ('Ocimene',     'C10H16', 136.23, 0.800),
    'linalool':           ('Linalool',    'C10H18O',154.25, 0.860),
    '1,8-cineole':        ('Cineol',      'C10H18O',154.25, 0.924),
    'eucalyptol':         ('Cineol',      'C10H18O',154.25, 0.924),
    'cineol':             ('Cineol',      'C10H18O',154.25, 0.924),
    'car-3-ene':          ('Car3ene',     'C10H16', 136.23, 0.858),
    '3-carene':           ('Car3ene',     'C10H16', 136.23, 0.858),
    'farnesene':          ('Farnesene',   'C15H24', 204.35, 0.841),
}


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
    # verify=False: FZJ uses an internal CA certificate not in the default bundle
    import urllib3
    urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
    resp = requests.get(url, timeout=15, verify=False)
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


def pulse_to_grid(pulses, grid_s, date_str):
    """
    Convert a list of (start_hhmm, stop_hhmm, value) injection pulse tuples
    to a summed array on the ENZ time grid.

    Each pulse contributes `value` during [start, stop); 0 outside.
    Multiple pulses for the same column are summed (additive injections).

    Parameters
    ----------
    pulses   : list of (start_hhmm: str, stop_hhmm: str, value: float)
    grid_s   : np.ndarray of seconds since ENZ_EPOCH
    date_str : str  e.g. '2024-07-09'
    """
    date_obj  = datetime.strptime(date_str, '%Y-%m-%d').replace(tzinfo=timezone.utc)
    day_s     = (date_obj - ENZ_EPOCH).total_seconds()
    result    = np.zeros(len(grid_s))

    for start_t, stop_t, val in pulses:
        try:
            sh, sm = map(int, start_t.strip().split(':'))
            eh, em = map(int, stop_t.strip().split(':'))
            t0 = day_s + sh * 3600 + sm * 60
            t1 = day_s + eh * 3600 + em * 60
            if t1 <= t0:
                continue
            mask = (grid_s >= t0) & (grid_s < t1)
            result[mask] += float(val)
        except (ValueError, AttributeError, TypeError):
            continue

    return result


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
        self._last_pasted_wiki_text = ""   # persists across re-opens of the paste dialog
        self._inj_table_updating   = False # re-entrance guard for cellChanged
        self._h2o_measured         = False # set True when H2O archive file found
        self._build_ui()

    # ── UI construction ────────────────────────────────────────────────────────

    def _build_ui(self):
        root = QVBoxLayout(self)
        root.setSpacing(6)

        # ── 1. Compact Experiment Setup (Data Source + Time Range in one bar) ─
        setup_group = QGroupBox("Experiment Setup")
        setup_outer = QHBoxLayout()
        setup_outer.setSpacing(10)

        # Left: radio-button source selectors
        src_grid = QGridLayout()
        src_grid.setColumnStretch(3, 1)
        src_grid.setHorizontalSpacing(6)
        src_grid.setVerticalSpacing(3)

        self.rb_archive = QRadioButton("Archive date:")
        self.rb_custom  = QRadioButton("Custom NC folder:")
        self.rb_manual  = QRadioButton("All manual (no NC)")
        self.rb_archive.setChecked(True)
        self.rb_archive.toggled.connect(self._on_source_changed)
        self.rb_custom.toggled.connect(self._on_source_changed)
        self.rb_manual.toggled.connect(self._on_source_changed)

        self.date_edit = QLineEdit("2024-07-09")
        self.date_edit.setFixedWidth(100)
        self.date_edit.setPlaceholderText("YYYY-MM-DD")

        self.archive_root_edit = QLineEdit(SAPHIR_ARCHIVE_ROOT)
        arc_browse = QPushButton("Browse"); arc_browse.setFixedWidth(65)
        arc_browse.clicked.connect(lambda: self._browse_dir(self.archive_root_edit))

        self.custom_nc_edit = QLineEdit()
        self.custom_nc_edit.setPlaceholderText(
            "Folder containing  YYYY-MM-DD.SAPHIR.*.nc  files"
        )
        cust_browse = QPushButton("Browse"); cust_browse.setFixedWidth(65)
        cust_browse.clicked.connect(lambda: self._browse_dir(self.custom_nc_edit))

        src_grid.addWidget(self.rb_archive,         0, 0)
        src_grid.addWidget(self.date_edit,          0, 1)
        src_grid.addWidget(QLabel("Archive root:"), 0, 2)
        src_grid.addWidget(self.archive_root_edit,  0, 3)
        src_grid.addWidget(arc_browse,              0, 4)
        src_grid.addWidget(self.rb_custom,          1, 0)
        src_grid.addWidget(self.custom_nc_edit,     1, 1, 1, 3)
        src_grid.addWidget(cust_browse,             1, 4)
        src_grid.addWidget(self.rb_manual,          2, 0, 1, 5)
        setup_outer.addLayout(src_grid, stretch=3)

        # Vertical separator
        sep = QFrame()
        sep.setFrameShape(QFrame.Shape.VLine)
        sep.setStyleSheet("color: #555;")
        setup_outer.addWidget(sep)

        # Right: time range (compact vertical form)
        time_grid = QGridLayout()
        time_grid.setVerticalSpacing(3)
        time_grid.setHorizontalSpacing(5)
        time_grid.addWidget(QLabel("Start (UTC):"), 0, 0)
        self.start_edit = QLineEdit("04:00"); self.start_edit.setFixedWidth(58)
        time_grid.addWidget(self.start_edit, 0, 1)
        time_grid.addWidget(QLabel("End:"), 1, 0)
        self.end_edit = QLineEdit("22:00"); self.end_edit.setFixedWidth(58)
        time_grid.addWidget(self.end_edit, 1, 1)
        time_grid.addWidget(QLabel("Step (s):"), 2, 0)
        self.step_edit = QLineEdit("60"); self.step_edit.setFixedWidth(58)
        time_grid.addWidget(self.step_edit, 2, 1)
        setup_outer.addLayout(time_grid, stretch=0)

        setup_group.setLayout(setup_outer)
        root.addWidget(setup_group)

        # ── 2. Main content splitter (Archive left | HV+Injections right) ─────
        main_splitter = QSplitter(Qt.Orientation.Horizontal)

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
        main_splitter.addWidget(arc_group)

        # ── Right: vertical splitter — HV (top, small) / Injections (bottom) ─
        v_splitter = QSplitter(Qt.Orientation.Vertical)

        # ── Roof / HV section ─────────────────────────────────────────────────
        hv_group = QGroupBox("Roof / UV lamp (HV)")
        hv_vbox = QVBoxLayout()
        hv_vbox.setSpacing(4)

        wiki_hbox = QHBoxLayout()
        self.wiki_fetch_btn = QPushButton("Fetch from Wiki")
        self.wiki_fetch_btn.setToolTip(
            "Read action log from saphirwiki/YYYY-MM-DD.\n"
            "Requires VPN connection to Forschungszentrum Jülich."
        )
        self.wiki_fetch_btn.clicked.connect(self._fetch_wiki)
        paste_btn = QPushButton("Paste Wiki Text")
        paste_btn.setToolTip(
            "Manually paste the action log text from the wiki page.\n"
            "Use this if the automatic fetch fails."
        )
        paste_btn.clicked.connect(self._paste_wiki_text)
        self.wiki_status_lbl = QLabel("Not fetched")
        self.wiki_status_lbl.setStyleSheet("color: #888; font-size: 10px;")
        wiki_hbox.addWidget(self.wiki_fetch_btn)
        wiki_hbox.addWidget(paste_btn)
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
        self.hv_table.verticalHeader().setDefaultSectionSize(22)
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
        v_splitter.addWidget(hv_group)

        # ── Injections / boundary conditions ─────────────────────────────────
        inj_group = QGroupBox("Injections & Boundary Conditions")
        inj_vbox = QVBoxLayout()
        inj_vbox.setSpacing(4)

        hint = QLabel(
            "All rows inactive until wiki-loaded or manually checked.  "
            "Fast injections (O\u2083, NO \u2026): 60 s default.  "
            "Slow (H\u2082O, seed): start/stop from wiki.  "
            "VOC rows added with compound name and estimated ppb when parseable."
        )
        hint.setStyleSheet("color: #888; font-size: 10px;")
        hint.setWordWrap(True)
        inj_vbox.addWidget(hint)

        # Columns: ☑(0) Name(1) Start(2) Stop(3) Dur(4) Vol.μl(5) Formula(6)
        #          Density(7) MW(8) Value_ppb(9) Notes(10)
        self.inj_table = QTableWidget(0, 11)
        self.inj_table.setHorizontalHeaderLabels([
            "☑", "Column", "Start\n(HH:MM)", "Stop\n(HH:MM)", "Dur.\n(s)",
            "Vol.\n(μl)", "Formula", "Density\n(g/mL)", "MW\n(g/mol)",
            "Value\n(ppb)", "Notes from wiki",
        ])
        self.inj_table.setColumnWidth(0, 28)
        self.inj_table.setColumnWidth(1, 130)
        self.inj_table.setColumnWidth(2, 80)
        self.inj_table.setColumnWidth(3, 80)
        self.inj_table.setColumnWidth(4, 60)
        self.inj_table.setColumnWidth(5, 58)   # Vol. (μl)
        self.inj_table.setColumnWidth(6, 68)   # Formula
        self.inj_table.setColumnWidth(7, 68)   # Density
        self.inj_table.setColumnWidth(8, 60)   # MW
        self.inj_table.setColumnWidth(9, 68)   # Value (ppb)
        self.inj_table.horizontalHeader().setSectionResizeMode(
            10, QHeaderView.ResizeMode.Stretch
        )
        self.inj_table.verticalHeader().setDefaultSectionSize(24)
        self.inj_table.cellChanged.connect(self._on_inj_cell_changed)

        for col_name, is_fast in DEFAULT_MANUAL_COLS:
            dur = str(FAST_DEFAULT_DURATION) if is_fast else ""
            self._add_inj_row(col_name, "", "", dur)  # all unchecked, no vol/formula/density/mw

        inj_btn_hbox = QHBoxLayout()
        add_inj_btn = QPushButton("+ Add row"); add_inj_btn.setFixedWidth(80)
        add_inj_btn.clicked.connect(
            lambda: self._add_inj_row("", "", "", str(FAST_DEFAULT_DURATION))
        )
        del_inj_btn = QPushButton("Remove"); del_inj_btn.setFixedWidth(80)
        del_inj_btn.clicked.connect(self._remove_selected_inj)
        inj_btn_hbox.addWidget(add_inj_btn)
        inj_btn_hbox.addWidget(del_inj_btn)
        inj_btn_hbox.addStretch()

        inj_vbox.addWidget(self.inj_table)
        inj_vbox.addLayout(inj_btn_hbox)
        inj_group.setLayout(inj_vbox)
        v_splitter.addWidget(inj_group)

        v_splitter.setSizes([140, 520])
        main_splitter.addWidget(v_splitter)
        main_splitter.setSizes([370, 730])
        root.addWidget(main_splitter, stretch=1)

        # ── 3. Project Root + Generate ────────────────────────────────────────
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
        gen_btn.setFixedHeight(38)
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

    # ── Static time helpers ────────────────────────────────────────────────────

    @staticmethod
    def _hhmm_to_s(hhmm):
        h, m = map(int, hhmm.strip().split(':'))
        return h * 3600 + m * 60

    @staticmethod
    def _s_to_hhmm(secs):
        secs = max(0, int(secs))
        return f"{secs // 3600:02d}:{(secs % 3600) // 60:02d}"

    # ── Injection table row management ────────────────────────────────────────

    def _add_inj_row(self, col_name, start_t, stop_t, duration_s,
                     vol_ul="", formula="", density="", mw="",
                     value="", notes="", checked=False):
        """
        Add one row to the injection table (11 columns).
        Cols: ☑(0) Name(1) Start(2) Stop(3) Dur(4) Vol.μl(5) Formula(6)
              Density(7) MW(8) Value_ppb(9) Notes(10)
        Duration is auto-calculated from Start+Stop if missing;
        Stop is auto-calculated from Start+Duration if Stop is missing.
        """
        self._inj_table_updating = True
        try:
            r = self.inj_table.rowCount()
            self.inj_table.insertRow(r)

            chk = QTableWidgetItem()
            chk.setFlags(Qt.ItemFlag.ItemIsUserCheckable | Qt.ItemFlag.ItemIsEnabled)
            chk.setCheckState(Qt.CheckState.Checked if checked else Qt.CheckState.Unchecked)
            self.inj_table.setItem(r, 0, chk)
            self.inj_table.setItem(r, 1, QTableWidgetItem(col_name))
            self.inj_table.setItem(r, 2, QTableWidgetItem(start_t))
            self.inj_table.setItem(r, 3, QTableWidgetItem(stop_t))

            # Auto-derive duration or stop if partial info given
            dur_str = str(duration_s)
            if start_t and stop_t and not duration_s:
                try:
                    dur_str = str(self._hhmm_to_s(stop_t) - self._hhmm_to_s(start_t))
                except ValueError:
                    pass
            elif start_t and duration_s and not stop_t:
                try:
                    new_stop = self._s_to_hhmm(self._hhmm_to_s(start_t) + int(duration_s))
                    self.inj_table.setItem(r, 3, QTableWidgetItem(new_stop))
                except ValueError:
                    pass

            self.inj_table.setItem(r, 4, QTableWidgetItem(dur_str))
            self.inj_table.setItem(r, 5, QTableWidgetItem(str(vol_ul)))
            self.inj_table.setItem(r, 6, QTableWidgetItem(str(formula)))
            self.inj_table.setItem(r, 7, QTableWidgetItem(str(density)))
            self.inj_table.setItem(r, 8, QTableWidgetItem(str(mw)))
            self.inj_table.setItem(r, 9, QTableWidgetItem(str(value)))

            notes_item = QTableWidgetItem(str(notes))
            notes_item.setForeground(QColor("#aaa"))
            self.inj_table.setItem(r, 10, notes_item)
        finally:
            self._inj_table_updating = False

    def _on_inj_cell_changed(self, row, col):
        """
        Auto-calculations triggered by cell edits:
        · Start or Stop changed  → recompute Duration
        · Duration changed       → recompute Stop
        · Vol (μl), Density, or MW changed → recompute Value (ppb)
        """
        if self._inj_table_updating:
            return
        self._inj_table_updating = True
        try:
            def txt(c):
                it = self.inj_table.item(row, c)
                return it.text().strip() if it else ""

            if col in (2, 3):
                start_t, stop_t = txt(2), txt(3)
                if start_t and stop_t:
                    try:
                        dur = self._hhmm_to_s(stop_t) - self._hhmm_to_s(start_t)
                        it  = self.inj_table.item(row, 4)
                        if it:
                            it.setText(str(dur))
                    except ValueError:
                        pass
            elif col == 4:
                start_t, dur_s = txt(2), txt(4)
                if start_t and dur_s:
                    try:
                        new_stop = self._s_to_hhmm(self._hhmm_to_s(start_t) + int(dur_s))
                        it = self.inj_table.item(row, 3)
                        if it:
                            it.setText(new_stop)
                    except ValueError:
                        pass
            elif col in (5, 7, 8):
                # Auto-calculate ppb from Vol (μl) + Density (g/mL) + MW (g/mol)
                vol_s, dens_s, mw_s = txt(5), txt(7), txt(8)
                if vol_s and dens_s and mw_s:
                    try:
                        vol_mL = float(vol_s) * 1e-3       # μl → mL
                        dens   = float(dens_s)
                        mw     = float(mw_s)
                        # n_air at standard conditions (T=298 K, P=1013 hPa, V=276 m³)
                        n_air  = (101_300.0 * 276.0) / (8.314 * 298.0)
                        n_voc  = (vol_mL * dens) / mw
                        ppb    = (n_voc / n_air) * 1e9
                        it = self.inj_table.item(row, 9)
                        if it:
                            it.setText(f'{ppb:.3f}')
                        else:
                            self.inj_table.setItem(row, 9, QTableWidgetItem(f'{ppb:.3f}'))
                    except (ValueError, ZeroDivisionError):
                        pass
        finally:
            self._inj_table_updating = False

    def _remove_selected_inj(self):
        rows = sorted({i.row() for i in self.inj_table.selectedItems()}, reverse=True)
        for r in rows:
            self.inj_table.removeRow(r)

    def _activate_inj_row(self, row):
        """Check the checkbox in the given injection row."""
        chk = self.inj_table.item(row, 0)
        if chk:
            chk.setCheckState(Qt.CheckState.Checked)

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

        # ── H2O measurement flag ───────────────────────────────────────────────
        h2o_row_idx = next((i for i, col in enumerate(ARCHIVE_MAP) if col == 'H2O'), None)
        h2o_measured = False
        if h2o_row_idx is not None:
            dot = self.archive_table.item(h2o_row_idx, 0)
            h2o_measured = (dot is not None
                            and dot.foreground().color().name() == '#4caf50')
        self._h2o_measured = h2o_measured

        for r in range(self.inj_table.rowCount()):
            col_it = self.inj_table.item(r, 1)
            if col_it and col_it.text() == 'H2O':
                note = ("⚠ H\u2082O measured — injection timescale for reference/plotting only"
                        if h2o_measured else "")
                col_it.setForeground(QColor("#aaa") if h2o_measured else QColor("white"))
                col_it.setToolTip(note)
                notes_it = self.inj_table.item(r, 10)
                if notes_it and h2o_measured and note not in notes_it.text():
                    notes_it.setText((notes_it.text() + "  " + note).strip())

        QMessageBox.information(
            self, "Probe complete",
            "● Green  = file & variable found\n"
            "● Orange = file found but variable name mismatch (edit the NC variable column)\n"
            "● Red    = file not found"
            + ("\n\nℹ H\u2082O measurement found — QH2O injection rows are reference/plotting only."
               if h2o_measured else "")
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
            "Make sure you are connected to the Forschungszentrum Jülich VPN.\n\n"
            "Alternatively, open the wiki page in a browser, copy the action log "
            "text, and use 'Paste Wiki Text'."
        )

    def _paste_wiki_text(self):
        """Open a text input where the user can paste the wiki action log manually."""
        from PyQt6.QtWidgets import QDialog, QDialogButtonBox, QPlainTextEdit

        dlg = QDialog(self)
        dlg.setWindowTitle("Paste Wiki Action Log")
        dlg.resize(600, 400)
        vbox = QVBoxLayout(dlg)
        vbox.addWidget(QLabel(
            "Open the wiki page in your browser, copy the action log text,\n"
            "and paste it below. Lines not matching  HH:MM <code>  are ignored."
        ))
        text_edit = QPlainTextEdit()
        text_edit.setPlaceholderText(
            "06:45 $ F.R. STATUS: Roof closed\n"
            "09:00 R F.R. CHANGE: open roof\n"
            "14:30 r F.R. CHANGE: close roof\n"
            "08:30 M M.S. CHANGE: NO injection ..."
        )
        if self._last_pasted_wiki_text:
            text_edit.setPlainText(self._last_pasted_wiki_text)
        vbox.addWidget(text_edit)
        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        buttons.accepted.connect(dlg.accept)
        buttons.rejected.connect(dlg.reject)
        vbox.addWidget(buttons)

        if dlg.exec() != QDialog.DialogCode.Accepted:
            return

        raw_text = text_edit.toPlainText()
        if not raw_text.strip():
            return

        self._last_pasted_wiki_text = raw_text   # remember for next re-open

        # Re-use the same regex parser as fetch_wiki_events
        pattern = re.compile(r'^\s*(\d{2}:\d{2})\s+([A-Za-z$])\b(.*)', re.MULTILINE)
        events = []
        for m in pattern.finditer(raw_text):
            desc = re.sub(r'\s{2,}', ' ', m.group(3).strip())
            events.append({'time': m.group(1), 'code': m.group(2), 'description': desc})

        if not events:
            QMessageBox.warning(self, "No events found",
                                "No lines matching  HH:MM <code>  were found in the pasted text.")
            return

        self.wiki_status_lbl.setText(f"Parsed {len(events)} events from pasted text — review below")
        self.wiki_status_lbl.setStyleSheet("color: #4caf50;")
        self._populate_from_wiki(events)

    # ── Injection table helpers ────────────────────────────────────────────────

    def _find_inj_row(self, col_name):
        """Return the first row index whose Column cell matches col_name, or None."""
        for r in range(self.inj_table.rowCount()):
            it = self.inj_table.item(r, 1)
            if it and it.text() == col_name:
                return r
        return None

    def _inj_row_is_empty(self, row):
        """True if Start, Stop, and Value (ppb) cells are all blank."""
        for c in (2, 3, 9):
            it = self.inj_table.item(row, c)
            if it and it.text().strip():
                return False
        return True

    def _set_inj_row(self, row, start_t, stop_t, dur_s, value, notes):
        """Overwrite an existing injection row's data cells (timing + value + notes)."""
        self._inj_table_updating = True
        try:
            for col_idx, val in ((2, start_t), (3, stop_t), (4, dur_s), (9, value)):
                it = self.inj_table.item(row, col_idx)
                if it:
                    it.setText(val)
                else:
                    self.inj_table.setItem(row, col_idx, QTableWidgetItem(val))
            notes_it = self.inj_table.item(row, 10)
            if notes_it:
                notes_it.setText(notes)
                notes_it.setForeground(QColor("#aaa"))
        finally:
            self._inj_table_updating = False

    def _mark_value_needed(self, row):
        """Colour the Value (ppb) cell orange to signal the user must fill in the amount."""
        it = self.inj_table.item(row, 9)
        if it:
            it.setForeground(QColor("#f0a500"))
            it.setToolTip("Enter the injection amount in ppb")

    def _parse_voc_description(self, desc):
        """
        Parse a VOC injection wiki description line.

        Returns
        -------
        (col_name, formula, density_str, mw_str, vol_ul_str, value_ppb_str, notes_str)
        All strings.  vol_ul_str is the raw extracted volume in μl (or '').
        value_ppb_str is computed ppb (or '' when conversion not possible).
        notes_str is the original wiki text + a brief conversion note.
        """
        desc_lower = desc.lower()

        # --- Identify compound ---
        col_suffix  = 'VOC'
        formula     = ''
        molar_mass  = None
        density     = None
        for key, (suf, form, mw, dens) in VOC_KNOWN.items():
            if key in desc_lower:
                col_suffix = suf
                formula    = form
                molar_mass = mw
                density    = dens
                break

        col_name    = col_suffix
        density_str = f'{density:.4f}' if density is not None else ''
        mw_str      = f'{molar_mass:.2f}' if molar_mass is not None else ''

        # --- Identify injection amount and unit ---
        unit_patterns = [
            (r'(\d+\.?\d*)\s*(\u03bcl|µl|ul)\b',  '\u03bcl'),
            (r'(\d+\.?\d*)\s*(ml|milliliter)\b',   'ml'),
            (r'(\d+\.?\d*)\s*(nl|nanoliter)\b',    'nl'),
            (r'(\d+\.?\d*)\s*ppm\b',               'ppm'),
            (r'(\d+\.?\d*)\s*ppb\b',               'ppb'),
        ]

        vol_ul_str    = ''
        value_ppb_str = ''
        note_suffix   = ''

        for pattern, unit in unit_patterns:
            m = re.search(pattern, desc_lower)
            if not m:
                continue
            val = float(m.group(1))
            if unit == 'ppb':
                value_ppb_str = f'{val:.3f}'
                note_suffix   = f'{val} ppb (direct)'
            elif unit == 'ppm':
                value_ppb_str = f'{val * 1000:.3f}'
                note_suffix   = f'{val} ppm \u2192 {val * 1000:.1f} ppb'
            elif unit in ('\u03bcl', 'ml', 'nl'):
                # Convert to μl for the Vol. column
                factor_to_ul = {'\u03bcl': 1.0, 'ml': 1000.0, 'nl': 1e-3}[unit]
                vol_ul        = val * factor_to_ul
                vol_ul_str    = f'{vol_ul:.4f}'
                vol_mL        = vol_ul * 1e-3
                if molar_mass and density:
                    n_air = (101_300.0 * 276.0) / (8.314 * 298.0)
                    n_voc = (vol_mL * density) / molar_mass
                    ppb   = (n_voc / n_air) * 1e9
                    value_ppb_str = f'{ppb:.3f}'
                    note_suffix   = (f'{val} {unit} \u2192 {ppb:.2f} ppb '
                                     f'(std 298K/1013hPa/276m\u00b3)')
                else:
                    note_suffix = (f'{val} {unit} \u2014 fill Density & MW '
                                   f'to auto-compute ppb')
            break

        notes_str = desc[:70] + (f'  |  {note_suffix}' if note_suffix else '')
        return col_name, formula, density_str, mw_str, vol_ul_str, value_ppb_str, notes_str

    def _populate_from_wiki(self, events):
        """
        Parse fetched/pasted wiki events and pre-fill HV table and injection rows.
        - Rows that receive data are checked (activated).
        - Fast injections (NO, O3 …): start = event time, 60 s default duration.
        - Slow injections (H2O, seed): start from uppercase code, stop from lowercase.
        - VOC (code X): compound name parsed, column named Injection_<Name>,
          ppb estimated from μl/ml volume if possible.
        All auto-filled values remain fully editable.
        """
        hv_opens  = []
        hv_closes = []

        fast_events        = {}   # col → [(start_t, desc)]
        slow_pairs         = {}   # col → [(start_t, stop_t, notes)]
        pending_slow_start = {}   # col → (start_t, desc)
        voc_events         = []   # [(start_t, col_name, formula, value_ppb_str, notes)]

        for ev in events:
            code = ev['code']
            t    = ev['time']
            desc = ev['description']

            if code == 'R':
                hv_opens.append((t, desc))
            elif code == 'r':
                hv_closes.append((t, desc))
            elif code == 'X':
                # VOC injection — parse compound name and convert value
                col_name, formula, density_str, mw_str, vol_ul_str, value_ppb_str, notes_str = \
                    self._parse_voc_description(desc)
                voc_events.append(
                    (t, col_name, formula, density_str, mw_str,
                     vol_ul_str, value_ppb_str, notes_str)
                )
            elif code in WIKI_FAST_CODE_TO_COL:
                col = WIKI_FAST_CODE_TO_COL[code]
                fast_events.setdefault(col, []).append((t, desc))
            elif code in WIKI_SLOW_START_CODE:
                col = WIKI_SLOW_START_CODE[code]
                pending_slow_start[col] = (t, desc)
            elif code in WIKI_SLOW_STOP_CODE:
                col = WIKI_SLOW_STOP_CODE[code]
                if col in pending_slow_start:
                    start_t, start_desc = pending_slow_start.pop(col)
                    notes = f"Start: {start_desc[:50]} | Stop: {desc[:50]}"
                    slow_pairs.setdefault(col, []).append((start_t, t, notes))
                else:
                    slow_pairs.setdefault(col, []).append(("", t, f"Stop: {desc[:60]}"))

        for col, (start_t, desc) in pending_slow_start.items():
            slow_pairs.setdefault(col, []).append(
                (start_t, "", f"Start: {desc[:60]} (stop not found)")
            )

        # ── Update HV table ────────────────────────────────────────────────────
        if hv_opens or hv_closes:
            self.hv_table.setRowCount(0)
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

        # ── Fast injections ────────────────────────────────────────────────────
        for col, entries in fast_events.items():
            for start_t, desc in entries:
                stop_t = self._s_to_hhmm(
                    self._hhmm_to_s(start_t) + FAST_DEFAULT_DURATION
                ) if start_t else ""
                notes_str = desc[:80]
                existing = self._find_inj_row(col)
                if existing is not None and self._inj_row_is_empty(existing):
                    self._set_inj_row(existing, start_t, stop_t,
                                      str(FAST_DEFAULT_DURATION), "", notes_str)
                    self._activate_inj_row(existing)
                    target_row = existing
                else:
                    self._add_inj_row(col, start_t, stop_t,
                                      str(FAST_DEFAULT_DURATION),
                                      notes=notes_str, checked=True)
                    target_row = self.inj_table.rowCount() - 1
                self._mark_value_needed(target_row)

        # ── Slow injections ────────────────────────────────────────────────────
        for col, pairs in slow_pairs.items():
            for start_t, stop_t, notes_str in pairs:
                dur_str = ""
                if start_t and stop_t:
                    try:
                        dur_str = str(self._hhmm_to_s(stop_t) - self._hhmm_to_s(start_t))
                    except ValueError:
                        pass
                existing = self._find_inj_row(col)
                if existing is not None and self._inj_row_is_empty(existing):
                    self._set_inj_row(existing, start_t, stop_t, dur_str, "", notes_str)
                    self._activate_inj_row(existing)
                    target_row = existing
                else:
                    self._add_inj_row(col, start_t, stop_t, dur_str,
                                      notes=notes_str, checked=True)
                    target_row = self.inj_table.rowCount() - 1
                self._mark_value_needed(target_row)

        # ── VOC injections ─────────────────────────────────────────────────────
        for start_t, col_name, formula, density_str, mw_str, vol_ul_str, value_ppb_str, notes_str \
                in voc_events:
            stop_t = self._s_to_hhmm(
                self._hhmm_to_s(start_t) + FAST_DEFAULT_DURATION
            ) if start_t else ""
            self._add_inj_row(col_name, start_t, stop_t,
                              str(FAST_DEFAULT_DURATION),
                              vol_ul=vol_ul_str, formula=formula,
                              density=density_str, mw=mw_str,
                              value=value_ppb_str, notes=notes_str,
                              checked=True)
            if not value_ppb_str:
                self._mark_value_needed(self.inj_table.rowCount() - 1)

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
        # Collect (start, stop, value) pulses grouped by column name
        pulses_by_col = {}
        for r in range(self.inj_table.rowCount()):
            chk_item   = self.inj_table.item(r, 0)
            if chk_item and chk_item.checkState() != Qt.CheckState.Checked:
                continue
            col_item   = self.inj_table.item(r, 1)
            start_item = self.inj_table.item(r, 2)
            stop_item  = self.inj_table.item(r, 3)
            val_item   = self.inj_table.item(r, 9)   # Value (ppb) column
            if not col_item:
                continue
            col     = col_item.text().strip()
            start_t = start_item.text().strip() if start_item else ""
            stop_t  = stop_item.text().strip()  if stop_item  else ""
            val_str = val_item.text().strip()   if val_item   else "0"
            if not col:
                continue
            
            # Prepend "Injection_" to the column name for the ENZ file
            full_col_name = f"Injection_{col}"

            if not start_t or not stop_t:
                if full_col_name not in column_data:
                    errors.append(f"{full_col_name}: missing Start or Stop — row skipped (set to 0).")
                    column_data.setdefault(full_col_name, np.zeros(len(grid_s)))
                continue
            try:
                if col in ('H2O', 'seed'):
                    val = 1.0
                else:
                    val = float(val_str) if val_str else 0.0
            except ValueError:
                errors.append(f"{full_col_name}: non-numeric value '{val_str}' — treated as 0.")
                val = 0.0
            pulses_by_col.setdefault(full_col_name, []).append((start_t, stop_t, val))

        for col, pulses in pulses_by_col.items():
            arr = pulse_to_grid(pulses, grid_s, date_str)
            # Add on top of archive data if present (e.g. QH2O reference on H2O)
            column_data[col] = column_data.get(col, np.zeros(len(grid_s))) + arr

        # ── Assemble in ENZ column order ───────────────────────────────────────
        # Collect any extra columns not in the standard order
        standard = [c for c in ENZ_COLUMN_ORDER if c in column_data]
        extra    = [c for c in column_data if c not in ENZ_COLUMN_ORDER]
        ordered  = {c: column_data[c] for c in standard + extra}

        # ── Preview — user confirms before writing ─────────────────────────────
        preview = _EnzPreviewDialog(
            grid_s, ordered, output_path, date_str, step_s, errors, self
        )
        if preview.exec() != QDialog.DialogCode.Accepted:
            return

        # ── Write ──────────────────────────────────────────────────────────────
        try:
            os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
            write_enz(output_path, grid_s, ordered)
        except Exception as exc:
            QMessageBox.critical(self, "Write error", str(exc))
            return

        n_cols = len(ordered)
        summary = (
            f"ENZ file written successfully.\n\n"
            f"  Path     : {output_path}\n"
            f"  Date     : {date_str}\n"
            f"  Time     : {self.start_edit.text()} – {self.end_edit.text()} UTC\n"
            f"  Timestep : {step_s} s\n"
            f"  Rows     : {len(grid_s)}\n"
            f"  Columns  : {n_cols + 1}  (TIME + {n_cols})\n"
        )
        if errors:
            summary += (
                f"\n⚠  {len(errors)} column(s) used fill values:\n"
                + "\n".join(f"    • {e}" for e in errors)
            )
            QMessageBox.warning(self, "ENZ written (with warnings)", summary)
        else:
            QMessageBox.information(self, "ENZ written", summary)

        # Notify main GUI so it can auto-load this folder as Project Root (Step 2)
        self.project_root_selected.emit(project_root)


# ══════════════════════════════════════════════════════════════════════════════
# ENZ pre-save preview dialog
# ══════════════════════════════════════════════════════════════════════════════

class _EnzPreviewDialog(QDialog):
    """
    4-panel matplotlib overview of the computed ENZ data shown before writing.
    User reviews the data then clicks "Confirm & Save ENZ" or "Cancel".
    """

    def __init__(self, grid_s, column_data, output_path,
                 date_str, step_s, errors, parent=None):
        super().__init__(parent)
        self.setWindowTitle("ENZ Preview — Review before saving")
        self.resize(1150, 680)
        self._grid_s      = grid_s
        self._col         = column_data
        self._output_path = output_path
        self._date_str    = date_str
        self._step_s      = step_s
        self._errors      = errors
        self._build_ui()

    def _build_ui(self):
        vbox = QVBoxLayout(self)
        vbox.setSpacing(6)

        # ── Info bar ──────────────────────────────────────────────────────────
        n_cols = len(self._col)
        info   = (f"  Output: {self._output_path}    │    "
                  f"Date: {self._date_str}    │    "
                  f"Step: {self._step_s} s    │    "
                  f"Rows: {len(self._grid_s)}    │    "
                  f"Columns: {n_cols + 1}  (TIME + {n_cols})")
        info_lbl = QLabel(info)
        info_lbl.setStyleSheet("font-size: 10px; color: #aaa;")
        vbox.addWidget(info_lbl)

        if self._errors:
            err_txt = (f"⚠  {len(self._errors)} column(s) used fill values:  "
                       + ";  ".join(self._errors[:4])
                       + ("  …" if len(self._errors) > 4 else ""))
            err_lbl = QLabel(err_txt)
            err_lbl.setStyleSheet("color: #f0a500; font-size: 10px;")
            err_lbl.setWordWrap(True)
            vbox.addWidget(err_lbl)

        # ── Matplotlib canvas ─────────────────────────────────────────────────
        try:
            from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
            from matplotlib.backends.backend_qtagg import NavigationToolbar2QT
            from matplotlib.figure import Figure

            fig    = Figure(tight_layout=True)
            canvas = FigureCanvasQTAgg(fig)
            tb     = NavigationToolbar2QT(canvas, self)
            self._draw_panels(fig)
            vbox.addWidget(tb)
            vbox.addWidget(canvas, stretch=1)
        except ImportError:
            vbox.addWidget(QLabel(
                "matplotlib not available — cannot show preview.\n"
                "Click 'Confirm & Save ENZ' to write the file directly."
            ))

        # ── Buttons ───────────────────────────────────────────────────────────
        btn_hbox = QHBoxLayout()
        btn_hbox.addStretch()
        cancel_btn = QPushButton("Cancel")
        cancel_btn.setFixedWidth(90)
        cancel_btn.clicked.connect(self.reject)
        save_btn = QPushButton("  Confirm & Save ENZ  ")
        save_btn.setFixedHeight(34)
        save_btn.setStyleSheet(
            "font-weight: bold; background-color: #1e7e34; color: white; padding: 0 16px;"
        )
        save_btn.clicked.connect(self.accept)
        btn_hbox.addWidget(cancel_btn)
        btn_hbox.addWidget(save_btn)
        vbox.addLayout(btn_hbox)

    # ── Plot ──────────────────────────────────────────────────────────────────

    def _draw_panels(self, fig):
        import numpy as np

        cd = self._col

        # Time axis: UTC hours from midnight
        date_obj  = datetime.strptime(self._date_str, '%Y-%m-%d').replace(tzinfo=timezone.utc)
        day_s     = (date_obj - ENZ_EPOCH).total_seconds()
        t_h       = (self._grid_s - day_s) / 3600.0

        axes = fig.subplots(4, 1, sharex=True)
        axes[-1].set_xlabel("Time UTC (h)", fontsize=9)

        def _hv_shade(ax):
            """Shade background yellow when roof is open (hv=1)."""
            if 'hv' in cd:
                ax.fill_between(
                    t_h, 0, 1,
                    where=cd['hv'] > 0.5,
                    transform=ax.get_xaxis_transform(),
                    alpha=0.10, color='yellow', zorder=0, label='HV on'
                )

        # ── Panel 0: Temperature ──────────────────────────────────────────────
        ax = axes[0]
        _hv_shade(ax)
        if 't' in cd:
            ax.plot(t_h, cd['t'], color='#e07b39', lw=1.3, label='T (K)')
            ax.set_ylabel('T (K)', fontsize=8, color='#e07b39')
            ax.tick_params(axis='y', colors='#e07b39', labelsize=7)
        if 'H2O' in cd:
            ax2 = ax.twinx()
            ax2.plot(t_h, cd['H2O'], color='steelblue', lw=0.9,
                     alpha=0.7, label='H₂O (ppb)')
            ax2.set_ylabel('H₂O (ppb)', fontsize=7, color='steelblue')
            ax2.tick_params(axis='y', colors='steelblue', labelsize=6)
        ax.set_title('Temperature & H₂O', fontsize=8)

        # ── Panel 1: NOx + O3 ────────────────────────────────────────────────
        ax = axes[1]
        _hv_shade(ax)
        _colors = {'NO': '#44aa44', 'NO2': '#dd8833', 'O3': '#9955cc'}
        for spec, col in _colors.items():
            if spec in cd:
                ax.plot(t_h, cd[spec], lw=1.2, color=col, label=spec)
        ax.set_ylabel('ppb', fontsize=8)
        ax.set_title('NO / NO₂ / O₃', fontsize=8)
        if any(s in cd for s in _colors):
            ax.legend(fontsize=7, loc='upper right', framealpha=0.5)

        # ── Panel 2: J-values ────────────────────────────────────────────────
        ax = axes[2]
        _hv_shade(ax)
        j_cols  = [c for c in cd if c.startswith('j')]
        ls_cycle = ['-', '--', '-.', ':', (0, (3, 1, 1, 1)), (0, (5, 2))]
        for i, jc in enumerate(j_cols):
            ax.plot(t_h, cd[jc], lw=0.9, ls=ls_cycle[i % len(ls_cycle)], label=jc)
        ax.set_ylabel('J (s⁻¹)', fontsize=8)
        ax.set_title('J-values', fontsize=8)
        if j_cols:
            ax.legend(fontsize=6, loc='upper right', ncol=3, framealpha=0.5)

        # ── Panel 3: Injections ───────────────────────────────────────────────
        ax = axes[3]
        _hv_shade(ax)
        inj_cols = [
            c for c in cd
            if c.startswith('Injection_')
            and np.any(cd[c] != 0)
        ]
        for i, ic in enumerate(inj_cols):
            ax.plot(t_h, cd[ic], lw=1.2, ls=ls_cycle[i % len(ls_cycle)], label=ic)
        ax.set_ylabel('ppb / a.u.', fontsize=8)
        ax.set_title('Injections & Boundary Conditions', fontsize=8)
        if inj_cols:
            ax.legend(fontsize=7, loc='upper right', framealpha=0.5)
        else:
            ax.text(0.5, 0.5, 'No active injections',
                    transform=ax.transAxes, ha='center', va='center',
                    color='#888', fontsize=9)

        for a in axes:
            a.tick_params(labelsize=7)
            a.grid(True, alpha=0.18, lw=0.5)
            a.set_xlim(t_h[0], t_h[-1])
