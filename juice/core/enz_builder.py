"""
enz_builder.py — the SAPHIR conditions as ENZ columns: archive files, solar zenith angle, the ENZ writer.

    compute_solar_zenith(grid_s, lat, lon)    solar zenith angle at the chamber's site
    read_nc_variable(path, name)              one variable of a SAPHIR archive NetCDF file, on ENZ seconds
    interpolate_to_grid(t, values, grid_s)    onto the ENZ time grid
    write_enz(path, grid_s, columns)          a .enz file in the IDL/FACSIMILE layout
    ARCHIVE_MAP, ARCHIVE_MAP_ALTS             which archive file and variable holds each ENZ column

ENZ times are seconds since 1 January 2000 00:00 UTC (ENZ_EPOCH).  Used by experiment_spec, field_scenario,
measured_inputs and the viewer.  The ENZ dialog of the early SAPHIR GUI that also lived here (EnzBuilderDialog,
with its wiki fetch and step, roof and pulse columns) was removed on 25 Nov 2026: nothing opened it any more; the
experiment window builds those columns in experiment_spec.
"""

import os
import re
import math
from datetime import datetime, timedelta, timezone

import numpy as np

# ── SAPHIR fixed constants ─────────────────────────────────────────────────────
SAPHIR_LAT = 50.9086       # degrees N, Jülich Research Centre
SAPHIR_LON = 6.3639        # degrees E

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


def compute_solar_zenith(grid_seconds, lat=SAPHIR_LAT, lon=SAPHIR_LON):
    """
    Compute solar zenith angle (degrees) at each time step for a site (default SAPHIR, Jülich).
    Uses Spencer (1971) / Iqbal (1983) algorithm. No external dependencies.

    Parameters
    ----------
    grid_seconds : array-like
        Seconds since ENZ_EPOCH (2000-01-01 00:00:00 UTC).
    lat, lon : float
        Site latitude (degrees N) and longitude (degrees E); from the chamber profile.

    Returns
    -------
    np.ndarray of SZA in degrees (0 = overhead, 90 = horizon).
    """
    lat_rad = math.radians(lat)
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
        solar_time = utc_min + eot + 4.0 * lon

        # Hour angle (radians): 0 at solar noon, ±π at midnight
        hour_angle = math.radians((solar_time / 4.0) - 180.0)

        cos_sza = (math.sin(lat_rad) * math.sin(dec)
                   + math.cos(lat_rad) * math.cos(dec) * math.cos(hour_angle))
        cos_sza = max(-1.0, min(1.0, cos_sza))
        results.append(math.degrees(math.acos(cos_sza)))

    return np.array(results)


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
