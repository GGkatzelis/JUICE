"""
measured_inputs.py — measured chamber conditions from a file the user provides (no Qt).

JUICE needs the physical conditions of the chamber over time: temperature, pressure, humidity,
dilution flow, photolysis frequencies, optionally NO / NO2 / O3 to hold them to measurements.
At SAPHIR they come from the archive; anywhere else the user can give one table:

    CSV (or Excel, first sheet): one time column + any of the standard columns below.
    Column names are matched case-insensitively, with common aliases; unknown columns are listed
    and ignored.  Time: ISO 'YYYY-MM-DD HH:MM[:SS]' (UTC), 'HH:MM[:SS]' of the experiment day,
    hours of day as a number, unix seconds, or seconds since 2000-01-01 (the ENZ epoch).

    standard name   aliases                       unit in the file        notes
    T               temp, temperature, TK, t_K    K  (°C if values < 150) chamber air temperature
    p               press, pressure, P_hPa        hPa (Pa if > 5000)
    H2O             h2o, water, H2O_percent       vol-%  (ppm if > 100)   or give RH:
    RH              rh, humidity                  %                       converted to H2O with T and p
    flow            dilute, total_flow            m3/h                    chamber dilution / replenishment flow
    roof            hv, lights, open              0 / 1                   roof open or lights on
    jNO2 jO1D jHONO jH2O2 jHCHO_r jHCHO_m          1/s                     photolysis frequencies
    NO NO2 O3       —                             ppb                     only needed to hold species to measurements

The template written by write_template() documents this in the file header.
"""
from __future__ import annotations

import datetime as _dt
import math
import os
import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np

ENZ_EPOCH_UNIX = 946684800.0

# ENZ column name → (aliases, unit text, kind)
COLUMNS: Dict[str, Tuple[List[str], str, str]] = {
    "t":       (["t", "temp", "temperature", "tk", "t_k", "tk_m_saphir", "temperature_k"], "K", "temperature"),
    "p":       (["p", "press", "pressure", "p_hpa", "pressure_hpa"], "hPa", "pressure"),
    "H2O":     (["h2o", "water", "h2o_percent", "h2o_vol", "water_vapour", "h2o_%"], "vol-%", "humidity"),
    "RH":      (["rh", "humidity", "relhum", "rel_hum", "rh_%"], "%", "humidity"),
    "dilute":  (["flow", "dilute", "total_flow", "dilution", "replenishment", "flow_m3h"], "m3/h", "flow"),
    "hv":      (["roof", "hv", "lights", "open", "roof_open", "light"], "0/1", "roof"),
    "jNO2":    (["jno2", "j_no2", "jno2_s-1"], "1/s", "photolysis"),
    "jO1D":    (["jo1d", "j_o1d", "jo(1d)"], "1/s", "photolysis"),
    "jHONO":   (["jhono", "j_hono"], "1/s", "photolysis"),
    "jH2O2":   (["jh2o2", "j_h2o2"], "1/s", "photolysis"),
    "jHCHO_r": (["jhcho_r", "jhcho_rad", "j_hcho_r"], "1/s", "photolysis"),
    "jHCHO_m": (["jhcho_m", "jhcho_mol", "j_hcho_m"], "1/s", "photolysis"),
    "NO":      (["no", "no_ppb"], "ppb", "species"),
    "NO2":     (["no2", "no2_ppb"], "ppb", "species"),
    "O3":      (["o3", "ozone", "o3_ppb"], "ppb", "species"),
}
TIME_ALIASES = ["time", "utc", "datetime", "date_time", "timestamp", "hour", "hours", "time_utc", "t_utc"]

TEMPLATE = """# JUICE measured inputs — the chamber conditions of one experiment, one row per time, UTC.
# Keep the header line (column names are matched ignoring case; see the aliases below), delete the
# columns you do not have, keep the rest.  Rows need not be evenly spaced: JUICE interpolates to its
# own time grid and holds the first / last value outside your rows.  Comma, semicolon or tab separated.
#
# REQUIRED
#   time      'YYYY-MM-DD HH:MM' (UTC), or 'HH:MM' of the experiment day, or hours of day as a number
#             (also: unix seconds, or seconds since 2000-01-01).  If no column is called time/utc/datetime,
#             the first column is taken.
# OPTIONAL — physical conditions (a missing column falls back to the constant typed in 'Where and when')
#   T         chamber air temperature, K            (aliases temp, temperature; °C is recognised when < 150)
#   p         pressure, hPa                         (aliases press, pressure; Pa is recognised when > 5000)
#   H2O       water vapour, vol-%                   (aliases water; ppm is recognised when > 100)
#   RH        relative humidity, %  — instead of H2O; converted with T and p (needs the T column)
#   flow      dilution / replenishment flow, m3/h   (aliases dilute, total_flow)
#   roof      1 = roof open / lights on, 0 = dark   (aliases hv, lights; values are rounded to 0/1)
#             — without this column the roof times typed in 'Where and when' are used
# OPTIONAL — light.  Give jNO2 whenever you have it: with jNO2 in the file every photolysis frequency the
#   file lacks is modelled from it (the mechanism's clear-sky spectrum scaled by your jNO2).  Without jNO2
#   the run uses the clear-sky model for the date, time and the Jülich sun.
#   jNO2 jO1D jHONO jH2O2 jHCHO_r jHCHO_m   photolysis frequencies, 1/s   (aliases j_no2, jo(1d), jhcho_rad, jhcho_mol …)
# OPTIONAL — measured trace gases, ppb: only used when a species is held to its measurement
#   NO NO2 O3                                                                 (aliases ozone, no_ppb …)
#
# Example rows (delete them):
time,T,p,H2O,flow,roof,jNO2,jO1D,jHONO,NO,NO2,O3
2024-07-09 06:00,295.6,1006.5,1.33,6.0,0,0,0,0,0.05,0.02,1.2
2024-07-09 06:01,295.6,1006.5,1.33,6.0,0,0,0,0,0.05,0.02,1.2
2024-07-09 12:20,303.1,1006.1,1.52,6.0,1,0.0071,2.1e-5,0.0013,0.30,1.9,35.0
"""


@dataclass
class InputsFile:
    path: str
    times_s: np.ndarray                                   # seconds since 2000-01-01 (ENZ epoch)
    cols: Dict[str, np.ndarray] = field(default_factory=dict)   # ENZ column name → values in ENZ units
    recognised: List[str] = field(default_factory=list)   # "T ← temperature (K)"
    ignored: List[str] = field(default_factory=list)
    notes: List[str] = field(default_factory=list)

    def summary(self) -> str:
        s = f"{len(self.times_s)} rows, {self.hours[0]:.2f}–{self.hours[-1]:.2f} h; " if len(self.times_s) else "no rows; "
        s += "recognised: " + (", ".join(self.cols) or "nothing")
        if self.ignored:
            s += "; ignored: " + ", ".join(self.ignored)
        return s

    @property
    def hours(self) -> np.ndarray:
        if not len(self.times_s):
            return self.times_s
        day0 = math.floor((self.times_s[0] + ENZ_EPOCH_UNIX) / 86400.0) * 86400.0 - ENZ_EPOCH_UNIX
        return (self.times_s - day0) / 3600.0


def write_template(path: str) -> str:
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write(TEMPLATE)
    return path


def _load_table(path: str):
    import pandas as pd
    ext = os.path.splitext(path)[1].lower()
    if ext in (".xlsx", ".xls", ".xlsm"):
        return pd.read_excel(path)
    with open(path, encoding="utf-8-sig", errors="replace") as f:
        head = [l for l in (f.readline() for _ in range(5)) if l]
    sep = ";" if sum(l.count(";") for l in head) > sum(l.count(",") for l in head) else ","
    if sep == "," and sum(l.count("\t") for l in head) > sum(l.count(",") for l in head):
        sep = "\t"
    return pd.read_csv(path, sep=sep, comment="#", encoding="utf-8-sig", skipinitialspace=True)


def _parse_time(series, date: str) -> np.ndarray:
    """Any of the accepted time forms → seconds since 2000-01-01."""
    import pandas as pd
    day = _dt.datetime.strptime(date, "%Y-%m-%d").replace(tzinfo=_dt.timezone.utc) if date else None
    if pd.api.types.is_numeric_dtype(series):
        v = series.to_numpy(dtype=float)
        med = float(np.nanmedian(v))
        if med > 1e9:                                   # unix seconds
            return v - ENZ_EPOCH_UNIX
        if med > 1e8:                                   # ENZ seconds since 2000
            return v
        if day is None:
            raise ValueError("hours of day need the experiment date")
        return day.timestamp() - ENZ_EPOCH_UNIX + v * 3600.0        # hours of the day
    txt = series.astype(str).str.strip()
    if txt.str.fullmatch(r"\d{1,2}:\d{2}(:\d{2})?").all():
        if day is None:
            raise ValueError("'HH:MM' times need the experiment date")
        parts = txt.str.split(":")
        secs = parts.map(lambda p: int(p[0]) * 3600 + int(p[1]) * 60 + (int(p[2]) if len(p) > 2 else 0)).to_numpy(dtype=float)
        return day.timestamp() - ENZ_EPOCH_UNIX + secs
    ts = pd.to_datetime(txt, utc=True, errors="coerce")
    if ts.isna().all():
        raise ValueError("time column not understood (use 'YYYY-MM-DD HH:MM', 'HH:MM', hours of day or seconds)")
    unix = ts.map(lambda x: x.timestamp() if x is not pd.NaT else float("nan")).to_numpy(dtype=float)
    return unix - ENZ_EPOCH_UNIX


def _es_hpa(T: np.ndarray) -> np.ndarray:
    """Saturation vapour pressure over water, hPa (Magnus)."""
    tc = T - 273.15
    return 6.112 * np.exp(17.62 * tc / (243.12 + tc))


def read_inputs_file(path: str, date: str = "", default_p_hpa: float = 1006.5) -> InputsFile:
    """Read the user's table; recognised columns come back in ENZ units (T K, p hPa, H2O vol-%, flow m3/h, J 1/s, ppb)."""
    df = _load_table(path)
    names = {str(c).strip().lower(): c for c in df.columns}
    tcol = next((names[a] for a in TIME_ALIASES if a in names), None)
    if tcol is None:
        tcol = df.columns[0]
    times = _parse_time(df[tcol], date)
    order = np.argsort(times); times = times[order]
    ok = np.isfinite(times)
    out = InputsFile(path=path, times_s=times[ok])
    used = {tcol}
    for enz, (aliases, unit, kind) in COLUMNS.items():
        src = next((names[a] for a in aliases if a in names), None)
        if src is None:
            continue
        used.add(src)
        v = np.asarray(df[src].to_numpy(dtype=float)[order][ok], dtype=float)
        med = float(np.nanmedian(v)) if np.isfinite(v).any() else float("nan")
        if enz == "t" and med < 150:
            v = v + 273.15; out.notes.append("T given in °C → K")
        if enz == "p" and med > 5000:
            v = v / 100.0; out.notes.append("p given in Pa → hPa")
        if enz == "H2O" and med > 100:
            v = v / 1e4; out.notes.append("H2O given in ppm → vol-%")
        out.cols[enz] = v
        out.recognised.append(f"{enz} ← {src} ({unit})")
    if "H2O" not in out.cols and "RH" in out.cols:
        T = out.cols.get("t"); p = out.cols.get("p", np.full_like(out.cols["RH"], default_p_hpa))
        if T is not None:
            out.cols["H2O"] = out.cols["RH"] / 100.0 * _es_hpa(T) / p * 100.0
            out.notes.append("H2O computed from RH, T and p")
        else:
            out.notes.append("RH given but no temperature column: humidity cannot be converted")
    out.cols.pop("RH", None)
    out.ignored = [str(c) for c in df.columns if c not in used]
    return out


def column_on_grid(inp: InputsFile, enz_name: str, grid_s: np.ndarray, fill: float) -> Optional[np.ndarray]:
    """The file's column interpolated to the ENZ grid (None when the file lacks it)."""
    if enz_name not in inp.cols or not len(inp.times_s):
        return None
    from enz_builder import interpolate_to_grid
    return np.asarray(interpolate_to_grid(inp.times_s, inp.cols[enz_name], grid_s, fill), dtype=float)


_CACHE: Dict[Tuple[str, float, str], InputsFile] = {}


def cached_inputs(path: str, date: str, default_p_hpa: float = 1006.5) -> Optional[InputsFile]:
    """read_inputs_file once per file version (path + mtime + date)."""
    if not path or not os.path.isfile(path):
        return None
    key = (os.path.abspath(path), os.path.getmtime(path), date)
    if key not in _CACHE:
        _CACHE[key] = read_inputs_file(path, date, default_p_hpa)
    return _CACHE[key]
