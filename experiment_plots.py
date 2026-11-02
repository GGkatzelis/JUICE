"""
experiment_plots.py — what to show for an experiment, computed once, drawn anywhere.

Pure functions (no Qt) shared by the step window (inputs plots) and the results viewer
(Overview tab):

  spec_for_run(cs_path)            the .experiment.json next to a result file, if any
  injected_species(spec)           MCM species that were injected (modelled products)
  chamber_sources(spec, cols)      HONO/HCHO wall sources, dark HONO [ppb/h], dilution [%/h]
  OBS_MAP / load_observations()    SAPHIR archive observations for a date, in model units;
                                   extracted once into <folder>/<date>.observations.csv (cache_dir)
  landing_panels(spec, columns)    the panels of the Overview tab for this run

SAPHIR-specific knowledge (instrument files, variable names, units) lives in OBS_MAP; another
chamber gets its own map.
"""
from __future__ import annotations

import os
import re
from typing import Dict, List, Optional, Tuple

import numpy as np

from experiment_spec import ExperimentSpec, load_spec
import injection_calc as ic

# ── observations: key → (file suffix, variable, unit of the file, factor to model unit, label, kind) ──
# kinds: 'ppb' (model: number density ÷ M × 1e9), 'cm-3' (model: number density), 's-1' (model: kOH from the .pd file)
OBS_MAP: Dict[str, Tuple[str, str, str, float, str, str]] = {
    "NO":   ("SAPHIR.TR43.NO_NO2_NOx.nc",                    "NO",      "ppb",  1.0,    "NO (TR43)",        "ppb"),
    "NO2":  ("SAPHIR.TR43.NO_NO2_NOx.nc",                    "NO2",     "ppb",  1.0,    "NO2 (TR43)",       "ppb"),
    "O3":   ("SAPHIR.TwoBTech_1014.O3.nc",                   "O3",      "ppb",  1.0,    "O3 (2B)",          "ppb"),
    "CO":   ("SAPHIR.PICARRO_CFKADS2270.CO_CO2_CH4_H2O.nc", "CO",      "ppm",  1000.0, "CO (Picarro)",     "ppb"),
    "HONO": ("SAPHIR.ICAD01062.HONO.nc",                     "HONO",    "ppb",  1.0,    "HONO (ICAD)",      "ppb"),
    "OH":   ("SAPHIR.LIF.ROx.nc",                            "OH_LIF",  "cm-3", 1.0,    "OH (LIF)",         "cm-3"),
    "HO2":  ("SAPHIR.LIF.ROx.nc",                            "HO2_LIF", "cm-3", 1.0,    "HO2 (LIF)",        "cm-3"),
    "RO2":  ("SAPHIR.LIF.ROx.nc",                            "RO2_LIF", "cm-3", 1.0,    "RO2 (LIF)",        "cm-3"),
    "kOH":  ("SAPHIR.LP_LIF.kOH.nc",                         "kOH",     "s-1",  1.0,    "kOH (LP-LIF)",     "s-1"),
    "jNO2": ("SAPHIR.SR.jNO2.nc",                            "jNO2",    "s-1",  1.0,    "jNO2 (SR)",        "s-1"),
}
# Older archive naming (e.g. 2017): tried after the OBS_MAP file; same variable unit as there.
OBS_OLD_NAMES: Dict[str, List[Tuple[str, str]]] = {
    "NO":  [("NOxContainer.no.nc", "NO")],
    "NO2": [("NOxContainer.no2.nc", "NO2")],
    "O3":  [("NOxContainer.o3_uv.nc", "O3")],
    "CO":  [("NOxContainer_PICARRO.CO.nc", "CO")],
    "OH":  [("LIF_SAPHIR.ROx.nc", "OH_LIF")],
    "HO2": [("LIF_SAPHIR.ROx.nc", "HO2_LIF")],
    "RO2": [("LIF_SAPHIR.ROx.nc", "RO2_LIF")],
    "kOH": [("LIF_SAPHIR.kOH.nc", "kOH")],
    "jNO2": [("saphir_photolysis.jno2.nc", "jNO2")],
}
ENZ_EPOCH_UNIX = 946684800.0            # 2000-01-01 00:00 UTC


def spec_for_run(cs_path: str) -> Optional[ExperimentSpec]:
    """The experiment description written next to a result file by the step window."""
    stem = os.path.splitext(os.path.basename(cs_path))[0]
    p = os.path.join(os.path.dirname(os.path.abspath(cs_path)), stem + ".experiment.json")
    if not os.path.isfile(p):
        return None
    try:
        return load_spec(p)
    except Exception:
        return None


def injected_species(spec: ExperimentSpec) -> List[str]:
    """MCM species that enter through injections (products and mixture components), in order."""
    out: List[str] = []
    for inj in spec.injections:
        for term in (inj.products or "").split("+"):
            name = term.split("*")[-1].strip().upper()
            if name and name not in ("XNOY", "HV") and name not in out and re.fullmatch(r"[A-Z][A-Z0-9]*", name):
                out.append(name)
    return out


def _const(expr: str, default: float = 0.0) -> float:
    try:
        return float(eval(expr.replace("@", "**").replace("D", "e").replace("d", "e"), {"__builtins__": {}}, {}))
    except Exception:
        return default


def chamber_sources(spec: ExperimentSpec, cols: Dict[str, np.ndarray], n: int) -> Dict[str, np.ndarray]:
    """Chamber contributions per hour from the ENZ columns and the chamber settings:
    'HONO wall source', 'HCHO wall source', 'dark HONO' [ppb/h]; 'dilution' [%/h]; plus 'jNO2' used."""
    c, ch = spec.conditions, spec.chamber
    T = np.asarray(cols.get("t", np.full(n, c.temperature_K)), dtype=float); P = c.pressure_hPa
    M = 7.2428e18 * P / np.maximum(T, 1.0)                     # molecules/cm³ = P/(k_B T)
    h2o = np.asarray(cols.get("H2O", np.full(n, c.h2o_percent)), dtype=float)   # vol-%
    hv = np.asarray(cols.get("hv", np.zeros(n)), dtype=float)
    j = cols.get("jNO2", cols.get("jno2"))
    j = np.asarray(j, dtype=float) if j is not None else np.zeros(n)
    if not np.any(j > 0) and ic.chamber(c.chamber_name).light_type == "sun":   # manual / clear-sky run: jNO2 = PHTFKT(1.165e-2, 0.244, -0.267); lamp/dark chambers carry their jNO2 column
        sza = np.asarray(cols.get("sz", np.full(n, 90.0)), dtype=float)
        cz = np.maximum(np.cos(np.radians(np.minimum(sza, 89.9))), 1e-3)
        j = np.where(sza < 89.9, 1.165e-2 * cz ** 0.244 * np.exp(-0.267 / cz), 0.0) * (hv > 0.5)
    to_ppb_h = 1e9 * 3600.0 / M
    rh = h2o * P / np.exp(21.36469 - 5339.66 / np.maximum(T, 1.0))
    out: Dict[str, np.ndarray] = {"jNO2": j}
    if ch.hono:
        out["HONO wall source"] = _const(ch.hono_const, 1.82) * 4.7e13 * (1 + (rh / 11.6) ** 2) * np.exp(-3950.0 / T) * j * to_ppb_h
    if ch.hcho:
        out["HCHO wall source"] = _const(ch.hcho_const, 1.0) * 2.8e6 * hv * to_ppb_h
    if ch.dark_hono:
        out["dark HONO"] = _const(ch.dark_hono_rate, 1e-13) * h2o * 1e-2 * M * to_ppb_h
    if ch.dilution:
        flow = np.asarray(cols.get("dilute", np.full(n, c.flow_m3h)), dtype=float)
        prof = ic.chamber(c.chamber_name)
        out["dilution"] = (flow + prof.leak_flow_m3h) / (prof.volume_m3 * 3600.0) * 3600.0 * 100.0
    return out


def archive_day(archive_root: str, date: str) -> str:
    try:
        y, m, d = date.split("-")
        p = os.path.join(archive_root, y, m, d)
        return p if os.path.isdir(p) else ""
    except Exception:
        return ""


def _midnight_of(date: str) -> float:
    import datetime as _dt
    y, m, d = (int(x) for x in date.split("-"))
    return _dt.datetime(y, m, d, tzinfo=_dt.timezone.utc).timestamp()


ObsDict = Dict[str, Tuple[np.ndarray, np.ndarray, str, str]]


def extract_observations(archive_root: str, date: str, keys: Optional[List[str]] = None) -> Tuple[ObsDict, List[str]]:
    """Read the mapped instrument files of that day straight from the archive.
    Returns ({key: (hours since the day's midnight UTC, values in model unit, label, kind)}, keys looked for)."""
    day = archive_day(archive_root, date)
    out: ObsDict = {}
    keys = list(keys or OBS_MAP)
    if not day:
        return out, keys
    try:
        from enz_builder import read_nc_variable
    except Exception:
        return out, keys
    midnight = _midnight_of(date)
    for key in keys:
        suffix, var, _unit, factor, label, kind = OBS_MAP[key]
        t_s = vals = None
        for suf, vname in [(suffix, var)] + OBS_OLD_NAMES.get(key, []):      # current name, then the older archive naming
            p = os.path.join(day, f"{date}.{suf}")
            if not os.path.isfile(p):
                continue
            try:
                t_s, vals = read_nc_variable(p, vname)           # seconds since 2000-01-01 (ENZ epoch)
                break
            except Exception:
                continue
        if t_s is None:
            continue
        hours = (np.asarray(t_s, dtype=float) + ENZ_EPOCH_UNIX - midnight) / 3600.0
        v = np.asarray(vals, dtype=float) * factor
        ok = np.isfinite(v)
        if ok.sum() == 0:
            continue
        out[key] = (hours[ok], v[ok], label, kind)
    return out, keys


def expected_koh(spec: ExperimentSpec, hours: np.ndarray, cols: Dict[str, np.ndarray],
                 obs_hours: Optional[np.ndarray] = None, obs_koh: Optional[np.ndarray] = None) -> Tuple[np.ndarray, List[Tuple[str, float, float]], float, str]:
    """First-order expectation of the measured OH reactivity from the injections alone:
    background (median of the measurement in the 15 min — or the last points of the hour — before the first
    injection adding ≥ 0.05 s⁻¹, else 0)
    + Σ kOH of each injection, ramped in linearly over its pulse and diluted afterwards with the
    chamber flow of the ENZ ('dilute' m³/h ÷ chamber volume).  No chemistry: products, OH loss
    of the injected compounds and chamber sources are ignored.
    Returns (expected kOH on `hours`, [(name, start hour, kOH s-1)], background, note)."""
    import oh_reactivity as R
    import injection_calc as ic
    c = spec.conditions
    T = float(np.nanmean(cols["t"])) if "t" in cols and np.isfinite(cols["t"]).any() else float(c.temperature_K or 298.0)
    p = float(c.pressure_hPa or 1013.25)
    vol = getattr(ic.chamber(c.chamber_name), "volume_m3", 276.0) or 276.0
    events = []
    for ev in spec.injections:
        if not ev.start or not ev.products.strip():
            continue
        k = 0.0
        for prod, ppb in [(ev.products, ev.total_ppb)] + list(ev.extra_channels):
            v, _missing = R.injection_koh(prod, ppb, T, p)
            k += v or 0.0
        if k > 0:
            hh, mm = ev.start.split(":")[:2]
            events.append((ev.name, int(hh) + int(mm) / 60.0, k, max(float(ev.duration_min or 1.0), 1e-3) / 60.0))
    # dilution: d ln C / dt = -flow / V   (flow in m³/h → per hour)
    flow = np.asarray(cols.get("dilute", np.zeros(len(hours))), dtype=float)
    rate_h = np.where(np.isfinite(flow), flow, 0.0) / vol
    dt = np.diff(hours, prepend=hours[0])
    cum = np.cumsum(rate_h * dt)                        # ∫ dilution dt from the start, dimensionless
    exp_k = np.zeros(len(hours))
    for _name, h0, k, dur in events:
        ramp = np.clip((hours - h0) / dur, 0.0, 1.0)
        i0 = int(np.searchsorted(hours, h0 + dur))
        decay = np.exp(-(cum - (cum[min(i0, len(cum) - 1)] if len(cum) else 0.0)))
        exp_k += k * ramp * np.where(hours >= h0 + dur, decay, 1.0)
    background, note = 0.0, "background 0 (no measurement before the first injection)"
    if obs_hours is not None and obs_koh is not None and events:
        big = [h for _n, h, k, _d in events if k >= 0.05]          # injections that change kOH noticeably (not humidification impurities)
        first = min(big or [h for _n, h, _k, _d in events])
        ok = np.isfinite(obs_koh) & (obs_hours < first)
        sel = ok & (obs_hours >= first - 0.25)
        window = "15 min"
        if sel.sum() < 3:                                # sparse measurement (e.g. every 10 min): the last points within the hour
            idx = np.where(ok & (obs_hours >= first - 1.0))[0][-3:]
            sel = np.zeros(len(obs_hours), bool); sel[idx] = True; window = "hour"
        if sel.any():
            background = float(np.median(obs_koh[sel]))
            note = (f"background {background:.2g} s⁻¹ = median of {int(sel.sum())} measurement(s) in the {window} before "
                    f"{int(first):02d}:{int(round(first % 1 * 60)):02d}")
    return exp_k + background, [(n, h, k) for n, h, k, _d in events], background, note


# ── one file per experiment day: everything the plots compare against, extracted once ──
OBS_FILE_SUFFIX = ".observations.csv"


def observation_file(folder: str, date: str) -> str:
    """`<folder>/<date>.observations.csv` — the extracted observations kept next to the results."""
    return os.path.join(folder, f"{date}{OBS_FILE_SUFFIX}")


def write_observation_file(path: str, obs: ObsDict, keys_checked: List[str], date: str, source: str) -> str:
    """Long table `hour_utc,key,value` with a comment header that records where each series came from.
    Plain CSV so it opens in Excel/pandas and can travel with the experiment folder."""
    import datetime as _dt
    lines = [f"# JUICE observations for {date}, extracted {_dt.datetime.now():%Y-%m-%d %H:%M} from {source}",
             "# hour_utc = hours since 00:00 UTC of that day; value in the model unit of the kind column",
             "# checked: " + ",".join(keys_checked)]
    for key in keys_checked:
        if key in obs:
            suffix, var, unit, factor, label, kind = OBS_MAP[key]
            lines.append(f"# series: {key}|{label}|{kind}|{date}.{suffix}|{var}|{unit}|{factor:g}|{len(obs[key][0])}")
    lines.append("hour_utc,key,value")
    for key in keys_checked:
        if key not in obs:
            continue
        h, v, _l, _k = obs[key]
        lines.extend(f"{hh:.6f},{key},{vv:.8g}" for hh, vv in zip(h, v))
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(lines) + "\n")
    return path


def read_observation_file(path: str) -> Tuple[ObsDict, List[str], str]:
    """Inverse of write_observation_file → (obs, keys checked, header line)."""
    import csv
    checked: List[str] = []; meta: Dict[str, Tuple[str, str]] = {}; header = ""
    data: Dict[str, Tuple[List[float], List[float]]] = {}
    with open(path, encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                if not header:
                    header = line[1:].strip()
                if line.startswith("# checked:"):
                    checked = [k for k in line.split(":", 1)[1].strip().split(",") if k]
                elif line.startswith("# series:"):
                    parts = line.split(":", 1)[1].strip().split("|")
                    if len(parts) >= 3:
                        meta[parts[0]] = (parts[1], parts[2])
                continue
            if line.startswith("hour_utc"):
                continue
            try:
                hs, key, vs = line.rstrip("\n").split(",", 2)
                data.setdefault(key, ([], []))
                data[key][0].append(float(hs)); data[key][1].append(float(vs))
            except ValueError:
                continue
    obs: ObsDict = {}
    for key, (h, v) in data.items():
        label, kind = meta.get(key, (OBS_MAP[key][4] if key in OBS_MAP else key, OBS_MAP[key][5] if key in OBS_MAP else "ppb"))
        obs[key] = (np.asarray(h), np.asarray(v), label, kind)
    return obs, checked, header


def load_observations(archive_root: str, date: str, keys: Optional[List[str]] = None,
                      midnight_unix: Optional[float] = None, cache_dir: str = "",
                      status: Optional[List[str]] = None) -> ObsDict:
    """{key: (hours since `midnight_unix`, values in model unit, label, kind)} for the keys available on that day.

    With `cache_dir`, the day's observation file there is used when it exists and covers the keys asked
    for; otherwise the archive is read and the file written (the archive is only needed once per day,
    and the plots work without the Z: drive afterwards).  `status` collects a one-line description."""
    keys = list(keys or OBS_MAP)
    obs: ObsDict = {}
    note = ""
    cache = observation_file(cache_dir, date) if cache_dir and date else ""
    if cache and os.path.isfile(cache):
        try:
            cached, checked, header = read_observation_file(cache)
            if all(k in checked for k in keys):
                obs = {k: v for k, v in cached.items() if k in keys}
                note = f"from {os.path.basename(cache)} ({header.split(', extracted')[0].replace('JUICE observations for ', '')})"
            else:
                note = "observation file lacks " + ", ".join(k for k in keys if k not in checked) + " — re-extracted"
        except Exception as e:
            note = f"observation file unreadable ({e}) — re-extracted"
    if not obs and not note.startswith("from "):
        day = archive_day(archive_root, date) if date else ""
        if day:
            obs, checked = extract_observations(archive_root, date, keys if not cache else sorted(set(keys) | set(OBS_MAP)))
            if cache:
                try:
                    write_observation_file(cache, obs, checked, date, day)
                    note = (note + "; " if note else "") + f"archive read, saved {os.path.basename(cache)}"
                except Exception as e:
                    note = (note + "; " if note else "") + f"archive read (could not save the observation file: {e})"
            else:
                note = "archive read"
            obs = {k: v for k, v in obs.items() if k in keys}
        else:
            note = f"no archive folder for {date or 'this day'}"
    if midnight_unix is not None and date and obs:
        shift = (_midnight_of(date) - midnight_unix) / 3600.0
        if abs(shift) > 1e-9:
            obs = {k: (h + shift, v, l, kd) for k, (h, v, l, kd) in obs.items()}
    if status is not None:
        status.append(note)
    return obs


def cs_species(path: str) -> List[str]:
    """Column names of a result file (header only, nothing else is read)."""
    from column_file_reader import _parse_column_defs
    return list(_parse_column_defs(path))


def read_columns(path: str, names: List[str]) -> Tuple[np.ndarray, Dict[str, np.ndarray]]:
    """A few named columns of a FACSIMILE-style column file (.cs/.pd/.kv/.enz) without parsing the
    rest: (time, {name: values}).  Comma- or space-separated rows, missing-E exponents repaired.
    A 47 MB .pd with 3800 columns gives its kOH column in ~0.3 s instead of 2 s."""
    import re as _re
    from column_file_reader import _parse_column_defs
    cols = _parse_column_defs(path)
    idx = {n: i for i, n in enumerate(cols)}
    want = [(n, idx[n]) for n in names if n in idx]
    if not want:
        return np.zeros(0), {}
    t: List[float] = []; data: Dict[str, List[float]] = {n: [] for n, _i in want}
    fix = _re.compile(r"(\d)([-+])(\d)")

    def num(tok: str) -> float:
        try:
            return float(tok)
        except ValueError:
            return float(fix.sub(r"\1E\2\3", tok))

    in_data = False
    with open(path, encoding="latin-1") as fh:
        for line in fh:
            s = line.strip()
            if not in_data:
                in_data = s.startswith("&&&")
                continue
            if not s:
                continue
            toks = (s.replace(",", " ") if "." in s else s.replace(",", ".")).split()
            try:
                row = [num(toks[i]) for _n, i in want]
            except (IndexError, ValueError):
                continue                               # incomplete row (e.g. a file still being written)
            t.append(num(toks[0]))
            for (n, _i), v in zip(want, row):
                data[n].append(v)
    tt = np.asarray(t, dtype=float)
    keep = np.r_[tt[1:] != tt[:-1], True] if len(tt) else np.zeros(0, bool)   # duplicate boundary rows: keep last
    return tt[keep], {n: np.asarray(v, dtype=float)[keep] for n, v in data.items()}


def enz_next_to(cs_path: str) -> str:
    """The .enz the run was driven with (same stem next to the result), '' if absent."""
    p = os.path.splitext(cs_path)[0] + ".enz"
    return p if os.path.isfile(p) else ""


def landing_panels(spec: Optional[ExperimentSpec], columns: List[str]) -> List[dict]:
    """Panels of the Overview tab: title, model species (present in `columns`), unit kind,
    observation keys, log axis flag.  Built from the experiment description when there is one."""
    cols = set(columns)
    have = lambda names: [s for s in names if s in cols]
    panels: List[dict] = []
    inj = have(injected_species(spec)) if spec else []
    if inj:
        panels.append(dict(title="Injected compounds", species=inj[:8], kind="ppb", obs=[], log=False))
    panels.append(dict(title="NOx and ozone", species=have(["NO", "NO2", "O3"]), kind="ppb", obs=["NO", "NO2", "O3"], log=False))
    panels.append(dict(title="Radicals", species=have(["HO2", "RO2", "OH"]), kind="cm-3", obs=["HO2", "RO2", "OH"], log=False, twin="OH"))
    panels.append(dict(title="HONO and CO", species=have(["HONO", "CO"]), kind="ppb", obs=["HONO", "CO"], log=False, twin="CO"))
    panels.append(dict(title="OH reactivity", species=[], kind="s-1", obs=["kOH"], log=False, koh=True))
    panels.append(dict(title="Chamber sources and conditions", species=[], kind="sources", obs=["jNO2"], log=False, sources=True))
    return [p for p in panels if p["species"] or p.get("sources") or p.get("koh")]
