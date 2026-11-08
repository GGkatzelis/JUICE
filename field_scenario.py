"""
field_scenario.py — a field case → ENZ + EASY → scenario JSON, the way experiment_spec does it
for a chamber experiment.  Everything the field adds to the chemistry is written as ordinary
EASY statements, so FACSIMILE can run the same files:

    conditions     one ENZ row per output step (step_s), conditions interpolated linearly, the
                   emissions held from one table row to the next; T, P, H2O, sz, hv from the table
                   or computed (sz from lat/lon per row; hv = 1, the sun is always "on")
    photolysis     measured jNO2 (and jO1D, jHONO, jH2O2, jHCHO) as in a chamber: the clear-sky
                   frequencies of the mechanism are scaled by the measured jNO2; without jNO2 the
                   clear-sky values for time and position
    emission       QE_X = INPUT(QE_X)   QE_X --> X   k = CONST(1.)        column in molec cm-3 s-1
    dilution       the mechanism's X + DIL --> reactions with DILUTE = CONST(KDIL + k_generic),
                   KDIL = INPUT(KDIL) in 1/s from field_spec.dilution_rate
    background     B_X = INPUT(B_X)     B_X --> X    k = CONST(KDIL*1e-9*M)      column in ppb
    held species   C_X = INPUT(C_X)     X = CONST(C_X*1e-9*M)                     column in ppb
                   (idl_parser / scenario_runner hold such a species per ENZ row, as FACSIMILE does)
    deposition     X -->   k = CONST(v/(HPBL*100.))  or  CONST(1/tau)
    spin-up        constrained site: the table is repeated before the start, period = whole days

    res = write_case(spec, table, out_dir)        # <name>.enz, <name>.easy, chemistry/, <name>.field.json
    path = prepare(spec, table, out_dir, include_dir, toml_path, solver)   # + scenario JSON
    python field_scenario.py case.field.json [--out DIR] [--run] [--solver fast|precise]
"""
from __future__ import annotations

import argparse
import math
import os
import subprocess
import sys
from typing import Dict, List, Optional, Tuple

import numpy as np

import field_spec as fs
import mechanism_registry as mreg
from easy_builder import J_INPUTS, MeasuredRow, default_form, render_form
from enz_builder import compute_solar_zenith, write_enz
from experiment_spec import CLEAR_SKY_JNO2, ExperimentSpec, ensure_mechanism, J_ENZ_NAME
from lamp_photolysis import MEASURED_J

HERE = os.path.dirname(os.path.abspath(__file__))
# table photolysis name (field_spec strips the leading j) → EASY variable of the chamber block
J_TABLE_TO_EASY = {"NO2": "jno2", "O1D": "jo1d", "HONO": "jhono", "H2O2": "jh2o2", "HCHO_r": "jhcho_r", "HCHO_m": "jhcho_m",
                   "HCHO_R": "jhcho_r", "HCHO_M": "jhcho_m"}
COL_KDIL, COL_HPBL = "KDIL", "HPBL"


# ── time grid and interpolation ──────────────────────────────────────────────────────────────────
def spin_up_period(table: fs.FieldTable) -> float:
    """The table is repeated with a period of whole days (a diurnal table of 24 h repeats daily)."""
    span = float(table.times_s[-1] - table.times_s[0])
    return max(1.0, math.ceil(span / 86400.0 - 1e-9)) * 86400.0


def field_grid(spec: fs.FieldSpec, table: fs.FieldTable) -> Tuple[np.ndarray, np.ndarray]:
    """(grid in ENZ seconds, the table time each grid row takes its values from).  Spin-up rows lie
    before the table and map back into it with the repeat period."""
    t0, t1 = float(table.times_s[0]), float(table.times_s[-1])
    step = float(spec.step_s)
    n_copies = 0
    period = spin_up_period(table)
    if spec.run_type == "constrained" and spec.spin_up_days > 0:
        n_copies = int(math.ceil(spec.spin_up_days * 86400.0 / period))
    start = t0 - n_copies * period
    n = int(math.floor((t1 - start) / step + 1e-9)) + 1
    grid = start + step * np.arange(n)
    src = np.where(grid < t0, t0 + np.mod(grid - t0, period), grid)
    return grid, src


def _linear(times: np.ndarray, values: np.ndarray, at: np.ndarray) -> np.ndarray:
    """Linear interpolation through the finite values, the ends held outside; NaN when nothing is finite."""
    ok = np.isfinite(values) & np.isfinite(times)
    if not ok.any():
        return np.full(len(at), np.nan)
    return np.interp(at, times[ok], values[ok])


def _step_hold(times: np.ndarray, values: np.ndarray, at: np.ndarray) -> np.ndarray:
    """The value of the latest table row at or before each time (the first row before the table)."""
    ok = np.isfinite(values) & np.isfinite(times)
    if not ok.any():
        return np.zeros(len(at))
    t, v = times[ok], values[ok]
    i = np.clip(np.searchsorted(t, at, side="right") - 1, 0, len(t) - 1)
    return v[i]


def gap_report(table: fs.FieldTable, name: str, values: np.ndarray) -> Optional[str]:
    """How empty cells of a column are bridged, for the notes."""
    bad = ~np.isfinite(values)
    if not bad.any():
        return None
    ok = np.flatnonzero(~bad)
    if not len(ok):
        return f"{name}: no values at all"
    longest = 0.0
    for a, b in zip(ok[:-1], ok[1:]):
        if b > a + 1:
            longest = max(longest, float(table.times_s[b] - table.times_s[a]))
    ends = int(ok[0]) + int(len(values) - 1 - ok[-1])
    return (f"{name}: {int(bad.sum())} empty cell(s) bridged by interpolation"
            + (f", longest gap {longest / 60.0:.0f} min" if longest else "") + (f", {ends} at the ends held" if ends else ""))


# ── ENZ columns ──────────────────────────────────────────────────────────────────────────────────
def enz_columns(spec: fs.FieldSpec, table: fs.FieldTable) -> Tuple[np.ndarray, Dict[str, np.ndarray], List[str]]:
    grid, src = field_grid(spec, table)
    tt = table.times_s
    c = table.conditions
    notes: List[str] = []
    cols: Dict[str, np.ndarray] = {}

    def lin(name: str, values: np.ndarray) -> np.ndarray:
        g = gap_report(table, name, values)
        if g:
            notes.append(g)
        return _linear(tt, values, src)

    if "sz" in c and np.isfinite(c["sz"]).all():
        cols["sz"] = lin("sz", c["sz"])
    else:
        if "lat" in c and "lon" in c:
            lat, lon = lin("lat", c["lat"]), lin("lon", c["lon"])
        else:
            lat, lon = np.full(len(grid), float(spec.lat)), np.full(len(grid), float(spec.lon))
        cols["sz"] = np.array([compute_solar_zenith([g], la, lo)[0] for g, la, lo in zip(grid, lat, lon)])
    cols["hv"] = np.ones(len(grid))
    cols["t"] = lin("T", c["t"])
    cols["P"] = lin("p", c["p"])
    cols["H2O"] = lin("H2O", c["H2O"])
    kd = fs.dilution_rate(table, spec)
    cols[COL_KDIL] = _step_hold(tt, kd, src) if spec.dilution in ("pbl", "tracer") else lin("kdil", kd)
    if "h" in c:
        cols[COL_HPBL] = lin("h", c["h"])
    for sp, rate in fs.emission_rates(table, spec).items():
        cols[f"QE_{sp}"] = _step_hold(tt, rate, src)
    for sp, v in table.background.items():
        cols[f"B_{sp}"] = lin(f"B_{sp}", v)
    for sp, v in table.constraints.items():
        cols[f"C_{sp}"] = np.maximum(lin(f"C_{sp}", v), 0.0)
    for name, v in table.photolysis.items():
        easy = J_TABLE_TO_EASY.get(name)
        if easy is None:
            notes.append(f"j{name}: no photolysis of that name in the chamber block, column written but not used")
            cols[f"j{name}"] = np.maximum(lin(f"j{name}", v), 0.0)
        else:
            cols[J_ENZ_NAME[easy]] = np.maximum(lin(f"j{name}", v), 0.0)
    if spec.run_type == "constrained" and spec.spin_up_days > 0:
        n_pre = int((grid < tt[0]).sum())
        notes.append(f"spin-up: {n_pre} rows ({n_pre * spec.step_s / 3600.0:.1f} h) before the table, repeating it with a "
                     f"period of {spin_up_period(table) / 86400.0:.0f} day(s)")
    return grid, cols, notes


# ── EASY ─────────────────────────────────────────────────────────────────────────────────────────
def chemistry_spec(spec: fs.FieldSpec) -> ExperimentSpec:
    """The chamber generator's description of the chemistry (its mechanism subset code is reused)."""
    x = ExperimentSpec(name=spec.name or "field_case")
    for k, v in (spec.chemistry or {}).items():
        if hasattr(x, k):
            setattr(x, k, v)
    x.reaction_adjustments = [tuple(a) for a in (x.reaction_adjustments or [])]
    return x


def field_statements(spec: fs.FieldSpec, table: fs.FieldTable, cols: Dict[str, np.ndarray]) -> List[str]:
    k_generic = 1.0 / (spec.generic_lifetime_h * 3600.0) if spec.run_type == "constrained" else 0.0
    out = [f"; FIELD CASE ({fs.RUN_TYPES[spec.run_type]}): what the field adds to the chemistry, written for FACSIMILE and JUICE alike",
           f"{COL_KDIL} = INPUT({COL_KDIL})   ; first-order dilution of every species, 1/s ({fs.DILUTION_SOURCES[spec.dilution]})"]
    if COL_HPBL in cols:
        out.append(f"{COL_HPBL} = INPUT({COL_HPBL})   ; boundary-layer height, m")
    if k_generic:
        out.append(f"; every species not held to a measurement is lost with a lifetime of {spec.generic_lifetime_h:g} h (in DILUTE)")
    emitted = sorted(k[3:] for k in cols if k.startswith("QE_"))
    if emitted:
        out.append("; emissions, molec cm-3 s-1 (one ENZ column per species)")
        for sp in emitted:
            out += [f"QE_{sp} = INPUT(QE_{sp})", f"QE_{sp} --> {sp}", f"k[QE_{sp} --> {sp}] = CONST(1.)"]
    backgrounds = sorted(k[2:] for k in cols if k.startswith("B_"))
    if backgrounds:
        out.append("; background mixing ratios (ppb) entrained with the dilution: + KDIL * background")
        for sp in backgrounds:
            out += [f"B_{sp} = INPUT(B_{sp})", f"B_{sp} --> {sp}", f"k[B_{sp} --> {sp}] = CONST({COL_KDIL}*1e-9*M)"]
    held = sorted(k[2:] for k in cols if k.startswith("C_"))
    if held:
        out.append("; species held to the measurement (ppb), set at every ENZ row")
        for sp in held:
            out += [f"C_{sp} = INPUT(C_{sp})", f"{sp} = CONST(C_{sp}*1e-9*M)"]
    if spec.deposition_velocity_cm_s or spec.deposition_lifetime_h:
        out.append("; deposition as a first-order loss")
        for sp, v in spec.deposition_velocity_cm_s.items():
            out += [f"{sp} -->", f"k[{sp} -->] = CONST({v:g}/({COL_HPBL}*100.))   ; {v:g} cm/s over the boundary layer"]
        for sp, tau in spec.deposition_lifetime_h.items():
            out += [f"{sp} -->", f"k[{sp} -->] = CONST(1./{tau * 3600.0:g})   ; lifetime {tau:g} h"]
    return out


def easy_text(spec: fs.FieldSpec, table: fs.FieldTable, chem: ExperimentSpec, cols: Dict[str, np.ndarray]) -> str:
    f = default_form()
    f.title = f"{spec.name}: {fs.RUN_TYPES[spec.run_type]}" + (f", {spec.site}" if spec.site else "")
    f.header_comments = [f"JUICE field case {spec.name} ({fs.RUN_TYPES[spec.run_type]}), written by field_scenario.py",
                         f"input table: {os.path.basename(spec.table)}; mechanism: {fs.chemistry_label(spec)}"]
    f.mech_stem, f.mech_incl_suffix = chem.mechanism_stem, chem.mech_incl_suffix
    f.use_anorg, f.use_j_clear_sky = chem.use_anorg, chem.use_j_clear_sky_incl
    f.includes_other = list(chem.extra_includes)
    # nothing a chamber does on its own
    f.walldep = f.hono_source = f.hcho_source = f.bgreac = f.dark_hono = f.light_source = False
    f.tracer_dilution = []
    k_generic = 1.0 / (spec.generic_lifetime_h * 3600.0) if spec.run_type == "constrained" else 0.0
    f.dilution, f.dilute_expr, f.dil_value = True, f"CONST({COL_KDIL}+{k_generic:.6e})", "CONST(1.)"
    # photolysis: measured jNO2 scales the clear-sky frequencies; otherwise clear sky
    # No saphir.jmeasured_02.incl: that chamber include re-parametrises a few photolysis reactions for the
    # roof-open hours and its NO3CH2CHO fit (exp(1.39/cos)) runs away at twilight, when the sun passes the
    # horizon during a field run.  The mechanism's own J<n> = PHTFKT(...) are used instead: j_clear_sky.incl
    # scales them with the measured jNO2 (ratio capped at 1), and a measured frequency replaces its J<n>.
    f.jmeasured_include = ""
    f.j_overrides = {}
    j_lines: List[str] = []
    measured = [(j, J_ENZ_NAME[j]) for j, _col, _r in J_INPUTS if J_ENZ_NAME[j] in cols]
    if any(j == "jno2" for j, _c in measured):
        f.j_rows = [MeasuredRow(j, "enz", col) for j, col in measured]
        j_lines = ["; measured photolysis frequencies replace the mechanism's J<n>; every other J<n> is the clear-sky value "
                   "scaled by jNO2 / jNO2(clear sky), at most 1"] + [f"J{MEASURED_J[j]} = CONST({j})   ; {col}" for j, col in measured]
    else:
        f.j_rows = [MeasuredRow("jno2", "const", "", CLEAR_SKY_JNO2)]
        j_lines = ["; no measured jNO2: clear-sky photolysis for time and position (jno2 = J4 clear sky, so the scaling is 1)"]
    f.setting_rows = [
        MeasuredRow("T", "enz", "t", comment="temperature in K"),
        MeasuredRow("P", "input", "P", comment="pressure in hPa"),
        MeasuredRow("SZA", "enz", "sz", comment="solar zenith angle"),
        MeasuredRow("SZ", "enz", "sz", comment="solar zenith angle"),
        MeasuredRow("HV", "enz", "hv", comment="the sun: always on, the zenith angle does the rest"),
        MeasuredRow("FLOW", "const", "", "0.", "no chamber flow"),
        MeasuredRow("H2O", "enz", "H2O", comment="mixing ratio of water"),
    ]
    f.enz_file = f"{spec.name}.enz"
    f.hmax, f.outstep = chem.hmax, chem.outstep
    f.initials = [("*", f"CONST({chem.initial_all})", "")] + [(sp, f"CONST(M*1e-9*{ppb:g})", "") for sp, ppb in spec.initials.items()]
    f.injections, f.chain_reactions = [], []
    f.extra_statements = []
    try:
        f.extra_statements += mreg.provenance_lines(chem.base_mechanism, list(chem.addons))
    except Exception as e:
        f.extra_statements.append(f"; CHEMISTRY: provenance unavailable ({e})")
    f.extra_statements += j_lines
    f.extra_statements += field_statements(spec, table, cols)
    for rxn, rate, cmt in chem.reaction_adjustments:
        if str(rxn).strip():
            f.extra_statements.append(f"k[{str(rxn).strip()}] = CONST({(rate or '0').strip()})" + (f"   ;{cmt.strip()}" if cmt.strip() else ""))
    return render_form(f)


# ── writing and running ──────────────────────────────────────────────────────────────────────────
def write_case(spec: fs.FieldSpec, table: fs.FieldTable, out_dir: str) -> dict:
    """<name>.enz, <name>.easy, the mechanism subset in chemistry/ and a copy of the case.
    Returns {'enz', 'easy', 'mechanism', 'notes', 'spec_path', 'chem'} like experiment_spec.write_experiment."""
    os.makedirs(out_dir, exist_ok=True)
    chem = chemistry_spec(spec)
    grid, cols, notes = enz_columns(spec, table)
    chem.conditions.j_source = "measured" if "jNO2" in cols else "clear_sky"
    res = {"mechanism": ensure_mechanism(chem, out_dir), "chem": chem}
    enz = os.path.join(out_dir, f"{spec.name}.enz")
    write_enz(enz, grid, cols)
    easy = os.path.join(out_dir, f"{spec.name}.easy")
    with open(easy, "w", encoding="latin-1", errors="replace", newline="\n") as fh:
        fh.write(easy_text(spec, table, chem, cols))
    res.update({"enz": enz, "easy": easy, "notes": notes, "grid": grid, "cols": cols,
                "spec_path": fs.save_spec(spec, os.path.join(out_dir, spec.name + fs.SPEC_SUFFIX))})
    return res


def prepare(spec: fs.FieldSpec, table: fs.FieldTable, out_dir: str, include_dir: str, toml_path: str = "",
            solver: Optional[dict] = None, export_kv_pd: bool = True) -> Tuple[str, dict]:
    """Write the case and its scenario JSON.  Returns (json path, write_case result)."""
    from scenario_prep import prepare_scenario, solver_from_toml
    res = write_case(spec, table, out_dir)
    extra = list((spec.chemistry or {}).get("include_dirs") or [])
    path = prepare_scenario(res, include_dir, export_kv_pd, toml_path, extra_dirs=extra, easy_version="",
                            solver=solver if solver is not None else solver_from_toml(toml_path))
    return path, res


def run(json_path: str, toml_path: str = "", log_path: str = "") -> int:
    """Run scenario_runner on a prepared case in its own process; the output goes to log_path or the console."""
    args = [sys.executable, "-u", os.path.join(HERE, "scenario_runner.py"), json_path] + (["--toml", toml_path] if toml_path else [])
    env = {**os.environ, "PYTHONIOENCODING": "utf-8", "PYTHONUTF8": "1"}
    if log_path:
        with open(log_path, "w", encoding="utf-8") as fh:
            return subprocess.run(args, stdout=fh, stderr=subprocess.STDOUT, env=env, cwd=HERE).returncode
    return subprocess.run(args, env=env, cwd=HERE).returncode


def main(argv=None) -> int:
    from scenario_prep import load_paths, pick_toml, solver_from_toml
    ap = argparse.ArgumentParser(description="field case → ENZ + EASY + scenario JSON (+ run)")
    ap.add_argument("case", help="<name>.field.json")
    ap.add_argument("--out", default="", help="output folder (default: <output_root>/<name>)")
    ap.add_argument("--run", action="store_true", help="also run the solver")
    ap.add_argument("--solver", choices=["fast", "precise"], default=None)
    ap.add_argument("--no-kv-pd", action="store_true", help="skip the reaction-rate export")
    ap.add_argument("-t", "--toml", default="")
    a = ap.parse_args(argv)
    toml = pick_toml(["-t", a.toml] if a.toml else [])
    paths = load_paths(toml)
    spec = fs.load_spec(a.case)
    table, rep = fs.check(spec)
    print(rep.text())
    if not rep.ok or table is None:
        return 2
    out = a.out or os.path.join(paths["output_root"], spec.name)
    path, res = prepare(spec, table, out, paths["include_dir"], toml, solver_from_toml(toml, a.solver), not a.no_kv_pd)
    for n in res["notes"]:
        print("note:", n)
    print("prepared", path)
    return run(path, toml) if a.run else 0


if __name__ == "__main__":
    sys.exit(main())
