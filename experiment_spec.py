"""
experiment_spec.py — experiment description → ENZ + EASY (Goal 2, step 3)
========================================================================

The user-facing layer of the JUICE scenario builder.  An `ExperimentSpec`
describes an experiment the way an experimenter thinks about it:

    conditions   : where and when (date, time window, roof schedule, T, P,
                   humidity, flow, J-value source), from the SAPHIR archive
                   or as manual constants
    chemistry    : the precursor compounds (MCM names) — the mechanism subset
                   is extracted from the full MCM by mcm_converter — or an
                   existing include set (mech stem)
    injections   : one entry per injection event: what, how much (total ppb),
                   when, how long; optional impurity channels and chains
    chamber      : what the SAPHIR chamber adds on its own (HONO/HCHO wall
                   sources, background reactivity, dilution, wall deposition,
                   dark HONO, NO2 light source), all with the IDL defaults
    initials, FACSIMILE parameters, expert overrides

From one spec two files are generated consistently (nothing is typed twice):

    <name>.enz   — time grid, sz, hv, dilute, t, H2O, one normalised pulse
                   column per injection (value 1/n over the n injection rows,
                   integral 1), optional measured columns
    <name>.easy  — through easy_builder.EasyForm → easy_scene.write_scene,
                   i.e. exactly the IDL layout; the injection rate holds the
                   TOTAL injected ppb: k[q --> X] = CONST(M*1e-9/60*ppb)

Golden tests (`python experiment_spec.py --selftest`) rebuild scene_05 and
scene_09 from specs and require statement-equivalent EASY files and identical
injection / roof columns in the ENZ.
"""

from __future__ import annotations

import os
import re
import sys
import math
from dataclasses import dataclass, field, asdict
from datetime import datetime, timedelta, timezone
from typing import Dict, List, Optional, Tuple

import numpy as np

from easy_builder import (EasyForm, MeasuredRow, InjectionRow, default_form, render_form,
                          J_INPUTS, HONO_RATE, HCHO_RATE, Y_RATE, DILUTE_DEF)
from easy_scene import statements_equivalent

ENZ_EPOCH = datetime(2000, 1, 1, 0, 0, 0, tzinfo=timezone.utc)
CLEAR_SKY_JNO2 = "PHTFKT_CS(1.165D-02,0.244,-0.267)"

# default ENZ column per injected compound (IDL convention); others → Q_<name>
STANDARD_COLUMNS = {"NO": "QNO", "NO2": "QNO2", "O3": "QO3", "H2O": "QH2O", "CO2": "QCO2",
                    "CL2": "QCl2", "SEED": "Q_seed", "CO": "QCO"}
# humidification carries chamber impurities (ppb per humidification, scene_05/09 defaults)
H2O_IMPURITY_DEFAULTS = [("HONO+XNOy", "0.002"), ("CO", "0.5"), ("NO2+XNOy", "0.00"), ("NO+XNOy", "0.000")]


# ─────────────────────────────────────────────────────────────────────────────
# Specification
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class Conditions:
    date: str = "2024-07-09"              # YYYY-MM-DD (UTC)
    start: str = "04:00"                  # HH:MM
    end: str = "21:59"
    step_s: int = 60
    source: str = "manual"                # manual | archive
    archive_root: str = r"Z:\IEK8-SAPHIR\data"
    temperature_K: float = 298.0
    temperature_mode: str = "enz"         # enz: T = CONST(T_in) | const: T = CONST(<temperature_K>)
    pressure_hPa: float = 1006.5
    h2o_percent: float = 1.0              # ENZ H2O column (vol-%); EASY: H2O = CONST(1e-2*M*H2O_in)
    flow_m3h: float = 6.0                 # ENZ dilute column (m^3/h); EASY: FLOW = INPUT(dilute)
    roof_open: List[Tuple[str, str]] = field(default_factory=list)   # [(open HH:MM, close HH:MM)]
    j_source: str = "measured"            # measured (ENZ j columns) | clear_sky (jno2 = PHTFKT_CS, rest from mechanism)
    constrain_species: List[str] = field(default_factory=list)       # e.g. ["NO2", "O3"] → held to ENZ columns


@dataclass
class InjectionEvent:
    name: str                             # label: "NO", "ethanol", "VCP2", ...
    products: str                         # EASY product string: "NO", "C2H5OH", "0.9*VOCin1+0.1*VOCin2", mixture;
                                          # empty → event is recorded in the ENZ only (not modelled, no EASY block)
    total_ppb: str = "0"                  # total injected amount (number or numeric expression)
    start: str = ""                       # HH:MM (empty → no ENZ pulse, column stays 0)
    duration_min: float = 1.0
    column: str = ""                      # ENZ switch column (default from name)
    switch: str = ""                      # EASY switch variable (default from column)
    extra_channels: List[Tuple[str, str]] = field(default_factory=list)   # (products, total_ppb): impurities / second product
    chain: List[Tuple[str, str, str]] = field(default_factory=list)       # follow-up reactions (lhs, rhs, rate expr)
    comment: str = ""

    def enz_column(self) -> str:
        if self.column:
            return self.column
        return STANDARD_COLUMNS.get(self.name.upper(), f"Q_{self.name}")

    def easy_switch(self) -> str:
        if self.switch:
            return self.switch
        col = self.enz_column()
        if re.fullmatch(r"Q[A-Z0-9]+", col):          # QNO → qNO, QH2O → qH2O
            return "q" + col[1:]
        return col                                     # Q_ethanol / qVCP2 → same name


@dataclass
class ChamberEffects:
    walldep: bool = True
    walldep_value: str = "1."
    hono: bool = True
    hono_lhs: str = "SHONO+HV"
    hono_rhs: str = "HONO+XNOy"
    hono_const: str = "2.6*0.7"
    hono_rate: str = HONO_RATE
    hcho: bool = True
    hcho_lhs: str = "qSHCHO+HV"
    hcho_rhs: str = "HCHO"
    hcho_const_name: str = "SHCHO"
    hcho_const: str = "1"
    hcho_rate: str = HCHO_RATE
    bgreac: bool = True
    y_name: str = "qY"
    y_product: str = "CH3O2"
    y_ppb: str = "250"
    y_expr: str = ""                      # full expression overrides y_ppb
    y_rate: str = Y_RATE
    dilution: bool = True
    dilute_expr: str = DILUTE_DEF
    dil_value: str = "CONST(1.)"
    dil_flow_const: str = ""
    tracer_dilution: List[str] = field(default_factory=lambda: ["XNOy", "NA"])
    dark_hono: bool = True
    dark_hono_rhs: str = "HONO+XNOy"
    dark_hono_rate: str = "1e-012/10"
    light_source: bool = False
    light_source_rhs: str = "NO2+XNOy"
    light_source_rate: str = "jno2/5e-3/3600*HV_in/1.6*0.55"
    light_source_initial: str = "M*1e-9*3.5*0"


@dataclass
class ExperimentSpec:
    name: str = "experiment"
    title: str = ""
    conditions: Conditions = field(default_factory=Conditions)
    precursors: List[str] = field(default_factory=list)      # MCM names → mechanism subset (backend)
    mechanism_stem: str = ""                                  # existing include set, e.g. "ethanol.fac"
    mech_incl_suffix: str = ".saphir.incl"
    use_anorg: bool = True
    use_j_clear_sky_incl: bool = True
    injections: List[InjectionEvent] = field(default_factory=list)
    chamber: ChamberEffects = field(default_factory=ChamberEffects)
    initial_all: str = "1."
    initials: List[Tuple[str, str]] = field(default_factory=list)   # (species, expression inside CONST)
    hmax: str = "30"
    outstep: str = "100000"
    extra_includes: List[str] = field(default_factory=list)
    overrides: List[str] = field(default_factory=list)              # expert statements, verbatim
    template: str = ""


# ─────────────────────────────────────────────────────────────────────────────
# Time grid and ENZ
# ─────────────────────────────────────────────────────────────────────────────

def _day_seconds(date: str) -> float:
    d = datetime.strptime(date, "%Y-%m-%d").replace(tzinfo=timezone.utc)
    return (d - ENZ_EPOCH).total_seconds()


def _hhmm_s(hhmm: str) -> int:
    h, m = hhmm.strip().split(":")
    return int(h) * 3600 + int(m) * 60


def time_grid(c: Conditions) -> np.ndarray:
    day = _day_seconds(c.date)
    return np.arange(day + _hhmm_s(c.start), day + _hhmm_s(c.end) + 1, c.step_s, dtype=float)


def pulse_column(grid: np.ndarray, c: Conditions, start: str, duration_min: float) -> np.ndarray:
    """Normalised injection pulse: value 1/n on the n rows from `start` (integral 1)."""
    col = np.zeros(len(grid))
    if not start:
        return col
    t0 = _day_seconds(c.date) + _hhmm_s(start)
    n = max(1, int(round(duration_min * 60.0 / c.step_s)))
    i0 = int(np.searchsorted(grid, t0 - 1e-6))
    i1 = min(i0 + n, len(grid))
    if i0 < len(grid):
        col[i0:i1] = 1.0 / n
    return col


def _solar_zenith(grid: np.ndarray) -> np.ndarray:
    try:
        from enz_builder import compute_solar_zenith
        return np.asarray(compute_solar_zenith(grid), dtype=float)
    except Exception:
        return np.full(len(grid), 45.0)


def _roof(grid: np.ndarray, c: Conditions) -> np.ndarray:
    hv = np.zeros(len(grid))
    day = _day_seconds(c.date)
    for o, cl in c.roof_open:
        t_o, t_c = day + _hhmm_s(o), day + _hhmm_s(cl)
        hv[(grid >= t_o) & (grid < t_c)] = 1.0
    return hv


def _archive_column(c: Conditions, enz_col: str, grid: np.ndarray, fill: float) -> Optional[np.ndarray]:
    """Read one standard ENZ column from the SAPHIR NetCDF archive (None if unavailable)."""
    try:
        from enz_builder import ARCHIVE_MAP, ARCHIVE_MAP_ALTS, read_nc_variable, interpolate_to_grid
    except Exception:
        return None
    if enz_col not in ARCHIVE_MAP:
        return None
    suffix, var, _f = ARCHIVE_MAP[enz_col]
    y, m, d = c.date.split("-")
    path = os.path.join(c.archive_root, y, m, d, f"{c.date}.{suffix}")
    if not os.path.isfile(path):
        return None
    for v in [var] + ARCHIVE_MAP_ALTS.get(var, []):
        try:
            t, vals = read_nc_variable(path, v)
            return np.asarray(interpolate_to_grid(t, vals, grid, fill), dtype=float)
        except Exception:
            continue
    return None


def enz_columns(spec: ExperimentSpec) -> Tuple[np.ndarray, Dict[str, np.ndarray], List[str]]:
    """Build the ENZ time grid and ordered column data for the spec.

    Returns (grid, columns, notes).  Column order follows the IDL files:
    sz, hv, dilute, t, H2O, injection switches, measured species, J values.
    """
    c = spec.conditions
    grid = time_grid(c)
    n = len(grid)
    notes: List[str] = []
    cols: Dict[str, np.ndarray] = {}
    cols["sz"] = _solar_zenith(grid)
    cols["hv"] = _roof(grid, c)

    def phys(name: str, const: float) -> np.ndarray:
        if c.source == "archive":
            arr = _archive_column(c, name, grid, const)
            if arr is not None:
                return arr
            notes.append(f"{name}: archive data not found, constant {const} used")
        return np.full(n, float(const))

    cols["dilute"] = phys("dilute", c.flow_m3h)
    cols["t"] = phys("t", c.temperature_K)
    cols["H2O"] = phys("H2O", c.h2o_percent)
    for inj in spec.injections:
        col = inj.enz_column()
        cols[col] = cols.get(col, np.zeros(n)) + pulse_column(grid, c, inj.start, inj.duration_min)
    if c.source == "archive":
        for name in ("NO2", "NO", "O3"):
            arr = _archive_column(c, name, grid, 0.0)
            if arr is not None:
                cols[name] = arr
    for sp in c.constrain_species:
        if sp not in cols:
            notes.append(f"constrained species {sp} has no ENZ column")
    if c.j_source == "measured":
        for jname, jcol, _rxn in J_INPUTS:
            enz_name = {"jno2": "jNO2", "jo1d": "jO1D", "jhono": "jHONO", "jh2o2": "jH2O2",
                        "jhcho_r": "jHCHO_r", "jhcho_m": "jHCHO_m"}[jname]
            arr = _archive_column(c, enz_name, grid, 0.0) if c.source == "archive" else None
            if arr is None:
                notes.append(f"{enz_name}: no measured J values, column filled with 0")
                arr = np.zeros(n)
            cols[enz_name] = arr
    return grid, cols, notes


def write_enz_file(path: str, grid: np.ndarray, cols: Dict[str, np.ndarray]) -> str:
    from enz_builder import write_enz
    write_enz(path, grid, cols)
    return path


# ─────────────────────────────────────────────────────────────────────────────
# EASY form
# ─────────────────────────────────────────────────────────────────────────────

def easy_form(spec: ExperimentSpec) -> EasyForm:
    c, ch = spec.conditions, spec.chamber
    f = default_form()
    f.title = spec.title or spec.name
    f.mech_stem = spec.mechanism_stem
    f.mech_incl_suffix = spec.mech_incl_suffix
    f.use_anorg = spec.use_anorg
    f.use_j_clear_sky = spec.use_j_clear_sky_incl
    f.includes_other = list(spec.extra_includes)

    f.walldep, f.walldep_value = ch.walldep, ch.walldep_value
    f.hono_source, f.hono_lhs, f.hono_rhs, f.hono_const, f.hono_rate = ch.hono, ch.hono_lhs, ch.hono_rhs, ch.hono_const, ch.hono_rate
    f.hcho_source, f.hcho_lhs, f.hcho_rhs, f.hcho_const, f.hcho_rate = ch.hcho, ch.hcho_lhs, ch.hcho_rhs, ch.hcho_const, ch.hcho_rate
    f.hcho_const_name = ch.hcho_const_name
    f.bgreac, f.y_name, f.y_product, f.y_ppb, f.y_const_expr, f.y_rate = ch.bgreac, ch.y_name, ch.y_product, ch.y_ppb, ch.y_expr, ch.y_rate
    f.dilution, f.dilute_expr, f.dil_value, f.dil_flow_const = ch.dilution, ch.dilute_expr, ch.dil_value, ch.dil_flow_const
    f.tracer_dilution = list(ch.tracer_dilution)
    f.dark_hono, f.dark_hono_rhs, f.dark_hono_rate = ch.dark_hono, ch.dark_hono_rhs, ch.dark_hono_rate
    f.light_source, f.light_source_rhs = ch.light_source, ch.light_source_rhs
    f.light_source_rate, f.light_source_initial = ch.light_source_rate, ch.light_source_initial

    # measured block
    if c.j_source == "measured":
        f.j_rows = [MeasuredRow(j, "enz", col) for j, col, _ in J_INPUTS]
        f.j_overrides = {j: rxn for j, _, rxn in J_INPUTS}
    else:
        f.j_rows = [MeasuredRow("jno2", "const", "", CLEAR_SKY_JNO2)]
        f.j_overrides = {"jno2": J_INPUTS[0][2]}
    t_row = (MeasuredRow("T", "enz", "t", comment="temperature in K") if c.temperature_mode == "enz"
             else MeasuredRow("T", "const", "t", f"{c.temperature_K:g}", "temperature in K"))
    f.setting_rows = [
        t_row,
        MeasuredRow("P", "const", "", f"{c.pressure_hPa:g}", "pressure in hPa"),
        MeasuredRow("SZA", "enz", "sz", comment="solar zenith angle"),
        MeasuredRow("SZ", "enz", "sz", comment="solar zenith angle"),
        MeasuredRow("HV", "enz", "hv", comment="switch for light"),
        MeasuredRow("FLOW", "enz", "dilute", comment="flow into the chamber in m^3/h"),
        MeasuredRow("H2O", "enz", "H2O", comment="mixing ratio of water"),
    ]
    f.enz_file = f"{spec.name}.enz"
    f.hmax, f.outstep = spec.hmax, spec.outstep
    f.initials = [("*", f"CONST({spec.initial_all})", "")] + [(sp, f"CONST({ex})", "") for sp, ex in spec.initials]

    f.injections, f.chain_reactions = [], []
    for inj in spec.injections:
        if not inj.products.strip() and not inj.extra_channels:
            continue                       # ENZ-only event (recorded, not modelled)
        sw, col = inj.easy_switch(), inj.enz_column()
        f.injections.append(InjectionRow(sw, col, inj.products, str(inj.total_ppb), inj.comment))
        for prod, ppb in inj.extra_channels:
            f.injections.append(InjectionRow(sw, col, prod, str(ppb)))
        f.chain_reactions.extend(inj.chain)
    f.extra_statements = []
    for sp in c.constrain_species:
        f.extra_statements += [f"{sp}_in=input({sp})", f"{sp}=CONST({sp}_in*1e-9*M)"]
    f.extra_statements += list(spec.overrides)
    return f


def render_easy(spec: ExperimentSpec) -> str:
    return render_form(easy_form(spec))


# ─────────────────────────────────────────────────────────────────────────────
# Backend: mechanism subset + file writing
# ─────────────────────────────────────────────────────────────────────────────

FULL_MCM_FAC = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            "MCM_database", "MCM_3.3.1", "mcm_3-3-1_facsimile_complete.fac")


def ensure_mechanism(spec: ExperimentSpec, out_dir: str, full_fac: str = FULL_MCM_FAC) -> dict:
    """Extract the MCM subset for spec.precursors into <out_dir>/chemistry and point
    the spec at it.  No-op when spec.mechanism_stem is already set."""
    if spec.mechanism_stem or not spec.precursors:
        return {}
    import io, contextlib
    from mcm_converter import convert_mcm_file
    chem = os.path.join(out_dir, "chemistry")
    with contextlib.redirect_stdout(io.StringIO()):
        res = convert_mcm_file(full_fac, chem, saphir=True,
                               measured=(spec.conditions.j_source == "measured"),
                               precursors=list(spec.precursors))
    spec.mechanism_stem = os.path.basename(str(res["mech"]))[:-len(".mech")]
    spec.mech_incl_suffix = ".incl"
    spec.use_anorg = False               # inorganic core is inlined in init mode
    return res


def write_experiment(spec: ExperimentSpec, out_dir: str) -> dict:
    """Write <name>.enz and <name>.easy (and the mechanism subset if needed)."""
    os.makedirs(out_dir, exist_ok=True)
    res = {"mechanism": ensure_mechanism(spec, out_dir)}
    grid, cols, notes = enz_columns(spec)
    res["enz"] = write_enz_file(os.path.join(out_dir, f"{spec.name}.enz"), grid, cols)
    text = render_easy(spec)
    easy = os.path.join(out_dir, f"{spec.name}.easy")
    with open(easy, "w", encoding="latin-1", errors="replace", newline="\n") as fh:
        fh.write(text)
    res["easy"] = easy
    res["notes"] = notes
    return res


# ─────────────────────────────────────────────────────────────────────────────
# Serialisation, templates, validation
# ─────────────────────────────────────────────────────────────────────────────

def spec_to_dict(spec: ExperimentSpec) -> dict:
    return asdict(spec)


def spec_from_dict(d: dict) -> ExperimentSpec:
    d = dict(d)
    cond = Conditions(**{k: v for k, v in d.pop("conditions", {}).items() if k in Conditions.__dataclass_fields__})
    cond.roof_open = [tuple(x) for x in cond.roof_open]
    ch = ChamberEffects(**{k: v for k, v in d.pop("chamber", {}).items() if k in ChamberEffects.__dataclass_fields__})
    inj = []
    for e in d.pop("injections", []):
        e = {k: v for k, v in e.items() if k in InjectionEvent.__dataclass_fields__}
        e["extra_channels"] = [tuple(x) for x in e.get("extra_channels", [])]
        e["chain"] = [tuple(x) for x in e.get("chain", [])]
        inj.append(InjectionEvent(**e))
    spec = ExperimentSpec(**{k: v for k, v in d.items() if k in ExperimentSpec.__dataclass_fields__})
    spec.conditions, spec.chamber, spec.injections = cond, ch, inj
    spec.initials = [tuple(x) for x in spec.initials]
    return spec


def save_spec(spec: ExperimentSpec, path: str) -> str:
    import json
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(spec_to_dict(spec), fh, indent=1, ensure_ascii=False)
    return path


def load_spec(path: str) -> ExperimentSpec:
    import json
    with open(path, encoding="utf-8") as fh:
        return spec_from_dict(json.load(fh))


TEMPLATE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "templates")


def list_templates() -> Dict[str, str]:
    """name → path of the experiment templates (templates/*.experiment.json)."""
    out = {}
    if os.path.isdir(TEMPLATE_DIR):
        for f in sorted(os.listdir(TEMPLATE_DIR)):
            if f.endswith(".experiment.json"):
                out[f[:-len(".experiment.json")]] = os.path.join(TEMPLATE_DIR, f)
    return out


def load_template(path: str) -> Tuple[str, str, ExperimentSpec]:
    """Return (title, description, spec) of a templates/*.experiment.json file.
    The template holds a partial spec; missing fields take the dataclass defaults."""
    import json
    with open(path, encoding="utf-8") as fh:
        d = json.load(fh)
    spec = spec_from_dict(d.get("spec", {}))
    spec.template = os.path.basename(path)[:-len(".experiment.json")]
    return d.get("title", spec.template), d.get("description", ""), spec


def mechanism_species(include_dir: str, stem: str) -> set:
    """Species names appearing in the mechanism files of <stem> (its .mech and
    .dil.mech) plus the inorganic anorg*.mech files and the standard chamber
    species; uppercase.  Used for validation warnings only."""
    names = {"XNOY", "NA", "SHONO", "SHCHO", "QSHCHO", "Y", "QY", "SNO2", "H2O", "HV", "DIL", "O2", "N2", "M"}
    if not include_dir or not stem:
        return names
    wanted = {f"{stem}.mech".lower(), f"{stem}.dil.mech".lower()}
    for root, _d, files in os.walk(include_dir):
        for f in files:
            fl = f.lower()
            if fl in wanted or (fl.startswith("anorg") and fl.endswith(".mech")):
                for line in open(os.path.join(root, f), encoding="latin-1"):
                    code = line.split(";", 1)[0]
                    if "-->" in code:
                        for tok in re.findall(r"[A-Za-z][A-Za-z0-9_]*", code.replace("-->", "+")):
                            names.add(tok.upper())
    return names


def validate(spec: ExperimentSpec, known_species: Optional[set] = None) -> List[Tuple[str, str]]:
    """Return (level, message) items; level is 'error' or 'warning'."""
    msgs: List[Tuple[str, str]] = []
    c = spec.conditions
    try:
        grid = time_grid(c)
        if len(grid) < 2:
            msgs.append(("error", "Time window is empty — check start and end time."))
    except Exception as e:
        msgs.append(("error", f"Date or time is not valid: {e}"))
    if not spec.mechanism_stem and not spec.precursors:
        msgs.append(("error", "No chemistry: choose compounds or an existing mechanism."))
    if c.j_source == "measured" and c.source != "archive":
        msgs.append(("warning", "Measured J values selected without archive data: the J columns will be 0 "
                                "(no photolysis). Use 'clear sky' or provide archive data."))
    if not c.roof_open:
        msgs.append(("warning", "Roof never opens (HV = 0 all day): dark experiment."))
    seen_cols = set()
    for inj in spec.injections:
        col = inj.enz_column()
        if inj.start:
            try:
                _hhmm_s(inj.start)
            except Exception:
                msgs.append(("error", f"Injection '{inj.name}': start time '{inj.start}' is not HH:MM."))
        if inj.products:
            try:
                float(eval(str(inj.total_ppb), {"__builtins__": {}}, {}))
            except Exception:
                msgs.append(("warning", f"Injection '{inj.name}': amount '{inj.total_ppb}' is an expression, "
                                        "it will be written verbatim."))
            if known_species:
                for term in inj.products.replace(" ", "").split("+"):
                    sp = re.sub(r"^[0-9.eE+\-]*\*", "", term).upper()
                    if sp and sp not in known_species:
                        msgs.append(("warning", f"Injection '{inj.name}': product '{sp}' is not a species of the mechanism."))
        elif not inj.start:
            msgs.append(("warning", f"Injection '{inj.name}': neither modelled nor timed — it has no effect."))
        seen_cols.add(col)
    if not spec.injections:
        msgs.append(("warning", "No injections defined."))
    return msgs


# ─────────────────────────────────────────────────────────────────────────────
# Golden specs: scene_05 (ethanol) and scene_09 (CHANEL)
# ─────────────────────────────────────────────────────────────────────────────

def spec_scene_05() -> ExperimentSpec:
    cond = Conditions(date="2024-07-09", start="04:00", end="21:59", temperature_K=300,
                      temperature_mode="const", roof_open=[("12:20", "17:20")])
    ch = ChamberEffects(hono_lhs="SHONO", hono_rhs="HONO", hono_const="HV_in*2.6*1.1*0.7",
                        hcho_lhs="SHCHO", hcho_rhs="HCHO", hcho_const_name="SHCHO", hcho_const="HV_in*1",
                        y_name="Y", y_product="HO2", y_ppb="230", tracer_dilution=[],
                        dark_hono_rhs="HONO", dark_hono_rate="1e-012",
                        light_source=True, light_source_rhs="NO2")
    inj = [
        InjectionEvent("NO", "NO", "8.1", "11:49", 5, extra_channels=[("NO2", "0.7")]),
        InjectionEvent("H2O", "HONO", "0.", "09:13", 66,
                       extra_channels=[("CO", "0.5"), ("NO2", "0.12"), ("NO", "0.0")]),
        InjectionEvent("O3", "O3", "113*0"),
        InjectionEvent("VCP2", "0.9*VOCin1+0.1*VOCin2", "575*180/175", "11:10", 2,
                       column="qVCP2", switch="qVOC",
                       chain=[("VOCin1", "VOCin_all", "1./60"), ("VOCin2", "VOCin_all", "1./1800"),
                              ("VOCin_all", "C2H5OH", "1./60")], comment="for ethanol"),
        # recorded in the chamber log but not modelled in this scene (ENZ columns only)
        InjectionEvent("gasoline", "", start="07:50", duration_min=2, column="qGASOLINE"),
        InjectionEvent("biogenic", "", start="07:52", duration_min=2, column="qBIOGENIC"),
        InjectionEvent("cooking", "", start="07:55", duration_min=2, column="qCOOKING"),
        InjectionEvent("diesel", "", start="08:12", duration_min=2, column="qDIESEL"),
        InjectionEvent("seed", "", start="10:20", duration_min=10, column="Q_seed"),
        InjectionEvent("isobutane", "", start="11:20", duration_min=2, column="q_isobutane"),
    ]
    return ExperimentSpec(name="scene_05", title="scene_05 rebuilt from ExperimentSpec", conditions=cond,
                          mechanism_stem="ethanol.fac", injections=inj, chamber=ch, initial_all="1")


def spec_scene_09() -> ExperimentSpec:
    cond = Conditions(date="2024-08-05", start="04:00", end="21:59", roof_open=[("09:33", "18:37")])
    ch = ChamberEffects(light_source=True)
    mix = {
        "VCP2": ("0.016042*LIMONENE+0.107333*PROPGLY+0.002273*APINENE+0.002273*BPINENE+0.469732*C2H5OH"
                 "+0.027435*ETHGLY+0.067800*NC12H26+0.209843*IPROPOL+0.026511*TOLUENE+0.035807*MXYL"
                 "+0.003172*BOX2PROL+0.008793*DIIPRETHER+0.004353*BUOX2ETOH", "117.813"),
        "gasoline": ("0.092669*CH3CHO+0.060555*BENZENE+0.242439*C2H5OH+0.146659*NC12H26+0.089599*TOLUENE"
                     "+0.044085*PXYL+0.032659*OXYL+0.067681*MXYL+0.021637*TM124B+0.009209*CPENT2ENE"
                     "+0.192807*IC5H12", "14.5878"),
        "biogenic": (" 0.201060*C5H8+0.205849*APINENE+0.065903*BPINENE+0.125067*LIMONENE+0.125067*X_limonene", "24.2329"),
        "cooking": ("0.268706*CH3CHO+0.338202*C2H5CHO+0.122994*C4H9CHO+0.257877*C4ALDB+0.012220*LIMONENE"
                    "+0.012220*X_limonene", "31.9793"),
    }
    inj = [
        InjectionEvent("H2O", "HONO+XNOy", "0.002", "05:26", 45,
                       extra_channels=[("CO", "0.5"), ("NO2+XNOy", "0.00"), ("NO+XNOy", "0.000")]),
        InjectionEvent("O3", "O3", "52", "09:36", 20),
        InjectionEvent("VCP2", mix["VCP2"][0], mix["VCP2"][1], "08:02", 2, column="qVCP2", switch="qVCP2"),
        InjectionEvent("gasoline", mix["gasoline"][0], mix["gasoline"][1], "07:55", 2, column="qGASOLINE", switch="Q_gasoline"),
        InjectionEvent("biogenic", mix["biogenic"][0], mix["biogenic"][1], "07:56", 2, column="qBIOGENIC", switch="Q_biogenic"),
        InjectionEvent("cooking", mix["cooking"][0], mix["cooking"][1], "07:59", 2, column="qCOOKING", switch="Q_cooking"),
        InjectionEvent("ethanol", "C2H5OH", "16.043", "07:58", 5, comment="11.5µL  1.29s-1"),
        InjectionEvent("ethylene", "C2H4", "7.3", "07:10", 5),
        InjectionEvent("propene", "C3H6", "3.7", "07:23", 5),
        InjectionEvent("13butadiene", "C4H6", "0.54", "07:49", 5),
        InjectionEvent("isobutene", "MEPROPENE", "1", "07:17", 5),
        InjectionEvent("isobutane", "IC4H10", "15.4", "07:27", 5),
        InjectionEvent("trans2butene", "TBUT2ENE", "0.22", "07:45", 5),
        InjectionEvent("seed", "", start="06:28", duration_min=10, column="Q_seed"),
    ]
    return ExperimentSpec(name="scene_09", title="scene_09 rebuilt from ExperimentSpec", conditions=cond,
                          mechanism_stem="chanel.all_mixes.fac", injections=inj, chamber=ch,
                          initials=[("CO", "M*1e-9*14.0"), ("CO2", "M*1e-6*0.25*0"), ("NO2", "M*1e-9*0.00")])


# ─────────────────────────────────────────────────────────────────────────────
# Self-test
# ─────────────────────────────────────────────────────────────────────────────

def _read_enz(path: str) -> Tuple[List[str], np.ndarray]:
    cols, rows, data = [], [], False
    for line in open(path, encoding="latin-1"):
        if data:
            s = line.replace(",", " ").split()
            if s:
                rows.append([float(v) for v in s])
        elif line.startswith("COLUMN"):
            cols.append(line.split("=", 1)[1].split("/")[0].strip())
        elif "&&&" in line:
            data = True
    return cols, np.array(rows)


def _selftest() -> int:
    import tempfile
    here = os.path.dirname(os.path.abspath(__file__))
    tmp = tempfile.mkdtemp(prefix="experiment_spec_")
    bad = 0
    for spec in (spec_scene_05(), spec_scene_09()):
        sc = spec.name
        res = write_experiment(spec, os.path.join(tmp, sc))
        orig_easy = os.path.join(here, "MCM_tested_scenarios", "IDL Input", sc, sc + ".easy")
        ok, diffs = statements_equivalent(orig_easy, res["easy"], ignore_duplicate_reactions=True)
        # ENZ: injection switch columns and roof must match the archived ENZ exactly
        orig_enz = os.path.join(here, "MCM_tested_scenarios", "JUICE archive", "JUICE", sc, sc + ".enz")
        enz_msgs = []
        if os.path.isfile(orig_enz):
            oc, oa = _read_enz(orig_enz)
            nc, na = _read_enz(res["enz"])
            if len(oa) != len(na) or not np.allclose(oa[:, 0], na[:, 0]):
                enz_msgs.append(f"time grid differs ({len(oa)} vs {len(na)} rows)")
            else:
                for col in [c for c in oc if c.upper().startswith("Q") or c == "hv"]:
                    if col not in nc:
                        if np.any(oa[:, oc.index(col)] != 0):
                            enz_msgs.append(f"column {col} missing in generated ENZ")
                        continue
                    if not np.allclose(oa[:, oc.index(col)], na[:, nc.index(col)], atol=1e-9):
                        enz_msgs.append(f"column {col} differs")
        else:
            enz_msgs.append("no archived ENZ to compare")
        good = ok and not enz_msgs
        bad += not good
        print(f"{'OK  ' if good else 'FAIL'} {sc}: easy-equivalent={ok} ({len(diffs)} diffs)  enz: "
              f"{'identical switch/roof columns' if not enz_msgs else '; '.join(enz_msgs)}")
        for d in diffs[:8]:
            print("       " + d[:150])
        if res["notes"]:
            print("       notes: " + "; ".join(res["notes"][:4]))
    print(f"{2 - bad}/2 golden specs pass  (files in {tmp})")
    return 1 if bad else 0


if __name__ == "__main__":
    if "--selftest" in sys.argv:
        sys.exit(_selftest())
    print(__doc__)
