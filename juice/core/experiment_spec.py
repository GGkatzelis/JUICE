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

Golden tests (`python -m juice.core.experiment_spec --selftest`) rebuild scene_05 and
scene_09 from specs and require statement-equivalent EASY files and identical
injection / roof columns in the ENZ.
"""

from __future__ import annotations
from juice._root import ROOT as _JUICE_ROOT

import os
import re
import sys
import math
from dataclasses import dataclass, field, asdict
from datetime import datetime, timedelta, timezone
from typing import Dict, List, Optional, Tuple

from juice.chemistry import mechanism_registry as mreg

import numpy as np

from juice.core.easy_builder import (EasyForm, MeasuredRow, InjectionRow, default_form, render_form,
                          J_INPUTS, HONO_RATE, HCHO_RATE, Y_RATE, DILUTE_DEF)
from juice.core.easy_scene import statements_equivalent

ENZ_EPOCH = datetime(2000, 1, 1, 0, 0, 0, tzinfo=timezone.utc)
CLEAR_SKY_JNO2 = "PHTFKT_CS(1.165D-02,0.244,-0.267)"

# default ENZ column per injected compound (IDL convention); others → Q_<name>
STANDARD_COLUMNS = {"NO": "QNO", "NO2": "QNO2", "O3": "QO3", "H2O": "QH2O", "CO2": "QCO2",
                    "CL2": "QCl2", "SEED": "Q_seed", "CO": "QCO"}
# humidification carries chamber impurities (ppb per humidification, scene_05/09 defaults)
H2O_IMPURITY_DEFAULTS = [("HONO+XNOy", "0.002"), ("CO", "0.5"), ("NO2+XNOy", "0.00"), ("NO+XNOy", "0.000")]

# Injection actions: what the operator did at the chamber.  Each action carries the
# defaults a new row gets in the wizard (name, product, total ppb, duration [min],
# ENZ column, EASY switch, comment, extra channels).  The action only steers defaults
# and labelling; the ENZ/EASY writers work from the InjectionEvent fields.
ACTIONS: Dict[str, dict] = {
    "humidification": dict(label="Humidification (water vapour + impurities)", name="H2O",
                           products=H2O_IMPURITY_DEFAULTS[0][0], total_ppb=H2O_IMPURITY_DEFAULTS[0][1],
                           duration_min=45, column="QH2O", switch="qH2O",
                           comment="humidification: HONO impurity (SAPHIR default)",
                           extra_channels=list(H2O_IMPURITY_DEFAULTS[1:])),
    "voc":            dict(label="VOC injection", name="VOC", products="", total_ppb="10", duration_min=2,
                           column="", switch="", comment="set the product to the MCM name and the total ppb"),
    "o3":             dict(label="Ozone injection", name="O3", products="O3", total_ppb="50", duration_min=20,
                           column="QO3", switch="qO3", comment=""),
    "no":             dict(label="NO injection", name="NO", products="NO", total_ppb="5", duration_min=1,
                           column="QNO", switch="qNO", comment=""),
    "no2":            dict(label="NO2 injection", name="NO2", products="NO2", total_ppb="5", duration_min=1,
                           column="QNO2", switch="qNO2", comment=""),
    "co":             dict(label="CO injection", name="CO", products="CO", total_ppb="100", duration_min=1,
                           column="QCO", switch="qCO", comment=""),
    "co2_tracer":     dict(label="CO2 tracer (dilution)", name="CO2", products="CO2", total_ppb="0", duration_min=1,
                           column="QCO2", switch="qCO2", comment="dilution tracer"),
    "seed":           dict(label="Seed aerosol (recorded only)", name="seed", products="", total_ppb="",
                           duration_min=10, column="Q_seed", switch="", comment="recorded only"),
    "mixture":        dict(label="Emission mixture (recorded only / surrogates)", name="mixture", products="",
                           total_ppb="", duration_min=2, column="", switch="", comment="recorded only"),
    "other":          dict(label="Other", name="", products="", total_ppb="0", duration_min=1,
                           column="", switch="", comment=""),
}
ACTION_ORDER = list(ACTIONS)

# wiki action-log codes (saphirwiki "HH:MM <code> <person> <type>: text") → action
WIKI_CODE_ACTION = {"H": "humidification", "h": "humidification", "X": "voc", "O": "o3", "M": "no", "D": "no2",
                    "c": "co", "C": "co2_tracer", "S": "seed", "s": "seed"}


def infer_action(name: str, products: str = "", column: str = "") -> str:
    """Best-effort action for specs written before the action field existed."""
    n = (name or "").strip().upper()
    if n == "H2O" or (column or "").upper() == "QH2O":
        return "humidification"
    if n in ("O3",):
        return "o3"
    if n == "NO":
        return "no"
    if n == "NO2":
        return "no2"
    if n == "CO":
        return "co"
    if n == "CO2":
        return "co2_tracer"
    if n == "SEED":
        return "seed"
    if n in ("GASOLINE", "DIESEL", "COOKING", "BIOGENIC", "VCP1", "VCP2", "MIXTURE"):
        return "mixture"
    return "voc" if (products or "").strip() else "other"


# ─────────────────────────────────────────────────────────────────────────────
# Specification
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class Conditions:
    date: str = "2024-07-09"              # YYYY-MM-DD (UTC)
    start: str = "04:00"                  # HH:MM
    end: str = "21:59"
    step_s: int = 60
    source: str = "manual"                # manual | archive | file
    archive_root: str = r"Z:\IEK8-SAPHIR\data"
    input_file: str = ""                  # source 'file': CSV/Excel (or an ENZ) with time + T, p, H2O/RH, flow, roof, sz, J values, NO/NO2/O3 (measured_inputs.py)
    temperature_K: float = 298.0
    temperature_mode: str = "enz"         # enz: T = CONST(T_in) | const: T = CONST(<temperature_K>)
    pressure_hPa: float = 1006.5
    h2o_percent: float = 1.0              # ENZ H2O column (vol-%); EASY: H2O = CONST(1e-2*M*H2O_in)
    flow_m3h: float = 6.0                 # ENZ dilute column (m^3/h); EASY: FLOW = INPUT(dilute)
    roof_open: List[Tuple[str, str]] = field(default_factory=list)   # [(open HH:MM, close HH:MM)]
    j_source: str = "measured"            # measured (ENZ j columns) | clear_sky (jno2 = PHTFKT_CS, rest from mechanism)
    constrain_species: List[str] = field(default_factory=list)       # e.g. ["NO2", "O3"] → held to ENZ columns
    chamber_name: str = "SAPHIR"          # injection_calc.CHAMBERS profile (volume, MFC/syringe standards)

    # In a lamp chamber "roof open" means "lights on": the same field under the name that fits.
    @property
    def lights_on(self) -> List[Tuple[str, str]]:
        return self.roof_open

    @lights_on.setter
    def lights_on(self, value: List[Tuple[str, str]]) -> None:
        self.roof_open = value


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
    action: str = ""                      # key of ACTIONS (humidification, voc, o3, ...); "" = infer from name
    delivery: Dict[str, str] = field(default_factory=dict)   # how the amount was derived (injection_calc.resolve_delivery)

    def action_key(self) -> str:
        return self.action if self.action in ACTIONS else infer_action(self.name, self.products, self.column)

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
    hcho_const_name: str = "qSHCHO"     # must match the reactant: scene_09 writes SHCHO here, which leaves qSHCHO an undefined species → source inert (FACSIMILE too)
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
    wall_loss: bool = False                                   # generic first-order wall losses of the chamber profile
    wall_loss_s1: Dict[str, float] = field(default_factory=dict)   # species → k [1/s]; empty for SAPHIR (saphir.walldep.incl instead)


@dataclass
class ExperimentSpec:
    name: str = "experiment"
    title: str = ""
    conditions: Conditions = field(default_factory=Conditions)
    precursors: List[str] = field(default_factory=list)      # MCM names → mechanism subset (backend)
    base_mechanism: str = mreg.DEFAULT_BASE                   # mechanism_registry base id (MCM download, unmodified)
    addons: List[str] = field(default_factory=list)          # mechanism_registry add-on ids, loaded after the base
    mechanism_stem: str = ""                                  # existing include set, e.g. "ethanol.fac"
    inorganic_only: bool = False                              # no VOC: inorganic core + chamber sources only (CO/NOx/O3/HOx tests)
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
    # (reaction "A + B --> C", new rate expression or "0" = remove, comment) → written exactly as EASY
    # does it: k[A + B --> C] = CONST(expr) after the mechanism includes (last definition wins)
    reaction_adjustments: List[Tuple[str, str, str]] = field(default_factory=list)
    import_notes: List[str] = field(default_factory=list)       # what the wiki / mixture import could not resolve (shown under 'Needs attention')
    template: str = ""
    include_dirs: List[str] = field(default_factory=list)       # more folders searched for #include files (an imported FACSIMILE experiment's own)
    easy_version: str = ""                                      # EASY that made the FACSIMILE run this came from ('5.96b'); '' = JUICE's rules (elog_generator.function_scope_for)
    reference_cs: str = ""                                      # FACSIMILE results of the same experiment: 'Compare with FACSIMILE' on the Run step
    gecko_core: str = "gecko"           # GECKO-A base: inorganic core 'gecko' (as generated) or 'mcm' (MCM v3.3.1 + aliases, measured J include)
    gecko_partitioning: bool = True     # GECKO-A base: load the gas <-> particle <-> wall partitioning include (Nannoolal C*)


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


def _solar_zenith(grid: np.ndarray, c: Optional[Conditions] = None) -> np.ndarray:
    """Solar zenith angle at the chamber's site (profile lat/lon; SAPHIR's when the profile has none)."""
    try:
        from juice.core import chambers
        from juice.core.enz_builder import compute_solar_zenith
        lat, lon = chambers.chamber(c.chamber_name if c else "").lat_lon
        return np.asarray(compute_solar_zenith(grid, lat, lon), dtype=float)
    except Exception:
        return np.full(len(grid), 45.0)


def chamber_defaults(profile_id: str, base: Optional[ChamberEffects] = None) -> ChamberEffects:
    """ChamberEffects for a chamber profile: its default switches and its DILUTE expression; every
    formula and strength stays as in `base` (or the SAPHIR defaults of the dataclass)."""
    from juice.core import chambers
    p = chambers.chamber(profile_id)
    ch = ChamberEffects(**{k: getattr(base, k) for k in ChamberEffects.__dataclass_fields__}) if base else ChamberEffects()
    for key in chambers.EFFECT_KEYS:
        setattr(ch, key, p.effect(key))
    ch.dilute_expr = chambers.render_dilute(p)
    ch.wall_loss_s1 = {str(k): float(v) for k, v in (p.wall_loss_s1 or {}).items()}
    ch.wall_loss = False                         # opt-in: irreversible first-order losses of inorganics only; organic
    return ch                                    # vapour-wall partitioning (reversible, volatility-dependent) comes with the SOA feature


def _profile(c: Conditions):
    from juice.core import chambers
    return chambers.chamber(c.chamber_name)


def apply_chamber(spec: ExperimentSpec, profile_id: str) -> ExperimentSpec:
    """Switch an experiment to another chamber: name, archive root, default chamber effects and
    DILUTE.  Called when the user picks a chamber, never when a saved experiment is loaded (what was
    saved is what runs).  For SAPHIR this reproduces the dataclass defaults exactly."""
    from juice.core import chambers
    p = chambers.chamber(profile_id)
    c = spec.conditions
    c.chamber_name = p.id
    if p.has_archive:
        c.archive_root = p.archive_root
    elif c.source == "archive":
        c.source = ""                               # nothing to read from: the user has to choose file or manual
    spec.chamber = chamber_defaults(p.id, spec.chamber)
    return spec


def _roof(grid: np.ndarray, c: Conditions) -> np.ndarray:
    hv = np.zeros(len(grid))
    day = _day_seconds(c.date)
    for o, cl in c.roof_open:
        t_o, t_c = day + _hhmm_s(o), day + _hhmm_s(cl)
        hv[(grid >= t_o) & (grid < t_c)] = 1.0
    return hv


# Older SAPHIR archive naming (e.g. 2017): same quantities, other file and variable names.
# Tried in order after the current name (enz_builder.ARCHIVE_MAP); units are the ENZ units already.
ARCHIVE_OLD_NAMES: Dict[str, List[Tuple[str, str]]] = {
    "t":       [("saphir_collected.selected_param.60sec.nc", "tk_m")],
    "dilute":  [("saphir_collected.selected_param.60sec.nc", "FIC401501601t")],
    "H2O":     [("NOxContainer_PICARRO.H2O.nc", "H2O")],
    "NO":      [("NOxContainer.no.nc", "NO")],
    "NO2":     [("NOxContainer.no2.nc", "NO2")],
    "O3":      [("NOxContainer.o3_uv.nc", "O3")],
    "jNO2":    [("saphir_photolysis.jno2.nc", "jNO2")],
    "jO1D":    [("saphir_photolysis.jo1d.nc", "jO1D")],
    "jHONO":   [("saphir_photolysis.jhono.nc", "jHONO")],
    "jH2O2":   [("saphir_photolysis.jh2o2.nc", "jH2O2")],
    "jHCHO_r": [("saphir_photolysis.jhchor.nc", "jHCHOR")],
    "jHCHO_m": [("saphir_photolysis.jhchom.nc", "jHCHOM")],
}


# The archive is on a network drive: a file-existence check costs ~50 ms and a NetCDF read 0.1–0.2 s.
# The live review refreshes after every edit, so every day folder is listed once and every
# (file, variable) read once per session; clear_archive_cache() forgets both (e.g. after new data arrived).
_DAY_LISTING: Dict[str, set] = {}
_NC_READS: Dict[Tuple[str, str], Optional[Tuple[np.ndarray, np.ndarray]]] = {}


def clear_archive_cache() -> None:
    _DAY_LISTING.clear(); _NC_READS.clear()


def _archive_has(path: str) -> bool:
    folder, name = os.path.split(path)
    if folder not in _DAY_LISTING:
        try:                                   # lower case: Windows file names ignore case (SAPHIR.SR.jNO2.nc vs .jno2.nc)
            _DAY_LISTING[folder] = {f.lower() for f in os.listdir(folder)}
        except OSError:
            _DAY_LISTING[folder] = set()
    return name.lower() in _DAY_LISTING[folder]


def _archive_read(path: str, var: str) -> Optional[Tuple[np.ndarray, np.ndarray]]:
    """(times in ENZ seconds, values) of one variable, read once per session; None if it is not there."""
    key = (path, var)
    if key not in _NC_READS:
        try:
            from juice.core.enz_builder import read_nc_variable
            t, vals = read_nc_variable(path, var)
            _NC_READS[key] = (np.asarray(t, dtype=float), np.asarray(vals, dtype=float))
        except Exception:
            _NC_READS[key] = None
    return _NC_READS[key]


def archive_sources(c: Conditions, enz_col: str) -> List[Tuple[str, List[str]]]:
    """(path, variable candidates) for one ENZ column: the current archive name first, then the older ones."""
    try:
        from juice.core.enz_builder import ARCHIVE_MAP, ARCHIVE_MAP_ALTS
    except Exception:
        return []
    y, m, d = c.date.split("-")
    day = os.path.join(c.archive_root, y, m, d)
    out = []
    if enz_col in ARCHIVE_MAP:
        suffix, var, _f = ARCHIVE_MAP[enz_col]
        out.append((os.path.join(day, f"{c.date}.{suffix}"), [var] + ARCHIVE_MAP_ALTS.get(var, [])))
    for suffix, var in ARCHIVE_OLD_NAMES.get(enz_col, []):
        out.append((os.path.join(day, f"{c.date}.{suffix}"), [var]))
    return out


def _archive_column(c: Conditions, enz_col: str, grid: np.ndarray, fill: float) -> Optional[np.ndarray]:
    """Read one standard ENZ column from the SAPHIR NetCDF archive (None if unavailable).
    Current file names first, then the older naming of the archive (ARCHIVE_OLD_NAMES)."""
    try:
        from juice.core.enz_builder import interpolate_to_grid
    except Exception:
        return None
    for path, candidates in archive_sources(c, enz_col):
        if not _archive_has(path):
            continue
        for v in candidates:
            got = _archive_read(path, v)
            if got is not None:
                return np.asarray(interpolate_to_grid(got[0], got[1], grid, fill), dtype=float)
    if enz_col == "H2O":
        return _archive_h2o_from_rh(c, grid, fill)
    return None


def _archive_h2o_from_rh(c: Conditions, grid: np.ndarray, fill: float) -> Optional[np.ndarray]:
    """Water vapour [vol-%] from the chamber's relative humidity, temperature and pressure in the
    collected file — for days without a Picarro water file (e.g. 2026-05-27).  Magnus formula, as in
    measured_inputs; negative RH readings of a dry chamber are taken as 0."""
    try:
        from juice.core.enz_builder import interpolate_to_grid
        from juice.core.measured_inputs import _es_hpa
    except Exception:
        return None
    y, m, d = c.date.split("-")
    path = os.path.join(c.archive_root, y, m, d, f"{c.date}.SAPHIR.collected.all_param.nc")
    if not _archive_has(path):
        return None

    def col(names, default=None):
        for v in names:
            got = _archive_read(path, v)
            if got is not None:
                return np.asarray(interpolate_to_grid(got[0], got[1], grid, np.nan), dtype=float)
        return default
    rh = col(["RH_SAPHIR_VAISALA2", "MT401B_RH_SAPHIR"])
    T = col(["TK_M_SAPHIR"])
    if rh is None or T is None:
        return None
    p = col(["p_SAPHIR"], np.full(len(grid), c.pressure_hPa))
    h2o = np.clip(rh, 0.0, None) / 100.0 * _es_hpa(T) / p * 100.0
    return np.where(np.isfinite(h2o), h2o, fill)


def _file_column(c: Conditions, enz_col: str, grid: np.ndarray, fill: float) -> Optional[np.ndarray]:
    """One standard ENZ column from the user's measured-inputs file (None if the file lacks it)."""
    try:
        from juice.core import measured_inputs as mi
        inp = mi.cached_inputs(c.input_file, c.date, c.pressure_hPa)
        if inp is None:
            return None
        return mi.column_on_grid(inp, enz_col, grid, fill)
    except Exception:
        return None


def _measured_column(c: Conditions, enz_col: str, grid: np.ndarray, fill: float) -> Optional[np.ndarray]:
    if c.source == "archive":
        return _archive_column(c, enz_col, grid, fill)
    if c.source == "file":
        return _file_column(c, enz_col, grid, fill)
    return None


# EASY photolysis variable → ENZ column written for it
J_ENZ_NAME = {"jno2": "jNO2", "jo1d": "jO1D", "jhono": "jHONO", "jh2o2": "jH2O2", "jhcho_r": "jHCHO_r", "jhcho_m": "jHCHO_m"}


def measured_j_rows(c: Conditions, grid: Optional[np.ndarray] = None) -> Tuple[List[Tuple[str, str, str]], List[str]]:
    """(the J_INPUTS rows the EASY file takes from measurements, notes).

    SAPHIR archive: always all six, as in the IDL/FACSIMILE workflow (a day without a J file gives
    a 0 column and a note — unchanged, so the golden scenarios stay byte-identical).
    User file: only the J values the file has.  A J value the file lacks is NOT overridden, so the
    mechanism's own value stays in force — the MCM clear-sky shape scaled by the measured jNO2
    (PHTFKT) — instead of being set to 0, which would switch that photolysis off.  Without jNO2
    nothing can be scaled: the run falls back to the clear-sky model.
    Lamp chamber: jNO2 always (the measured column, else jNO2 at full power × lights); the other
    five only when a file has them; everything else comes from the lamp J table (lamp_photolysis)."""
    p = _profile(c)
    if p.light_type == "lamps":
        grid = time_grid(c) if grid is None else grid
        have = [J_INPUTS[0]]
        jfull = float(p.light.get("jno2_full_s1", 0) or 0)
        notes: List[str] = []
        if c.source == "file" and _file_column(c, "jNO2", grid, 0.0) is not None:
            pass
        else:
            notes.append(f"jNO2: {'not in the input file — ' if c.source == 'file' else ''}jNO2 at full lamp power ({jfull:g} 1/s) × lights")
        if c.source == "file":
            have += [row for row in J_INPUTS[1:] if _file_column(c, J_ENZ_NAME[row[0]], grid, 0.0) is not None]
        missing = [J_ENZ_NAME[row[0]] for row in J_INPUTS if row not in have]
        if missing:
            notes.append(f"{', '.join(missing)}: from the lamp J table (jNO2 × ratio of this light source)")
        return have, notes
    if c.j_source != "measured" or c.source != "file":
        return list(J_INPUTS), []
    grid = time_grid(c) if grid is None else grid
    have = [row for row in J_INPUTS if _file_column(c, J_ENZ_NAME[row[0]], grid, 0.0) is not None]
    if not any(row[0] == "jno2" for row in have):
        return [], ["no jNO2 in the input file: photolysis from the clear-sky model (nothing to scale the mechanism's J values with)"]
    missing = [J_ENZ_NAME[row[0]] for row in J_INPUTS if row not in have]
    notes = [f"{', '.join(missing)}: not in the input file — modelled from the measured jNO2 (MCM clear-sky shape)"] if missing else []
    return have, notes


def input_sources(spec: ExperimentSpec) -> Dict[str, str]:
    """Where each physical input of the run comes from: 'measured (archive)', 'measured (file)',
    'constant <value>', 'clear-sky model', 'roof times' — the readiness report of the experiment."""
    c = spec.conditions
    where = {"archive": "measured (archive)", "file": "measured (file)"}.get(c.source, "")
    grid = time_grid(c) if (c.date and c.start and c.end) else np.zeros(2)
    out: Dict[str, str] = {}
    for enz, label, const, unit in (("t", "temperature", c.temperature_K, "K"), ("H2O", "humidity", c.h2o_percent, "%"),
                                    ("dilute", "flow", c.flow_m3h, "m³/h")):
        arr = _measured_column(c, enz, grid, const) if where else None
        out[label] = where if arr is not None else f"constant {const:g} {unit}"
    out["pressure"] = (where if (c.source == "file" and _file_column(c, "p", grid, c.pressure_hPa) is not None)
                       else f"constant {c.pressure_hPa:g} hPa")
    light = _profile(c).light_type
    if light == "lamps":
        from juice.chemistry import lamp_photolysis as lp
        rows, j_notes = measured_j_rows(c, grid)
        meas = [J_ENZ_NAME[r[0]] for r in rows]
        lj = lp.lamp_j(_profile(c))
        out["photolysis"] = (f"lamps — {where + ': ' if where else ''}{', '.join(meas)}; " if where else "lamps — jNO2 at full power × lights; ") \
            + f"other J values: {lp.summary(lj)}"
    elif light == "dark":
        out["photolysis"] = "none (dark chamber profile)"
    elif c.j_source == "measured" and where:
        have = [j for j in J_ENZ_NAME.values() if _measured_column(c, j, grid, 0.0) is not None]
        if c.source == "file":
            rows, _n = measured_j_rows(c, grid)
            if not rows:
                out["photolysis"] = "clear-sky model — the input file has no jNO2 (give at least jNO2 to use measured light)"
            else:
                missing = [j for j in J_ENZ_NAME.values() if j not in have]
                out["photolysis"] = f"{where}: " + ", ".join(have) + (f"; modelled from jNO2: {', '.join(missing)}" if missing else "")
        else:
            out["photolysis"] = (f"{where}: " + ", ".join(have)) if have else "NO measured J values found — columns would be 0"
    else:
        out["photolysis"] = "clear-sky model (jNO2 from the solar zenith angle, the rest from the mechanism)"
    key = "lights" if light != "sun" else "roof"
    out[key] = (f"{where} ({key} column)" if (c.source == "file" and _file_column(c, "hv", grid, 0.0) is not None)
                else (", ".join(f"{o}–{cl}" for o, cl in c.roof_open) or ("never open (dark)" if light == "sun" else "never on (dark)")))
    if c.constrain_species:
        missing = [s for s in c.constrain_species if _measured_column(c, s, grid, 0.0) is None]
        out["held to measurements"] = ", ".join(c.constrain_species) + (f" — MISSING: {', '.join(missing)}" if missing else "")
    return out


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
    cols["sz"] = _solar_zenith(grid, c)
    cols["hv"] = _roof(grid, c)
    light = _profile(c).light_type                            # sun: roof 0/1; lamps: lights 0–1; dark: 0

    def phys(name: str, const: float) -> np.ndarray:
        if c.source in ("archive", "file"):
            arr = _measured_column(c, name, grid, const)
            if arr is not None:
                return arr
            notes.append(f"{name}: not in the {'archive' if c.source == 'archive' else 'input file'}, constant {const} used")
        return np.full(n, float(const))

    if c.source == "file":                                 # a roof column in the file wins over typed roof times
        hv = _file_column(c, "hv", grid, 0.0)
        if hv is not None:
            cols["hv"] = (np.clip(np.asarray(hv, float), 0.0, 1.0) if light == "lamps"     # lamps may be dimmed: 0–1
                          else (np.asarray(hv) > 0.5).astype(float))
        sz = _file_column(c, "sz", grid, np.nan)           # a tabulated solar zenith angle (e.g. an ENZ) wins over the computed one
        if sz is not None and np.isfinite(sz).all():
            cols["sz"] = np.asarray(sz, float)
    if light == "dark" and np.any(cols["hv"] > 0):
        notes.append("dark chamber: the light switch is forced to 0 (profile light type 'dark')")
        cols["hv"] = np.zeros(n)
    cols["dilute"] = phys("dilute", c.flow_m3h)
    cols["t"] = phys("t", c.temperature_K)
    cols["H2O"] = phys("H2O", c.h2o_percent)
    for inj in spec.injections:
        col = inj.enz_column()
        cols[col] = cols.get(col, np.zeros(n)) + pulse_column(grid, c, inj.start, inj.duration_min)
    if c.source in ("archive", "file"):
        for name in ("NO2", "NO", "O3"):
            arr = _measured_column(c, name, grid, 0.0)
            if arr is not None:
                cols[name] = arr
    for sp in c.constrain_species:
        if sp not in cols:
            notes.append(f"constrained species {sp} has no ENZ column")
    if c.j_source == "measured" or light == "lamps":
        rows, j_notes = measured_j_rows(c, grid)
        notes.extend(j_notes)
        for jname, jcol, _rxn in rows:
            enz_name = J_ENZ_NAME[jname]
            arr = _measured_column(c, enz_name, grid, 0.0) if c.source in ("archive", "file") else None
            if arr is None and light == "lamps" and jname == "jno2":
                arr = float(_profile(c).light.get("jno2_full_s1", 0) or 0) * cols["hv"]   # noted by measured_j_rows
            elif arr is None:
                notes.append(f"{enz_name}: no measured J values, column filled with 0")
                arr = np.zeros(n)
            cols[enz_name] = arr
    return grid, cols, notes


def write_enz_file(path: str, grid: np.ndarray, cols: Dict[str, np.ndarray]) -> str:
    from juice.core.enz_builder import write_enz
    write_enz(path, grid, cols)
    return path


# ─────────────────────────────────────────────────────────────────────────────
# EASY form
# ─────────────────────────────────────────────────────────────────────────────

def is_gecko_base(base_id: str) -> bool:
    return mreg.is_gecko(base_id)


def gecko_includes(spec: ExperimentSpec) -> Tuple[List[str], List[str], str]:
    """(includes_pre, includes_other, dilution stem) of a GECKO-A base: the physics lines (M, O2, N2 - an MCM
    include set carries them in its .saphir.incl) and the alias file first, the GECKO-A inorganic core when it is
    the chosen core, then names, organic chemistry and (optionally) partitioning."""
    b = mreg.base(spec.base_mechanism)
    fl = b.files
    pre = ([fl["physics"]] if "physics" in fl else []) + [fl["aliases"]] + ([fl["anorg"]] if (spec.gecko_core or "gecko") == "gecko" else [])
    other = [fl["names"], fl["organic"]] + ([fl["soa"]] if spec.gecko_partitioning and "soa" in fl else [])
    return pre, other, b.id


# MCM clear-sky parameters (l, m, n of J = l*cos^m(SZA)*exp(n*sec SZA)) of the six measured-J channels GECKO-A's
# GPHOT(...) anchors on, as in anorg.mcm331.phot.rate / default_minimum.incl (MCM v3.3.1).
GECKO_ANCHOR_PHTFKT = {"jo1d": "6.073D-05,1.743,-0.474", "jno2": "1.165D-02,0.244,-0.267", "jh2o2": "1.041D-05,0.723,-0.279",
                       "jhono": "2.644D-03,0.261,-0.288", "jhcho_r": "4.642D-05,0.762,-0.353", "jhcho_m": "6.853D-05,0.477,-0.323"}


def gecko_form(f: EasyForm, spec: ExperimentSpec, lamps: bool = False) -> None:
    """Apply a GECKO-A base to a form whose photolysis rows are already set: include files, the .dil pair, the
    inorganic core choice (GECKO-A's own core drops the MCM anorg set and the SAPHIR measured-J include with its MCM
    rate overrides; j_clear_sky stays for PHTFKT), and an anchor for every measured-J name GECKO-A's GPHOT(...) refers
    to that this run does not measure: under the sun the MCM clear-sky PHTFKT(l,m,n) (j_clear_sky.incl: with the HV
    switch and the measured-jNO2 scaling), under lamps the lamp J<n>.  Never CONST(J<n>) under the sun: J1 etc. exist
    only with the MCM phot include, and an undefined name is not zero in every path."""
    from juice.chemistry import lamp_photolysis as lp
    pre, other, stem = gecko_includes(spec)
    ours = set(pre) | set(other)
    f.includes_pre = pre + [x for x in f.includes_pre if x not in ours]
    f.includes_other = other + [x for x in f.includes_other if x not in ours]
    f.mech_stem = ""
    f.dil_stem = stem
    core = spec.gecko_core or "gecko"
    f.use_anorg = core == "mcm"
    f.use_j_clear_sky = True
    if core == "gecko":
        f.jmeasured_include = ""
        f.j_overrides = {}
    have = {r.name for r in f.j_rows if r.mode != "off"}
    for j, n in lp.MEASURED_J.items():
        if j in have:
            continue
        if lamps:
            f.extra_statements.append(f"{j} = CONST(J{n})   ; GECKO-A photolysis anchor, not measured here: the lamp J{n}")
        else:
            f.extra_statements.append(f"{j} = PHTFKT({GECKO_ANCHOR_PHTFKT[j]})   ; GECKO-A photolysis anchor, not measured here: "
                                      f"MCM clear-sky J{n}, HV and measured-jNO2 scaling as in j_clear_sky.incl")


# Lamp chambers: F0AM channels without an MCM number (JN<k>) take the lamp J of the MCM channel that absorbs in
# the same region, times their clear-sky ratio at 30 deg (approximate; noted in the EASY file).
JN_LAMP_REFERENCE = {1: 18, 2: 18, 3: 22, 4: 41, 14: 41, 15: 41, 16: 41, 17: 41, 21: 3, 22: 3, 52: 18}


def jn_lamp_lines(channels) -> List[str]:
    """JN<k> under lamps: the lamp J of the MCM channel absorbing in the same region x the clear-sky ratio at 30 deg."""
    here = _JUICE_ROOT
    mcm = _phtfkt_table(os.path.join(here, "chemistry", "photolysis", "mcm_photolysis.csv"), "J")
    f0am = _phtfkt_table(os.path.join(here, "chemistry", "photolysis", "f0am_jn_clear_sky.csv"), "Jn")
    out = []
    for k in channels:
        if k not in f0am:
            out.append(f"JN{k} = CONST(0.0)   ; F0AM channel Jn{k} is zero in its clear-sky table"); continue
        ref = JN_LAMP_REFERENCE.get(k, 41)
        ratio = _clear_sky(f0am[k]) / _clear_sky(mcm[ref]) if ref in mcm else 0.0
        out.append(f"JN{k} = CONST(J{ref}*{ratio:.4e})   ; approximate: lamp J{ref} x clear-sky Jn{k}/J{ref} at 30 deg")
    return out


def _phtfkt_table(path: str, key_prefix: str) -> Dict[int, Tuple[float, float, float]]:
    """{number: (l, m, n)} of a J table (chemistry/photolysis/mcm_photolysis.csv: J<n>; f0am_jn_clear_sky.csv: Jn<k>)."""
    import csv
    out: Dict[int, Tuple[float, float, float]] = {}
    with open(path, encoding="utf-8") as fh:
        rows = csv.DictReader(l for l in fh if not l.startswith("#"))
        for r in rows:
            m = re.fullmatch(key_prefix + r"(\d+)", (r.get("J") or "").strip())
            if m and r.get("l"):
                out[int(m.group(1))] = (float(r["l"]), float(r["m"]), float(r["n"]))
    return out


def _clear_sky(p: Tuple[float, float, float], sza: float = 30.0) -> float:
    c = math.cos(math.radians(sza))
    return p[0] * c ** p[1] * math.exp(-p[2] / c)


def lumped_form(f: EasyForm, spec: ExperimentSpec, lamps: bool = False, lamp_lines: Optional[List[str]] = None) -> List[str]:
    """Apply a lumped base (a mechanism of the F0AM library, f0am_import.py) to a form whose photolysis rows are set:
    its own includes (physics, aliases, chemistry, J names) and dilution pair, no MCM inorganic set, no SAPHIR
    measured-J include.  Its J names refer to J<n> (MCM numbers) and JN<k> (F0AM channels), defined here for the
    light of the run: under the sun the MCM clear-sky PHTFKT(l,m,n) (HV switch and measured-jNO2 scaling of
    j_clear_sky.incl) or the measured J where the run has it; under lamps the lamp J<n> (the k[...] lines of
    the MCM inorganic reactions are dropped) and JN<k> scaled from a lamp J<n>.  Returns the lamp lines to use."""
    from juice.chemistry import lamp_photolysis as lp
    b = mreg.base(spec.base_mechanism)
    fl, info = b.files, b.lumped or {}
    pre, other = [fl["physics"], fl["aliases"]], [fl["chemistry"], fl["phot"]]
    ours = set(pre) | set(other)
    f.includes_pre = pre + [x for x in f.includes_pre if x not in ours]
    f.includes_other = other + [x for x in f.includes_other if x not in ours]
    f.mech_stem, f.dil_stem = "", b.id
    f.use_anorg, f.use_j_clear_sky = False, True
    f.jmeasured_include, f.j_overrides = "", {}
    here = _JUICE_ROOT
    mcm = _phtfkt_table(os.path.join(here, "chemistry", "photolysis", "mcm_photolysis.csv"), "J")
    f0am = _phtfkt_table(os.path.join(here, "chemistry", "photolysis", "f0am_jn_clear_sky.csv"), "Jn")
    lines = [f"; PHOTOLYSIS of {b.id}: J<n> = MCM v3.3.1 J number n, JN<k> = F0AM channel Jn<k> "
             f"({'lamps' if lamps else 'sun'})"]
    measured = {lp.MEASURED_J[r.name]: r.name for r in f.j_rows if r.mode == "enz" and r.name in lp.MEASURED_J}
    lamp_out = [l for l in (lamp_lines or []) if not l.lstrip().startswith("k[")]
    defined = {int(m.group(1)) for l in lamp_out for m in [re.match(r"\s*J(\d+)\s*=", l)] if m}
    for n in info.get("j_mcm", []):
        if lamps:
            if n not in defined:
                lines.append(f"J{n} = CONST(0.0)   ; no lamp J{n}: this light source gives no value for it")
        elif n in measured:
            lines.append(f"J{n} = CONST({measured[n]})   ; measured")
        elif n in mcm:
            l_, m_, n_ = mcm[n]
            lines.append(f"J{n} = PHTFKT({l_:.4e},{m_:g},{-n_:g})")
        else:
            lines.append(f"J{n} = CONST(0.0)   ; no MCM parameters for J{n}")
    if lamps:
        lines += jn_lamp_lines(info.get("j_f0am", []))
    else:
        for k in info.get("j_f0am", []):
            if k not in f0am:
                lines.append(f"JN{k} = CONST(0.0)   ; F0AM channel Jn{k} is zero in its clear-sky table"); continue
            l_, m_, n_ = f0am[k]
            lines.append(f"JN{k} = PHTFKT({l_:.4e},{m_:g},{-n_:g})   ; F0AM clear-sky table, fitted (chemistry/photolysis/f0am_jn_clear_sky.csv)")
    if lamps:                                          # after the lamp J<n> they refer to
        return lamp_out + lines
    f.extra_statements += lines
    return lamp_out


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
    p = _profile(c)
    lamps = p.light_type == "lamps"
    rows = measured_j_rows(c)[0] if (c.j_source == "measured" or lamps) else []
    if rows:
        f.j_rows = [MeasuredRow(j, "enz", col) for j, col, _ in rows]
        f.j_overrides = {j: rxn for j, _, rxn in rows}
    else:                                              # clear sky — also a user file without jNO2
        f.j_rows = [MeasuredRow("jno2", "const", "", CLEAR_SKY_JNO2)]
        f.j_overrides = {"jno2": J_INPUTS[0][2]}
    lamp_lines: List[str] = []
    if lamps:
        # no saphir.jmeasured_02.incl: it swaps NO2+hv-->NO+O for NO2+hv-->O3+NO and re-defines NO3/HNO3/
        # CH3CHO photolysis with the solar shape (PHOTMMSC).  The lamp block binds every J<n> to jNO2 ×
        # the ratio of this light source; the inorganic reactions keep their mechanism text.
        from juice.chemistry import lamp_photolysis as lp
        f.jmeasured_include = ""
        f.j_overrides["jno2"] = "NO2 + HV --> NO + O"
        lamp_lines = lp.easy_lines(lp.lamp_j(p), measured={j: lp.MEASURED_J[j] for j, _c, _r in rows if j in lp.MEASURED_J},
                                   profile_id=p.id)
    wall_lines: List[str] = []
    if ch.wall_loss and ch.wall_loss_s1:
        wall_lines = [f"; WALL LOSS (chamber profile {p.id}): first-order loss to the walls, k in 1/s"]
        for sp, k in ch.wall_loss_s1.items():
            wall_lines += [f"{sp} -->", f"k[{sp} --> ] = CONST({float(k):g})"]
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
    if is_gecko_base(spec.base_mechanism):
        gecko_form(f, spec, lamps=lamps)
    elif mreg.is_lumped(spec.base_mechanism):
        lamp_lines = lumped_form(f, spec, lamps=lamps, lamp_lines=lamp_lines)
    if lamps and spec.addons:                          # add-ons with F0AM channels: their JN<k> under the lamps
        try:
            jn = sorted({k for a in mreg.addons_for(spec.base_mechanism, list(spec.addons)) for k in a.j_f0am})
        except (KeyError, ValueError):
            jn = []
        if jn:
            lamp_lines = list(lamp_lines) + ["; add-on photolysis channels under the lamps (last definition wins)"] + jn_lamp_lines(jn)
    try:                                               # every generated EASY states its chemistry
        f.extra_statements += mreg.provenance_lines(spec.base_mechanism, list(spec.addons))
    except Exception as e:
        f.extra_statements.append(f"; CHEMISTRY: provenance unavailable ({e})")
    f.extra_statements += lamp_lines + wall_lines
    for sp in c.constrain_species:
        f.extra_statements += [f"{sp}_in=input({sp})", f"{sp}=CONST({sp}_in*1e-9*M)"]
    for rxn, rate, cmt in spec.reaction_adjustments:
        if rxn.strip():
            f.extra_statements.append(f"k[{rxn.strip()}] = CONST({(rate or '0').strip()})" + (f"   ;{cmt.strip()}" if cmt.strip() else ""))
    f.extra_statements += list(spec.overrides)
    return f


def render_easy(spec: ExperimentSpec) -> str:
    return render_form(easy_form(spec))


# ─────────────────────────────────────────────────────────────────────────────
# Backend: mechanism subset + file writing
# ─────────────────────────────────────────────────────────────────────────────


def base_fac(base_id: str = "") -> str:
    """Path of the base mechanism's .fac (registry; falls back to the old MCM_database layout)."""
    try:
        return mreg.base(base_id).fac
    except KeyError:
        return os.path.join(_JUICE_ROOT,
                            "MCM_database", "MCM_3.3.1", "mcm_3-3-1_facsimile_complete.fac")


FULL_MCM_FAC = base_fac()


def ensure_mechanism(spec: ExperimentSpec, out_dir: str, full_fac: str = "") -> dict:
    """Extract the base-mechanism subset for spec.precursors into <out_dir>/chemistry, copy the
    chosen add-ons next to it and register them as includes, and point the spec at it.
    No-op when spec.mechanism_stem is already set (existing include set)."""
    if is_gecko_base(spec.base_mechanism):
        return _ensure_gecko(spec, out_dir)
    if mreg.is_lumped(spec.base_mechanism):
        return _ensure_lumped(spec, out_dir)
    if spec.mechanism_stem or not (spec.precursors or spec.inorganic_only):
        return {}
    import io, contextlib, shutil
    from juice.chemistry.mcm_converter import convert_mcm_file
    chem = os.path.join(out_dir, "chemistry")
    addons = mreg.addons_for(spec.base_mechanism, list(spec.addons)) if spec.addons else []
    seed = list(spec.precursors) or ["CO"]             # inorganic only: any inorganic seed → init mode, the converter adds the whole inorganic core
    for a in addons:                                   # what the add-on consumes must be in the subset
        seed += [sp for sp in a.needs_base_species if sp not in seed]
    with contextlib.redirect_stdout(io.StringIO()):
        res = convert_mcm_file(full_fac or base_fac(spec.base_mechanism), chem, saphir=True,
                               measured=(spec.conditions.j_source == "measured"),
                               precursors=seed)
    spec.mechanism_stem = os.path.basename(str(res["mech"]))[:-len(".mech")]
    spec.mech_incl_suffix = ".incl"
    spec.use_anorg = False               # inorganic core is inlined in init mode
    _ensure_comprc(os.path.join(chem, spec.mechanism_stem + ".incl"))
    for a in addons:                                   # add-ons after the base: last definition wins
        for p in a.paths():
            shutil.copy2(p, os.path.join(chem, os.path.basename(p)))
            if os.path.basename(p) not in spec.extra_includes:
                spec.extra_includes.append(os.path.basename(p))
    ro2_new = [s for a in addons for s in a.ro2]
    if ro2_new:                                        # the add-ons' peroxy radicals count in RO2, as in their source
        ro2_incl = os.path.join(chem, spec.mechanism_stem + ".ro2.incl")
        if os.path.isfile(ro2_incl):
            text = open(ro2_incl, encoding="latin-1").read().rstrip()
            have = set(re.findall(r"[A-Za-z0-9_]+", text))
            extra = [s for s in ro2_new if s not in have]
            if extra and text.endswith(")"):
                with open(ro2_incl, "w", encoding="latin-1", newline="\n") as fh:
                    fh.write(text[:-1] + " + " + " + ".join(extra) + ")\n")
    res["addons"] = [a.id for a in addons]
    return res


def _ensure_gecko(spec: ExperimentSpec, out_dir: str) -> dict:
    """A GECKO-A base: its include files are copied into <out_dir>/chemistry (the run folder stays self-contained,
    FACSIMILE and JUICE find them next to the EASY); the form takes its include lines from gecko_includes()."""
    import shutil
    b = mreg.base(spec.base_mechanism)
    chem = os.path.join(out_dir, "chemistry")
    os.makedirs(chem, exist_ok=True)
    for role, fn in b.files.items():
        src = os.path.join(b.folder, fn)
        if os.path.isfile(src):
            shutil.copy2(src, os.path.join(chem, fn))
    spec.mechanism_stem = ""
    spec.inorganic_only = False
    voc = (b.gecko or {}).get("precursor", "")
    if voc and voc.upper() not in {p.upper() for p in spec.precursors}:
        spec.precursors = [voc] + list(spec.precursors)
    return {"mech": os.path.join(chem, b.files["organic"]), "gecko": b.id, "species": (b.gecko or {}).get("n_species")}


def _ensure_lumped(spec: ExperimentSpec, out_dir: str) -> dict:
    """A lumped base: its include files are copied into <out_dir>/chemistry; the whole mechanism is used."""
    import shutil
    b = mreg.base(spec.base_mechanism)
    chem = os.path.join(out_dir, "chemistry")
    os.makedirs(chem, exist_ok=True)
    for fn in b.files.values():
        src = os.path.join(b.folder, fn)
        if os.path.isfile(src):
            shutil.copy2(src, os.path.join(chem, fn))
    spec.mechanism_stem, spec.inorganic_only = "", False
    return {"mech": os.path.join(chem, b.files["chemistry"]), "lumped": b.id, "species": (b.lumped or {}).get("n_species")}


# FACSIMILE fall-off function; normally defined in anorg.mcm331.dark.rate, which the
# init-mode subset does not include (its inorganic core is inlined) — so the subset's
# .incl must carry the definition before the first KMT/KFPAN line that uses it.
COMPRC_DEF = """; Troe fall-off (copied from anorg.mcm331.dark.rate; init-mode subset has no anorg files)
COMPRC = FUNCTION[(%1)*(%2)/((%1)+(%2))*$
                  (%3)@(1/(1+(LOG10((%1)/(%2))/$
                  ((0.75-1.27*LOG10(%3))))@2))]

"""


def _ensure_comprc(incl_path: str) -> None:
    if not os.path.isfile(incl_path):
        return
    with open(incl_path, encoding="latin-1") as fh:
        lines = fh.readlines()
    if any(re.match(r"\s*COMPRC\s*=\s*FUNCTION", l) for l in lines):
        return
    first_use = next((i for i, l in enumerate(lines) if "COMPRC(" in l and not l.lstrip().startswith(";")), None)
    if first_use is None:
        return
    lines[first_use:first_use] = [COMPRC_DEF]
    with open(incl_path, "w", encoding="latin-1") as fh:
        fh.writelines(lines)


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
        e["delivery"] = dict(e.get("delivery") or {})
        ev = InjectionEvent(**e)
        if not ev.action:
            ev.action = ev.action_key()
        inj.append(ev)
    spec = ExperimentSpec(**{k: v for k, v in d.items() if k in ExperimentSpec.__dataclass_fields__})
    spec.conditions, spec.chamber, spec.injections = cond, ch, inj
    spec.initials = [tuple(x) for x in spec.initials]
    spec.reaction_adjustments = [tuple(x) for x in spec.reaction_adjustments]
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


TEMPLATE_DIR = os.path.join(_JUICE_ROOT, "setups", "templates")


def list_templates() -> Dict[str, str]:
    """name → path of the experiment templates (setups/templates/*.experiment.json)."""
    out = {}
    if os.path.isdir(TEMPLATE_DIR):
        for f in sorted(os.listdir(TEMPLATE_DIR)):
            if f.endswith(".experiment.json"):
                out[f[:-len(".experiment.json")]] = os.path.join(TEMPLATE_DIR, f)
    return out


def load_template(path: str) -> Tuple[str, str, ExperimentSpec]:
    """Return (title, description, spec) of a setups/templates/*.experiment.json file.
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


def injection_label(col: str, spec: Optional["ExperimentSpec"] = None) -> str:
    """Plot label for an injection switch column: 'QNO' → 'NO injection', 'Q_ethanol' →
    'ethanol injection'.  Only for display — files keep the Q column names.  With a spec, the
    name the user gave the injection is used when it maps to this column."""
    if spec is not None:
        names = [inj.name for inj in spec.injections if inj.name and inj.enz_column() == col]
        if names:
            return " / ".join(dict.fromkeys(names)) + " injection"
    base = col[2:] if col[:2].upper() == "Q_" else (col[1:] if col[:1].upper() == "Q" else col)
    return (base or col) + " injection"


LEVEL_ICON = {"error": "✖ ", "warning": "⚠ ", "info": "ℹ "}


def level_icon(level: str) -> str:
    """The mark a validation level is shown with in the GUIs."""
    return LEVEL_ICON.get(level, "⚠ ")


def validate(spec: ExperimentSpec, known_species: Optional[set] = None) -> List[Tuple[str, str]]:
    """Return (level, message) items; level is 'error' (must be fixed), 'warning' (check it) or
    'info' (a note on what the experiment is, nothing to fix — e.g. a dark experiment)."""
    msgs: List[Tuple[str, str]] = []
    c = spec.conditions
    if not (c.date and c.start and c.end):
        msgs.append(("error", "Enter the experiment date and the time window (Where and when)."))
    else:
        try:
            grid = time_grid(c)
            if len(grid) < 2:
                msgs.append(("error", "Time window is empty — check start and end time."))
        except Exception as e:
            msgs.append(("error", f"Date or time is not valid: {e}"))
    if c.source not in ("archive", "manual", "file"):
        msgs.append(("error", "Choose where the chamber conditions come from: measured (archive), measured (file) or manual."))
    if c.source == "file":
        if not c.input_file:
            msgs.append(("error", "Measured inputs from a file: choose the file (CSV or Excel; 'Template…' writes an example)."))
        elif not os.path.isfile(c.input_file):
            msgs.append(("error", f"Measured inputs file not found: {c.input_file}"))
        else:
            try:
                from juice.core import measured_inputs as mi
                inp = mi.cached_inputs(c.input_file, c.date, c.pressure_hPa)
                if inp is not None and not inp.cols:
                    msgs.append(("error", "Measured inputs file: no column recognised — see the template for the column names."))
                elif inp is not None:
                    for note in inp.notes:
                        if "cannot" in note:
                            msgs.append(("warning", "Measured inputs file: " + note))
            except Exception as e:
                msgs.append(("error", f"Measured inputs file cannot be read: {e}"))
    light = _profile(c).light_type
    dark_run = not c.roof_open            # roof never open / lights never on: photolysis settings do not matter
    try:
        src = input_sources(spec)
        if src.get("photolysis", "").startswith("NO measured") and not dark_run:
            msgs.append(("warning", "Photolysis: measured J values selected but none found — the J columns would be 0 (switch to clear sky or provide them)."))
        if c.source in ("archive", "file"):
            consts = [k for k in ("temperature", "humidity", "flow") if src.get(k, "").startswith("constant")]
            if consts:
                msgs.append(("warning", "Inputs: " + ", ".join(consts) + " not measured — constant value(s) used."))
        if "held to measurements" in src and "MISSING" in src["held to measurements"]:
            msgs.append(("error", "Held species without measurements: " + src["held to measurements"].split("MISSING: ")[1]))
    except Exception:
        pass
    if not spec.mechanism_stem and not spec.precursors and not spec.inorganic_only and not mreg.is_complete(spec.base_mechanism):
        msgs.append(("error", "No chemistry: choose compounds, an existing mechanism, or 'Inorganic chemistry only'."))
    if c.j_source == "measured" and c.source == "manual" and light == "sun" and not dark_run:
        msgs.append(("warning", "Measured J values selected without archive data: the J columns will be 0 "
                                "(no photolysis). Use 'clear sky' or provide archive data."))
    if light == "lamps" and not dark_run and not float(_profile(c).light.get("jno2_full_s1", 0) or 0) and c.source != "file":
        msgs.append(("error", "The chamber profile gives no jNO2 at full lamp power and no file with jNO2 is used: no photolysis."))
    if dark_run and light != "dark":
        msgs.append(("info", "Roof never opens: run as a dark experiment (no photolysis)." if light == "sun"
                     else "Lights never on: run as a dark experiment (no photolysis)."))
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
    for note in spec.import_notes:
        if note.lstrip("⚠✖ ").startswith("no roof open/close"):   # older wiki imports: covered by the dark-experiment note above
            continue
        msgs.append(("warning", "Import: " + note.lstrip("⚠✖ ").strip()))
    return msgs


def known_species_for(spec: ExperimentSpec, include_dir: str, out_dir: str = "") -> set:
    """Species an injection may produce without a warning: the standard chamber species, the species of
    an existing include set (mechanism_stem) or — compounds mode — the chosen compounds, the add-on
    species, and whatever a previously generated chemistry/ folder of the experiment contains."""
    names = mechanism_species(include_dir, spec.mechanism_stem) if spec.mechanism_stem else mechanism_species("", "")
    for d in (spec.include_dirs if spec.mechanism_stem else []):       # an imported experiment's own include folders
        if d and os.path.isdir(d):
            names |= mechanism_species(d, spec.mechanism_stem)
    names |= {p.upper() for p in spec.precursors}
    names |= {"NO", "NO2", "O3", "CO", "CO2", "HONO", "H2O2", "HNO3", "SO2", "H2", "CH4", "HCHO"}
    for inj in spec.injections:                        # holding species of a delayed release are defined by the chain itself
        for lhs, rhs, _rate in inj.chain:
            names |= {t.upper() for t in re.findall(r"[A-Za-z][A-Za-z0-9_]*", f"{lhs}+{rhs}")}
    try:
        names |= {s.upper() for s in mreg.addon_species(list(spec.addons))}
        names |= set(mreg.species_names(spec.base_mechanism))          # a GECKO-A base: its species table
    except Exception:
        pass
    chem = os.path.join(out_dir, "chemistry") if out_dir else ""
    if chem and os.path.isdir(chem):
        for f in os.listdir(chem):
            if f.lower().endswith((".mech", ".incl")):
                try:
                    for line in open(os.path.join(chem, f), encoding="latin-1"):
                        code = line.split(";", 1)[0]
                        if "-->" in code:
                            for tok in re.findall(r"[A-Za-z][A-Za-z0-9_]*", code.replace("-->", "+")):
                                names.add(tok.upper())
                except OSError:
                    pass
    return names


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
    ch = ChamberEffects(light_source=True, hcho_const_name="SHCHO")   # scene_09's mismatched constant name, kept so the golden stays statement-equivalent
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
    here = _JUICE_ROOT
    tmp = tempfile.mkdtemp(prefix="experiment_spec_")
    bad = 0
    for spec in (spec_scene_05(), spec_scene_09()):
        sc = spec.name
        res = write_experiment(spec, os.path.join(tmp, sc))
        orig_easy = os.path.join(here, "validation", "IDL Input", sc, sc + ".easy")
        ok, diffs = statements_equivalent(orig_easy, res["easy"], ignore_duplicate_reactions=True)
        # ENZ: injection switch columns and roof must match the archived ENZ exactly
        orig_enz = os.path.join(here, "validation", "IDL Input", sc, sc + ".enz")
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
