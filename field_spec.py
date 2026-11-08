"""
field_spec.py — a field case as data: what is modelled, and the one table that drives it (no Qt).

A field case is a small JSON file (<name>.field.json, a `FieldSpec`) plus one input table (CSV, or
the first sheet of an Excel file) with one row per time.  Three kinds of run share that table:

    constrained   a site: measured species are held to the data, radicals and products run free
    plume         a plume or transect: chemistry ages while the air is diluted towards background
    lagrangian    a parcel along a trajectory: emissions, boundary-layer growth, background
                  dC/dt = E(t)/(h(t)·A) + chemistry − k(t)·(C − C_background)     (Stockwell et al., 2024)

The table
    time         as in measured_inputs: 'YYYY-MM-DD HH:MM' (UTC), 'HH:MM' or hours of day (these two
                 need the case's date), unix seconds, seconds since 2000-01-01
    conditions   T  p  H2O or RH  sz  lat  lon  h  kdil
    photolysis   any jNO2-style column: j<name> or J_<name>, 1/s.  The names belong to the mechanism
    E_<species>  emission     unit in brackets, e.g. "E_NO [mol/s]"; else the case's default unit
    B_<species>  background   mixing ratio the air is diluted towards ("B_O3 [ppb]"; ppb by default)
    C_<species>  constraint   mixing ratio the species is held to     ("C_NO [ppt]"; ppb by default)

Species names are not checked against a mechanism here and carry no MCM assumption: the mechanism
is whatever the case's `chemistry` holds (the choices of the Chemistry step: an MCM subset, or an
existing mechanism include set such as RACM2B-VCP or GECKO-A).

    spec = load_spec("x.field.json")          # or FieldSpec(...)
    table, report = check(spec)               # never raises on a bad table: report.errors says why
    emission_rates(table, spec)               # {species: molec cm-3 s-1 per row}
    dilution_rate(table, spec)                # s-1 per row
    deposition_rates(table, spec)             # {species: s-1 per row}
    write_template("t.csv", "lagrangian")

This module only reads, converts and checks.  Turning a case into a run is a separate step, which
also decides how empty cells in the conditions and constraints are bridged.
`python field_spec.py --selftest` runs the module's checks.
"""
from __future__ import annotations

import json
import os
import re
import sys
from dataclasses import asdict, dataclass, field, fields
from typing import Dict, List, Optional, Tuple

import numpy as np

import measured_inputs as mi

SPEC_VERSION = 1
SPEC_SUFFIX = ".field.json"
EXAMPLES_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "field_examples")

RUN_TYPES = {
    "constrained": "Constrained site",
    "plume": "Plume or transect",
    "lagrangian": "Lagrangian air parcel",
}
DILUTION_SOURCES = {
    "pbl": "from the growth of the boundary layer (column h)",
    "table": "from the table (column kdil)",
    "constant": "one constant rate",
    "tracer": "from the decay of a measured tracer",
    "none": "no dilution",
}
# The steps of the field window; every check says which step fixes it.
STEP_CASE, STEP_CHEMISTRY, STEP_SOURCES, STEP_MIXING, STEP_RUN = range(5)
DEFAULT_DILUTION = {"constrained": "none", "plume": "table", "lagrangian": "pbl"}

NA = 6.02214076e23
# Air density in the same form as the DENSITY function of the scenario files: M = p[hPa] / (T * 1.379e-19).
DENSITY_CONST = 1.379e-19

# unit → (needs boundary-layer height, needs footprint area, factor to the base of that unit's form)
#   ppb/s            mixing-ratio rate              → multiplied by M·1e-9
#   molec/cm3/s      already a volume rate
#   flux units       molec cm-2 s-1 after the factor → divided by h [cm]
#   mol/s, molec/s   molec s-1 after the factor      → divided by h·A [cm3]
EMISSION_UNITS: Dict[str, Tuple[bool, bool, float]] = {
    "ppb/s":       (False, False, 1.0),
    "molec/cm3/s": (False, False, 1.0),
    "molec/cm2/s": (True,  False, 1.0),
    "mol/km2/s":   (True,  False, NA / 1e10),
    "mol/km2/h":   (True,  False, NA / 1e10 / 3600.0),
    "mol/s":       (True,  True,  NA),
    "molec/s":     (True,  True,  1.0),
}
MIXING_UNITS = {"ppb": 1.0, "ppt": 1e-3, "ppm": 1e3, "molec/cm3": None}

_CONDITIONS = {
    "h":    (["h", "pbl", "pblh", "blh", "pbl_m", "mixing_height", "boundary_layer_height"], "m"),
    "kdil": (["kdil", "k_dil", "dilution_rate", "kdil_s-1"], "1/s"),
    "lat":  (["lat", "latitude"], "deg"),
    "lon":  (["lon", "long", "longitude"], "deg"),
}
_FROM_MEASURED = ("t", "p", "H2O", "RH", "sz")
_J_STANDARD = {a: k[1:] for k, (al, _u, kind) in mi.COLUMNS.items() if kind == "photolysis" for a in al}
_PREFIXED = re.compile(r"^([EeBbCc])_([^\[\]\s]+)\s*(?:\[\s*([^\]]+?)\s*\])?$")
_J_NAME = re.compile(r"^[jJ]_?([A-Z0-9][A-Za-z0-9_()]*)$")


@dataclass
class FieldSpec:
    name: str = "field_case"
    description: str = ""
    run_type: str = "constrained"
    table: str = ""                          # input table; kept relative to the .field.json when saved next to it
    date: str = ""                           # YYYY-MM-DD, only needed for 'HH:MM' / hours-of-day times
    site: str = ""
    lat: Optional[float] = None              # used when the table has no lat / lon columns
    lon: Optional[float] = None
    chemistry: Dict = field(default_factory=dict)  # base_mechanism, addons, precursors, mechanism_stem, inorganic_only,
                                                   # include_dirs, reaction_adjustments: the fields of the Chemistry step
    step_s: int = 60                         # output step of the run; inputs may be coarser
    default_emission_unit: str = "molec/cm2/s"
    area_km2: float = 0.0                    # footprint A of E/(h·A); needed by the mol/s and molec/s units
    dilution: str = "none"                   # key of DILUTION_SOURCES
    k_dil_per_h: float = 0.0                 # dilution 'constant'
    k_collapse_per_h: float = 0.3            # dilution 'pbl': rate while the layer does not grow
    tracer: str = ""                         # dilution 'tracer': species with a C_ column
    spin_up_days: float = 0.0                # constrained: repeat the first day before the output
    generic_lifetime_h: float = 24.0         # constrained: first-order loss of everything that is not held
    initials: Dict[str, float] = field(default_factory=dict)    # ppb
    deposition_velocity_cm_s: Dict[str, float] = field(default_factory=dict)   # loss = v / h
    deposition_lifetime_h: Dict[str, float] = field(default_factory=dict)      # loss = 1 / lifetime
    version: int = SPEC_VERSION


@dataclass
class FieldTable:
    path: str
    times_s: np.ndarray
    conditions: Dict[str, np.ndarray] = field(default_factory=dict)   # t K, p hPa, H2O vol-%, sz deg, lat, lon, h m, kdil 1/s
    photolysis: Dict[str, np.ndarray] = field(default_factory=dict)   # mechanism's name → 1/s
    emissions: Dict[str, np.ndarray] = field(default_factory=dict)    # species → as given in the file
    emission_units: Dict[str, str] = field(default_factory=dict)      # species → unit
    background: Dict[str, np.ndarray] = field(default_factory=dict)   # species → ppb
    constraints: Dict[str, np.ndarray] = field(default_factory=dict)  # species → ppb
    ignored: List[str] = field(default_factory=list)
    notes: List[str] = field(default_factory=list)
    problems: List[Tuple[int, str]] = field(default_factory=list)     # (step, reading error found column by column)
    unreadable: Dict[str, int] = field(default_factory=dict)          # column → cells that were not numbers

    @property
    def n(self) -> int:
        return len(self.times_s)


@dataclass
class Report:
    errors: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    info: List[str] = field(default_factory=list)
    issues: List[Tuple[str, int, str]] = field(default_factory=list)   # ("error" | "warning", step, message)
    step: int = STEP_CASE                                              # the step the next message belongs to

    def at(self, step: int) -> None:
        self.step = step

    def err(self, message: str, step: Optional[int] = None) -> None:
        self.errors.append(message); self.issues.append(("error", self.step if step is None else step, message))

    def warn(self, message: str, step: Optional[int] = None) -> None:
        self.warnings.append(message); self.issues.append(("warning", self.step if step is None else step, message))

    def worst(self, step: int) -> str:
        """'error', 'warning' or '' for one step."""
        levels = {lvl for lvl, st, _m in self.issues if st == step}
        return "error" if "error" in levels else "warning" if "warning" in levels else ""

    @property
    def ok(self) -> bool:
        return not self.errors

    def text(self) -> str:
        return "\n".join([f"ERROR    {e}" for e in self.errors] + [f"WARNING  {w}" for w in self.warnings]
                         + [f"         {i}" for i in self.info])


# ── spec file ────────────────────────────────────────────────────────────────────────────────────
def save_spec(spec: FieldSpec, path: str) -> str:
    if not path.lower().endswith(SPEC_SUFFIX):
        path += SPEC_SUFFIX
    data = asdict(spec)
    if data["table"]:
        table, folder = os.path.abspath(data["table"]), os.path.dirname(os.path.abspath(path))
        try:
            rel = os.path.relpath(table, folder)
            data["table"] = (rel if not rel.startswith("..") else table).replace("\\", "/")
        except ValueError:                       # another drive
            data["table"] = table
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
        f.write("\n")
    return path


def load_spec(path: str) -> FieldSpec:
    """Raises ValueError when the file is not a field case this JUICE can read."""
    with open(path, encoding="utf-8") as f:
        data = json.load(f)
    if not isinstance(data, dict):
        raise ValueError(f"{os.path.basename(path)} is not a JUICE field case")
    if int(data.get("version", 1)) > SPEC_VERSION:
        raise ValueError(f"{os.path.basename(path)} was written by a newer JUICE (version {data['version']})")
    known = {f.name for f in fields(FieldSpec)}
    clean = {k: v for k, v in data.items() if k in known}
    try:
        for k in ("area_km2", "k_dil_per_h", "k_collapse_per_h", "spin_up_days", "generic_lifetime_h"):
            if k in clean:
                clean[k] = float(clean[k])
        for k in ("lat", "lon"):
            if clean.get(k) is not None:
                clean[k] = float(clean[k])
        if "step_s" in clean:
            clean["step_s"] = int(clean["step_s"])
        for k in ("initials", "deposition_velocity_cm_s", "deposition_lifetime_h"):
            if k in clean:
                clean[k] = {str(a).upper(): float(b) for a, b in dict(clean[k]).items()}
        if "chemistry" in clean:
            clean["chemistry"] = dict(clean["chemistry"])
    except (TypeError, ValueError) as e:
        raise ValueError(f"{os.path.basename(path)} has a value of the wrong kind: {e}") from None
    spec = FieldSpec(**clean)
    if spec.table and not os.path.isabs(spec.table):
        spec.table = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(path)), spec.table))
    return spec


def list_examples() -> Dict[str, Tuple[str, str]]:
    """{name: (path of the .field.json, one-line description)} of the synthetic examples in field_examples/."""
    out: Dict[str, Tuple[str, str]] = {}
    if not os.path.isdir(EXAMPLES_DIR):
        return out
    for d in sorted(os.listdir(EXAMPLES_DIR)):
        p = os.path.join(EXAMPLES_DIR, d, d + SPEC_SUFFIX)
        if os.path.isfile(p):
            try:
                spec = load_spec(p)
                out[d] = (p, f"{RUN_TYPES.get(spec.run_type, spec.run_type)}: {spec.description}")
            except (OSError, ValueError):
                continue
    return out


def is_example(path: str) -> bool:
    try:
        return os.path.commonpath([os.path.abspath(path), EXAMPLES_DIR]) == EXAMPLES_DIR
    except ValueError:
        return False


def copy_case(path: str, folder: str) -> str:
    """Copy a case and its table into <folder>/<name>/ and return the path of the copied .field.json.
    The copy refers to its own table, so the original (an example, say) is never touched."""
    import shutil
    spec = load_spec(path)
    dest = os.path.join(folder, spec.name)
    os.makedirs(dest, exist_ok=True)
    if spec.table and os.path.isfile(spec.table):
        target = os.path.join(dest, os.path.basename(spec.table))
        if os.path.abspath(target) != os.path.abspath(spec.table):
            shutil.copy2(spec.table, target)
        spec.table = target
    return save_spec(spec, os.path.join(dest, spec.name))


def parse_initials(text: str) -> Dict[str, float]:
    """'O3=20, NO=0.5' → {'O3': 20.0, 'NO': 0.5} (ppb).  Raises ValueError on anything else."""
    out: Dict[str, float] = {}
    for part in re.split(r"[,;\n]+", text or ""):
        part = part.strip()
        if not part:
            continue
        m = re.fullmatch(r"([^\s=]+)\s*=\s*(\S+)", part)
        if not m:
            raise ValueError(f"'{part}' is not of the form NAME=value")
        try:
            v = float(m.group(2))
        except ValueError:
            raise ValueError(f"'{m.group(2)}' is not a number (in '{part}')") from None
        if not np.isfinite(v) or v < 0:
            raise ValueError(f"'{part}': an initial mixing ratio must be zero or positive")
        out[m.group(1).upper()] = v
    return out


def format_initials(values: Dict[str, float]) -> str:
    return ", ".join(f"{k}={v:g}" for k, v in values.items())


def chemistry_chosen(spec: FieldSpec) -> bool:
    c = spec.chemistry or {}
    return bool(c.get("precursors") or c.get("mechanism_stem") or c.get("inorganic_only"))


def chemistry_label(spec: FieldSpec) -> str:
    c = spec.chemistry or {}
    if c.get("mechanism_stem"):
        return f"mechanism include set '{c['mechanism_stem']}'"
    base = c.get("base_mechanism", "")
    if c.get("precursors"):
        return f"{base} subset for {', '.join(c['precursors'])}"
    if c.get("inorganic_only"):
        return f"{base} inorganic chemistry only"
    return "not chosen yet"


# ── reading the table ────────────────────────────────────────────────────────────────────────────
def _load_dataframe(path: str):
    """CSV with '#' comment lines, comma, semicolon or tab separated, decimal point or decimal comma
    (a semicolon table from a German Excel).  Other formats go through measured_inputs."""
    import pandas as pd
    if os.path.splitext(path)[1].lower() in (".xlsx", ".xls", ".xlsm", ".enz"):
        return mi._load_table(path)
    with open(path, encoding="utf-8-sig", errors="replace") as f:
        lines = [l for l in f if l.strip() and not l.lstrip().startswith("#")][:6]
    if not lines:
        raise ValueError("the file has no header line")
    counts = {s: sum(l.count(s) for l in lines) for s in (",", ";", "\t")}
    sep = max(counts, key=counts.get)
    decimal = "," if sep != "," and any(re.search(r"\d,\d", l) for l in lines[1:]) else "."
    return pd.read_csv(path, sep=sep, decimal=decimal, comment="#", encoding="utf-8-sig", skipinitialspace=True)


def _numeric(df, col, order, table: FieldTable) -> np.ndarray:
    import pandas as pd
    raw = df[col]
    v = pd.to_numeric(raw, errors="coerce").to_numpy(dtype=float)
    bad = int((raw.notna().to_numpy() & ~np.isfinite(v)).sum())
    if bad:
        table.unreadable[str(col)] = bad
    return v[order]


def _to_ppb(v: np.ndarray, unit: str, M: Optional[np.ndarray]) -> np.ndarray:
    if unit == "molec/cm3":
        return v / M * 1e9
    return v * MIXING_UNITS[unit]


def air_density(table: FieldTable) -> Optional[np.ndarray]:
    """molec cm-3 per row, the same density as the scenario files; None without T and p."""
    if "t" not in table.conditions or "p" not in table.conditions:
        return None
    return table.conditions["p"] / (table.conditions["t"] * DENSITY_CONST)


def read_table(path: str, date: str = "") -> FieldTable:
    """Read the table.  A file or time column that cannot be read raises ValueError; problems with
    single columns are collected in `problems` so that all of them can be shown at once."""
    if not os.path.isfile(path):
        raise ValueError(f"file not found: {path}")
    try:
        df = _load_dataframe(path)
    except Exception as e:
        raise ValueError(f"could not read {os.path.basename(path)}: {e}") from None
    if df.shape[1] < 2 or df.shape[0] < 1:
        raise ValueError(f"{os.path.basename(path)} has no data rows (a header line and at least one row are needed)")
    df = df.rename(columns=lambda c: str(c).strip())
    if df.columns.duplicated().any():
        raise ValueError("repeated column name(s): " + ", ".join(sorted(set(df.columns[df.columns.duplicated()]))))
    names = {str(c).lower(): c for c in df.columns}
    tcol = next((names[a] for a in mi.TIME_ALIASES if a in names), df.columns[0])
    try:
        times = mi._parse_time(df[tcol], date)
    except ValueError:
        raise
    except Exception as e:
        raise ValueError(f"time column '{tcol}' cannot be read: {e}") from None
    order = np.argsort(times, kind="stable")
    t = FieldTable(path=path, times_s=times[order])
    bad_t = int((~np.isfinite(times)).sum())
    if bad_t:
        t.problems.append((STEP_CASE, f"{bad_t} row(s) have a time that cannot be read"))
    if np.any(np.diff(times[np.isfinite(times)]) < 0):
        t.notes.append("rows were not in time order and have been sorted")
    used = {tcol}

    for key in _FROM_MEASURED:
        aliases = mi.COLUMNS[key][0]
        src = next((names[a] for a in aliases if a in names), None)
        if src is None:
            continue
        used.add(src)
        v = _numeric(df, src, order, t)
        med = float(np.nanmedian(v)) if np.isfinite(v).any() else float("nan")
        if key == "t" and med < 150:
            v = v + 273.15; t.notes.append("T given in °C, converted to K")
        if key == "p" and med > 5000:
            v = v / 100.0; t.notes.append("p given in Pa, converted to hPa")
        if key == "H2O" and med > 100:
            v = v / 1e4; t.notes.append("H2O given in ppm, converted to vol-%")
        t.conditions[key] = v
    for key, (aliases, _unit) in _CONDITIONS.items():
        src = next((names[a] for a in aliases if a in names), None)
        if src is not None:
            used.add(src); t.conditions[key] = _numeric(df, src, order, t)
    if "H2O" not in t.conditions and "RH" in t.conditions:
        if "t" in t.conditions and "p" in t.conditions:
            t.conditions["H2O"] = t.conditions["RH"] / 100.0 * mi._es_hpa(t.conditions["t"]) / t.conditions["p"] * 100.0
            t.notes.append("H2O computed from RH, T and p")
        else:
            t.problems.append((STEP_CASE, "RH cannot be converted to H2O without the T and p columns"))

    M = air_density(t)
    seen: Dict[Tuple[str, str], str] = {}
    for col in df.columns:
        if col in used:
            continue
        m = _PREFIXED.match(col)
        if m:
            used.add(col)
            kind, species, unit = m.group(1).upper(), m.group(2).upper(), m.group(3)
            step = STEP_MIXING if kind == "B" else STEP_SOURCES
            if (kind, species) in seen:
                t.problems.append((step, f"{species} is given twice: columns '{seen[(kind, species)]}' and '{col}'"))
                continue
            seen[(kind, species)] = col
            v = _numeric(df, col, order, t)
            if kind == "E":
                if unit is not None and unit not in EMISSION_UNITS:
                    t.problems.append((step, f"column '{col}': unknown emission unit '{unit}' (known: {', '.join(EMISSION_UNITS)})"))
                    continue
                t.emissions[species] = v
                if unit:
                    t.emission_units[species] = unit
            else:
                u = unit or "ppb"
                if u not in MIXING_UNITS:
                    t.problems.append((step, f"column '{col}': unknown mixing-ratio unit '{u}' (known: {', '.join(MIXING_UNITS)})"))
                    continue
                if u == "molec/cm3" and M is None:
                    t.problems.append((step, f"column '{col}': molec/cm3 needs the T and p columns"))
                    continue
                (t.background if kind == "B" else t.constraints)[species] = _to_ppb(v, u, M)
            continue
        low = str(col).lower()
        j = _J_STANDARD.get(low)
        if j is None:
            jm = _J_NAME.match(str(col))
            j = jm.group(1) if jm else None
        if j is not None:
            if j in t.photolysis:
                t.problems.append((STEP_CASE, f"photolysis frequency {j} is given twice (column '{col}')"))
                continue
            used.add(col)
            t.photolysis[j] = _numeric(df, col, order, t)
    t.ignored = [str(c) for c in df.columns if c not in used]
    return t


def resolve_emission_unit(table: FieldTable, spec: FieldSpec, species: str) -> str:
    return table.emission_units.get(species, spec.default_emission_unit)


def emission_rates(table: FieldTable, spec: FieldSpec) -> Dict[str, np.ndarray]:
    """Every E_ column as a volume source in molec cm-3 s-1.  Assumes a table `check` has accepted."""
    out: Dict[str, np.ndarray] = {}
    M = air_density(table)
    for sp, v in table.emissions.items():
        unit = resolve_emission_unit(table, spec, sp)
        need_h, need_a, factor = EMISSION_UNITS[unit]
        r = v * factor
        if unit == "ppb/s":
            r = r * M * 1e-9
        if need_h:
            r = r / (table.conditions["h"] * 100.0)
        if need_a:
            r = r / (spec.area_km2 * 1e10)
        out[sp] = r
    return out


def dilution_rate(table: FieldTable, spec: FieldSpec) -> np.ndarray:
    """First-order dilution in s-1 per row.  Assumes a table `check` has accepted.

    'pbl' follows Stockwell et al. (2024): the logarithmic growth of h to the next row while the
    layer grows, otherwise the collapse rate; the last row repeats the last interval.
    'tracer' is the logarithmic decay of the tracer above its background, never below zero."""
    n, d = table.n, spec.dilution
    if d == "none":
        return np.zeros(n)
    if d == "constant":
        return np.full(n, spec.k_dil_per_h / 3600.0)
    if d == "table":
        return table.conditions["kdil"].copy()
    dt = np.diff(table.times_s)
    with np.errstate(divide="ignore", invalid="ignore"):
        if d == "pbl":
            h = table.conditions["h"]
            g = np.log(h[1:] / h[:-1]) / dt
            g = np.append(g, g[-1])
            return np.where(g > 0, g, spec.k_collapse_per_h / 3600.0)
        tracer = spec.tracer.upper()
        x = table.constraints[tracer] - table.background.get(tracer, np.zeros(n))
        k = np.log(x[:-1] / x[1:]) / dt
    k = np.append(k, k[-1])
    return np.maximum(np.nan_to_num(k, nan=0.0, posinf=0.0, neginf=0.0), 0.0)


def deposition_rates(table: FieldTable, spec: FieldSpec) -> Dict[str, np.ndarray]:
    """First-order deposition in s-1 per row and species.  Assumes a table `check` has accepted."""
    out = {sp: np.full(table.n, 1.0 / (tau * 3600.0)) for sp, tau in spec.deposition_lifetime_h.items()}
    for sp, v in spec.deposition_velocity_cm_s.items():
        out[sp] = v / (table.conditions["h"] * 100.0)
    return out


# ── checking ─────────────────────────────────────────────────────────────────────────────────────
_RANGES = {"t": (180.0, 330.0, "K"), "p": (100.0, 1100.0, "hPa"), "H2O": (0.0, 8.0, "vol-%"), "sz": (0.0, 180.0, "deg")}


def _check_spec(spec: FieldSpec, r: Report) -> None:
    if spec.run_type not in RUN_TYPES:
        r.err(f"unknown kind of run '{spec.run_type}' (known: {', '.join(RUN_TYPES)})")
    if not re.fullmatch(r"[A-Za-z0-9_.\-]+", spec.name or ""):
        r.err("the case name may only contain letters, digits, _ . and -")
    if spec.date:
        try:
            import datetime as _dt
            _dt.datetime.strptime(spec.date, "%Y-%m-%d")
        except ValueError:
            r.err(f"date '{spec.date}' is not of the form YYYY-MM-DD")
    if spec.dilution not in DILUTION_SOURCES:
        r.err(f"unknown dilution source '{spec.dilution}' (known: {', '.join(DILUTION_SOURCES)})", STEP_MIXING)
    if spec.default_emission_unit not in EMISSION_UNITS:
        r.err(f"unknown emission unit '{spec.default_emission_unit}' (known: {', '.join(EMISSION_UNITS)})", STEP_SOURCES)
    if not 1 <= spec.step_s <= 3600:
        r.err("the output step must be between 1 and 3600 s")
    elif spec.step_s > 60:
        r.warn(f"output step {spec.step_s} s: the RO2 sum is frozen within a step and has been verified at 60 s")
    if spec.lat is not None and not -90 <= spec.lat <= 90:
        r.err(f"latitude {spec.lat} is outside -90 to 90")
    if spec.lon is not None and not -180 <= spec.lon <= 360:
        r.err(f"longitude {spec.lon} is outside -180 to 360")
    for name, value, step in (("area", spec.area_km2, STEP_SOURCES), ("dilution rate", spec.k_dil_per_h, STEP_MIXING),
                              ("collapse rate", spec.k_collapse_per_h, STEP_MIXING), ("spin-up", spec.spin_up_days, STEP_CASE)):
        if value < 0:
            r.err(f"the {name} cannot be negative", step)
    if spec.run_type == "constrained" and spec.generic_lifetime_h <= 0:
        r.err("the lifetime of unconstrained species must be positive", STEP_MIXING)
    if spec.dilution == "constant" and spec.k_dil_per_h <= 0:
        r.err("dilution 'constant' needs a rate above zero", STEP_MIXING)
    if not chemistry_chosen(spec):
        r.err("no chemistry chosen yet", STEP_CHEMISTRY)
    r.at(STEP_MIXING)
    for sp in set(spec.deposition_velocity_cm_s) & set(spec.deposition_lifetime_h):
        r.err(f"deposition of {sp}: give a velocity or a lifetime, not both")
    for sp, v in {**spec.deposition_velocity_cm_s, **spec.deposition_lifetime_h}.items():
        if not np.isfinite(v) or v <= 0:
            r.err(f"deposition of {sp}: the value must be above zero")
    r.at(STEP_CASE)


def _finite_or_report(v: np.ndarray, label: str, r: Report, fatal: bool, step: Optional[int] = None) -> None:
    n = int((~np.isfinite(v)).sum())
    if n:
        (r.err if fatal else r.warn)(f"{label}: {n} of {len(v)} values are empty or not numbers", step)


def _check_table(t: FieldTable, spec: FieldSpec, r: Report) -> None:
    for step, message in t.problems:
        r.err(message, step)
    r.info.extend(t.notes)
    for col, n in t.unreadable.items():
        r.warn(f"column '{col}': {n} cell(s) are text, not numbers")
    if t.n < 2:
        r.err("the table needs at least two rows")
    ts = t.times_s[np.isfinite(t.times_s)]
    if len(ts) >= 2:
        dt = np.diff(ts)
        if np.any(dt <= 0):
            r.err(f"{int((dt <= 0).sum())} row(s) repeat an earlier time")
        pos = dt[dt > 0]
        if len(pos):
            med = float(np.median(pos))
            r.info.append(f"{t.n} rows, {(ts[-1] - ts[0]) / 3600.0:.2f} h, median spacing {med / 60.0:.1f} min")
            if float(pos.max()) > 5 * med:
                r.warn(f"the largest gap between rows is {pos.max() / 60.0:.0f} min, the median {med / 60.0:.1f} min")

    c = t.conditions
    for key, label in (("t", "temperature T"), ("p", "pressure p")):
        if key not in c:
            r.err(f"column {key.upper() if key == 't' else key} ({label}) is missing")
    if "H2O" not in c:
        r.err("humidity is missing: give H2O or RH")
    for key, (lo, hi, unit) in _RANGES.items():
        if key in c:
            _finite_or_report(c[key], key, r, False)
            v = c[key][np.isfinite(c[key])]
            if len(v) and (v.min() < lo or v.max() > hi):
                r.warn(f"{key} spans {v.min():.4g} to {v.max():.4g} {unit}, outside the usual {lo:g} to {hi:g}")
    if "lat" in c and np.any(np.abs(c["lat"][np.isfinite(c["lat"])]) > 90):
        r.err("the lat column has values outside -90 to 90")
    if "lon" in c:
        v = c["lon"][np.isfinite(c["lon"])]
        if len(v) and (v.min() < -180 or v.max() > 360):
            r.err("the lon column has values outside -180 to 360")
    if ("lat" in c) != ("lon" in c):
        r.err("give both lat and lon columns, or neither")

    positioned = ("lat" in c and "lon" in c) or (spec.lat is not None and spec.lon is not None)
    if not (positioned or "sz" in c):
        r.err("the sun position is unknown: give lat and lon (table or case), or an sz column")
    if not t.photolysis and not positioned:
        r.err("no photolysis: give j columns, or a position so that clear-sky values can be used")
    elif "NO2" not in t.photolysis:
        r.warn("no jNO2 column: photolysis would be clear sky, without clouds or aerosol")
    for name, v in t.photolysis.items():
        _finite_or_report(v, f"j{name}", r, False)
        if np.any(v[np.isfinite(v)] < 0):
            r.err(f"j{name} has negative values")

    if "h" in c:
        _finite_or_report(c["h"], "h", r, False)
        if np.any(c["h"][np.isfinite(c["h"])] <= 0):
            r.err("the boundary-layer height h must be above zero")
    if "kdil" in c and np.any(c["kdil"][np.isfinite(c["kdil"])] < 0):
        r.err("kdil has negative values", STEP_MIXING)

    r.at(STEP_SOURCES)
    for sp, v in t.emissions.items():
        _finite_or_report(v, f"E_{sp}", r, True)
        if np.any(v[np.isfinite(v)] < 0):
            r.err(f"E_{sp} has negative values")
        unit = resolve_emission_unit(t, spec, sp)
        if unit in EMISSION_UNITS:
            need_h, need_a, _f = EMISSION_UNITS[unit]
            if need_h and "h" not in c:
                r.err(f"E_{sp} is in {unit}: this needs the boundary-layer height column h")
            if need_a and spec.area_km2 <= 0:
                r.err(f"E_{sp} is in {unit}: this needs the footprint area of the case")
            if unit == "ppb/s" and air_density(t) is None:
                r.err(f"E_{sp} is in ppb/s: this needs the T and p columns")
    for sp, v in t.background.items():
        _finite_or_report(v, f"B_{sp}", r, True, STEP_MIXING)
        if np.any(v[np.isfinite(v)] < 0):
            r.err(f"B_{sp} has negative values", STEP_MIXING)
    for sp, v in t.constraints.items():
        _finite_or_report(v, f"C_{sp}", r, False)
        n = int((v[np.isfinite(v)] < 0).sum())
        if n:
            r.warn(f"C_{sp}: {n} negative value(s)")

    r.at(STEP_MIXING)
    if spec.deposition_velocity_cm_s and "h" not in c:
        r.err("deposition velocities need the boundary-layer height column h; give lifetimes instead if there is none")
    d = spec.dilution
    if d == "pbl" and "h" not in c:
        r.err("dilution from the boundary layer needs the column h")
    if d == "table" and "kdil" not in c:
        r.err("dilution from the table needs the column kdil")
    if d == "tracer":
        if not spec.tracer:
            r.err("dilution from a tracer: name the tracer")
        elif spec.tracer.upper() not in t.constraints:
            r.err(f"dilution from a tracer: the table has no column C_{spec.tracer.upper()}")

    r.at(STEP_SOURCES)
    rt = spec.run_type
    if rt == "constrained":
        if not t.constraints:
            r.err("a constrained run needs at least one measured species (a C_<species> column)")
        else:
            lacking = [s for s in ("NO", "NO2", "O3") if s not in t.constraints]
            if lacking:
                r.warn("usually held at a site, but not in the table: " + ", ".join(lacking))
        if t.emissions:
            r.warn("the table has emission columns; a constrained run normally has none")
    elif rt == "plume":
        if d == "none":
            r.warn("a plume without dilution: is this intended?", STEP_MIXING)
    elif rt == "lagrangian":
        if not t.emissions:
            r.err("a Lagrangian run needs at least one emission (an E_<species> column)")
        if "h" not in c:
            r.err("a Lagrangian run needs the boundary-layer height column h", STEP_CASE)
        if not ("lat" in c and "lon" in c):
            r.warn("no lat and lon columns: the parcel is treated as staying at one place", STEP_CASE)
    if d != "none" and not t.background:
        r.info.append("no background columns: air is diluted with clean air")
    if t.ignored:
        r.warn("columns not used: " + ", ".join(t.ignored)
               + ".  Species need a prefix: E_ emission, B_ background, C_ held to the data", STEP_CASE)
    r.info.append(f"{len(t.constraints)} species held, {len(t.emissions)} emitted, {len(t.background)} with a background, "
                  f"photolysis frequencies: {len(t.photolysis)}")


def check(spec: FieldSpec) -> Tuple[Optional[FieldTable], Report]:
    """The table of a case and everything wrong with the pair.  Bad input is reported, never raised."""
    r = Report()
    _check_spec(spec, r)
    if not spec.table:
        r.err("no input table chosen", STEP_CASE)
        return None, r
    try:
        t = read_table(spec.table, spec.date)
    except ValueError as e:
        r.err(str(e), STEP_CASE)
        return None, r
    _check_table(t, spec, r)
    return t, r


# ── template ─────────────────────────────────────────────────────────────────────────────────────
_TEMPLATE_HEAD = """# JUICE field case: {title}.  One row per time, UTC.  Comma, semicolon or tab separated.
# The rows below are invented examples.  Delete them and any columns you do not have.
#
# time   'YYYY-MM-DD HH:MM' (UTC), or 'HH:MM' / hours of day together with the date of the case,
#        or unix seconds, or seconds since 2000-01-01.
#
# CONDITIONS  (names are matched ignoring case)
#   T       air temperature, K (°C is recognised when below 150)
#   p       pressure, hPa (Pa is recognised when above 5000)
#   H2O     water vapour, vol-%, or RH relative humidity, %
#   lat lon position in degrees: per row along a trajectory, or once in the case for a fixed site
#   sz      solar zenith angle, degrees (otherwise computed from time and position)
#   h       boundary-layer height, m
#   kdil    first-order dilution rate, 1/s
# PHOTOLYSIS, 1/s
#   jNO2 jO1D jHONO ... or any j<name> / J_<name>, with the names your mechanism uses.
# SPECIES  (names as in your mechanism; the unit in brackets is optional)
#   E_<species> [unit]   emission.  Units: {units}
#                        molec/cm2/s, mol/km2/s and mol/km2/h are divided by the height h;
#                        mol/s and molec/s by h and the footprint area of the case.
#   B_<species> [unit]   background the air is diluted towards.  Units: ppb (default), ppt, ppm, molec/cm3
#   C_<species> [unit]   measured value the species is held to.  Same units
{extra}#
"""
_TEMPLATES = {
    "constrained": ("#   Hold everything you measured that is long-lived.  Leave a measured product such as HCHO out\n"
                    "#   to compare the model with it afterwards.\n",
                    "time,lat,lon,T,p,RH,jNO2,C_NO [ppb],C_NO2 [ppb],C_O3 [ppb],C_CO [ppb],C_C5H8 [ppt]\n"
                    "2021-08-15 18:00,34.14,-118.13,296.0,985.0,70.0,0.0030,1.2,8.5,32.0,180.0,250\n"
                    "2021-08-15 19:00,34.14,-118.13,298.5,985.0,62.0,0.0062,0.9,6.0,48.0,170.0,520\n"
                    "2021-08-15 20:00,34.14,-118.13,301.0,984.5,55.0,0.0078,0.7,4.5,66.0,160.0,780\n"),
    "plume": ("#   Each row is a time since the start of the plume.  kdil often comes from the decay of CO.\n",
              "time,lat,lon,T,p,RH,kdil,jNO2,B_O3 [ppb],B_CO [ppb]\n"
              "2021-08-15 20:00,34.14,-118.13,298.0,950.0,45.0,1.2e-4,0.0080,40,120\n"
              "2021-08-15 20:15,34.15,-118.05,298.5,950.0,44.0,1.0e-4,0.0081,40,120\n"
              "2021-08-15 20:30,34.16,-117.97,299.0,950.0,43.0,0.8e-4,0.0082,40,120\n"),
    "lagrangian": ("#   Each row is one point of the trajectory.  Set the footprint area in the case when you use mol/s.\n",
                   "time,lat,lon,T,p,RH,h,jNO2,E_NO [mol/s],E_NO2 [mol/s],E_CO [mol/s],E_C5H8 [mol/km2/s],B_O3 [ppb],B_CO [ppb]\n"
                   "2021-08-15 17:00,33.95,-118.43,294.0,1005.0,75.0,450,0.0065,3.1,0.4,28.0,0.002,30,160\n"
                   "2021-08-15 17:15,33.98,-118.38,295.0,1004.0,73.0,520,0.0069,4.0,0.5,35.0,0.003,30,160\n"
                   "2021-08-15 17:30,34.01,-118.33,296.0,1003.0,71.0,600,0.0072,5.2,0.6,44.0,0.004,30,160\n"),
}


def write_template(path: str, run_type: str) -> str:
    if run_type not in _TEMPLATES:
        raise ValueError(f"unknown kind of run '{run_type}'")
    extra, body = _TEMPLATES[run_type]
    head = _TEMPLATE_HEAD.format(title=RUN_TYPES[run_type], units=", ".join(EMISSION_UNITS), extra=extra)
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write(head + body)
    return path


def selftest() -> int:
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        for rt in RUN_TYPES:
            p = write_template(os.path.join(d, rt + ".csv"), rt)
            spec = FieldSpec(name=rt, run_type=rt, table=p, area_km2=64.0, chemistry={"mechanism_stem": "x"},
                             dilution=DEFAULT_DILUTION[rt])
            _t, rep = check(spec)
            print(rt, "ok" if rep.ok else rep.text())
            if not rep.ok:
                return 1
    return 0


if __name__ == "__main__":
    sys.exit(selftest() if "--selftest" in sys.argv else 0)
