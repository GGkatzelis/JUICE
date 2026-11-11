"""
juice.py — JUICE from code
==========================

Everything the JUICE windows do, as plain Python calls.  The windows themselves use the same
functions, so an experiment built here and one built by clicking produce the same files.

    import juice

    exp = juice.new_experiment("isoprene_day", chamber="SAPHIR")
    c = exp.conditions
    c.date, c.start, c.end = "2024-07-09", "06:00", "18:00"
    c.source = "manual"                      # manual | archive | file
    c.roof_open = [("08:00", "16:00")]
    exp.precursors = ["C5H8"]                # MCM names: the mechanism subset is cut from the full MCM
    exp.injections = [juice.Injection(name="NO", products="NO", total_ppb="5", start="07:30"),
                      juice.Injection(name="isoprene", products="C5H8", total_ppb="10", start="08:30")]

    run = juice.run(exp, r"D:\\runs\\isoprene_day", solver="fast")   # write -> prepare -> solve
    print(run.summary())
    o3 = run.results().ppb("O3")

The four steps can also be called one by one (what juice.run does):

    files   = juice.write(exp, out_dir)                    # <name>.enz, <name>.easy, chemistry/, <name>.experiment.json
    json    = juice.prepare(exp, out_dir, solver="fast")   # + EASY -> ELOG -> <name>.json (the solver input)
    cs      = juice.solve(json)                            # scenario_runner in its own process -> <name>.cs
    results = juice.results(cs)                            # the .cs as a table

What to choose from:

    juice.chambers()     chamber profiles (SAPHIR, generic_lamps, ...)
    juice.mechanisms()   base mechanisms and add-ons
    juice.solvers()      solver profiles with what each one costs
    juice.templates()    the experiment templates of the Start step
    juice.OPTIONS        the allowed values of the fields that take a keyword

Field cases (site, plume, Lagrangian parcel) work the same way with juice.new_field_case().

Every experiment the windows prepare also gets a <name>_run.py next to its files (juice.write_script):
the same experiment written out as code, which runs it from start to end when executed.
"""
from __future__ import annotations

import datetime as _dt
import os
import subprocess
import sys
from dataclasses import fields, is_dataclass
from typing import Dict, List, Optional, Tuple, Union

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

import experiment_spec as es                                    # noqa: E402
import field_spec as fs                                         # noqa: E402
import mechanism_registry as mreg                               # noqa: E402
import scenario_prep as sp                                      # noqa: E402

from experiment_spec import (ExperimentSpec, Conditions, ChamberEffects,  # noqa: E402,F401
                             InjectionEvent as Injection)
from field_spec import FieldSpec                                # noqa: E402,F401

Spec = Union[ExperimentSpec, FieldSpec]

# the fields that take one of a few keywords: shown in the generated scripts and in the tutorial
OPTIONS: Dict[str, Tuple[str, ...]] = {
    "solver": tuple(sp.SOLVER_PROFILES),
    "conditions.source": ("manual", "archive", "file"),
    "conditions.temperature_mode": ("enz", "const"),
    "conditions.j_source": ("measured", "clear_sky"),
    "gecko_core": ("gecko", "mcm"),
    "field.run_type": tuple(fs.RUN_TYPES),
    "field.dilution": tuple(fs.DILUTION_SOURCES),
}


# ── what to choose from ───────────────────────────────────────────────────────────────────────────
def chambers() -> Dict[str, str]:
    """Chamber profile id -> its name (chambers/*.chamber.json)."""
    import chambers as ch
    return {k: (getattr(p, "name", "") or k) for k, p in ch.load_dir().items()}


def mechanisms() -> Dict[str, Dict[str, str]]:
    """{'bases': {id: name}, 'addons': {id: name}} of the mechanism registry."""
    return {"bases": {k: getattr(b, "name", k) for k, b in mreg.load_bases().items()},
            "addons": {k: getattr(a, "name", k) for k, a in mreg.load_addons().items()}}


def solvers(toml: str = "") -> Dict[str, str]:
    """Solver profile -> one-line summary (with the numbers of the settings file applied)."""
    return {k: str(v.get("summary", "")) for k, v in sp.solver_profiles(toml or default_toml()).items()}


def templates() -> Dict[str, str]:
    """Template name -> path (templates/*.experiment.json)."""
    return es.list_templates()


def default_toml() -> str:
    """The settings file the windows use: GG_GUI.toml, else FR_GUI.toml next to JUICE ('' if none)."""
    return sp.pick_toml([])


def paths(toml: str = "") -> Dict[str, str]:
    """[paths] of the settings file: include_dir (the include files), output_root, easy_root."""
    return sp.load_paths(toml or default_toml())


# ── building an experiment ───────────────────────────────────────────────────────────────────────
def new_experiment(name: str = "experiment", chamber: str = "SAPHIR", template: str = "") -> ExperimentSpec:
    """A chamber experiment with the defaults of `chamber` (what picking the chamber on the Start step
    does), or a copy of a template (templates() lists them; the template's own chamber is kept)."""
    if template:
        path = templates().get(template, template)
        _t, _d, spec = es.load_template(path)
        spec.name = name
        return spec
    spec = ExperimentSpec(name=name)
    es.apply_chamber(spec, chamber)
    return spec


def new_field_case(name: str = "field_case", run_type: str = "constrained", table: str = "") -> FieldSpec:
    """A field case of one kind of run (constrained | plume | lagrangian) reading `table`, with the
    dilution that kind of run uses by default (what the first step of the field window does)."""
    if run_type not in fs.RUN_TYPES:
        raise ValueError(f"run_type must be one of {', '.join(fs.RUN_TYPES)}")
    return FieldSpec(name=name, run_type=run_type, table=os.path.abspath(table) if table else "",
                     dilution=fs.DEFAULT_DILUTION[run_type])


def load(path: str) -> Spec:
    """An experiment (.experiment.json) or field case (.field.json) saved by JUICE."""
    return fs.load_spec(path) if path.lower().endswith(fs.SPEC_SUFFIX) else es.load_spec(path)


def save(spec: Spec, path: str) -> str:
    return fs.save_spec(spec, path) if isinstance(spec, FieldSpec) else es.save_spec(spec, path)


def check(spec: Spec, out_dir: str = "", toml: str = "") -> List[Tuple[str, str]]:
    """[(level, message)], level 'error' | 'warning' | 'info': the checks of the review panel.
    juice.write refuses an experiment with errors."""
    if isinstance(spec, FieldSpec):
        _table, rep = fs.check(spec)
        return ([("error", m) for m in rep.errors] + [("warning", m) for m in rep.warnings]
                + [("info", m) for m in rep.info])
    known = es.known_species_for(spec, paths(toml)["include_dir"], out_dir)
    return list(es.validate(spec, known))


# ── the four steps ───────────────────────────────────────────────────────────────────────────────
def write(spec: Spec, out_dir: str, toml: str = "") -> dict:
    """Write the input files: <name>.enz, <name>.easy, the mechanism subset in chemistry/ and the saved
    spec (<name>.experiment.json or .field.json).  Returns {'enz', 'easy', 'spec_path', 'notes', ...}."""
    out_dir = os.path.abspath(out_dir)
    if isinstance(spec, FieldSpec):
        import field_scenario as fsc
        table, rep = fs.check(spec)
        if not rep.ok or table is None:
            raise ValueError("the field case has errors:\n" + rep.text())
        return fsc.write_case(spec, table, out_dir)
    errors = [m for lvl, m in check(spec, out_dir, toml) if lvl == "error"]
    if errors:
        raise ValueError("the experiment has errors:\n  " + "\n  ".join(errors))
    res = es.write_experiment(spec, out_dir)
    res["spec_path"] = es.save_spec(spec, os.path.join(out_dir, spec.name + ".experiment.json"))
    return res


def prepare(spec: Spec, out_dir: str, solver: Optional[str] = None, toml: str = "", kv_pd: bool = True,
            files: Optional[dict] = None) -> str:
    """Write the files (unless `files` from write() is given) and turn them into the solver input
    <name>.json: EASY -> ELOG -> scenario JSON, with the solver profile `solver` (None = the settings
    file's).  kv_pd: also export reaction rates and the k/v tables (slower, needed for Process Analysis).
    Returns the JSON path."""
    toml = toml or default_toml()
    res = files or write(spec, out_dir, toml)
    if isinstance(spec, FieldSpec):
        extra, version = list((spec.chemistry or {}).get("include_dirs") or []), ""
    else:
        extra, version = list(spec.include_dirs), (spec.easy_version or None)
    return sp.prepare_scenario(res, paths(toml)["include_dir"], kv_pd, toml, extra_dirs=extra,
                               easy_version=version, solver=sp.solver_from_toml(toml, solver))


def solve(json_path: str, toml: str = "", log: Optional[str] = None) -> str:
    """Run the solver (scenario_runner.py, in its own process, as the Run step does) on a prepared
    <name>.json.  log: a file for the solver output (None = <name>_run.log next to the JSON,
    '' = the console).  Returns the path of the results file <name>.cs; raises if the run fails."""
    toml = toml or default_toml()
    args = [sys.executable, "-u", os.path.join(HERE, "scenario_runner.py"), json_path] + (["--toml", toml] if toml else [])
    env = {**os.environ, "PYTHONIOENCODING": "utf-8", "PYTHONUTF8": "1"}
    if log is None:
        log = json_path[:-5] + "_run.log"
    if log:
        with open(log, "w", encoding="utf-8") as fh:
            rc = subprocess.run(args, stdout=fh, stderr=subprocess.STDOUT, env=env, cwd=HERE).returncode
    else:
        rc = subprocess.run(args, env=env, cwd=HERE).returncode
    cs = json_path[:-5] + ".cs"
    if rc != 0 or not os.path.isfile(cs):
        raise RuntimeError(f"the solver stopped (exit code {rc})" + (f"; see {log}" if log else ""))
    return cs


class Results:
    """A results file (.cs) as a table: .df has TIME (s since 2000-01-01 UTC) and one column per species
    in molecules cm-3; .utc and .hour give the time axis; .ppb(name) converts with the air density."""

    def __init__(self, cs_path: str, pressure_hPa: float = 1013.25, temperature_K: float = 298.0):
        from run_benchmark import parse_cs_file
        self.path = cs_path
        self.df = parse_cs_file(cs_path)
        self.M = pressure_hPa * 100.0 / (1.380649e-23 * temperature_K) * 1e-6     # molecules cm-3

    @property
    def species(self) -> List[str]:
        return [c for c in self.df.columns if c != "TIME"]

    @property
    def utc(self):
        epoch = _dt.datetime(2000, 1, 1, tzinfo=_dt.timezone.utc)
        return [epoch + _dt.timedelta(seconds=float(t)) for t in self.df["TIME"]]

    @property
    def hour(self):
        """Hours since midnight UTC of the first output time."""
        t = self.df["TIME"].to_numpy(float)
        return (t - (t[0] // 86400.0) * 86400.0) / 3600.0

    def __getitem__(self, name: str):
        return self.df[name].to_numpy(float)

    def ppb(self, name: str):
        """Mixing ratio in ppb, with the air density of the pressure and temperature given (constant)."""
        return self[name] / self.M * 1e9


def results(cs_path: str, spec: Optional[Spec] = None) -> Results:
    """Read a results file; with the experiment given, ppb() uses its pressure and temperature."""
    if isinstance(spec, ExperimentSpec):
        return Results(cs_path, spec.conditions.pressure_hPa, spec.conditions.temperature_K)
    return Results(cs_path)


class Run:
    """What juice.run returns: the files of one run and its results."""

    def __init__(self, spec: Spec, out_dir: str, files: dict, json_path: str, cs_path: str, seconds: float):
        self.spec, self.out_dir, self.files, self.json, self.cs, self.seconds = spec, out_dir, files, json_path, cs_path, seconds

    def results(self) -> Results:
        return results(self.cs, self.spec)

    def summary(self) -> str:
        r = self.results()
        lines = [f"{self.spec.name}: {len(r.df)} output times, {len(r.species)} species, {self.seconds:.0f} s",
                 f"  results  {self.cs}", f"  inputs   {os.path.dirname(self.json)}"]
        for n in ("O3", "NO", "NO2", "OH", "HO2"):
            if n in r.df.columns:
                v = r[n]
                lines.append(f"  {n:4s} start {v[0]:.3e}  end {v[-1]:.3e}  max {v.max():.3e} molec/cm3")
        return "\n".join(lines)


def run(spec: Spec, out_dir: str, solver: Optional[str] = None, toml: str = "", kv_pd: bool = True,
        log: Optional[str] = None, script: bool = True) -> Run:
    """write -> prepare -> solve.  script: also write <name>_run.py (this experiment as code).
    Returns a Run (files, JSON, results file, run time)."""
    import time
    toml = toml or default_toml()
    t0 = time.time()
    files = write(spec, out_dir, toml)
    if script:
        write_script(spec, os.path.join(os.path.abspath(out_dir), spec.name + "_run.py"), solver=solver, toml=toml, kv_pd=kv_pd)
    js = prepare(spec, out_dir, solver, toml, kv_pd, files=files)
    cs = solve(js, toml, log)
    return Run(spec, os.path.abspath(out_dir), files, js, cs, time.time() - t0)


# ── an experiment as code ────────────────────────────────────────────────────────────────────────
# <name>_run.py follows the steps of the window. Every field of the experiment is written out with what it
# means (and the allowed values), so the file documents the experiment completely and rebuilds it exactly.
# A layout entry is ("# text", None) for a comment line, ("", None) for an empty line, or (field, comment).

_RULE = "# " + "═" * 100

_CHAMBER_LAYOUT = {
    "conditions": [
        ("# Which day and window (UTC); with typed conditions the date only sets the sun's position.", None),
        ("date", "YYYY-MM-DD, UTC"),
        ("start", "HH:MM UTC"),
        ("end", "HH:MM UTC"),
        ("step_s", "output step in seconds: one ENZ row and one result row per step"),
        ("", None),
        ("# Where the physical conditions come from. 'file': a CSV/Excel table you prepare (columns T, p, H2O or RH,",
         None),
        ("#   flow, roof, jNO2/jO1D/..., NO/NO2/O3; missing ones fall back to the constants below). An empty, documented",
         None),
        ("#   table: import measured_inputs; measured_inputs.write_template('my_conditions.csv'). An ENZ file works too.",
         None),
        ("source", "manual | file | archive (the SAPHIR data archive of that day)"),
        ("input_file", "the table for source 'file' (full path)"),
        ("archive_root", "the archive drive for source 'archive'"),
        ("", None),
        ("# Constants: used for anything the file or archive does not give.", None),
        ("temperature_K", "K"),
        ("temperature_mode", "const = temperature_K | enz = the T column of the file/archive"),
        ("pressure_hPa", "hPa"),
        ("h2o_percent", "water vapour, vol-%"),
        ("flow_m3h", "replenishment flow, m3/h: sets the dilution"),
        ("", None),
        ("# Light: roof open (lamp chamber: lights on) as [(open, close), ...]; a 'roof' column of the file wins.",
         None),
        ("roof_open", "[(HH:MM, HH:MM), ...]"),
        ("j_source", "clear_sky = the mechanism's clear-sky J | measured = jNO2, jO1D, ... of the file/archive"),
        ("", None),
        ("# Held to measurements at every step (NO, NO2, O3; needs source 'file' or 'archive' with those columns).",
         None),
        ("constrain_species", "e.g. ['NO2', 'O3']"),
        ("chamber_name", "the chamber profile (set by new_experiment above)"),
    ],
    "chemistry": [
        ("base_mechanism", "juice.mechanisms()['bases']: MCM_v3.3.1 | CRACMM2 | CRACMM1_aq | RACM2 | a GECKO-A base"),
        ("precursors", "MCM names of the compounds; the MCM subset is cut from them (empty for a lumped base)"),
        ("addons", "juice.mechanisms()['addons'], e.g. ['NOAA_BB']"),
        ("inorganic_only", "True: inorganic chemistry and chamber sources only, no VOC"),
        ("mechanism_stem", "instead of precursors: an existing EASY include set (with include_dirs)"),
        ("include_dirs", "folders searched for #include files (your own chemistry, an imported experiment's)"),
        ("extra_includes", "your own EASY include files, loaded after the mechanism"),
        ("gecko_core", "GECKO-A base only: gecko (its own inorganic core) | mcm (MCM core + measured-J include)"),
        ("gecko_partitioning", "GECKO-A base only: load the gas-particle-wall partitioning include"),
        ("", None),
        ("# Initial values: everything starts at initial_all molecules cm-3; initials are expressions inside CONST()",
         None),
        ("#   (M = air density, so '20e-9*M' is 20 ppb).", None),
        ("initial_all", "molecules cm-3 for every species not listed in initials"),
        ("initials", "[(species, expression), ...], e.g. [('O3', '30e-9*M')]"),
        ("", None),
        ("# Rate changes: (reaction as the mechanism writes it, new rate in EASY notation or '0' to remove, comment).",
         None),
        ("reaction_adjustments", "e.g. [('OH + NO2 --> HNO3', '1.1*KMT08', 'test: +10 %')]"),
    ],
    "chamber": [
        ("# Each process: a switch and its EASY expressions (they use the ENZ inputs: HV_in = roof, T_in, jno2, FLOW).",
         None),
        ("hono", "HONO wall source under light"),
        ("hono_lhs", None), ("hono_rhs", None), ("hono_const", "strength"), ("hono_rate", None),
        ("hcho", "HCHO wall source under light"),
        ("hcho_lhs", None), ("hcho_rhs", None), ("hcho_const_name", None), ("hcho_const", "strength"), ("hcho_rate", None),
        ("bgreac", "background OH reactivity 'Y'"),
        ("y_name", None), ("y_product", "what Y + OH makes"), ("y_ppb", "CO-equivalents, ppb"),
        ("y_expr", "a full expression instead of y_ppb"), ("y_rate", None),
        ("dilution", "dilution by the replenishment flow"),
        ("dilute_expr", "the dilution rate, 1/s"), ("dil_value", None), ("dil_flow_const", "a constant flow instead of the ENZ column"),
        ("tracer_dilution", "tracers that only dilute"),
        ("walldep", "wall deposition of HNO3, HONO, H2O2, HO2, O3, N2O5 (saphir.walldep.incl)"),
        ("walldep_value", "factor on it"),
        ("dark_hono", "HONO source in the dark"),
        ("dark_hono_rhs", None), ("dark_hono_rate", None),
        ("light_source", "NO2 source under light"),
        ("light_source_rhs", None), ("light_source_rate", None), ("light_source_initial", None),
        ("wall_loss", "first-order wall losses of the chamber profile"),
        ("wall_loss_s1", "species -> k, 1/s"),
    ],
    "easy": [
        ("# Settings of the EASY file. Change only if you know why.", None),
        ("hmax", "longest internal time step, s"),
        ("outstep", None),
        ("mech_incl_suffix", None),
        ("use_anorg", "load the MCM inorganic include set"),
        ("use_j_clear_sky_incl", "load j_clear_sky.incl (PHTFKT)"),
        ("overrides", "EASY statements written verbatim after everything else (last definition wins)"),
        ("", None),
        ("# Provenance: what the experiment came from (kept so that this file rebuilds it exactly).", None),
        ("template", None), ("easy_version", "EASY of the FACSIMILE run it came from"),
        ("reference_cs", "FACSIMILE results to compare with"), ("import_notes", None),
    ],
}

_FIELD_LAYOUT = {
    "case": [
        ("# One input table, one row per time (UTC): T, p, RH or H2O, h, kdil, jNO2..., and species columns", None),
        ("#   E_<sp> emission, B_<sp> background, C_<sp> held, O_<sp> observed only (units in [brackets]).", None),
        ("#   An empty table with the columns of this kind of run: field_spec.write_template(path, run_type).", None),
        ("description", None),
        ("date", "YYYY-MM-DD, needed when the table gives HH:MM or hours of day"),
        ("site", None),
        ("lat", "degrees, when the table has no lat column"),
        ("lon", "degrees, when the table has no lon column"),
        ("step_s", "output step in seconds"),
        ("spin_up_days", "constrained site: repeat the first day this many times before the output"),
        ("generic_lifetime_h", "constrained site: first-order loss of everything neither held nor deposited"),
    ],
    "chemistry": [
        ("# The Chemistry step: base_mechanism, precursors, addons, inorganic_only, mechanism_stem, include_dirs,", None),
        ("#   reaction_adjustments (the same meaning as for a chamber experiment).", None),
        ("chemistry", None),
    ],
    "sources": [
        ("default_emission_unit", "unit of E_ columns without [unit]: ppb/s, molec/cm3/s, molec/cm2/s, mol/km2/s, mol/km2/h, mol/s, molec/s"),
        ("area_km2", "footprint A of E/(h*A); needed by mol/s and molec/s"),
        ("emission_scale", "factor on the E_ columns: {species | 'NOX' | 'VOC' | 'ALL': factor}"),
        ("initials", "ppb at the start: {species: ppb}"),
    ],
    "mixing": [
        ("dilution", " | ".join(fs.DILUTION_SOURCES) + " (pbl: growth of h; table: kdil column; tracer: decay of a held tracer)"),
        ("k_dil_per_h", "dilution 'constant', 1/h"),
        ("k_collapse_per_h", "dilution 'pbl': rate while the layer does not grow, 1/h"),
        ("tracer", "dilution 'tracer': a species with a C_ column"),
        ("deposition_velocity_cm_s", "{species: cm/s}; loss = v/h"),
        ("deposition_lifetime_h", "{species: h}; loss = 1/lifetime"),
    ],
}


def _diff(obj, base) -> List[str]:
    """Field names of a dataclass whose value differs from `base`."""
    return [f.name for f in fields(obj) if getattr(obj, f.name) != getattr(base, f.name)]


def _emit(target: str, obj, layout, done: set) -> List[str]:
    L = []
    for name, cmt in layout:
        if name.startswith("#") or name == "":
            L.append(name); continue
        if not hasattr(obj, name):
            continue
        done.add(name)
        text = f"{target}.{name} = {getattr(obj, name)!r}"
        L.append(text + (f"    # {cmt}" if cmt else ""))
    return L


def _rest(target: str, obj, done: set) -> List[str]:
    """Fields the layout does not name (a newer JUICE may have more): written so the rebuild stays exact."""
    return [f"{target}.{f.name} = {getattr(obj, f.name)!r}" for f in fields(obj) if f.name not in done
            and f.name not in ("name", "conditions", "chamber", "injections")]


def _step(n: int, title: str, what: str) -> List[str]:
    return ["", _RULE, f"# STEP {n} — {title}", f"# {what}", _RULE]


def _injection_code(ev: Injection) -> List[str]:
    base = Injection(name=ev.name, products=ev.products)
    args = [f"name={ev.name!r}", f"products={ev.products!r}"] + [f"{k}={getattr(ev, k)!r}" for k in _diff(ev, base)]
    return ["    juice.Injection(" + ",\n                    ".join(args) + "),"]


def experiment_code(spec: ExperimentSpec) -> List[str]:
    """The body of build(): code that rebuilds `spec` exactly, step by step as in the window."""
    chamber = spec.conditions.chamber_name
    L = _step(1, "EXPERIMENT: which chamber, when, and where the conditions come from",
              "The chamber profile sets the volume, the dilution and the chamber processes; juice.chambers() lists them.")
    L += [f"exp = juice.new_experiment({spec.name!r}, chamber={chamber!r})",
          f"exp.title = {spec.title!r}", "", "c = exp.conditions"]
    done_c: set = set()
    L += _emit("c", spec.conditions, _CHAMBER_LAYOUT["conditions"], done_c) + _rest("c", spec.conditions, done_c)

    L += _step(2, "CHEMISTRY: which mechanism, which compounds",
               "juice.mechanisms() lists the bases and the add-ons.")
    done_s = {"title"}
    L += _emit("exp", spec, _CHAMBER_LAYOUT["chemistry"], done_s)

    L += _step(3, "INJECTIONS: what goes into the chamber, how much, when",
               "One Injection per event: total_ppb = the whole amount, spread over duration_min from start (HH:MM UTC).")
    L += ["# products: a species or a split ('0.9*X + 0.1*Y'); empty = recorded in the ENZ only. Further arguments:",
          "#   column / switch (ENZ column, EASY switch), extra_channels=[(products, ppb)] impurities,",
          "#   chain=[(lhs, rhs, rate)] a delayed release, comment, action, delivery (how the amount was derived).",
          "# From a liquid volume or a flow: injection_calc.ppb_liquid / ppb_mfc / ppb_gas_syringe.",
          "exp.injections = ["]
    for ev in spec.injections:
        L += _injection_code(ev)
    L += ["]"]

    L += _step(4, "CHAMBER: what the chamber adds on its own",
               f"The {chamber} profile set these; switch a process off or change its expressions here.")
    L += ["ch = exp.chamber"]
    done_ch: set = set()
    L += _emit("ch", spec.chamber, _CHAMBER_LAYOUT["chamber"], done_ch) + _rest("ch", spec.chamber, done_ch)

    L += ["", _RULE, "# EASY SETTINGS AND PROVENANCE", _RULE]
    L += _emit("exp", spec, _CHAMBER_LAYOUT["easy"], done_s) + _rest("exp", spec, done_s)
    return L


def field_code(spec: FieldSpec) -> List[str]:
    """The body of build() for a field case, step by step as in the field window."""
    L = _step(1, "FIELD CASE: which kind of run, which table",
              f"run_type: {' | '.join(fs.RUN_TYPES)} (constrained site, plume or transect, Lagrangian parcel).")
    L += [f"case = juice.new_field_case({spec.name!r}, run_type={spec.run_type!r},",
          f"                            table={spec.table!r})"]
    done: set = {"name", "run_type", "table"}
    L += _emit("case", spec, _FIELD_LAYOUT["case"], done)
    L += _step(2, "CHEMISTRY", "Any base mechanism: the MCM subset, a lumped mechanism or a GECKO-A base.")
    L += _emit("case", spec, _FIELD_LAYOUT["chemistry"], done)
    L += _step(3, "EMISSIONS AND CONSTRAINTS", "Emissions are the E_ columns, held species the C_ columns of the table.")
    L += _emit("case", spec, _FIELD_LAYOUT["sources"], done)
    L += _step(4, "MIXING AND DEPOSITION", "How the air mass is diluted, and what deposits.")
    L += _emit("case", spec, _FIELD_LAYOUT["mixing"], done)
    rest = _rest("case", spec, done)
    if rest:
        L += ["", "# other fields"] + rest
    return L


_SCRIPT = '''"""
{name} — this JUICE {what} as code.

Written by JUICE {version} on {when} from {origin}. It follows the steps of the window; every field
is written out with what it means. Running this file rebuilds the {what}, writes its input files into
OUT, runs the solver and prints a summary:

    python {fname}

Change anything in build() and run it again. juice.chambers(), juice.mechanisms(), juice.solvers() and
juice.OPTIONS list what can be chosen; docs/juice_from_code.md explains the rest. Note: running it in
this folder overwrites the run that is here; set OUT to another folder to keep both.
"""
import os
import sys

sys.path.insert(0, {here!r})          # the JUICE folder
import juice                                                           # noqa: E402

OUT = os.path.dirname(os.path.abspath(__file__))   # where the input files and results go
TOML = {toml!r}         # settings file (paths, solver numbers); '' = JUICE's default

# The run files in OUT: <name>.enz (time table), <name>.easy (FACSIMILE master file), chemistry/ (mechanism
# include files), <name>.{suffix} (this {what} as data), <name>.json (solver input), <name>.cs (results),
# <name>_run.log (solver output). The .enz and .easy run in FACSIMILE unchanged.


def build():
{body}
    return {var}


def main():
    {var} = build()

    {rule}
    # STEP 5 — RUN
    {rule}
    # The review panel's checks: errors stop the run, warnings and notes are worth reading.
    for level, message in juice.check({var}, OUT, TOML):
        if level != "info":
            print(level + ":", message)

    # The solver (juice.solvers() says what each costs):
    #   "fast"     BDF with the exact Jacobian (the default)
    #   "precise"  the same at a tighter tolerance: the closest agreement with FACSIMILE
    #   "twostep"  Jacobian-free: very large mechanisms (GECKO-A) and the complete MCM
    #   None       the profile of the settings file
    solver = {solver!r}
    kv_pd = {kv_pd!r}        # also write reaction rates (.kv) and production/destruction (.pd) for Process Analysis

    # The four steps, one by one, so that you can stop in between:
    files = juice.write({var}, OUT, TOML)                                    # 1  input files
    #   <-- here you may edit the .easy or .enz by hand; the next step reads them from disk
    json_path = juice.prepare({var}, OUT, solver=solver, toml=TOML, kv_pd=kv_pd, files=files)   # 2  solver input
    cs_path = juice.solve(json_path, TOML)                                   # 3  the solver (log: <name>_run.log)
    r = juice.results(cs_path, {var})                                         # 4  results as a table

    # r["OH"] molecules cm-3, r.ppb("O3"), r.hour (hours of the day, UTC), r.species, r.df (pandas).
    # The viewer: python results_viewer.py <cs_path>
    print(f"{{len(r.df)}} output times, {{len(r.species)}} species -> {{cs_path}}")
    for n in ("O3", "NO", "NO2", "OH", "HO2"):
        if n in r.species:
            print(f"  {{n:4s}} start {{r[n][0]:.3e}}  end {{r[n][-1]:.3e}}  max {{r[n].max():.3e}} molec/cm3")


if __name__ == "__main__":
    main()
'''


def script_text(spec: Spec, solver: Optional[str] = None, toml: str = "", kv_pd: bool = True, origin: str = "the JUICE window") -> str:
    try:
        from version import __version__ as ver
    except Exception:
        ver = ""
    field = isinstance(spec, FieldSpec)
    body = field_code(spec) if field else experiment_code(spec)
    return _SCRIPT.format(name=spec.name, what="field case" if field else "experiment", version=ver,
                          when=_dt.datetime.now().strftime("%Y-%m-%d %H:%M"), origin=origin,
                          fname=f"{spec.name}_run.py", here=HERE, solver=solver, toml=toml, kv_pd=bool(kv_pd),
                          suffix="field.json" if field else "experiment.json", rule=_RULE,
                          body="\n".join(("    " + l) if l else "" for l in body),
                          var="case" if field else "exp")


def write_script(spec: Spec, path: str, solver: Optional[str] = None, toml: str = "", kv_pd: bool = True,
                 origin: str = "the JUICE window") -> str:
    """Write `spec` as a runnable Python file (<name>_run.py).  build() in that file returns a spec equal
    to `spec`, field by field; running the file writes, prepares and solves it."""
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(script_text(spec, solver, toml, kv_pd, origin))
    return path


def build_from_script(path: str) -> Spec:
    """build() of a <name>_run.py, without running it."""
    import importlib.util
    s = importlib.util.spec_from_file_location("_juice_script", path)
    mod = importlib.util.module_from_spec(s)
    s.loader.exec_module(mod)
    return mod.build()


if __name__ == "__main__":
    print(__doc__)
