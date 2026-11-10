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
# fields always written out (with the allowed values as a comment), even when they have the default
_SHOW_CONDITIONS = ("date", "start", "end", "step_s", "source", "temperature_K", "temperature_mode",
                    "pressure_hPa", "h2o_percent", "flow_m3h", "roof_open", "j_source")
_SHOW_SPEC = ("title", "base_mechanism", "addons", "precursors", "initials")
_SHOW_FIELD = ("run_type", "table", "date", "chemistry", "dilution", "step_s")
_COMMENTS = {
    "source": "manual | archive (SAPHIR data archive) | file (CSV/Excel with T, p, H2O, J ...)",
    "temperature_mode": "enz (T from the ENZ column) | const (temperature_K)",
    "j_source": "measured (ENZ J columns) | clear_sky (clear-sky J of the mechanism)",
    "roof_open": "[(open HH:MM, close HH:MM), ...]; lamp chambers: lights on/off",
    "base_mechanism": "juice.mechanisms()['bases']",
    "addons": "juice.mechanisms()['addons']",
    "precursors": "MCM names; the mechanism subset is cut from the full MCM",
    "initials": "[(species, expression inside CONST), ...], e.g. ('NO2', '20e-9*M')",
    "run_type": " | ".join(fs.RUN_TYPES),
    "dilution": " | ".join(fs.DILUTION_SOURCES),
}


def _diff(obj, base) -> List[str]:
    """Field names of a dataclass whose value differs from `base`."""
    return [f.name for f in fields(obj) if getattr(obj, f.name) != getattr(base, f.name)]


def _line(target: str, name: str, value, show_comment: bool = True) -> str:
    text = f"{target}.{name} = {value!r}"
    cmt = _COMMENTS.get(name) if show_comment else None
    return text + (f"    # {cmt}" if cmt else "")


def _injection_code(ev: Injection) -> str:
    base = Injection(name=ev.name, products=ev.products)
    args = [f"name={ev.name!r}", f"products={ev.products!r}"] + [f"{k}={getattr(ev, k)!r}" for k in _diff(ev, base)]
    return "    juice.Injection(" + ",\n                    ".join(args) + "),"


def experiment_code(spec: ExperimentSpec) -> List[str]:
    """The body of build(): code that rebuilds `spec` exactly."""
    chamber = spec.conditions.chamber_name
    base = new_experiment(spec.name, chamber)
    L = [f"exp = juice.new_experiment({spec.name!r}, chamber={chamber!r})    # chamber: juice.chambers()",
         "", "# when and under which conditions", "c = exp.conditions"]
    for k in _SHOW_CONDITIONS:
        L.append(_line("c", k, getattr(spec.conditions, k)))
    for k in _diff(spec.conditions, base.conditions):
        if k not in _SHOW_CONDITIONS:
            L.append(_line("c", k, getattr(spec.conditions, k)))
    L += ["", "# chemistry"]
    for k in _SHOW_SPEC:
        L.append(_line("exp", k, getattr(spec, k)))
    skip = set(_SHOW_SPEC) | {"name", "conditions", "chamber", "injections"}
    for k in _diff(spec, base):
        if k not in skip:
            L.append(_line("exp", k, getattr(spec, k)))
    L += ["", "# injections: one entry per event (total_ppb = the whole amount, start HH:MM, duration_min)",
          "exp.injections = ["] + [_injection_code(ev) for ev in spec.injections] + ["]"]
    ch = _diff(spec.chamber, base.chamber)
    L += ["", f"# chamber effects: the {chamber} defaults" + (", except:" if ch else " (exp.chamber.<name> to change one)")]
    for k in ch:
        L.append(_line("exp.chamber", k, getattr(spec.chamber, k), False))
    return L


def field_code(spec: FieldSpec) -> List[str]:
    base = new_field_case(spec.name, spec.run_type, spec.table)
    L = [f"case = juice.new_field_case({spec.name!r}, run_type={spec.run_type!r},",
         f"                            table={spec.table!r})    # run_type: {_COMMENTS['run_type']}", ""]
    for k in _SHOW_FIELD:
        if k not in ("run_type", "table"):
            L.append(_line("case", k, getattr(spec, k)))
    for k in _diff(spec, base):
        if k not in _SHOW_FIELD and k != "name":
            L.append(_line("case", k, getattr(spec, k)))
    return L


_SCRIPT = '''"""
{name} — this JUICE {what} as code.

Written by JUICE {version} on {when} from {origin}.  Running this file
rebuilds the experiment, writes its input files into OUT, runs the solver and prints a summary:

    python {fname}

Change anything below and run it again: this is the same experiment the window holds, field by field.
juice.chambers(), juice.mechanisms(), juice.solvers() and juice.OPTIONS list what can be chosen.
"""
import os
import sys

sys.path.insert(0, {here!r})          # the JUICE folder
import juice

OUT = os.path.dirname(os.path.abspath(__file__))   # where the input files and results go
SOLVER = {solver!r}       # {solvers}  (None = the settings file's profile)
TOML = {toml!r}         # settings file (paths, solver numbers); '' = JUICE's default
KV_PD = {kv_pd!r}         # also export reaction rates (Process Analysis); False is faster


def build():
{body}
    return {var}


if __name__ == "__main__":
    {var} = build()
    for level, message in juice.check({var}, OUT, TOML):
        if level != "info":
            print(level + ":", message)
    run = juice.run({var}, OUT, solver=SOLVER, toml=TOML, kv_pd=KV_PD, script=False)
    print(run.summary())
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
                          fname=f"{spec.name}_run.py", here=HERE, solver=solver,
                          solvers=" | ".join(sp.SOLVER_PROFILES), toml=toml, kv_pd=bool(kv_pd),
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
