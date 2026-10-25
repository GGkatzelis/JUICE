"""
study.py — many runs from existing experiments (no Qt).

A *study* = base experiments × one variant × an optional sweep grid:

    base experiments   .experiment.json files written by the step GUI
    variant            modifications applied identically to every base experiment
                       (reaction rate, rate × factor, extra add-on chemistry, chamber parameter,
                       injected amounts, initial values, a condition)
    sweep              one or two axes, e.g. initial NO 0.1 … 20 ppb (log) × injected VOC × 0.5 … 4
    metrics            numbers taken from every finished .cs (O3 maximum, Ox gain while lit, OH mean, …)

Every (experiment × grid point) is one run in its own folder <study folder>/runs/<name>/ with its
own .easy / .enz / scenario JSON / .cs, so the base experiments are never touched and any run opens
in the results viewer.  The modifications are applied on the ExperimentSpec, so they reach the
EASY the way the GUI writes them (k[...] = CONST(...) after the includes, #include of add-ons,
injection ppb in the EASY rate) — the FACSIMILE "last definition wins" convention, nothing is
patched into scenario JSONs.  Each run's .easy carries '; STUDY:' provenance lines.

Execution: one scenario_runner process per run, N workers (BLAS/numba threads pinned to 1 per
worker), idempotent — a run with a finished .cs is skipped, so a study can be stopped, resumed and
extended.  Runtime: the first run is measured alone ("pilot") and the estimate for the rest is
n × t_pilot / workers, refined from the median of finished runs.

    python study_runner.py my.study.json --workers 4        # CLI; --dry writes inputs only
"""
from __future__ import annotations

import copy
import csv
import json
import math
import os
import re
import shutil
import subprocess
import sys
import threading
import time
from dataclasses import dataclass, field, asdict
from typing import Callable, Dict, List, Optional, Tuple

import numpy as np

from experiment_spec import (ExperimentSpec, InjectionEvent, load_spec, save_spec, write_experiment,
                             ensure_mechanism, base_fac, ACTIONS)
import mechanism_registry as mreg
from scenario_prep import HERE, load_paths, solver_from_toml, prepare_scenario

PPB = 1e-9

# ─────────────────────────────────────────────────────────────────────────────
# Model
# ─────────────────────────────────────────────────────────────────────────────

MOD_KINDS: Dict[str, Tuple[str, str, str]] = {   # kind → (label, what 'target' is, what 'value' is)
    "rate":            ("Reaction rate = expression", "reaction as in the mechanism, e.g. C5H8 + OH --> CISOPA", "rate expression, 0 switches the reaction off"),
    "rate_factor":     ("Reaction rate × factor", "reaction as in the mechanism", "factor, e.g. 0.8"),
    "addon":           ("Extra add-on chemistry", "add-on id from mechanisms/addons", "(none)"),
    "chamber":         ("Chamber parameter", "field of ChamberEffects, e.g. hono_const, walldep_value, hcho_const", "new value / expression"),
    "injection_scale": ("Injected amount × factor", "injection name, action (voc, no, no2, o3, mixture, …) or 'all'", "factor"),
    "injection_set":   ("Injected amount = ppb", "injection name", "total ppb"),
    "chamber_scale":   ("Chamber parameter × factor", "field of ChamberEffects, e.g. hono_const (chamber HONO source), dark_hono_rate", "factor"),
    "initial":         ("Initial value at t0", "species", "ppb"),
    "condition":       ("Condition", "field of Conditions, e.g. temperature_K, h2o_percent", "new value"),
    "base_mechanism":  ("Base mechanism (MCM version) — rebuilds the chemistry", "base id from mechanisms/base, e.g. MCM_v3.3.1", "(none)"),
}
CHEMISTRY_KINDS = ("addon", "base_mechanism")      # these rebuild the subset once per experiment
AXIS_KINDS: Dict[str, str] = {
    "initial":         "initial value of a species at t0 [ppb]",
    "injection_scale": "× factor on the injected amount of an injection / action / 'all'",
    "injection_set":   "injected amount of an injection [ppb]",
    "rate_factor":     "× factor on a reaction rate",
    "chamber":         "chamber parameter value",
    "chamber_scale":   "× factor on a chamber parameter (e.g. hono_const = chamber HONO source strength)",
    "condition":       "condition value (temperature_K, h2o_percent, …)",
}


@dataclass
class Modification:
    kind: str
    target: str = ""
    value: str = ""
    comment: str = ""


@dataclass
class Axis:
    kind: str
    target: str
    values: List[float] = field(default_factory=list)
    label: str = ""                         # column name in the results; default from kind and target
    also: List[List[str]] = field(default_factory=list)   # further [kind, target] pairs varied together with the same value

    def members(self) -> List[Tuple[str, str]]:
        return [(self.kind, self.target)] + [(k, t) for k, t in self.also]

    def name(self) -> str:
        if self.label:
            return self.label
        if self.also:
            return _safe("+".join(t for _k, t in self.members())) + ("_x" if self.kind in ("injection_scale", "rate_factor", "chamber_scale") else "")
        if self.kind == "initial":
            return f"{self.target}_t0_ppb"
        if self.kind == "injection_set":
            return f"{self.target}_ppb"
        if self.kind in ("injection_scale", "rate_factor"):
            return f"{self.target}_x"
        return f"{self.target}"


@dataclass
class Study:
    name: str = "study"
    folder: str = ""                                    # where runs/ and study_results.csv go
    experiments: List[str] = field(default_factory=list)   # .experiment.json paths
    modifications: List[Modification] = field(default_factory=list)
    axes: List[Axis] = field(default_factory=list)      # 0, 1 or 2 axes
    metrics: List[str] = field(default_factory=lambda: ["o3_max", "ox_lit_delta", "oh_mean_lit", "nox_end", "hcho_max"])
    workers: int = max(1, (os.cpu_count() or 2) - 1)
    export_kv_pd: bool = False
    toml: str = ""
    solver: Dict[str, object] = field(default_factory=dict)   # frozen at creation: identical for every run
    notes: str = ""
    created: str = ""


def save_study(st: Study, path: str) -> str:
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(asdict(st), fh, indent=1, ensure_ascii=False)
    return path


def load_study(path: str) -> Study:
    with open(path, encoding="utf-8") as fh:
        d = json.load(fh)
    d["modifications"] = [Modification(**m) for m in d.get("modifications", [])]
    d["axes"] = [Axis(**a) for a in d.get("axes", [])]
    st = Study(**{k: v for k, v in d.items() if k in Study.__dataclass_fields__})
    if not st.folder:
        st.folder = os.path.dirname(os.path.abspath(path))
    return st


def new_study(name: str, folder: str, toml: str) -> Study:
    st = Study(name=name, folder=folder, toml=toml, solver=solver_from_toml(toml),
               created=time.strftime("%Y-%m-%d %H:%M"))
    return st


def parse_values(text: str) -> List[float]:
    """'0.5, 1, 2'  |  'lin 0 20 5'  |  'log 0.1 20 6'  → list of values."""
    t = text.strip()
    m = re.match(r"^(lin|log)\s+([-\d.eE+]+)\s+([-\d.eE+]+)\s+(\d+)$", t, re.I)
    if m:
        a, b, n = float(m.group(2)), float(m.group(3)), int(m.group(4))
        if m.group(1).lower() == "log":
            if a <= 0 or b <= 0:
                raise ValueError("log range needs positive limits")
            vals = np.logspace(math.log10(a), math.log10(b), n)
        else:
            vals = np.linspace(a, b, n)
        return [float(f"{v:.6g}") for v in vals]
    return [float(x) for x in re.split(r"[,\s;]+", t) if x]


# ─────────────────────────────────────────────────────────────────────────────
# Applying modifications to a spec
# ─────────────────────────────────────────────────────────────────────────────

def _norm_rxn(s: str) -> str:
    return re.sub(r"\s+", "", s).upper()


def _fmt(v) -> str:
    return f"{float(v):.6g}" if not isinstance(v, str) else v


def _coerce(current, text: str):
    """The new value in the type of the field it replaces (bool / int / float / str)."""
    if isinstance(current, bool):
        return str(text).strip().lower() in ("1", "true", "yes", "on")
    if isinstance(current, int):
        return int(float(text))
    if isinstance(current, float):
        return float(text)
    return str(text)


def _scale_expr(expr: str, factor) -> str:
    expr = (expr or "0").strip()
    try:
        return f"{float(expr) * float(factor):.6g}"
    except ValueError:
        return f"({expr})*{_fmt(factor)}"


def _events_matching(spec: ExperimentSpec, target: str) -> List[InjectionEvent]:
    t = target.strip()
    if "," in t:                                   # several names / classes varied together
        out = []
        for part in t.split(","):
            out += [e for e in _events_matching(spec, part) if e not in out]
        return out
    if t.lower() == "all":
        return list(spec.injections)
    by_name = [e for e in spec.injections if e.name.upper() == t.upper()]
    if by_name:
        return by_name
    if t.lower() in ACTIONS:
        return [e for e in spec.injections if e.action_key() == t.lower()]
    if t.lower() == "voc":
        return [e for e in spec.injections if e.action_key() in ("voc", "mixture")]
    if t.lower() == "nox":
        return [e for e in spec.injections if e.action_key() in ("no", "no2")]
    return [e for e in spec.injections if e.products.strip().upper() == t.upper()]


def mcm_rate_of(reaction: str, fac_path: str = "") -> Optional[str]:
    """The rate expression (EASY form) of a reaction of the base MCM, or None."""
    from mcm_species_db import load_fac_reactions
    want = _norm_rxn(reaction)
    for r in load_fac_reactions(fac_path or base_fac()):
        if _norm_rxn(r["easy"]) == want:
            return r.get("easy_rate") or r.get("rate")
    return None


def apply_modification(spec: ExperimentSpec, mod: Modification, fac_path: str = "") -> List[str]:
    """Mutates spec; returns notes (warnings) for the validation report."""
    notes: List[str] = []
    k, t, v = mod.kind, mod.target.strip(), str(mod.value).strip()
    spec.reaction_adjustments = [tuple(x) for x in (spec.reaction_adjustments or [])]
    spec.addons = list(spec.addons or [])
    spec.initials = [tuple(x) for x in (spec.initials or [])]
    if k == "rate":
        spec.reaction_adjustments = [a for a in spec.reaction_adjustments if _norm_rxn(a[0]) != _norm_rxn(t)]
        spec.reaction_adjustments.append((t, v, mod.comment or "study variant"))
    elif k == "rate_factor":
        base_rate = next((a[1] for a in spec.reaction_adjustments if _norm_rxn(a[0]) == _norm_rxn(t)), None)
        if base_rate is None:
            base_rate = mcm_rate_of(t, fac_path)
        if base_rate is None:
            notes.append(f"rate × factor: reaction '{t}' not found in the base mechanism — skipped")
        else:
            spec.reaction_adjustments = [a for a in spec.reaction_adjustments if _norm_rxn(a[0]) != _norm_rxn(t)]
            spec.reaction_adjustments.append((t, f"({base_rate})*{v}", mod.comment or f"study: × {v}"))
    elif k == "addon":
        if t not in spec.addons:
            spec.addons.append(t)
        spec.mechanism_stem = ""                  # subset must be rebuilt with the add-on's seed species
    elif k == "base_mechanism":
        spec.base_mechanism = t
        spec.mechanism_stem = ""                  # subset rebuilt from the other MCM version
    elif k == "chamber_scale":
        if not hasattr(spec.chamber, t):
            notes.append(f"chamber parameter '{t}' unknown — skipped")
        else:
            setattr(spec.chamber, t, _scale_expr(str(getattr(spec.chamber, t)), v))
    elif k in ("chamber", "condition"):
        obj = spec.chamber if k == "chamber" else spec.conditions
        if not hasattr(obj, t):
            notes.append(f"{k} parameter '{t}' unknown — skipped")
        else:
            setattr(obj, t, _coerce(getattr(obj, t), v))
    elif k in ("injection_scale", "injection_set"):
        evs = _events_matching(spec, t)
        if not evs:
            notes.append(f"no injection matches '{t}' — skipped")
        for e in evs:
            if k == "injection_scale":
                e.total_ppb = _scale_expr(e.total_ppb, v)
                e.extra_channels = [(p, _scale_expr(a, v)) for p, a in e.extra_channels]
            else:
                e.total_ppb = _fmt(float(v))
    elif k == "initial":
        spec.initials = [x for x in spec.initials if x[0].upper() != t.upper()]
        spec.initials.append((t.upper(), f"M*1e-9*{_fmt(float(v))}"))
    else:
        notes.append(f"unknown modification '{k}'")
    return notes


def apply_point(spec: ExperimentSpec, axes: List[Axis], point: Dict[str, float], fac_path: str = "") -> List[str]:
    notes = []
    for ax in axes:
        for kind, target in ax.members():
            notes += apply_modification(spec, Modification(kind, target, _fmt(point[ax.name()])), fac_path)
    return notes


def grid_points(axes: List[Axis]) -> List[Dict[str, float]]:
    pts: List[Dict[str, float]] = [{}]
    for ax in axes:
        pts = [dict(p, **{ax.name(): v}) for p in pts for v in ax.values]
    return pts


# ─────────────────────────────────────────────────────────────────────────────
# Plan: one RunPlan per (experiment × point)
# ─────────────────────────────────────────────────────────────────────────────

def _safe(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.+-]+", "_", s).strip("_")


def run_name(base_name: str, point: Dict[str, float]) -> str:
    if not point:
        return base_name
    return base_name + "__" + "__".join(f"{_safe(k)}-{v:g}" for k, v in point.items())


@dataclass
class RunPlan:
    base_path: str
    base_name: str
    point: Dict[str, float]
    name: str
    folder: str
    status: str = "planned"        # planned | ready | running | done | failed | skipped
    elapsed: float = float("nan")
    progress: float = 0.0
    message: str = ""

    @property
    def json_path(self) -> str:
        return os.path.join(self.folder, self.name + ".json")

    @property
    def cs_path(self) -> str:
        return os.path.join(self.folder, self.name + ".cs")

    @property
    def log_path(self) -> str:
        return os.path.join(self.folder, self.name + ".log")

    @property
    def info_path(self) -> str:
        return os.path.join(self.folder, self.name + ".run.json")


def plan_runs(st: Study) -> List[RunPlan]:
    plans = []
    for exp in st.experiments:
        base = os.path.basename(exp)[:-len(".experiment.json")] if exp.endswith(".experiment.json") else os.path.splitext(os.path.basename(exp))[0]
        for pt in grid_points(st.axes):
            name = run_name(base, pt)
            plans.append(RunPlan(exp, base, pt, name, os.path.join(st.folder, "runs", name)))
    return plans


def check_guards(n_runs: int) -> Tuple[str, str]:
    if n_runs > 200:
        return "confirm", f"{n_runs} runs — above 200, this needs an explicit confirmation."
    if n_runs > 50:
        return "warn", f"{n_runs} runs — above 50; check the estimate before starting."
    return "ok", f"{n_runs} runs"


# ─────────────────────────────────────────────────────────────────────────────
# Validation (before anything is written)
# ─────────────────────────────────────────────────────────────────────────────

def _subset_reactions(exp_path: str, stem: str) -> set:
    chem = os.path.join(os.path.dirname(exp_path), "chemistry")
    out = set()
    for f in (f"{stem}.mech", f"{stem}.dil.mech"):
        p = os.path.join(chem, f)
        if os.path.isfile(p):
            for line in open(p, encoding="latin-1"):
                line = line.split(";")[0].strip()
                if "-->" in line and not line.startswith("k["):
                    out.add(_norm_rxn(line))
    return out


def validate_study(st: Study) -> List[Tuple[str, str]]:
    msgs: List[Tuple[str, str]] = []
    if not st.experiments:
        msgs.append(("error", "No base experiment chosen."))
    if not st.folder:
        msgs.append(("error", "No study folder."))
    for exp in st.experiments:
        if not os.path.isfile(exp):
            msgs.append(("error", f"Experiment file not found: {exp}")); continue
        try:
            spec = load_spec(exp)
        except Exception as e:
            msgs.append(("error", f"{os.path.basename(exp)}: cannot read ({e})")); continue
        if not (spec.precursors or spec.mechanism_stem or spec.inorganic_only):
            msgs.append(("error", f"{spec.name}: no chemistry in this experiment"))
        if spec.mechanism_stem and not os.path.isdir(os.path.join(os.path.dirname(exp), "chemistry")) and not any(m.kind in CHEMISTRY_KINDS for m in st.modifications):
            msgs.append(("warning", f"{spec.name}: no chemistry/ folder next to it — the subset will be rebuilt"))
        try:
            mcm = {_norm_rxn(r["easy"]) for r in __import__("mcm_species_db").load_fac_reactions(base_fac(spec.base_mechanism))}
        except Exception:
            mcm = set()
        sub = _subset_reactions(exp, spec.mechanism_stem) if spec.mechanism_stem else set()
        for m in st.modifications + [Modification(k, t, "") for a in st.axes for k, t in a.members()]:
            if m.kind in ("rate", "rate_factor"):
                r = _norm_rxn(m.target)
                if mcm and r not in mcm and r not in sub:
                    msgs.append(("error", f"{spec.name}: reaction '{m.target}' is not in the mechanism"))
                elif sub and r not in sub and not any(x.kind in CHEMISTRY_KINDS for x in st.modifications):
                    msgs.append(("warning", f"{spec.name}: reaction '{m.target}' is in the MCM but not in this experiment's subset (not reachable from its compounds) — the adjustment has no effect there"))
            elif m.kind in ("injection_scale", "injection_set") and not _events_matching(spec, m.target):
                msgs.append(("warning", f"{spec.name}: no injection matches '{m.target}'"))
            elif m.kind == "addon":
                try:
                    ids = {a.id for a in mreg.addons_for(spec.base_mechanism, [m.target])}
                except Exception:
                    ids = set()
                if m.target not in ids:
                    msgs.append(("error", f"add-on '{m.target}' not installed or not compatible with {spec.base_mechanism}"))
            elif m.kind in ("chamber", "chamber_scale") and not hasattr(spec.chamber, m.target):
                msgs.append(("error", f"chamber parameter '{m.target}' unknown"))
            elif m.kind == "base_mechanism" and m.target not in mreg.load_bases():
                msgs.append(("error", f"base mechanism '{m.target}' is not installed (mechanisms/base)"))
            elif m.kind == "condition" and not hasattr(spec.conditions, m.target):
                msgs.append(("error", f"condition '{m.target}' unknown"))
    for ax in st.axes:
        if not ax.values:
            msgs.append(("error", f"axis '{ax.name()}' has no values"))
    n = len(st.experiments) * max(1, len(grid_points(st.axes)))
    lvl, txt = check_guards(n)
    if lvl != "ok":
        msgs.append(("warning", txt))
    return msgs


# ─────────────────────────────────────────────────────────────────────────────
# Materialisation: spec → run folder (easy, enz, json)
# ─────────────────────────────────────────────────────────────────────────────

def _chemistry_for(st: Study, exp_path: str, spec_base: ExperimentSpec) -> Tuple[str, dict]:
    """The chemistry folder every run of this base experiment shares, and the spec fields that
    point at it.  Rebuilt once per study in <study>/chemistry/<base>/ when the variant changes
    the chemistry (add-on); otherwise the base experiment's own chemistry/ folder."""
    chem_mods = [m for m in st.modifications if m.kind in CHEMISTRY_KINDS]
    base_dir = os.path.dirname(os.path.abspath(exp_path))
    if not chem_mods and spec_base.mechanism_stem and os.path.isdir(os.path.join(base_dir, "chemistry")):
        return os.path.join(base_dir, "chemistry"), {}
    out = os.path.join(st.folder, "chemistry", spec_base.name)
    marker = os.path.join(out, "chemistry.json")
    if os.path.isfile(marker):
        return os.path.join(out, "chemistry"), json.load(open(marker, encoding="utf-8"))
    spec = copy.deepcopy(spec_base)
    for m in chem_mods:
        apply_modification(spec, m)
    if not spec.mechanism_stem or chem_mods:
        spec.mechanism_stem = "" if spec.precursors or spec.inorganic_only else spec.mechanism_stem
    ensure_mechanism(spec, out)
    fields = {"mechanism_stem": spec.mechanism_stem, "mech_incl_suffix": spec.mech_incl_suffix,
              "use_anorg": spec.use_anorg, "extra_includes": list(spec.extra_includes), "addons": list(spec.addons)}
    os.makedirs(out, exist_ok=True)
    json.dump(fields, open(marker, "w", encoding="utf-8"), indent=1)
    return os.path.join(out, "chemistry"), fields


def materialise(st: Study, plan: RunPlan, include_dir: str, force: bool = False) -> str:
    """Write the run's spec, EASY, ENZ and scenario JSON.  Returns the JSON path."""
    if not force and os.path.isfile(plan.json_path) and os.path.isfile(plan.json_path[:-5] + ".easy"):
        plan.status = "ready"; return plan.json_path
    spec = load_spec(plan.base_path)
    spec.addons = list(spec.addons or []); spec.reaction_adjustments = list(spec.reaction_adjustments or [])
    spec.initials = list(spec.initials or []); spec.overrides = list(spec.overrides or [])
    base_name = spec.name
    chem_dir, chem_fields = _chemistry_for(st, plan.base_path, spec)
    notes: List[str] = []
    for m in st.modifications:
        if m.kind not in CHEMISTRY_KINDS:
            notes += apply_modification(spec, m)
    for k, v in chem_fields.items():
        setattr(spec, k, v)
    notes += apply_point(spec, st.axes, plan.point)
    if spec.conditions.source == "archive" and spec.conditions.date:
        pass                                   # ENZ rebuilt from the archive as the base was; the base .enz is the fallback below
    spec.name = plan.name
    spec.overrides = [f"; STUDY: {st.name}  base experiment: {base_name}  run: {plan.name}",
                      f"; STUDY: variant: " + ("; ".join(f"{m.kind} {m.target} {m.value}".strip() for m in st.modifications) or "none"),
                      f"; STUDY: point: " + (", ".join(f"{k} = {v:g}" for k, v in plan.point.items()) or "base"),
                      f"; STUDY: solver: " + ", ".join(f"{k}={v}" for k, v in st.solver.items())] + spec.overrides
    os.makedirs(plan.folder, exist_ok=True)
    base_enz = os.path.join(os.path.dirname(plan.base_path), base_name + ".enz")
    try:
        res = write_experiment(spec, plan.folder)
    except Exception as e:
        # the archive may be unreachable (Z: not mounted): the base experiment's ENZ is reused as it is
        # (its switch columns are the same pulses; only the EASY amounts change)
        if spec.conditions.source == "archive" and os.path.isfile(base_enz):
            from experiment_spec import render_easy
            enz = shutil.copy2(base_enz, os.path.join(plan.folder, plan.name + ".enz"))
            easy = os.path.join(plan.folder, plan.name + ".easy")
            with open(easy, "w", encoding="latin-1", errors="replace", newline="\n") as fh:
                fh.write(render_easy(spec))
            res = {"enz": enz, "easy": easy, "notes": []}
            notes.append(f"archive not readable ({e}); ENZ copied from {os.path.basename(base_enz)}")
        else:
            raise
    save_spec(spec, os.path.join(plan.folder, plan.name + ".experiment.json"))
    js = prepare_scenario(res, include_dir, st.export_kv_pd, st.toml, extra_dirs=[chem_dir], solver=st.solver or None)
    info = _read_info(plan); info.update({"base": base_name, "point": plan.point, "notes": notes, "materialised": time.strftime("%Y-%m-%d %H:%M:%S")})
    _write_info(plan, info)
    plan.status = "ready"; plan.message = "; ".join(notes)
    return js


def _read_info(plan: RunPlan) -> dict:
    try:
        return json.load(open(plan.info_path, encoding="utf-8"))
    except Exception:
        return {}


def _write_info(plan: RunPlan, info: dict) -> None:
    os.makedirs(plan.folder, exist_ok=True)
    json.dump(info, open(plan.info_path, "w", encoding="utf-8"), indent=1, default=str)


def is_done(plan: RunPlan) -> bool:
    if not os.path.isfile(plan.cs_path) or not os.path.isfile(plan.json_path):
        return False
    if os.path.getmtime(plan.cs_path) < os.path.getmtime(plan.json_path):
        return False
    return _read_info(plan).get("exit", 1) == 0


# ─────────────────────────────────────────────────────────────────────────────
# Metrics from a finished .cs
# ─────────────────────────────────────────────────────────────────────────────

METRICS: Dict[str, Tuple[str, str]] = {        # key → (label, unit)
    "o3_max":        ("O3 maximum", "ppb"),
    "o3_end":        ("O3 at the end", "ppb"),
    "o3_lit_delta":  ("O3 gain while lit", "ppb"),
    "ox_lit_delta":  ("Ox (O3+NO2) gain while lit", "ppb"),
    "o3_time_max":   ("Time of the O3 maximum", "h UTC"),
    "oh_mean_lit":   ("OH mean while lit", "cm-3"),
    "oh_max":        ("OH maximum", "cm-3"),
    "ho2_mean_lit":  ("HO2 mean while lit", "cm-3"),
    "ro2_mean_lit":  ("RO2 mean while lit", "cm-3"),
    "nox_end":       ("NOx at the end", "ppb"),
    "no2_max":       ("NO2 maximum", "ppb"),
    "hono_max":      ("HONO maximum", "ppb"),
    "hcho_max":      ("HCHO maximum", "ppb"),
    "h2o2_max":      ("H2O2 maximum", "ppb"),
    "pan_max":       ("PAN maximum", "ppb"),
    "hno3_end":      ("HNO3 at the end", "ppb"),
}


def compute_metrics(cs_path: str, keys: List[str]) -> Dict[str, float]:
    from experiment_plots import read_columns
    t, d = read_columns(cs_path, ["M", "HV", "O3", "NO", "NO2", "OH", "HO2", "RO2", "HONO", "HCHO", "H2O2", "PAN", "HNO3"])
    t = np.asarray(t, float)
    M = np.asarray(d.get("M"), float) if "M" in d else None

    def ppb(name):
        if name not in d or M is None:
            return None
        return np.asarray(d[name], float) / M / PPB

    def raw(name):
        return np.asarray(d[name], float) if name in d else None
    lit = np.asarray(d["HV"], float) > 0.5 if "HV" in d else np.ones(len(t), bool)
    if not lit.any():
        lit = np.ones(len(t), bool)
    day0 = math.floor((t[0] + 946684800.0) / 86400.0) * 86400.0 - 946684800.0   # midnight of the first day, ENZ seconds
    hours = (t - day0) / 3600.0
    o3, no, no2 = ppb("O3"), ppb("NO"), ppb("NO2")
    out: Dict[str, float] = {}
    for k in keys:
        v = float("nan")
        try:
            if k == "o3_max" and o3 is not None: v = float(np.nanmax(o3))
            elif k == "o3_end" and o3 is not None: v = float(o3[-1])
            elif k == "o3_time_max" and o3 is not None: v = float(hours[int(np.nanargmax(o3))])
            elif k == "o3_lit_delta" and o3 is not None:
                i = np.where(lit)[0]; v = float(o3[i[-1]] - o3[i[0]])
            elif k == "ox_lit_delta" and o3 is not None and no2 is not None:
                i = np.where(lit)[0]; ox = o3 + no2; v = float(ox[i[-1]] - ox[i[0]])
            elif k == "oh_mean_lit" and raw("OH") is not None: v = float(np.nanmean(raw("OH")[lit]))
            elif k == "oh_max" and raw("OH") is not None: v = float(np.nanmax(raw("OH")))
            elif k == "ho2_mean_lit" and raw("HO2") is not None: v = float(np.nanmean(raw("HO2")[lit]))
            elif k == "ro2_mean_lit" and raw("RO2") is not None: v = float(np.nanmean(raw("RO2")[lit]))
            elif k == "nox_end" and no is not None and no2 is not None: v = float(no[-1] + no2[-1])
            elif k == "no2_max" and no2 is not None: v = float(np.nanmax(no2))
            elif k == "hono_max" and ppb("HONO") is not None: v = float(np.nanmax(ppb("HONO")))
            elif k == "hcho_max" and ppb("HCHO") is not None: v = float(np.nanmax(ppb("HCHO")))
            elif k == "h2o2_max" and ppb("H2O2") is not None: v = float(np.nanmax(ppb("H2O2")))
            elif k == "pan_max" and ppb("PAN") is not None: v = float(np.nanmax(ppb("PAN")))
            elif k == "hno3_end" and ppb("HNO3") is not None: v = float(ppb("HNO3")[-1])
        except Exception:
            v = float("nan")
        out[k] = v
    return out


# ─────────────────────────────────────────────────────────────────────────────
# Execution
# ─────────────────────────────────────────────────────────────────────────────

_FRACTION = re.compile(r"Fraction Done:\s*([\d.]+)")


def _worker_env() -> dict:
    env = dict(os.environ)
    env.update({"OMP_NUM_THREADS": "1", "MKL_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "NUMBA_NUM_THREADS": "1",
                "PYTHONIOENCODING": "utf-8", "PYTHONUTF8": "1"})
    return env


def _tail_fraction(log_path: str) -> float:
    try:
        with open(log_path, "rb") as fh:
            fh.seek(0, 2); size = fh.tell(); fh.seek(max(0, size - 4000))
            txt = fh.read().decode("utf-8", errors="replace")
        m = _FRACTION.findall(txt)
        return float(m[-1]) if m else 0.0
    except Exception:
        return 0.0


def prior_seconds(st: Study) -> float:
    """Before a pilot: the last run time of the base experiments (<name>.run.json written by the step
    window), median over the experiments; 0 when none has run yet."""
    times = []
    for exp in st.experiments:
        p = exp[:-len(".experiment.json")] + ".run.json" if exp.endswith(".experiment.json") else ""
        try:
            t = float(json.load(open(p, encoding="utf-8")).get("elapsed", 0))
            if t > 0:
                times.append(t)
        except Exception:
            pass
    return float(np.median(times)) if times else 0.0


def estimate_seconds(n_remaining: int, per_run: float, workers: int) -> float:
    if not n_remaining or not per_run or math.isnan(per_run):
        return 0.0
    return math.ceil(n_remaining / max(1, workers)) * per_run


def fmt_duration(s: float) -> str:
    if not s or math.isnan(s):
        return "—"
    s = int(round(s)); h, r = divmod(s, 3600); m, sec = divmod(r, 60)
    return f"{h} h {m:02d} min" if h else (f"{m} min {sec:02d} s" if m else f"{sec} s")


class StudyRunner:
    """Materialises and runs the plans of a study with N worker processes.

    on_event(kind, plan, runner) is called from the runner thread with kind in
    'materialised' | 'started' | 'progress' | 'finished' | 'failed' | 'skipped' | 'all_done'.
    """

    def __init__(self, st: Study, include_dir: str = "", python: str = sys.executable):
        self.study = st
        self.include_dir = include_dir or load_paths(st.toml)["include_dir"]
        self.python = python
        self.plans: List[RunPlan] = plan_runs(st)
        self.stop_event = threading.Event()
        self.elapsed: List[float] = []
        self._thread: Optional[threading.Thread] = None

    # ── inputs ──
    def materialise_all(self, on_event: Callable = None, force: bool = False) -> List[str]:
        errors = []
        for p in self.plans:
            try:
                materialise(self.study, p, self.include_dir, force=force)
                if on_event:
                    on_event("materialised", p, self)
            except Exception as e:
                p.status = "failed"; p.message = f"inputs: {e}"; errors.append(f"{p.name}: {e}")
                if on_event:
                    on_event("failed", p, self)
        return errors

    # ── one run (blocking) ──
    def run_one(self, p: RunPlan, on_event: Callable = None) -> float:
        if is_done(p):
            p.status = "done"; p.elapsed = _read_info(p).get("elapsed", float("nan")); p.progress = 1.0
            if on_event:
                on_event("skipped", p, self)
            return p.elapsed
        if not os.path.isfile(p.json_path):
            materialise(self.study, p, self.include_dir)
        args = [self.python, "-u", os.path.join(HERE, "scenario_runner.py"), p.json_path]
        if self.study.toml:
            args += ["--toml", self.study.toml]
        t0 = time.monotonic(); p.status = "running"; p.progress = 0.0
        if on_event:
            on_event("started", p, self)
        with open(p.log_path, "w", encoding="utf-8") as log:
            proc = subprocess.Popen(args, cwd=HERE, env=_worker_env(), stdout=log, stderr=subprocess.STDOUT)
            while proc.poll() is None:
                if self.stop_event.is_set():
                    proc.kill(); p.status = "failed"; p.message = "stopped"; return float("nan")
                time.sleep(0.5)
                fr = _tail_fraction(p.log_path)
                if fr > p.progress:
                    p.progress = fr
                    if on_event:
                        on_event("progress", p, self)
        p.elapsed = time.monotonic() - t0
        ok = proc.returncode == 0 and os.path.isfile(p.cs_path)
        info = _read_info(p); info.update({"exit": proc.returncode, "elapsed": p.elapsed, "finished": time.strftime("%Y-%m-%d %H:%M:%S")})
        if ok:
            try:
                info["metrics"] = compute_metrics(p.cs_path, list(METRICS))
            except Exception as e:
                info["metrics_error"] = str(e)
            try:
                info["cs_bytes"] = os.path.getsize(p.cs_path)
            except OSError:
                pass
        _write_info(p, info)
        p.status = "done" if ok else "failed"; p.progress = 1.0 if ok else p.progress
        if not ok:
            p.message = f"exit {proc.returncode}, see {os.path.basename(p.log_path)}"
        else:
            self.elapsed.append(p.elapsed)
        if on_event:
            on_event("finished" if ok else "failed", p, self)
        return p.elapsed

    # ── all remaining runs with N workers (blocking; call from a thread) ──
    def run_all(self, on_event: Callable = None, workers: int = 0, plans: List[RunPlan] = None) -> None:
        workers = max(1, workers or self.study.workers)
        todo = [p for p in (plans or self.plans) if not is_done(p)]
        for p in (plans or self.plans):
            if is_done(p) and p.status != "done":
                p.status = "done"; p.elapsed = _read_info(p).get("elapsed", float("nan")); p.progress = 1.0
                if on_event:
                    on_event("skipped", p, self)
        lock = threading.Lock(); queue = list(todo)

        def work():
            while not self.stop_event.is_set():
                with lock:
                    if not queue:
                        return
                    p = queue.pop(0)
                self.run_one(p, on_event)
        threads = [threading.Thread(target=work, daemon=True) for _ in range(min(workers, len(queue)))]
        for th in threads:
            th.start()
        for th in threads:
            th.join()
        if on_event:
            on_event("all_done", None, self)

    def start(self, on_event: Callable = None, workers: int = 0, plans: List[RunPlan] = None) -> threading.Thread:
        self.stop_event.clear()
        self._thread = threading.Thread(target=self.run_all, args=(on_event, workers, plans), daemon=True)
        self._thread.start()
        return self._thread

    def stop(self):
        self.stop_event.set()

    # ── estimates ──
    def per_run_seconds(self) -> float:
        done = [x for x in self.elapsed if x and not math.isnan(x)]
        if not done:
            done = [_read_info(p).get("elapsed") for p in self.plans if is_done(p)]
            done = [x for x in done if x]
        return float(np.median(done)) if done else float("nan")

    def remaining(self) -> int:
        return sum(1 for p in self.plans if not is_done(p) and p.status not in ("done", "running"))

    def estimate(self, workers: int = 0) -> Tuple[float, str]:
        s = estimate_seconds(self.remaining(), self.per_run_seconds(), workers or self.study.workers)
        return s, fmt_duration(s)

    def disk_estimate(self) -> float:
        sizes = [_read_info(p).get("cs_bytes", 0) for p in self.plans if is_done(p)]
        return float(np.median(sizes)) * self.remaining() if sizes else float("nan")


# ─────────────────────────────────────────────────────────────────────────────
# Results table
# ─────────────────────────────────────────────────────────────────────────────

def collect_results(st: Study, plans: List[RunPlan] = None) -> List[dict]:
    rows = []
    for p in (plans or plan_runs(st)):
        info = _read_info(p)
        row = {"run": p.name, "experiment": p.base_name, **{k: v for k, v in p.point.items()},
               "status": "done" if is_done(p) else p.status, "elapsed_s": info.get("elapsed"), "cs": p.cs_path if is_done(p) else ""}
        row.update({k: info.get("metrics", {}).get(k) for k in st.metrics})
        rows.append(row)
    return rows


def write_results(st: Study, rows: List[dict]) -> str:
    path = os.path.join(st.folder, f"{st.name}_results.csv")
    if not rows:
        return path
    keys = list(rows[0].keys())
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=keys); w.writeheader(); w.writerows(rows)
    return path


def normalised_sensitivity(xs: List[float], ys: List[float]) -> List[float]:
    """d ln y / d ln x between neighbouring grid points (percent change of the metric per percent
    change of the parameter); NaN where a log is undefined."""
    out = []
    for i in range(1, len(xs)):
        try:
            out.append(math.log(ys[i] / ys[i - 1]) / math.log(xs[i] / xs[i - 1]))
        except (ValueError, ZeroDivisionError, TypeError):
            out.append(float("nan"))
    return out
