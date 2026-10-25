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
    "injection_new":   ("New injection (not in the experiment)", "species|HH:MM|minutes, e.g. NO|11:00|5", "total ppb (0 = none)"),
    "condition":       ("Condition", "field of Conditions, e.g. temperature_K, h2o_percent", "new value"),
    "base_mechanism":  ("Base mechanism (MCM version) — rebuilds the chemistry", "base id from mechanisms/base, e.g. MCM_v3.3.1", "(none)"),
}
CHEMISTRY_KINDS = ("addon", "base_mechanism")      # these rebuild the subset once per experiment
AXIS_KINDS: Dict[str, str] = {
    "initial":         "initial value of a species at t0 [ppb]",
    "injection_scale": "× factor on the injected amount of an injection / action / 'all'",
    "injection_set":   "injected amount of an injection [ppb]",
    "injection_new":   "amount of an injection added to the experiment [ppb]; target species|HH:MM|minutes",
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
        if self.kind == "injection_new":
            sp, start, _dur = parse_new_injection(self.target)
            return f"{sp}_added_{start.replace(':', '')}_ppb"
        return f"{self.target}"


@dataclass
class Study:
    name: str = "study"
    folder: str = ""                                    # where runs/ and study_results.csv go
    experiments: List[str] = field(default_factory=list)   # .experiment.json paths
    modifications: List[Modification] = field(default_factory=list)
    axes: List[Axis] = field(default_factory=list)      # 0, 1 or 2 axes
    metrics: List[str] = field(default_factory=lambda: ["po3_lit_int", "o3_max", "ox_lit_delta", "po3_mean_lit", "ope", "regime", "nox_end", "ho2_ro2"])
    workers: int = max(1, (os.cpu_count() or 2) - 1)
    export_kv_pd: bool = False
    toml: str = ""
    solver: Dict[str, object] = field(default_factory=dict)   # frozen at creation: identical for every run
    notes: str = ""
    created: str = ""
    mode: str = ""                                      # "sweep" (one experiment, its own knobs) | "batch" (several experiments)
    focus: str = "ozone"                                # key of FOCUS: which results are ticked by default


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


def parse_new_injection(target: str) -> Tuple[str, str, float]:
    """'NO|11:00|5' → ('NO', '11:00', 5.0); duration defaults to the SAPHIR action default, else 1 min."""
    parts = [x.strip() for x in target.split("|")]
    sp = parts[0].upper() if parts else ""
    start = parts[1] if len(parts) > 1 and re.fullmatch(r"\d{1,2}:\d{2}", parts[1]) else ""
    try:
        dur = float(parts[2]) if len(parts) > 2 and parts[2] else None
    except ValueError:
        dur = None
    if dur is None:
        key = {"NO": "no", "NO2": "no2", "O3": "o3", "CO": "co"}.get(sp)
        dur = float(ACTIONS[key]["duration_min"]) if key else 1.0
    return sp, start, dur


def add_injection(spec: ExperimentSpec, target: str, ppb: str) -> List[str]:
    """Add (or set) an injection the experiment does not have: species at HH:MM over minutes, total ppb.
    NO / NO2 / O3 / CO use the SAPHIR names and switch columns (QNO, …), other species Q_<name>.  If the
    experiment already injects that species, the new one is named <species>_add with its own column."""
    sp, start, dur = parse_new_injection(target)
    if not sp or not start:
        return [f"new injection '{target}': give species|HH:MM|minutes"]
    key = {"NO": "no", "NO2": "no2", "O3": "o3", "CO": "co"}.get(sp, "voc")
    tag = f"{sp}@{start}"
    ev = next((e for e in spec.injections if e.comment.startswith(f"study: added {tag}")), None)
    if ev is None:
        taken = {e.name.upper() for e in spec.injections}
        name = sp if sp not in taken else f"{sp}_add"
        col, sw = "", ""
        if name == sp and key != "voc":
            col, sw = ACTIONS[key]["column"], ACTIONS[key]["switch"]
        else:
            col, sw = f"Q_{name}", f"q_{name}"
        ev = InjectionEvent(name, ACTIONS[key]["products"] or sp if key != "voc" else sp, _fmt(float(ppb)), start, dur, col, sw,
                            comment=f"study: added {tag} ({dur:g} min)", action=key)
        spec.injections.append(ev)
    else:
        ev.total_ppb = _fmt(float(ppb))
    if key == "voc" and sp not in [p.upper() for p in (spec.precursors or [])] and not spec.mechanism_stem:
        return [f"new injection of {sp}: it is not among the experiment's compounds — add it on the Chemistry step if it should react"]
    return []


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
    elif k == "injection_new":
        notes += add_injection(spec, t, v or "0")
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
    "po3_mean_lit":  ("P(O3) mean, lit", "ppb/h"),
    "po3_lit_int":   ("P(O3) integral, lit", "ppb"),
    "ope":           ("OPE ΔOx/ΔNOz", "-"),
    "regime":        ("Regime P(H2O2)/P(HNO3)", "-"),
    "ho2_ro2":       ("HO2/RO2 while lit", "-"),
    "oh_exposure":   ("OH exposure ∫OH dt", "cm-3 s"),
    "voc_consumed":  ("Injected VOC consumed", "%"),
}

# A focus ticks a default set of results; the first one is plotted first.
FOCUS: Dict[str, Tuple[str, List[str]]] = {
    "ozone":    ("Ozone production / isopleth", ["po3_lit_int", "o3_max", "ox_lit_delta", "po3_mean_lit", "ope", "regime", "nox_end", "ho2_ro2"]),
    "radicals": ("Radicals (HOx / ROx)", ["oh_mean_lit", "ho2_mean_lit", "ro2_mean_lit", "ho2_ro2", "oh_exposure", "po3_mean_lit", "hono_max"]),
    "nox":      ("NOx fate", ["nox_end", "hno3_end", "pan_max", "ope", "no2_max", "regime"]),
    "voc":      ("VOC oxidation", ["voc_consumed", "oh_exposure", "oh_mean_lit", "hcho_max", "h2o2_max"]),
}
METRIC_HELP: Dict[str, str] = {
    "po3_mean_lit": "gross ozone production k(HO2+NO)[HO2][NO] + k(RO2+NO)[RO2][NO] (MCM: 3.45e-12·e^(270/T), KRO2NO 2.7e-12·e^(360/T)), mean while lit",
    "po3_lit_int": "the same production rate integrated over the lit period — the ozone made, before its losses",
    "ope": "Ox (O3+NO2) gained per NOx turned into products while lit; NOz = all nitrogen species except NO and NO2 (HONO, HNO3, N2O5, PANs, organic nitrates, particulate nitrate), counted by their N atoms",
    "regime": "radical loss via peroxides vs via nitric acid: P(H2O2) = k(HO2+HO2)[HO2]², P(HNO3) = KMT08[OH][NO2], ratio of their lit-period integrals; high = NOx-limited, low = VOC-limited, transition typically ~0.3–0.5",
    "ho2_ro2": "mean HO2 / mean RO2 while lit (RO2 = the mechanism's RO2 sum)",
    "oh_exposure": "∫ OH dt over the whole run",
    "voc_consumed": "1 − Σ precursors at the end / Σ precursors at their maximum (the experiment's chosen compounds)",
}


def compute_metrics(cs_path: str, keys: List[str]) -> Dict[str, float]:
    from experiment_plots import read_columns
    base = ["M", "T", "H2O", "HV", "O3", "NO", "NO2", "OH", "HO2", "RO2", "HONO", "HCHO", "H2O2", "PAN", "HNO3"]
    extra = []
    if any(k in keys for k in ("ope",)):
        extra += _nitrogen_species(cs_path)
    if "voc_consumed" in keys:
        extra += _precursors_for(cs_path)
    t, d = read_columns(cs_path, base + [x for x in extra if x not in base])
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
            elif k in ("po3_mean_lit", "po3_lit_int", "regime", "ho2_ro2", "oh_exposure", "ope", "voc_consumed"):
                v = _derived(k, t, d, lit, M, cs_path)
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


# Runs in parallel slow each other down (memory bandwidth, CPU boost clock): measured 2026-10-25 on a
# 10-core office PC, a run that takes 49 s alone took ~135 s with 11 in parallel → ~0.18 per extra run.
PARALLEL_SLOWDOWN = 0.18


def parallel_factor(workers: int) -> float:
    return 1.0 + PARALLEL_SLOWDOWN * max(0, workers - 1)


def estimate_seconds(n_remaining: int, per_run: float, workers: int, alone: bool = True) -> float:
    """Wall time for n runs.  `alone`: per_run was measured with no other run going (the pilot, or the
    base experiment) — then each run is assumed slower by parallel_factor(workers)."""
    if not n_remaining or not per_run or math.isnan(per_run):
        return 0.0
    w = max(1, min(workers, n_remaining))
    return math.ceil(n_remaining / w) * per_run * (parallel_factor(w) if alone else 1.0)


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
    def measured_alone(self) -> bool:
        """True while the only finished runs ran one at a time (the pilot)."""
        return len([x for x in self.elapsed if x and not math.isnan(x)]) <= 1 and             sum(1 for p in self.plans if is_done(p)) <= 1

    def per_run_seconds(self) -> float:
        done = [x for x in self.elapsed if x and not math.isnan(x)]
        if not done:
            done = [_read_info(p).get("elapsed") for p in self.plans if is_done(p)]
            done = [x for x in done if x]
        return float(np.median(done)) if done else float("nan")

    def remaining(self) -> int:
        return sum(1 for p in self.plans if not is_done(p) and p.status not in ("done", "running"))

    def estimate(self, workers: int = 0) -> Tuple[float, str]:
        s = estimate_seconds(self.remaining(), self.per_run_seconds(), workers or self.study.workers, alone=self.measured_alone())
        return s, fmt_duration(s)

    def disk_estimate(self) -> float:
        sizes = [_read_info(p).get("cs_bytes", 0) for p in self.plans if is_done(p)]
        return float(np.median(sizes)) * self.remaining() if sizes else float("nan")


# ─────────────────────────────────────────────────────────────────────────────
# Results table
# ─────────────────────────────────────────────────────────────────────────────

# What each result is and what to expect from it across a sweep (shown under the plot).
METRIC_MEANING: Dict[str, Tuple[str, str]] = {
    "o3_max": ("highest O3 of the run",
               "rises with VOC everywhere; with NOx it rises in the NOx-sensitive part and falls beyond the ridge (NOx-saturated, titration at very high NO)."),
    "o3_end": ("O3 at the last time step",
               "like the maximum, lowered by dilution and by NO titration late in the run."),
    "o3_lit_delta": ("O3 at the end of the lit period minus at its start",
                     "net ozone made in the light; can be negative when NO injected in the light titrates O3 faster than it is produced."),
    "ox_lit_delta": ("Ox = O3 + NO2 gained while lit",
                     "net photochemical ozone without the NO–O3 titration artefact; the cleanest single number for an isopleth."),
    "o3_time_max": ("hour (UTC) of the O3 maximum",
                    "later when NOx is high (titration delays the build-up) or when VOC is low (slow production); early when NOx runs out."),
    "oh_mean_lit": ("mean OH while lit",
                    "highest at intermediate NOx; falls at high NOx (OH + NO2 → HNO3) and at high VOC (OH consumed by the VOC)."),
    "oh_max": ("highest OH of the run",
               "follows the mean; spikes right after an NO injection (HO2 + NO → OH) or HONO photolysis at roof opening."),
    "ho2_mean_lit": ("mean HO2 while lit",
                     "rises with VOC; falls strongly with NOx (HO2 + NO)."),
    "ro2_mean_lit": ("mean RO2 sum while lit",
                     "rises with VOC; falls with NOx (RO2 + NO); large values at low NOx mean RO2 self/cross reactions matter."),
    "nox_end": ("NO + NO2 at the end",
                "close to 0 where the run became NOx-limited; grows with the NOx injected and with the chamber HONO source."),
    "no2_max": ("highest NO2 of the run",
                "rises with the NOx injected; NO → NO2 conversion is faster with more VOC (more peroxy radicals)."),
    "hono_max": ("highest HONO of the run",
                 "set by the chamber HONO source and humidity; barely depends on the injections."),
    "hcho_max": ("highest HCHO of the run",
                 "a first-generation product: rises with the VOC oxidised, i.e. with VOC and with OH."),
    "h2o2_max": ("highest H2O2 of the run",
                 "a NOx-limited marker: large at low NOx (HO2 + HO2), suppressed at high NOx."),
    "pan_max": ("highest PAN of the run",
                "needs both acyl peroxy radicals (VOC) and NO2: largest where both are high; thermally unstable at high T."),
    "hno3_end": ("HNO3 at the end",
                 "a NOx-saturated marker: rises with NOx (OH + NO2), a sink of both radicals and NOx."),
    "po3_mean_lit": ("mean gross ozone production rate while lit, k[HO2][NO] + k[RO2][NO]",
                     "the isopleth of production itself: ridge where adding NOx stops helping; NOx-sensitive below/left, VOC-sensitive beyond."),
    "po3_lit_int": ("ozone produced over the lit period (gross, before its losses)",
                    "larger than the Ox gain by the O3 lost to photolysis, OH/HO2 and dilution; same regime pattern as the rate."),
    "ope": ("Ox made per NOx turned into products (HNO3, PANs, nitrates, …)",
            "high (10–20+) at low NOx where each NOx cycles many times; low (< 5) at high NOx where NOx is removed quickly."),
    "regime": ("radical loss to peroxides / to HNO3",
               "high (> ~0.5) = NOx-sensitive, low (< ~0.2) = VOC-sensitive / NOx-saturated; transition near 0.35 (red dotted line)."),
    "ho2_ro2": ("mean HO2 / mean RO2 while lit",
                "≈ 1 for small VOC; below 1 for large VOC with many RO2 steps; rises when NO converts RO2 quickly to HO2."),
    "oh_exposure": ("∫ OH dt over the run",
                    "how much oxidation the air saw; with k_OH of a VOC it gives the fraction reacted (1 − e^(−k·exposure))."),
    "voc_consumed": ("share of the injected compounds reacted by the end",
                     "rises with OH exposure; near 100 % for reactive VOC (isoprene, monoterpenes) in a full sunny day."),
}


def _nitrogen_species(cs_path: str) -> List[str]:
    """Species of the run with nitrogen (MCM formula), except NO and NO2, with NA (particulate nitrate)."""
    try:
        from mcm_species_db import MCM_DB
        from experiment_plots import cs_species
        names = cs_species(cs_path)
    except Exception:
        return []
    out = []
    for n in names:
        sp = MCM_DB.get(n.upper())
        if (sp is not None and sp.n_N > 0 and n.upper() not in ("NO", "NO2")) or n.upper() in ("NA", "HONO", "N2O5", "NO3", "HNO3"):
            out.append(n)
    return out


def _n_atoms(name: str) -> int:
    if name.upper() in ("NA", "HONO", "NO3", "HNO3"):
        return 1
    if name.upper() == "N2O5":
        return 2
    try:
        from mcm_species_db import MCM_DB
        sp = MCM_DB.get(name.upper())
        return int(sp.n_N) if sp is not None else 0
    except Exception:
        return 0


def _precursors_for(cs_path: str) -> List[str]:
    try:
        from experiment_plots import spec_for_run
        spec = spec_for_run(cs_path)
        return [p.upper() for p in (spec.precursors or [])] if spec else []
    except Exception:
        return []


def _kmt08(T: np.ndarray, M: np.ndarray) -> np.ndarray:
    """OH + NO2 (+M) → HNO3, MCM v3.3.1 KMT08."""
    k0 = 3.2e-30 * M * (T / 300.0) ** -4.5; ki = 3.0e-11; fc = 0.41
    nc = 0.75 - 1.27 * np.log10(fc)
    f = 10 ** (np.log10(fc) / (1 + (np.log10(k0 / ki) / nc) ** 2))
    return (k0 * ki) * f / (k0 + ki)


def _derived(key: str, t: np.ndarray, d: Dict[str, np.ndarray], lit: np.ndarray, M: Optional[np.ndarray], cs_path: str) -> float:
    """Production rates, efficiency, regime indicator, exposure, VOC consumption (see METRIC_HELP)."""
    g = lambda n: np.asarray(d[n], dtype=float) if n in d else None
    T = g("T") if g("T") is not None else np.full(len(t), 298.0)
    dt = np.diff(t, prepend=t[0])
    if key == "oh_exposure":
        return float(np.nansum(g("OH") * dt)) if g("OH") is not None else float("nan")
    if key == "voc_consumed":
        names = [n for n in _precursors_for(cs_path) if n in d]
        if not names:
            return float("nan")
        tot = sum(np.asarray(d[n], dtype=float) for n in names)
        return float(100.0 * (1.0 - tot[-1] / np.nanmax(tot))) if np.nanmax(tot) > 0 else float("nan")
    if M is None:
        return float("nan")
    ho2, ro2, no, oh, no2 = g("HO2"), g("RO2"), g("NO"), g("OH"), g("NO2")
    if key == "ho2_ro2":
        return float(np.nanmean(ho2[lit]) / np.nanmean(ro2[lit])) if ho2 is not None and ro2 is not None and np.nanmean(ro2[lit]) > 0 else float("nan")
    if key in ("po3_mean_lit", "po3_lit_int"):
        if ho2 is None or no is None:
            return float("nan")
        p = 3.45e-12 * np.exp(270.0 / T) * ho2 * no + (2.7e-12 * np.exp(360.0 / T) * ro2 * no if ro2 is not None else 0.0)
        p_ppb_s = p / M * 1e9
        if key == "po3_mean_lit":
            return float(np.nanmean(p_ppb_s[lit]) * 3600.0)
        return float(np.nansum((p_ppb_s * dt)[lit]))
    if key == "regime":
        if ho2 is None or oh is None or no2 is None:
            return float("nan")
        h2o = g("H2O") if g("H2O") is not None else 0.01 * M
        kmt06 = 1 + 1.40e-21 * np.exp(2200.0 / T) * h2o
        k_ho2 = 2.20e-13 * kmt06 * np.exp(600.0 / T) + 1.90e-33 * M * kmt06 * np.exp(980.0 / T)
        p_h2o2 = np.nansum((k_ho2 * ho2 ** 2 * dt)[lit]); p_hno3 = np.nansum((_kmt08(T, M) * oh * no2 * dt)[lit])
        return float(p_h2o2 / p_hno3) if p_hno3 > 0 else float("nan")
    if key == "ope":
        o3 = g("O3")
        if o3 is None or no2 is None:
            return float("nan")
        noz = np.zeros(len(t))
        for n in _nitrogen_species(cs_path):
            if n in d:
                noz += _n_atoms(n) * np.asarray(d[n], dtype=float)
        i = np.where(lit)[0]
        if not len(i):
            return float("nan")
        d_ox = (o3 + no2)[i[-1]] - (o3 + no2)[i[0]]; d_noz = noz[i[-1]] - noz[i[0]]
        return float(d_ox / d_noz) if d_noz > 0 else float("nan")
    return float("nan")


def collect_results(st: Study, plans: List[RunPlan] = None) -> List[dict]:
    rows = []
    for p in (plans or plan_runs(st)):
        info = _read_info(p)
        have = info.get("metrics", {})
        if is_done(p) and any(k not in have for k in st.metrics):     # results added after the run: computed from its .cs once
            try:
                have = dict(have); have.update(compute_metrics(p.cs_path, [k for k in METRICS if k not in have]))
                info["metrics"] = have; _write_info(p, info)
            except Exception:
                pass
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


NOX_SPECIES = ("NO", "NO2", "NOX")


def nox_axis(axes: List[Axis]) -> int:
    """Index of the axis that varies NOx (an NO/NO2 injection or group, an initial NO/NO2, a new NO/NO2 injection); -1 if none."""
    for i, ax in enumerate(axes):
        for kind, target in ax.members():
            names = [x.strip().upper() for x in re.split(r"[,|]", target)]
            if kind in ("initial", "injection_set", "injection_scale", "injection_new") and names and \
               (names[0] in NOX_SPECIES or names[0].startswith(("NO_", "NO2_")) or all(n in NOX_SPECIES or n.startswith(("NO_", "NO2_")) for n in names if n)):
                return i
        if ax.name().upper().startswith(("NOX_", "NO_", "NO2_", "NO+")):
            return i
    return -1


def ozone_ridge(rows: List[dict], nox_axis_name: str, other_axis_name: str, key: str, log: bool) -> List[dict]:
    """Where ozone (result `key`) peaks along the NOx axis, for every value of the other axis.
    Per level: the three sampled NOx values around the maximum, a parabola through them (in log NOx
    when `log`), its vertex = ridge.  'inside' when the maximum is interior; 'above' when ozone is
    still rising at the highest NOx sampled (ridge beyond the range → the whole level is NOx-sensitive);
    'below' when it falls from the lowest NOx on (whole level VOC-sensitive).
    Returns [{"level": v, "ridge": x or None, "status": inside|above|below, "n": points}] sorted by level."""
    out = []
    for lvl in sorted({r[other_axis_name] for r in rows}):
        pts = sorted((r[nox_axis_name], r[key]) for r in rows if r[other_axis_name] == lvl and r.get(key) is not None and np.isfinite(r[key]))
        if len(pts) < 2:
            continue
        xs = np.array([p[0] for p in pts], float); zs = np.array([p[1] for p in pts], float)
        xt = np.log10(xs) if log else xs
        i = int(np.argmax(zs))
        if i == len(xs) - 1:
            out.append({"level": lvl, "ridge": None, "status": "above", "n": len(xs)}); continue
        if i == 0:
            out.append({"level": lvl, "ridge": None, "status": "below", "n": len(xs)}); continue
        a, b, c = np.polyfit(xt[i - 1:i + 2], zs[i - 1:i + 2], 2)
        xv = -b / (2 * a) if a < 0 else xt[i]
        xv = float(np.clip(xv, xt[i - 1], xt[i + 1]))
        out.append({"level": lvl, "ridge": (10 ** xv if log else xv), "status": "inside", "n": len(xs)})
    return out


def regime_map(x: np.ndarray, y: np.ndarray, z: np.ndarray, nox_on_x: bool, n: int = 80):
    """Interpolate a result onto a regular grid and its derivative along the NOx axis.
    Returns (X, Y, Z, dZ/dNOx) on the grid; NaN outside the sampled area.  The ridge (dZ/dNOx = 0)
    separates NOx-sensitive (dZ/dNOx > 0) from VOC-sensitive / NOx-saturated (dZ/dNOx < 0)."""
    from scipy.interpolate import griddata
    gx = np.linspace(np.nanmin(x), np.nanmax(x), n); gy = np.linspace(np.nanmin(y), np.nanmax(y), n)
    X, Y = np.meshgrid(gx, gy)
    Z = griddata((x, y), z, (X, Y), method="linear")
    dZdy, dZdx = np.gradient(Z, gy, gx)
    return X, Y, Z, (dZdx if nox_on_x else dZdy)


def axis_description(st: Study, ax: Axis) -> Tuple[str, float, str]:
    """(readable axis label with unit, total base ppb of the scaled injections (0 if not a scale axis), short name)
    for the plots: a factor axis says what it multiplies and what 1 means; amount axes carry their unit."""
    spec = None
    try:
        spec = load_spec(st.experiments[0]) if st.experiments else None
    except Exception:
        spec = None
    base_ppb, names = 0.0, []
    for kind, target in ax.members():
        if kind == "injection_scale" and spec is not None:
            for e in _events_matching(spec, target):
                names.append(e.name)
                try:
                    base_ppb += float(e.total_ppb) + sum(float(p) for _pr, p in e.extra_channels if _pr.strip().upper() == e.products.strip().upper())
                except ValueError:
                    pass
    k, t = ax.kind, ax.target
    if k == "injection_scale":
        what = ("injections " + ", ".join(names)) if names else f"injections {t}"
        lab = f"{what} × factor (1 = as in the experiment" + (f", {base_ppb:.3g} ppb in total)" if base_ppb else ")")
        short = "× factor"
    elif k == "injection_new":
        sp, start, dur = parse_new_injection(t); lab = f"{sp} added at {start} over {dur:g} min [ppb]"; short = f"{sp} added [ppb]"
    elif k == "injection_set":
        lab = f"{t} injected [ppb]"; short = lab
    elif k == "initial":
        lab = f"{t} at the start [ppb]"; short = lab
    elif k == "rate_factor":
        lab = f"k({t}) × factor (1 = MCM)"; short = "× factor"
    elif k == "chamber_scale":
        lab = f"chamber {t} × factor (1 = as in the experiment)"; short = "× factor"
    else:
        lab = f"{k} {t}"; short = lab
    if ax.also:
        lab += " — together with " + ", ".join(f"{kk} {tt}" for kk, tt in ax.also)
    return lab, base_ppb, short


# Quantities of the isopleth video: key → (label, unit)
TS_QUANTITIES: Dict[str, Tuple[str, str]] = {
    "O3": ("O3", "ppb"), "Ox": ("Ox = O3 + NO2", "ppb"), "NO": ("NO", "ppb"), "NO2": ("NO2", "ppb"), "NOx": ("NOx", "ppb"),
    "OH": ("OH", "cm-3"), "HO2": ("HO2", "cm-3"), "RO2": ("RO2", "cm-3"),
    "PO3": ("P(O3) instantaneous", "ppb/h"), "REGIME": ("P(H2O2)/P(HNO3) instantaneous", "-"),
}


def time_series(st: Study, plans: List[RunPlan], step_min: float = 5.0, progress: Callable = None) -> Optional[dict]:
    """Every finished run's time series on one time grid (minutes step), for the isopleth video.
    Kept in <study>/<name>.timeseries.npz and reused while the runs are unchanged.
    Returns {"hours": (nt,), "runs": [names], "axes": {axis name: (nrun,)}, "hv": (nt,), key: (nrun, nt) …}."""
    done = [p for p in plans if is_done(p)]
    if not done:
        return None
    sig = json.dumps([[p.name, round(os.path.getmtime(p.cs_path))] for p in done] + [step_min])
    cache = os.path.join(st.folder, f"{st.name}.timeseries.npz")
    if os.path.isfile(cache):
        try:
            z = np.load(cache, allow_pickle=False)
            if str(z["sig"]) == sig:
                out = {k: z[k] for k in z.files if k not in ("sig", "runs", "axis_names")}
                out["runs"] = list(z["runs"]); out["axes"] = {n: z["axis_" + n] for n in z["axis_names"]}
                return out
        except Exception:
            pass
    from experiment_plots import read_columns
    names = ["M", "T", "H2O", "HV", "O3", "NO", "NO2", "OH", "HO2", "RO2"]
    grid = None; series = {k: [] for k in TS_QUANTITIES}; hv = None
    for i, p in enumerate(done):
        if progress:
            progress(i, len(done), p.name)
        t, d = read_columns(p.cs_path, names)
        t = np.asarray(t, float)
        if grid is None:
            grid = np.arange(t[0], t[-1] + 1, step_min * 60.0)
            hv = np.interp(grid, t, np.asarray(d.get("HV", np.zeros(len(t))), float))
        on = lambda a: np.interp(grid, t, np.asarray(a, float))
        M = on(d["M"]); T = on(d["T"]) if "T" in d else np.full(len(grid), 298.0)
        ppb = lambda n: on(d[n]) / M * 1e9 if n in d else np.full(len(grid), np.nan)
        raw = lambda n: on(d[n]) if n in d else np.full(len(grid), np.nan)
        o3, no, no2, oh, ho2, ro2 = ppb("O3"), ppb("NO"), ppb("NO2"), raw("OH"), raw("HO2"), raw("RO2")
        no_c = raw("NO")
        po3 = (3.45e-12 * np.exp(270.0 / T) * ho2 * no_c + 2.7e-12 * np.exp(360.0 / T) * ro2 * no_c) / M * 1e9 * 3600.0
        h2o = raw("H2O") if "H2O" in d else 0.01 * M
        kmt06 = 1 + 1.40e-21 * np.exp(2200.0 / T) * h2o
        p_h2o2 = (2.20e-13 * kmt06 * np.exp(600.0 / T) + 1.90e-33 * M * kmt06 * np.exp(980.0 / T)) * ho2 ** 2
        p_hno3 = _kmt08(T, M) * oh * raw("NO2")
        with np.errstate(divide="ignore", invalid="ignore"):
            regime = np.where(p_hno3 > 0, p_h2o2 / p_hno3, np.nan)
        for k, v in (("O3", o3), ("Ox", o3 + no2), ("NO", no), ("NO2", no2), ("NOx", no + no2), ("OH", oh), ("HO2", ho2),
                     ("RO2", ro2), ("PO3", po3), ("REGIME", regime)):
            series[k].append(v)
    day0 = math.floor((grid[0] + 946684800.0) / 86400.0) * 86400.0 - 946684800.0
    out = {"hours": (grid - day0) / 3600.0, "hv": hv, "runs": [p.name for p in done],
           "axes": {ax.name(): np.array([p.point.get(ax.name(), np.nan) for p in done], float) for ax in st.axes}}
    out.update({k: np.vstack(v) for k, v in series.items()})
    try:
        np.savez_compressed(cache, sig=sig, runs=np.array(out["runs"]), axis_names=np.array(list(out["axes"])),
                            **{"axis_" + n: v for n, v in out["axes"].items()},
                            **{k: v for k, v in out.items() if k not in ("runs", "axes")})
    except Exception:
        pass
    return out


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
