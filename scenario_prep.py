"""
scenario_prep.py — from the written experiment files to the scenario JSON the runner reads (no Qt).

Shared by the step-by-step GUI (juice_gui), the study runner (many runs) and the command line:

    paths  = load_paths(toml_path)                 # [paths] of the settings TOML (+ .local overlay)
    solver = solver_from_toml(toml_path)           # the validated [solver] keys the runner expects in the JSON
    json   = prepare_scenario(res, include_dir, export_kv_pd, toml_path, extra_dirs=[...])

`res` is write_experiment's result dict ('easy', 'enz').  The JSON stem = the EASY stem, so the
runner writes <stem>.cs next to it.
"""
from __future__ import annotations

import json
import os
import tomllib
from typing import Dict, Iterable

HERE = os.path.dirname(os.path.abspath(__file__))


def load_paths(toml_path: str) -> Dict[str, str]:
    """[paths] of the TOML (+ .local.toml overlay); defaults relative to this file."""
    base = {"easy_root": os.path.join(HERE, "MCM_tested_scenarios", "IDL Input"),
            "include_dir": os.path.join(HERE, "MCM_tested_scenarios", "IDL Input", "include files"),
            "output_root": os.path.join(HERE, "MCM_tested_scenarios", "JUICE archive", "JUICE"),
            "last_easy_dir": ""}
    for p in (toml_path, toml_path[:-5] + ".local.toml" if toml_path.endswith(".toml") else ""):
        if p and os.path.isfile(p):
            try:
                with open(p, "rb") as fh:
                    base.update({k: str(v) for k, v in tomllib.load(fh).get("paths", {}).items()})
            except Exception as e:
                print(f"[scenario_prep] could not read {p}: {e}")
    return base


def pick_toml(argv) -> str:
    """-t/--toml <file> from argv, else GG_GUI.toml, else FR_GUI.toml next to this file, else ''."""
    argv = list(argv)
    for i, a in enumerate(argv):
        if a in ("-t", "--toml") and i + 1 < len(argv):
            p = argv[i + 1]
            return p if os.path.isabs(p) else os.path.join(HERE, p)
    for name in ("GG_GUI.toml", "FR_GUI.toml"):
        if os.path.isfile(os.path.join(HERE, name)):
            return os.path.join(HERE, name)
    return ""


# ── Solver profiles ──────────────────────────────────────────────────────────
# Two validated ways to run the same chemistry.  Both agree with FACSIMILE far inside the 2 % benchmark
# tolerance for every species that reaches 1 ppt (benchmark 2026-11-04, scenes 04-09); they differ in
# run time and in how tightly they track FACSIMILE.  The TOML may override the numbers in
# [solver.profiles.<name>]; the names and these defaults are the reference.
SOLVER_PROFILES: Dict[str, dict] = {
    "fast": {
        "rtol": 1e-4, "perstep_warm": True,
        "label": "Fast (rtol 1e-4, warm continuation)",
        "summary": "About half the run time of 'precise'. Every species above 1 ppt within 0.6 % of FACSIMILE, "
                   "O3, NOx, OH and HO2 within 0.25 %; trace species below 1e-6 ppt at most 7 % off. "
                   "For O3, SOA and product studies, sweeps and isopleths.",
    },
    "precise": {
        "rtol": 1e-6, "perstep_warm": False,
        "label": "Precise (rtol 1e-6, cold restart per ENZ row)",
        "summary": "About twice the run time. In large mechanisms (scene_06-09) O3, NOx, OH and HO2 within 0.03 % "
                   "of FACSIMILE (5-10x tighter than 'fast'), every species above 1 ppt within 0.35 % (scene_09: 0.06 %). "
                   "Rare large relative excursions (up to 5x) in species below 0.01 ppt. "
                   "For validation against FACSIMILE and reference runs.",
    },
}
DEFAULT_PROFILE = "fast"


def _tobool(v) -> bool:
    return v if isinstance(v, bool) else str(v).strip().lower() in ("1", "true", "yes", "on")


def _read_solver_section(toml_path: str) -> dict:
    sol: dict = {}
    for p in (toml_path, toml_path[:-5] + ".local.toml" if toml_path and toml_path.endswith(".toml") else ""):
        if p and os.path.isfile(p):
            try:
                with open(p, "rb") as fh:
                    sol.update(tomllib.load(fh).get("solver", {}))
            except Exception as e:
                print(f"[scenario_prep] could not read {p}: {e}")
    return sol


def resolve_profile(sol: dict, profile: str | None = None) -> dict:
    """{'solver_profile', 'rtol', 'perstep_warm'} from a [solver] section.

    `profile` (or the section's `profile` key) names one of SOLVER_PROFILES, with the numbers of
    [solver.profiles.<name>] on top.  A section without a `profile` key and with rtol_default is a
    settings file from before profiles existed: its rtol_default / perstep_warm are used as written
    (profile 'custom'), so an old TOML reproduces exactly what it did."""
    name = (profile or sol.get("profile") or "").strip().lower()
    if not name and "rtol_default" in sol:
        return {"solver_profile": "custom", "rtol": float(sol["rtol_default"]),
                "perstep_warm": _tobool(sol.get("perstep_warm", False))}
    name = name or DEFAULT_PROFILE
    if name not in SOLVER_PROFILES:
        raise ValueError(f"unknown solver profile '{name}' (choose {', '.join(SOLVER_PROFILES)})")
    p = dict(SOLVER_PROFILES[name])
    p.update((sol.get("profiles", {}) or {}).get(name, {}))
    return {"solver_profile": name, "rtol": float(p["rtol"]), "perstep_warm": _tobool(p["perstep_warm"])}


def solver_profiles(toml_path: str = "") -> Dict[str, dict]:
    """The profiles with the TOML's numbers applied: {name: {rtol, perstep_warm, label, summary}}."""
    sol = _read_solver_section(toml_path)
    out = {}
    for name, p in SOLVER_PROFILES.items():
        out[name] = {**p, **resolve_profile({"profiles": sol.get("profiles", {})}, name)}
    return out


def default_profile(toml_path: str = "") -> str:
    """The profile a run uses when nobody chooses: [solver] profile, else 'custom' for a pre-profile TOML."""
    return resolve_profile(_read_solver_section(toml_path))["solver_profile"]


def solver_from_toml(toml_path: str, profile: str | None = None) -> dict:
    """The solver keys scenario_runner reads from the scenario JSON (not from the TOML),
    taken from the TOML's [solver] section so every run uses a validated set.
    `profile` ('fast' | 'precise') overrides the TOML's [solver] profile.
    Fallbacks are the validated values, never the runner's own (python / linear / 1e-2)."""
    sol = _read_solver_section(toml_path)

    def _f(key, default):
        try:
            return float(sol.get(key, default))
        except (TypeError, ValueError):
            return default
    return {**resolve_profile(sol, profile),
            "ode_backend": str(sol.get("ode_backend", "numba")),
            "solver_method": str(sol.get("solver_method", "BDF")),
            "jac_mode": str(sol.get("jac_mode", "analytic")),
            "interpolation": str(sol.get("interpolation", "piecewise")),
            "first_step": _f("first_step", 1e-6)}


def prepare_scenario(res: dict, include_dir: str, export_kv_pd: bool = True, toml_path: str = "",
                     extra_dirs: Iterable[str] = (), solver: dict | None = None,
                     easy_version: str | None = None) -> str:
    """EASY + ENZ (from write_experiment) → canonical elog → scenario JSON next to them.
    `extra_dirs` are searched for include files after the run folder and before include_dir
    (a study run uses its base experiment's chemistry folder).  `solver` replaces the TOML's
    [solver] keys when given (a study freezes them once for all its runs).
    `easy_version` is the EASY that made the FACSIMILE run to reproduce ('5.96b', '6.00'); it sets
    how redefined FUNCTION macros apply.  None = read it from a .fac/.out in the run folder or
    extra_dirs; '' or not found = EASY 6.00 rules (JUICE's own experiments).  Returns the JSON path."""
    from elog_generator import generate_both, easy_version as _easy_version, function_scope_for
    from idl_parser import parse_idl_easy
    out = os.path.dirname(res["easy"])
    extra = [d for d in extra_dirs if d]
    dirs = [out, os.path.join(out, "chemistry")] + extra + [include_dir]
    if easy_version is None:
        stem = os.path.splitext(os.path.basename(res["easy"]))[0]
        easy_version = _easy_version(out, *extra, stem=stem)
    scope = function_scope_for(easy_version)
    print(f"[scenario_prep] EASY version {easy_version or '(none: JUICE rules = 6.00)'}: "
          f"FUNCTION scope '{scope}'")
    raw, canon = generate_both(res["easy"], search_dirs=dirs, output_dir=out, function_scope=scope)
    cfg = parse_idl_easy(canon, enz_path=res["enz"], chemistry_folder=out, include_fallbacks=dirs[1:])
    cfg.update({"folder": out, "enz_path": res["enz"], "elog_path": canon, "export_kv_pd": bool(export_kv_pd),
                "easy_version": easy_version or "", "function_scope": scope})
    cfg.update(solver if solver is not None else solver_from_toml(toml_path))
    name = os.path.splitext(os.path.basename(res["easy"]))[0]
    path = os.path.join(out, name + ".json")
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(cfg, fh, indent=1, default=str)
    return path
