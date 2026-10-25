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


def solver_from_toml(toml_path: str) -> dict:
    """The solver keys scenario_runner reads from the scenario JSON (not from the TOML),
    taken from the TOML's [solver] section so every run uses the validated set.
    Fallbacks are the validated values, never the runner's own (python / linear / 1e-2)."""
    sol: dict = {}
    for p in (toml_path, toml_path[:-5] + ".local.toml" if toml_path.endswith(".toml") else ""):
        if p and os.path.isfile(p):
            try:
                with open(p, "rb") as fh:
                    sol.update(tomllib.load(fh).get("solver", {}))
            except Exception as e:
                print(f"[scenario_prep] could not read {p}: {e}")

    def _f(key, default):
        try:
            return float(sol.get(key, default))
        except (TypeError, ValueError):
            return default
    return {"rtol": _f("rtol_default", 1e-4),
            "ode_backend": str(sol.get("ode_backend", "numba")),
            "solver_method": str(sol.get("solver_method", "BDF")),
            "jac_mode": str(sol.get("jac_mode", "analytic")),
            "interpolation": str(sol.get("interpolation", "piecewise")),
            "first_step": _f("first_step", 1e-6)}


def prepare_scenario(res: dict, include_dir: str, export_kv_pd: bool = True, toml_path: str = "",
                     extra_dirs: Iterable[str] = (), solver: dict | None = None) -> str:
    """EASY + ENZ (from write_experiment) → canonical elog → scenario JSON next to them.
    `extra_dirs` are searched for include files after the run folder and before include_dir
    (a study run uses its base experiment's chemistry folder).  `solver` replaces the TOML's
    [solver] keys when given (a study freezes them once for all its runs).  Returns the JSON path."""
    from elog_generator import generate_both
    from idl_parser import parse_idl_easy
    out = os.path.dirname(res["easy"])
    dirs = [out, os.path.join(out, "chemistry")] + [d for d in extra_dirs if d] + [include_dir]
    raw, canon = generate_both(res["easy"], search_dirs=dirs, output_dir=out)
    cfg = parse_idl_easy(canon, enz_path=res["enz"], chemistry_folder=out, include_fallbacks=dirs[1:])
    cfg.update({"folder": out, "enz_path": res["enz"], "elog_path": canon, "export_kv_pd": bool(export_kv_pd)})
    cfg.update(solver if solver is not None else solver_from_toml(toml_path))
    name = os.path.splitext(os.path.basename(res["easy"]))[0]
    path = os.path.join(out, name + ".json")
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(cfg, fh, indent=1, default=str)
    return path
