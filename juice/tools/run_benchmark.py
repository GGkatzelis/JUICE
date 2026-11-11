"""
run_benchmark.py — JUICE vs FACSIMILE on the reference scenarios.

For every scene in validation/IDL Input that has <scene>.easy, <scene>.enz and the
FACSIMILE reference <scene>.cs (scene_01 … scene_09):

    1. copy the .easy and .enz to <work>/<scene>/  (the reference folder is never written to)
    2. prepare the scenario JSON exactly as the step GUI does (scenario_prep.prepare_scenario,
       solver settings from settings.toml [solver], profile fast or precise)
    3. run scenario_runner.py in its own process (N in parallel with --workers)
    4. compare the JUICE .cs with the FACSIMILE .cs (compare_juice_facsimile.compare)

and print one table: settings and concentrations within --tol, the worst species, run time.
The same table is written to <work>/benchmark_results.md, the numbers to benchmark_results.json,
each scene's full comparison to <work>/<scene>/<scene>_difference_JUICE_Facsimile.log.

    python -m juice.tools.run_benchmark                          # all scenes, work dir = output_root/benchmark
    python -m juice.tools.run_benchmark scene_03 scene_04        # some scenes
    python -m juice.tools.run_benchmark --workers 3              # parallel solver processes
    python -m juice.tools.run_benchmark --compare-only           # re-compare existing results, no runs
    python -m juice.tools.run_benchmark --kv-pd                  # also export and compare rate constants (k, v)
    python -m juice.tools.run_benchmark --solver precise         # the slower, tighter solver profile (default: the
                                                     # TOML's [solver] profile); work dir .../benchmark/<profile>

Exit code 0 when every scene's settings and concentrations are within --tol (default 2 %).

`parse_cs_file` is kept for check_lamp_chamber.py.
"""
from __future__ import annotations
from juice._root import ROOT as _JUICE_ROOT

import argparse
import contextlib
import datetime as _dt
import io
import json
import os
import re
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor

import pandas as pd

HERE = _JUICE_ROOT
REF_ROOT = os.path.join(HERE, "validation", "IDL Input")
TOML = os.path.join(HERE, "settings.toml")
PROFILE: str | None = None          # --solver: fast | precise (None = the TOML's [solver] profile)
PROFILE_DESC = ""                   # the resolved profile, for the report header
S_LABEL = "1. s-values (settings)"
C_LABEL = "2. c-values (concentrations)"


def parse_cs_file(cs_path):
    """Parse a FACSIMILE .cs output file into a DataFrame (TIME + one column per output)."""
    with open(cs_path, encoding='latin-1') as f:
        lines = f.readlines()
    col_names = ['TIME']
    data_start = 0
    for i, line in enumerate(lines):
        m = re.match(r'\s*COLUMN\s+\d+=(\S+)', line)
        if m:
            name = m.group(1).split('<')[0]
            if name != 'TIME':
                col_names.append(name)
        if i > 5 and line.strip():
            try:
                vals = line.strip().split()
                if len(vals) > 5:
                    float(vals[0]); float(vals[1])
                    data_start = i
                    break
            except ValueError:
                pass
    rows = []
    for line in lines[data_start:]:
        line = re.sub(r'(?<=\d)[dD](?=[+-]?\d)', 'E', line.strip())
        if not line:
            continue
        try:
            vals = [float(x) for x in line.split()]
            if len(vals) == len(col_names):
                rows.append(vals)
        except ValueError:
            pass
    return pd.DataFrame(rows, columns=col_names)


def reference_scenes() -> list[str]:
    """Scenes with the three files the benchmark needs, in name order."""
    out = []
    for name in sorted(os.listdir(REF_ROOT)):
        d = os.path.join(REF_ROOT, name)
        if os.path.isdir(d) and all(os.path.isfile(os.path.join(d, f"{name}{ext}"))
                                    for ext in (".easy", ".enz", ".cs")):
            out.append(name)
    return out


def _worker_env(parallel: bool) -> dict:
    env = dict(os.environ)
    env.update({"PYTHONIOENCODING": "utf-8", "PYTHONUTF8": "1"})
    if parallel:      # one BLAS / numba thread per solver process, as study.py does
        env.update({"OMP_NUM_THREADS": "1", "MKL_NUM_THREADS": "1",
                    "OPENBLAS_NUM_THREADS": "1", "NUMBA_NUM_THREADS": "1"})
    return env


def prepare(scene: str, work: str, kv_pd: bool) -> str:
    """Copy the scene's inputs to the work folder and write the scenario JSON; returns its path."""
    from juice.core import scenario_prep as sp
    src = os.path.join(REF_ROOT, scene)
    dst = os.path.join(work, scene)
    os.makedirs(dst, exist_ok=True)
    for ext in (".easy", ".enz"):
        shutil.copy2(os.path.join(src, scene + ext), dst)
    res = {"easy": os.path.join(dst, scene + ".easy"), "enz": os.path.join(dst, scene + ".enz")}
    include_dir = sp.load_paths(TOML)["include_dir"]
    # the reference folder (and its chemistry/ subfolder) is searched for include files first
    return sp.prepare_scenario(res, include_dir, export_kv_pd=kv_pd, toml_path=TOML,
                               extra_dirs=[src, os.path.join(src, "chemistry")],
                               solver=sp.solver_from_toml(TOML, PROFILE))


def prepare_quietly(scene: str, work: str, kv_pd: bool) -> tuple[str, str]:
    """prepare() with its parser report written to <scene>_prepare.log; returns (json path, error)."""
    os.makedirs(os.path.join(work, scene), exist_ok=True)
    buf = io.StringIO()
    try:
        with contextlib.redirect_stdout(buf):
            path, err = prepare(scene, work, kv_pd), ""
    except Exception as e:
        path, err = "", f"preparation failed — {type(e).__name__}: {e}"
    with open(os.path.join(work, scene, f"{scene}_prepare.log"), "w", encoding="utf-8") as fh:
        fh.write(buf.getvalue())
    return path, err


def run(scene: str, work: str, json_path: str, parallel: bool) -> dict:
    """Solve one prepared scene in its own process; returns {'scene', 'seconds', 'error'}."""
    t0 = time.monotonic()
    log_path = os.path.join(work, scene, f"{scene}_run.log")
    try:
        with open(log_path, "w", encoding="utf-8") as log:
            proc = subprocess.run([sys.executable, "-u", "-m", "juice.core.scenario_runner",
                                   json_path, "--toml", TOML],
                                  cwd=HERE, env=_worker_env(parallel), stdout=log, stderr=subprocess.STDOUT)
        err = "" if proc.returncode == 0 else f"solver exit code {proc.returncode} (see {log_path})"
        if not err and not os.path.isfile(os.path.join(work, scene, scene + ".cs")):
            err = f"solver wrote no result (see {log_path})"
    except Exception as e:
        err = f"{type(e).__name__}: {e}"
    return {"scene": scene, "seconds": time.monotonic() - t0, "error": err}


def compare(scene: str, work: str, tol: float, floor: float) -> dict:
    """Run the comparison quietly (full report into the scene's log); returns its summary."""
    from juice.tools import compare_juice_facsimile as cjf
    log_path = os.path.join(work, scene, f"{scene}_difference_JUICE_Facsimile.log")
    buf = io.StringIO()
    try:
        with contextlib.redirect_stdout(buf):
            summary = cjf.compare(scene, tol, no_stop=True, floor=floor,
                                  juice_dir_override=os.path.join(work, scene),
                                  fac_dir_override=os.path.join(REF_ROOT, scene))
    except Exception as e:
        summary = {"ok": False, "groups": {}, "error": f"{type(e).__name__}: {e}"}
    with open(log_path, "w", encoding="utf-8") as fh:
        fh.write(buf.getvalue())
    return summary


def _group_cells(g: dict | None, tol: float) -> tuple[str, str]:
    if not g:
        return "—", "—"
    within = g["n"] - g["n_exceed"]
    worst = g["worst"][0] if g["worst"] else (0.0, "")
    return f"{within}/{g['n']}", (f"{worst[1]} {worst[0]*100:.2g} %" if worst[1] else "—")


def report(rows: list[dict], tol: float, floor: float) -> str:
    """The results as a Markdown table (also printed)."""
    lines = [f"JUICE vs FACSIMILE — {_dt.datetime.now():%Y-%m-%d %H:%M}, tolerance ±{tol*100:g} % "
             f"(max relative deviation over the run; concentrations measured against "
             f"max(max|FAC|, {floor:g} molec/cm³))"] + ([f"solver profile: {PROFILE_DESC}"] if PROFILE_DESC else []) + ["",
             "| scene | settings within | concentrations within | worst concentration | "
             "within 0.1 % | run time | verdict |",
             "|---|---|---|---|---|---|---|"]
    for r in rows:
        s = r["summary"]["groups"].get(S_LABEL)
        c = r["summary"]["groups"].get(C_LABEL)
        s_in, _ = _group_cells(s, tol)
        c_in, c_worst = _group_cells(c, tol)
        fine = (f"{sum(1 for rel, _ in c['all'] if rel <= 1e-3)}/{c['n']}" if c else "—")
        t = r["run"]["seconds"] if r.get("run") else None
        t_str = "—" if t is None else (f"{t:.0f} s" if t < 120 else f"{t/60:.1f} min")
        problem = (r.get("run") or {}).get("error") or r["summary"].get("error", "")
        ok = (not problem) and s is not None and c is not None and not s["n_exceed"] and not c["n_exceed"]
        verdict = "pass" if ok else ("**FAIL** " + problem if problem else "**FAIL**")
        lines.append(f"| {r['scene']} | {s_in} | {c_in} | {c_worst} | {fine} | {t_str} | {verdict} |")
    return "\n".join(lines)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="JUICE vs FACSIMILE on the reference scenarios.")
    ap.add_argument("scenes", nargs="*", help="scene names (default: every reference scene)")
    ap.add_argument("--work", default="", help="work folder (default: output_root/benchmark)")
    ap.add_argument("--workers", type=int, default=1, help="solver processes in parallel (default 1)")
    ap.add_argument("--tol", type=float, default=0.02, help="max relative deviation (default 0.02 = 2 %%)")
    ap.add_argument("--floor", type=float, default=100.0,
                    help="concentration floor, molec/cm³ (default 100 = the solver's atol_floor)")
    ap.add_argument("--compare-only", action="store_true", help="compare existing results, do not run")
    ap.add_argument("--kv-pd", action="store_true", help="export and compare k/v values too (slower)")
    ap.add_argument("--toml", default="", help="settings file for the solver (default settings.toml), "
                                               "e.g. to compare two solver sets")
    ap.add_argument("--solver", choices=["fast", "precise", "twostep"], default=None,
                    help="solver profile (default: the TOML's [solver] profile). fast = rtol 1e-4, warm; "
                         "precise = rtol 1e-6, cold restart per ENZ row, ~2x slower, tightest agreement")
    a = ap.parse_args(argv)
    global TOML, PROFILE, PROFILE_DESC
    if a.toml:
        TOML = os.path.abspath(a.toml)
    PROFILE = a.solver

    sys.path.insert(0, HERE)
    from juice.core import scenario_prep as sp
    profile = sp.solver_from_toml(TOML, PROFILE)
    PROFILE_DESC = (f"{profile['solver_profile']} ({profile['solver_method']}, rtol {profile['rtol']:g}, "
                    f"{'warm continuation' if profile['perstep_warm'] else 'cold restart per ENZ row'})")
    print(f"solver profile: {PROFILE_DESC}")
    available = reference_scenes()
    scenes = a.scenes or available
    unknown = [s for s in scenes if s not in available]
    if unknown:
        print(f"not a reference scene (needs .easy, .enz and .cs in {REF_ROOT}): {', '.join(unknown)}")
        return 2
    work = a.work or os.path.join(sp.load_paths(TOML)["output_root"], "benchmark", profile["solver_profile"])
    os.makedirs(work, exist_ok=True)
    print(f"benchmark: {len(scenes)} scene(s) → {work}")

    runs: dict[str, dict] = {}
    if not a.compare_only:
        parallel = a.workers > 1
        prepared = {}
        for s in scenes:                                   # sequential: the parsers print a lot
            json_path, err = prepare_quietly(s, work, a.kv_pd)
            if err:
                runs[s] = {"scene": s, "seconds": 0.0, "error": err}
                print(f"  {s}: {err}")
            else:
                prepared[s] = json_path
        with ThreadPoolExecutor(max_workers=max(1, a.workers)) as pool:
            futures = {s: pool.submit(run, s, work, j, parallel) for s, j in prepared.items()}
            for s in prepared:
                runs[s] = futures[s].result()
                r = runs[s]
                print(f"  {s}: {'solved' if not r['error'] else 'FAILED'} in {r['seconds']:.0f} s"
                      + (f" — {r['error']}" if r["error"] else ""))

    rows = []
    for s in scenes:
        summary = compare(s, work, a.tol, a.floor)
        rows.append({"scene": s, "run": runs.get(s), "summary": summary})

    table = report(rows, a.tol, a.floor)
    print("\n" + table + "\n")
    print(f"full comparisons: {work}\\<scene>\\<scene>_difference_JUICE_Facsimile.log")
    with open(os.path.join(work, "benchmark_results.md"), "w", encoding="utf-8") as fh:
        fh.write(table + "\n")
    with open(os.path.join(work, "benchmark_results.json"), "w", encoding="utf-8") as fh:
        json.dump(rows, fh, indent=1, default=str)

    ok = all(r["summary"].get("ok") and not (r.get("run") or {}).get("error") for r in rows)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
