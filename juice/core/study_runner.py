"""
study_runner.py — run a JUICE study (many experiments, a variant, a sweep) from the command line.

    python -m juice.core.study_runner <name>.study.json [--workers N] [--dry] [--yes] [--force]

    --dry     write every run's input files (easy, enz, json), run nothing
    --yes     no questions: skip the pilot confirmation (for servers / overnight)
    --force   rewrite inputs even when they exist (after editing the study)

The first run is executed alone as a pilot; its time gives the estimate for the rest, which is
shown before the remaining runs start.  Finished runs are skipped, so the command can be repeated.
"""
from __future__ import annotations

import os
import sys
import time

from juice.core import study as S


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")   # Windows consoles default to cp1252
    except Exception:
        pass
    if not argv or argv[0] in ("-h", "--help"):
        print(__doc__); return 0
    path = argv[0]
    workers = int(argv[argv.index("--workers") + 1]) if "--workers" in argv else 0
    dry, yes, force = "--dry" in argv, "--yes" in argv, "--force" in argv
    st = S.load_study(path)
    what = f"{len(st.tables)} input table(s) with the field case {os.path.basename(st.field_case)}" if st.field_case else \
        f"{len(st.experiments)} experiment(s), {len(st.modifications)} modification(s)"
    print(f"study '{st.name}' in {st.folder}: {what}, {len(st.axes)} axis/axes → {len(S.plan_runs(st))} runs")
    msgs = S.validate_study(st)
    for lvl, m in msgs:
        print(f"  [{lvl}] {m}")
    if any(lvl == "error" for lvl, _ in msgs):
        print("fix the errors first"); return 2
    r = S.StudyRunner(st)
    t0 = time.monotonic()
    errs = r.materialise_all(lambda k, p, _r: print(f"  inputs: {p.name}" + (f"  ({p.message})" if p.message else "")), force=force)
    print(f"inputs written in {S.fmt_duration(time.monotonic() - t0)}" + (f"; {len(errs)} failed" if errs else ""))
    if dry:
        return 1 if errs else 0

    def ev(kind, p, rr):
        if kind == "started":
            print(f"  ▶ {p.name}")
        elif kind in ("finished", "skipped"):
            print(f"  ✔ {p.name}  {S.fmt_duration(p.elapsed)}" + ("  (already done)" if kind == "skipped" else ""))
        elif kind == "failed":
            print(f"  ✖ {p.name}  {p.message}")
        elif kind == "all_done":
            pass
    todo = [p for p in r.plans if not S.is_done(p)]
    if todo:
        pilot = todo[0]
        print(f"pilot run: {pilot.name}")
        r.run_one(pilot, ev)
        rest = [p for p in r.plans if not S.is_done(p)]
        if rest:
            secs, txt = r.estimate(workers)
            wk = workers or st.workers
            print(f"{len(rest)} runs remaining; pilot {S.fmt_duration(r.per_run_seconds())} alone, about "
                  f"{S.parallel_factor(min(wk, len(rest))):.1f}× slower with {wk} in parallel ≈ {txt}")
            lvl, gtxt = S.check_guards(len(rest))
            if not yes:
                ans = input(f"{gtxt}. Continue? [y/N] ").strip().lower()
                if ans not in ("y", "yes"):
                    print("stopped after the pilot; rerun to continue"); return 0
            r.run_all(ev, workers, rest)
    rows = S.collect_results(st, r.plans)
    out = S.write_results(st, rows)
    ok = sum(1 for x in rows if x["status"] == "done")
    print(f"{ok}/{len(rows)} runs done in {S.fmt_duration(time.monotonic() - t0)}; results: {out}")
    return 0 if ok == len(rows) else 1


if __name__ == "__main__":
    sys.exit(main())
