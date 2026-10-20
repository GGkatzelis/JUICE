"""
check_easy_roundtrip.py — regression check for easy_scene.py (Goal 2: EASY writer)
=================================================================================

For every scene under MCM_tested_scenarios/IDL Input/scene_XX:

  1. document level : parse + render must be byte-identical to the source file
  2. scene level    : extract_scene + write_scene must be statement-equivalent
                      (same active statements, same last-wins definitions)
  3. pipeline level : elog_generator + idl_parser on the original and on the
                      regenerated file must give identical scenario JSON

Usage:
    python check_easy_roundtrip.py                # all scenes
    python check_easy_roundtrip.py scene_05 scene_09
Regenerated files and elogs are written to a temporary directory.
"""

import sys
import os
import json
import shutil
import io
import tempfile
import contextlib

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
os.chdir(HERE)

from easy_scene import roundtrip_report, statements_equivalent          # noqa: E402
from elog_generator import generate_both                                # noqa: E402
from idl_parser import parse_idl_easy                                   # noqa: E402

IDL_INPUT = os.path.join(HERE, "MCM_tested_scenarios", "IDL Input")
INC = os.path.join(IDL_INPUT, "include files")
SKIP = {"folder", "enz_path", "elog_path", "chemistry_folder", "mechanism_files",
        "skipped_notes", "parser_notes"}


def _pipeline(easy_path, work, search):
    os.makedirs(work, exist_ok=True)
    with contextlib.redirect_stdout(io.StringIO()):
        _raw, canon = generate_both(easy_path, search_dirs=search, output_dir=work)
        d = os.path.dirname(easy_path)
        enz = next((os.path.join(d, f) for f in os.listdir(d) if f.lower().endswith(".enz")), None)
        cfg = parse_idl_easy(canon, enz_path=enz, chemistry_folder=d, include_fallbacks=[INC])
    return canon, {k: v for k, v in cfg.items() if k not in SKIP}


def check(scene, tmp):
    orig = os.path.join(IDL_INPUT, scene, scene + ".easy")
    if not os.path.isfile(orig):
        return f"{scene}: no .easy file"
    # Work on copies: elog_generator removes stale *_juice*.elog files next to
    # the .easy it processes, which must never touch the committed scene folders.
    rdir = os.path.join(tmp, scene)
    odir = os.path.join(rdir, "orig_src")
    os.makedirs(odir, exist_ok=True)
    rep = roundtrip_report(orig, os.path.join(rdir, scene + ".regen.easy"))
    orig_copy = os.path.join(odir, scene + ".easy")
    shutil.copy(orig, orig_copy)
    regen = os.path.join(rdir, scene + ".easy")
    shutil.copy(rep["regenerated"], regen)
    for f in os.listdir(os.path.dirname(orig)):
        if f.lower().endswith(".enz"):
            shutil.copy(os.path.join(os.path.dirname(orig), f), rdir)
            shutil.copy(os.path.join(os.path.dirname(orig), f), odir)
    search = [os.path.dirname(orig), INC]
    c1, j1 = _pipeline(orig_copy, os.path.join(rdir, "orig"), search)
    c2, j2 = _pipeline(regen, os.path.join(rdir, "regen"), search)
    ok_e, d_e = statements_equivalent(c1, c2)
    kd = [k for k in sorted(set(j1) | set(j2))
          if json.dumps(j1.get(k), sort_keys=True, default=str) != json.dumps(j2.get(k), sort_keys=True, default=str)]
    ok = rep["level1_byte_identical"] and rep["level2_equivalent"] and ok_e and not kd
    line = (f"{'OK  ' if ok else 'FAIL'} {scene}: byte-identical={rep['level1_byte_identical']} "
            f"statement-equivalent={rep['level2_equivalent']} elog-equivalent={ok_e} "
            f"json-identical={not kd} unstructured={rep['unstructured']}")
    for d in (rep["diffs"] + d_e)[:5]:
        line += "\n       " + d[:150]
    if kd:
        line += "\n       JSON keys differing: " + ", ".join(kd)
    return line


if __name__ == "__main__":
    scenes = sys.argv[1:] or sorted(d for d in os.listdir(IDL_INPUT)
                                     if d.startswith("scene_") and os.path.isdir(os.path.join(IDL_INPUT, d)))
    tmp = tempfile.mkdtemp(prefix="easy_roundtrip_")
    bad = 0
    for sc in scenes:
        msg = check(sc, tmp)
        bad += msg.startswith("FAIL")
        print(msg, flush=True)
    print(f"\n{len(scenes) - bad}/{len(scenes)} scenes pass  (work dir: {tmp})")
    sys.exit(1 if bad else 0)
