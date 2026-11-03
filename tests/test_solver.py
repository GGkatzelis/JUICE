"""Solver side: JUICE against FACSIMILE, and the UV-lamp chamber check.

Marked slow — they run the solver.  Quick suite only:  pytest -m "not slow"
The full nine-scene comparison is `python run_benchmark.py` (minutes); here one scene keeps the
FACSIMILE agreement under test on every change.
"""
import contextlib
import io
import os

import pytest

from conftest import IDL_INPUT


def test_compare_reports_missing_result(tmp_path):
    """No JUICE result → a clear FileNotFoundError naming the file (not a silent pass)."""
    import compare_juice_facsimile as cjf
    with contextlib.redirect_stdout(io.StringIO()), pytest.raises(FileNotFoundError, match="JUICE .cs not found"):
        cjf.compare("scene_04", 0.02, juice_dir_override=str(tmp_path),
                    fac_dir_override=os.path.join(IDL_INPUT, "scene_04"))


@pytest.mark.slow
def test_scene_04_matches_facsimile(tmp_path):
    """scene_04 (inorganic + CH4 chemistry, SAPHIR day): every setting and species within 2 %."""
    import run_benchmark
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        rc = run_benchmark.main(["scene_04", "--work", str(tmp_path)])
    assert rc == 0, out.getvalue()


@pytest.mark.slow
def test_uv_lamp_chamber():
    """Generic UV-lamp chamber: photostationary state, OH and Ox under the lamps, dark chemistry,
    and the pinned reference run."""
    import check_lamp_chamber
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        rc = check_lamp_chamber.main([])
    assert rc == 0, out.getvalue()
