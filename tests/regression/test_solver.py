"""Solver side: JUICE against FACSIMILE, and the UV-lamp chamber check.

Marked slow — they run the solver.  Quick suite only:  pytest -m "not slow"
The full nine-scene comparison is `python -m juice.tools.run_benchmark` (minutes); here one scene keeps the
FACSIMILE agreement under test on every change.
"""
import contextlib
import io
import os

import pytest

from conftest import IDL_INPUT


def test_compare_reports_missing_result(tmp_path):
    """No JUICE result → a clear FileNotFoundError naming the file (not a silent pass)."""
    from juice.tools import compare_juice_facsimile as cjf
    with contextlib.redirect_stdout(io.StringIO()), pytest.raises(FileNotFoundError, match="JUICE .cs not found"):
        cjf.compare("scene_04", 0.02, juice_dir_override=str(tmp_path),
                    fac_dir_override=os.path.join(IDL_INPUT, "scene_04"))


# Concentration floor for this test: 1e-3 ppt (≈ 2.5e4 molec/cm³), a thousand times below the 1 ppt
# 'trace' tier.  The benchmark's own default (100 molec/cm³ = the solver's atol) is kept for detailed
# comparisons, but at that level species that stay at a few molecules/cm³ are rounding noise and
# depend on the CPU: SA in scene_04 (mean ≈ 3 molec/cm³, ~1e-7 ppt) is 0.03 % off FACSIMILE on one
# machine and 3.3 % on the GitHub runner — a difference of ~3 molecules/cm³.  Settings (J, k) are
# always compared without a floor.
TEST_FLOOR = 2.5e4


@pytest.mark.slow
def test_scene_04_matches_facsimile(tmp_path):
    """scene_04 (inorganic + CH4 chemistry, SAPHIR day): every setting within 2 %, and every species
    that rises above 1e-3 ppt within 2 %."""
    from juice.tools import run_benchmark
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        rc = run_benchmark.main(["scene_04", "--work", str(tmp_path), "--floor", str(TEST_FLOOR)])
    assert rc == 0, out.getvalue()


@pytest.mark.slow
def test_uv_lamp_chamber():
    """Generic UV-lamp chamber: photostationary state, OH and Ox under the lamps, dark chemistry,
    and the pinned reference run."""
    from juice.tools import check_lamp_chamber
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        rc = check_lamp_chamber.main([])
    assert rc == 0, out.getvalue()
