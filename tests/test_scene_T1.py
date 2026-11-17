"""F. Rohrer's scene_T1 (validation/scene_T1): an ENZ file that starts at t = 0 with a temperature step at
3600 s. k = c*X follows X, and the density step at 3600 s scales every concentration by 280/300.  Before
17 Nov 2026 the ENZ rows were read one line late when the times start near 0, the step was never applied
and X ended 3.7 % high."""
import contextlib
import io
import math
import os
import shutil
import subprocess
import sys

import pytest

from juice._root import ROOT
from juice.core import scenario_prep as sp
from juice.core.column_file_reader import load_column_file

SRC = os.path.join(ROOT, "validation", "scene_T1")
TOML = os.path.join(ROOT, "settings.toml")
C, X0 = 19375.0 / (421.0 * 14400.0 * 1e18), 1e9


def exact():
    """X, XX at 3600 s (before the density step) and at 7200 s, from the header of scene_T1.easy."""
    a, r = 4 * C * X0 ** 2 * 3600, 280.0 / 300.0
    x1 = X0 / math.sqrt(1 + a)
    xx1 = (X0 - x1) / 2
    x2 = r * X0 / math.sqrt(1 + a + a * r * r)
    xx2 = r * xx1 + (r * x1 - x2) / 2
    return (x1, xx1), (x2, xx2)


def run(tmp, profile):
    d = os.path.join(str(tmp), profile)
    os.makedirs(d)
    for ext in (".easy", ".enz"):
        shutil.copy2(os.path.join(SRC, "scene_T1" + ext), d)
    with contextlib.redirect_stdout(io.StringIO()):
        js = sp.prepare_scenario({"easy": os.path.join(d, "scene_T1.easy"), "enz": os.path.join(d, "scene_T1.enz")},
                                 sp.load_paths(TOML)["include_dir"], export_kv_pd=False, toml_path=TOML,
                                 solver=sp.solver_from_toml(TOML, profile))
    subprocess.run([sys.executable, "-m", "juice.core.scenario_runner", js, "--toml", TOML], capture_output=True, check=True)
    df = load_column_file(os.path.join(d, "scene_T1.cs"))
    df = df[0] if isinstance(df, tuple) else df
    cols = {c.upper(): c for c in df.columns}
    return {k: df[cols[k]].to_numpy(dtype=float) for k in ("X", "XX", "T")}


def test_the_header_values():
    (x1, xx1), (x2, xx2) = exact()
    assert x1 == pytest.approx(1.458318e8, rel=1e-6) and xx1 == pytest.approx(4.270841e8, rel=1e-6)
    assert x2 == pytest.approx(1.0e8, rel=1e-6) and xx2 == pytest.approx(4.166667e8, rel=1e-6)


@pytest.mark.parametrize("profile,tol", [("precise", 1e-4), ("fast", 5e-4), ("twostep", 1e-3)])
def test_the_temperature_step_is_applied(tmp_path, profile, tol):
    r, ((x1, xx1), (x2, xx2)) = run(tmp_path, profile), exact()
    assert list(r["T"]) == [280.0, 280.0, 300.0]
    assert r["X"][1] == pytest.approx(x1, rel=tol) and r["XX"][1] == pytest.approx(xx1, rel=tol)
    assert r["X"][2] == pytest.approx(x2, rel=tol) and r["XX"][2] == pytest.approx(xx2, rel=tol)
