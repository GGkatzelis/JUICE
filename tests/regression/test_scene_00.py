"""F. Rohrer's clock scenario (validation/scene_00_clock): settings that depend on a species must be live."""
import contextlib
import io
import os
import shutil
import subprocess
import sys

import numpy as np
import pytest

from juice._root import ROOT
from juice.core import scenario_prep as sp
from juice.core.column_file_reader import load_column_file

SRC = os.path.join(ROOT, "validation", "scene_00_clock")
TOML = os.path.join(ROOT, "settings.toml")
T_END, K, T0, A_T = 86400.0, 1e-5, 280.0, 100.0 / 86400.0


def closed_form():
    """(A) of the header of scene_00.easy: zeit, T, X, XX at 86 400 s."""
    from scipy.special import lambertw
    L = T0 / A_T
    u = float(np.real(-lambertw(-2 * np.exp(-(2 + T_END / L)), -1) / 2))
    x = 1 / (1 + 2 * K * L * (2 * np.log(u) + 1 / u - 1))
    return {"ZEIT": L * (u - 1), "T": T0 * u, "X": x / u, "XX": (1 - x) / (2 * u)}


def run(tmp, profile, live):
    d = os.path.join(str(tmp), f"{profile}_{live}")
    os.makedirs(d)
    for ext in (".easy", ".enz"):
        shutil.copy2(os.path.join(SRC, "scene_00" + ext), d)
    with contextlib.redirect_stdout(io.StringIO()):
        js = sp.prepare_scenario({"easy": os.path.join(d, "scene_00.easy"), "enz": os.path.join(d, "scene_00.enz")},
                                 sp.load_paths(TOML)["include_dir"], export_kv_pd=False, toml_path=TOML,
                                 solver={**sp.solver_from_toml(TOML, profile), "live_settings": live})
    subprocess.run([sys.executable, "-m", "juice.core.scenario_runner", js, "--toml", TOML], capture_output=True, check=True)
    df = load_column_file(os.path.join(d, "scene_00.cs"))
    df = df[0] if isinstance(df, tuple) else df
    return {c.upper(): float(df[c].values[-1]) for c in df.columns}


def test_closed_form_is_the_solution_of_the_equations():
    from scipy.integrate import solve_ivp

    def rhs(_t, y):
        z, x, xx = y
        dln = -A_T / (T0 + 2 * A_T * z)                 # d ln M / dt, with the clock itself diluted
        return [1.0 + z * dln, -2 * K * x * x + x * dln, K * x * x + xx * dln]
    y = solve_ivp(rhs, (0, T_END), [0.0, 1.0, 0.0], method="LSODA", rtol=1e-11, atol=1e-14).y[:, -1]
    a = closed_form()
    assert y[0] == pytest.approx(a["ZEIT"], rel=1e-7) and y[1] == pytest.approx(a["X"], rel=1e-7)
    assert a["X"] == pytest.approx(0.3046754, abs=1e-7) and a["XX"] == pytest.approx(0.2276622, abs=1e-7)


@pytest.mark.parametrize("profile,tol", [("twostep", 1e-4), ("fast", 2e-3)])
def test_live_settings_give_the_correct_solution(tmp_path, profile, tol):
    r, a = run(tmp_path, profile, True), closed_form()
    for key in ("ZEIT", "T", "X", "XX"):
        assert r[key] == pytest.approx(a[key], rel=tol), (profile, key)   # fast: atol_floor 100 >> species of order 1


def test_frozen_settings_give_the_uncorrected_solution(tmp_path):
    r = run(tmp_path, "twostep", False)
    assert r["ZEIT"] == pytest.approx(86400.0) and r["T"] == pytest.approx(380.0)
    assert r["X"] == pytest.approx(1 / (1 + 2 * K * T_END), rel=1e-4)
