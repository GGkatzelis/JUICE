"""TWOSTEP (Verwer 1994) as the third solver profile: the integrator against a tight scipy reference on the
mini mechanism, held species in production/loss form, and what the profile resolves to."""
import os

import numpy as np
import pytest

import scenario_prep as sp
import twostep_solver as ts

TOML = os.path.join(sp.HERE, "FR_GUI.toml")


def _reference(sysm, y0, rows, env):
    from scipy.integrate import solve_ivp
    y, out = y0.copy(), []
    for ta, tb, hv in rows:
        env[0] = hv
        sol = solve_ivp(sysm.rhs, (ta, tb), y, method="BDF", rtol=1e-9, atol=1e-3, max_step=60.0)
        y = np.maximum(sol.y[:, -1], ts.SMALLC)
        out.append(y.copy())
    return out


def test_twostep_matches_a_tight_bdf_reference_on_the_mini_mechanism():
    species, k, _env_names, arrays, y0, _M = ts._mini_mechanism()
    env = np.array([1.0, 1.0])
    sysm = ts.PLSystem(len(species), k, env, arrays)
    rows = [(0.0, 1800.0, 1.0), (1800.0, 3600.0, 0.0)]          # light, then dark (roof closure)
    refs = _reference(sysm, y0, rows, env)
    y = y0.copy()
    stats = {}
    for (ta, tb, hv), yr in zip(rows, refs):
        env[0] = hv
        y = ts.twostep(sysm, y, ta, tb, dtmin=1e-6, dtmax=60.0, atol=100.0, rtol=1e-3, numit=3, stats=stats)
        mask = yr > 1e6                                         # species above ~0.04 ppt
        assert np.max(np.abs(y[mask] / yr[mask] - 1.0)) < 1e-3, species
    assert stats["naccpt"] > 100 and stats["nrejec"] == 0


def test_twostep_integrates_relative_to_the_row_start():
    """JUICE times are ~8e8 s since 2000; an absolute end test would stop a 60-s row ~1 s early."""
    species, k, _env_names, arrays, y0, _M = ts._mini_mechanism()
    env = np.array([0.0, 1.0])
    sysm = ts.PLSystem(len(species), k, env, arrays)
    t0 = 8.2e8
    y_abs = ts.twostep(sysm, y0, t0, t0 + 60.0, dtmin=1e-6, dtmax=60.0, atol=100.0, rtol=1e-4, numit=3)
    y_rel = ts.twostep(sysm, y0, 0.0, 60.0, dtmin=1e-6, dtmax=60.0, atol=100.0, rtol=1e-4, numit=3)
    assert np.allclose(y_abs, y_rel, rtol=1e-12, atol=0.0)


def test_held_species_stay_exactly_where_the_row_started():
    """The runner holds a species by zeroing its P and L: the BDF2 history combination then reproduces the
    value (a1 + a2 = 1) in the implicit-Euler start, the predictor and every step, up to rounding; the runner
    restores the row-start value afterwards, so in a run they are exact."""
    species, k, _env_names, arrays, y0, _M = ts._mini_mechanism()
    env = np.array([1.0, 1.0])
    held = np.array([species.index("NO"), species.index("O3")])

    def extra_P(t, y, P):
        P[held] = 0.0

    def extra_L(t, y, L):
        L[held] = 0.0

    sysm = ts.PLSystem(len(species), k, env, arrays, extra_P=extra_P, extra_L=extra_L)
    y = ts.twostep(sysm, y0, 0.0, 1800.0, dtmin=1e-6, dtmax=60.0, atol=100.0, rtol=1e-3, numit=3)
    assert np.allclose(y[held], y0[held], rtol=1e-12, atol=0.0)
    assert y[species.index("OH")] != y0[species.index("OH")]   # the rest did move


def test_twostep_profile_resolves_to_the_verwer_solver():
    sol = sp.solver_from_toml(TOML, "twostep")
    assert sol["solver_profile"] == "twostep" and sol["solver_method"] == "TWOSTEP"
    assert (sol["rtol"], sol["perstep_warm"], sol["twostep_numit"], sol["twostep_dtmin"]) == (1e-3, False, 3, 1e-4)
    assert (sol["twostep_numit_max"], sol["twostep_conv"]) == (10, 0.1)
    # the BDF profiles keep the TOML's method and carry no TWOSTEP keys
    for name in ("fast", "precise"):
        sol = sp.solver_from_toml(TOML, name)
        assert sol["solver_method"] == "BDF" and "twostep_numit" not in sol
        assert "solver_method" not in sp.resolve_profile({}, name)     # so the runner's --solver can fall back
    assert sp.solver_profiles(TOML)["twostep"]["solver_method"] == "TWOSTEP"
    # a TOML override of the numbers is honoured, strings included
    r = sp.resolve_profile({"profiles": {"twostep": {"twostep_numit": "4", "twostep_dtmin": "0.01", "rtol": "1e-4"}}}, "twostep")
    assert (r["twostep_numit"], r["twostep_dtmin"], r["rtol"]) == (4, 0.01, 1e-4)
