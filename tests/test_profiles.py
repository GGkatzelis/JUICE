"""Solver profiles (fast | precise): what each resolves to, and that old settings files keep their meaning."""
import os

import pytest

import scenario_prep as sp

TOML = os.path.join(sp.HERE, "FR_GUI.toml")


def test_repository_toml_default_is_fast():
    sol = sp.solver_from_toml(TOML)
    assert (sol["solver_profile"], sol["rtol"], sol["perstep_warm"]) == ("fast", 1e-4, True)


def test_precise_profile():
    sol = sp.solver_from_toml(TOML, "precise")
    assert (sol["solver_profile"], sol["rtol"], sol["perstep_warm"]) == ("precise", 1e-6, False)


def test_toml_numbers_match_the_reference_profiles():
    """FR_GUI.toml [solver.profiles.*] must not drift from the benchmarked SOLVER_PROFILES."""
    for name, p in sp.solver_profiles(TOML).items():
        ref = sp.SOLVER_PROFILES[name]
        assert (p["rtol"], p["perstep_warm"]) == (ref["rtol"], ref["perstep_warm"]), name


def test_never_rtol_1e6_with_warm_continuation():
    """rtol 1e-6 + warm continuation diverges in trace species (HCOCO3 x1400, scene_06)."""
    for p in sp.solver_profiles(TOML).values():
        assert not (p["rtol"] < 1e-5 and p["perstep_warm"])


def test_pre_profile_toml_is_used_as_written():
    """A settings file from before profiles (rtol_default, no profile key) reproduces what it did."""
    assert sp.resolve_profile({"rtol_default": "1e-06"}) == \
        {"solver_profile": "custom", "rtol": 1e-6, "perstep_warm": False}
    assert sp.resolve_profile({"rtol_default": "1e-04", "perstep_warm": "true"}) == \
        {"solver_profile": "custom", "rtol": 1e-4, "perstep_warm": True}


def test_explicit_profile_beats_the_toml():
    assert sp.resolve_profile({"profile": "precise"}, "fast")["solver_profile"] == "fast"
    assert sp.resolve_profile({"rtol_default": "1e-06"}, "fast")["rtol"] == 1e-4


def test_unknown_profile_is_an_error():
    with pytest.raises(ValueError):
        sp.resolve_profile({}, "turbo")
