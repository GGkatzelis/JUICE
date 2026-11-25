"""The exact-answer experiment (tests/exact_small.py) under every solver option a user can set in settings.toml or
for one run: each must still give the closed-form answer, and the multi-core reaction kernels must give exactly
the result of one core.  The options are the keys of the [solver] section that the scenario file carries."""
import numpy as np
import pytest

from exact_small import prepare, run, worst_deviation

pytestmark = pytest.mark.slow

UNTOUCHED = list("ABCDEFG")      # pieces that do not depend on a setting from the ENZ file (light, flow, temperature)


@pytest.mark.parametrize("options, tol", [
    ({"jac_mode": "fd"}, 6e-3),                          # finite-difference Jacobian instead of the analytic one
    ({"ode_backend": "python"}, 6e-3),                   # the plain Python kernels instead of the compiled ones
    ({"perstep_snap": False}, 6e-3),                     # continuous mode: settings not held per row
    ({"live_settings": False}, 6e-3),                    # settings held per output step, as EASY up to 6.00
    ({"atol_floor": 1.0}, 6e-3),                         # a lower absolute tolerance floor
    # the fast tolerance with a cold restart at every row (no profile uses it; a settings file can): the restarts
    # add up on the species that decays 250-fold, 0.85 % against 0.46 % with warm continuation (25 Nov 2026)
    ({"perstep_warm": False}, 1e-2),
])
def test_every_option_still_gives_the_closed_form(tmp_path, options, tol):
    df = run(prepare(tmp_path, "fast", options))
    dev = worst_deviation(df)
    assert dev[0] < tol, (options, dev)


def test_linear_inputs_change_only_what_follows_the_inputs(tmp_path):
    """With 'linear' the ENZ values ramp from one row to the next instead of stepping: the light, flow and
    temperature switches take a minute instead of an instant.  The pieces without an ENZ setting are unchanged;
    the others move by about that minute (the temperature-step decay most: its rate is seven times higher after)."""
    df = run(prepare(tmp_path, "fast", {"interpolation": "linear"}))
    assert worst_deviation(df, UNTOUCHED)[0] < 2e-3
    dev = worst_deviation(df)
    assert dev[0] < 0.05, dev


@pytest.mark.parametrize("profile", ["fast", "twostep"])
def test_multi_core_kernels_give_the_result_of_one_core(tmp_path, profile):
    """parallel_kernels.py runs automatically from ode_threads_min_reactions (20 000) reactions on, i.e. in every
    GECKO-A run; lowering the threshold to 1 runs it on the seven reactions here."""
    one, log1 = run(prepare(tmp_path / "one", profile, {"ode_threads": 1}), with_log=True)
    two, log2 = run(prepare(tmp_path / "two", profile, {"ode_threads": 2, "ode_threads_min_reactions": 1}), with_log=True)
    assert "reaction kernels on 2 cores" in log2 and "reaction kernels on 2 cores" not in log1
    assert list(one.columns) == list(two.columns) and np.array_equal(one.index, two.index)
    assert np.array_equal(one.to_numpy(), two.to_numpy())               # bit-identical, not just close
