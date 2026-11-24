"""A small experiment whose every result can be written down: seven independent pieces of chemistry in one EASY
file, run through the whole chain (EASY + ENZ -> parser -> scenario -> solver -> .cs) with each solver profile.

    A --> B                      first-order decay
    C <--> D                     relaxation to an equilibrium
    E + F --> G                  second order, unequal amounts
    H --> I       k = HV * 5e-4  only while the light is on (HV from the ENZ file: on from 1800 s to 5400 s)
    QJ --> J                     an injection of 8 ppb over the two rows at 600 s and 660 s
    L + DIL -->                  dilution, flow on from 3600 s
    N --> O       k = 1e-2*exp(-1000/T), T from the ENZ file: 250 K, then 500 K from 3600 s (with P following T,
                                 so that the air density and every concentration stay as they are)

Two hours, one ENZ row per minute.  This is the test of the solver and of the two parsers from outside: a wrong
row convention, a setting read one row late, a lost reaction or a wrong unit shows up as a wrong number here."""
import json
import math

import numpy as np
import pytest

from exact_small import END, M, exact, prepare, run


@pytest.fixture(scope="module")
def prepared(tmp_path_factory):
    """The scenario file of every profile, prepared once."""
    return {profile: prepare(tmp_path_factory.mktemp(profile), profile) for profile in ("precise", "fast", "twostep")}


def test_closed_forms_at_the_moments_that_matter():
    assert exact(0)["A"] == 1e9 and exact(0)["J"] == 0.0 and exact(0)["G"] == 0.0
    assert exact(END)["A"] == pytest.approx(1e9 * math.exp(-0.72))
    assert exact(1e6)["C"] == pytest.approx(1e9) and exact(1e6)["F"] == pytest.approx(0.0, abs=1.0)
    assert exact(1800)["H"] == 1e9 and exact(5400)["H"] == exact(END)["H"] == pytest.approx(1e9 * math.exp(-1.8))
    assert exact(600)["J"] == 0.0 and exact(660)["J"] == pytest.approx(4e-9 * M) and exact(720)["J"] == exact(END)["J"]
    assert exact(3600)["L"] == 1e9 and exact(END)["L"] == pytest.approx(1e9 * math.exp(-0.1))


def test_the_scenario_holds_what_the_easy_file_says(prepared):
    scenario = json.load(open(prepared["fast"], encoding="utf-8"))
    reactions = {(r["reactants"], r["products"], r["rate"]) for r in scenario["custom_reactions"]}
    assert reactions == {("A", "B", "1.E-4"), ("C", "D", "2.E-4"), ("D", "C", "1.E-4"), ("E + F", "G", "1.E-13"),
                         ("H", "I", "HV*5.E-4"), ("L + DIL", "", "(FLOW/(276.*3600.))"),
                         ("N", "O", "1.E-2*EXP(-1000./T)")}
    (injection,) = scenario["injection_reactions"]
    assert (injection["reactants"], injection["products"]) == ("QJ", "J")
    (pulse,) = scenario["burst_injections"]                   # read back from the rate and the ENZ column
    assert (pulse["target"], float(pulse["ppb"]), float(pulse["duration_s"])) == ("J", 8.0, 120.0)
    assert pulse["start_utc"].endswith("00:10:00")
    assert set(scenario["initials"]) == set("ACEFHLN")
    assert scenario["idl_parameters"]["M"] == "((P)/((T)*1.379E-19))" and scenario["idl_parameters"]["DIL"] == "1."
    assert (scenario["outstep"], scenario["hmax"]) == (600.0, 30.0)


@pytest.mark.slow
# fast is 0.5 % off on N, a species that has decayed fifty-fold by then (relative tolerance 1e-4 per step);
# precise and twostep are within 0.03 % and 0.3 % on everything.
@pytest.mark.parametrize("profile, tol", [("precise", 4e-4), ("fast", 6e-3), ("twostep", 3e-3)])
def test_every_piece_follows_its_closed_form(prepared, profile, tol):
    df = run(prepared[profile])
    times = np.asarray(df.index, dtype=float)
    assert times[0] == 0.0 and times[-1] == END
    worst = []
    for t in times:
        if t % 600:
            continue
        want = exact(t)
        for name, value in want.items():
            got = float(df.loc[t, name])
            scale = max(abs(value), 1e7)                 # 1 % of the starting amounts: below that the floor decides
            worst.append((abs(got - value) / scale, name, t, got, value))
    worst.sort(reverse=True)
    assert worst[0][0] < tol, worst[:5]


@pytest.mark.slow
def test_settings_in_the_result_are_those_of_the_enz_rows(prepared):
    df = run(prepared["fast"])
    # a result line shows the settings of the row that ends at its time: the light goes on at 1800 s and off at 5400 s
    hv = [float(df.loc[t, "HV"]) for t in (0.0, 1800.0, 2400.0, 5400.0, 6000.0)]
    assert hv == [0.0, 0.0, 1.0, 1.0, 0.0]
    assert [float(df.loc[t, "T"]) for t in (3600.0, 4200.0)] == [250.0, 500.0]
    assert [float(df.loc[t, "FLOW"]) for t in (3600.0, 4200.0)] == [0.0, 27.6]
    assert float(df.loc[0.0, "M"]) == pytest.approx(M, rel=1e-6) and float(df.loc[END, "M"]) == pytest.approx(M, rel=1e-6)
