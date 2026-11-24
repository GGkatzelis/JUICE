"""The budget engine behind Process Analysis and Reactivity (results_viewer.BudgetEngine), on the run of the
exact-answer experiment (tests/exact_small.py): every production and loss term there has a value that can be
written down.  Rates are in ppb per hour, reactivities in 1/s."""
import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import numpy as np
import pytest

from exact_small import M
from juice.gui import results_viewer as rv

PPB_H = 1e9 * 3600.0 / M              # molecules cm-3 s-1 -> ppb per hour


@pytest.fixture(scope="module")
def engine_and_run(exact_small_run):
    js, df = exact_small_run
    engine = rv.BudgetEngine()
    engine.load(js)
    return engine, df


def at(series, df, t):
    return float(series[list(df.index).index(float(t))])


def test_the_engine_holds_each_reaction_once(engine_and_run):
    engine, _df = engine_and_run
    labels = sorted(engine._make_label(r) for r in engine.mechanism.reactions)
    assert labels == ["A → B", "C → D", "D → C", "E + F → G", "H → I", "L + DIL → ", "N → O", "QJ → J"]


def test_first_order_loss_and_the_matching_production(engine_and_run):
    engine, df = engine_and_run
    prod, dest = engine.compute("A", df)
    assert prod == {} and list(dest) == ["A → B"]
    assert dest["A → B"] == pytest.approx(1e-4 * df["A"].to_numpy() * PPB_H, rel=1e-6)
    prod_b, dest_b = engine.compute("B", df)
    assert dest_b == {} and prod_b["A → B"] == pytest.approx(dest["A → B"])


def test_both_directions_of_an_equilibrium(engine_and_run):
    engine, df = engine_and_run
    prod, dest = engine.compute("C", df)
    assert at(dest["C → D"], df, 0) == pytest.approx(2e-4 * 3e9 * PPB_H)        # all of C at the start
    assert at(prod["D → C"], df, 0) == 0.0
    assert at(prod["D → C"], df, 7200) == pytest.approx(1e-4 * at(df["D"].to_numpy(), df, 7200) * PPB_H)


def test_second_order_reaction_counts_for_both_reactants(engine_and_run):
    engine, df = engine_and_run
    rate = 1e-13 * df["E"].to_numpy() * df["F"].to_numpy() * PPB_H
    assert engine.compute("E", df)[1]["E + F → G"] == pytest.approx(rate, rel=1e-6)
    assert engine.compute("F", df)[1]["E + F → G"] == pytest.approx(rate, rel=1e-6)
    assert engine.compute("G", df)[0]["E + F → G"] == pytest.approx(rate, rel=1e-6)


def test_photolysis_like_loss_only_in_the_rows_with_light(engine_and_run):
    engine, df = engine_and_run
    loss = engine.compute("H", df)[1]["H → I"]
    assert at(loss, df, 1200) == 0.0 and at(loss, df, 6600) == 0.0
    assert at(loss, df, 3000) == pytest.approx(5e-4 * at(df["H"].to_numpy(), df, 3000) * PPB_H)


def test_injection_and_dilution_terms(engine_and_run):
    engine, df = engine_and_run
    prod, _dest = engine.compute("J", df)
    injection = prod["QJ → J"]
    assert injection.max() == pytest.approx(8.0 / 60.0 * 0.5 * 3600.0)           # 8 ppb over two minutes = 240 ppb/h
    assert at(injection, df, 3000) == 0.0
    loss = engine.compute("L", df)[1]["L + DIL → "]
    assert at(loss, df, 3000) == 0.0                                            # no flow yet
    assert at(loss, df, 6000) == pytest.approx(27.6 / (276.0 * 3600.0) * at(df["L"].to_numpy(), df, 6000) * PPB_H)


def test_rate_constant_follows_the_temperature_of_the_row(engine_and_run):
    engine, df = engine_and_run
    loss = engine.compute("N", df)[1]["N → O"]
    n = df["N"].to_numpy()
    assert at(loss, df, 3000) == pytest.approx(1e-2 * np.exp(-1000.0 / 250.0) * at(n, df, 3000) * PPB_H)
    assert at(loss, df, 6000) == pytest.approx(1e-2 * np.exp(-1000.0 / 500.0) * at(n, df, 6000) * PPB_H)


def test_reactivity_is_the_rate_without_the_species_itself(engine_and_run):
    engine, df = engine_and_run
    assert engine.compute_reactivity("A", df)["A → B"] == pytest.approx(np.full(len(df), 1e-4))
    assert engine.compute_reactivity("E", df)["E + F → G"] == pytest.approx(1e-13 * df["F"].to_numpy())
    assert engine.compute_reactivity("B", df) == {}                              # B reacts with nothing


def test_species_that_is_not_in_the_run_and_an_engine_that_is_not_loaded(engine_and_run):
    engine, df = engine_and_run
    assert engine.compute("NOSUCH", df) == ({}, {})
    empty = rv.BudgetEngine()
    assert empty.compute("A", df) == ({}, {}) and empty.compute_reactivity("A", df) == {}


def test_pathway_source_gives_the_producers_of_a_species(engine_and_run):
    engine, df = engine_and_run
    producers, hours = engine.pathway_source(df)
    assert len(hours) == len(df)
    (from_e_and_f,) = producers("g")
    assert from_e_and_f.label == "E + F → G" and from_e_and_f.rate_ppb_h == pytest.approx(engine.compute("G", df)[0]["E + F → G"])
    assert (from_e_and_f.reactants, from_e_and_f.yield_, from_e_and_f.photolysis) == ([(1.0, "E"), (1.0, "F")], 1.0, False)
    assert producers("A") == []
