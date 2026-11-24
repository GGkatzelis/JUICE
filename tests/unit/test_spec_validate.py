"""experiment_spec.validate: what the review of an experiment says, one condition at a time.  The start is an
inorganic experiment with manual conditions, about which the review has nothing to object."""
import pytest

from juice.core import experiment_spec as es
from juice.core.experiment_spec import InjectionEvent


def experiment(**conditions):
    spec = es.ExperimentSpec()
    spec.inorganic_only = True
    c = spec.conditions
    c.date, c.start, c.end, c.source, c.j_source = "2024-07-09", "08:00", "12:00", "manual", "clear_sky"
    c.roof_open = [("09:00", "11:00")]
    for k, v in conditions.items():
        setattr(c, k, v)
    spec.injections = [InjectionEvent("NO", "NO", "2", "09:30", 1.0, "QNO", "qNO", action="no")]
    return spec


def said(spec, level, known_species=None):
    return [text for lvl, text in es.validate(spec, known_species) if lvl == level]


def has(messages, *parts):
    return any(all(p in m for p in parts) for m in messages)


def test_the_starting_experiment_has_nothing_to_fix():
    assert es.validate(experiment()) == []


@pytest.mark.parametrize("conditions, part", [
    (dict(date=""), "Enter the experiment date"),
    (dict(end=""), "Enter the experiment date"),
    (dict(start="12:00", end="08:00"), "Time window is empty"),
    (dict(start="12:00", end="12:00"), "Time window is empty"),
    (dict(date="9 July"), "Date or time is not valid"),
    (dict(start="eight"), "Date or time is not valid"),
    (dict(source="somewhere"), "Choose where the chamber conditions come from"),
    (dict(source="file"), "choose the file"),
    (dict(source="file", input_file="no_such_file.csv"), "Measured inputs file not found"),
])
def test_conditions_that_must_be_fixed(conditions, part):
    assert has(said(experiment(**conditions), "error"), part)


def test_an_experiment_without_chemistry_is_an_error():
    spec = experiment()
    spec.inorganic_only = False
    assert has(said(spec, "error"), "No chemistry")
    spec.precursors = ["C2H5OH"]
    assert not has(said(spec, "error"), "No chemistry")


def test_measured_photolysis_without_data_is_a_warning_only_when_the_roof_opens():
    spec = experiment(j_source="measured")
    assert has(said(spec, "warning"), "Measured J values selected without archive data")
    spec.conditions.roof_open = []
    assert not has(said(spec, "warning"), "Measured J values")
    assert said(spec, "info") == ["Roof never opens: run as a dark experiment (no photolysis)."]


def test_injection_times_and_amounts():
    spec = experiment()
    spec.injections = [InjectionEvent("NO", "NO", "2", "half past nine", 1.0, "QNO", "qNO", action="no"),
                       InjectionEvent("ethanol", "C2H5OH", "ppb_today", "10:00", 2.0, action="voc"),
                       InjectionEvent("fan", "", "0", "", 1.0),
                       InjectionEvent("humidification", "", "0", "08:30", 30.0)]
    assert has(said(spec, "error"), "Injection 'NO'", "is not HH:MM")
    warnings = said(spec, "warning")
    assert has(warnings, "Injection 'ethanol'", "is an expression", "written verbatim")
    assert has(warnings, "Injection 'fan'", "neither modelled nor timed")
    assert not has(warnings, "humidification")                          # timed, not modelled: an ENZ column only
    spec.injections = []
    assert said(spec, "warning") == ["No injections defined."]


def test_injected_species_must_be_in_the_mechanism():
    spec = experiment()
    spec.injections = [InjectionEvent("mix", "0.6*LIMONENE+0.4*XYZ123", "10", "10:00", 2.0, action="mixture")]
    warnings = said(spec, "warning", known_species={"LIMONENE", "NO", "NO2"})
    assert has(warnings, "Injection 'mix'", "product 'XYZ123' is not a species of the mechanism")
    assert not has(warnings, "LIMONENE")
    assert not has(said(spec, "warning"), "is not a species")           # without the species list: not checked


def test_partitioning_needs_compounds_and_a_known_volatility_method():
    spec = experiment()
    spec.partitioning = True
    spec.inorganic_only, spec.mechanism_stem = False, "ethanol.fac"    # an existing include set: no compounds chosen
    assert has(said(spec, "warning"), "partitioning is switched on", "this run has none")
    spec.mechanism_stem, spec.precursors, spec.volatility_method = "", ["APINENE"], "guesswork"
    assert has(said(spec, "error"), "unknown volatility method 'guesswork'")


def test_particle_file_must_exist():
    spec = experiment()
    spec.inorganic_only, spec.precursors, spec.partitioning = False, ["APINENE"], True
    spec.particles, spec.particle_file = "file", "no_such_distribution.csv"
    assert has(said(spec, "error"), "size-distribution file is not found")
