"""The small functions of the results viewer that need no window: species classes of the nitrogen and carbon
budgets, atom counts, colours, the time axis, stoichiometry and labels."""
import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import pytest

from juice.core.mechanism import Reaction
from juice.gui import results_viewer as rv


@pytest.mark.parametrize("species, group", [
    ("NO", "NOx"), ("no2", "NOx"), ("NO3", "NOx"),
    ("HNO3", "HNO3"), ("HO2NO2", "HNO3"),
    ("HONO", "HONO"), ("N2O5", "N2O5"),
    ("PAN", "PANs"), ("MPAN", "PANs"), ("C5PAN5", "PANs"), ("PPN", "PANs"), ("CH3CO3NO2", "PANs"),
    ("CH3NO3", "Org. Nitrates"), ("ISOPBNO3", "Org. Nitrates"), ("CH3O2NO2", "Org. Nitrates"), ("C2H5ONO2", "Org. Nitrates"),
    ("CLNO2", "Other N"), ("CLONO2", "Other N"), ("NA", "Other N"), ("NH3", "Other N"),
])
def test_nitrogen_budget_groups(species, group):
    assert rv._nitrogen_class(species) == group


def test_pernitric_acid_and_halogen_nitrates_are_not_organic_nitrates():
    """Bug of 24 Nov 2026: HO2NO2 (the MCM's name of HNO4) and ClNO2 / ClONO2 were drawn as organic nitrates."""
    assert "Org. Nitrates" not in {rv._nitrogen_class(s) for s in ("HO2NO2", "CLNO2", "CLONO2", "BRONO2")}


@pytest.mark.parametrize("species, group", [
    ("C5H8", "Parent VOC"), ("CO", "CO"), ("CO2", "Other"), ("HNO3", "Other"), ("HO2NO2", "Other"), ("O3", "Other"),
    ("PAN", "PANs"), ("CH3CO3NO2", "PANs"), ("CH3NO3", "Org. Nitrates"), ("CH3O2NO2", "Org. Nitrates"),
    ("CH3CO3H", "Hydroperoxides"), ("CH3CO3", "Acyl Peroxy"), ("CH3O2", "Peroxy Radicals"), ("ISOPBO2", "Peroxy Radicals"),
    ("CH3OOH", "Hydroperoxides"), ("CH3CO2H", "Carboxylic Acids"), ("HCOOH", "Carboxylic Acids"),
    ("HCHO", "Aldehydes"), ("CH3CHO", "Aldehydes"), ("PINAL", "Aldehydes"),
    ("CH3COCH3", "Ketones"), ("NOPINONE", "Ketones"), ("C2H5OH", "Alcohols"), ("ETHGLY", "Alcohols"), ("CH4", "Other"),
])
def test_carbon_budget_groups(species, group):
    assert rv._carbon_class(species, parent_voc={"C5H8"}) == group


def test_formic_acid_is_an_acid_and_the_cooh_names_of_the_mcm_are_hydroperoxides():
    """Bug of 24 Nov 2026: HCOOH was drawn with the hydroperoxides, because names ending in OOH were tested first."""
    assert rv._carbon_class("HCOOH") == "Carboxylic Acids"
    assert rv._carbon_class("APINCOOH") == "Hydroperoxides" and rv._carbon_class("BZEMUCOOH") == "Hydroperoxides"


def test_an_injected_compound_is_a_parent_only_when_it_is_named_as_one():
    assert rv._carbon_class("C2H5OH", parent_voc={"C2H5OH"}) == "Parent VOC"
    assert rv._carbon_class("C2H5OH") == "Alcohols"
    assert rv._carbon_class("HNO3", parent_voc={"HNO3"}) == "Other"            # an inorganic species never is


def test_atom_counts_from_the_species_table_or_from_the_name(monkeypatch):
    if rv.MCM_DB.get("APINENE") is not None:
        assert rv._carbon_count("APINENE") == 10 and rv._nitrogen_count("PAN") == 1 and rv._carbon_count("O3") == 0
    monkeypatch.setattr(rv, "MCM_DB", {})                                       # a species the table does not have
    assert rv._carbon_count("C5H8") == 5 and rv._carbon_count("CH3CHO") == 2 and rv._carbon_count("Cl2") == 0
    assert rv._nitrogen_count("N2O5") == 2 and rv._nitrogen_count("CH3CO3NO2") == 1 and rv._nitrogen_count("HONO") == 1
    assert rv._nitrogen_count("XONO2NO3") == 2 and rv._nitrogen_count("NO") == 1 and rv._nitrogen_count("CH3NO") == 1
    # bug of 24 Nov 2026: the letters NO inside a name were counted as a nitrogen atom
    assert rv._nitrogen_count("ETHANOL") == 0 and rv._nitrogen_count("NOPINONE") == 0


@pytest.mark.parametrize("hours, text", [(0.0, "00:00"), (9.5, "09:30"), (13.999, "14:00"), (24.25, "00:15"), (-0.5, "23:30")])
def test_time_axis_labels(hours, text):
    assert rv._fmt_hhmm(hours, None) == text


def test_colours_are_distinct_and_other_is_grey():
    labels = [f"reaction {i}" for i in range(30)] + ["Other (12)"]
    colours = rv._distinct_colors(labels)
    assert colours[-1] == "#9a9a9a" and len({tuple(c) for c in colours[:-1]}) == 30
    assert rv._distinct_colors(["other", "a"])[1] == rv._distinct_colors(["a"])[0]     # 'other' does not use up a colour


def test_shades_of_one_colour_run_from_light_to_dark():
    base = (0.2, 0.4, 0.8)
    assert rv._shade(base, 0, 1) == pytest.approx(base)                         # a single member keeps the colour
    shades = [rv._shade(base, k, 5) for k in range(5)]
    brightness = [sum(s) for s in shades]
    assert brightness == sorted(brightness, reverse=True) and shades[2] == pytest.approx(base)
    assert all(0.0 <= v <= 1.0 for s in shades for v in s)


def test_juice_result_is_told_from_a_facsimile_one(tmp_path):
    juice = tmp_path / "a.cs"
    juice.write_text("# Calculated by JUICE v4.0\nCOLUMN 1=TIME\n", encoding="latin-1")
    facsimile = tmp_path / "b.cs"
    facsimile.write_text(" * Facsimile-Output / EASY Ver. 6.00\n", encoding="latin-1")
    assert rv._is_juice_output(str(juice)) and not rv._is_juice_output(str(facsimile))
    assert not rv._is_juice_output(str(tmp_path / "missing.cs"))


def test_stoichiometry_and_reaction_labels():
    assert rv.BudgetEngine._parse_stoich(["oh", " 2*HO2 ", "0.5*ch3o2", "C2H5OH*X"]) == \
        [(1.0, "OH"), (2.0, "HO2"), (0.5, "CH3O2"), (1.0, "C2H5OH")]
    assert rv.BudgetEngine._make_label(Reaction(["OH", "2*HO2"], ["H2O", "O2"], "k")) == "OH + HO2 → H2O + O2"
    assert rv.BudgetEngine._make_label(Reaction([], ["HONO"], "k")) == "hv/wall → HONO"
    long = rv.BudgetEngine._make_label(Reaction(["A", "B", "C"], ["D", "E", "F", "G"], "k"))
    assert long == "A + B + ... → D + E + F + ..."
