"""injection_calc.py, function by function: injected amount → ppb for each delivery mode, compound names → MCM
species, and mixtures.  The test chamber holds exactly 100 mol of air at 300 K and 1000 hPa, so every expected
mixing ratio can be worked out by hand: 1e-4 mol of a compound is 1000 ppb."""
from types import SimpleNamespace as NS

import pytest

from juice.core import injection_calc as ic
from juice.core.injection_calc import Mixture, MixtureComponent

T, P = 300.0, 1000.0
V = ic.R_GAS * T / 1000.0                       # m³ that hold 100 mol at T, P  (2.494 m³)
CH = NS(volume_m3=V, syringe_T_K=300.0, mfc_std_T_K=300.0, mfc_std_P_hPa=1000.0)


def test_the_test_chamber_holds_100_mol():
    assert ic.moles_air(V, T, P) == pytest.approx(100.0)
    assert ic.moles_air(2 * V, T, P) == pytest.approx(200.0)          # twice the volume
    assert ic.moles_air(V, 2 * T, P) == pytest.approx(50.0)           # twice the temperature
    assert ic.moles_air(V, T, P / 4) == pytest.approx(25.0)           # a quarter of the pressure


# ── delivery modes ───────────────────────────────────────────────────────────

def test_liquid_syringe():
    # 10 µl × 0.5 g/ml = 5 mg; ÷ 50 g/mol = 1e-4 mol; in 100 mol of air = 1000 ppb
    ppb, text = ic.ppb_liquid(10.0, 50.0, 0.5, T, P, CH)
    assert ppb == pytest.approx(1000.0)
    assert "10 µl" in text and "1e+03 ppb" in text


def test_gas_syringe_is_a_volume_ratio_when_filled_at_chamber_conditions():
    # 1 ml of pure gas at the chamber's T and P is 1e-6 m³ of the chamber's V
    ppb, _ = ic.ppb_gas_syringe(1.0, 1.0, T, P, CH)
    assert ppb == pytest.approx(1e-6 / V * 1e9)
    diluted, _ = ic.ppb_gas_syringe(1.0, 0.01, T, P, CH)
    assert diluted == pytest.approx(ppb / 100.0)


def test_gas_syringe_filled_colder_holds_more():
    cold = NS(volume_m3=V, syringe_T_K=150.0, mfc_std_T_K=300.0, mfc_std_P_hPa=1000.0)
    assert ic.ppb_gas_syringe(1.0, 1.0, T, P, cold)[0] == pytest.approx(2 * ic.ppb_gas_syringe(1.0, 1.0, T, P, CH)[0])


def test_mass_flow_controller():
    # 60 sccm for 60 s = 60 standard cm³; the standard is the chamber's T and P here, so again a volume ratio
    ppb, text = ic.ppb_mfc(60.0, 60.0, 1.0, T, P, CH)
    assert ppb == pytest.approx(60e-6 / V * 1e9)
    assert "60 std-cm³" in text


def test_mass_flow_controller_uses_its_own_standard_conditions():
    std = NS(volume_m3=V, syringe_T_K=300.0, mfc_std_T_K=273.15, mfc_std_P_hPa=1013.25)
    ratio = ic.ppb_mfc(60.0, 60.0, 1.0, T, P, std)[0] / ic.ppb_mfc(60.0, 60.0, 1.0, T, P, CH)[0]
    assert ratio == pytest.approx((1013.25 / 1000.0) * (300.0 / 273.15))


@pytest.mark.parametrize("value, unit, fraction", [
    (100.0, "ppm", 1e-4), (5.0, "ppb", 5e-9), (1.0, "%", 0.01), (2.0, "percent", 0.02), (0.3, "pure", 1.0),
    (0.25, "fraction", 0.25), (0.25, "", 0.25), (0.25, None, 0.25), (100.0, "PPM", 1e-4),
])
def test_bottle_fraction(value, unit, fraction):
    assert ic.bottle_fraction(value, unit) == pytest.approx(fraction)


def test_delivery_dict_reaches_the_right_calculation():
    assert ic.resolve_delivery({"mode": "ppb", "value": "7.5"}, T, P, CH) == (7.5, "7.5 ppb as stated")
    assert ic.resolve_delivery({"value": "3"}, T, P, CH)[0] == 3.0                         # no mode: ppb
    liquid = {"mode": "liquid", "vol_ul": "10", "mw": "50", "density": "0.5"}
    assert ic.resolve_delivery(liquid, T, P, CH)[0] == pytest.approx(1000.0)
    syringe = {"mode": "gas_syringe", "vol_ml": "1", "bottle_value": "1", "bottle_unit": "%"}
    assert ic.resolve_delivery(syringe, T, P, CH)[0] == pytest.approx(1e-6 / V * 1e9 / 100.0)
    assert ic.resolve_delivery({"mode": "gas_syringe", "vol_ml": "1"}, T, P, CH)[0] == pytest.approx(1e-6 / V * 1e9)


def test_delivery_flow_in_sccm_or_slpm():
    sccm = {"mode": "mfc", "flow_sccm": "60", "duration_s": "60", "bottle_value": "100", "bottle_unit": "ppm"}
    slpm = {"mode": "mfc", "flow_slpm": "0.06", "duration_s": "60", "bottle_value": "100", "bottle_unit": "ppm"}
    expected = 60e-6 / V * 1e9 * 1e-4
    assert ic.resolve_delivery(sccm, T, P, CH)[0] == pytest.approx(expected)
    assert ic.resolve_delivery(slpm, T, P, CH)[0] == pytest.approx(expected)


def test_delivery_with_a_missing_or_unreadable_field_gives_none_and_says_why():
    ppb, why = ic.resolve_delivery({"mode": "liquid", "vol_ul": "10", "mw": "50"}, T, P, CH)
    assert ppb is None and "density" in why
    ppb, why = ic.resolve_delivery({"mode": "ppb", "value": "seven"}, T, P, CH)
    assert ppb is None and "invalid" in why
    ppb, why = ic.resolve_delivery({"mode": "bucket"}, T, P, CH)
    assert ppb is None and "unknown delivery mode 'bucket'" in why


# ── compounds ────────────────────────────────────────────────────────────────

def test_same_formula_ignores_order_and_explicit_ones():
    assert ic._same_formula("C5H6O", "C5H6O1") and ic._same_formula("C5H6O", "OC5H6")
    assert not ic._same_formula("C5H6O", "C5H6O2") and not ic._same_formula("C5H6O", "C5H8O")
    assert ic._same_formula("C10H30O5Si5", "Si5O5C10H30")            # two-letter element


def test_common_name_inside_a_sentence():
    rc = ic.resolve_compound("20 µl isoprene by syringe")
    assert (rc.mcm, rc.status, rc.mw, rc.density) == ("C5H8", "mcm", 68.12, 0.681)


def test_longer_name_wins_over_the_name_it_contains():
    assert ic.resolve_compound("octadecane").status == "no_mcm"       # not 'decane'
    assert ic.resolve_compound("n-decane").mcm == "NC10H22"
    assert ic.resolve_compound("beta-pinene").mcm == "BPINENE"


def test_compound_without_mcm_species_carries_its_surrogate_hint_and_formula():
    rc = ic.resolve_compound("Texanol")
    assert (rc.mcm, rc.status, rc.candidates, rc.formula) == ("", "no_mcm", ["NC12H26"], "C12H24O3")
    assert rc.mw == 216.32 and "NC12H26" in rc.note
    d5 = ic.resolve_compound("D5 siloxane")
    assert d5.status == "no_mcm" and d5.candidates == [] and "no MCM analogue" in d5.note


def test_typed_mcm_name_and_nothing_at_all():
    assert ic.resolve_compound("nc12h26").mcm == "NC12H26"
    empty = ic.resolve_compound("   ")
    assert (empty.status, empty.note) == ("unresolved", "empty")
    foo = ic.resolve_compound("foo")
    assert foo.status == "unresolved" and "not found" in foo.note
    assert ic.resolve_compound(None).status == "unresolved"            # never raises


def test_formula_with_several_mcm_species_is_ambiguous():
    rc = ic.resolve_compound("C10H16")
    assert rc.status == "ambiguous" and {"APINENE", "BPINENE", "LIMONENE"} <= set(rc.candidates)


def test_addon_species_take_the_compound_only_by_its_own_name():
    addon = {"FURFURAL": "C5H4O2", "MEFURFURAL": "C6H6O2", "M2F": "C5H6O"}
    assert ic.resolve_compound("5-methyl-furfural", addon).mcm == "MEFURFURAL"
    assert ic.resolve_compound("2-methylfuran", addon).mcm == "M2F"                   # by the add-on's short name
    only_furfural = ic.resolve_compound("5-methyl-furfural", {"FURFURAL": "C5H4O2"})
    assert only_furfural.mcm == "" and only_furfural.status == "no_mcm"               # never the name it contains
    assert ic.resolve_compound("m2f", addon).mcm == "M2F"                             # typed species name


def test_addon_species_of_another_formula_is_not_the_compound():
    assert ic.resolve_compound("2-methylfuran", {"M2F": "C6H8O"}).mcm == ""


# ── mixtures ─────────────────────────────────────────────────────────────────

def test_safe_name():
    assert ic.safe_name(" VCP mix #2 ") == "VCP_mix_2" and ic.safe_name("***") == "mixture"


def test_mixture_file_round_trip(tmp_path):
    m = Mixture("Test mix", "two things", ["TM"], [MixtureComponent("limonene", 5.0, "ul", "LIMONENE", "", 136.23, 0.841),
                                                   MixtureComponent("texanol", 2.0, "mg", "", "NC12H26", 216.32, 0.95)],
                "bench notes", 20.0)
    path = ic.save_mixture(m, str(tmp_path / "t.mixture.json"))
    assert ic.load_mixture(path) == m


def test_mixture_folder_skips_unreadable_and_other_files(tmp_path):
    ic.save_mixture(Mixture("A"), str(tmp_path / "a.mixture.json"))
    (tmp_path / "broken.mixture.json").write_text("{not json", encoding="utf-8")
    (tmp_path / "notes.txt").write_text("x", encoding="utf-8")
    assert list(ic.load_mixtures(str(tmp_path))) == ["A"]
    assert ic.load_mixtures(str(tmp_path / "missing")) == {}


def test_mixture_is_found_by_name_or_alias_as_a_whole_word():
    mixtures = {"Cooking": Mixture("Cooking"), "VCP2": Mixture("VCP2", aliases=["VCP"])}
    assert ic.find_mixture("10:15 injected 20 µl VCP solution", mixtures).name == "VCP2"
    assert ic.find_mixture("cooking emissions", mixtures).name == "Cooking"
    assert ic.find_mixture("VCPX bottle", mixtures) is None
    assert ic.find_mixture("", mixtures) is None


def test_recipe_lines():
    comps = ic.parse_mixture_text("* Limonene 5.61 μl\nethanol, 59,19, ul\ntoluene\t6.08\nwater 2 %\n\nno number here\n")
    assert [(c.name, c.amount, c.unit, c.mcm) for c in comps] == [
        ("Limonene", 5.61, "ul", "LIMONENE"), ("ethanol", 59.19, "ul", "C2H5OH"), ("toluene", 6.08, "ul", "TOLUENE"),
        ("water", 2.0, "vol%", "")]
    assert comps[0].mw == 136.23 and comps[0].density == 0.841
    assert ic.parse_mixture_text("texanol 3", default_unit="mg")[0].unit == "mg"


def two_liquids():
    # 10 µl and 30 µl of density 1 g/ml and 100 g/mol: 1e-4 and 3e-4 mol
    return Mixture("M", components=[MixtureComponent("limonene", 10.0, "ul", "LIMONENE", "", 100.0, 1.0),
                                    MixtureComponent("ethanol", 30.0, "ul", "C2H5OH", "", 100.0, 1.0)])


def test_mixture_moles_as_written_and_scaled_to_the_injected_volume():
    moles, v_recipe, text = ic.mixture_moles(two_liquids(), None)
    assert moles == pytest.approx({"limonene": 1e-4, "ethanol": 3e-4}) and v_recipe == 40.0 and "used as is" in text
    half, _, text = ic.mixture_moles(two_liquids(), 20.0)
    assert half == pytest.approx({"limonene": 0.5e-4, "ethanol": 1.5e-4}) and "scale 0.5" in text


@pytest.mark.parametrize("amount, unit, mol", [
    (2.0, "ml", 2e-2),          # 2 ml × 1 g/ml ÷ 100 g/mol
    (50.0, "mg", 5e-4),
    (0.05, "g", 5e-4),
    (50.0, "wt%", 5e-4),        # relative recipe: read as mg
    (40.0, "vol%", 4e-4),       # relative recipe: read as µl
    (3.0, "mol%", 3e-3),        # relative moles: 3 × 100 g/mol = 300 mg
])
def test_mixture_units(amount, unit, mol):
    m = Mixture("M", components=[MixtureComponent("x", amount, unit, "X", "", 100.0, 1.0)])
    assert ic.mixture_moles(m, None)[0] == pytest.approx({"x": mol})


def test_mole_fractions_of_a_mol_percent_recipe_survive_the_scaling():
    m = Mixture("M", components=[MixtureComponent("a", 25.0, "mol%", "A", "", 50.0, 0.5),
                                 MixtureComponent("b", 75.0, "mol%", "B", "", 200.0, 2.0)])
    moles, _, _ = ic.mixture_moles(m, 10.0)
    assert moles["b"] / moles["a"] == pytest.approx(3.0)


def test_components_without_molar_mass_or_given_as_ppb_per_ul_are_left_out_of_the_moles():
    m = Mixture("M", components=[MixtureComponent("known", 10.0, "ul", "X", "", 100.0, 1.0),
                                 MixtureComponent("no mass", 10.0, "ul", "Y", "", None, 1.0),
                                 MixtureComponent("direct", 0.2, "ppb/ul", "Z"),
                                 MixtureComponent("odd unit", 1.0, "drops", "W", "", 100.0, 1.0)])
    moles, v_recipe, _ = ic.mixture_moles(m, None)
    assert list(moles) == ["known"] and v_recipe == 10.0


def test_mixture_event_carries_mole_fractions_and_the_modelled_total():
    ev, species, alerts = ic.mixture_event(two_liquids(), None, "10:00", 2.0, T, P, CH)
    assert ev.products == "0.75*C2H5OH+0.25*LIMONENE"        # largest first
    assert float(ev.total_ppb) == pytest.approx(4000.0)      # 1000 + 3000 ppb in 100 mol of air
    assert (ev.name, ev.column, ev.start, ev.duration_min, ev.action) == ("M", "Q_M", "10:00", 2.0, "mixture")
    assert set(species) == {"LIMONENE", "C2H5OH"} and alerts == []


def test_mixture_event_scales_with_the_injected_volume_and_falls_back_to_the_reference_volume():
    ev, _, _ = ic.mixture_event(two_liquids(), 20.0, "10:00", 1.0, T, P, CH)
    assert float(ev.total_ppb) == pytest.approx(2000.0) and ev.delivery["injected_ul"] == "20"
    m = two_liquids()
    m.reference_ul = 80.0
    assert float(ic.mixture_event(m, None, "10:00", 1.0, T, P, CH)[0].total_ppb) == pytest.approx(8000.0)


def test_two_components_on_one_species_are_summed_and_ppb_per_ul_is_taken_directly():
    m = two_liquids()
    m.components.append(MixtureComponent("more ethanol", 10.0, "ul", "C2H5OH", "", 100.0, 1.0))
    m.components.append(MixtureComponent("toluene", 2.0, "ppb/ul", "TOLUENE", "", 92.14))
    ev, species, _ = ic.mixture_event(m, 50.0, "10:00", 1.0, T, P, CH)       # recipe is 50 µl: scale 1
    # ethanol 3000 + 1000, limonene 1000, toluene 2 ppb/µl × 50 µl = 100 → 5100 ppb
    assert float(ev.total_ppb) == pytest.approx(5100.0)
    assert ev.products == "0.7843*C2H5OH+0.1961*LIMONENE+0.01961*TOLUENE"
    assert species == ["LIMONENE", "C2H5OH", "TOLUENE"]


def test_mixture_with_nothing_to_model_is_an_error_and_has_no_products():
    m = Mixture("Empty", components=[MixtureComponent("mystery", 10.0, "ul", "", "NC12H26", None, None)])
    ev, species, alerts = ic.mixture_event(m, None, "10:00", 1.0, T, P, CH)
    assert ev.products == "" and species == []
    levels = [level for level, _ in alerts]
    assert "error" in levels and any("no component can be modelled" in text for _, text in alerts)


@pytest.mark.parametrize("line", ["ethanol, 59.19, ul", "ethanol;59.19;ul", "ethanol\t59.19\tul", "ethanol 59,19 ul",
                                  "- ethanol: 59.19ul", "ethanol, 59.19 ul,"])
def test_recipe_line_separators(line):
    """Bug of 24 Nov 2026: a separator between number and unit ('ethanol, 59.19, ul', the form the docstring
    gives) dropped the line without a word."""
    (c,) = ic.parse_mixture_text(line)
    assert (c.name, c.amount, c.unit) == ("ethanol", 59.19, "ul")
