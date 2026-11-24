"""lamp_photolysis.py: reading a lamp spectrum, the built-in stand-ins, the clear-sky ratios, the convolution with
the absorption data, and which of the profile's entries wins."""
import math
from types import SimpleNamespace as NS

import numpy as np
import pytest

from juice.chemistry import lamp_photolysis as lp

needs_table = pytest.mark.skipif(not lp.sigma_phi_available(), reason="chemistry/photolysis/sigma_phi_298K.csv is not there")


def spectrum_file(tmp_path, text, name="lamp.csv"):
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return str(path)


ROWS = "".join(f"{lam},{value}\n" for lam, value in ((360, 4.0), (340, 2.0), (350, 3.0), (370, -1.0), (380, 1.0)))


# ── a spectrum file ──────────────────────────────────────────────────────────

def test_energy_spectrum_is_turned_into_photons_sorted_and_cleared_of_negative_values(tmp_path):
    lam, F, note = lp.spectrum_from_file(spectrum_file(tmp_path, "wavelength_nm,intensity\n" + ROWS))
    assert list(lam) == [340, 350, 360, 370, 380]
    assert list(F) == [2.0 * 340, 3.0 * 350, 4.0 * 360, 0.0, 1.0 * 380]          # energy × wavelength
    assert "340-380 nm, 5 points, read as energy units" in note and "energy assumed" in note


def test_photon_spectrum_is_taken_as_it_is(tmp_path):
    lam, F, note = lp.spectrum_from_file(spectrum_file(tmp_path, "# units: photons\n; a comment\n\n" + ROWS))
    assert list(F) == [2.0, 3.0, 4.0, 0.0, 1.0] and "read as photons units" in note and "assumed" not in note


def test_units_given_by_the_caller_win_over_the_header(tmp_path):
    path = spectrum_file(tmp_path, "# units: photons\n" + ROWS)
    assert list(lp.spectrum_from_file(path, units="energy")[1])[:2] == [680.0, 1050.0]
    path = spectrum_file(tmp_path, "# Units = W m-2 nm-1\n" + ROWS, "w.csv")
    assert list(lp.spectrum_from_file(path)[1])[:2] == [680.0, 1050.0]


def test_blank_or_tab_separated_columns(tmp_path):
    text = "# units: photons\n" + ROWS.replace(",", "\t")
    assert list(lp.spectrum_from_file(spectrum_file(tmp_path, text))[1]) == [2.0, 3.0, 4.0, 0.0, 1.0]


def test_a_file_with_fewer_than_five_rows_is_refused(tmp_path):
    with pytest.raises(ValueError, match="fewer than 5"):
        lp.spectrum_from_file(spectrum_file(tmp_path, "340,1\n350,2\nnm,counts\n"))


# ── built-in stand-ins and clear sky ─────────────────────────────────────────

@pytest.mark.parametrize("name", list(lp.BUILTIN_SPECTRA))
def test_stand_in_is_a_gaussian_of_the_stated_centre_and_width(name):
    centre, fwhm, _ = lp.BUILTIN_SPECTRA[name]
    lam, F, note = lp.builtin_spectrum(name)
    assert lam[int(F.argmax())] == pytest.approx(centre, abs=0.5) and F.max() == pytest.approx(1.0, abs=1e-3)
    assert np.interp(centre + fwhm / 2, lam, F) == pytest.approx(0.5, abs=0.01)
    assert np.interp(centre - fwhm / 2, lam, F) == pytest.approx(0.5, abs=0.01)
    assert name in note


def test_narrower_short_side_changes_only_the_short_side():
    name = next(iter(lp.BUILTIN_SPECTRA))
    centre, fwhm, _ = lp.BUILTIN_SPECTRA[name]
    lam, wide, _ = lp.builtin_spectrum(name)
    _, narrow, _ = lp.builtin_spectrum(name, sigma_short=0.5 * fwhm / 2.3548)
    assert np.array_equal(wide[lam >= centre], narrow[lam >= centre])
    assert np.all(narrow[lam < centre - 1] < wide[lam < centre - 1])


def test_clear_sky_ratios_are_relative_to_no2_and_follow_the_mcm_expression():
    ratios = lp.clear_sky_ratios(30.0)
    assert ratios[lp.J_NO2] == 1.0
    c = math.cos(math.radians(30.0))
    l, m, n = lp.MCM_LMN[1]
    l4, m4, n4 = lp.MCM_LMN[lp.J_NO2]
    assert ratios[1] == pytest.approx(l * c ** m * math.exp(-n / c) / (l4 * c ** m4 * math.exp(-n4 / c)))
    assert lp.clear_sky_ratios(60.0)[1] < ratios[1]                    # j(O1D) falls faster than jNO2 with the sun
    assert lp.clear_sky_ratios(95.0) == lp.clear_sky_ratios(120.0)     # below the horizon: the same floor


def test_missing_j_numbers_are_scaled_from_a_covered_one_by_the_clear_sky_ratio():
    cs = lp.clear_sky_ratios(30.0)
    ratios = {lp.J_NO2: 1.0, 1: 5e-4}
    filled = lp._fill_missing(ratios)
    assert set(ratios) == set(lp.MCM_LMN) and set(filled) == set(lp.MCM_LMN) - {lp.J_NO2, 1}
    for j, ref in filled.items():
        assert ref in (lp.J_NO2, 1)
        assert ratios[j] == pytest.approx({lp.J_NO2: 1.0, 1: 5e-4}[ref] * cs[j] / cs[ref])


# ── the convolution ──────────────────────────────────────────────────────────

@needs_table
def test_ratios_are_normalised_to_no2_and_a_bluer_lamp_photolyses_more_ozone():
    lam = np.arange(250.0, 500.0, 1.0)
    at_365, notes = lp.convolve_ratios(lam, np.exp(-0.5 * ((lam - 365.0) / 10.0) ** 2))
    at_310, _ = lp.convolve_ratios(lam, np.exp(-0.5 * ((lam - 310.0) / 10.0) ** 2))
    assert at_365[lp.J_NO2] == 1.0 and at_310[lp.J_NO2] == 1.0 and notes == []
    assert at_310[1] > 100 * at_365[1]                                  # O3 -> O1D needs light below about 320 nm
    assert at_365[1] < 1e-4


@needs_table
def test_scaling_the_spectrum_changes_nothing_and_a_narrow_one_is_noted():
    lam = np.arange(320.0, 381.0, 1.0)
    F = np.exp(-0.5 * ((lam - 350.0) / 10.0) ** 2)
    one, notes = lp.convolve_ratios(lam, F)
    thousand, _ = lp.convolve_ratios(lam, 1000.0 * F)
    assert thousand == pytest.approx(one)
    assert notes == ["spectrum covers only 320-380 nm; photolysis outside is taken as zero"]


@needs_table
def test_a_spectrum_that_no2_does_not_absorb_cannot_be_used():
    lam = np.arange(600.0, 700.0, 1.0)
    with pytest.raises(ValueError, match="no overlap with the NO2 absorption"):
        lp.convolve_ratios(lam, np.ones(len(lam)))


# ── what the profile says ────────────────────────────────────────────────────

def test_profile_without_anything_falls_back_to_clear_sky_and_says_so():
    lj = lp.lamp_j(NS(light={}))
    assert lj.method == "clear-sky fallback" and lj.ratios == pytest.approx(lp.clear_sky_ratios(30.0))
    assert "NO lamp spectrum" in lj.source and any("give the profile a lamp spectrum" in n for n in lj.notes)
    assert lp.lamp_j(NS()).method == "clear-sky fallback"              # a profile without a light entry at all


def test_solar_like_source_uses_its_effective_zenith_angle():
    lj = lp.lamp_j(NS(light={"solar_like": True, "effective_sza": 45}))
    assert lj.method == "solar-like" and "45 deg" in lj.source
    assert lj.ratios[1] == pytest.approx(lp.clear_sky_ratios(45.0)[1])


@needs_table
def test_built_in_spectrum_is_flagged_as_synthetic_with_its_sensitivity():
    name = next(iter(lp.BUILTIN_SPECTRA))
    lj = lp.lamp_j(NS(light={"spectrum_builtin": name}))
    assert lj.method == "builtin spectrum" and "SYNTHETIC" in lj.source and lj.ratios[lp.J_NO2] == 1.0
    assert any("stand-in sensitivity" in n for n in lj.notes)
    assert set(lj.ratios) >= set(lp.MCM_LMN)                            # every MCM J has a value, filled where needed


@needs_table
def test_spectrum_file_of_the_profile_is_used_and_wins_over_the_stand_in(tmp_path):
    lam = np.arange(300.0, 420.0, 2.0)
    text = "# units: photons\n" + "".join(f"{x},{math.exp(-0.5 * ((x - 365.0) / 8.0) ** 2)}\n" for x in lam)
    spectrum_file(tmp_path, text, "my_lamp.csv")
    lj = lp.lamp_j(NS(light={"spectrum": "my_lamp.csv", "spectrum_builtin": next(iter(lp.BUILTIN_SPECTRA))}), str(tmp_path))
    assert lj.method == "spectrum" and "my_lamp.csv" in lj.source and "read as photons units" in lj.notes[0]


def test_spectrum_file_that_is_missing_is_said_and_clear_sky_is_used():
    lj = lp.lamp_j(NS(light={"spectrum": "no_such_lamp.csv"}), "no_such_folder")
    assert lj.method == "clear-sky fallback" and any("lamp spectrum not usable" in n for n in lj.notes)


def test_measured_ratios_override_and_no2_stays_one():
    lj = lp.lamp_j(NS(light={"j_ratios": {"J1": 1e-5, "3": 2e-4, "J4": 7.0, "J11": None, "jx": 1.0}}))
    assert lj.ratios[1] == 1e-5 and lj.ratios[3] == 2e-4 and lj.ratios[lp.J_NO2] == 1.0
    assert lj.ratios[11] == pytest.approx(lp.clear_sky_ratios(30.0)[11])      # no value given: not overridden
    assert "measured ratios for J1, J3, J4" in lj.source and lj.method == "clear-sky fallback"
    complete = lp.lamp_j(NS(light={"j_ratios": {f"J{j}": 1e-4 for j in (1, 3, 7, 11, 12)}}))
    assert complete.method == "measured ratios"                         # the five that matter are all measured
