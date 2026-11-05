"""A roof that never opens (lights never on) is a dark experiment: a note, not a problem."""
import experiment_spec as es


def _spec(chamber: str, roof_open):
    s = es.ExperimentSpec(name="dark_check")
    es.apply_chamber(s, chamber)
    c = s.conditions
    c.source, c.date, c.start, c.end, c.step_s = "manual", "2024-07-09", "08:00", "12:00", 60
    c.temperature_K, c.temperature_mode, c.pressure_hPa, c.h2o_percent, c.flow_m3h = 298.0, "const", 1013.0, 1.0, 0.5
    c.roof_open = roof_open
    s.inorganic_only = True
    return s


def _dark_msgs(msgs):
    return [(lvl, m) for lvl, m in msgs if "never" in m and "dark" in m]


import pytest


@pytest.mark.parametrize("chamber", ["SAPHIR", "EUPHORE", "HELIOS", "generic_lamps"])
def test_closed_roof_or_lights_off_is_a_note(chamber):
    msgs = es.validate(_spec(chamber, []))
    dark = _dark_msgs(msgs)
    assert dark and all(lvl == "info" for lvl, _ in dark)


@pytest.mark.parametrize("chamber", ["SAPHIR", "EUPHORE"])
def test_dark_run_has_no_photolysis_warnings(chamber):
    s = _spec(chamber, [])
    s.conditions.j_source = "measured"          # irrelevant in the dark: must not be flagged
    assert not [m for lvl, m in es.validate(s) if lvl != "info" and ("photolysis" in m.lower() or "j values" in m.lower())]


def test_lit_run_still_warns_about_missing_j():
    s = _spec("SAPHIR", [("09:00", "11:00")])
    s.conditions.j_source = "measured"
    assert any("j values" in m.lower() for lvl, m in es.validate(s) if lvl == "warning")


def test_open_roof_gives_no_dark_note():
    assert not _dark_msgs(es.validate(_spec("SAPHIR", [("09:00", "11:00")])))


def test_old_wiki_note_is_not_repeated_as_a_warning():
    s = _spec("SAPHIR", [])
    s.import_notes = ["no roof open/close (R/r) entries: dark experiment"]
    assert not [m for lvl, m in es.validate(s) if lvl == "warning" and "roof" in m.lower()]


def test_level_icons():
    assert es.level_icon("info").startswith("ℹ") and es.level_icon("error").startswith("✖")
