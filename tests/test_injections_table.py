"""Injected compounds against the mechanism: add-on hints, the compounds table of the Injections step, and the
compound search of complete bases (GECKO-A starting compounds)."""
import os

import pytest

import juice
from juice.chemistry import mechanism_registry as mreg
from juice.core import injection_calc as ic
from juice.core.experiment_spec import InjectionEvent

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")


@pytest.mark.parametrize("text, species", [("furan", "FURAN"), ("2-methylfuran", "M2F"), ("guaiacol", "GUAIACOL"),
                                           ("C4H4O", "FURAN")])
def test_addon_compounds_are_found(text, species):
    hit = ic.find_in_addons(text)
    assert hit and hit[:2] == (species, "NOAA_BB")


def test_mcm_compounds_are_not_addon_hits():
    assert ic.find_in_addons("isoprene") is None and ic.find_in_addons("toluene") is None


def test_mixture_component_in_an_addon():
    m = ic.Mixture("bb", components=[ic.MixtureComponent("furan", 1.0, "ul")])
    alerts = ic.resolve_mixture(m)
    assert m.components[0].mcm == "FURAN"
    assert any("add-on NOAA_BB" in msg for _lvl, msg in alerts)
    assert not any(lvl == "error" and "not an MCM species" in msg for lvl, msg in alerts)


def test_check_names_the_addon_to_tick():
    exp = juice.new_experiment("probe", "SAPHIR")
    exp.precursors = ["C5H8"]
    exp.injections = [InjectionEvent("furan", "FURAN", "5", "10:00", 2)]
    msgs = [m for _lvl, m in juice.check(exp) if "is not a species" in m]
    assert msgs and "add-on NOAA_BB" in msgs[0]
    exp.addons = ["NOAA_BB"]; exp.precursors = ["C5H8", "FURAN"]
    assert not [m for _lvl, m in juice.check(exp) if "is not a species" in m]


def test_starting_compounds():
    assert [n for n, _f, _note in mreg.starting_compounds("GECKOA_ethanol")] == ["C2H5OH"]
    assert ("APINENE", "C10H16") in [(n, f) for n, f, _note in mreg.starting_compounds("RACM2")]
    assert mreg.starting_compounds("MCM_v3.3.1") == []


@pytest.fixture(scope="module")
def app():
    from PyQt6.QtWidgets import QApplication
    yield QApplication.instance() or QApplication([])


def _rows(inj):
    return [[inj.comp.item(r, c).text() for c in range(6)] for r in range(inj.comp.rowCount())]


def test_compounds_table_marks_and_replaces(app):
    from juice.gui import juice_gui
    from juice.gui.experiment_wizard import _inj_add_row
    w = juice_gui.JuiceWindow("")
    chem, inj = w.page_chem, w.page_inj
    chem.base_combo.setCurrentIndex(chem.base_combo.findData("MCM_v3.3.1"))
    chem.chosen.addItem("C5H8")
    _inj_add_row(inj.table, "mixture", ["BBmix", "0.5*C5H8 + 0.3*TOLUENE + 0.2*FURAN", "20", "10:00", "2", "", "", ""])
    inj.refresh_compounds()
    rows = _rows(inj)
    assert [r[2] for r in rows] == ["C5H8", "TOLUENE", "FURAN"]
    assert [r[4] for r in rows] == ["10", "6", "4"]                                  # share x total
    assert rows[0][5].startswith("✔") and "chosen compounds" in rows[1][5] and "NOAA_BB" in rows[2][5]
    from PyQt6.QtCore import Qt
    orange = lambda r: inj.comp.item(r, 0).background().style() != Qt.BrushStyle.NoBrush
    assert [orange(r) for r in range(3)] == [False, True, True]
    assert inj.table.item(0, 2).background().color().alpha() > 0                      # the product cell too
    inj.comp.item(1, 6).setText("BENZENE")                                           # 'Use instead'
    assert inj.table.item(0, 2).text() == "0.5*C5H8 + 0.3*BENZENE + 0.2*FURAN"
    w.close()


def test_gecko_search_offers_its_starting_compounds(app):
    from juice.gui import juice_gui
    w = juice_gui.JuiceWindow("")
    chem = w.page_chem
    chem.base_combo.setCurrentIndex(chem.base_combo.findData("GECKOA_ethanol"))
    items = [chem.candidates.item(i).text() for i in range(chem.candidates.count())]
    assert len(items) == 1 and items[0].startswith("C2H5OH ")
    assert [chem.chosen.item(i).text() for i in range(chem.chosen.count())] == ["C2H5OH"]
    assert chem.tree_tabs.count() == 1                                                # its first-generation chemistry
    spec = w.collect_spec()
    assert spec.precursors == ["C2H5OH"]
    w.close()


def test_log_import_points_to_the_addon():
    from juice.importers import saphir_wiki as sw
    imp = sw.events_to_import(sw.parse_action_log("10:00 X furan injection 20 ppb\n10:30 X isoprene 5 ppb"), "2024-07-09")
    assert [(e.name, e.products) for e in imp.injections] == [("FURAN", "FURAN"), ("C5H8", "C5H8")]
    assert any("add-on NOAA_BB" in m for _lvl, m in imp.checks)


def test_chosen_compound_outside_the_mcm_is_flagged():
    exp = juice.new_experiment("probe", "SAPHIR")
    exp.precursors = ["C5H8", "FURAN"]                       # FURAN chosen, NOAA_BB not ticked
    exp.injections = [InjectionEvent("FURAN", "FURAN", "20", "09:10", 2)]
    assert any("add-on NOAA_BB" in m for _lvl, m in juice.check(exp))
