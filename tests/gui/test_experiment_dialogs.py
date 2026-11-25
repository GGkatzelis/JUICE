"""The dialogs of the experiment window, opened without a screen: the amount calculator, the mixture injection, the
SAPHIR wiki import, the chamber editor and the GECKO-A import.  Message boxes are replaced by a recorder, and the
two dialogs that write into the repository (recipes into setups/mixtures/, chamber profiles into setups/chambers/)
write into a temporary folder instead."""
import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import pytest
from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import QApplication, QDialog, QFileDialog, QInputDialog, QMessageBox

from juice.core import chambers as C
from juice.core import injection_calc as ic
from juice.core.experiment_spec import Conditions
from juice.gui import chamber_editor as ce
from juice.gui import experiment_wizard as ew

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication([])


@pytest.fixture
def boxes(monkeypatch):
    """Every message box: recorded as (kind, text); a question is answered with `boxes.answer` (Yes by default)."""
    class Seen(list):
        answer = QMessageBox.StandardButton.Yes
    seen = Seen()

    def record(kind):
        def box(*args, **kwargs):
            seen.append((kind, args[2] if len(args) > 2 else ""))
            return seen.answer if kind == "question" else QMessageBox.StandardButton.Ok
        return box
    for kind in ("warning", "information", "critical", "question"):
        monkeypatch.setattr(QMessageBox, kind, staticmethod(record(kind)))
    return seen


@pytest.fixture
def mixture_dir(tmp_path, monkeypatch):
    folder = tmp_path / "mixtures"
    folder.mkdir()
    monkeypatch.setattr(ic, "MIXTURE_DIR", str(folder))
    return folder


def said(boxes, kind, part):
    return any(k == kind and part in text for k, text in boxes)


# ── amount calculator ────────────────────────────────────────────────────────

COND = Conditions()                       # SAPHIR, 298 K, 1006.5 hPa
SAPHIR = ic.chamber("SAPHIR")


def mode(dlg, key):
    dlg.mode.setCurrentIndex([m[0] for m in ew.AmountDialog.MODES].index(key))


def test_liquid_injection_with_the_molar_mass_and_density_looked_up(app, boxes):
    dlg = ew.AmountDialog(compound="20 µl isoprene", cond=COND)
    assert (dlg.fields["mw"].text(), dlg.fields["density"].text()) == ("68.12", "0.681")      # looked up on opening
    mode(dlg, "liquid"); dlg.fields["vol_ul"].setText("20")
    expected = ic.ppb_liquid(20.0, 68.12, 0.681, COND.temperature_K, COND.pressure_hPa, SAPHIR)[0]
    assert dlg.result_ppb == pytest.approx(expected) and dlg.result.toPlainText().startswith(f"= {expected:.4g} ppb")
    assert [k for k, e in dlg.fields.items() if e.isEnabled()] == ["vol_ul", "mw", "density"]
    assert not dlg.bottle_unit.isEnabled()
    dlg._accept()
    assert QDialog.result(dlg) == QDialog.DialogCode.Accepted      # the dialog's own 'result' is its text box
    assert dlg.result_delivery["mode"] == "liquid" and dlg.result_delivery["vol_ul"] == "20"


def test_mass_flow_controller_and_gas_syringe(app, boxes):
    dlg = ew.AmountDialog(cond=COND)
    mode(dlg, "mfc")
    for key, value in (("flow_sccm", "5000"), ("duration_s", "50"), ("bottle_value", "100")):
        dlg.fields[key].setText(value)
    expected = ic.ppb_mfc(5000.0, 50.0, 1e-4, COND.temperature_K, COND.pressure_hPa, SAPHIR)[0]
    assert dlg.result_ppb == pytest.approx(expected) and dlg.bottle_unit.isEnabled()
    mode(dlg, "gas_syringe"); dlg.fields["vol_ml"].setText("3"); dlg.fields["bottle_value"].setText("1"); dlg.bottle_unit.setCurrentText("%")
    assert dlg.result_ppb == pytest.approx(ic.ppb_gas_syringe(3.0, 0.01, COND.temperature_K, COND.pressure_hPa, SAPHIR)[0])


def test_an_earlier_delivery_is_shown_again(app, boxes):
    dlg = ew.AmountDialog(delivery={"mode": "ppb", "value": "7.5"}, cond=COND)
    assert ew.AmountDialog.MODES[dlg.mode.currentIndex()][0] == "ppb" and dlg.result_ppb == 7.5
    assert [k for k, e in dlg.fields.items() if e.isEnabled()] == ["value"]


def test_nothing_to_compute_is_not_accepted_and_an_unknown_compound_is_said(app, boxes):
    dlg = ew.AmountDialog(cond=COND)
    mode(dlg, "liquid")
    dlg._accept()
    assert QDialog.result(dlg) != QDialog.DialogCode.Accepted and dlg.result_delivery == {}
    assert said(boxes, "warning", "missing or invalid field")
    dlg.compound.setText("foo")
    dlg._lookup()
    assert said(boxes, "information", "'foo' not found")


# ── mixture injection ────────────────────────────────────────────────────────

RECIPE = "* Limonene 10 μl\nethanol, 30, ul\ntexanol 5 mg\n"


def table(dlg):
    return [[dlg.table.item(r, c).text() if dlg.table.item(r, c) else "" for c in range(8)] for r in range(dlg.table.rowCount())]


def test_a_pasted_recipe_becomes_one_injection_of_mole_fractions(app, boxes, mixture_dir):
    dlg = ew.MixtureDialog(cond=COND)
    dlg.name.setText("Cleaner"); dlg.aliases.setText("cleaner, spray")
    dlg.paste.setPlainText(RECIPE); dlg._parse()
    rows = table(dlg)
    assert [(r[0], r[1], r[2], r[3]) for r in rows] == [("Limonene", "10", "ul", "LIMONENE"), ("ethanol", "30", "ul", "C2H5OH"),
                                                         ("texanol", "5", "mg", "")]
    preview = dlg.preview.toPlainText()
    assert "products: " in preview and "LIMONENE" in preview and "texanol: no MCM species" in preview
    dlg.injected.setText("80"); dlg.start.setText("10:30"); dlg.dur.setText("3")
    dlg._accept()
    assert dlg.result() == QDialog.DialogCode.Accepted
    ev = dlg.result_event
    assert (ev.name, ev.start, ev.duration_min, ev.action) == ("Cleaner", "10:30", 3.0, "mixture")
    assert set(dlg.result_precursors) == {"LIMONENE", "C2H5OH"} and "C2H5OH" in ev.products
    saved = ic.load_mixtures(str(mixture_dir))
    assert list(saved) == ["Cleaner"] and saved["Cleaner"].aliases == ["cleaner", "spray"]


def test_the_injected_volume_scales_the_amount(app, boxes, mixture_dir):
    dlg = ew.MixtureDialog(cond=COND)
    dlg.paste.setPlainText(RECIPE); dlg._parse()
    total = lambda: float(dlg._compute()[0].total_ppb)
    dlg.injected.setText("40")                              # the recipe's own volume (10 + 30 µl; texanol is weighed)
    at_40 = total()
    dlg.injected.setText("80")
    assert total() == pytest.approx(2 * at_40, rel=1e-3)


def test_a_surrogate_is_found_by_formula_and_used(app, boxes, mixture_dir):
    if not ic.MCM_DB:
        pytest.skip("the MCM species table is not there")
    dlg = ew.MixtureDialog(cond=COND)
    dlg.paste.setPlainText(RECIPE); dlg._parse()
    dlg.table.setCurrentCell(2, 0)                          # texanol: no MCM species, so the finder starts
    assert dlg.sur_query.text() == "C12H24O3" and dlg.sur_list.count() > 0
    names = [dlg.sur_list.item(i).data(Qt.ItemDataRole.UserRole) for i in range(dlg.sur_list.count())]
    dlg._use_surrogate(dlg.sur_list.item(names.index(next(n for n in names if n))))
    assert table(dlg)[2][4] == next(n for n in names if n)
    assert "NOT modelled" not in dlg.preview.toPlainText()


def test_a_recipe_without_anything_to_model_is_refused(app, boxes, mixture_dir):
    dlg = ew.MixtureDialog(cond=COND)
    dlg.paste.setPlainText("texanol 5 mg\n"); dlg._parse()
    dlg._accept()
    assert dlg.result() != QDialog.DialogCode.Accepted and said(boxes, "warning", "Fix the ✖ alerts first")
    assert list(ic.load_mixtures(str(mixture_dir))) == []              # nothing saved either


def test_saved_recipes_are_offered_and_loaded(app, boxes, mixture_dir):
    ic.save_mixture(ic.Mixture("Kitchen", aliases=["cooking"], components=ic.parse_mixture_text("ethanol 5 ul\n")))
    dlg = ew.MixtureDialog(cond=COND)
    assert [dlg.saved.itemText(i) for i in range(dlg.saved.count())] == ["— saved mixtures —", "Kitchen"]
    dlg.saved.setCurrentText("Kitchen")
    assert dlg.name.text() == "Kitchen" and dlg.aliases.text() == "cooking" and table(dlg)[0][:4] == ["ethanol", "5", "ul", "C2H5OH"]


def test_preload_shows_a_recipe_with_its_notes(app, boxes, mixture_dir):
    dlg = ew.MixtureDialog(cond=COND)
    dlg.preload(ic.Mixture("VCP", aliases=["VCP solution"], components=ic.parse_mixture_text(RECIPE)), "120", "09:15", "4")
    assert (dlg.name.text(), dlg.injected.text(), dlg.start.text(), dlg.dur.text()) == ("VCP", "120", "09:15", "4")
    assert "not an MCM species" in table(dlg)[2][7]


# ── SAPHIR wiki import ───────────────────────────────────────────────────────

LOG = """2024-07-09 action log
06:00 H F.R. CHANGE: humidification start
06:40 h F.R. CHANGE: humidification stop
08:00 R F.R. CHANGE: open roof
09:00 M M.S. CHANGE: 2 ppb NO
09:30 $ M.S. instrument check
10:00 X M.S. CHANGE: 20 µl isoprene injected
14:00 r F.R. CHANGE: close roof
"""


@pytest.fixture
def wiki(app, boxes, mixture_dir):
    return ew.WikiImportDialog(date="2024-07-09", archive_root=os.path.join(HERE, "no_archive_here"))


def ticks(dlg):
    return [dlg.table.item(r, 0).checkState() == Qt.CheckState.Checked for r in range(dlg.table.rowCount())]


def test_a_pasted_log_is_read_and_imported(wiki):
    wiki.paste.setPlainText(LOG); wiki._parse_paste()
    assert wiki.table.rowCount() == 7 and ticks(wiki) == [True, True, True, True, False, True, True]      # notes unticked
    assert wiki.status.text().startswith("7 action-log lines from pasted text")
    preview = wiki.preview.toPlainText()
    assert preview.startswith("Date 2024-07-09;  roof open: 08:00–14:00") and "day folder not found" in preview
    assert "Precursors (MCM): C5H8" in preview
    wiki._accept()
    imp = wiki.result_import
    assert wiki.result() == QDialog.DialogCode.Accepted and imp.roof_open == [("08:00", "14:00")]
    assert [e.name for e in imp.injections] == ["H2O", "NO", "C5H8"] and imp.precursors == ["C5H8"]


def test_an_unticked_line_is_left_out(wiki):
    wiki.paste.setPlainText(LOG); wiki._parse_paste()
    wiki.table.item(3, 0).setCheckState(Qt.CheckState.Unchecked)       # the NO injection
    assert [e.name for e in wiki._current_import().injections] == ["H2O", "C5H8"]


def test_the_dialog_reopens_as_it_was(wiki, app):
    wiki.paste.setPlainText(LOG); wiki._parse_paste()
    wiki.table.item(3, 0).setCheckState(Qt.CheckState.Unchecked); wiki.replace_box.setChecked(False); wiki.user.setText("me")
    memory = wiki.memory()
    assert "pw" not in memory and "password" not in str(memory).lower()
    again = ew.WikiImportDialog()
    again.restore(memory)
    assert again.date.text() == "2024-07-09" and again.user.text() == "me" and not again.replace_box.isChecked()
    assert ticks(again) == ticks(wiki) and again.status.text().endswith("restored from the last import of this experiment")
    ew.WikiImportDialog().restore({})                                   # nothing to restore: nothing happens


def test_a_saved_page_is_loaded_and_its_date_taken_from_the_name(app, boxes, mixture_dir, tmp_path, monkeypatch):
    page = tmp_path / "2023-05-17.html"
    page.write_text("<html><p>08:00 R open roof</p><p>12:00 r close roof</p></html>", encoding="utf-8")
    monkeypatch.setattr(QFileDialog, "getOpenFileName", lambda *a, **k: (str(page), ""))
    dlg = ew.WikiImportDialog()
    dlg._load_file()
    assert dlg.date.text() == "2023-05-17" and dlg.table.rowCount() == 2 and "2023-05-17.html" in dlg.status.text()


def test_fetching_needs_a_date_and_reports_a_failure(wiki, boxes):
    wiki.date.setText("9 July")
    wiki._fetch()
    assert said(boxes, "warning", "YYYY-MM-DD")
    wiki._fetch_failed("The wiki page needs a login")
    assert wiki.status.text() == "Fetch failed: The wiki page needs a login" and wiki.fetch_btn.isEnabled()
    wiki._fetched(LOG)
    assert wiki.paste.toPlainText() == LOG and wiki.table.rowCount() == 7


def test_selecting_a_row_marks_its_line_in_the_original_text(wiki):
    wiki.paste.setPlainText(LOG); wiki._parse_paste()
    wiki._show_source_line(5)
    (selection,) = wiki.paste.extraSelections()
    assert selection.cursor.selectedText() == "10:00 X M.S. CHANGE: 20 µl isoprene injected"
    wiki._show_source_line(99)                                          # outside the table: nothing happens


def test_nothing_parsed_cannot_be_used_and_a_mixture_offers_its_surrogates(wiki, boxes, mixture_dir):
    wiki._accept()
    assert wiki.result_import is None and said(boxes, "information", "Nothing parsed yet")
    ic.save_mixture(ic.Mixture("VCP2", aliases=["VCP"], components=ic.parse_mixture_text("limonene 10 ul\nethanol 30 ul\n")))
    wiki.paste.setPlainText("10:00 X M.S. CHANGE: VCP solution (20 ul)\n"); wiki._parse_paste()
    assert wiki.sur_btn.isEnabled() and "Mixture recipe(s) used: VCP2" in wiki.preview.toPlainText()


# ── chamber editor ───────────────────────────────────────────────────────────

@pytest.fixture
def chamber_dir(tmp_path, monkeypatch):
    """save_profile writes into a temporary folder and does not reload the installed profiles."""
    folder = tmp_path / "chambers"
    written = []

    def save(p, directory=str(folder)):
        os.makedirs(directory, exist_ok=True)
        path = os.path.join(directory, p.id + C.SUFFIX)
        with open(path, "w", encoding="utf-8") as fh:
            import json
            json.dump(C.profile_to_dict(p), fh)
        written.append(p)
        return path
    monkeypatch.setattr(C, "save_profile", save)
    return written


def test_saphir_opens_read_only(app, boxes):
    dlg = ce.ChamberDialog(C.chamber("SAPHIR"))
    assert dlg.locked and not dlg.save_btn.isEnabled() and not dlg.vol.isEnabled() and dlg.save_as_btn.isEnabled()


def test_a_new_lamp_chamber_is_made_from_a_copy(app, boxes, chamber_dir):
    dlg = ce.ChamberDialog(C.chamber("SAPHIR"), as_new=True)
    assert dlg.id.isEnabled() and dlg.id.text() == "" and dlg.save_as_btn.isHidden()
    dlg.id.setText("my lab box"); dlg.name.setText("My lab box"); dlg.vol.setText("2,5"); dlg.leak.setText("0.1")
    assert dlg.dil_preview.text().startswith("DILUTE = ")
    dlg.ltype.setCurrentIndex(C.LIGHT_TYPES.index("lamps"))
    assert dlg.lat.isHidden() and not dlg.jno2.isHidden() and not dlg.builtin.isHidden()
    dlg.jno2.setText("2e-3"); dlg.builtin.setCurrentIndex(1)
    dlg.jratios.insertRow(0)
    from PyQt6.QtWidgets import QTableWidgetItem
    dlg.jratios.setItem(0, 0, QTableWidgetItem("1")); dlg.jratios.setItem(0, 1, QTableWidgetItem("5e-4"))
    dlg._add_wall("O3", "2.1e-6")
    dlg._save()
    (saved,) = chamber_dir
    assert dlg.result() == QDialog.DialogCode.Accepted and dlg.saved_id == saved.id == C.safe_id("my lab box")
    assert (saved.name, saved.volume_m3, saved.leak_flow_m3h, saved.light_type) == ("My lab box", 2.5, 0.1, "lamps")
    assert saved.light["jno2_full_s1"] == 2e-3 and saved.light["spectrum_builtin"] == "blacklight_350"
    assert saved.light["j_ratios"] == {"J1": 5e-4} and saved.wall_loss_s1["O3"] == 2.1e-6


def test_what_a_profile_needs_before_it_is_saved(app, boxes, chamber_dir):
    dlg = ce.ChamberDialog(C.chamber("SAPHIR"), as_new=True)
    dlg.id.setText("box")
    dlg._save()
    # bug of 25 Nov 2026: a new chamber without a name took the name of the chamber it was copied from ("SAPHIR")
    assert said(boxes, "warning", "Give the chamber a name") and chamber_dir == []
    dlg.id.setText("")
    dlg.name.setText("Box"); dlg.vol.setText("0")
    dlg._save()
    assert said(boxes, "warning", "volume must be a positive number")
    dlg.vol.setText("ten")
    dlg._save()
    assert said(boxes, "warning", "A number could not be read") and dlg.dil_preview.text().endswith("is not a number)")
    dlg.vol.setText("10"); dlg._save()
    assert said(boxes, "warning", "Give the new chamber an id")
    dlg.id.setText("SAPHIR"); dlg._save()
    assert said(boxes, "warning", "The id SAPHIR is reserved")
    dlg.id.setText("box"); dlg.ltype.setCurrentIndex(C.LIGHT_TYPES.index("lamps"))
    boxes.answer = QMessageBox.StandardButton.No
    dlg._save()
    assert said(boxes, "question", "No jNO2 at full lamp power") and chamber_dir == []


def test_saphir_is_saved_as_a_new_unverified_chamber(app, boxes, chamber_dir, monkeypatch):
    monkeypatch.setattr(QInputDialog, "getText", staticmethod(lambda *a, **k: ("saphir_like", True)))
    dlg = ce.ChamberDialog(C.chamber("SAPHIR"))
    dlg._save_as()
    (saved,) = chamber_dir
    assert saved.id == "saphir_like" and saved.name == "saphir_like" and not saved.verified
    assert "Derived from the SAPHIR profile" in saved.note


def test_the_j_ratios_of_a_light_source_can_be_previewed(app, boxes):
    dlg = ce.ChamberDialog(C.chamber("SAPHIR"), as_new=True)
    dlg.ltype.setCurrentIndex(C.LIGHT_TYPES.index("lamps")); dlg.solar_like.setChecked(True)
    dlg._preview_ratios()
    assert dlg.ratio_lab.text() and not dlg.ratio_lab.text().startswith("cannot compute")


# ── GECKO-A import ───────────────────────────────────────────────────────────

ETHANOL = os.path.join(HERE, "fixtures", "gecko_ethanol")


def test_a_generator_folder_is_described_and_named(app, boxes):
    dlg = ew.GeckoImportDialog(start_dir=ETHANOL)
    assert "106 species, 366 reactions, generator maxgen 6; primary O02000 CH3CH2(OH)" == dlg.summary.text()
    assert dlg.voc.text() == "C2H5OH" and dlg.base_id.text() == "GECKOA_C2H5OH"
    dlg.voc.setText("ETOH")
    assert dlg.base_id.text() == "GECKOA_ETOH"
    dlg.folder.setText(HERE)
    assert dlg.summary.text().startswith("not a generator output folder")


def test_the_import_starts_only_with_a_folder_and_a_name_and_reports_its_result(app, boxes, monkeypatch):
    started = []

    class FakeThread:
        def __init__(self, *args):
            started.append(args)
            from types import SimpleNamespace as NS
            self.message = NS(connect=lambda f: None)
            self.done = NS(connect=lambda f: setattr(self, "_done", f))
            self.failed = NS(connect=lambda f: setattr(self, "_failed", f))

        def start(self):
            self._done("GECKOA_C2H5OH")
    monkeypatch.setattr(ew, "_GeckoImportThread", FakeThread)
    dlg = ew.GeckoImportDialog()
    dlg._start()
    assert said(boxes, "warning", "Choose a generator output folder") and not started
    dlg.folder.setText(ETHANOL); dlg.voc.setText("")
    dlg._start()
    assert said(boxes, "warning", "Give the primary species a name") and not started
    dlg.voc.setText("C2H5OH"); dlg.maxgen.setValue(3)
    dlg._start()
    (args,) = started
    assert args == (ETHANOL, "GECKOA_C2H5OH", "C2H5OH", 3, None, "krechmer")
    assert dlg.result_id == "GECKOA_C2H5OH" and dlg.result() == QDialog.DialogCode.Accepted
    dlg2 = ew.GeckoImportDialog(start_dir=ETHANOL)
    dlg2.buttons.setEnabled(False)
    dlg2._failed("ValueError: broken file")
    assert dlg2.buttons.isEnabled() and "failed: ValueError: broken file" in dlg2.log.toPlainText()


# ── a wiki import put into the experiment window ─────────────────────────────

def test_a_wiki_import_fills_the_experiment_window(app, boxes, mixture_dir, tmp_path):
    from juice.gui import juice_gui
    from juice.importers import saphir_wiki
    w = juice_gui.JuiceWindow("")
    try:
        w.page_start.name_edit.setText("my_experiment")
        cond = w.page_cond.get(); cond.archive_root = str(tmp_path / "archive"); w.page_cond.set_from(cond)
        imp = saphir_wiki.events_to_import(saphir_wiki.parse_action_log(LOG), "2024-07-09", mixtures={})
        w.apply_wiki_import(imp)
        spec = w.collect_spec()
        assert spec.name == "saphir_20240709"                          # the placeholder name is replaced by the date
        assert spec.conditions.date == "2024-07-09" and spec.conditions.roof_open == [("08:00", "14:00")]
        assert (spec.conditions.source, spec.conditions.j_source) == ("manual", "clear_sky")     # no day folder in the archive
        assert "C5H8" in spec.precursors and [e.name for e in spec.injections] == ["H2O", "NO", "C5H8"]
        assert w.import_notes == [m for _l, m in imp.checks]
        more = saphir_wiki.events_to_import(saphir_wiki.parse_action_log("11:00 O F.R. CHANGE: 50 ppb ozone\n"), "2024-07-09", mixtures={})
        w.apply_wiki_import(more, replace=False)                        # appended, not replacing
        assert [e.name for e in w.collect_spec().injections] == ["H2O", "NO", "C5H8", "O3"]
        assert w.collect_spec().name == "saphir_20240709"               # a name already set is kept
        (tmp_path / "archive" / "2024" / "07" / "09").mkdir(parents=True)
        w.apply_wiki_import(imp)
        assert (w.collect_spec().conditions.source, w.collect_spec().conditions.j_source) == ("archive", "measured")
    finally:
        w.dirty = False
        w.close()
