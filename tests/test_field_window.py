"""The first start-menu question and the field window, run without a screen."""
import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import pytest
from PyQt6.QtWidgets import QApplication, QMessageBox

import field_spec as fs
import field_window as fw
import juice_start as js


@pytest.fixture(scope="module")
def app():
    yield QApplication.instance() or QApplication([])


@pytest.fixture
def calls(monkeypatch):
    seen = []
    monkeypatch.setattr(js, "launch", lambda kind, args=(), toml="": seen.append((kind, list(args))) or True)
    return seen


def _titles(menu):
    from juice_start import _Choice
    return [c.findChildren(js.QLabel)[0].text() for c in menu.findChildren(_Choice)]


def test_first_choice_leads_to_todays_chamber_menu(app, calls):
    first = js.FirstChoice()
    assert _titles(first) == ["Chamber experiment", "Field observations"]
    menu = first.open_menu(js.StartWindow)
    assert _titles(menu) == ["New experiment", "Open a saved experiment", "New study", "Import a FACSIMILE experiment",
                             "Results of one experiment", "Results of a study"]
    menu.close()


def test_chamber_menu_without_back_has_no_back_button(app):
    labels = [b.text() for b in js.StartWindow().findChildren(js.QPushButton)]
    assert labels == ["Close"]


def test_field_menu_starts_a_field_window(app, calls):
    menu = js.FieldMenu("", on_back=lambda: None)
    assert _titles(menu)[:2] == ["New field case", "Open a saved field case"]
    from juice_start import _Choice
    menu.findChildren(_Choice)[0].clicked.emit()
    assert calls == [("field", [])]
    assert "Back" in [b.text() for b in menu.findChildren(js.QPushButton)]


def test_field_file_opens_the_field_window(calls):
    assert js.open_file(r"C:\x\la.field.json") is True
    assert calls == [("field", [r"C:\x\la.field.json"])]


def test_other_files_still_open_their_windows(calls):
    js.open_file("a.experiment.json"); js.open_file("a.study.json"); js.open_file("a.cs")
    assert [c[0] for c in calls] == ["experiment", "study", "results"]


def _window(tmp_path, run_type="lagrangian", area=64.0):
    table = fs.write_template(str(tmp_path / "t.csv"), run_type)
    w = fw.FieldWindow()
    w.run_buttons[run_type].setChecked(True)
    w.name.setText("case1"); w.table.setText(table); w.page_chem.chosen.addItem("C5H8"); w.area.setValue(area)
    w.refresh_all()
    return w, table


def test_window_has_the_five_steps(app):
    w = fw.FieldWindow()
    assert [w.steps.item(i).text() for i in range(w.steps.count())] == \
           ["1   Field case", "2   Chemistry", "3   Emissions and constraints", "4   Mixing and deposition", "5   Run"]
    assert w.stack.count() == 5 and not w.page_run.run_btn.isEnabled()
    assert isinstance(w.page_chem, fw.ChemistryPage)
    w._step(+1); assert w.steps.currentRow() == 1 and w.stack.currentIndex() == 1
    w._step(-1); w._step(-1); assert w.steps.currentRow() == 0
    w.dirty = False


def test_window_checks_a_template(app, tmp_path):
    w, _table = _window(tmp_path, area=0.0)
    assert w.dilution.currentData() == "pbl"
    w.steps.setCurrentRow(fs.STEP_SOURCES)
    assert not w.plot_sources.fig.axes[0].lines                  # no plot from a division by a missing area
    w.steps.setCurrentRow(fs.STEP_SOURCES)
    assert not w.plot_sources.fig.axes[0].lines                  # no plot from a division by a missing area
    rep = w.case[2]
    assert any("footprint area" in e for e in rep.errors)
    assert w.steps.item(fs.STEP_SOURCES).text().endswith("✖")
    w.area.setValue(64.0)
    rep = w.refresh_all()
    assert rep.ok, rep.text()
    assert w.steps.item(fs.STEP_SOURCES).text().endswith("✔")
    roles = [w.review.preview.item(r, 0).text() for r in range(w.review.preview.rowCount())]
    assert "emission" in roles and "background" in roles
    assert w.emit_table.rowCount() == 4 and w.back_table.rowCount() == 2
    w.dirty = False


def test_clicking_a_problem_opens_its_step(app, tmp_path):
    w, _table = _window(tmp_path, area=0.0)
    item = next(w.review.checks.item(i) for i in range(w.review.checks.count()) if "footprint" in w.review.checks.item(i).text())
    w.review.checks.itemClicked.emit(item)
    assert w.steps.currentRow() == fs.STEP_SOURCES and w.stack.currentIndex() == fs.STEP_SOURCES
    w.dirty = False


def test_blank_window_marks_nothing(app):
    w = fw.FieldWindow()
    assert not w.dirty and w.steps.item(0).text() == "1   Field case"
    assert "Nothing entered" in w.review.summary.text()


def test_steps_show_what_belongs_to_the_kind_of_run(app):
    w = fw.FieldWindow()
    assert w.emit_box.isHidden() and not w.life_box.isHidden() and w.init_box.isHidden()
    w.run_buttons["lagrangian"].setChecked(True)
    assert not w.emit_box.isHidden() and w.life_box.isHidden() and not w.init_box.isHidden()
    assert w.mix_form.isRowVisible(w.k_collapse) and not w.mix_form.isRowVisible(w.tracer)
    w.dilution.setCurrentIndex(w.dilution.findData("tracer"))
    assert w.mix_form.isRowVisible(w.tracer) and not w.mix_form.isRowVisible(w.k_collapse)
    assert w.kind_box.title() == "What this run needs: Lagrangian air parcel"
    w.dirty = False


@pytest.mark.parametrize("run_type", list(fs.RUN_TYPES))
def test_every_plot_draws(app, tmp_path, run_type):
    w, _table = _window(tmp_path, run_type)
    for step, plot in ((0, w.plot_conditions), (2, w.plot_sources), (3, w.plot_mixing)):
        w.steps.setCurrentRow(step)
        assert plot.plotted and plot.fig.axes, plot.mode
    w.dirty = False


def test_plots_survive_a_table_they_cannot_use(app, tmp_path):
    w, _table = _window(tmp_path, "lagrangian")
    p = tmp_path / "broken.csv"
    p.write_text("time,T,p,RH,E_NO [mol/s]\n2021-08-15 12:00,298,1000,50,1\n2021-08-15 13:00,298,1000,50,1\n", encoding="utf-8")
    w.table.setText(str(p))
    for step in (0, 2, 3):
        w.steps.setCurrentRow(step)
    assert w.plot_sources.plotted and not w.refresh_all().ok
    w.table.setText(""); w.steps.setCurrentRow(0)
    assert w.plot_conditions.plotted
    w.dirty = False


def test_save_and_reopen_restores_every_step(app, tmp_path):
    table = fs.write_template(str(tmp_path / "t.csv"), "plume")
    w = fw.FieldWindow()
    w.run_buttons["plume"].setChecked(True)
    w.name.setText("plume1"); w.table.setText(table); w.lat.setText("34,14"); w.lon.setText("-118.13"); w.step_s.setValue(30)
    w.page_chem.chosen.addItem("C5H8")
    fw._add_row(w.initials, ["o3", "20"]); fw._add_row(w.deposition, ["O3", "0.5", ""]); fw._add_row(w.deposition, ["HNO3", "", "2"])
    path = str(tmp_path / "plume1" / "plume1.field.json")
    assert w.write_to(path) and not w.dirty
    again = fw.FieldWindow(path)
    spec, problems = again.spec_from_form()
    assert not problems and not again.dirty
    assert (spec.run_type, spec.lat, spec.lon, spec.initials, spec.step_s, spec.dilution) == ("plume", 34.14, -118.13, {"O3": 20.0}, 30, "table")
    assert spec.deposition_velocity_cm_s == {"O3": 0.5} and spec.deposition_lifetime_h == {"HNO3": 2.0}
    assert spec.chemistry["precursors"] == ["C5H8"] and os.path.samefile(spec.table, table)
    assert again.folder.text() == str(tmp_path)


def test_existing_mechanism_round_trip(app, tmp_path):
    w, _table = _window(tmp_path, "constrained")
    w.page_chem.mode_existing.setChecked(True); w.page_chem.mech_combo.setCurrentText("RACM2B_VCP")
    path = str(tmp_path / "c.field.json")
    assert w.write_to(path)
    again = fw.FieldWindow(path)
    assert again.page_chem.mode_existing.isChecked() and again.spec_from_form()[0].chemistry["mechanism_stem"] == "RACM2B_VCP"


def test_save_uses_the_project_folder(app, tmp_path):
    w, _table = _window(tmp_path, "constrained")
    w.folder.setText(str(tmp_path / "project"))
    assert w.save()
    assert os.path.isfile(tmp_path / "project" / "case1" / "case1.field.json") and not w.dirty


def test_form_problems_block_saving_and_are_listed(app, tmp_path, monkeypatch):
    shown = []
    monkeypatch.setattr(QMessageBox, "warning", lambda *a: shown.append(a[2]))
    w = fw.FieldWindow()
    w.lat.setText("north"); fw._add_row(w.deposition, ["O3", "1", "2"]); fw._add_row(w.initials, ["O3", "x"])
    rep = w.refresh_all()
    text = rep.text()
    assert "Latitude" in text and "one of the two" in text and "not a number" in text
    assert not w.write_to(str(tmp_path / "x.field.json")) and "Latitude" in shown[0]
    assert not (tmp_path / "x.field.json").exists()
    w.dirty = False


def test_save_without_a_name_asks_for_one(app, monkeypatch):
    shown = []
    monkeypatch.setattr(QMessageBox, "warning", lambda *a: shown.append(a[2]))
    w = fw.FieldWindow()
    assert not w.save() and "name" in shown[0]


def test_broken_case_file_is_reported(app, tmp_path, monkeypatch):
    shown = []
    monkeypatch.setattr(QMessageBox, "warning", lambda *a: shown.append(a[2]))
    p = tmp_path / "x.field.json"
    p.write_text("not json", encoding="utf-8")
    w = fw.FieldWindow()
    assert not w.open_path(str(p)) and shown and w.path == ""
