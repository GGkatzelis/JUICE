"""The experiment window (juice_gui.JuiceWindow) opened without a screen, on a small real experiment (the UV-lamp
chamber with ethanol) whose project folder is a temporary one: the review panel and its guidance, the input
plots, the step marks, and the Run step from 'Prepare' to a finished run and the comparison with a reference.
Starting another program (viewer, study window, file browser) is recorded instead of done."""
import json
import os
import shutil

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import numpy as np
import pytest
from PyQt6.QtCore import QProcess
from PyQt6.QtWidgets import QApplication, QDialog, QLabel, QMessageBox

import juice
from juice.core.experiment_spec import InjectionEvent
from juice.gui import juice_gui as jg

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
TOML = os.path.join(REPO, "settings.toml")


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication([])


@pytest.fixture
def started(monkeypatch):
    """Programs the window would start: (program, arguments)."""
    calls = []
    monkeypatch.setattr(QProcess, "start", lambda self, prog, args=(): calls.append((prog, list(args))))
    monkeypatch.setattr(jg, "_open_path", lambda path: calls.append(("open", [path])))
    return calls


def experiment(name="uv_ethanol"):
    spec = juice.new_experiment(name, template="uv_chamber_lamps")
    c = spec.conditions
    c.start, c.end = "08:00", "11:00"
    c.archive_root = ""                       # no archive: this machine may see the SAPHIR drive, the tests must not
    spec.precursors = ["C2H5OH"]
    spec.initials = list(spec.initials) + [("NO2", "M*1e-9*20")]
    spec.injections = [InjectionEvent("ethanol", "C2H5OH", "200", "09:30", 2.0, action="voc")]
    return spec


@pytest.fixture
def window(app, tmp_path):
    w = jg.JuiceWindow(TOML)
    w.load_spec(experiment())
    w.page_start.out_edit.setText(str(tmp_path))
    w.refresh_all()
    yield w
    w.dirty = False
    w.close()


# ── small helpers of the module ──────────────────────────────────────────────

def test_a_blank_experiment_and_what_counts_as_entered():
    blank = jg.blank_spec()
    assert jg.is_blank(blank) and blank.name == "" and blank.conditions.roof_open == []
    blank.name = "named only"
    assert jg.is_blank(blank)                                           # a name alone is not an experiment
    blank.conditions.date = "2024-07-09"
    assert not jg.is_blank(blank)


@pytest.mark.parametrize("text, value", [("2.6*0.7", 1.82), ("1.0D-3", 1e-3), ("2@3", 8.0), ("TEMP*2", 5.0), ("", 5.0)])
def test_small_constant_expressions(text, value):
    assert jg._const(text, 5.0) == pytest.approx(value)


def test_the_list_of_what_goes_into_the_files():
    spec = experiment()
    spec.reaction_adjustments = [("C2H5OH + OH --> X", "0", ""), ("A --> B", "1e-3", "")]
    bullets = jg.spec_bullets(spec)
    assert bullets[0].startswith("Chemistry: MCM_v3.3.1 subset for C2H5OH")
    assert "Reaction adjusted: C2H5OH + OH --> X  →  switched off" in bullets and "Reaction adjusted: A --> B  →  k = 1e-3" in bullets
    assert any(b.startswith("Time: 2024-07-09 08:00–11:00 UTC in 60 s steps; roof open 09:00–13:00") for b in bullets)
    assert any(b.startswith("Injection 09:30 VOC injection: 200 ppb of C2H5OH") for b in bullets)
    assert bullets[-1] == "Initial values: 2 species set, all others at the floor value"
    spec.precursors, spec.inorganic_only = [], True
    assert jg.spec_bullets(spec)[0].startswith("Chemistry: inorganic only")
    spec.inorganic_only = False
    assert jg.spec_bullets(spec)[0] == "Chemistry: not chosen yet"


# ── review panel ─────────────────────────────────────────────────────────────

def items(lst):
    return [lst.item(i).text() for i in range(lst.count())]


def test_the_review_of_a_complete_experiment(window):
    rev = window.review
    assert rev.summary.text().startswith("<b>uv_ethanol</b> — 2024-07-09 08:00–11:00 UTC, manual conditions, step 60 s")
    assert rev.checks_title.text() == "<b>Needs attention</b> — nothing" and items(rev.checks)[0] == "✔ No problems found"
    assert any("C2H5OH" in b for b in items(rev.bullets))
    assert [window.steps.item(i).text() for i in range(4)] == ["1   Experiment  ✔", "2   Chemistry  ✔", "3   Injections  ✔", "4   Chamber  ✔"]


def test_the_review_of_a_blank_window(app):
    w = jg.JuiceWindow(TOML)
    try:
        assert w.review.summary.text().startswith("<i>Nothing entered yet.</i>") and w.review.checks.count() == 0
        assert w.steps.item(0).text() == "1   Experiment"                  # no mark before anything is entered
    finally:
        w.dirty = False; w.close()


def test_problems_are_listed_and_mark_their_step(window):
    window.page_cond.date.setText("")
    window.page_chem.chosen.clear()
    window.refresh_all()
    rev = window.review
    assert "to fix before running" in rev.checks_title.text()
    texts = items(rev.checks)
    assert texts[0].startswith("✖") and any("Enter the experiment date" in t for t in texts) and any("No chemistry" in t for t in texts)
    assert window.steps.item(0).text().endswith("✖") and window.steps.item(1).text().endswith("✖")


@pytest.mark.parametrize("message, step", [
    ("Enter the experiment date and the time window (Where and when).", 0),
    ("Time window is empty — check start and end time.", 0),
    ("Measured inputs from a file: choose the file (CSV or Excel; 'Template…' writes an example).", 0),
    ("Photolysis: measured J values selected but none found — the J columns would be 0.", 0),
    ("Inputs: temperature not measured — constant value(s) used.", 0),
    ("Held species without measurements: NO", 0),
    ("Choose where the chamber conditions come from: measured (archive), measured (file) or manual.", 0),
    ("Measured J values selected without archive data: the J columns will be 0.", 0),
    ("Roof never opens: run as a dark experiment (no photolysis).", 0),
    ("No chemistry: choose compounds, an existing mechanism, or 'Inorganic chemistry only'.", 1),
    ("Injection 'ethanol': product 'XYZ' is not a species of the mechanism.", 1),
    ("Injection 'ethanol': product 'XYZ' is not a species of the mechanism GECKOA_x (used as a whole, so it cannot be added).", 1),
    ("Injection 'ethanol': start time 'soon' is not HH:MM.", 2),
    ("Injection 'ethanol': amount '2*x' is an expression, it will be written verbatim.", 2),
    ("Injection 'fan': neither modelled nor timed — it has no effect.", 2),
    ("No injections defined.", 2),
    ("Import: 10:00: texanol: no MCM species and no surrogate — injected but NOT modelled", 2),
    ("Import: 10:00: roof closed without a preceding open", 2),
])
def test_every_message_leads_to_its_step(window, message, step):
    got, hint, _act = window.review._guide(message)
    assert got == step and hint


def test_clicking_a_message_opens_its_step_and_says_how_to_fix_it(window):
    from PyQt6.QtWidgets import QListWidgetItem
    from PyQt6.QtCore import Qt
    item = QListWidgetItem("✖ x"); item.setData(Qt.ItemDataRole.UserRole, ("error", "Injection 'ethanol': start time 'soon' is not HH:MM."))
    window.review._guide_clicked(item)
    assert window.steps.currentRow() == 2 and window.review.fix.text().startswith("<b>How to fix</b> (Injections step)")
    assert window.page_inj.table.currentRow() == 0                       # the ethanol row is selected
    item.setData(Qt.ItemDataRole.UserRole, ("warning", "something nobody wrote guidance for"))
    window.review._guide_clicked(item)
    assert window.review.fix.text() == "No guidance for this message yet."


def test_a_species_from_an_add_on_is_ticked_by_the_guidance(window):
    from juice.chemistry import mechanism_registry as mreg
    ids = list(mreg.load_addons())
    if not ids:
        pytest.skip("no add-on installed")
    step, _hint, act = window.review._guide(f"Injection 'x': product 'FURAN' is not a species of the mechanism; it is in the add-on {ids[0]}: tick it")
    assert step == 1
    act()
    assert ids[0] in window.page_chem.addon_ids()


# ── input plots ──────────────────────────────────────────────────────────────

@pytest.mark.parametrize("plot, n_axes", [("plot_cond", 4), ("plot_inj", 2), ("plot_chamber", 2)])
def test_each_step_plot_draws(window, plot, n_axes):
    p = getattr(window, plot)
    p.plot()
    assert p.plotted and len(p.fig.axes) == n_axes
    assert all(ax.get_xlabel() == "hour (UTC)" for ax in p.fig.axes if ax.get_xlabel())


def test_the_review_plot_stacks_conditions_and_injections(window):
    window.review.tabs.setCurrentIndex(1)                               # the Inputs plot tab draws on its first visit
    assert window.review.plot.plotted and len(window.review.plot.fig.axes) >= 3


def test_the_injection_plot_says_when_there_is_no_measured_reactivity(window):
    window.plot_inj.plot()
    assert "no measured kOH for this day" in window.plot_inj.note.text()


def test_measured_reactivity_is_drawn_against_the_expected_one(window, monkeypatch):
    hours = np.arange(8.0, 11.0, 1 / 60)
    monkeypatch.setattr(jg.InputsPlot, "_measured_koh", lambda self, spec: (hours, np.full(len(hours), 3.0)))
    window.plot_inj.plot()
    (ax,) = [a for a in window.plot_inj.fig.axes if a.get_title().startswith("OH reactivity: measured vs expected")]
    assert [ln.get_label() for ln in ax.lines][:2] == ["measured (LP-LIF)", "expected from the injections"]


def test_nothing_to_plot_is_said_in_the_plot(app):
    w = jg.JuiceWindow(TOML)
    try:
        w.plot_cond.plot()
        assert w.plot_cond.fig.axes[0].texts[0].get_text().startswith("nothing to plot yet")
    finally:
        w.dirty = False; w.close()


# ── run step ─────────────────────────────────────────────────────────────────

def test_prepare_needs_an_experiment_without_errors(window):
    window.page_cond.date.setText("")
    assert window.page_run.prepare() == ""
    assert "Cannot prepare" in window.page_run.log.toPlainText() and window.page_run.bar.format() == "fix the errors first"


@pytest.fixture
def prepared(window):
    js = window.page_run.prepare()
    assert js, window.page_run.log.toPlainText()
    return window, js


def test_prepare_writes_the_input_files_and_the_script(prepared):
    window, js = prepared
    run = window.page_run
    out = os.path.dirname(js)
    assert out == window.out_dir() and os.path.isfile(js)
    for name in ("uv_ethanol.enz", "uv_ethanol.easy", "uv_ethanol.experiment.json", "uv_ethanol_run.py"):
        assert os.path.isfile(os.path.join(out, name)), name
    assert run.files.text().startswith("Written: uv_ethanol.enz, uv_ethanol.easy, uv_ethanol.experiment.json, uv_ethanol.json")
    assert run.run_btn.isEnabled() and not run.view_btn.isEnabled() and run.bar.format() == "inputs ready"
    assert window.steps.item(jg.RUN_STEP).text().endswith("✔")
    files = items(window.review.files)
    assert "uv_ethanol.easy   — model input shared with FACSIMILE (EASY)" in files
    assert any(f.startswith("chemistry") and f.endswith("chemistry include files") for f in files)


def test_choosing_another_solver_asks_for_prepare_again(prepared):
    window, _js = prepared
    run = window.page_run
    other = next(i for i in range(run.solver.count()) if i != run.solver.currentIndex())
    run.solver.setCurrentIndex(other)
    assert run.json_path == "" and not run.run_btn.isEnabled() and run.solver_note.text()


@pytest.mark.slow
def test_a_run_from_the_window_and_the_comparison_with_a_reference(prepared, tmp_path):
    window, js = prepared
    run = window.page_run
    run.run()
    assert run._proc.waitForFinished(300_000), run.log.toPlainText()[-2000:]
    QApplication.processEvents()
    assert run.bar.format().startswith("done") and run.view_btn.isEnabled() and run.run_btn.isEnabled()
    info = json.load(open(js[:-5] + ".run.json", encoding="utf-8"))
    assert info["exit"] == 0 and info["cs_bytes"] > 0 and info["elapsed"] > 0
    assert "$ python -m juice.core.scenario_runner uv_ethanol.json" in run.log.toPlainText()
    # an imported FACSIMILE experiment keeps its FACSIMILE result as a reference: here the run itself stands in
    ref = tmp_path / "reference" / "uv_ethanol.cs"
    ref.parent.mkdir()
    shutil.copy2(run.cs_path(), ref)
    window._base.reference_cs = str(ref)
    run.refresh()
    assert run.cmp_btn.isEnabled()
    run.compare()
    log = run.log.toPlainText()
    assert "JUICE vs FACSIMILE (uv_ethanol.cs)" in log and "✔ the two runs agree" in log
    assert run.bar.format() == "agree with FACSIMILE"
    assert os.path.isfile(os.path.join(os.path.dirname(js), "uv_ethanol_difference_JUICE_Facsimile.log"))


def test_results_studies_and_the_folder_are_opened_as_programs_of_their_own(prepared, started, monkeypatch):
    window, js = prepared
    run = window.page_run
    open(run.cs_path(), "w").close()
    run.open_results()
    assert started[-1][1][:2] == ["-m", "juice.gui.results_viewer"] and started[-1][1][2] == run.cs_path()
    run._study_from_here()
    spec_path = os.path.join(window.out_dir(), "uv_ethanol.experiment.json")
    assert started[-1][1][:2] == ["-m", "juice.gui.study_window"]
    assert started[-1][1][-4:] == ["--mode", "sweep", "--experiment", spec_path]
    window._open_results_file(run.cs_path())
    assert started[-1][1][-2:] == ["--toml", TOML]
    window._open_folder()
    assert started[-1] == ("open", [window.out_dir()])
    run.stop()                                                          # nothing running: nothing happens


def test_sweeping_before_the_files_exist_asks_first(window, started, monkeypatch):
    monkeypatch.setattr(QMessageBox, "question", staticmethod(lambda *a, **k: QMessageBox.StandardButton.No))
    window.page_run._study_from_here()
    assert started == [] and not os.path.isdir(window.out_dir())


# ── file menu ────────────────────────────────────────────────────────────────

def test_save_open_and_new(window, tmp_path):
    window._save()
    path = os.path.join(window.out_dir(), "uv_ethanol.experiment.json")
    assert os.path.isfile(path)
    other = jg.JuiceWindow(TOML)
    try:
        other._open(path)
        spec = other.collect_spec()
        assert spec.name == "uv_ethanol" and spec.precursors == ["C2H5OH"] and spec.injections[0].total_ppb == "200"
        assert other.page_start.out_edit.text() == str(tmp_path)
        other._new()
        assert jg.is_blank(other.collect_spec())
    finally:
        other.dirty = False; other.close()


def test_the_settings_dialog_lists_where_juice_reads_and_writes(window, monkeypatch):
    shown = []

    def capture(dialog):
        shown.append([lb.text() for lb in dialog.findChildren(QLabel)])
        return 0
    monkeypatch.setattr(QDialog, "exec", capture)
    window._settings()
    (labels,) = shown
    assert TOML in labels and window.include_dir in labels
    assert any(t.startswith("BDF (analytic Jacobian") for t in labels)
    assert any("rtol" in t and "(the Run step can pick the other profile per run)" in t for t in labels)


def test_back_and_next(window):
    window.steps.setCurrentRow(0)
    window._step(+1)
    assert window.steps.currentRow() == 1 and window.back_btn.isEnabled()
    for _ in range(5):
        window._step(+1)
    assert window.steps.currentRow() == jg.RUN_STEP and not window.next_btn.isEnabled()
    assert window.page_run.folder.text() == window.out_dir()
    window._step(-1)
    assert window.steps.currentRow() == jg.RUN_STEP - 1
