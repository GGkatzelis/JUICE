"""The expert editor (easy_builder.EasyBuilderDialog; Review step -> "Open expert editor ..."), opened without a screen.

Every golden scene is loaded into the dialog's widgets and read back: the EASY written from what the widgets hold must
be statement-equivalent to the scene, so no field of the form is lost on its way through a table or a text box.  Then
the actions: an ENZ chosen fills FILES[ENZ] and the output path, its Q columns become injection rows, and writing puts
the .easy and a copy of the ENZ in the output folder."""
import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import pytest
from PyQt6.QtWidgets import QApplication, QFileDialog, QMessageBox

from juice.core import easy_builder as eb
from juice.core.easy_scene import statements_equivalent

HERE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
IDL = os.path.join(HERE, "validation", "IDL Input")
INCLUDES = os.path.join(IDL, "include files")
SCENES = sorted(d for d in os.listdir(IDL) if d.startswith("scene_") and os.path.isfile(os.path.join(IDL, d, d + ".easy")))


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication([])


@pytest.fixture
def boxes(monkeypatch):
    seen = []
    for kind in ("warning", "information", "critical"):
        monkeypatch.setattr(QMessageBox, kind, staticmethod(
            lambda *a, _k=kind, **kw: seen.append((_k, a[1] if len(a) > 1 else "")) or QMessageBox.StandardButton.Ok))
    return seen


@pytest.fixture
def dlg(app, boxes):
    d = eb.EasyBuilderDialog(include_dir=INCLUDES)
    d.boxes = boxes
    return d


@pytest.mark.parametrize("scene", SCENES)
def test_a_scene_survives_the_dialogs_widgets(dlg, scene, tmp_path):
    src = os.path.join(IDL, scene, scene + ".easy")
    dlg.set_form(eb.load_form(src))
    out = tmp_path / (scene + ".easy")
    out.write_text(eb.render_form(dlg.get_form()), encoding="latin-1", errors="replace")   # as _write does
    ok, diffs = statements_equivalent(src, str(out))
    assert ok, diffs[:6]


def test_choosing_an_enz_fills_the_file_names_and_its_q_columns(dlg):
    enz = os.path.join(IDL, "scene_09", "scene_09.enz")
    dlg.enz_edit.setText(enz)
    dlg._enz_changed()
    assert dlg.files_enz.text() == "scene_09.enz"
    assert dlg.out_edit.text() == os.path.join(IDL, "scene_09", "scene_09.easy")
    q_cols = [c for c in dlg.enz_columns if c.upper().startswith("Q")]
    assert q_cols
    before = dlg.inj_table.rowCount()
    dlg._fill_from_enz()
    rows = eb.EasyBuilderDialog._get_rows(dlg.inj_table)
    assert dlg.inj_table.rowCount() == before + len(q_cols)
    assert {r[1] for r in rows[before:]} == set(q_cols) and all(r[4] == "from ENZ" for r in rows[before:])
    dlg._fill_from_enz()                                          # every column is used now: nothing added
    assert dlg.inj_table.rowCount() == before + len(q_cols)


def test_loading_previewing_and_writing(dlg, tmp_path, monkeypatch):
    src = os.path.join(IDL, "scene_05", "scene_05.easy")
    monkeypatch.setattr(QFileDialog, "getOpenFileName", staticmethod(lambda *a, **k: (src, "")))
    dlg._load_easy()
    assert dlg.enz_edit.text() == os.path.join(IDL, "scene_05", "scene_05.enz")       # FILES[ENZ] found next to it
    shown = tmp_path / "preview.easy"
    shown.write_text(dlg.preview.toPlainText(), encoding="latin-1", errors="replace")
    ok, diffs = statements_equivalent(src, str(shown))                               # comments may be re-spaced
    assert ok, diffs[:6]
    out = tmp_path / "written" / "copy.easy"
    dlg.out_edit.setText(str(out))
    written = []
    dlg.easy_written.connect(written.append)
    dlg._write(open_after=False)
    assert dlg.boxes[-1][0] == "information"
    assert (tmp_path / "written" / "scene_05.enz").is_file()                           # the ENZ goes with it
    ok, diffs = statements_equivalent(src, str(out))
    assert ok, diffs[:6]
    dlg._write(open_after=True)
    assert written == [str(out)]


def test_writing_without_a_path_is_refused(dlg):
    dlg.out_edit.setText("")
    dlg._write(open_after=False)
    assert dlg.boxes == [("warning", "No output path")]


def test_a_file_that_cannot_be_read_is_said(dlg, tmp_path, monkeypatch):
    bad = tmp_path / "bad.easy"
    bad.write_bytes(b"\xff\xfe\x00 not an easy file {{{")
    monkeypatch.setattr(eb, "load_form", lambda p: (_ for _ in ()).throw(ValueError("cannot parse")))
    monkeypatch.setattr(QFileDialog, "getOpenFileName", staticmethod(lambda *a, **k: (str(bad), "")))
    dlg._load_easy()
    assert dlg.boxes == [("critical", "Load failed")]
