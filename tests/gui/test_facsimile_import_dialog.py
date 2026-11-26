"""The FACSIMILE import dialog, opened without a screen, on the golden scenes: choosing the file finds the ENZ, the
EASY version and the FACSIMILE results; reading shows what the import made of them; accepting copies the ENZ into
the experiment folder.  The dialog's settings go to a temporary file, not to the user's registry.  Also the import's
self-check (fidelity) and the helpers that find the files."""
import os
import shutil

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import pytest
from PyQt6.QtCore import QSettings
from PyQt6.QtWidgets import QApplication, QDialogButtonBox, QMessageBox

from juice.importers import facsimile_import as fi

HERE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
IDL = os.path.join(HERE, "validation", "IDL Input")


def _scene(name):
    return os.path.join(IDL, name, name + ".easy"), os.path.join(IDL, name, name + ".enz")


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication([])


@pytest.fixture
def boxes(monkeypatch):
    seen = []
    for kind in ("warning", "information", "critical"):
        monkeypatch.setattr(QMessageBox, kind, staticmethod(
            lambda *a, _k=kind, **kw: seen.append((_k, a[2] if len(a) > 2 else "")) or QMessageBox.StandardButton.Ok))
    return seen


@pytest.fixture
def dialog(app, tmp_path, monkeypatch, boxes):
    ini = str(tmp_path / "settings.ini")
    monkeypatch.setattr(fi, "QSettings", lambda *_a: QSettings(ini, QSettings.Format.IniFormat))
    project = tmp_path / "project"
    project.mkdir()

    def make(path=""):
        d = fi.ImportDialog(path=path, project_root=str(project), include_dir=os.path.join(IDL, "include files"))
        d.boxes = boxes
        return d
    return make


def _ok(d):
    return d.buttons.button(QDialogButtonBox.StandardButton.Ok)


def test_choosing_a_scene_finds_its_files(dialog):
    easy, enz = _scene("scene_09")
    d = dialog(easy)
    assert os.path.samefile(d.enz.text(), enz) and d.name.text() == "scene_09"
    found = d.found.text()
    assert "EASY 6.00" in found and "6.00 rule" in found and "FACSIMILE results (.cs) found" in found
    assert not _ok(d).isEnabled()


def test_reading_and_accepting_copies_the_enz(dialog, tmp_path):
    easy, enz = _scene("scene_09")
    d = dialog(easy)
    d.name.setText("my scene 09!")
    d.read()
    assert _ok(d).isEnabled() and d.result_spec is not None
    report = d.report.toPlainText()
    assert report.startswith("my_scene_09") or "2024-08-05 04:00" in report
    assert "Roof open: 09:33–18:37" in report and "Injections (" in report and "O3" in report
    assert "include file(s) not found" not in report
    d._accept()
    out = tmp_path / "project" / "my_scene_09"
    assert d.result_spec.name == "my_scene_09"
    assert d.result_spec.conditions.input_file == str(out / "my_scene_09.measured.enz")
    assert (out / "my_scene_09.measured.enz").read_bytes() == open(enz, "rb").read()
    assert os.path.isfile(enz)                                   # the original stays where it was


def test_missing_files_are_said(dialog, tmp_path):
    d = dialog()
    d.read()
    assert d.boxes[-1] == ("warning", "Choose the .easy of a FACSIMILE experiment first.")
    easy, _enz = _scene("scene_05")
    folder = tmp_path / "alone"
    folder.mkdir()
    shutil.copy(easy, folder / "scene_05.easy")                  # no ENZ, no .fac, no .cs, no include files
    d._file_changed(str(folder / "scene_05.easy"))
    assert d.enz.text() == ""
    assert "EASY version unknown" in d.found.text() and "no FACSIMILE results" in d.found.text()
    d.read()
    assert d.boxes[-1] == ("warning", "The ENZ file was not found: choose it.")
    shutil.copy(_scene("scene_05")[1], folder / "scene_05.enz")
    d.enz.setText(str(folder / "scene_05.enz"))
    d.include_dir = str(folder)                                  # JUICE's include folder not there either
    d.read()
    report = d.report.toPlainText()
    assert "✖ include file(s) not found in any folder:" in report


def test_only_the_elog_is_left(dialog, tmp_path):
    """An .elog without its .easy: the master file and its includes are rebuilt next to it."""
    folder = tmp_path / "elog_only"
    folder.mkdir()
    for ext in (".elog", ".enz", ".fac"):
        shutil.copy(os.path.join(IDL, "scene_04", "scene_04" + ext), folder / ("scene_04" + ext))
    d = dialog(str(folder / "scene_04.elog"))
    assert "rebuilt with its" in d.found.text()
    assert d.easy_path == str(folder / "scene_04.from_elog" / "scene_04.easy")
    d.read()
    assert _ok(d).isEnabled() and "include file(s) not found" not in d.report.toPlainText()


def test_include_folders_are_remembered(dialog, tmp_path, monkeypatch):
    from PyQt6.QtWidgets import QFileDialog
    d = dialog()
    monkeypatch.setattr(QFileDialog, "getExistingDirectory", staticmethod(lambda *a, **k: str(tmp_path)))
    d._add_dir()
    d._add_dir()                                                 # the same folder twice: once
    assert d.dirs.count() == 1
    assert dialog().dirs.count() == 1                            # a new dialog reads it back
    d.dirs.item(0).setSelected(True)
    d._remove_dir()
    assert dialog().dirs.count() == 0


# ── the self-check and the helpers ───────────────────────────────────────────────────────────────────────────────

def test_fidelity_is_empty_for_a_faithful_import_and_names_a_change():
    easy, enz = _scene("scene_09")
    spec, _notes = fi.import_facsimile(easy, enz)
    assert fi.fidelity(spec, easy, enz) == []
    o3 = next(i for i in spec.injections if i.name == "O3")
    o3.start = "10:36"                                           # an hour later than the ENZ pulse
    notes = fi.fidelity(spec, easy, enz)
    assert any(n.startswith("ENZ: column") and "differs" in n for n in notes), notes
    spec.conditions.end = "20:59"
    assert any("time grid differs" in n for n in fi.fidelity(spec, easy, enz))
    assert any("could not be checked" in n for n in fi.fidelity(spec, easy, os.path.join(IDL, "missing.enz")))


def test_finding_the_enz_and_names(tmp_path):
    (tmp_path / "a.easy").write_text("; FILES[ENZ] = old.enz\nFILES[ENZ] = run.enz ; the day\n")
    assert fi.enz_name_in(str(tmp_path / "a.easy")) == "run.enz"
    assert fi.enz_name_in(str(tmp_path / "missing.easy")) == ""
    assert fi.find_enz(str(tmp_path), "a", "run.enz") == ""                # named, but not there; no other ENZ
    (tmp_path / "other.enz").write_text("")
    assert fi.find_enz(str(tmp_path), "a", "run.enz") == str(tmp_path / "other.enz")   # the only one
    (tmp_path / "a.enz").write_text("")
    assert fi.find_enz(str(tmp_path), "a", "run.enz") == str(tmp_path / "a.enz")
    (tmp_path / "run.enz").write_text("")
    assert fi.find_enz(str(tmp_path), "a", "run.enz") == str(tmp_path / "run.enz")
    assert fi.find_enz(str(tmp_path / "nope"), "a", "") == ""
    assert fi.safe_name("Day 3 (rerun)") == "Day_3_rerun" and fi.safe_name("??") == "experiment"


def test_elog_rebuilt_in_a_temporary_folder_when_its_own_cannot_be_written(tmp_path, monkeypatch):
    from juice.core import elog_generator as eg
    real = eg.split_elog
    tried = []

    def split(elog, d, stem):
        tried.append(d)
        if len(tried) == 1:
            raise PermissionError(d)
        return real(elog, d, stem)
    monkeypatch.setattr(eg, "split_elog", split)
    easy, incs = fi.easy_from_elog(os.path.join(IDL, "scene_04", "scene_04.elog"))
    assert tried[0].endswith("scene_04.from_elog") and os.path.basename(os.path.dirname(easy)).startswith("juice_scene_04_")
    assert os.path.isfile(easy) and incs
    monkeypatch.setattr(eg, "split_elog", lambda *a: (_ for _ in ()).throw(PermissionError("no")))
    with pytest.raises(PermissionError, match="no folder to write"):
        fi.easy_from_elog(os.path.join(IDL, "scene_04", "scene_04.elog"))
