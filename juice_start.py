"""
juice_start.py — the JUICE start menu: choose what you need, then the matching window opens.

    python juice_start.py                    # the menu (what the JUICE desktop icon opens)
    python juice_start.py <file> [-t toml]   # a file dropped on the icon opens straight away:
                                             #   .cs → results viewer, .study.json → study window,
                                             #   .experiment.json → experiment window

Model an experiment:  new experiment · open a saved experiment · new study (sweep or several experiments)
Look at results:      one experiment's results (.cs) · a study's results (.study.json)

Every choice starts its own, independent window (its own process) and closes the menu, so clicking
the icon again always gives a fresh menu, whatever is already open.  `choose_results()` is the same
"one experiment or a study?" question for the results viewer when it is started without a file.
"""
from __future__ import annotations

import os
import sys
from typing import Optional, Tuple

from PyQt6.QtCore import Qt, QProcess, pyqtSignal
from PyQt6.QtGui import QPixmap
from PyQt6.QtWidgets import (QApplication, QDialog, QFileDialog, QFrame, QGroupBox, QHBoxLayout, QLabel, QMessageBox,
                             QPushButton, QVBoxLayout, QWidget)

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = {"experiment": "juice_gui.py", "study": "study_window.py", "results": "results_viewer.py"}


def _toml_args(toml: str) -> list:
    return ["--toml", toml] if toml else []


def launch(kind: str, args: list = (), toml: str = "") -> bool:
    """Start one JUICE window in its own process (independent of this one)."""
    ok, _pid = QProcess.startDetached(sys.executable, [os.path.join(HERE, SCRIPTS[kind])] + list(args) + _toml_args(toml), HERE)
    return bool(ok)


def open_file(path: str, toml: str = "") -> bool:
    """The window that belongs to a file: .cs → viewer, .study.json → study, .experiment.json → experiment."""
    p = path.lower()
    if p.endswith(".study.json"):
        return launch("study", [path], toml)
    if p.endswith(".experiment.json"):
        return launch("experiment", [path], toml)
    if p.endswith(".cs") or p.endswith("_cs.csv"):
        return launch("results", [path], toml)
    return False


def _start_dir(toml: str) -> str:
    try:
        from scenario_prep import load_paths
        paths = load_paths(toml)
        return next((d for d in (paths.get("output_root"), paths.get("easy_root")) if d and os.path.isdir(d)), HERE)
    except Exception:
        return HERE


def pick_result(parent, toml: str) -> str:
    p, _ = QFileDialog.getOpenFileName(parent, "Results of one experiment", _start_dir(toml),
                                       "JUICE / FACSIMILE results (*.cs *_cs.csv);;All files (*)")
    return p


def pick_study(parent, toml: str) -> str:
    p, _ = QFileDialog.getOpenFileName(parent, "Results of a study", _start_dir(toml), "Study (*.study.json)")
    return p


def pick_experiment(parent, toml: str) -> str:
    p, _ = QFileDialog.getOpenFileName(parent, "Open a saved experiment", _start_dir(toml), "Experiment (*.experiment.json)")
    return p


class _Choice(QFrame):
    """A large clickable card: a title line and a wrapped description (a frame, so it grows with the text)."""
    clicked = pyqtSignal()

    def __init__(self, title: str, text: str, primary: bool = False):
        super().__init__()
        self.setCursor(Qt.CursorShape.PointingHandCursor)
        self.setObjectName("choicePrimary" if primary else "choice")
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        bg, hover, fg, sub = (("#6b3e85", "#4a3e85", "#ffffff", "#f2ecf7") if primary else
                              ("#f3f3f3", "#efe6f7", "#4a3e85", "#555555"))
        self.setStyleSheet(f"QFrame#{self.objectName()} {{ background: {bg}; border: 1px solid "
                           f"{'#6b3e85' if primary else '#c8c8c8'}; border-radius: 6px; }}"
                           f"QFrame#{self.objectName()}:hover {{ background: {hover}; border-color: #6b3e85; }}")
        lay = QVBoxLayout(self); lay.setContentsMargins(14, 10, 14, 10); lay.setSpacing(3)
        t = QLabel(title); d = QLabel(text); d.setWordWrap(True)
        t.setStyleSheet(f"font-size: 14px; font-weight: bold; color: {fg}; background: transparent;")
        d.setStyleSheet(f"color: {sub}; background: transparent;")
        for w in (t, d):
            w.setAttribute(Qt.WidgetAttribute.WA_TransparentForMouseEvents)
            lay.addWidget(w)
        self.setToolTip(text)

    def mouseReleaseEvent(self, ev):
        if ev.button() == Qt.MouseButton.LeftButton and self.rect().contains(ev.position().toPoint()):
            self.clicked.emit()
        super().mouseReleaseEvent(ev)

    def keyPressEvent(self, ev):
        if ev.key() in (Qt.Key.Key_Return, Qt.Key.Key_Enter, Qt.Key.Key_Space):
            self.clicked.emit()
        else:
            super().keyPressEvent(ev)


class ResultsChoice(QDialog):
    """'One experiment or a study?' — the results viewer's first question when it has no file."""

    def __init__(self, toml: str = "", parent=None):
        super().__init__(parent)
        self.toml, self.result_kind, self.path = toml, "", ""
        self.setWindowTitle("JUICE — Results"); self.setMinimumWidth(520)
        v = QVBoxLayout(self)
        v.addWidget(QLabel("Which results do you want to look at?"))
        b1 = _Choice("One experiment", "The time series, budgets and reactivity of one run (a .cs file, "
                     "JUICE or FACSIMILE).", primary=True)
        b2 = _Choice("A study", "A sweep or a set of experiments run with one change: its results table, "
                     "isopleths and every run (a .study.json file).")
        b1.clicked.connect(lambda: self._pick("results")); b2.clicked.connect(lambda: self._pick("study"))
        v.addWidget(b1); v.addWidget(b2)
        cancel = QPushButton("Cancel"); cancel.clicked.connect(self.reject)
        h = QHBoxLayout(); h.addStretch(); h.addWidget(cancel); v.addLayout(h)

    def _pick(self, kind: str):
        p = pick_result(self, self.toml) if kind == "results" else pick_study(self, self.toml)
        if p:
            self.result_kind, self.path = kind, p
            self.accept()


def choose_results(parent=None, toml: str = "") -> Optional[Tuple[str, str]]:
    """('results', <.cs>) or ('study', <.study.json>), or None when cancelled."""
    d = ResultsChoice(toml, parent)
    return (d.result_kind, d.path) if d.exec() == QDialog.DialogCode.Accepted else None


class StartWindow(QWidget):
    def __init__(self, toml: str = ""):
        super().__init__()
        import juice_theme as theme
        self.toml = toml
        self.setWindowTitle("JUICE"); self.setWindowIcon(theme.app_icon()); self.setMinimumWidth(820)
        v = QVBoxLayout(self); v.setContentsMargins(22, 18, 22, 16); v.setSpacing(12)

        head = QHBoxLayout()
        if os.path.isfile(theme.LOGO):
            logo = QLabel(); logo.setPixmap(QPixmap(theme.LOGO).scaledToHeight(64, Qt.TransformationMode.SmoothTransformation))
            head.addWidget(logo)
        tl = QVBoxLayout()
        t = QLabel("JUICE"); t.setObjectName("title"); tl.addWidget(t)
        s = QLabel("Jülich Utility for Integrating Chemical Equations — what do you want to do?"); s.setObjectName("subtitle")
        tl.addWidget(s); head.addLayout(tl); head.addStretch(); v.addLayout(head)

        cols = QHBoxLayout(); cols.setSpacing(14)
        g = QGroupBox("Model an experiment"); gl = QVBoxLayout(g); gl.setSpacing(8)
        self._add(gl, "New experiment", "Start from scratch: chamber, conditions, chemistry, injections, then run the "
                  "solver. Templates, the SAPHIR archive and wiki logs are offered on the first step.",
                  lambda: launch("experiment", [], self.toml), primary=True)
        self._add(gl, "Open a saved experiment", "Change or rerun an experiment described earlier (.experiment.json).",
                  self._open_experiment)
        self._add(gl, "New study", "Many runs at once: sweep one experiment's injections or conditions over a grid "
                  "(isopleths, sensitivities), or run several experiments with one change.",
                  lambda: launch("study", [], self.toml))
        gl.addStretch(); cols.addWidget(g, 1)

        g = QGroupBox("Look at results"); gl = QVBoxLayout(g); gl.setSpacing(8)
        self._add(gl, "Results of one experiment", "A single run in the results viewer: time series, budgets, "
                  "reactivity, comparison with measurements (.cs, JUICE or FACSIMILE).", self._open_result)
        self._add(gl, "Results of a study", "A finished or partly finished study: results table, isopleths, "
                  "isopleth in time, and each run (.study.json). Nothing is run again.", self._open_study)
        gl.addStretch(); cols.addWidget(g, 1)
        v.addLayout(cols)

        foot = QHBoxLayout()
        n = QLabel("Each choice opens its own window. Click the JUICE icon again for another one at any time.")
        n.setObjectName("subtitle"); foot.addWidget(n); foot.addStretch()
        q = QPushButton("Close"); q.clicked.connect(self.close); foot.addWidget(q)
        v.addLayout(foot)

    def _add(self, layout, title, text, slot, primary=False):
        b = _Choice(title, text, primary); b.clicked.connect(lambda: self._go(slot)); layout.addWidget(b)

    def _go(self, slot):
        """Run a choice; close the menu once a window has been started."""
        started = slot()
        if started:
            self.close()
        elif started is False:
            QMessageBox.warning(self, "JUICE", "The window could not be started (see the console).")

    def _open_experiment(self):
        p = pick_experiment(self, self.toml)
        return launch("experiment", [p], self.toml) if p else None

    def _open_result(self):
        p = pick_result(self, self.toml)
        return launch("results", [p], self.toml) if p else None

    def _open_study(self):
        p = pick_study(self, self.toml)
        return launch("study", [p], self.toml) if p else None


def main(argv=None) -> int:
    argv = list(sys.argv if argv is None else argv)
    sys.path.insert(0, HERE)
    import juice_theme as theme
    from scenario_prep import pick_toml
    toml = pick_toml(argv)
    skip = {i + 1 for i, a in enumerate(argv) if a in ("-t", "--toml")}
    files = [a for i, a in enumerate(argv[1:], 1) if i not in skip and not a.startswith("-") and os.path.isfile(a)]
    theme.taskbar_identity()
    app = QApplication(argv); theme.apply_theme(app)
    if files:                                     # dropped on the icon: straight to the right window
        if open_file(files[0], toml):
            return 0
        QMessageBox.warning(None, "JUICE", f"JUICE does not know what to open {os.path.basename(files[0])} with.\n"
                                           "It opens .cs results, .study.json studies and .experiment.json experiments.")
    w = StartWindow(toml); w.show(); theme.force_taskbar_icon(w)
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
