"""
juice_gui.py — JUICE, step by step.

One window, six steps down the left, the active step in the middle, a live review
on the right, and the solver run as the last step.  No pop-ups on the main path.

    1 Experiment        name, folder, template, wiki import, date/time, conditions, roof
    2 Chemistry         precursors (MCM subset) or an existing mechanism
    3 Injections        actions, amounts, mixtures
    4 Chamber           wall deposition, sources, background reactivity, dilution
    5 Initial values    concentrations at the start
    6 Run               preset, prepare the input files, run the solver, open the results

The only thing a user edits is the experiment description (<name>.experiment.json).
The ENZ, EASY and scenario JSON files are generated from it in the project folder
(IDL/FACSIMILE naming: QNO, QO3, QH2O, Q_<name> columns; q<NAME> switches) and are
never hand-edited here; FR_GUI.py remains the expert tool for that.

    python juice_gui.py                 # paths from GG_GUI.toml if present, else FR_GUI.toml
    python juice_gui.py -t FR_GUI.toml  # another TOML

The step pages are the ones of experiment_wizard.py (shared through ExperimentHost),
so the QWizard behind FR_GUI's "New experiment" button and this window stay identical.
"""
from __future__ import annotations

import json
import os
import re
import sys
import time
import tomllib
from typing import Dict, List, Optional

from PyQt6.QtCore import Qt, QProcess
from PyQt6.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QGridLayout, QLabel,
                             QListWidget, QListWidgetItem, QStackedWidget, QPushButton, QGroupBox, QScrollArea,
                             QComboBox, QLineEdit, QCheckBox, QProgressBar, QTextEdit, QPlainTextEdit, QMessageBox,
                             QFileDialog, QSplitter, QSizePolicy)

from experiment_spec import load_spec, validate, mechanism_species
from experiment_wizard import ExperimentHost, _info
from elog_generator import generate_both
from idl_parser import parse_idl_easy

HERE = os.path.dirname(os.path.abspath(__file__))
STEPS = ["1  Experiment", "2  Chemistry", "3  Injections", "4  Chamber", "5  Initial values", "6  Run"]

# Solver presets → keys written into the scenario JSON (everything else comes from the TOML)
PRESETS = {
    "Standard (validated against FACSIMILE)": dict(perstep_snap=True, perstep_warm=True, interpolation="piecewise",
                                                   rtol=1e-4, first_step=1e-6, export_kv_pd=True),
    "Quick look (no rate/budget files)":      dict(perstep_snap=True, perstep_warm=True, interpolation="piecewise",
                                                   rtol=1e-4, first_step=1e-6, export_kv_pd=False),
    "Custom":                                 None,
}


# ─────────────────────────────────────────────────────────────────────────────
def load_paths(toml_path: str) -> Dict[str, str]:
    """[paths] of the TOML (+ .local.toml overlay); defaults relative to this file."""
    base = {"easy_root": os.path.join(HERE, "MCM_tested_scenarios", "IDL Input"),
            "include_dir": os.path.join(HERE, "MCM_tested_scenarios", "IDL Input", "include files"),
            "output_root": os.path.join(HERE, "MCM_tested_scenarios", "JUICE archive", "JUICE"),
            "last_easy_dir": ""}
    for p in (toml_path, toml_path[:-5] + ".local.toml" if toml_path.endswith(".toml") else ""):
        if p and os.path.isfile(p):
            try:
                with open(p, "rb") as fh:
                    base.update({k: str(v) for k, v in tomllib.load(fh).get("paths", {}).items()})
            except Exception as e:
                print(f"[juice_gui] could not read {p}: {e}")
    return base


def pick_toml(argv: List[str]) -> str:
    for i, a in enumerate(argv):
        if a in ("-t", "--toml") and i + 1 < len(argv):
            p = argv[i + 1]
            return p if os.path.isabs(p) else os.path.join(HERE, p)
    for name in ("GG_GUI.toml", "FR_GUI.toml"):
        if os.path.isfile(os.path.join(HERE, name)):
            return os.path.join(HERE, name)
    return ""


def prepare_scenario(res: dict, include_dir: str, solver: dict) -> str:
    """EASY + ENZ (from write_all) → canonical elog → scenario JSON next to them.
    Returns the JSON path.  The JSON stem = experiment name, so the runner writes <name>.cs."""
    out = os.path.dirname(res["easy"])
    raw, canon = generate_both(res["easy"], search_dirs=[out, os.path.join(out, "chemistry"), include_dir], output_dir=out)
    cfg = parse_idl_easy(canon, enz_path=res["enz"], chemistry_folder=out,
                         include_fallbacks=[os.path.join(out, "chemistry"), include_dir])
    cfg.update({"folder": out, "enz_path": res["enz"], "elog_path": canon,
                "ode_backend": "numba", "solver_method": "BDF", "jac_mode": "analytic"})
    cfg.update(solver)
    name = os.path.splitext(os.path.basename(res["easy"]))[0]
    path = os.path.join(out, name + ".json")
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(cfg, fh, indent=1, default=str)
    return path


def _boxed(page: QWidget, title: str = "") -> QWidget:
    """A wizard page shown as a plain section: title as group box, page inside."""
    box = QGroupBox(title or page.title())
    lay = QVBoxLayout(box); lay.setContentsMargins(8, 6, 8, 6)
    sub = page.subTitle()
    if sub:
        lay.addWidget(_info(sub))
    lay.addWidget(page)
    return box


def _scroll(*widgets: QWidget) -> QScrollArea:
    inner = QWidget(); lay = QVBoxLayout(inner); lay.setContentsMargins(0, 0, 0, 0)
    for w in widgets:
        lay.addWidget(w)
    lay.addStretch()
    sa = QScrollArea(); sa.setWidgetResizable(True); sa.setWidget(inner); sa.setFrameShape(QScrollArea.Shape.NoFrame)
    return sa


# ─────────────────────────────────────────────────────────────────────────────
class RunPage(QWidget):
    """Step 6: preset → prepare the input files → run scenario_runner → open the results."""

    def __init__(self, win: "JuiceWindow"):
        super().__init__()
        self.win = win
        self.json_path = ""
        self._proc: Optional[QProcess] = None
        self._t0 = 0.0
        lay = QVBoxLayout(self)
        g = QGridLayout(); r = 0
        g.addWidget(QLabel("Project folder"), r, 0); self.folder = QLabel(""); self.folder.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        g.addWidget(self.folder, r, 1, 1, 3); r += 1
        g.addWidget(QLabel("Solver preset"), r, 0); self.preset = QComboBox(); self.preset.addItems(list(PRESETS))
        self.preset.currentIndexChanged.connect(self._preset_changed); g.addWidget(self.preset, r, 1, 1, 3); r += 1
        self.custom = QWidget(); cg = QGridLayout(self.custom); cg.setContentsMargins(0, 0, 0, 0)
        cg.addWidget(QLabel("rtol"), 0, 0); self.rtol = QLineEdit("1e-4"); cg.addWidget(self.rtol, 0, 1)
        cg.addWidget(QLabel("first_step"), 0, 2); self.first_step = QLineEdit("1e-6"); cg.addWidget(self.first_step, 0, 3)
        self.warm = QCheckBox("Warm restart across ENZ rows"); self.warm.setChecked(True); cg.addWidget(self.warm, 1, 0, 1, 2)
        self.kvpd = QCheckBox("Write rate (.kv) and budget (.pd) files"); self.kvpd.setChecked(True); cg.addWidget(self.kvpd, 1, 2, 1, 2)
        g.addWidget(self.custom, r, 0, 1, 4); r += 1
        lay.addLayout(g)
        lay.addWidget(_info("Standard = the settings validated against FACSIMILE (piecewise ENZ, warm restart, rtol 1e-4, "
                            "rate and budget files on). Quick look only skips the large rate/budget files. "
                            f"Everything else (atol floor, Newton control, hmax) comes from {os.path.basename(self.win.toml_path) or 'the TOML'}."))
        hb = QHBoxLayout()
        self.prep_btn = QPushButton("1. Prepare input files"); self.prep_btn.clicked.connect(self.prepare); hb.addWidget(self.prep_btn)
        self.run_btn = QPushButton("2. Run simulation"); self.run_btn.setEnabled(False); self.run_btn.clicked.connect(self.run); hb.addWidget(self.run_btn)
        self.stop_btn = QPushButton("Stop"); self.stop_btn.setEnabled(False); self.stop_btn.clicked.connect(self.stop); hb.addWidget(self.stop_btn)
        self.view_btn = QPushButton("3. Open results"); self.view_btn.setEnabled(False); self.view_btn.clicked.connect(self.open_results); hb.addWidget(self.view_btn)
        hb.addStretch(); lay.addLayout(hb)
        self.files = QLabel(""); self.files.setWordWrap(True); lay.addWidget(self.files)
        self.bar = QProgressBar(); self.bar.setRange(0, 100); self.bar.setTextVisible(True); self.bar.setFormat("not started")
        self.bar.setStyleSheet("QProgressBar { text-align: center; } QProgressBar::chunk { background-color: #2ecc40; }")
        lay.addWidget(self.bar)
        self.log = QTextEdit(); self.log.setReadOnly(True)
        self.log.setStyleSheet("font-family: monospace; font-size: 11px; background: #1e1e1e; color: #d4d4d4;")
        lay.addWidget(self.log, stretch=1)
        self._preset_changed()

    def _preset_changed(self):
        self.custom.setVisible(PRESETS[self.preset.currentText()] is None)

    def solver_settings(self) -> dict:
        p = PRESETS[self.preset.currentText()]
        if p is not None:
            return dict(p)
        return dict(perstep_snap=True, perstep_warm=self.warm.isChecked(), interpolation="piecewise",
                    rtol=float(self.rtol.text() or 1e-4), first_step=float(self.first_step.text() or 1e-6),
                    export_kv_pd=self.kvpd.isChecked())

    def refresh(self):
        self.folder.setText(self.win.out_dir())

    # ── actions ──
    def prepare(self) -> str:
        """Write ENZ/EASY/experiment.json and the scenario JSON.  Returns the JSON path ('' on error)."""
        self.refresh()
        spec = self.win.collect_spec()
        known = mechanism_species(self.win.include_dir, spec.mechanism_stem) if spec.mechanism_stem else None
        errors = [m for lvl, m in validate(spec, known) if lvl == "error"]
        if errors:
            self.log.append("Cannot prepare:\n  " + "\n  ".join(errors)); self.bar.setFormat("fix the errors first"); return ""
        try:
            res = self.win.write_all()
            self.json_path = prepare_scenario(res, self.win.include_dir, self.solver_settings())
        except Exception as e:
            self.log.append(f"Prepare failed: {e}"); self.bar.setFormat("prepare failed"); return ""
        out = os.path.dirname(self.json_path)
        self.files.setText("Written: " + ", ".join(os.path.basename(p) for p in (res["enz"], res["easy"], res["spec_path"], self.json_path))
                           + (f"  |  chemistry: {len(res['mechanism'])} include files" if isinstance(res.get("mechanism"), (list, tuple)) else ""))
        for n in res.get("notes", []):
            self.log.append("note: " + n)
        self.log.append(f"Inputs ready in {out}")
        self.bar.setFormat("inputs ready"); self.run_btn.setEnabled(True)
        self.view_btn.setEnabled(os.path.isfile(self.cs_path()))
        self.win.set_step_state(5, "ok")
        return self.json_path

    def cs_path(self) -> str:
        return self.json_path[:-5] + ".cs" if self.json_path else ""

    def run(self):
        if not self.json_path:
            return
        # no --gui: that flag makes the runner spawn the viewer itself; here "Open results" does it
        args = ["-u", os.path.join(HERE, "scenario_runner.py"), self.json_path]
        if self.win.toml_path:
            args += ["--toml", self.win.toml_path]
        self._proc = QProcess(self)
        self._proc.setProcessChannelMode(QProcess.ProcessChannelMode.MergedChannels)
        self._proc.readyReadStandardOutput.connect(self._on_output)
        self._proc.finished.connect(self._on_finished)
        self._t0 = time.monotonic()
        self.bar.setValue(0); self.bar.setFormat("running   0%")
        self.run_btn.setEnabled(False); self.prep_btn.setEnabled(False); self.stop_btn.setEnabled(True); self.view_btn.setEnabled(False)
        self.log.append(f"$ python {' '.join(os.path.basename(a) if os.path.isabs(a) else a for a in args)}")
        self._proc.start(sys.executable, args)

    def stop(self):
        if self._proc and self._proc.state() != QProcess.ProcessState.NotRunning:
            self._proc.kill()

    def _on_output(self):
        raw = bytes(self._proc.readAllStandardOutput()).decode("utf-8", errors="replace")
        self.log.moveCursor(self.log.textCursor().MoveOperation.End); self.log.insertPlainText(raw)
        self.log.moveCursor(self.log.textCursor().MoveOperation.End)
        for m in re.finditer(r"Fraction Done:\s*([\d.]+)", raw):
            pct = min(int(float(m.group(1)) * 100), 99)
            el = int(time.monotonic() - self._t0); mm, ss = divmod(el, 60)
            self.bar.setValue(pct); self.bar.setFormat(f"running   {mm}m {ss:02d}s   {pct}%")

    def _on_finished(self, code, _status):
        el = int(time.monotonic() - self._t0); mm, ss = divmod(el, 60)
        ok = code == 0 and os.path.isfile(self.cs_path())
        self.bar.setValue(100 if ok else self.bar.value())
        self.bar.setFormat(f"done   {mm}m {ss:02d}s" if ok else f"failed (exit {code})   {mm}m {ss:02d}s")
        self.run_btn.setEnabled(True); self.prep_btn.setEnabled(True); self.stop_btn.setEnabled(False); self.view_btn.setEnabled(ok)
        self.win.set_step_state(5, "ok" if ok else "error")

    def open_results(self):
        cs = self.cs_path()
        if not os.path.isfile(cs):
            return
        p = QProcess(self); p.start(sys.executable, [os.path.join(HERE, "results_viewer.py"), cs])
        self._viewer = p


# ─────────────────────────────────────────────────────────────────────────────
class JuiceWindow(QMainWindow, ExperimentHost):
    def __init__(self, toml_path: str = ""):
        super().__init__()
        self.toml_path = toml_path
        paths = load_paths(toml_path)
        self.output_root = paths["output_root"]
        self._init_host(paths["include_dir"], paths["output_root"] if os.path.isdir(paths["output_root"]) else paths["easy_root"])
        self.setWindowTitle("JUICE — Jülich Utility for Integrating Chemical Equations")
        self.resize(1500, 900)

        # ── left: steps ──
        self.steps = QListWidget(); self.steps.setFixedWidth(190)
        self.steps.setStyleSheet("QListWidget { font-size: 14px; } QListWidget::item { padding: 10px 6px; }")
        for s in STEPS:
            self.steps.addItem(QListWidgetItem(s))
        self._state = ["todo"] * len(STEPS)
        self.steps.currentRowChanged.connect(self._go)

        # ── middle: pages ──
        self.stack = QStackedWidget()
        self.stack.addWidget(_scroll(_boxed(self.page_start, "Experiment"), _boxed(self.page_cond, "Where and when")))
        self.stack.addWidget(_scroll(_boxed(self.page_chem, "Chemistry")))
        self.stack.addWidget(_scroll(_boxed(self.page_inj, "Injections")))
        self.stack.addWidget(_scroll(_boxed(self.page_chamber, "What the chamber adds on its own")))
        self.stack.addWidget(_scroll(_boxed(self.page_init, "Initial values")))
        self.page_run = RunPage(self)
        self.stack.addWidget(_scroll(_boxed(self.page_run, "Run") if False else self.page_run))
        nav = QHBoxLayout()
        self.back_btn = QPushButton("◀ Back"); self.back_btn.clicked.connect(lambda: self._step(-1)); nav.addWidget(self.back_btn)
        nav.addStretch()
        self.next_btn = QPushButton("Next ▶"); self.next_btn.clicked.connect(lambda: self._step(+1)); nav.addWidget(self.next_btn)
        mid = QWidget(); ml = QVBoxLayout(mid); ml.setContentsMargins(0, 0, 0, 0); ml.addWidget(self.stack, stretch=1); ml.addLayout(nav)

        # ── right: review (the wizard's review page, always visible) ──
        self.page_review.setMinimumWidth(380)
        self.show_files = QCheckBox("Show the generated EASY text"); self.show_files.setChecked(False)
        self.show_files.toggled.connect(self._toggle_preview)
        right = QWidget(); rl = QVBoxLayout(right); rl.setContentsMargins(4, 4, 4, 4)
        rl.addWidget(QLabel("<b>Review</b> — updates when you change step"))
        rl.addWidget(self.page_review, stretch=1); rl.addWidget(self.show_files)
        self._toggle_preview(False)

        split = QSplitter(); split.addWidget(self.steps); split.addWidget(mid); split.addWidget(right)
        split.setStretchFactor(1, 3); split.setStretchFactor(2, 2)
        self.setCentralWidget(split)

        # ── menu ──
        m = self.menuBar().addMenu("File")
        m.addAction("New experiment (SAPHIR defaults)", self._new)
        m.addAction("Open experiment (.experiment.json) …", self._open)
        m.addAction("Save experiment", self._save)
        m.addSeparator()
        m.addAction("Expert tool (FR_GUI) …", self._expert)
        m.addSeparator()
        m.addAction("Quit", self.close)
        self.statusBar().showMessage(f"paths from {os.path.basename(toml_path) if toml_path else 'defaults'}  |  include files: {self.include_dir}")

        self.load_spec(self._base, keep_name=True)
        self.steps.setCurrentRow(0)

    # ── navigation ──
    def _go(self, row: int):
        if row < 0:
            return
        self.stack.setCurrentIndex(row)
        self.back_btn.setEnabled(row > 0); self.next_btn.setEnabled(row < len(STEPS) - 1)
        if row == len(STEPS) - 1:
            self.page_run.refresh()
        self.page_review.refresh()
        self._update_states()

    def _step(self, d: int):
        row = self.steps.currentRow()
        if d > 0:
            page = (self.page_cond if row == 0 else self.page_chem if row == 1 else None)
            if row == 0 and not self.page_start.validatePage():
                return
            if page is not None and not page.validatePage():
                return
        self.steps.setCurrentRow(max(0, min(len(STEPS) - 1, row + d)))

    def set_step_state(self, i: int, state: str):
        self._state[i] = state
        self._paint_states()

    def _update_states(self):
        """Green/orange/red marks from the spec validation; the current step is at least 'seen'."""
        try:
            spec = self.collect_spec()
            msgs = validate(spec, None)
        except Exception:
            msgs = [("error", "spec")]
        has_err = any(l == "error" for l, _ in msgs)
        txt = " ".join(m for _, m in msgs).lower()
        self._state[0] = "error" if any(k in txt for k in ("date", "time window")) else ("warn" if "j values" in txt or "roof" in txt else "ok")
        self._state[1] = "error" if "no chemistry" in txt else "ok"
        self._state[2] = "warn" if "injection" in txt else "ok"
        self._state[3] = "ok"; self._state[4] = "ok"
        if self._state[5] == "todo" and not has_err:
            self._state[5] = "todo"
        self._paint_states()

    def _paint_states(self):
        marks = {"todo": "", "ok": "  ✔", "warn": "  ⚠", "error": "  ✖"}
        for i, s in enumerate(STEPS):
            self.steps.item(i).setText(s + marks.get(self._state[i], ""))

    def _toggle_preview(self, on: bool):
        self.page_review.preview.setVisible(on)
        for w in (self.page_review.overrides,):
            w.setVisible(on)

    # ── file menu ──
    def _new(self):
        from experiment_spec import ExperimentSpec
        self.load_spec(ExperimentSpec()); self.page_start.name_edit.setText("my_experiment"); self.steps.setCurrentRow(0)

    def _open(self):
        p, _ = QFileDialog.getOpenFileName(self, "Open experiment", self.easy_root, "Experiment (*.experiment.json)")
        if not p:
            return
        try:
            spec = load_spec(p)
        except Exception as e:
            QMessageBox.critical(self, "Open failed", str(e)); return
        self.load_spec(spec)
        self.page_start.out_edit.setText(os.path.dirname(os.path.dirname(p)) or self.easy_root)
        self.steps.setCurrentRow(0)

    def _save(self):
        try:
            res = self.write_all()
            self.statusBar().showMessage(f"saved {res['spec_path']}")
        except Exception as e:
            QMessageBox.critical(self, "Save failed", str(e))

    def _expert(self):
        p = QProcess(self); args = [os.path.join(HERE, "FR_GUI.py"), "-timeout", "0"]
        if self.toml_path:
            args += ["-t", self.toml_path]
        p.start(sys.executable, args); self._expert_proc = p


def main(argv=None):
    argv = list(sys.argv if argv is None else argv)
    app = QApplication(argv)
    w = JuiceWindow(pick_toml(argv))
    w.show()
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
