"""
study_window.py — run several experiments at once, with one consistent change, or sweep a
parameter grid of one experiment (isopleths, sensitivities).  Separate window, same look as the
step GUI and the results viewer; the engine is study.py (no Qt), also driven by study_runner.py.

    python study_window.py [<name>.study.json] [-t <toml>]
"""
from __future__ import annotations

import math
import os
import sys
import threading
from typing import List, Optional

from PyQt6.QtCore import Qt, QProcess, pyqtSignal, QObject
from PyQt6.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QGridLayout, QLabel, QLineEdit,
                             QPushButton, QListWidget, QTableWidget, QTableWidgetItem, QComboBox, QSpinBox, QCheckBox,
                             QGroupBox, QTabWidget, QFileDialog, QMessageBox, QSplitter, QScrollArea, QHeaderView,
                             QProgressBar, QTextEdit, QAbstractItemView, QSizePolicy)
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import numpy as np

import juice_theme as theme
import study as S
from scenario_prep import pick_toml, load_paths, HERE
from experiment_wizard import _info, _more_info

F = theme.PLOT_FONT


def _table(headers: List[str]) -> QTableWidget:
    t = QTableWidget(0, len(headers)); t.setHorizontalHeaderLabels(headers)
    t.horizontalHeader().setStretchLastSection(True); t.verticalHeader().setVisible(False)
    t.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
    return t


class _Events(QObject):
    event = pyqtSignal(str, object)       # kind, RunPlan (or None)


class StudyWindow(QMainWindow):
    def __init__(self, toml_path: str = "", study_path: str = "", experiment: str = ""):
        super().__init__()
        self.toml_path = toml_path
        self.paths = load_paths(toml_path)
        self.study = S.new_study("study", "", toml_path)
        self.study_path = ""
        self.runner: Optional[S.StudyRunner] = None
        self._events = _Events(); self._events.event.connect(self._on_event)
        self._pilot: Optional[S.RunPlan] = None
        self._viewers: list = []
        self.setWindowTitle("JUICE — Study: several experiments or a sweep")
        self.setWindowIcon(theme.app_icon()); self.resize(1600, 940)
        self._build()
        m = self.menuBar().addMenu("File")
        m.addAction("New study", self._new); m.addAction("Open study …", self._open); m.addAction("Save study", self._save)
        m.addAction("Save study as …", self._save_as); m.addSeparator()
        m.addAction("Open study folder", lambda: self._open_path(self.study.folder)); m.addSeparator(); m.addAction("Quit", self.close)
        if study_path:
            self._load(study_path)
        else:
            if experiment and os.path.isfile(experiment):      # opened from the step window: start from that experiment
                self.study.experiments = [experiment]
                base = os.path.basename(experiment)[:-len(".experiment.json")]
                self.study.name = f"{base}_study"
                self.study.folder = os.path.join(os.path.dirname(experiment), "studies", self.study.name)
            self._from_study()
        self._readiness()

    # ── layout ──
    def _build(self):
        left = QWidget(); ll = QVBoxLayout(left); ll.setContentsMargins(6, 6, 6, 6)
        ll.addWidget(_more_info("A study reruns saved experiments many times without touching them: pick the experiments, "
                                "optionally one change for all of them, optionally a sweep — then pilot, estimate, run.", [
            "<b>One experiment, a sweep</b>: add the experiment, leave the variant empty, switch on one or two axes "
            "(e.g. initial NO log 0.1 … 20 ppb × injected VOC × 0.5 … 4) → an ozone isopleth or a sensitivity curve.",
            "<b>Many experiments, one change</b>: add the experiments, put the change in the variant table "
            "(a reaction rate, an add-on, a chamber parameter, scaled injections) → every experiment is rerun with it; "
            "compare against the original runs in the viewer.",
            "<b>Both</b>: the sweep is applied to every experiment; the plot shows the first one, the table all of them.",
            "<b>Nothing changed</b>: just reruns the experiments with this study's solver settings (useful after a mechanism update).",
            "Yellow fields still need input. 'Check' lists what the study would do and what it cannot do (e.g. a reaction that is "
            "not in an experiment's subset). 'Run pilot, then estimate' runs one case alone, measures it and asks before the rest.",
            "Every run lives in &lt;study folder&gt;/runs/&lt;experiment&gt;__&lt;axis&gt;-&lt;value&gt;/ with provenance lines in its .easy; "
            "finished runs are skipped when the study is repeated or extended."], "ⓘ How a study works"))
        g = QGroupBox("Study"); gl = QGridLayout(g)
        gl.addWidget(QLabel("Name"), 0, 0); self.name = QLineEdit(); gl.addWidget(self.name, 0, 1, 1, 2)
        gl.addWidget(QLabel("Folder"), 1, 0); self.folder = QLineEdit(); gl.addWidget(self.folder, 1, 1)
        b = QPushButton("Browse …"); b.clicked.connect(self._browse_folder); gl.addWidget(b, 1, 2)
        gl.addWidget(QLabel("Notes"), 2, 0); self.notes = QLineEdit(); gl.addWidget(self.notes, 2, 1, 1, 2)
        ll.addWidget(g)

        g = QGroupBox("Base experiments (.experiment.json written by the step window)"); gl = QVBoxLayout(g)
        self.exps = QListWidget(); self.exps.setMaximumHeight(110); gl.addWidget(self.exps)
        self.exps.model().rowsInserted.connect(self._readiness); self.exps.model().rowsRemoved.connect(self._readiness)
        hb = QHBoxLayout()
        for label, fn in (("Add …", self._add_exps), ("Add all in a folder …", self._add_folder), ("Remove", self._remove_exp)):
            b = QPushButton(label); b.clicked.connect(fn); hb.addWidget(b)
        hb.addStretch(); gl.addLayout(hb); ll.addWidget(g)

        g = QGroupBox("Variant — applied identically to every experiment (leave empty to run them as they are)"); gl = QVBoxLayout(g)
        self.mods = _table(["Change", "Target", "Value", "Comment"]); self.mods.setMaximumHeight(150)
        self.mods.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeMode.ResizeToContents); gl.addWidget(self.mods)
        hb = QHBoxLayout()
        b = QPushButton("+ change"); b.clicked.connect(lambda: self._add_mod()); hb.addWidget(b)
        b = QPushButton("− row"); b.clicked.connect(lambda: self.mods.removeRow(self.mods.currentRow()) if self.mods.currentRow() >= 0 else None); hb.addWidget(b)
        hb.addStretch(); gl.addLayout(hb)
        gl.addWidget(_more_info("Reaction rate (expression or × factor), extra add-on, chamber parameter, injected amounts, initial values, a condition.",
                                [f"<b>{v[0]}</b>: target = {v[1]}; value = {v[2]}" for v in S.MOD_KINDS.values()] +
                                ["A rate change is written as k[reaction] = CONST(expression) after the mechanism — exactly the EASY way (last definition wins).",
                                 "An add-on rebuilds the MCM subset once per experiment in &lt;study&gt;/chemistry/."], "ⓘ What can be changed"))
        ll.addWidget(g)

        g = QGroupBox("Sweep — one or two axes (every experiment × every grid point)"); gl = QGridLayout(g)
        self.axes_ui = []
        for i in range(2):
            on = QCheckBox(f"Axis {i + 1}"); kind = QComboBox(); kind.addItems(list(S.AXIS_KINDS)); kind.setToolTip("\n".join(f"{k}: {v}" for k, v in S.AXIS_KINDS.items()))
            target = QLineEdit(); target.setPlaceholderText("NO  |  voc  |  C5H8 + OH --> CISOPA  |  hono_const")
            values = QLineEdit(); values.setPlaceholderText("0.5, 1, 2   |   lin 0 20 5   |   log 0.1 20 6")
            gl.addWidget(on, i, 0); gl.addWidget(kind, i, 1); gl.addWidget(target, i, 2); gl.addWidget(values, i, 3)
            for w in (on, kind, target, values):
                (w.toggled if isinstance(w, QCheckBox) else w.currentIndexChanged if isinstance(w, QComboBox) else w.textChanged).connect(self._count)
            self.axes_ui.append((on, kind, target, values))
        gl.setColumnStretch(2, 2); gl.setColumnStretch(3, 3)
        self.count = QLabel(""); self.count.setStyleSheet("color:#555"); gl.addWidget(self.count, 2, 0, 1, 4)
        gl.addWidget(_info("Target: a species for 'initial' (values in ppb); an injection name, an action, 'voc', 'nox' or 'all' for the "
                           "injection kinds; a reaction for 'rate_factor'. Values: a list, 'lin a b n' or 'log a b n'. "
                           "Two axes = a grid of all combinations."), 3, 0, 1, 4)
        ll.addWidget(g)

        g = QGroupBox("Metrics taken from every run"); gl = QGridLayout(g)
        self.metric_boxes = {}
        for i, (k, (label, unit)) in enumerate(S.METRICS.items()):
            cb = QCheckBox(f"{label} [{unit}]"); self.metric_boxes[k] = cb; gl.addWidget(cb, i // 3, i % 3)
        ll.addWidget(g)

        g = QGroupBox("Execution"); gl = QGridLayout(g)
        gl.addWidget(QLabel("Parallel runs"), 0, 0); self.workers = QSpinBox(); self.workers.setRange(1, max(1, os.cpu_count() or 1))
        self.workers.setToolTip("Each run is its own process pinned to one core (BLAS/numba threads = 1).  Leave one core free for yourself."); gl.addWidget(self.workers, 0, 1)
        self.kvpd = QCheckBox("write kv/pd files (large; off for sweeps)"); gl.addWidget(self.kvpd, 0, 2, 1, 2)
        hb = QHBoxLayout()
        self.check_btn = QPushButton("Check"); self.check_btn.clicked.connect(self._check); hb.addWidget(self.check_btn)
        self.dry_btn = QPushButton("Write inputs only"); self.dry_btn.clicked.connect(lambda: self._materialise(run_after=None)); hb.addWidget(self.dry_btn)
        self.pilot_btn = QPushButton("Run pilot, then estimate"); self.pilot_btn.setObjectName("primary"); self.pilot_btn.clicked.connect(self._run_pilot); hb.addWidget(self.pilot_btn)
        self.run_btn = QPushButton("Run all"); self.run_btn.clicked.connect(self._run_all); hb.addWidget(self.run_btn)
        self.stop_btn = QPushButton("Stop"); self.stop_btn.setEnabled(False); self.stop_btn.clicked.connect(self._stop); hb.addWidget(self.stop_btn)
        hb.addStretch(); gl.addLayout(hb, 1, 0, 1, 4)
        self.bar = QProgressBar(); self.bar.setFormat("%v / %m runs"); gl.addWidget(self.bar, 2, 0, 1, 4)
        self.status = QLabel("—"); self.status.setWordWrap(True); self.status.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Preferred); gl.addWidget(self.status, 3, 0, 1, 4)
        self.ready = QLabel(""); self.ready.setWordWrap(True); self.ready.setStyleSheet("color:#6a4c93; font-style: italic")
        self.ready.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Preferred); gl.addWidget(self.ready, 4, 0, 1, 4)
        ll.addWidget(g)
        self.name.textChanged.connect(self._readiness); self.folder.textChanged.connect(self._readiness)
        self.mods.itemChanged.connect(self._readiness)
        ll.addStretch()
        scroll = QScrollArea(); scroll.setWidget(left); scroll.setWidgetResizable(True); scroll.setMinimumWidth(560)

        # right: runs / results / plot / messages
        self.tabs = QTabWidget()
        self.runs = _table(["Run", "Experiment", "Point", "State", "Progress", "Time", "Note"])
        self.runs.itemDoubleClicked.connect(lambda _it: self._open_run())
        rw = QWidget(); rl = QVBoxLayout(rw); rl.setContentsMargins(0, 0, 0, 0); rl.addWidget(self.runs)
        hb = QHBoxLayout()
        for label, fn in (("Open in viewer", self._open_run), ("Open run folder", lambda: self._open_path(self._sel_plan().folder if self._sel_plan() else "")),
                          ("Show log", self._show_log)):
            b = QPushButton(label); b.clicked.connect(fn); hb.addWidget(b)
        hb.addStretch(); rl.addLayout(hb); self.tabs.addTab(rw, "Runs")
        self.results = _table(["Run"]); self.results.itemDoubleClicked.connect(lambda _it: self._open_run(self.results))
        rw = QWidget(); rl = QVBoxLayout(rw); rl.setContentsMargins(0, 0, 0, 0); rl.addWidget(self.results)
        hb = QHBoxLayout()
        b = QPushButton("Refresh"); b.clicked.connect(self._refresh_results); hb.addWidget(b)
        b = QPushButton("Export CSV …"); b.clicked.connect(self._export); hb.addWidget(b)
        hb.addStretch(); rl.addLayout(hb); self.tabs.addTab(rw, "Results")
        pw = QWidget(); pl = QVBoxLayout(pw); pl.setContentsMargins(0, 0, 0, 0)
        hb = QHBoxLayout(); hb.addWidget(QLabel("Metric")); self.plot_metric = QComboBox(); self.plot_metric.currentIndexChanged.connect(self._plot); hb.addWidget(self.plot_metric)
        self.plot_log = QCheckBox("log axes"); self.plot_log.toggled.connect(self._plot); hb.addWidget(self.plot_log)
        b = QPushButton("Save image …"); b.clicked.connect(self._save_image); hb.addWidget(b); hb.addStretch(); pl.addLayout(hb)
        self.fig = Figure(figsize=(8, 5.5)); self.canvas = FigureCanvas(self.fig); pl.addWidget(self.canvas, stretch=1)
        self.plot_note = QLabel(""); self.plot_note.setWordWrap(True); self.plot_note.setStyleSheet("color:#555"); pl.addWidget(self.plot_note)
        self.tabs.addTab(pw, "Plot")
        self.msgs = QTextEdit(); self.msgs.setReadOnly(True); self.tabs.addTab(self.msgs, "Messages")

        split = QSplitter(); split.addWidget(scroll); split.addWidget(self.tabs); split.setStretchFactor(0, 2); split.setStretchFactor(1, 3)
        split.setSizes([620, 980]); self.setCentralWidget(split)

    # ── study ↔ widgets ──
    def _from_study(self):
        # snapshot first: filling the widgets fires change signals that rebuild self.study from the widgets
        st = self.study
        name, folder, notes, exps = st.name, st.folder, st.notes, list(st.experiments)
        mods, axes, metrics, workers, kvpd = list(st.modifications), list(st.axes), list(st.metrics), st.workers, st.export_kv_pd
        self._loading = True
        try:
            self.name.setText(name); self.folder.setText(folder); self.notes.setText(notes)
            self.exps.clear(); self.exps.addItems(exps)
            self.mods.setRowCount(0)
            for m in mods:
                self._add_mod(m)
            for k, cb in self.metric_boxes.items():
                cb.setChecked(k in metrics)
            for i, (on, kind, target, values) in enumerate(self.axes_ui):
                ax = axes[i] if i < len(axes) else None
                on.setChecked(ax is not None)
                if ax:
                    kind.setCurrentText(ax.kind); target.setText(ax.target); values.setText(", ".join(f"{v:g}" for v in ax.values))
                else:
                    target.clear(); values.clear()
            self.workers.setValue(max(1, min(self.workers.maximum(), workers))); self.kvpd.setChecked(kvpd)
        finally:
            self._loading = False
        self.setWindowTitle(f"JUICE — Study · {name}")
        self._count(); self._fill_runs()

    def _to_study(self) -> S.Study:
        st = self.study
        st.name = self.name.text().strip() or "study"; st.folder = self.folder.text().strip(); st.notes = self.notes.text()
        st.experiments = [self.exps.item(i).text() for i in range(self.exps.count())]
        st.modifications = []
        for r in range(self.mods.rowCount()):
            kind = self.mods.cellWidget(r, 0).currentData()
            cell = lambda c: (self.mods.item(r, c).text().strip() if self.mods.item(r, c) else "")
            if cell(1) or kind == "addon":
                st.modifications.append(S.Modification(kind, cell(1), cell(2), cell(3)))
        st.axes = []
        for on, kind, target, values in self.axes_ui:
            if on.isChecked() and target.text().strip():
                try:
                    vals = S.parse_values(values.text())
                except ValueError:
                    vals = []
                st.axes.append(S.Axis(kind.currentText(), target.text().strip(), vals))
        st.metrics = [k for k, cb in self.metric_boxes.items() if cb.isChecked()]
        st.workers = self.workers.value(); st.export_kv_pd = self.kvpd.isChecked()
        if not st.toml:
            st.toml = self.toml_path
        if not st.solver:
            st.solver = S.solver_from_toml(st.toml)
        return st

    def _add_mod(self, m: S.Modification = None):
        r = self.mods.rowCount(); self.mods.insertRow(r)
        cb = QComboBox()
        for k, (label, _t, _v) in S.MOD_KINDS.items():
            cb.addItem(label, k)
        if m:
            cb.setCurrentIndex(max(0, list(S.MOD_KINDS).index(m.kind))) if m.kind in S.MOD_KINDS else None
        self.mods.setCellWidget(r, 0, cb)
        for c, v in enumerate((m.target if m else "", m.value if m else "", m.comment if m else ""), start=1):
            self.mods.setItem(r, c, QTableWidgetItem(v))

    def _readiness(self, *_a):
        """Yellow on what still needs input, and one line saying what the study would do."""
        if getattr(self, "_loading", False):
            return
        missing = []
        theme.mark_missing(self.name, not self.name.text().strip())
        if not self.name.text().strip():
            missing.append("a study name")
        theme.mark_missing(self.exps, self.exps.count() == 0)
        if self.exps.count() == 0:
            missing.append("at least one base experiment")
        n_axes = 0
        for i, (on, kind, target, values) in enumerate(self.axes_ui):
            if on.isChecked():
                bad_t = not target.text().strip()
                try:
                    bad_v = not S.parse_values(values.text())
                except ValueError:
                    bad_v = True
                theme.mark_missing(target, bad_t); theme.mark_missing(values, bad_v)
                if bad_t or bad_v:
                    missing.append(f"{'the target' if bad_t else 'the values'} of axis {i + 1}")
                else:
                    n_axes += 1
            else:
                theme.mark_missing(target, False); theme.mark_missing(values, False)
        n_mods = sum(1 for r in range(self.mods.rowCount()) if self.mods.item(r, 1) and self.mods.item(r, 1).text().strip())
        if self.mods.rowCount() > n_mods:
            missing.append("the target of a change in the variant table")
        if missing:
            self.ready.setText("Still needed: " + ", ".join(missing) + ".")
        else:
            what = []
            if n_mods:
                what.append(f"{n_mods} change(s) applied to every experiment")
            if n_axes:
                what.append(f"a sweep over {n_axes} axis/axes")
            self.ready.setText("Ready: " + (" and ".join(what) if what else "the experiments rerun unchanged") +
                               f" — {self.exps.count()} experiment(s). Next: Check, then 'Run pilot, then estimate'.")

    def _count(self, *_a):
        if getattr(self, "_loading", False):
            return
        self._readiness()
        try:
            st = self._to_study(); n_pts = len(S.grid_points(st.axes)); n = len(st.experiments) * n_pts
            txt = f"{len(st.experiments)} experiment(s) × {n_pts} grid point(s) = {n} runs"
            prior = S.prior_seconds(st)
            if prior:
                txt += f";  last run time of the base experiment(s) ≈ {S.fmt_duration(prior)} → about {S.fmt_duration(S.estimate_seconds(n, prior, st.workers))} with {st.workers} workers (before a pilot)"
            lvl, g = S.check_guards(n)
            if lvl != "ok":
                txt += f"   ⚠ {g}"
            self.count.setText(txt)
        except Exception as e:
            self.count.setText(str(e))

    # ── experiments ──
    def _add_exps(self):
        start = self.paths.get("output_root", "")
        files, _ = QFileDialog.getOpenFileNames(self, "Base experiments", start, "Experiment (*.experiment.json)")
        for f in files:
            if f not in [self.exps.item(i).text() for i in range(self.exps.count())]:
                self.exps.addItem(f)
        self._count()

    def _add_folder(self):
        d = QFileDialog.getExistingDirectory(self, "Folder with experiments", self.paths.get("output_root", ""))
        if not d:
            return
        for root, _dirs, files in os.walk(d):
            if os.path.basename(root) == "runs" or "studies" in root:
                continue
            for f in files:
                if f.endswith(".experiment.json"):
                    p = os.path.join(root, f)
                    if p not in [self.exps.item(i).text() for i in range(self.exps.count())]:
                        self.exps.addItem(p)
        self._count()

    def _remove_exp(self):
        for it in self.exps.selectedItems():
            self.exps.takeItem(self.exps.row(it))
        self._count()

    def _browse_folder(self):
        d = QFileDialog.getExistingDirectory(self, "Study folder", self.folder.text() or self.paths.get("output_root", ""))
        if d:
            self.folder.setText(d)

    # ── file menu ──
    def _new(self):
        self.study = S.new_study("study", "", self.toml_path); self.study_path = ""; self.runner = None; self._from_study()

    def _open(self):
        p, _ = QFileDialog.getOpenFileName(self, "Open study", self.paths.get("output_root", ""), "Study (*.study.json)")
        if p:
            self._load(p)

    def _load(self, p: str):
        try:
            self.study = S.load_study(p); self.study_path = p; self.runner = None
        except Exception as e:
            QMessageBox.critical(self, "Open failed", str(e)); return
        self._from_study(); self._refresh_results(); self._log(f"opened {p}")

    def _save(self, path: str = "") -> str:
        st = self._to_study()
        if not st.folder:
            st.folder = os.path.join(self.paths.get("output_root", HERE), "studies", st.name); self.folder.setText(st.folder)
        os.makedirs(st.folder, exist_ok=True)
        path = path or self.study_path or os.path.join(st.folder, st.name + ".study.json")
        S.save_study(st, path); self.study_path = path; self.setWindowTitle(f"JUICE — Study · {st.name}")
        self.statusBar().showMessage(f"saved {path}", 6000); return path

    def _save_as(self):
        st = self._to_study()
        p, _ = QFileDialog.getSaveFileName(self, "Save study as", os.path.join(st.folder or self.paths.get("output_root", ""), st.name + ".study.json"), "Study (*.study.json)")
        if p:
            self._save(p if p.endswith(".study.json") else os.path.splitext(p)[0] + ".study.json")

    # ── checks, inputs, runs ──
    def _check(self) -> bool:
        st = self._to_study(); msgs = S.validate_study(st)
        self.msgs.clear()
        for lvl, m in msgs:
            self._log(("✖ " if lvl == "error" else "⚠ ") + m)
        if not msgs:
            self._log("✔ nothing to fix")
        self.tabs.setCurrentWidget(self.msgs) if msgs else None
        return not any(lvl == "error" for lvl, _ in msgs)

    def _ensure_runner(self) -> Optional[S.StudyRunner]:
        if not self._check():
            return None
        self._save()
        self.runner = S.StudyRunner(self.study, self.paths.get("include_dir", ""))
        self._fill_runs(); return self.runner

    def _materialise(self, run_after):
        r = self._ensure_runner()
        if r is None:
            return
        self._busy(True); self.status.setText("writing inputs …")

        def work():
            errs = r.materialise_all(lambda k, p, _r: self._events.event.emit(k, p))
            self._events.event.emit("inputs_done", errs)
            if run_after and not errs:
                run_after()
        threading.Thread(target=work, daemon=True).start()

    def _run_pilot(self):
        def after():
            todo = [p for p in self.runner.plans if not S.is_done(p)]
            if not todo:
                self._events.event.emit("all_done", None); return
            self._pilot = todo[0]
            self._events.event.emit("pilot_started", self._pilot)
            self.runner.run_one(self._pilot, lambda k, p, _r: self._events.event.emit(k, p))
            self._events.event.emit("pilot_done", self._pilot)
        self._materialise(after)

    def _run_all(self):
        def after():
            self.runner.start(lambda k, p, _r: self._events.event.emit(k, p), self.workers.value())
        r = self._ensure_runner()
        if r is None:
            return
        n = r.remaining(); lvl, g = S.check_guards(n)
        prior = r.per_run_seconds() or S.prior_seconds(self.study)
        est = f"≈ {S.fmt_duration(S.estimate_seconds(n, prior, self.workers.value()))} with {self.workers.value()} workers" if prior else "no time estimate yet (run a pilot first for one)"
        if lvl != "ok" or not prior:
            if QMessageBox.question(self, "Run all", f"{g}\n{est}\n\nStart?") != QMessageBox.StandardButton.Yes:
                return
        self._materialise(after)

    def _stop(self):
        if self.runner:
            self.runner.stop(); self.status.setText("stopping after the current runs …")

    def _busy(self, on: bool):
        for b in (self.check_btn, self.dry_btn, self.pilot_btn, self.run_btn):
            b.setEnabled(not on)
        self.stop_btn.setEnabled(on)

    # ── events from the runner thread ──
    def _on_event(self, kind: str, p):
        if kind == "inputs_done":
            errs = p
            self._busy(False); self._fill_runs()
            self.status.setText(f"inputs written ({len(errs)} failed)" if errs else "inputs written")
            for e in errs:
                self._log("✖ " + e)
            return
        if kind == "pilot_started":
            self._busy(True); self.status.setText(f"pilot: {p.name} running alone to measure the time …"); return
        if kind == "pilot_done":
            self._busy(False); self._update_row(p); self._refresh_results()
            rest = [x for x in self.runner.plans if not S.is_done(x)]
            if p.status != "done":
                self.status.setText(f"pilot failed: {p.message}"); return
            if not rest:
                self.status.setText("all runs done"); return
            secs, txt = self.runner.estimate(self.workers.value()); lvl, g = S.check_guards(len(rest))
            disk = self.runner.disk_estimate()
            msg = (f"Pilot {p.name}: {S.fmt_duration(p.elapsed)}.\n\n{len(rest)} runs remaining × {S.fmt_duration(self.runner.per_run_seconds())} "
                   f"/ {self.workers.value()} workers ≈ {txt}" + (f"; about {disk / 1e6:.0f} MB of results" if disk and not math.isnan(disk) else "") +
                   (f"\n\n⚠ {g}" if lvl != "ok" else "") + "\n\nRun the remaining runs now?")
            if QMessageBox.question(self, "Estimate", msg) == QMessageBox.StandardButton.Yes:
                self._busy(True); self.runner.start(lambda k, pp, _r: self._events.event.emit(k, pp), self.workers.value(), rest)
            else:
                self.status.setText(f"stopped after the pilot ({txt} estimated for the rest); 'Run all' continues")
            return
        if kind == "all_done":
            self._busy(False); self._refresh_results()
            done = sum(1 for x in self.runner.plans if S.is_done(x)); self.bar.setValue(done)
            self.status.setText(f"{done} / {len(self.runner.plans)} runs done; results in {self.study.folder}")
            S.write_results(self.study, S.collect_results(self.study, self.runner.plans)); return
        if p is not None:
            self._update_row(p)
            if kind in ("finished", "failed", "skipped"):
                done = sum(1 for x in self.runner.plans if x.status == "done"); self.bar.setValue(done)
                secs, txt = self.runner.estimate(self.workers.value())
                self.status.setText(f"{done} / {len(self.runner.plans)} done; per run ≈ {S.fmt_duration(self.runner.per_run_seconds())}; remaining ≈ {txt}")
                if kind == "failed":
                    self._log(f"✖ {p.name}: {p.message}")

    # ── runs table ──
    def _fill_runs(self):
        plans = self.runner.plans if self.runner else S.plan_runs(self._to_study())
        self.runs.setRowCount(0); self.bar.setMaximum(max(1, len(plans))); self.bar.setValue(0)
        for p in plans:
            if S.is_done(p) and p.status == "planned":
                p.status = "done"; p.progress = 1.0; p.elapsed = S._read_info(p).get("elapsed", float("nan"))
            r = self.runs.rowCount(); self.runs.insertRow(r)
            for c, v in enumerate((p.name, p.base_name, ", ".join(f"{k} = {v:g}" for k, v in p.point.items()), p.status, "", "", p.message)):
                it = QTableWidgetItem(v); it.setData(Qt.ItemDataRole.UserRole, p.name); self.runs.setItem(r, c, it)
            bar = QProgressBar(); bar.setRange(0, 100); bar.setValue(int(p.progress * 100)); self.runs.setCellWidget(r, 4, bar)
            self.runs.item(r, 5).setText(S.fmt_duration(p.elapsed) if p.elapsed == p.elapsed else "")
        self.bar.setValue(sum(1 for p in plans if p.status == "done"))
        self.runs.resizeColumnsToContents()

    def _update_row(self, p: S.RunPlan):
        for r in range(self.runs.rowCount()):
            if self.runs.item(r, 0).data(Qt.ItemDataRole.UserRole) == p.name:
                self.runs.item(r, 3).setText(p.status); self.runs.cellWidget(r, 4).setValue(int(p.progress * 100))
                self.runs.item(r, 5).setText(S.fmt_duration(p.elapsed) if p.elapsed == p.elapsed else ""); self.runs.item(r, 6).setText(p.message)
                return

    def _sel_plan(self, table: QTableWidget = None) -> Optional[S.RunPlan]:
        table = table or self.runs
        r = table.currentRow()
        if r < 0 or not table.item(r, 0):
            return None
        name = table.item(r, 0).data(Qt.ItemDataRole.UserRole) or table.item(r, 0).text()
        plans = self.runner.plans if self.runner else S.plan_runs(self._to_study())
        return next((p for p in plans if p.name == name), None)

    def _open_run(self, table: QTableWidget = None):
        p = self._sel_plan(table)
        if p is None:
            return
        if not os.path.isfile(p.cs_path):
            QMessageBox.information(self, "Open", f"{p.name} has no results yet."); return
        args = [os.path.join(HERE, "results_viewer.py"), p.cs_path] + (["--toml", self.toml_path] if self.toml_path else [])
        proc = QProcess(self); proc.start(sys.executable, args); self._viewers.append(proc)

    def _show_log(self):
        p = self._sel_plan()
        if p and os.path.isfile(p.log_path):
            self.msgs.setPlainText(open(p.log_path, encoding="utf-8", errors="replace").read()[-20000:]); self.tabs.setCurrentWidget(self.msgs)

    def _open_path(self, path: str):
        if path and os.path.exists(path):
            QProcess.startDetached("explorer", [os.path.normpath(path)])

    def _log(self, text: str):
        self.msgs.append(text)

    # ── results ──
    def _rows(self) -> List[dict]:
        st = self._to_study()
        plans = self.runner.plans if self.runner else S.plan_runs(st)
        return S.collect_results(st, plans)

    def _refresh_results(self):
        rows = self._rows(); st = self.study
        axis_cols = [a.name() for a in st.axes]
        cols = ["run", "experiment"] + axis_cols + ["status", "elapsed_s"] + st.metrics
        self.results.setColumnCount(len(cols)); self.results.setHorizontalHeaderLabels([S.METRICS.get(c, (c, ""))[0] if c in S.METRICS else c for c in cols])
        self.results.setRowCount(0)
        for row in rows:
            r = self.results.rowCount(); self.results.insertRow(r)
            for c, k in enumerate(cols):
                v = row.get(k)
                txt = "" if v is None or (isinstance(v, float) and math.isnan(v)) else (f"{v:.4g}" if isinstance(v, float) else str(v))
                it = QTableWidgetItem(txt); it.setData(Qt.ItemDataRole.UserRole, row["run"]); self.results.setItem(r, c, it)
        self.results.resizeColumnsToContents()
        cur = self.plot_metric.currentData()
        self.plot_metric.blockSignals(True); self.plot_metric.clear()
        for k in st.metrics:
            self.plot_metric.addItem(f"{S.METRICS[k][0]} [{S.METRICS[k][1]}]", k)
        i = self.plot_metric.findData(cur); self.plot_metric.setCurrentIndex(i if i >= 0 else 0); self.plot_metric.blockSignals(False)
        self._plot()

    def _export(self):
        p, _ = QFileDialog.getSaveFileName(self, "Export results", os.path.join(self.study.folder, f"{self.study.name}_results.csv"), "CSV (*.csv)")
        if p:
            rows = self._rows(); st = self.study
            import csv
            with open(p, "w", newline="", encoding="utf-8") as fh:
                w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()) if rows else ["run"]); w.writeheader(); w.writerows(rows)
            self.statusBar().showMessage(f"exported {p}", 6000)

    def _plot(self, *_a):
        st = self.study; rows = [r for r in self._rows() if r["status"] == "done"]
        key = self.plot_metric.currentData()
        self.fig.clear()
        if not rows or not key:
            self.canvas.draw(); self.plot_note.setText("no finished runs yet"); return
        label, unit = S.METRICS[key]; ax = self.fig.add_subplot(111)
        axes = st.axes; log = self.plot_log.isChecked()
        exps = sorted({r["experiment"] for r in rows})
        if len(axes) >= 2:
            a1, a2 = axes[0].name(), axes[1].name()
            exp = exps[0]
            sub = [r for r in rows if r["experiment"] == exp and r.get(key) is not None and not math.isnan(r[key])]
            x = np.array([r[a1] for r in sub]); y = np.array([r[a2] for r in sub]); z = np.array([r[key] for r in sub])
            if len(sub) >= 4:
                xs, ys = (np.log10(x), np.log10(y)) if log else (x, y)
                try:
                    cs = ax.tricontourf(xs, ys, z, levels=14, cmap="viridis"); self.fig.colorbar(cs, ax=ax, label=f"{label} [{unit}]")
                    ax.tricontour(xs, ys, z, levels=7, colors="white", linewidths=0.6, alpha=0.7)
                except Exception:
                    sc = ax.scatter(xs, ys, c=z, cmap="viridis", s=60); self.fig.colorbar(sc, ax=ax, label=f"{label} [{unit}]")
                ax.scatter(xs, ys, c="k", s=8)
                for xi, yi, zi in zip(xs, ys, z):
                    ax.annotate(f"{zi:.3g}", (xi, yi), fontsize=F["small"], color="white", xytext=(3, 3), textcoords="offset points")
            else:
                sc = ax.scatter(x, y, c=z, cmap="viridis", s=80); self.fig.colorbar(sc, ax=ax, label=f"{label} [{unit}]")
            ax.set_xlabel(a1 + (" (log10)" if log else "")); ax.set_ylabel(a2 + (" (log10)" if log else ""))
            ax.set_title(f"{label} — {exp}" + (f"  (+{len(exps) - 1} more experiments: see Results)" if len(exps) > 1 else ""), loc="left")
            self.plot_note.setText(f"isopleth of {label} over {a1} × {a2}; {len(sub)} grid points, interpolated between them")
        elif len(axes) == 1:
            a1 = axes[0].name(); notes = []
            for exp in exps:
                sub = sorted([r for r in rows if r["experiment"] == exp and r.get(key) is not None], key=lambda r: r[a1])
                x = [r[a1] for r in sub]; y = [r[key] for r in sub]
                ax.plot(x, y, "o-", label=exp)
                sens = S.normalised_sensitivity(x, y)
                if sens:
                    notes.append(f"{exp}: % change of {label} per % change of {a1} between neighbouring points = " + ", ".join("—" if math.isnan(s) else f"{s:+.2f}" for s in sens))
            if log:
                ax.set_xscale("log"); ax.set_yscale("log")
            ax.set_xlabel(a1); ax.set_ylabel(f"{label} [{unit}]"); ax.grid(alpha=0.3); ax.legend(fontsize=F["legend"]); ax.set_title(f"{label} vs {a1}", loc="left")
            self.plot_note.setText("\n".join(notes))
        else:
            ys = [next((r[key] for r in rows if r["experiment"] == e), float("nan")) for e in exps]
            ax.bar(range(len(exps)), ys, color="#6a4c93"); ax.set_xticks(range(len(exps))); ax.set_xticklabels(exps, rotation=20, ha="right", fontsize=F["tick"])
            ax.set_ylabel(f"{label} [{unit}]"); ax.set_title(f"{label} per experiment (variant applied)", loc="left"); ax.grid(axis="y", alpha=0.3)
            self.plot_note.setText("one value per experiment; compare with the base runs in the viewer's Comparison tab")
        self.fig.tight_layout(); self.canvas.draw()

    def _save_image(self):
        p, _ = QFileDialog.getSaveFileName(self, "Save image", os.path.join(self.study.folder or "", f"{self.study.name}_{self.plot_metric.currentData() or 'plot'}.png"), "PNG (*.png);;SVG (*.svg)")
        if p:
            self.fig.savefig(p, dpi=200, bbox_inches="tight")


def main(argv=None):
    argv = list(sys.argv if argv is None else argv)
    theme.taskbar_identity()
    app = QApplication(argv); theme.apply_theme(app)
    files = [a for a in argv[1:] if a.endswith(".study.json")]
    exp = argv[argv.index("--experiment") + 1] if "--experiment" in argv and argv.index("--experiment") + 1 < len(argv) else ""
    w = StudyWindow(pick_toml(argv), files[0] if files else "", exp); w.show(); theme.force_taskbar_icon(w)
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
