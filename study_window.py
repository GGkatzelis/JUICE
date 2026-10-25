"""
study_window.py — two kinds of multi-run studies, one window:

  • Sweep one experiment: its injections, chamber sources and initial values are listed; tick what
    goes on axis 1 / axis 2 (× factor or = value), give the values → a grid of runs (isopleths,
    sensitivities).
  • Several experiments, one change: a list of saved experiments and one generic change applied to
    all of them (reaction rate, add-on, MCM version, chamber parameter, injections by class).

The engine is study.py (no Qt); study_runner.py drives the same study files from the command line.

    python study_window.py [<name>.study.json] [--mode sweep|batch] [--experiment <path> …] [-t <toml>]
"""
from __future__ import annotations

import csv
import re
import math
import os
import sys
import threading
from typing import List, Optional, Tuple

from PyQt6.QtCore import Qt, QProcess, pyqtSignal, QObject
from PyQt6.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QGridLayout, QLabel, QLineEdit,
                             QPushButton, QListWidget, QTableWidget, QTableWidgetItem, QComboBox, QSpinBox, QCheckBox,
                             QGroupBox, QTabWidget, QFileDialog, QMessageBox, QSplitter, QScrollArea, QHeaderView,
                             QProgressBar, QTextEdit, QAbstractItemView, QSizePolicy, QRadioButton, QInputDialog)
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import numpy as np

import juice_theme as theme
import study as S
from scenario_prep import pick_toml, load_paths, HERE
from experiment_spec import load_spec, ACTIONS
from experiment_wizard import _info, _more_info

F = theme.PLOT_FONT
AXIS_CHOICES = ["—", "Axis 1", "Axis 2"]

# rows of the "what to vary" table that come from the chamber, not from an injection
CHAMBER_ROWS = [  # (label, field, base-value field, scale kind, set kind)
    ("Chamber HONO source (light, SHONO)", "hono_const", "chamber_scale", "chamber"),
    ("Dark HONO source rate", "dark_hono_rate", "chamber_scale", "chamber"),
    ("Chamber HCHO source (light)", "hcho_const", "chamber_scale", "chamber"),
    ("Wall loss factor (WALLDEP)", "walldep_value", "chamber_scale", "chamber"),
]


def _table(headers: List[str]) -> QTableWidget:
    t = QTableWidget(0, len(headers)); t.setHorizontalHeaderLabels(headers)
    t.horizontalHeader().setStretchLastSection(True); t.verticalHeader().setVisible(False)
    t.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
    return t


class _Events(QObject):
    event = pyqtSignal(str, object)       # kind, RunPlan (or payload)


class StudyWindow(QMainWindow):
    def __init__(self, toml_path: str = "", study_path: str = "", experiments: List[str] = None, mode: str = ""):
        super().__init__()
        self.toml_path = toml_path
        self.paths = load_paths(toml_path)
        self.study = S.new_study("study", "", toml_path)
        self.study_path = ""
        self.runner: Optional[S.StudyRunner] = None
        self._events = _Events(); self._events.event.connect(self._on_event)
        self._viewers: list = []
        self._vary_rows: List[dict] = []          # sweep mode: one dict per table row (kinds, target, widgets)
        self._loading = False
        self.setWindowTitle("JUICE — Study")
        self.setWindowIcon(theme.app_icon()); self.resize(1640, 960)
        self._build()
        m = self.menuBar().addMenu("File")
        m.addAction("New study", self._new); m.addAction("Open study …", self._open); m.addAction("Save study", self._save)
        m.addAction("Save study as …", self._save_as); m.addSeparator()
        m.addAction("Open study folder", lambda: self._open_path(self.study.folder)); m.addSeparator(); m.addAction("Quit", self.close)
        experiments = [e for e in (experiments or []) if e and os.path.isfile(e)]
        if study_path:
            self._load(study_path)
        else:
            if experiments:
                self.study.experiments = experiments
                base = os.path.basename(experiments[0])[:-len(".experiment.json")]
                self.study.name = f"{base}_sweep" if (mode == "sweep" or len(experiments) == 1 and mode != "batch") else "batch_" + base
                self.study.folder = os.path.join(os.path.dirname(experiments[0]), "studies", self.study.name)
            self._from_study(mode or ("sweep" if len(experiments) <= 1 else "batch"))
        self._readiness()

    # ── layout ──
    def _build(self):
        left = QWidget(); ll = QVBoxLayout(left); ll.setContentsMargins(6, 6, 6, 6)
        ll.addWidget(_more_info("A study reruns saved experiments many times without touching them. Choose the kind of study, fill "
                                "the yellow fields, then Check → pilot → run.", [
            "<b>Sweep one experiment</b>: its injections, chamber sources and initial values are listed. Put what you want to vary on "
            "axis 1 or axis 2 (× factor or = value), give the axis values → a grid of runs: an ozone isopleth or a sensitivity curve. "
            "Several rows on the same axis are varied together (e.g. all VOC injections × the same factor).",
            "<b>Several experiments, one change</b>: a list of saved experiments and one change applied to every one of them — a reaction "
            "rate, an add-on, another MCM version, a chamber parameter, injections by class (voc, nox, …). Every experiment keeps its own "
            "steps; compare with the original runs in the viewer. A generic sweep can be added on top.",
            "<b>No change at all</b> is also possible: the experiments are rerun with their own chemistry files and this study's solver "
            "settings — a reproduction, not a mechanism update. A new MCM version is a chemistry change: use 'Base mechanism' in the variant.",
            "'Check' lists what the study will do and what it cannot do (a reaction that is not in an experiment's subset, an injection "
            "class that does not exist there). 'Run pilot, then estimate' runs one case alone, measures it and asks before the rest.",
            "Every run lives in &lt;study folder&gt;/runs/&lt;experiment&gt;__&lt;axis&gt;-&lt;value&gt;/ with provenance lines in its .easy; "
            "finished runs are skipped when a study is repeated or extended."], "ⓘ How a study works"))

        g = QGroupBox("Study"); gl = QGridLayout(g)
        self.mode_sweep = QRadioButton("Sweep one experiment"); self.mode_batch = QRadioButton("Several experiments, one change")
        self.mode_sweep.setChecked(True); self.mode_sweep.toggled.connect(self._mode_changed)
        hb = QHBoxLayout(); hb.addWidget(self.mode_sweep); hb.addWidget(self.mode_batch); hb.addStretch(); gl.addLayout(hb, 0, 0, 1, 3)
        gl.addWidget(QLabel("Name"), 1, 0); self.name = QLineEdit(); gl.addWidget(self.name, 1, 1, 1, 2)
        gl.addWidget(QLabel("Folder"), 2, 0); self.folder = QLineEdit(); gl.addWidget(self.folder, 2, 1)
        b = QPushButton("Browse …"); b.clicked.connect(self._browse_folder); gl.addWidget(b, 2, 2)
        gl.addWidget(QLabel("Notes"), 3, 0); self.notes = QLineEdit(); gl.addWidget(self.notes, 3, 1, 1, 2)
        ll.addWidget(g)

        # ── sweep mode ──
        self.sweep_box = QGroupBox("What to vary in this experiment"); sl = QVBoxLayout(self.sweep_box)
        hb = QHBoxLayout(); hb.addWidget(QLabel("Experiment")); self.exp_path = QLineEdit(); self.exp_path.setReadOnly(True); hb.addWidget(self.exp_path, stretch=1)
        b = QPushButton("Choose …"); b.clicked.connect(self._choose_exp); hb.addWidget(b); sl.addLayout(hb)
        self.exp_status = _info("no experiment chosen"); sl.addWidget(self.exp_status)
        self.vary = _table(["Vary on", "What", "Kind", "In the experiment", "How"])
        self.vary.setMinimumHeight(230)
        for c in (0, 2, 4):
            self.vary.horizontalHeader().setSectionResizeMode(c, QHeaderView.ResizeMode.ResizeToContents)
        sl.addWidget(self.vary)
        hb = QHBoxLayout()
        b = QPushButton("+ initial value of a species …"); b.clicked.connect(self._add_initial_row); hb.addWidget(b)
        b = QPushButton("+ reaction rate × factor …"); b.clicked.connect(self._add_rate_row); hb.addWidget(b)
        b = QPushButton("+ new injection …"); b.setToolTip("An injection the experiment does not have (e.g. NO at 11:00 over 5 min); its amount goes on an axis, e.g. 0, 2, 5, 10 ppb.")
        b.clicked.connect(self._add_new_injection_row); hb.addWidget(b)
        b = QPushButton("− custom row"); b.clicked.connect(self._remove_custom_row); hb.addWidget(b)
        hb.addStretch(); sl.addLayout(hb)
        vg = QGridLayout()
        self.axis_values = []
        for i in range(2):
            vg.addWidget(QLabel(f"Axis {i + 1} values"), 2 * i, 0)
            e = QLineEdit(); e.setPlaceholderText("a list: 0.5, 1, 2   |   lin 0 20 5 = 5 values from 0 to 20   |   log 0.1 20 6 = 6 log-spaced values from 0.1 to 20")
            e.setMinimumWidth(420); e.textChanged.connect(self._count)
            vg.addWidget(e, 2 * i, 1); self.axis_values.append(e)
            lb = QLabel(""); lb.setStyleSheet("color:#555"); lb.setWordWrap(True); lb.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Preferred)
            vg.addWidget(lb, 2 * i + 1, 1); setattr(self, f"axis_desc{i}", lb)
        vg.setColumnStretch(1, 1); sl.addLayout(vg)
        sl.addWidget(_info("× factor multiplies what the experiment has (1 = unchanged); = value replaces it (ppb for injections and "
                           "initial values). 'All VOC injections together' keeps the VOC mix and scales it as a whole. "
                           "The chamber HONO source is the main dark/light NO source besides injections — vary it here, not as an injection."))
        ll.addWidget(self.sweep_box)

        # ── batch mode ──
        self.batch_box = QGroupBox("Experiments and the change applied to all of them"); bl = QVBoxLayout(self.batch_box)
        self.exps = QListWidget(); self.exps.setMaximumHeight(110); bl.addWidget(self.exps)
        self.exps.model().rowsInserted.connect(self._readiness); self.exps.model().rowsRemoved.connect(self._readiness)
        hb = QHBoxLayout()
        for label, fn in (("Add …", self._add_exps), ("Add all in a folder …", self._add_folder), ("Remove", self._remove_exp)):
            b = QPushButton(label); b.clicked.connect(fn); hb.addWidget(b)
        hb.addStretch(); bl.addLayout(hb)
        ll.addWidget(self.batch_box)

        # ── the change table: in batch mode the point of the study, in sweep mode an optional fixed change on every grid point ──
        self.change_box = QGroupBox(""); bl = QVBoxLayout(self.change_box)
        self.change_box.toggled.connect(lambda on: (self._mods_widgets_visible(on), self._count()))
        self.mods = _table(["Change", "Target", "Value", "Comment"]); self.mods.setMaximumHeight(150)
        self.mods.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeMode.ResizeToContents); bl.addWidget(self.mods)
        self.mods.itemChanged.connect(self._readiness)
        hb = QHBoxLayout()
        b = QPushButton("+ change"); b.clicked.connect(lambda: (self._add_mod(), self._readiness())); hb.addWidget(b)
        b = QPushButton("− row"); b.clicked.connect(lambda: (self.mods.removeRow(self.mods.currentRow()) if self.mods.currentRow() >= 0 else None, self._readiness())); hb.addWidget(b)
        hb.addStretch(); bl.addLayout(hb)
        bl.addWidget(_more_info("Targets are generic because every experiment has its own steps: a reaction, an add-on id, a chamber "
                                "field, an injection class (voc, nox, no, o3, mixture, all) or an injection name present in all of them.",
                                [f"<b>{v[0]}</b>: target = {v[1]}; value = {v[2]}" for v in S.MOD_KINDS.values()] +
                                ["A rate change is written as k[reaction] = CONST(expression) after the mechanism — the EASY way (last definition wins).",
                                 "An add-on or another base mechanism rebuilds the MCM subset once per experiment in &lt;study&gt;/chemistry/."],
                                "ⓘ What can be changed"))
        self._mods_widgets = [bl.itemAt(i).widget() or bl.itemAt(i).layout() for i in range(bl.count())]
        ll.addWidget(self.change_box)
        self.gen_sweep = QGroupBox("Also sweep a generic parameter in every experiment"); self.gen_sweep.setCheckable(True); self.gen_sweep.setChecked(False)
        gl = QGridLayout(self.gen_sweep); self.axes_ui = []
        for i in range(2):
            on = QCheckBox(f"Axis {i + 1}"); kind = QComboBox(); kind.addItems(list(S.AXIS_KINDS)); kind.setToolTip("\n".join(f"{k}: {v}" for k, v in S.AXIS_KINDS.items()))
            target = QLineEdit(); target.setPlaceholderText("NO  |  voc  |  C5H8 + OH --> CISOPA  |  hono_const")
            values = QLineEdit(); values.setPlaceholderText("0.5, 1, 2   |   lin 0 20 5   |   log 0.1 20 6")
            gl.addWidget(on, i, 0); gl.addWidget(kind, i, 1); gl.addWidget(target, i, 2); gl.addWidget(values, i, 3)
            for w in (on, kind, target, values):
                (w.toggled if isinstance(w, QCheckBox) else w.currentIndexChanged if isinstance(w, QComboBox) else w.textChanged).connect(self._count)
            self.axes_ui.append((on, kind, target, values))
        gl.setColumnStretch(2, 2); gl.setColumnStretch(3, 3)
        self.gen_sweep.toggled.connect(self._count); ll.addWidget(self.gen_sweep)

        self.count = QLabel(""); self.count.setStyleSheet("color:#555"); self.count.setWordWrap(True); ll.addWidget(self.count)

        g = QGroupBox("Results read from every run"); gl = QGridLayout(g)
        note = _info("Outputs, not inputs: after each run these numbers are read from its result file and collected in the "
                     "Results table, the CSV and the isopleth / sensitivity plot. They change nothing in the model.")
        note.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Preferred); gl.addWidget(note, 0, 0, 1, 2)
        self.metric_boxes = {}
        for i, (k, (label, unit)) in enumerate(S.METRICS.items()):
            cb = QCheckBox(f"{label} [{unit}]"); self.metric_boxes[k] = cb; gl.addWidget(cb, 1 + i // 2, i % 2)
        ll.addWidget(g)

        g = QGroupBox("Execution"); gl = QGridLayout(g)
        gl.addWidget(QLabel("Parallel runs"), 0, 0); self.workers = QSpinBox(); self.workers.setRange(1, max(1, os.cpu_count() or 1))
        self.workers.setToolTip("Each run is its own process pinned to one core (BLAS/numba threads = 1).  Leave one core free for yourself."); gl.addWidget(self.workers, 0, 1)
        self.kvpd = QCheckBox("write kv/pd files (large; off for sweeps)"); gl.addWidget(self.kvpd, 0, 2, 1, 2)
        hb = QHBoxLayout()
        self.check_btn = QPushButton("Check"); self.check_btn.clicked.connect(self._check); hb.addWidget(self.check_btn)
        self.check_btn.setToolTip("Lists what the study will run and anything it cannot do. Writes and runs nothing.")
        self.dry_btn = QPushButton("Write inputs only"); self.dry_btn.clicked.connect(lambda: self._materialise(run_after=None)); hb.addWidget(self.dry_btn)
        self.dry_btn.setToolTip("Writes every run's .easy, .enz and scenario JSON into its run folder, runs nothing.\n"
                                "To inspect a run before spending CPU, or to run the study elsewhere with study_runner.py.")
        self.pilot_btn = QPushButton("Run pilot, then estimate"); self.pilot_btn.setObjectName("primary"); self.pilot_btn.clicked.connect(self._run_pilot); hb.addWidget(self.pilot_btn)
        self.pilot_btn.setToolTip("Runs the first case alone, measures its time and disk size, then asks before starting the rest.")
        self.run_btn = QPushButton("Run all"); self.run_btn.clicked.connect(self._run_all); hb.addWidget(self.run_btn)
        self.run_btn.setToolTip("Runs every case that has no finished result yet, in parallel.")
        self.stop_btn = QPushButton("Stop"); self.stop_btn.setEnabled(False); self.stop_btn.clicked.connect(self._stop); hb.addWidget(self.stop_btn)
        self.stop_btn.setToolTip("Stops the running cases; finished ones are kept and skipped next time.")
        hb.addStretch(); gl.addLayout(hb, 1, 0, 1, 4)
        exec_help = _more_info("Usual order: Check → Run pilot, then estimate → yes. Finished runs are never redone.", [
            "<b>Check</b>: lists what will run (axes, changes) and what cannot work, e.g. a reaction not in an experiment's subset. Writes nothing.",
            "<b>Write inputs only</b>: writes every run's .easy, .enz and scenario JSON into &lt;study&gt;/runs/&lt;run&gt;/ without running. "
            "Use it to open one .easy and confirm the change is where you expect, to see the disk space a big sweep needs, or to copy the "
            "study to a faster machine and run it there: <code>python study_runner.py &lt;name&gt;.study.json --workers 8 --yes</code>.",
            "<b>Run pilot, then estimate</b>: runs the first case alone, then shows <i>remaining runs × pilot time ÷ parallel runs</i> and the "
            "expected disk space, and asks before starting the rest.",
            "<b>Run all</b>: runs everything not yet finished. Without a pilot the estimate comes from the base experiment's last run, if known.",
            "<b>Stop</b>: kills the running cases; finished ones stay. Running the study again continues where it stopped, and adding grid "
            "points later only runs the new ones.",
            "<b>Parallel runs</b>: each run is its own process on one core. Leave a core free if you keep working on the machine.",
            "<b>kv/pd files</b>: rates and production/loss budgets for the process analysis in the viewer; large, so off for sweeps."], "ⓘ What the buttons do")
        exec_help.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Preferred); gl.addWidget(exec_help, 5, 0, 1, 4)
        self.bar = QProgressBar(); self.bar.setFormat("%v / %m runs"); gl.addWidget(self.bar, 2, 0, 1, 4)
        self.status = QLabel("—"); self.status.setWordWrap(True); self.status.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Preferred); gl.addWidget(self.status, 3, 0, 1, 4)
        self.ready = QLabel(""); self.ready.setWordWrap(True); self.ready.setStyleSheet("color:#6a4c93; font-style: italic")
        self.ready.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Preferred); gl.addWidget(self.ready, 4, 0, 1, 4)
        ll.addWidget(g)
        ll.addStretch()
        self.name.textChanged.connect(self._readiness); self.folder.textChanged.connect(self._readiness)
        scroll = QScrollArea(); scroll.setWidget(left); scroll.setWidgetResizable(True); scroll.setMinimumWidth(600)

        # right: runs / results / plot / messages
        self.tabs = QTabWidget()
        self.runs = _table(["Run", "Experiment", "Point", "State", "Progress", "Time", "Note"])
        self.runs.itemDoubleClicked.connect(lambda _it: self._open_run())
        rw = QWidget(); rl = QVBoxLayout(rw); rl.setContentsMargins(0, 0, 0, 0)
        rl.addWidget(_more_info("One row per run: experiment × grid point. Double-click a finished run to open it in the viewer.", [
            "<b>State</b>: planned (nothing written) → ready (inputs written) → running → done | failed. Done = a result newer than its inputs "
            "and a clean exit; such runs are skipped when the study runs again.",
            "<b>Point</b>: the axis values of this run, e.g. VOC_x = 2, NOx_x = 0.5 (× factors) or NO_t0_ppb = 5 (a value).",
            "<b>Time</b>: wall-clock time of the run; the median of finished runs drives the remaining-time estimate.",
            "<b>Note</b>: what could not be applied (e.g. archive unreachable → the base experiment's .enz reused) or why a run failed.",
            "<b>Show log</b>: the solver output of the selected run; <b>Open run folder</b>: its .easy, .enz, JSON and results."], "ⓘ About this table"))
        rl.addWidget(self.runs)
        hb = QHBoxLayout()
        for label, fn in (("Open in viewer", self._open_run), ("Open run folder", lambda: self._open_path(self._sel_plan().folder if self._sel_plan() else "")),
                          ("Show log", self._show_log)):
            b = QPushButton(label); b.clicked.connect(fn); hb.addWidget(b)
        hb.addStretch(); rl.addLayout(hb); self.tabs.addTab(rw, "Runs")
        self.results = _table(["Run"]); self.results.itemDoubleClicked.connect(lambda _it: self._open_run(self.results))
        rw = QWidget(); rl = QVBoxLayout(rw); rl.setContentsMargins(0, 0, 0, 0)
        rl.addWidget(_more_info("The numbers read from every finished run, one row per run. Double-click a row to open that run.", [
            "Columns: the run, its experiment, the axis values, then the results ticked under 'Results read from every run'.",
            "Concentrations in ppb come from number density ÷ M; radicals (OH, HO2, RO2) stay in molecules cm⁻³.",
            "'While lit' = the period with the roof open or the lamps on (HV > 0.5); 'gain while lit' = value at the end of that "
            "period minus value at its start. Ox = O3 + NO2, so it does not move when NO titrates O3.",
            "<b>Export CSV</b> writes the table; the same file is saved as &lt;study&gt;_results.csv after every complete run of the study.",
            "Changed which results are ticked? Press <b>Refresh</b>: the numbers are recomputed from the stored values, no rerun needed."], "ⓘ About this table"))
        rl.addWidget(self.results)
        hb = QHBoxLayout()
        b = QPushButton("Refresh"); b.clicked.connect(self._refresh_results); hb.addWidget(b)
        b = QPushButton("Export CSV …"); b.clicked.connect(self._export); hb.addWidget(b)
        hb.addStretch(); rl.addLayout(hb); self.tabs.addTab(rw, "Results")
        pw = QWidget(); pl = QVBoxLayout(pw); pl.setContentsMargins(0, 0, 0, 0)
        pl.addWidget(_more_info("What is drawn depends on the number of axes: two → isopleth, one → curves, none → one bar per experiment.", [
            "<b>Two axes</b>: filled contours of the chosen result over axis 1 × axis 2 for the first experiment, interpolated between the "
            "grid points (black dots, value printed next to each). Few points give coarse contours: 4 × 4 or more for a usable isopleth.",
            "<b>One axis</b>: the result against the axis value, one line per experiment. The note below gives the normalised sensitivity "
            "between neighbouring points: % change of the result per % change of the parameter (+1 proportional, 0 insensitive, "
            "negative = the result falls as the parameter rises).",
            "<b>No axis</b> (several experiments, one change): one bar per experiment; compare each with its original run in the viewer.",
            "<b>log axes</b> for log-spaced sweeps (log 0.1 20 6). <b>Save image</b> writes PNG or SVG."], "ⓘ About this plot"))
        hb = QHBoxLayout(); hb.addWidget(QLabel("Metric")); self.plot_metric = QComboBox(); self.plot_metric.currentIndexChanged.connect(self._plot); hb.addWidget(self.plot_metric)
        self.plot_log = QCheckBox("log axes"); self.plot_log.toggled.connect(self._plot); hb.addWidget(self.plot_log)
        b = QPushButton("Save image …"); b.clicked.connect(self._save_image); hb.addWidget(b); hb.addStretch(); pl.addLayout(hb)
        self.fig = Figure(figsize=(8, 5.5)); self.canvas = FigureCanvas(self.fig); pl.addWidget(self.canvas, stretch=1)
        self.plot_note = QLabel(""); self.plot_note.setWordWrap(True); self.plot_note.setStyleSheet("color:#555"); pl.addWidget(self.plot_note)
        self.tabs.addTab(pw, "Plot")
        self.msgs = QTextEdit(); self.msgs.setReadOnly(True); self.tabs.addTab(self.msgs, "Messages")

        split = QSplitter(); split.addWidget(scroll); split.addWidget(self.tabs); split.setStretchFactor(0, 65); split.setStretchFactor(1, 35)
        split.setSizes([1066, 574]); self.setCentralWidget(split)
        self._mode_changed()

    def _mode_changed(self, *_a):
        sweep = self.mode_sweep.isChecked()
        self.sweep_box.setVisible(sweep); self.batch_box.setVisible(not sweep); self.gen_sweep.setVisible(not sweep)
        if sweep:
            self.change_box.setTitle("Also apply a fixed change to every run (optional)")
            self.change_box.setCheckable(True)
            self.change_box.setChecked(self.mods.rowCount() > 0)
        else:
            self.change_box.setCheckable(False)
            self.change_box.setTitle("Change applied to every experiment")
        self._mods_widgets_visible(self.change_box.isChecked() if sweep else True)
        if not self._loading:
            self._count()

    def _mods_widgets_visible(self, on: bool):
        for w in getattr(self, "_mods_widgets", []):
            if hasattr(w, "setVisible"):
                w.setVisible(on)
            else:
                for j in range(w.count()):
                    if w.itemAt(j).widget():
                        w.itemAt(j).widget().setVisible(on)

    def _mods_active(self) -> bool:
        return (not self.mode_sweep.isChecked()) or self.change_box.isChecked()

    # ── sweep mode: the experiment's own knobs ──
    def _choose_exp(self):
        p, _ = QFileDialog.getOpenFileName(self, "Experiment to sweep", self.paths.get("output_root", ""), "Experiment (*.experiment.json)")
        if p:
            self._load_exp(p)

    def _load_exp(self, path: str, assignments: dict = None):
        """Fill the 'what to vary' table from the experiment; assignments = {(kind, target): (axis index, how)} to restore."""
        self.exp_path.setText(path); self.vary.setRowCount(0); self._vary_rows = []
        try:
            spec = load_spec(path)
        except Exception as e:
            self.exp_status.setText(f"cannot read: {e}"); return
        self._spec = spec
        inj = [e for e in spec.injections if e.products.strip()]
        amount = lambda e: f"{e.name} {e.total_ppb} ppb at {e.start or '—'}"
        # the two groups first: true injections only (VOC = voc + mixture actions, NOx = NO and NO2 injections incl. their
        # impurity channels); chamber sources are listed separately below and are never part of a group
        vocs = [e for e in inj if e.action_key() in ("voc", "mixture")]
        if vocs:
            self._vary_add("All VOC injections together: " + ", ".join(e.name for e in vocs), "group of injections", "; ".join(amount(e) for e in vocs),
                           ("injection_scale", ",".join(e.name for e in vocs)), None, hows=("× factor",), label="VOC")
        nox = [e for e in inj if e.action_key() in ("no", "no2")]
        if nox:
            self._vary_add("All NOx injections together: " + ", ".join(e.name for e in nox), "group of injections", "; ".join(amount(e) for e in nox),
                           ("injection_scale", ",".join(e.name for e in nox)), None, hows=("× factor",), label="NOx")
        for e in inj:
            self._vary_add(e.name, ACTIONS.get(e.action_key(), {}).get("label", e.action_key()), f"{e.total_ppb} ppb at {e.start or '—'}" + (f", {e.duration_min:g} min" if e.duration_min else ""),
                           ("injection_scale", e.name), ("injection_set", e.name), hows=("× factor", "= ppb"))
        ch = spec.chamber
        for label, fld, k_scale, k_set in CHAMBER_ROWS:
            if fld == "hono_const" and not ch.hono:
                continue
            if fld == "hcho_const" and not ch.hcho:
                continue
            self._vary_add(label, "chamber", str(getattr(ch, fld)), (k_scale, fld), (k_set, fld), hows=("× factor", "= value"))
        for sp, expr in spec.initials:
            self._vary_add(f"Initial {sp}", "initial value", expr, None, ("initial", sp), hows=("= ppb",), custom=True)
        self.exp_status.setText(f"{spec.name}: {len(inj)} injection(s)" + (", chamber sources" if ch.hono or ch.hcho else "") +
                                f"; {'measured conditions' if spec.conditions.source in ('archive', 'file') else 'manual conditions'}; "
                                f"{len(spec.precursors)} compound(s)" if spec.precursors else f"{spec.name}: {len(inj)} injection(s)")
        placed = set()                                   # a saved member goes to the first row that offers it
        for row in self._vary_rows:
            for key in ("scale", "set"):
                kt = row.get(key)
                if assignments and kt in assignments and kt not in placed:
                    ax, how = assignments[kt]
                    row["axis"].setCurrentIndex(ax); row["how"].setCurrentText(how); placed.add(kt)
                    break
        self._count()

    def _vary_add(self, what, kind, base, scale: Optional[Tuple[str, str]], setk: Optional[Tuple[str, str]], hows, label: str = "", custom: bool = False):
        r = self.vary.rowCount(); self.vary.insertRow(r)
        axis = QComboBox(); axis.addItems(AXIS_CHOICES); axis.currentIndexChanged.connect(self._count); self.vary.setCellWidget(r, 0, axis)
        for c, v in ((1, what), (2, kind), (3, base)):
            it = QTableWidgetItem(v); it.setFlags(it.flags() & ~Qt.ItemFlag.ItemIsEditable); self.vary.setItem(r, c, it)
        how = QComboBox(); how.addItems(list(hows)); how.currentIndexChanged.connect(self._count); self.vary.setCellWidget(r, 4, how)
        self._vary_rows.append({"what": what, "scale": scale, "set": setk, "axis": axis, "how": how, "label": label, "custom": custom})

    def _add_initial_row(self):
        sp, ok = QInputDialog.getText(self, "Initial value", "Species (MCM name), its value at t0 goes on an axis in ppb:")
        if ok and sp.strip():
            self._vary_add(f"Initial {sp.strip().upper()}", "initial value", "0 (not set)", None, ("initial", sp.strip().upper()), hows=("= ppb",), custom=True); self._count()

    def _add_rate_row(self):
        rx, ok = QInputDialog.getText(self, "Reaction rate", "Reaction as in the mechanism, e.g.  C5H8 + OH --> CISOPA  (its MCM rate × the axis factor):")
        if ok and rx.strip():
            self._vary_add(rx.strip(), "reaction rate", "MCM rate", ("rate_factor", rx.strip()), None, hows=("× factor",), custom=True); self._count()

    def _add_new_injection_row(self):
        sp, ok = QInputDialog.getText(self, "New injection", "Compound to inject (NO, NO2, O3, CO or an MCM name):", text="NO")
        if not ok or not sp.strip():
            return
        start, ok = QInputDialog.getText(self, "New injection", "Start time HH:MM (UTC):",
                                         text=(self._spec.conditions.roof_open[0][0] if getattr(self, "_spec", None) and self._spec.conditions.roof_open else "11:00"))
        if not ok or not re.fullmatch(r"\d{1,2}:\d{2}", start.strip()):
            return
        _sp, _st, default_dur = S.parse_new_injection(f"{sp}|{start}|")
        dur, ok = QInputDialog.getDouble(self, "New injection", "Duration of the pulse [min]:", default_dur, 0.1, 600, 1)
        if not ok:
            return
        target = f"{sp.strip().upper()}|{start.strip()}|{dur:g}"
        self._vary_add(f"New injection: {sp.strip().upper()} at {start.strip()} over {dur:g} min", "new injection", "not in the experiment",
                       None, ("injection_new", target), hows=("= ppb",), custom=True)
        self._vary_rows[-1]["axis"].setCurrentIndex(1 if not any(r["axis"].currentIndex() == 1 for r in self._vary_rows[:-1]) else 2)
        if not self.axis_values[self._vary_rows[-1]["axis"].currentIndex() - 1].text().strip():
            self.axis_values[self._vary_rows[-1]["axis"].currentIndex() - 1].setText("0, 2, 5, 10, 20")
        self._count()

    def _remove_custom_row(self):
        r = self.vary.currentRow()
        if r >= 0 and self._vary_rows[r]["custom"]:
            self.vary.removeRow(r); del self._vary_rows[r]; self._count()

    def _sweep_axes(self) -> List[S.Axis]:
        axes = []
        for i in range(2):
            members = []
            for row in self._vary_rows:
                if row["axis"].currentIndex() == i + 1:
                    how = row["how"].currentText()
                    kt = row["scale"] if how.startswith("×") else row["set"]
                    if kt and kt not in [m[0] for m in members]:
                        members.append((kt, row["label"] or row["what"]))
            if not members:
                continue
            try:
                vals = S.parse_values(self.axis_values[i].text())
            except ValueError:
                vals = []
            (k0, t0), lab0 = members[0]
            scale = k0 in ("injection_scale", "chamber_scale", "rate_factor")
            label = (S._safe(lab0) if len(members) == 1 and self._vary_rows and any(r["label"] for r in self._vary_rows if (r["scale"] == (k0, t0) or r["set"] == (k0, t0))) else S._safe(t0.split(",")[0]) if len(members) == 1 else S._safe("+".join(S._safe(t.split(",")[0]) for (_k, t), _l in members)))
            if k0 == "injection_new" and len(members) == 1:
                sp_, st_, _du = S.parse_new_injection(t0)
                label = f"{sp_}_added_{st_.replace(':', '')}"
            label += "_x" if scale else ("_t0_ppb" if k0 == "initial" else "_ppb" if k0 in ("injection_set", "injection_new") else "")
            axes.append(S.Axis(k0, t0, vals, label, [[k, t] for (k, t), _l in members[1:]]))
        return axes

    # ── study ↔ widgets ──
    def _from_study(self, mode: str = ""):
        st = self.study
        name, folder, notes, exps = st.name, st.folder, st.notes, list(st.experiments)
        mods, axes, metrics, workers, kvpd = list(st.modifications), list(st.axes), list(st.metrics), st.workers, st.export_kv_pd
        mode = mode or st.mode
        sweep = mode == "sweep" if mode else (len(exps) == 1 and not mods)
        self._loading = True
        try:
            (self.mode_sweep if sweep else self.mode_batch).setChecked(True)
            self.name.setText(name); self.folder.setText(folder); self.notes.setText(notes)
            self.exps.clear(); self.exps.addItems(exps)
            self.mods.setRowCount(0)
            for m in mods:
                self._add_mod(m)
            for k, cb in self.metric_boxes.items():
                cb.setChecked(k in metrics)
            if sweep:
                assignments = {}
                for i, ax in enumerate(axes[:2]):
                    for k, t in ax.members():
                        assignments[(k, t)] = (i + 1, "× factor" if k in ("injection_scale", "chamber_scale", "rate_factor") else ("= ppb" if k in ("injection_set", "initial", "injection_new") else "= value"))
                    self.axis_values[i].setText(", ".join(f"{v:g}" for v in ax.values))
                if exps:
                    self._load_exp(exps[0], assignments)
                    for (k, t), (i, how) in assignments.items():      # members with no row in the table (e.g. an initial value) → custom rows
                        if any(r["scale"] == (k, t) or r["set"] == (k, t) for r in self._vary_rows):
                            continue
                        if not any(r["scale"] == (k, t) or r["set"] == (k, t) for r in self._vary_rows):
                            if k == "initial":
                                self._vary_add(f"Initial {t}", "initial value", "0 (not set)", None, (k, t), hows=("= ppb",), custom=True)
                            elif k == "rate_factor":
                                self._vary_add(t, "reaction rate", "MCM rate", (k, t), None, hows=("× factor",), custom=True)
                            elif k in ("injection_scale", "injection_set"):     # generic class from a CLI-made study (voc, nox, all, …)
                                self._vary_add(f"Injections: {t}", "injection class", "", ("injection_scale", t), ("injection_set", t), hows=("× factor", "= ppb"), custom=True)
                            elif k == "injection_new":
                                sp_, st_, du_ = S.parse_new_injection(t)
                                self._vary_add(f"New injection: {sp_} at {st_} over {du_:g} min", "new injection", "not in the experiment",
                                               None, (k, t), hows=("= ppb",), custom=True)
                            elif k in ("chamber_scale", "chamber"):
                                self._vary_add(f"Chamber {t}", "chamber", "", ("chamber_scale", t), ("chamber", t), hows=("× factor", "= value"), custom=True)
                            else:
                                continue
                            self._vary_rows[-1]["axis"].setCurrentIndex(i); self._vary_rows[-1]["how"].setCurrentText(how)
                for i in range(2):
                    on, kind, target, values = self.axes_ui[i]; on.setChecked(False); target.clear(); values.clear()
            else:
                self.gen_sweep.setChecked(bool(axes))
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
        self._mode_changed(); self._count(); self._fill_runs()

    def _to_study(self) -> S.Study:
        st = self.study
        st.name = self.name.text().strip() or "study"; st.folder = self.folder.text().strip(); st.notes = self.notes.text()
        st.modifications = []
        if self._mods_active():
            for r in range(self.mods.rowCount()):
                kind = self.mods.cellWidget(r, 0).currentData()
                cell = lambda c: (self.mods.item(r, c).text().strip() if self.mods.item(r, c) else "")
                if cell(1):
                    st.modifications.append(S.Modification(kind, cell(1), cell(2), cell(3)))
        if self.mode_sweep.isChecked():
            st.mode = "sweep"
            st.experiments = [self.exp_path.text()] if self.exp_path.text() else []
            st.axes = self._sweep_axes()
        else:
            st.mode = "batch"
            st.experiments = [self.exps.item(i).text() for i in range(self.exps.count())]
            st.axes = []
            if self.gen_sweep.isChecked():
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
        if m and m.kind in S.MOD_KINDS:
            cb.setCurrentIndex(list(S.MOD_KINDS).index(m.kind))
        self.mods.setCellWidget(r, 0, cb)
        for c, v in enumerate((m.target if m else "", m.value if m else "", m.comment if m else ""), start=1):
            self.mods.setItem(r, c, QTableWidgetItem(v))

    # ── readiness: yellow marks + one honest line ──
    def _readiness(self, *_a):
        if self._loading:
            return
        missing = []
        theme.mark_missing(self.name, not self.name.text().strip())
        if not self.name.text().strip():
            missing.append("a study name")
        sweep = self.mode_sweep.isChecked()
        n_axes = n_mods = 0
        if sweep:
            theme.mark_missing(self.exp_path, not self.exp_path.text())
            if not self.exp_path.text():
                missing.append("the experiment to sweep")
            used = [i for i in range(2) if any(r["axis"].currentIndex() == i + 1 for r in self._vary_rows)]
            theme.mark_missing(self.vary, self.exp_path.text() != "" and not used)
            if self.exp_path.text() and not used:
                missing.append("at least one row put on axis 1")
            for i in range(2):
                need = i in used
                try:
                    bad = need and not S.parse_values(self.axis_values[i].text())
                except ValueError:
                    bad = True
                theme.mark_missing(self.axis_values[i], bad)
                if bad:
                    missing.append(f"the values of axis {i + 1}")
                elif need:
                    n_axes += 1
                rows = [r for r in self._vary_rows if r["axis"].currentIndex() == i + 1]
                desc = ""
                if rows:
                    hows = {r["how"].currentText() for r in rows}
                    try:
                        vals = S.parse_values(self.axis_values[i].text())
                    except ValueError:
                        vals = []
                    unit = "× the experiment's amount" if hows == {"× factor"} else ("ppb" if hows <= {"= ppb"} else "")
                    shown = ", ".join(f"{v:.3g}" for v in vals[:12]) + (f", … ({len(vals)} values)" if len(vals) > 12 else "")
                    desc = ("varies " + ", ".join(r["what"] for r in rows) +
                            (f"  →  {len(vals)} runs at {shown} {unit}" if vals else "") +
                            ("   (1 = as in the experiment)" if hows == {"× factor"} and vals else ""))
                getattr(self, f"axis_desc{i}").setText(desc)
        else:
            theme.mark_missing(self.exps, self.exps.count() == 0)
            if self.exps.count() == 0:
                missing.append("at least one experiment")

            if self.gen_sweep.isChecked():
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
        if self._mods_active():
            n_mods = sum(1 for r in range(self.mods.rowCount()) if self.mods.item(r, 1) and self.mods.item(r, 1).text().strip())
            if self.mods.rowCount() > n_mods:
                missing.append("the target of a change in the change table")
        if missing:
            self.ready.setText("Still needed: " + ", ".join(missing) + ".")
        else:
            if sweep:
                what = f"a sweep over {n_axes} axis/axes of {os.path.basename(self.exp_path.text())[:-16]}" + \
                       (f", with {n_mods} fixed change(s) on every run" if n_mods else "")
            else:
                parts = ([f"{n_mods} change(s) applied to every experiment"] if n_mods else []) + ([f"a sweep over {n_axes} axis/axes"] if n_axes else [])
                what = (" and ".join(parts) if parts else "the experiments rerun as they are (no change; their chemistry files, this study's solver settings)") + f" — {self.exps.count()} experiment(s)"
            self.ready.setText("Ready: " + what + ". Next: Check, then 'Run pilot, then estimate'.")

    def _count(self, *_a):
        if self._loading:
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

    # ── batch mode: experiments ──
    def _add_exps(self):
        files, _ = QFileDialog.getOpenFileNames(self, "Experiments", self.paths.get("output_root", ""), "Experiment (*.experiment.json)")
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
        self.study = S.new_study("study", "", self.toml_path); self.study_path = ""; self.runner = None
        self.exp_path.clear(); self.vary.setRowCount(0); self._vary_rows = []; self.exp_status.setText("no experiment chosen")
        for e in self.axis_values:
            e.clear()
        self._from_study("sweep")

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
            base = os.path.dirname(st.experiments[0]) if st.experiments else self.paths.get("output_root", HERE)
            st.folder = os.path.join(base, "studies", st.name); self.folder.setText(st.folder)
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
        self._log("This study will run: " + self.ready.text().replace("Ready: ", "").split(". Next")[0])
        for ax in st.axes:
            self._log(f"  axis '{ax.name()}': " + "; ".join(f"{k} {t}" for k, t in ax.members()) + f" over {', '.join(f'{v:g}' for v in ax.values)}")
        for m in st.modifications:
            self._log(f"  change: {m.kind} {m.target} {m.value}")
        for lvl, m in msgs:
            self._log(("✖ " if lvl == "error" else "⚠ ") + m)
        if not msgs:
            self._log("✔ nothing to fix")
        self.tabs.setCurrentWidget(self.msgs)
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
            pilot = todo[0]
            self._events.event.emit("pilot_started", pilot)
            self.runner.run_one(pilot, lambda k, p, _r: self._events.event.emit(k, p))
            self._events.event.emit("pilot_done", pilot)
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
        self.results.setColumnCount(len(cols)); self.results.setHorizontalHeaderLabels([S.METRICS[c][0] if c in S.METRICS else c for c in cols])
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
            rows = self._rows()
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
            ax.set_ylabel(f"{label} [{unit}]"); ax.set_title(f"{label} per experiment (change applied)", loc="left"); ax.grid(axis="y", alpha=0.3)
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
    exps = [argv[i + 1] for i, a in enumerate(argv) if a == "--experiment" and i + 1 < len(argv)]
    mode = argv[argv.index("--mode") + 1] if "--mode" in argv and argv.index("--mode") + 1 < len(argv) else ""
    w = StudyWindow(pick_toml(argv), files[0] if files else "", exps, mode); w.show(); theme.force_taskbar_icon(w)
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
