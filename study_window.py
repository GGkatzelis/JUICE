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
        self._loading = True                            # building fires change signals: they must not read empty widgets back into the study
        self._build()
        self._loading = False
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
        b = QPushButton("+ initial value …"); b.setToolTip("A species' value at the start of the run goes on an axis (ppb)."); b.clicked.connect(self._add_initial_row); hb.addWidget(b)
        b = QPushButton("+ rate × factor …"); b.setToolTip("A reaction's MCM rate times the axis factor."); b.clicked.connect(self._add_rate_row); hb.addWidget(b)
        b = QPushButton("+ new injection …"); b.setToolTip("An injection the experiment does not have (e.g. NO at 11:00 over 5 min); its amount goes on an axis, e.g. 0, 2, 5, 10 ppb.")
        b.clicked.connect(self._add_new_injection_row); hb.addWidget(b)
        b = QPushButton("− row"); b.setToolTip("Removes a row added with the + buttons."); b.clicked.connect(self._remove_custom_row); hb.addWidget(b)
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
        fh = QHBoxLayout(); fh.addWidget(QLabel("Focus"))
        self.focus = QComboBox()
        for k, (label, _keys) in S.FOCUS.items():
            self.focus.addItem(label, k)
        self.focus.addItem("Custom (tick below)", "custom")
        self.focus.currentIndexChanged.connect(self._focus_changed); fh.addWidget(self.focus); fh.addStretch()
        gl.addLayout(fh, 10, 0, 1, 2)
        note = _info("Outputs, not inputs: after each run these numbers are read from its result file and collected in the "
                     "Results table, the CSV and the isopleth / sensitivity plot. They change nothing in the model.")
        note.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Preferred); gl.addWidget(note, 0, 0, 1, 2)
        self.metric_boxes = {}
        for i, (k, (label, unit)) in enumerate(S.METRICS.items()):
            cb = QCheckBox(label + ("" if unit == "-" else f" [{unit}]")); self.metric_boxes[k] = cb; gl.addWidget(cb, 11 + i // 2, i % 2)
            if k in S.METRIC_HELP:
                cb.setToolTip(S.METRIC_HELP[k])
            cb.toggled.connect(self._metric_toggled)
        mh = _more_info("Ozone focus: production rate, efficiency and regime come from each run's HO2, RO2, OH and NOx, so they show "
                        "why O3 changes across an isopleth, not only how much.",
                        [f"<b>{S.METRICS[k][0]}</b>: {v}" for k, v in S.METRIC_HELP.items()], "ⓘ Definitions")
        mh.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Preferred); gl.addWidget(mh, 40, 0, 1, 2)
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
            "<b>Regimes</b> (two axes, one of them NOx): <i>NOx-sensitive</i> = the result rises with NOx (little NOx; radicals end as "
            "peroxides); <i>VOC-sensitive / NOx-saturated</i> = it rises with VOC and falls with more NOx (radicals end as HNO3). The dashed "
            "ridge (dashed) is where ozone production peaks along the NOx axis: for every VOC level a parabola through the three sampled "
            "NOx values around the maximum; where ozone is still rising at the highest NOx sampled the level is NOx-sensitive throughout "
            "and the note says so. The ridge comes from ozone production, so it is the same on every plot of the study, whatever is "
            "coloured. Hatched = VOC-sensitive side. Red dotted: P(H2O2)/P(HNO3) = 0.35, the radical-budget indicator of the same "
            "split (high = NOx-sensitive). Four NOx values give a rough ridge — six or more per axis for a clean one.",
            "<b>log axes</b> for log-spaced sweeps (log 0.1 20 6). <b>Save image</b> writes PNG or SVG."], "ⓘ About this plot"))
        hb = QHBoxLayout(); hb.addWidget(QLabel("Metric")); self.plot_metric = QComboBox(); self.plot_metric.currentIndexChanged.connect(self._plot); hb.addWidget(self.plot_metric)
        self.plot_log = QCheckBox("log axes"); self.plot_log.toggled.connect(self._plot); hb.addWidget(self.plot_log)
        self.show_regimes = QCheckBox("regimes"); self.show_regimes.setChecked(True); self.show_regimes.toggled.connect(self._plot)
        self.show_regimes.setToolTip("Isopleths with a NOx axis: the ridge where the result stops rising with NOx, the NOx-sensitive and\n"
                                     "VOC-sensitive (NOx-saturated) sides, and the 0.35 contour of the regime indicator when it was computed.")
        hb.addWidget(self.show_regimes)
        hb.addWidget(QLabel("Legend"))
        self.legend_pos = QComboBox()
        self.legend_pos.addItems(["where dragged / automatic", "upper left", "upper right", "lower left", "lower right", "below the plot", "hidden"])
        self.legend_pos.setToolTip("Where the legend goes. You can also drag it with the mouse; double-click it to put it back.")
        self.legend_pos.currentIndexChanged.connect(self._legend_pos_changed); hb.addWidget(self.legend_pos)
        b = QPushButton("Save image …"); b.clicked.connect(self._save_image); hb.addWidget(b); hb.addStretch(); pl.addLayout(hb)
        self.fig = Figure(figsize=(8, 5.5)); self.canvas = FigureCanvas(self.fig); pl.addWidget(self.canvas, stretch=1)
        import plot_interaction                       # as in the results viewer: wheel = zoom, drag = pan, double-click = reset,
        plot_interaction.attach(self.canvas, hover=False)   # and the edges of the plot can be dragged to resize it (e.g. the x axis up)
        plot_interaction.attach_dividers(self.canvas)
        self.canvas.setToolTip("Wheel: zoom the x axis - Ctrl+wheel or right-drag: zoom the values - drag: pan - double-click: reset\n"
                               "Drag an edge of the plot (resize cursor) to move it, e.g. the bottom axis further up.")
        self.plot_note = QLabel(""); self.plot_note.setWordWrap(True); self.plot_note.setStyleSheet("color:#555"); pl.addWidget(self.plot_note)
        self.tabs.addTab(pw, "Plot")
        self.tabs.addTab(self._build_video_tab(), "Isopleth in time")
        self.msgs = QTextEdit(); self.msgs.setReadOnly(True); self.tabs.addTab(self.msgs, "Messages")

        split = QSplitter(); split.addWidget(scroll); split.addWidget(self.tabs); split.setStretchFactor(0, 65); split.setStretchFactor(1, 35)
        split.setSizes([1066, 574]); self.setCentralWidget(split)
        self._mode_changed()

    def _focus_changed(self, *_a):
        key = self.focus.currentData()
        self._plot_first = True                          # the plot shows the focus's first result
        if key in S.FOCUS and not self._loading:
            self._loading = True
            try:
                for k, cb in self.metric_boxes.items():
                    cb.setChecked(k in S.FOCUS[key][1])
            finally:
                self._loading = False
            self._count()

    def _metric_toggled(self, *_a):
        if self._loading:
            return
        ticked = {k for k, cb in self.metric_boxes.items() if cb.isChecked()}
        match = next((k for k, (_l, keys) in S.FOCUS.items() if set(keys) == ticked), "custom")
        self.focus.blockSignals(True); self.focus.setCurrentIndex(self.focus.findData(match)); self.focus.blockSignals(False)

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
            match = next((k for k, (_l, keys) in S.FOCUS.items() if set(keys) == set(metrics)), "custom")
            self.focus.setCurrentIndex(self.focus.findData(match))
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
        focus = self.focus.currentData()
        order = S.FOCUS[focus][1] if focus in S.FOCUS else []
        ticked = [k for k, cb in self.metric_boxes.items() if cb.isChecked()]
        st.metrics = [k for k in order if k in ticked] + [k for k in ticked if k not in order]     # the focus order: first = plotted first
        st.focus = focus
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
            done = sum(1 for p in S.plan_runs(st) if S.is_done(p)) if st.folder and os.path.isdir(os.path.join(st.folder, "runs")) else 0
            todo = n - done
            if done:
                txt += f" — {done} already done (kept), {todo} to run" if todo else f" — all {done} already done"
            prior = S.prior_seconds(st)
            if prior and todo:
                txt += f";  last run time of the base experiment(s) ≈ {S.fmt_duration(prior)} → about {S.fmt_duration(S.estimate_seconds(todo, prior, st.workers))} with {st.workers} workers (before a pilot)"
            lvl, g = S.check_guards(todo or n)
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
            wk = self.workers.value(); fac = S.parallel_factor(min(wk, len(rest)))
            best = min(range(1, wk + 1), key=lambda k: S.estimate_seconds(len(rest), p.elapsed, k))
            msg = (f"Pilot {p.name}: {S.fmt_duration(p.elapsed)} alone.\n\n{len(rest)} runs remaining with {wk} in parallel "
                   f"(each about {fac:.1f}× slower when sharing the machine) ≈ {txt}" +
                   (f"\n{best} in parallel would take about as long or less ({S.fmt_duration(S.estimate_seconds(len(rest), p.elapsed, best))})"
                    if best < wk else "") + (f"; about {disk / 1e6:.0f} MB of results" if disk and not math.isnan(disk) else "") +
                   (f"\n\n⚠ {g}" if lvl != "ok" else "") + "\n\nRun the remaining runs now?")
            if QMessageBox.question(self, "Estimate", msg) == QMessageBox.StandardButton.Yes:
                self._busy(True); self.runner.start(lambda k, pp, _r: self._events.event.emit(k, pp), self.workers.value(), rest)
            else:
                self.status.setText(f"stopped after the pilot ({txt} estimated for the rest); 'Run all' continues")
            return
        if kind == "all_done":
            self._busy(False); self._refresh_results()
            if self.tabs.currentWidget() is getattr(self, "_v_tab", None):
                self._v_refresh_if_stale()
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
        rows = S.collect_results(st, plans)
        need = [k for k in ("regime", "po3_lit_int") if k not in st.metrics]
        if need and len(st.axes) >= 2:          # the isopleth's ridge (ozone) and indicator contour need them, even when not ticked
            extra = {r["run"]: r for r in S.collect_results(S.Study(**{**st.__dict__, "metrics": need}), plans)}
            for r in rows:
                for k in need:
                    r[k] = extra.get(r["run"], {}).get(k)
        return rows

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
        i = -1 if getattr(self, "_plot_first", False) else self.plot_metric.findData(cur); self._plot_first = False
        self.plot_metric.setCurrentIndex(i if i >= 0 else 0); self.plot_metric.blockSignals(False)
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
            dropped = 0
            if log:                                          # no logarithm of 0: those runs cannot be placed on log axes
                keep = [r for r in sub if r[a1] > 0 and r[a2] > 0]; dropped = len(sub) - len(keep); sub = keep
            x = np.array([r[a1] for r in sub]); y = np.array([r[a2] for r in sub]); z = np.array([r[key] for r in sub])
            if len(sub) >= 4:
                xs, ys = (np.log10(x), np.log10(y)) if log else (x, y)
                cs = None
                try:
                    cs = ax.tricontourf(xs, ys, z, levels=14, cmap="viridis"); self.fig.colorbar(cs, ax=ax, label=f"{label} [{unit}]")
                    ax.tricontour(xs, ys, z, levels=7, colors="white", linewidths=0.6, alpha=0.7)
                except Exception:
                    sc = ax.scatter(xs, ys, c=z, cmap="viridis", s=60); self.fig.colorbar(sc, ax=ax, label=f"{label} [{unit}]")
                self._run_markers(ax, xs, ys, z, cs)
                for xi, yi, zi in zip(xs, ys, z):
                    ax.annotate(f"{zi:.3g}", (xi, yi), fontsize=F["small"], color="white", xytext=(3, 3), textcoords="offset points")
                regime_note = self._regimes(ax, st, sub, a1, a2, xs, ys, z, key, log) if self.show_regimes.isChecked() else ""
            else:
                sc = ax.scatter(x, y, c=z, cmap="viridis", s=80); self.fig.colorbar(sc, ax=ax, label=f"{label} [{unit}]")
            (lab1, tot1, _s1), (lab2, tot2, _s2) = S.axis_description(st, axes[0]), S.axis_description(st, axes[1])
            ax.set_xlabel(lab1 + (" — log10 scale" if log else ""), fontsize=F["label"] - 1)
            ax.set_ylabel(lab2 + (" — log10 scale" if log else ""), fontsize=F["label"] - 1)
            self._ppb_scale(ax, "x", tot1, log); self._ppb_scale(ax, "y", tot2, log)
            ax.set_title(f"{label} — {exp}" + (f"  (+{len(exps) - 1} more experiments: see Results)" if len(exps) > 1 else ""), loc="left")
            self.plot_note.setText(self._meaning(key, "Colour") + f"\nisopleth over {a1} × {a2}; {len(sub)} grid points, interpolated between them"
                                   + (f"; {dropped} run(s) at 0 left out — no logarithm of 0 (untick 'log axes' to see them)" if dropped else "")
                                   + (("\n" + regime_note) if len(sub) >= 4 and regime_note else ""))
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
            lab1, tot1, _s1 = S.axis_description(st, axes[0])
            ax.set_xlabel(lab1, fontsize=F["label"] - 1); ax.set_ylabel(f"{label} [{unit}]"); ax.grid(alpha=0.3); ax.legend(fontsize=F["legend"]); ax.set_title(f"{label} vs {a1}", loc="left")
            if tot1 and not log:
                self._ppb_scale(ax, "x", tot1, False)
            self.plot_note.setText("\n".join([self._meaning(key, "Lines")] + notes))
        else:
            ys = [next((r[key] for r in rows if r["experiment"] == e), float("nan")) for e in exps]
            ax.bar(range(len(exps)), ys, color="#6a4c93"); ax.set_xticks(range(len(exps))); ax.set_xticklabels(exps, rotation=20, ha="right", fontsize=F["tick"])
            ax.set_ylabel(f"{label} [{unit}]"); ax.set_title(f"{label} per experiment (change applied)", loc="left"); ax.grid(axis="y", alpha=0.3)
            self.plot_note.setText(self._meaning(key, "Bars") + "\none value per experiment; compare with the base runs in the viewer's Comparison tab")
        self._apply_legend_pos(self.fig)
        self.fig.tight_layout(); self.canvas.draw()

    # ── isopleth in time (video) ──
    def _build_video_tab(self) -> QWidget:
        vw = QWidget(); vl = QVBoxLayout(vw); vl.setContentsMargins(0, 0, 0, 0)
        vl.addWidget(_more_info("The isopleth at every moment of the day: one frame per time step, the same colour scale for the whole "
                                "day, so the evolution through the injections and the lit period is visible.", [
            "Needs two axes. The first opening reads every run's time series once and keeps them in the study folder "
            "(&lt;study&gt;.timeseries.npz); after that frames are instant.",
            "<b>P(O3) instantaneous</b> = k(HO2+NO)[HO2][NO] + k(RO2+NO)[RO2][NO] at that time; <b>P(H2O2)/P(HNO3) instantaneous</b> "
            "= the regime indicator at that time (high = NOx-sensitive, ~0.35 transition).",
            "<b>regimes</b>: the ozone ridge of that moment, from the instantaneous P(O3) along the NOx axis (dashed), VOC-sensitive "
            "side hatched — it moves through the day as NOx is consumed and radicals change.",
            "<b>Export video</b> writes MP4 (or GIF) of the frames between the From/To times at the chosen frames per second."],
            "ⓘ About the video"))
        hb = QHBoxLayout()
        hb.addWidget(QLabel("Quantity")); self.v_qty = QComboBox()
        for k, (lab, unit) in S.TS_QUANTITIES.items():
            self.v_qty.addItem(f"{lab} [{unit}]" if unit != "-" else lab, k)
        self.v_qty.setCurrentIndex(self.v_qty.findData("O3")); self.v_qty.currentIndexChanged.connect(lambda *_: self._v_draw(rescale=True)); hb.addWidget(self.v_qty)
        self.v_log = QCheckBox("log axes"); self.v_log.toggled.connect(lambda *_: self._v_draw()); hb.addWidget(self.v_log)
        self.v_fixed = QCheckBox("one colour scale for the day"); self.v_fixed.setChecked(True); self.v_fixed.toggled.connect(lambda *_: self._v_draw()); hb.addWidget(self.v_fixed)
        self.v_reg = QCheckBox("regimes"); self.v_reg.setChecked(True); self.v_reg.toggled.connect(lambda *_: self._v_draw()); hb.addWidget(self.v_reg)
        hb.addWidget(QLabel("Legend")); self.v_legend = QComboBox()
        self.v_legend.addItems(["where dragged / automatic", "upper left", "upper right", "lower left", "lower right", "below the plot", "hidden"])
        self.v_legend.currentIndexChanged.connect(lambda *_: self._v_draw()); hb.addWidget(self.v_legend)
        hb.addStretch(); vl.addLayout(hb)
        hb = QHBoxLayout()
        self.v_play = QPushButton("▶ Play"); self.v_play.setCheckable(True); self.v_play.toggled.connect(self._v_toggle_play); hb.addWidget(self.v_play)
        from PyQt6.QtWidgets import QSlider
        self.v_slider = QSlider(Qt.Orientation.Horizontal); self.v_slider.valueChanged.connect(lambda *_: self._v_draw()); hb.addWidget(self.v_slider, stretch=1)
        self.v_time = QLabel("--:--"); self.v_time.setMinimumWidth(110); hb.addWidget(self.v_time)
        hb.addWidget(QLabel("frames/s")); self.v_fps = QSpinBox(); self.v_fps.setRange(1, 30); self.v_fps.setValue(6); hb.addWidget(self.v_fps)
        hb.addWidget(QLabel("From")); self.v_from = QLineEdit(); self.v_from.setPlaceholderText("HH:MM"); self.v_from.setMaximumWidth(70); hb.addWidget(self.v_from)
        hb.addWidget(QLabel("To")); self.v_to = QLineEdit(); self.v_to.setPlaceholderText("HH:MM"); self.v_to.setMaximumWidth(70); hb.addWidget(self.v_to)
        b = QPushButton("Export video …"); b.clicked.connect(self._v_export); hb.addWidget(b)
        b = QPushButton("Load / refresh"); b.setToolTip("Reads the time series of the finished runs (once; then from the study folder)."); b.clicked.connect(lambda: self._v_load(force=True)); hb.addWidget(b)
        vl.addLayout(hb)
        self.v_fig = Figure(figsize=(8, 5.5)); self.v_canvas = FigureCanvas(self.v_fig); vl.addWidget(self.v_canvas, stretch=1)
        import plot_interaction
        plot_interaction.attach(self.v_canvas, hover=False); plot_interaction.attach_dividers(self.v_canvas)
        self.v_note = QLabel(""); self.v_note.setWordWrap(True); self.v_note.setStyleSheet("color:#555"); vl.addWidget(self.v_note)
        self._v_data = None; self._v_scale = {}
        from PyQt6.QtCore import QTimer
        self._v_timer = QTimer(self); self._v_timer.timeout.connect(self._v_step)
        self._v_tab = vw
        self.tabs.currentChanged.connect(lambda i: self._v_refresh_if_stale() if self.tabs.widget(i) is vw else None)
        return vw

    def _v_refresh_if_stale(self):
        """Re-read when the finished runs differ from what the video holds (runs added or finished since)."""
        st = self._to_study()
        plans = self.runner.plans if self.runner else S.plan_runs(st)
        done = sorted(p.name for p in plans if S.is_done(p))
        if self._v_data is None or done != sorted(self._v_data.get("runs", [])):
            self._v_load()

    def _v_load(self, force: bool = False):
        st = self._to_study()
        if len(st.axes) < 2:
            self.v_note.setText("The video needs a study with two axes (an isopleth)."); return
        plans = self.runner.plans if self.runner else S.plan_runs(st)
        from PyQt6.QtWidgets import QProgressDialog
        dlg = QProgressDialog("Reading the runs' time series …", "Cancel", 0, len(plans), self); dlg.setMinimumDuration(300)

        def prog(i, n, name):
            dlg.setMaximum(n); dlg.setValue(i); dlg.setLabelText(f"Reading {i + 1}/{n}: {name}"); QApplication.processEvents()
        if force:
            try:
                os.remove(os.path.join(st.folder, f"{st.name}.timeseries.npz"))
            except OSError:
                pass
        self._v_data = S.time_series(st, plans, 5.0, prog); dlg.close()
        if not self._v_data:
            self.v_note.setText("No finished runs yet."); return
        n = len(self._v_data["hours"]); self.v_slider.setRange(0, n - 1)
        lit = np.where(self._v_data["hv"] > 0.5)[0]
        self.v_slider.setValue(int(lit[len(lit) // 2]) if len(lit) else n // 2)
        self._v_scale = {}; self._v_draw(rescale=True)

    def _v_frame_values(self, key: str, i: int):
        d = self._v_data; st = self.study
        a1, a2 = st.axes[0].name(), st.axes[1].name()
        x, y, z = d["axes"][a1], d["axes"][a2], d[key][:, i]
        ok = np.isfinite(z) & np.isfinite(x) & np.isfinite(y)
        if self.v_log.isChecked():
            ok &= (x > 0) & (y > 0)
        return a1, a2, x[ok], y[ok], z[ok], d["PO3"][ok, i], int((~ok).sum())

    def _v_draw(self, rescale: bool = False, fig=None):
        if not self._v_data:
            return
        fig = fig or self.v_fig; i = self.v_slider.value(); key = self.v_qty.currentData(); d = self._v_data; st = self.study
        log = self.v_log.isChecked()
        a1, a2, x, y, z, po3, dropped = self._v_frame_values(key, i)
        hh = d["hours"][i]; lit = d["hv"][i] > 0.5
        self.v_time.setText(f"{int(hh):02d}:{int(round(hh % 1 * 60)) % 60:02d} UTC" + (" ☀" if lit else ""))
        fig.clear(); ax = fig.add_subplot(111)
        lab, unit = S.TS_QUANTITIES[key]
        if len(z) < 4:
            ax.text(0.5, 0.5, "fewer than 4 runs with values at this time", ha="center", transform=ax.transAxes)
            self.v_canvas.draw(); return
        xs, ys = (np.log10(x), np.log10(y)) if log else (x, y)
        levels = 14
        if self.v_fixed.isChecked():
            if key not in self._v_scale or rescale:
                allv = d[key][np.isfinite(d[key])]
                lo, hi = (np.nanpercentile(allv, 1), np.nanpercentile(allv, 99)) if allv.size else (0, 1)
                self._v_scale[key] = (lo, hi if hi > lo else lo + 1)
            lo, hi = self._v_scale[key]; levels = np.linspace(lo, hi, 15)
        cs = None
        try:                                   # values beyond the day's scale: end colours (arrows on the colour bar), never holes
            cs = ax.tricontourf(xs, ys, z, levels=levels, cmap="viridis", extend="both" if not isinstance(levels, int) else "neither")
            fig.colorbar(cs, ax=ax, label=f"{lab} [{unit}]" if unit != "-" else lab)
        except Exception:
            sc = ax.scatter(xs, ys, c=z, cmap="viridis", s=60); fig.colorbar(sc, ax=ax)
        self._run_markers(ax, xs, ys, z, cs)
        note = ""
        if self.v_reg.isChecked():
            note = self._v_ridge(ax, st, a1, a2, x, y, po3, log)
        (l1, t1, _), (l2, t2, _) = S.axis_description(st, st.axes[0]), S.axis_description(st, st.axes[1])
        ax.set_xlabel(l1 + (" — log10 scale" if log else ""), fontsize=F["label"] - 1); ax.set_ylabel(l2 + (" — log10 scale" if log else ""), fontsize=F["label"] - 1)
        self._ppb_scale(ax, "x", t1, log); self._ppb_scale(ax, "y", t2, log)
        ax.set_title(f"{lab} at {int(hh):02d}:{int(round(hh % 1 * 60)) % 60:02d} UTC" + ("  (roof open)" if lit else "  (dark)"), loc="left")
        self._apply_legend_pos(fig, self.v_legend)
        fig.tight_layout()
        if fig is self.v_fig:
            self.v_canvas.draw()
            self.v_note.setText(self._v_meaning(key) + (f"\n{note}" if note else "") +
                                (f"\n{dropped} run(s) at 0 left out on log axes" if dropped else ""))

    def _v_meaning(self, key: str) -> str:
        return {"PO3": "Colour = gross ozone production at this moment (k[HO2][NO] + k[RO2][NO]); it peaks where NO meets enough peroxy radicals.",
                "REGIME": "Colour = regime indicator at this moment: high = NOx-sensitive, low = VOC-sensitive, transition ~0.35.",
                "O3": "Colour = O3 at this moment: the result of production so far, titration by NO and dilution.",
                "Ox": "Colour = Ox = O3 + NO2 at this moment: ozone made so far without the NO titration effect."}.get(
            key, f"Colour = {S.TS_QUANTITIES[key][0]} at this moment.")

    def _v_ridge(self, ax, st, a1, a2, x, y, po3, log) -> str:
        i_nox = S.nox_axis(st.axes[:2])
        if i_nox < 0:
            return "regimes: no axis varies NOx"
        if np.nanmax(po3) <= 0:
            return "regimes: no ozone production at this time (dark or no NO)"
        rows = [{a1: xi, a2: yi, "po3": pi} for xi, yi, pi in zip(x, y, po3)]
        a_nox, a_oth = (a1, a2) if i_nox == 0 else (a2, a1)
        ridge = S.ozone_ridge(rows, a_nox, a_oth, "po3", log)
        tr = (lambda v: np.log10(v)) if log else (lambda v: v)
        nmax = tr(max(r[a_nox] for r in rows)); nmin = tr(min(r[a_nox] for r in rows))
        pts = [(tr(d["level"]), tr(d["ridge"]) if d["status"] == "inside" else (nmax if d["status"] == "above" else nmin), d["status"]) for d in ridge]
        if not pts:
            return ""
        lv = np.array([p[0] for p in pts]); xr = np.array([p[1] for p in pts])
        if any(p[2] == "inside" for p in pts):
            ax.plot(*((xr, lv) if i_nox == 0 else (lv, xr)), "k--", lw=1.6, label="ozone ridge (now)")
        if np.any(xr < nmax - 1e-9):
            (ax.fill_betweenx if i_nox == 0 else ax.fill_between)(lv, xr, nmax, facecolor="none", hatch="///", edgecolor="k", alpha=0.35,
                                                                  linewidth=0, label="VOC-sensitive now")
        if ax.get_legend_handles_labels()[0]:
            ax.legend(fontsize=F["small"], loc="upper left", framealpha=0.85)
        inside = [d for d in ridge if d["status"] == "inside"]
        return ("ridge of the instantaneous P(O3): " + "; ".join(f"{a_oth} {d['level']:g} → {d['ridge']:.3g}" for d in inside)) if inside \
            else "no ozone ridge inside the sampled NOx range at this time (NOx-sensitive everywhere)" if all(d["status"] == "above" for d in ridge) else ""

    def _v_toggle_play(self, on: bool):
        self.v_play.setText("❚❚ Pause" if on else "▶ Play")
        if on and self._v_data:
            self._v_timer.start(int(1000 / self.v_fps.value()))
        else:
            self._v_timer.stop()

    def _v_step(self):
        v = self.v_slider.value() + 1
        self.v_slider.setValue(v if v <= self.v_slider.maximum() else 0)

    def _v_range(self):
        h = self._v_data["hours"]

        def idx(txt, default):
            try:
                hh, mm = txt.strip().split(":"); return int(np.argmin(abs(h - (int(hh) + int(mm) / 60.0))))
            except Exception:
                return default
        return idx(self.v_from.text(), 0), idx(self.v_to.text(), len(h) - 1)

    def _v_export(self):
        if not self._v_data:
            self._v_load()
            if not self._v_data:
                return
        try:
            import imageio.v2 as imageio
        except ImportError:
            QMessageBox.warning(self, "Video", "Saving a video needs the imageio package (pip install imageio imageio-ffmpeg)."); return
        key = self.v_qty.currentData()
        p, _ = QFileDialog.getSaveFileName(self, "Export video", os.path.join(self.study.folder, f"{self.study.name}_{key}_isopleth.mp4"),
                                           "MP4 video (*.mp4);;GIF animation (*.gif)")
        if not p:
            return
        i0, i1 = self._v_range()
        frames = list(range(min(i0, i1), max(i0, i1) + 1))
        fps = self.v_fps.value()
        try:
            writer = imageio.get_writer(p, mode="I", duration=1.0 / fps, loop=0) if p.lower().endswith(".gif") else \
                imageio.get_writer(p, fps=fps, codec="libx264", quality=8, pixelformat="yuv420p", macro_block_size=2)
        except Exception as e:
            QMessageBox.critical(self, "Video", f"Could not open {p}: {e}"); return
        from PyQt6.QtWidgets import QProgressDialog
        from matplotlib.backends.backend_agg import FigureCanvasAgg
        dlg = QProgressDialog("Writing frames …", "Cancel", 0, len(frames), self); dlg.setMinimumDuration(0)
        fig = Figure(figsize=(10, 6.6), dpi=110); FigureCanvasAgg(fig)
        keep = self.v_slider.value()
        try:
            for n, i in enumerate(frames):
                if dlg.wasCanceled():
                    break
                self.v_slider.blockSignals(True); self.v_slider.setValue(i); self.v_slider.blockSignals(False)
                self._v_draw(fig=fig); fig.canvas.draw()
                frame = np.asarray(fig.canvas.buffer_rgba())[:, :, :3]
                writer.append_data(frame.copy()); dlg.setValue(n + 1); QApplication.processEvents()
        finally:
            writer.close(); dlg.close()
            self.v_slider.setValue(keep); self._v_draw()
        self.statusBar().showMessage(f"video written: {p} ({len(frames)} frames at {fps} frames/s)", 10000)

    def _legend_pos_changed(self, *_a):
        lg = getattr(self.canvas, "_juice_legends", None)
        if lg is not None and self.legend_pos.currentIndex() > 0:
            lg.ops.clear(); lg._orig.clear()          # a chosen position replaces a dragged one
        self._plot()

    def _apply_legend_pos(self, fig, combo=None):
        combo = combo or self.legend_pos
        choice = combo.currentText()
        if combo.currentIndex() == 0:
            return
        for ax in fig.axes:
            lg = ax.get_legend()
            if lg is None:
                continue
            if choice == "hidden":
                lg.set_visible(False); continue
            if choice == "below the plot":
                lg.set_bbox_to_anchor((0.5, -0.16), transform=ax.transAxes); lg.set_loc("upper center")
            else:
                lg.set_bbox_to_anchor(None); lg.set_loc(choice)

    def _run_markers(self, ax, xs, ys, z, cs):
        """One circle per run, filled with that run's own value on the contour's colour scale, thin white edge:
        a circle that stands out from its surroundings shows where the contour between runs is interpolated."""
        if cs is not None:
            ax.scatter(xs, ys, c=z, cmap=cs.cmap, norm=cs.norm, s=34, edgecolors="white", linewidths=0.8, zorder=5)
        else:
            ax.scatter(xs, ys, c="k", s=8, zorder=5)

    def _ppb_scale(self, ax, which: str, total_ppb: float, log: bool):
        """Second scale on the opposite side of a factor axis: the total injected ppb that factor means."""
        if not total_ppb:
            return
        if log:
            fwd = lambda v: np.power(10.0, v) * total_ppb
            inv = lambda p: np.log10(np.clip(p, 1e-12, None) / total_ppb)
        else:
            fwd = lambda v: v * total_ppb
            inv = lambda p: p / total_ppb
        try:
            if which == "x":
                sec = ax.secondary_xaxis("top", functions=(fwd, inv)); sec.set_xlabel("total injected [ppb]", fontsize=F["small"] + 1)
            else:
                sec = ax.secondary_yaxis("right", functions=(fwd, inv)); sec.set_ylabel("total injected [ppb]", fontsize=F["small"] + 1)
            sec.tick_params(labelsize=F["small"])
        except Exception:
            pass

    def _meaning(self, key: str, what: str) -> str:
        label, unit = S.METRICS.get(key, (key, ""))
        is_, expect = S.METRIC_MEANING.get(key, ("", ""))
        return f"{what} = {label}" + (f" [{unit}]" if unit and unit != "-" else "") + (f": {is_}." if is_ else ".") + \
               (f"  What to expect: {expect}" if expect else "")

    def _regimes(self, ax, st, sub, a1, a2, xs, ys, z, key, log) -> str:
        """Ozone ridge, VOC-sensitive side (hatched) and the indicator contour on a two-axis isopleth with a NOx axis.
        The ridge comes from ozone production (not from the coloured result), so it is the same on every plot of the study."""
        from matplotlib.patches import Patch
        from matplotlib.lines import Line2D
        i_nox = S.nox_axis(st.axes[:2])
        if i_nox < 0:
            return "regimes: no axis varies NOx (NO/NO2 injection or initial value), so the regimes are not drawn"
        okey = next((k for k in ("po3_lit_int", "ox_lit_delta", "o3_max") if all(r.get(k) is not None and np.isfinite(r[k]) for r in sub)), None)
        if okey is None:
            return "regimes: no ozone result available for the runs (P(O3), Ox gain or O3 maximum)"
        a_nox, a_oth = (a1, a2) if i_nox == 0 else (a2, a1)
        ridge = S.ozone_ridge(sub, a_nox, a_oth, okey, log)
        if not ridge:
            return "regimes: too few grid points"
        tr = (lambda v: np.log10(v)) if log else (lambda v: v)
        nox_vals = sorted({r[a_nox] for r in sub}); nmin, nmax = tr(nox_vals[0]), tr(nox_vals[-1])
        handles = []
        # ridge through the per-level maxima (interior ones); levels with the maximum at an edge extend the line to that edge
        pts = []
        for d in ridge:
            xr = tr(d["ridge"]) if d["status"] == "inside" else (nmax if d["status"] == "above" else nmin)
            pts.append((tr(d["level"]), xr, d["status"]))
        inside = [p for p in pts if p[2] == "inside"]
        if inside:
            lv = [p[0] for p in pts]; xr = [p[1] for p in pts]
            line = (xr, lv) if i_nox == 0 else (lv, xr)
            ax.plot(*line, "k--", lw=1.8)
            handles.append(Line2D([], [], color="k", ls="--", lw=1.8, label="ozone ridge: NOx of maximum ozone production"))
        # VOC-sensitive side = between the ridge and the highest NOx sampled
        lv = np.array([p[0] for p in pts]); xr = np.array([p[1] for p in pts])
        if np.any(xr < nmax - 1e-9):
            if i_nox == 0:
                ax.fill_betweenx(lv, xr, nmax, facecolor="none", hatch="///", edgecolor="k", alpha=0.35, linewidth=0)
            else:
                ax.fill_between(lv, xr, nmax, facecolor="none", hatch="///", edgecolor="k", alpha=0.35, linewidth=0)
            handles.append(Patch(facecolor="none", hatch="///", edgecolor="k", alpha=0.5, label="VOC-sensitive (NOx-saturated): ozone falls with more NOx"))
        handles.append(Patch(facecolor="none", edgecolor="none", label="elsewhere NOx-sensitive: ozone rises with NOx"))
        # chemical indicator contour, with its own legend entry
        rows_r = [r for r in sub if r.get("regime") is not None and np.isfinite(r["regime"])]
        ind = ""
        if len(rows_r) >= 4:
            xr_ = np.array([tr(r[a1]) for r in rows_r]); yr_ = np.array([tr(r[a2]) for r in rows_r]); zr = np.array([r["regime"] for r in rows_r], float)
            if np.nanmin(zr) < 0.35 < np.nanmax(zr):
                try:
                    ax.tricontour(xr_, yr_, zr, levels=[0.35], colors="#d62728", linewidths=1.3, linestyles=":")
                    handles.append(Line2D([], [], color="#d62728", ls=":", lw=1.3, label="P(H2O2)/P(HNO3) = 0.35 (radical-budget indicator)"))
                    ind = "; red dotted = indicator 0.35"
                except Exception:
                    pass
            else:
                ind = f"; indicator {np.nanmin(zr):.2g}–{np.nanmax(zr):.2g} everywhere, no 0.35 contour"
        ax.legend(handles=handles, fontsize=F["small"], loc="upper left", framealpha=0.85)
        above = [f"{d['level']:g}" for d in ridge if d["status"] == "above"]; below = [f"{d['level']:g}" for d in ridge if d["status"] == "below"]
        fmt = lambda d: f"{a_oth} = {d['level']:g}: {d['ridge']:.3g}"
        note = (f"regimes from {S.METRICS[okey][0]} along {a_nox} (same on every plot of this study); ridge at " +
                ("; ".join(fmt(d) for d in ridge if d["status"] == "inside") or "no interior maximum") +
                (f"; still rising at the highest {a_nox} sampled for {a_oth} = {', '.join(above)} (NOx-sensitive throughout)" if above else "") +
                (f"; falling from the lowest {a_nox} on for {a_oth} = {', '.join(below)} (VOC-sensitive throughout)" if below else "") + ind +
                (f"; only {ridge[0]['n']} NOx values per level — the ridge is a parabola through the three points around the maximum, extend the axis or add points to sharpen it" if ridge[0]['n'] <= 4 else ""))
        return note

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
