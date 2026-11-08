"""
field_window.py — JUICE for field observations, step by step, in the layout of the chamber window.

One window, five steps down the left, the active step in the middle, a review panel on the right
that always shows what the case contains and what still needs attention.

    1 Field case                 kind of run, name, site, date, input table — with the conditions
    2 Chemistry                  the chamber window's Chemistry page: MCM subset or existing mechanism
    3 Emissions and constraints  what enters the air or is held to data — with their time series
    4 Mixing and deposition      dilution, background, deposition, initial values — with loss per hour
    5 Run                        write the input files (field_scenario), run the solver, open the results

The only thing a user edits is the case (<name>.field.json, field_spec.FieldSpec).  The input files
(ENZ, EASY, chemistry, scenario JSON) are generated from it in <project folder>/<name>/.

    python field_window.py [case.field.json] [-t file.toml]
"""
from __future__ import annotations

import io
import contextlib
import json
import os
import re
import sys
import time
from typing import Dict, List, Optional, Tuple

import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from PyQt6.QtCore import Qt, QProcess, QProcessEnvironment, QTimer
from PyQt6.QtGui import QKeySequence, QPixmap
from PyQt6.QtWidgets import (QApplication, QButtonGroup, QCheckBox, QComboBox, QDoubleSpinBox, QFileDialog, QFormLayout,
                             QGridLayout, QGroupBox, QHBoxLayout, QHeaderView, QLabel, QLineEdit, QListWidget, QListWidgetItem,
                             QMainWindow, QMessageBox, QProgressBar, QPushButton, QRadioButton, QSpinBox, QSplitter,
                             QStackedWidget, QTableWidget, QTableWidgetItem, QTabWidget, QTextEdit, QVBoxLayout, QWidget)

import field_scenario as fsc
import field_spec as fs
import juice_theme as theme
from experiment_spec import ExperimentSpec
from experiment_wizard import ChemistryPage, _add_row, _del_row, _info, _rows, _set_rows, _table
from juice_gui import _boxed, _open_path, _scroll
from scenario_prep import default_profile, load_paths, pick_toml, solver_from_toml, solver_profiles

STEPS = ["Field case", "Chemistry", "Emissions and constraints", "Mixing and deposition", "Run"]
RUN_STEP = fs.STEP_RUN
CASE_FILTER = "Field case (*.field.json)"
TABLE_FILTER = "Input table (*.csv *.txt *.tsv *.xlsx *.xls);;All files (*)"
CHEM_KEYS = ("base_mechanism", "addons", "precursors", "mechanism_stem", "inorganic_only", "include_dirs", "reaction_adjustments")
RUN_TEXT = {
    "constrained": "A measurement site. Measured species are held to the data; radicals and products run free. "
                   "No emissions or transport are needed.",
    "plume": "A plume or transect. Chemistry ages while the air is diluted towards background.",
    "lagrangian": "An air parcel along a trajectory, with emissions, a growing boundary layer and background air.",
}
NEEDS = {
    "constrained": ["T, p and H2O or RH", "one C_<species> column per measured species you hold",
                    "jNO2 if you have it, and the position of the site"],
    "plume": ["T, p and H2O or RH", "a dilution rate: the column kdil, or C_ of a tracer such as CO",
              "B_<species> columns for the air it mixes with", "initial values on the Mixing and deposition step"],
    "lagrangian": ["T, p and H2O or RH, and lat and lon per row", "the boundary-layer height column h",
                   "one E_<species> column per emitted species", "B_<species> columns for the air it mixes with"],
}
CONDITION_UNITS = {"t": "K", "p": "hPa", "H2O": "vol-%", "RH": "%", "sz": "deg", "lat": "deg", "lon": "deg", "h": "m", "kdil": "1/s"}
COLORS = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#000000", "#F0E442"]   # fixed order, colour-blind safe


def _span(v: np.ndarray) -> Tuple[str, str]:
    finite = v[np.isfinite(v)]
    return (f"{finite.min():.4g}", f"{finite.max():.4g}") if len(finite) else ("-", "-")


def preview_rows(t: fs.FieldTable, spec: fs.FieldSpec) -> List[Tuple[str, str, str, str, str]]:
    """What the table was understood to contain: role, name, unit, min, max."""
    rows = [("condition", k, CONDITION_UNITS.get(k, ""), *_span(v)) for k, v in t.conditions.items()]
    rows += [("photolysis", f"j{k}", "1/s", *_span(v)) for k, v in t.photolysis.items()]
    rows += [("emission", k, fs.resolve_emission_unit(t, spec, k), *_span(v)) for k, v in t.emissions.items()]
    rows += [("background", k, "ppb", *_span(v)) for k, v in t.background.items()]
    rows += [("held to data", k, "ppb", *_span(v)) for k, v in t.constraints.items()]
    return rows


def _readonly_table(headers: List[str]) -> QTableWidget:
    t = QTableWidget(0, len(headers)); t.setHorizontalHeaderLabels(headers)
    t.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch); t.verticalHeader().setVisible(False)
    t.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers); t.setMinimumHeight(110)
    return t


def _fill(t: QTableWidget, rows: List[Tuple[str, ...]]) -> None:
    t.setRowCount(len(rows))
    for i, row in enumerate(rows):
        for j, text in enumerate(row):
            t.setItem(i, j, QTableWidgetItem(text))


class FieldPlot(QWidget):
    """The inputs of one step as the model will see them: 'conditions', 'sources' or 'mixing'."""

    NOTES = {"conditions": "Conditions of the input table against time.",
             "sources": "Emissions as a volume source, and the species held to the data.",
             "mixing": "First-order loss per hour from dilution and deposition."}

    def __init__(self, win: "FieldWindow", mode: str, height: float = 2.8):
        super().__init__()
        self.win, self.mode, self.plotted = win, mode, False
        lay = QVBoxLayout(self); lay.setContentsMargins(0, 0, 0, 0)
        self.fig = Figure(figsize=(6, height), tight_layout=True); self.canvas = FigureCanvas(self.fig)
        self.canvas.setMinimumHeight(int(height * 95)); lay.addWidget(self.canvas, stretch=1)
        try:
            from plot_interaction import attach
            attach(self.canvas)
        except Exception as e:
            print(f"[field] plot interaction not attached: {e}")
        self.note = _info(self.NOTES[mode]); lay.addWidget(self.note)

    def _message(self, text: str) -> None:
        ax = self.fig.add_subplot(111); ax.text(0.5, 0.5, text, ha="center", va="center", wrap=True, fontsize=9); ax.axis("off")
        self.canvas.draw(); self.plotted = True

    def plot(self) -> None:
        self.fig.clear()
        spec, table, rep = self.win.case
        if table is None:
            self._message("nothing to plot yet: choose an input table on the Field case step")
            return
        step = {"sources": fs.STEP_SOURCES, "mixing": fs.STEP_MIXING}.get(self.mode)
        if step is not None and rep.worst(step) == "error":
            self._message("nothing to plot yet: fix the problems of this step first")
            return
        try:
            drawn = getattr(self, "_" + self.mode)(spec, table, (table.times_s - table.times_s[0]) / 3600.0)
        except (KeyError, ValueError, IndexError, FloatingPointError):
            self.fig.clear()
            self._message("nothing to plot yet: see Needs attention on the right")
            return
        if not drawn:
            self._message({"sources": "no emissions and no species held to the data",
                           "mixing": "no dilution or deposition chosen", "conditions": "no conditions in the table"}[self.mode])
            return
        for ax in self.fig.axes:
            ax.tick_params(labelsize=7); ax.set_xlabel("hours since the first row", fontsize=8); ax.grid(alpha=0.25)
        self.canvas.draw(); self.plotted = True

    def _series(self, ax, hours, series: Dict[str, np.ndarray], ylabel: str) -> None:
        for i, (name, v) in enumerate(list(series.items())[:len(COLORS)]):
            ax.plot(hours, v, color=COLORS[i], lw=1.3, marker="o", ms=2.5, label=name)
        ax.set_ylabel(ylabel, fontsize=8)
        if len(series) >= 2:
            ax.legend(fontsize=7, frameon=False)
        if len(series) > len(COLORS):
            self.note.setText(f"{self.NOTES[self.mode]} The first {len(COLORS)} of {len(series)} are shown.")

    def _conditions(self, spec, t, hours) -> bool:
        panels = [(k, label) for k, label in (("t", "T [K]"), ("h", "boundary layer h [m]"), ("H2O", "H2O [vol-%]")) if k in t.conditions]
        extra = [("jNO2 [1/s]", t.photolysis["NO2"])] if "NO2" in t.photolysis else []
        n = len(panels) + len(extra)
        for i, (k, label) in enumerate(panels[:3]):
            self._series(self.fig.add_subplot(1, n, i + 1), hours, {k: t.conditions[k]}, label)
        for label, v in extra:
            self._series(self.fig.add_subplot(1, n, n), hours, {"jNO2": v}, label)
        return n > 0

    def _sources(self, spec, t, hours) -> bool:
        panels = []
        if t.emissions:
            panels.append((fs.emission_rates(t, spec), "emission [molec cm⁻³ s⁻¹]"))
        if t.constraints:
            panels.append((t.constraints, "held to data [ppb]"))
        for i, (series, label) in enumerate(panels):
            self._series(self.fig.add_subplot(1, len(panels), i + 1), hours, series, label)
        return bool(panels)

    def _mixing(self, spec, t, hours) -> bool:
        series: Dict[str, np.ndarray] = {}
        if spec.dilution != "none":
            series["dilution"] = fs.dilution_rate(t, spec) * 3600.0
        if spec.run_type == "constrained":
            series["species not held"] = np.full(t.n, 1.0 / spec.generic_lifetime_h)
        for sp, k in fs.deposition_rates(t, spec).items():
            series[f"deposition {sp}"] = k * 3600.0
        if series:
            self._series(self.fig.add_subplot(111), hours, series, "loss [1/h]")
        return bool(series)


class FieldReview(QWidget):
    """Always visible: summary, what needs attention (a click opens its step), what goes into the model, the table."""

    def __init__(self, win: "FieldWindow"):
        super().__init__()
        self.win = win
        lay = QVBoxLayout(self); lay.setContentsMargins(4, 4, 4, 4)
        lay.addWidget(QLabel("<b>Review</b> — follows your entries"))
        self.tabs = QTabWidget(); lay.addWidget(self.tabs, stretch=1)
        st = QWidget(); sl = QVBoxLayout(st); sl.setSpacing(4)
        self.summary = QLabel(""); self.summary.setWordWrap(True); sl.addWidget(self.summary)
        self.checks_title = QLabel("<b>Needs attention</b>"); sl.addWidget(self.checks_title)
        self.checks = QListWidget(); self.checks.setWordWrap(True)
        self.checks.setToolTip("click an entry: the step it belongs to opens")
        self.checks.itemClicked.connect(lambda it: self.win.steps.setCurrentRow(it.data(Qt.ItemDataRole.UserRole)))
        sl.addWidget(self.checks, stretch=1)
        sl.addWidget(QLabel("<b>What goes into the model</b>"))
        self.bullets = QListWidget(); self.bullets.setWordWrap(True); sl.addWidget(self.bullets, stretch=1)
        self.tabs.addTab(st, "Summary")
        self.preview = _readonly_table(["Role", "Name", "Unit", "Min", "Max"])
        self.tabs.addTab(self.preview, "Input table")

    def refresh(self, spec: fs.FieldSpec, table: Optional[fs.FieldTable], rep: fs.Report, blank: bool) -> None:
        self.checks.clear(); self.bullets.clear()
        _fill(self.preview, preview_rows(table, spec) if table is not None else [])
        if blank:
            self.summary.setText("<i>Nothing entered yet.</i><br>Give the case a name, choose the kind of run and its input table.")
            self.checks_title.setText("<b>Needs attention</b>")
            return
        self.summary.setText(f"<b>{spec.name or '(no name)'}</b> — {fs.RUN_TYPES[spec.run_type]}" + (f", {spec.site}" if spec.site else ""))
        n_err, n_warn = len(rep.errors), len(rep.warnings)
        self.checks_title.setText("<b>Needs attention</b>" + (f" — {n_err} to fix, {n_warn} to check" if n_err + n_warn else " — nothing"))
        if not n_err + n_warn:
            self.checks.addItem(QListWidgetItem("✔ No problems found"))
        for level, step, message in sorted(rep.issues, key=lambda i: (i[0] != "error", i[1])):
            it = QListWidgetItem(("✖ " if level == "error" else "⚠ ") + message); it.setData(Qt.ItemDataRole.UserRole, step)
            self.checks.addItem(it)
        for b in self._bullets(spec, table, rep):
            self.bullets.addItem(QListWidgetItem("• " + b))

    @staticmethod
    def _bullets(spec: fs.FieldSpec, table: Optional[fs.FieldTable], rep: fs.Report) -> List[str]:
        out = [f"Chemistry: {fs.chemistry_label(spec)}"]
        out += [f"Input table: {i}" for i in rep.info]
        out.append(f"Dilution: {fs.DILUTION_SOURCES.get(spec.dilution, spec.dilution)}")
        n_dep = len(spec.deposition_velocity_cm_s) + len(spec.deposition_lifetime_h)
        out.append(f"Deposition: {n_dep} species" if n_dep else "Deposition: none")
        if spec.run_type == "constrained":
            out.append(f"Species not held: lifetime {spec.generic_lifetime_h:g} h; spin-up {spec.spin_up_days:g} days")
        else:
            out.append(f"Initial values: {len(spec.initials)} species set")
        out.append(f"Output step: {spec.step_s} s")
        return out


class FieldRunPage(QWidget):
    """Last step: write the input files (field_scenario) → run scenario_runner → open the results."""

    def __init__(self, win: "FieldWindow"):
        super().__init__()
        self.win, self.json_path, self._proc, self._t0 = win, "", None, 0.0
        lay = QVBoxLayout(self)
        self.state = QLabel(""); self.state.setWordWrap(True); lay.addWidget(self.state)
        g = QGridLayout(); r = 0
        g.addWidget(QLabel("Case folder"), r, 0); self.folder = QLabel(""); self.folder.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        g.addWidget(self.folder, r, 1, 1, 3); r += 1
        self.kvpd = QCheckBox("Also write the reaction rates and production/loss budgets (large files; needed for the process "
                              "analysis and reactivity views in the viewer)")
        self.kvpd.setChecked(True); g.addWidget(self.kvpd, r, 0, 1, 4); r += 1
        g.addWidget(QLabel("Solver"), r, 0); self.solver = QComboBox()
        self._profiles = solver_profiles(self.win.toml)
        for name, p in self._profiles.items():
            self.solver.addItem(p["label"], name)
        dflt = default_profile(self.win.toml)
        if dflt in self._profiles:
            self.solver.setCurrentIndex(self.solver.findData(dflt))
        self.solver.setToolTip("\n\n".join(f"{p['label']}:\n{p['summary']}" for p in self._profiles.values()))
        self.solver.currentIndexChanged.connect(self._solver_changed)
        g.addWidget(self.solver, r, 1, 1, 3); r += 1
        lay.addLayout(g)
        lay.addWidget(_info("What the field adds (emissions, dilution towards the background, species held to the data, deposition) "
                            "is written as ordinary reactions and inputs, so the same files run in FACSIMILE. "
                            "Every input table row is interpolated to the output step; emissions are held from one row to the next."))
        hb = QHBoxLayout()
        self.prep_btn = QPushButton("Prepare input files"); self.prep_btn.clicked.connect(self.prepare); hb.addWidget(self.prep_btn)
        self.run_btn = QPushButton("Run"); self.run_btn.setObjectName("primary"); self.run_btn.setEnabled(False); self.run_btn.clicked.connect(self.run)
        hb.addWidget(self.run_btn)
        self.stop_btn = QPushButton("Stop"); self.stop_btn.setEnabled(False); self.stop_btn.clicked.connect(self.stop); hb.addWidget(self.stop_btn)
        self.view_btn = QPushButton("Open results"); self.view_btn.setEnabled(False); self.view_btn.clicked.connect(self.open_results); hb.addWidget(self.view_btn)
        hb.addStretch(); lay.addLayout(hb)
        self.bar = QProgressBar(); self.bar.setRange(0, 100); self.bar.setValue(0); self.bar.setFormat(""); lay.addWidget(self.bar)
        self.files = QLabel(""); self.files.setWordWrap(True); self.files.setStyleSheet("color: #555555;"); lay.addWidget(self.files)
        self.log = QTextEdit(); self.log.setReadOnly(True); self.log.setMinimumHeight(260)
        self.log.setStyleSheet("font-family: Consolas, monospace; font-size: 10px;"); lay.addWidget(self.log, stretch=1)

    def cs_path(self) -> str:
        return os.path.join(self.win.out_dir(), self.win.name.text().strip() + ".cs")

    def refresh(self):
        self.folder.setText(self.win.out_dir())
        rep = self.win.case[2]
        self.state.setText("The case is complete." if rep.ok else
                           f"{len(rep.errors)} problem(s) to fix first; see Needs attention on the right.")
        self.prep_btn.setEnabled(rep.ok and self._proc is None)
        self.view_btn.setEnabled(os.path.isfile(self.cs_path()))

    def _solver_changed(self, *_a):
        self.json_path = ""; self.run_btn.setEnabled(False)       # the prepared JSON carries the solver: prepare again

    def prepare(self) -> str:
        """Write ENZ, EASY, chemistry and the scenario JSON into the case folder.  Returns the JSON path ('' on error)."""
        self.win.refresh_all()
        spec, table, rep = self.win.case
        if not rep.ok or table is None:
            self.log.append("Cannot prepare:\n  " + "\n  ".join(rep.errors)); self.bar.setFormat("fix the problems first"); return ""
        self.bar.setValue(0); self.bar.setFormat("preparing …"); QApplication.processEvents()
        buf = io.StringIO()
        try:
            with contextlib.redirect_stdout(buf):
                self.json_path, res = fsc.prepare(spec, table, self.win.out_dir(), self.win.include_dir, self.win.toml,
                                                  solver_from_toml(self.win.toml, self.solver.currentData()), self.kvpd.isChecked())
        except Exception as e:
            self.log.append(f"Prepare failed: {e}"); self.bar.setFormat("prepare failed"); return ""
        for line in buf.getvalue().splitlines():
            if "[scenario_prep]" in line or "WARNING" in line:
                self.log.append(line)
        for n in res.get("notes", []):
            self.log.append("note: " + n)
        out = os.path.dirname(self.json_path)
        missing = []
        try:
            with open(os.path.join(out, spec.name + "_juice.elog"), encoding="utf-8", errors="replace") as fh:
                missing = sorted({m.group(1).strip() for m in re.finditer(r"^;#include (.+?) ! NOT FOUND", fh.read(), re.M)})
        except OSError:
            pass
        if missing:
            self.log.append(f"✖ {len(missing)} include file(s) not found: " + ", ".join(missing)
                            + "\n  Add the folder that holds them under 'Include folders' on the Chemistry step.")
        self.files.setText("Written: " + ", ".join(os.path.basename(p) for p in (res["enz"], res["easy"], res["spec_path"], self.json_path)))
        self.log.append(f"Inputs ready in {out}")
        self.bar.setValue(100); self.bar.setFormat("inputs ready" if not missing else f"inputs ready — {len(missing)} include(s) missing")
        self.run_btn.setEnabled(True); self.view_btn.setEnabled(os.path.isfile(self.cs_path()))
        return self.json_path

    def run(self):
        if not self.json_path:
            return
        args = ["-u", os.path.join(fsc.HERE, "scenario_runner.py"), self.json_path] + (["--toml", self.win.toml] if self.win.toml else [])
        self._proc = QProcess(self)
        env = QProcessEnvironment.systemEnvironment(); env.insert("PYTHONIOENCODING", "utf-8"); env.insert("PYTHONUTF8", "1")
        self._proc.setProcessEnvironment(env); self._proc.setProcessChannelMode(QProcess.ProcessChannelMode.MergedChannels)
        self._proc.readyReadStandardOutput.connect(self._on_output); self._proc.finished.connect(self._on_finished)
        self._t0 = time.monotonic()
        self.bar.setValue(0); self.bar.setFormat("running   0%")
        for b, on in ((self.run_btn, False), (self.prep_btn, False), (self.stop_btn, True), (self.view_btn, False)):
            b.setEnabled(on)
        self.log.append("$ python scenario_runner.py " + os.path.basename(self.json_path))
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
        try:
            info = {"exit": code, "elapsed": time.monotonic() - self._t0, "finished": time.strftime("%Y-%m-%d %H:%M:%S"),
                    "cs_bytes": os.path.getsize(self.cs_path()) if ok else 0}
            with open(self.json_path[:-5] + ".run.json", "w", encoding="utf-8") as fh:
                json.dump(info, fh, indent=1)
        except Exception:
            pass
        self._proc = None
        self.bar.setValue(100 if ok else 0)
        self.bar.setFormat(f"done   {mm}m {ss:02d}s" if ok else f"failed (exit {code})   {mm}m {ss:02d}s")
        for b, on in ((self.run_btn, True), (self.prep_btn, True), (self.stop_btn, False), (self.view_btn, ok)):
            b.setEnabled(on)
        self.win.set_run_state("ok" if ok else "error")

    def open_results(self):
        cs = self.cs_path()
        if os.path.isfile(cs):
            p = QProcess(self); p.start(sys.executable, [os.path.join(fsc.HERE, "results_viewer.py"), cs]); self._viewer = p


class FieldWindow(QMainWindow):
    def __init__(self, path: str = "", toml: str = ""):
        super().__init__()
        self.toml, self.path, self.dirty, self._loading = toml, "", False, False
        paths = load_paths(toml)
        self.include_dir, self.output_root = paths["include_dir"], paths["output_root"]
        self.easy_root = self.output_root if os.path.isdir(self.output_root) else paths["easy_root"]
        self.setWindowIcon(theme.app_icon()); self.resize(1550, 920)
        self._timer = QTimer(self); self._timer.setSingleShot(True); self._timer.setInterval(400)
        self._timer.timeout.connect(self.refresh_all)

        left = QWidget(); ll = QVBoxLayout(left); ll.setContentsMargins(6, 6, 6, 6); left.setFixedWidth(270)
        if os.path.isfile(theme.LOGO):
            lg = QLabel(); lg.setPixmap(QPixmap(theme.LOGO).scaledToWidth(150, Qt.TransformationMode.SmoothTransformation))
            lg.setAlignment(Qt.AlignmentFlag.AlignCenter); ll.addWidget(lg)
        t = QLabel("JUICE"); t.setObjectName("title"); t.setAlignment(Qt.AlignmentFlag.AlignCenter); ll.addWidget(t)
        st = QLabel("field observations → model run"); st.setObjectName("subtitle"); st.setAlignment(Qt.AlignmentFlag.AlignCenter); ll.addWidget(st)
        self.steps = QListWidget(); self.steps.setObjectName("steps")
        for i, s in enumerate(STEPS):
            self.steps.addItem(QListWidgetItem(f"{i + 1}   {s}"))
        self._state = ["todo"] * len(STEPS)
        ll.addWidget(self.steps, stretch=1)

        self.plot_conditions = FieldPlot(self, "conditions"); self.plot_sources = FieldPlot(self, "sources", 3.2)
        self.plot_mixing = FieldPlot(self, "mixing")
        self.stack = QStackedWidget()
        self.stack.addWidget(_scroll(*self._case_page(), _boxed(self.plot_conditions, "Conditions the model will see")))
        self.page_chem = self._chemistry_page()
        self.stack.addWidget(_scroll(_boxed(self.page_chem, "Chemistry")))
        self.stack.addWidget(_scroll(*self._sources_page(), _boxed(self.plot_sources, "Emissions and constraints against time")))
        self.stack.addWidget(_scroll(*self._mixing_page(), _boxed(self.plot_mixing, "Loss per hour by process")))
        self.page_run = FieldRunPage(self)
        self.stack.addWidget(_scroll(_boxed(self.page_run, "Run the model")))
        nav = QHBoxLayout()
        self.back_btn = QPushButton("◀ Back"); self.back_btn.clicked.connect(lambda: self._step(-1)); nav.addWidget(self.back_btn)
        nav.addStretch()
        self.next_btn = QPushButton("Next ▶"); self.next_btn.setObjectName("primary"); self.next_btn.clicked.connect(lambda: self._step(+1))
        nav.addWidget(self.next_btn)
        mid = QWidget(); ml = QVBoxLayout(mid); ml.setContentsMargins(0, 0, 0, 0); ml.addWidget(self.stack, stretch=1); ml.addLayout(nav)

        self.review = FieldReview(self); self.review.setMinimumWidth(360); self.review.setMaximumWidth(560)
        split = QSplitter(); split.addWidget(left); split.addWidget(mid); split.addWidget(self.review)
        split.setStretchFactor(0, 0); split.setStretchFactor(1, 5); split.setStretchFactor(2, 2)
        split.setSizes([270, 900, 420]); split.setChildrenCollapsible(False)
        self.setCentralWidget(split)

        m = self.menuBar().addMenu("File")
        m.addAction("New field case", QKeySequence("Ctrl+N"), self.new_case)
        m.addAction("Open field case (.field.json) …", QKeySequence("Ctrl+O"), self._open_dialog)
        m.addAction("Save field case", QKeySequence("Ctrl+S"), self.save)
        m.addAction("Save field case as …", self.save_as)
        m.addAction("Open case folder", self._open_folder)
        m.addSeparator()
        m.addAction("Quit", self.close)

        self._wire()
        self.case: Tuple[fs.FieldSpec, Optional[fs.FieldTable], fs.Report] = (fs.FieldSpec(), None, fs.Report())
        self.steps.currentRowChanged.connect(self._go)
        self._apply(self._blank())
        self.steps.setCurrentRow(0)
        if path:
            self.open_path(path)
        self._update_title()

    @staticmethod
    def _blank() -> fs.FieldSpec:
        return fs.FieldSpec(name="", dilution=fs.DEFAULT_DILUTION["constrained"])

    # ── pages ────────────────────────────────────────────────────────────────────────────────────
    def _case_page(self) -> List[QWidget]:
        g1 = QGroupBox("What kind of field run?"); gv = QVBoxLayout(g1); gv.setSpacing(2)
        self.run_group = QButtonGroup(self); self.run_buttons = {}
        for key, title in fs.RUN_TYPES.items():
            b = QRadioButton(title); d = _info(RUN_TEXT[key]); d.setContentsMargins(22, 0, 0, 6)
            self.run_group.addButton(b); self.run_buttons[key] = b
            gv.addWidget(b); gv.addWidget(d)
        g2 = QGroupBox("The case"); f = QFormLayout(g2)
        self.name = self._line("e.g. pasadena_2021-08-15"); self.site = self._line("e.g. Pasadena, CalTech")
        self.folder = self._line(""); self.folder.setText(self.easy_root)
        self.date = self._line("YYYY-MM-DD, only for 'HH:MM' or hours-of-day times")
        self.lat = self._line("degrees north, when the table has no lat column")
        self.lon = self._line("degrees east, when the table has no lon column")
        self.description = self._line("")
        self.step_s = QSpinBox(); self.step_s.setRange(1, 3600); self.step_s.setSuffix(" s")
        fh = QHBoxLayout(); fh.addWidget(self.folder, 1); b = QPushButton("Choose…"); b.clicked.connect(self._choose_folder); fh.addWidget(b)
        f.addRow("Name", self.name); f.addRow("Project folder", fh)
        for label, w in (("Site", self.site), ("Date", self.date), ("Latitude", self.lat), ("Longitude", self.lon),
                         ("Notes", self.description), ("Output step", self.step_s)):
            f.addRow(label, w)
        g3 = QGroupBox("Input table"); gv = QVBoxLayout(g3)
        h = QHBoxLayout(); self.table = self._line("CSV or Excel, one row per time"); h.addWidget(self.table, 1)
        b = QPushButton("Choose…"); b.clicked.connect(self._choose_table); h.addWidget(b); gv.addLayout(h)
        h = QHBoxLayout(); b = QPushButton("Write a template for this kind of run…"); b.clicked.connect(self._write_template); h.addWidget(b); h.addStretch()
        gv.addLayout(h)
        self.kind_box = QGroupBox(""); self.kind_form = QFormLayout(self.kind_box)
        self.kind_needs = QLabel(""); self.kind_needs.setWordWrap(True); self.kind_form.addRow(self.kind_needs)
        self.spin_up = self._spin(0.0, 60.0, 1, " days")
        self.spin_up.setToolTip("Repeat the first day before the output, so that unmeasured products reach a steady cycle.")
        self.kind_form.addRow("Spin-up", self.spin_up)
        return [g1, g2, g3, self.kind_box]

    def _chemistry_page(self) -> ChemistryPage:
        page = ChemistryPage(self)
        page.setSubTitle("Which chemistry does this case use? The matching part of the MCM is extracted in the background, "
                         "or pick an existing mechanism include set, which may be any mechanism written as include files.")
        page.mode_inorganic.setToolTip("NOx, O3, HOx, CO, SO2 chemistry and photolysis only.")
        return page

    def _sources_page(self) -> List[QWidget]:
        self.emit_box = QGroupBox("Emissions"); v = QVBoxLayout(self.emit_box)
        f = QFormLayout(); self.unit = QComboBox(); self.unit.addItems(list(fs.EMISSION_UNITS))
        self.unit.setToolTip("Used for emission columns that have no unit in brackets.")
        self.area = self._spin(0.0, 1e6, 2, " km²"); self.area.setToolTip("Footprint A in E/(h·A), needed for the mol/s and molec/s units.")
        f.addRow("Default emission unit", self.unit); f.addRow("Footprint area", self.area); v.addLayout(f)
        self.emit_table = _readonly_table(["Species", "Unit", "Min", "Max"]); v.addWidget(self.emit_table)
        v.addWidget(_info("Emissions come from the E_<species> columns of the input table."))
        self.held_box = QGroupBox("Held to measurements"); v = QVBoxLayout(self.held_box)
        self.held_table = _readonly_table(["Species", "Min [ppb]", "Max [ppb]", "Empty rows"]); v.addWidget(self.held_table)
        v.addWidget(_info("Species held come from the C_<species> columns of the input table. "
                          "Leave out a measured product to compare the model with it afterwards."))
        return [self.emit_box, self.held_box]

    def _mixing_page(self) -> List[QWidget]:
        g = QGroupBox("Dilution and background"); v = QVBoxLayout(g); self.mix_form = f = QFormLayout()
        self.dilution = QComboBox()
        for k, text in fs.DILUTION_SOURCES.items():
            self.dilution.addItem(text, k)
        self.k_dil = self._spin(0.0, 100.0, 3, " 1/h"); self.k_collapse = self._spin(0.0, 100.0, 3, " 1/h")
        self.k_collapse.setToolTip("Dilution while the boundary layer does not grow. Stockwell et al. (2024) use 0.3 per hour.")
        self.tracer = self._line("a species with a C_ column, e.g. CO")
        for label, w in (("Dilution", self.dilution), ("Dilution rate", self.k_dil), ("Rate when layer stops growing", self.k_collapse), ("Tracer", self.tracer)):
            f.addRow(label, w)
        v.addLayout(f)
        self.back_table = _readonly_table(["Species", "Min [ppb]", "Max [ppb]"]); v.addWidget(self.back_table)
        v.addWidget(_info("Background comes from the B_<species> columns of the input table. Species without one are diluted with clean air."))
        d = QGroupBox("Deposition"); v = QVBoxLayout(d)
        self.deposition = _table(["Species", "Velocity [cm/s]", "or lifetime [h]"]); self.deposition.setMinimumHeight(110); v.addWidget(self.deposition)
        v.addLayout(self._row_buttons(self.deposition))
        v.addWidget(_info("A velocity is divided by the boundary-layer height, so it needs the column h. A lifetime is a constant first-order loss."))
        self.life_box = QGroupBox("Species not held to the data"); f2 = QFormLayout(self.life_box)
        self.lifetime = self._spin(0.1, 1e4, 1, " h")
        self.lifetime.setToolTip("First-order loss of every species that is not held to the data, standing in for mixing and deposition.")
        f2.addRow("Lifetime", self.lifetime)
        self.init_box = QGroupBox("Initial values"); v = QVBoxLayout(self.init_box)
        self.initials = _table(["Species", "Initial [ppb]"]); self.initials.setMinimumHeight(110); v.addWidget(self.initials)
        v.addLayout(self._row_buttons(self.initials))
        return [g, d, self.life_box, self.init_box]

    def _row_buttons(self, table: QTableWidget) -> QHBoxLayout:
        h = QHBoxLayout()
        a = QPushButton("+ row"); a.clicked.connect(lambda: _add_row(table, [""] * table.columnCount())); h.addWidget(a)
        b = QPushButton("− row"); b.clicked.connect(lambda: _del_row(table)); h.addWidget(b); h.addStretch()
        return h

    def _line(self, placeholder: str) -> QLineEdit:
        w = QLineEdit(); w.setPlaceholderText(placeholder)
        return w

    def _spin(self, lo: float, hi: float, decimals: int, suffix: str) -> QDoubleSpinBox:
        w = QDoubleSpinBox(); w.setRange(lo, hi); w.setDecimals(decimals); w.setSuffix(suffix)
        return w

    def _wire(self):
        for b in self.run_buttons.values():
            b.toggled.connect(self._run_type_toggled)
        for w in (self.name, self.site, self.folder, self.date, self.lat, self.lon, self.description, self.table, self.tracer):
            w.textChanged.connect(self._changed)
        for w in (self.step_s, self.spin_up, self.area, self.k_dil, self.k_collapse, self.lifetime):
            w.valueChanged.connect(self._changed)
        for w in (self.unit, self.dilution):
            w.currentIndexChanged.connect(self._changed)
        for t in (self.deposition, self.initials):
            t.itemChanged.connect(self._changed); t.model().rowsRemoved.connect(self._changed)
        c = self.page_chem
        for w in (c.base_combo, c.mech_combo):
            w.currentIndexChanged.connect(self._changed)
        c.mech_combo.currentTextChanged.connect(self._changed)
        c.addon_combo.changed.connect(self._changed)
        for w in (c.mode_compounds, c.mode_existing, c.mode_inorganic):
            w.toggled.connect(self._changed)
        for lst in (c.chosen, c.include_dirs):
            lst.model().rowsInserted.connect(self._changed); lst.model().rowsRemoved.connect(self._changed)
        c.adjust.itemChanged.connect(self._changed); c.adjust.model().rowsRemoved.connect(self._changed)

    # ── form ⇄ spec ──────────────────────────────────────────────────────────────────────────────
    def run_type(self) -> str:
        return next(k for k, b in self.run_buttons.items() if b.isChecked())

    def _species_values(self, rows: List[List[str]], label: str, problems: List[Tuple[int, str]]) -> Dict[str, float]:
        out: Dict[str, float] = {}
        for name, text in rows:
            species = name.strip().upper()
            try:
                value = float(text.replace(",", "."))
            except ValueError:
                problems.append((fs.STEP_MIXING, f"{label}: '{text}' is not a number (for {name or 'a row without a species'})"))
                continue
            if not species:
                problems.append((fs.STEP_MIXING, f"{label}: a row has a value but no species"))
            elif species in out:
                problems.append((fs.STEP_MIXING, f"{label}: {species} is listed twice"))
            else:
                out[species] = value
        return out

    def spec_from_form(self) -> Tuple[fs.FieldSpec, List[Tuple[int, str]]]:
        """The case as entered, plus (step, what could not be understood in the form)."""
        problems: List[Tuple[int, str]] = []

        def number(w: QLineEdit, label: str) -> Optional[float]:
            s = w.text().strip().replace(",", ".")
            if not s:
                return None
            try:
                return float(s)
            except ValueError:
                problems.append((fs.STEP_CASE, f"{label} '{w.text()}' is not a number"))
                return None
        initials = self._species_values(_rows(self.initials), "Initial values", problems)
        for sp, v in initials.items():
            if v < 0:
                problems.append((fs.STEP_MIXING, f"Initial values: {sp} cannot be negative"))
        dep = _rows(self.deposition)
        for name, vel, life in dep:
            if bool(vel) == bool(life):
                problems.append((fs.STEP_MIXING, f"Deposition of {name or 'a row without a species'}: give a velocity or a lifetime, one of the two"))
        velocity = self._species_values([[n, v] for n, v, _l in dep if v and not _l], "Deposition", problems)
        lifetime = self._species_values([[n, l] for n, _v, l in dep if l and not _v], "Deposition", problems)
        tmp = ExperimentSpec(name="x"); self.page_chem.apply(tmp)
        chemistry = {k: getattr(tmp, k) for k in CHEM_KEYS}
        chemistry["reaction_adjustments"] = [list(a) for a in chemistry["reaction_adjustments"]]
        spec = fs.FieldSpec(name=self.name.text().strip(), description=self.description.text().strip(), run_type=self.run_type(),
                            table=self.table.text().strip(), date=self.date.text().strip(), site=self.site.text().strip(),
                            lat=number(self.lat, "Latitude"), lon=number(self.lon, "Longitude"), chemistry=chemistry,
                            step_s=self.step_s.value(), default_emission_unit=self.unit.currentText(), area_km2=self.area.value(),
                            dilution=self.dilution.currentData(), k_dil_per_h=self.k_dil.value(), k_collapse_per_h=self.k_collapse.value(),
                            tracer=self.tracer.text().strip(), spin_up_days=self.spin_up.value(), generic_lifetime_h=self.lifetime.value(),
                            initials=initials, deposition_velocity_cm_s=velocity, deposition_lifetime_h=lifetime)
        return spec, problems

    def _apply(self, spec: fs.FieldSpec):
        self._loading = True
        try:
            self.run_buttons.get(spec.run_type, self.run_buttons["constrained"]).setChecked(True)
            for w, v in ((self.name, spec.name), (self.site, spec.site), (self.date, spec.date), (self.description, spec.description),
                         (self.table, spec.table), (self.tracer, spec.tracer), (self.lat, "" if spec.lat is None else f"{spec.lat:g}"),
                         (self.lon, "" if spec.lon is None else f"{spec.lon:g}")):
                w.setText(v)
            self.unit.setCurrentText(spec.default_emission_unit)
            self.dilution.setCurrentIndex(max(0, self.dilution.findData(spec.dilution)))
            for w, v in ((self.area, spec.area_km2), (self.k_dil, spec.k_dil_per_h), (self.k_collapse, spec.k_collapse_per_h),
                         (self.spin_up, spec.spin_up_days), (self.lifetime, spec.generic_lifetime_h)):
                w.setValue(v)
            self.step_s.setValue(spec.step_s)
            _set_rows(self.initials, [[k, f"{v:g}"] for k, v in spec.initials.items()])
            _set_rows(self.deposition, [[k, f"{v:g}", ""] for k, v in spec.deposition_velocity_cm_s.items()]
                      + [[k, "", f"{v:g}"] for k, v in spec.deposition_lifetime_h.items()])
            tmp = ExperimentSpec(name="x")
            for k in CHEM_KEYS:
                if k in (spec.chemistry or {}):
                    setattr(tmp, k, spec.chemistry[k])
            self.page_chem.set_from(tmp)
        finally:
            self._loading = False
        self._timer.stop()
        self._sync(); self.refresh_all()

    def _run_type_toggled(self, checked: bool):
        if not checked:
            return
        if not self._loading:
            self.dilution.setCurrentIndex(self.dilution.findData(fs.DEFAULT_DILUTION[self.run_type()]))
        self._changed()

    def _sync(self):
        """Show what belongs to the kind of run and the dilution source."""
        rt, d = self.run_type(), self.dilution.currentData()
        self.kind_box.setTitle(f"What this run needs: {fs.RUN_TYPES[rt]}")
        self.kind_needs.setText("In the input table:<br>" + "<br>".join("• " + n for n in NEEDS[rt]))
        self.kind_form.setRowVisible(self.spin_up, rt == "constrained")
        self.emit_box.setVisible(rt != "constrained")
        self.life_box.setVisible(rt == "constrained")
        self.init_box.setVisible(rt != "constrained")
        for w, on in ((self.k_dil, d == "constant"), (self.k_collapse, d == "pbl"), (self.tracer, d == "tracer")):
            self.mix_form.setRowVisible(w, on)

    def _changed(self, *_):
        if self._loading:
            return
        self.dirty = True
        self._sync(); self._update_title(); self._timer.start()

    # ── review ───────────────────────────────────────────────────────────────────────────────────
    def is_blank(self, spec: fs.FieldSpec) -> bool:
        return not (spec.name or spec.table or fs.chemistry_chosen(spec))

    def refresh_all(self) -> fs.Report:
        spec, problems = self.spec_from_form()
        table, rep = fs.check(spec)
        for step, message in problems:
            rep.err(message, step)
        self.case = (spec, table, rep)
        blank = self.is_blank(spec)
        self.review.refresh(spec, table, rep, blank)
        self._fill_step_tables(spec, table)
        for i in range(RUN_STEP):
            level = "todo" if blank else (rep.worst(i) or "ok")
            self._state[i] = level
        self._paint_states()
        theme.mark_missing(self.name, not spec.name)
        theme.mark_missing(self.table, not spec.table)
        if hasattr(self, "page_run"):
            self.page_run.refresh()
        plot = {0: self.plot_conditions, 2: self.plot_sources, 3: self.plot_mixing}.get(self.steps.currentRow())
        if plot is not None:
            plot.plot()
        return rep

    def _fill_step_tables(self, spec: fs.FieldSpec, t: Optional[fs.FieldTable]):
        if t is None:
            for w in (self.emit_table, self.held_table, self.back_table):
                w.setRowCount(0)
            self.held_box.setVisible(spec.run_type == "constrained")
            return
        _fill(self.emit_table, [(k, fs.resolve_emission_unit(t, spec, k), *_span(v)) for k, v in t.emissions.items()])
        _fill(self.held_table, [(k, *_span(v), str(int((~np.isfinite(v)).sum()))) for k, v in t.constraints.items()])
        _fill(self.back_table, [(k, *_span(v)) for k, v in t.background.items()])
        self.held_box.setVisible(spec.run_type == "constrained" or bool(t.constraints))

    def set_run_state(self, state: str):
        self._state[RUN_STEP] = state
        self._paint_states()

    def _paint_states(self):
        marks = {"todo": "", "ok": "  ✔", "warning": "  ⚠", "error": "  ✖"}
        for i, s in enumerate(STEPS):
            self.steps.item(i).setText(f"{i + 1}   {s}" + marks.get(self._state[i], ""))

    def _go(self, row: int):
        if row < 0:
            return
        self.stack.setCurrentIndex(row)
        self.back_btn.setEnabled(row > 0); self.next_btn.setEnabled(row < RUN_STEP)
        self.refresh_all()

    def _step(self, d: int):
        self.steps.setCurrentRow(max(0, min(RUN_STEP, self.steps.currentRow() + d)))

    # ── files ────────────────────────────────────────────────────────────────────────────────────
    def out_dir(self) -> str:
        return os.path.join(self.folder.text().strip() or self.easy_root, self.name.text().strip() or "field_case")

    def _choose_folder(self):
        p = QFileDialog.getExistingDirectory(self, "Project folder", self.folder.text().strip() or self.easy_root)
        if p:
            self.folder.setText(p)

    def _dialog_dir(self) -> str:
        return os.path.dirname(self.path) if self.path else (self.folder.text().strip() or self.easy_root)

    def _choose_table(self):
        p, _ = QFileDialog.getOpenFileName(self, "Input table", self._dialog_dir(), TABLE_FILTER)
        if p:
            self.table.setText(p)

    def _write_template(self):
        rt = self.run_type()
        p, _ = QFileDialog.getSaveFileName(self, "Write a template", os.path.join(self._dialog_dir(), f"{rt}_template.csv"), "CSV (*.csv)")
        if not p:
            return
        try:
            fs.write_template(p, rt)
        except OSError as e:
            QMessageBox.warning(self, "JUICE", f"The template could not be written:\n{e}")
            return
        self.table.setText(p)

    def _open_folder(self):
        d = os.path.dirname(self.path) if self.path else self.out_dir()
        if os.path.isdir(d):
            _open_path(d)
        else:
            QMessageBox.information(self, "JUICE", "The case has no folder yet. Save it first.")

    def new_case(self):
        if self._confirm_discard():
            self.path = ""
            self._apply(self._blank())
            self.dirty = False; self.steps.setCurrentRow(0); self._update_title()

    def _open_dialog(self):
        if not self._confirm_discard():
            return
        p, _ = QFileDialog.getOpenFileName(self, "Open a field case", self._dialog_dir(), CASE_FILTER)
        if p:
            self.open_path(p)

    def open_path(self, path: str) -> bool:
        try:
            spec = fs.load_spec(path)
        except (OSError, ValueError) as e:
            QMessageBox.warning(self, "JUICE", f"This field case could not be opened:\n{e}")
            return False
        self.path = os.path.abspath(path)
        self._loading = True
        self.folder.setText(os.path.dirname(os.path.dirname(self.path)))
        self._loading = False
        self._apply(spec)
        self.dirty = False; self.steps.setCurrentRow(0); self._update_title()
        return True

    def write_to(self, path: str) -> bool:
        spec, problems = self.spec_from_form()
        if problems:
            QMessageBox.warning(self, "JUICE", "Please correct this first:\n" + "\n".join(m for _s, m in problems))
            return False
        try:
            os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
            self.path = os.path.abspath(fs.save_spec(spec, path))
        except OSError as e:
            QMessageBox.warning(self, "JUICE", f"The case could not be saved:\n{e}")
            return False
        self.dirty = False; self._update_title()
        self.statusBar().showMessage(f"saved {self.path}", 8000)
        return True

    def save(self) -> bool:
        if self.path:
            return self.write_to(self.path)
        name = self.name.text().strip()
        if not name:
            QMessageBox.warning(self, "JUICE", "Give the case a name first.")
            return False
        return self.write_to(os.path.join(self.out_dir(), name + fs.SPEC_SUFFIX))

    def save_as(self) -> bool:
        name = (self.name.text().strip() or "field_case") + fs.SPEC_SUFFIX
        p, _ = QFileDialog.getSaveFileName(self, "Save the field case", os.path.join(self.out_dir(), name), CASE_FILTER)
        return self.write_to(p) if p else False

    def _update_title(self):
        name = os.path.basename(self.path) if self.path else "new field case"
        self.setWindowTitle(f"JUICE — {name}{' *' if self.dirty else ''}")

    def _confirm_discard(self) -> bool:
        if not self.dirty:
            return True
        ans = QMessageBox.question(self, "JUICE", "Save the changes to this field case?",
                                   QMessageBox.StandardButton.Save | QMessageBox.StandardButton.Discard | QMessageBox.StandardButton.Cancel)
        if ans == QMessageBox.StandardButton.Save:
            return self.save()
        return ans == QMessageBox.StandardButton.Discard

    def closeEvent(self, ev):
        if self._confirm_discard():
            ev.accept()
        else:
            ev.ignore()


def main(argv=None) -> int:
    argv = list(sys.argv if argv is None else argv)
    toml = pick_toml(argv)
    skip = {i + 1 for i, a in enumerate(argv) if a in ("-t", "--toml")}
    paths = [a for i, a in enumerate(argv[1:], 1) if i not in skip and not a.startswith("-") and os.path.isfile(a)]
    theme.taskbar_identity()
    app = QApplication(argv); theme.apply_theme(app)
    w = FieldWindow(paths[0] if paths else "", toml); w.show(); theme.force_taskbar_icon(w)
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
