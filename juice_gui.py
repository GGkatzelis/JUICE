"""
juice_gui.py — JUICE, step by step.

One window, five steps down the left, the active step in the middle, a review panel on
the right that always shows what the experiment contains, and the solver run as the
last step.  No pop-ups on the main path.

    1 Experiment   name, folder, load (saved / wiki / pasted notes), template, date,
                   conditions, roof  — with a plot of the conditions the model will see
    2 Chemistry    precursors (MCM subset) or an existing mechanism, first-generation
                   chemistry tree, reaction adjustments (EASY k[...] overrides)
    3 Injections   actions, amounts, mixtures, delayed release — with the injection timeline
    4 Chamber      wall deposition, sources, background reactivity, dilution, initial
                   values — with a plot of what the chamber adds per hour
    5 Run          prepare the input files, run the solver, open the results

The only thing a user edits is the experiment description (<name>.experiment.json).
The input files (ENZ, EASY, scenario JSON) are generated from it in the project
folder with the IDL/FACSIMILE naming (QNO, QO3, QH2O, Q_<name>; q<NAME> switches)
and can be opened from the review panel, but never need to be.  (The former expert
tool for file-level editing, FR_GUI.py, is not part of the JUICE repository.)

    python juice_gui.py                 # paths from GG_GUI.toml if present, else FR_GUI.toml
    python juice_gui.py -t FR_GUI.toml  # another TOML

The step pages are the ones of experiment_wizard.py (shared through ExperimentHost).
Solver numerics come from the TOML (solver_from_toml), never from this window.
"""
from __future__ import annotations

import contextlib
import io
import json
import os
import re
import subprocess
import sys
import time
import tomllib
from typing import Dict, List, Optional, Tuple

import numpy as np
from PyQt6.QtCore import Qt, QProcess, QProcessEnvironment, QTimer
from PyQt6.QtGui import QPixmap, QIcon
from PyQt6.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QGridLayout, QLabel,
                             QListWidget, QListWidgetItem, QStackedWidget, QPushButton, QGroupBox, QScrollArea,
                             QCheckBox, QProgressBar, QTextEdit, QMessageBox, QFileDialog, QSplitter, QTabWidget, QFrame,
                             QDialog, QFormLayout, QDialogButtonBox, QComboBox)
from PyQt6.QtGui import QKeySequence, QColor
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

from experiment_spec import load_spec, save_spec, validate, mechanism_species, enz_columns, ExperimentSpec, ACTIONS, known_species_for, input_sources, level_icon, injection_label
from experiment_wizard import ExperimentHost, ReviewPage, CHAMBER_ITEMS, _info
from elog_generator import generate_both
from idl_parser import parse_idl_easy
import injection_calc as ic
import juice_theme as theme
import mechanism_registry as mreg

COMPARE_FLOOR = 2.5e4      # molec/cm³ ≈ 1e-3 ppt, as tests/test_solver.py: below this CPU rounding dominates
COMPARE_TOL = 0.02         # the benchmark tolerance

HERE = os.path.dirname(os.path.abspath(__file__))
LOGO = theme.LOGO
app_icon = theme.app_icon
STEPS = ["Experiment", "Chemistry", "Injections", "Chamber", "Run"]
RUN_STEP = len(STEPS) - 1
ACCENT, ACCENT_DARK, ACCENT_SOFT = theme.ACCENT, theme.ACCENT_DARK, theme.ACCENT_SOFT


# ─────────────────────────────────────────────────────────────────────────────
from scenario_prep import (load_paths, pick_toml, solver_from_toml, prepare_scenario,   # no-Qt helpers shared with study_runner
                           solver_profiles, default_profile)


def blank_spec() -> ExperimentSpec:
    """What a fresh window shows: nothing chosen, no defaults.  Defaults arrive with a
    template or an import."""
    spec = ExperimentSpec(name="")
    c = spec.conditions
    c.date = c.start = c.end = ""; c.source = ""; c.roof_open = []
    spec.injections = []
    return spec


def is_blank(spec: ExperimentSpec) -> bool:
    c = spec.conditions
    # the name alone does not count: collect_spec() substitutes "experiment" when the field is empty
    return not (c.date or c.start or c.source or spec.injections or spec.precursors or spec.mechanism_stem)


def spec_bullets(spec: ExperimentSpec) -> List[str]:
    """Plain-language list of what goes into the input files (no file syntax)."""
    c = spec.conditions
    ch = ic.chamber(c.chamber_name)
    out = []
    if spec.precursors:
        out.append(f"Chemistry: {spec.base_mechanism} subset for {', '.join(spec.precursors)} + the complete inorganic core and photolysis"
                   + (f"; add-ons: {', '.join(spec.addons)}" if spec.addons else ""))
    elif spec.mechanism_stem:
        out.append(f"Chemistry: mechanism include set '{spec.mechanism_stem}' + inorganic core and photolysis")
    elif spec.inorganic_only:
        out.append(f"Chemistry: inorganic only ({spec.base_mechanism} NOx/O3/HOx/CO/SO2 core and photolysis, chamber HONO/HCHO sources) — no VOC")
    else:
        out.append("Chemistry: not chosen yet")
    for rxn, rate, _cmt in spec.reaction_adjustments:
        out.append(f"Reaction adjusted: {rxn}  →  " + ("switched off" if (rate or "0").strip() in ("0", "0.", "0.0") else f"k = {rate}"))
    try:                                   # readiness: where every physical input comes from
        src = input_sources(spec)
        head = {"archive": f"Conditions from the SAPHIR archive for {c.date}", "file": f"Conditions from {os.path.basename(c.input_file) or 'a file'}",
                "manual": "Conditions typed as constants"}.get(c.source, "Conditions: source not chosen")
        out.append(head + " — " + "; ".join(f"{k}: {v}" for k, v in src.items()))
    except Exception:
        if c.source == "archive":
            out.append(f"Conditions: measured T, P, humidity, flow and J values from the archive for {c.date}")
        elif c.source == "manual":
            out.append(f"Conditions: constant T {c.temperature_K:g} K, P {c.pressure_hPa:g} hPa, H2O {c.h2o_percent:g} %, flow {c.flow_m3h:g} m³/h")
    out.append(f"Time: {c.date} {c.start}–{c.end} UTC in {c.step_s} s steps; roof open "
               + (", ".join(f"{o}–{cl}" for o, cl in c.roof_open) or "never (dark)"))
    out.append(f"Chamber: {ch.name}, {ch.volume_m3:g} m³" + (f"; {', '.join(c.constrain_species)} held to measurements" if c.constrain_species else ""))
    for inj in spec.injections:
        label = ACTIONS.get(inj.action_key(), ACTIONS["other"])["label"]
        what = f"{inj.total_ppb} ppb of {inj.products[:50]}" if inj.products else "recorded only"
        extra = f" + {len(inj.extra_channels)} impurity channel(s)" if inj.extra_channels else ""
        mode = inj.delivery.get("mode", "")
        out.append(f"Injection {inj.start or '--:--'} {label}: {what}{extra}" + (f" [from {mode}]" if mode and mode != "ppb" else "")
                   + (" (delayed release)" if inj.chain else ""))
    on = [t for k, t, _d in CHAMBER_ITEMS if getattr(spec.chamber, k)]
    out.append("Chamber effects: " + (", ".join(t.split(" (")[0] for t in on) or "none"))
    out.append(f"Initial values: {len(spec.initials)} species set, all others at the floor value")
    if spec.overrides:
        out.append(f"Expert overrides: {len(spec.overrides)} statement(s) appended")
    return out


def _boxed(page: QWidget, title: str = "") -> QWidget:
    """A wizard page shown as a plain section: title as group box, page inside."""
    box = QGroupBox(title or (page.title() if hasattr(page, "title") else ""))
    lay = QVBoxLayout(box); lay.setContentsMargins(8, 10, 8, 6)
    sub = page.subTitle() if hasattr(page, "subTitle") else ""
    if sub:
        lay.addWidget(_info(sub))
    lay.addWidget(page)
    return box


def _scroll(*widgets: QWidget) -> QScrollArea:
    inner = QWidget(); lay = QVBoxLayout(inner); lay.setContentsMargins(0, 0, 0, 0)
    for w in widgets:
        lay.addWidget(w)
    lay.addStretch()
    sa = QScrollArea(); sa.setWidgetResizable(True); sa.setWidget(inner); sa.setFrameShape(QFrame.Shape.NoFrame)
    sa.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)      # pages fit the column width
    return sa


def _open_path(path: str):
    try:
        os.startfile(path)                      # Windows
    except AttributeError:                      # pragma: no cover
        subprocess.Popen(["xdg-open", path])


def _const(expr: str, default: float = 0.0) -> float:
    """Numeric value of a small EASY constant expression such as '2.6*0.7' (no names)."""
    try:
        return float(eval(expr.replace("@", "**").replace("D", "e").replace("d", "e"), {"__builtins__": {}}, {}))
    except Exception:
        return default


# ─────────────────────────────────────────────────────────────────────────────
class InputsPlot(QWidget):
    """Matplotlib view of the inputs as the model will see them.  Modes:
    'conditions' (T + roof, jNO2 + H2O), 'injections' (pulse timeline),
    'chamber' (what the chamber adds per hour), 'all' (the three stacked)."""

    def __init__(self, win: "JuiceWindow", mode: str, height: float = 3.2):
        super().__init__()
        self.win, self.mode = win, mode
        lay = QVBoxLayout(self); lay.setContentsMargins(0, 0, 0, 0)
        self.fig = Figure(figsize=(6, height), tight_layout=True); self.canvas = FigureCanvas(self.fig)
        self.canvas.setMinimumHeight(int(height * 95)); lay.addWidget(self.canvas, stretch=1)
        try:                                   # as in the results viewer: wheel zoom, drag pan, double-click reset, right-click menu
            from plot_interaction import attach
            attach(self.canvas)
            self.canvas.setToolTip("wheel: zoom the time · Ctrl+wheel or right-drag: zoom the values · drag: pan · double-click: reset · "
                                   "right-click: zoom / log scale menu")
        except Exception as e:
            print(f"[gui] plot interaction not attached: {e}")
        hb = QHBoxLayout(); b = QPushButton("Update plot"); b.clicked.connect(self.plot); hb.addWidget(b)
        self.note = _info({"conditions": "Temperature, roof, photolysis and humidity as the model will see them (archive data is read on update).",
                           "injections": "When each injection pulse is on; the model spreads the total ppb over the pulse.",
                           "chamber": "Chamber contributions per hour from the chosen effects and the conditions.",
                           "all": "Conditions and injections as the model will see them."}[mode])
        hb.addWidget(self.note, stretch=1); lay.addLayout(hb)
        self.plotted = False

    def _inputs(self):
        spec = self.win.collect_spec()
        grid, cols, notes = self.win.inputs_for(spec)
        c = spec.conditions
        h0 = float(c.start.split(":")[0]) + float(c.start.split(":")[1]) / 60.0
        hours = (grid - grid[0]) / 3600.0 + h0
        return spec, hours, cols, notes

    def plot(self):
        self.fig.clear()
        try:
            spec, hours, cols, notes = self._inputs()
        except Exception as e:
            ax = self.fig.add_subplot(111); ax.text(0.5, 0.5, f"nothing to plot yet:\n{e}", ha="center", va="center", wrap=True, fontsize=8)
            self.canvas.draw(); return
        if self.mode == "all":
            axes = [self.fig.add_subplot(311), self.fig.add_subplot(312), self.fig.add_subplot(313)]
            self._conditions(axes[0], axes[1], hours, cols); self._injections(axes[2], hours, cols, spec)
        elif self.mode == "conditions":
            self._conditions(self.fig.add_subplot(121), self.fig.add_subplot(122), hours, cols)
        elif self.mode == "injections":
            obs = self._measured_koh(spec)
            if obs is not None and len(obs[0]):
                gs = self.fig.add_gridspec(2, 1, height_ratios=[1, 1.6])
                self._injections(self.fig.add_subplot(gs[0]), hours, cols, spec)
                koh_note = self._koh_panel(self.fig.add_subplot(gs[1]), spec, hours, cols, obs)
                notes = list(notes) + [koh_note]
            else:
                self._injections(self.fig.add_subplot(111), hours, cols, spec)
                notes = list(notes) + (["no measured kOH for this day — only the pulses are shown"] if spec.conditions.date else [])
        else:
            self._chamber(self.fig.add_subplot(111), spec, hours, cols)
        for ax in self.fig.axes:
            ax.tick_params(labelsize=7)
        if notes:
            self.note.setText("; ".join(n[:110] for n in notes[:3]))
        self.canvas.draw(); self.plotted = True

    def _conditions(self, ax1, ax2, hours, cols):
        if "t" in cols:
            ax1.plot(hours, cols["t"], color=ACCENT, lw=1.2); ax1.set_ylabel("T [K]", color=ACCENT, fontsize=8)
        if "hv" in cols:
            axb = ax1.twinx(); axb.fill_between(hours, 0, cols["hv"], color="#ffd966", alpha=0.4, step="post"); axb.set_ylim(0, 1.05); axb.set_yticks([])
        ax1.set_title("temperature, light on (yellow)", fontsize=8); ax1.set_xlabel("hour (UTC)", fontsize=8)
        jn = [k for k in cols if k.lower() == "jno2"] or [k for k in cols if k.lower().startswith("j")][:1]
        for k in jn:
            ax2.plot(hours, cols[k], lw=1, color="#3b6fc4", label=k)
        ax2.set_ylabel("jNO2 [1/s]", color="#3b6fc4", fontsize=8)
        if "H2O" in cols:
            axc = ax2.twinx(); axc.plot(hours, cols["H2O"], color="#2e9e4f", lw=1, ls="--"); axc.set_ylabel("H2O [%]", color="#2e9e4f", fontsize=8)
        ax2.set_title("photolysis, humidity (dashed)", fontsize=8); ax2.set_xlabel("hour (UTC)", fontsize=8)

    _KOH_OBS: Dict[Tuple[str, str], object] = {}        # (archive root, date) → measured kOH (hours, values) or None, per session

    def _measured_koh(self, spec):
        c = spec.conditions
        if not c.date:
            return None
        key = (c.archive_root, c.date)
        if key not in InputsPlot._KOH_OBS:
            got = None
            try:
                import experiment_plots as ep
                out = self.win.out_dir()
                obs = ep.load_observations(c.archive_root, c.date, keys=["kOH"], cache_dir=out if os.path.isdir(out) else "")
                if "kOH" in obs:
                    got = (np.asarray(obs["kOH"][0], dtype=float), np.asarray(obs["kOH"][1], dtype=float))
            except Exception:
                got = None
            InputsPlot._KOH_OBS[key] = got
        return InputsPlot._KOH_OBS[key]

    def _koh_panel(self, ax, spec, hours, cols, obs):
        import experiment_plots as ep
        oh_h, oh_k = obs
        exp_k, events, bg, note = ep.expected_koh(spec, hours, cols, oh_h, oh_k)
        lo, hi = hours[0], hours[-1]
        sel = (oh_h >= lo) & (oh_h <= hi)
        ax.plot(oh_h[sel], oh_k[sel], ".", ms=2.5, color="#555555", alpha=0.7, label="measured (LP-LIF)")
        ax.plot(hours, exp_k, color=ACCENT, lw=1.6, label="expected from the injections")
        ymax = max(np.nanmax(exp_k) if len(exp_k) else 0, np.nanpercentile(oh_k[sel], 99) if sel.any() else 0) * 1.15 or 1.0
        for name, h0, k in events:
            ax.axvline(h0, color=ACCENT, lw=0.6, alpha=0.35)
            ax.annotate(f"{name[:14]} +{k:.2g}", (h0, ymax * 0.97), rotation=90, fontsize=6, va="top", ha="right", color="#444")
        if "hv" in cols:
            axb = ax.twinx(); axb.fill_between(hours, 0, cols["hv"], color="#ffd966", alpha=0.2, step="post"); axb.set_ylim(0, 1.05); axb.set_yticks([])
        ax.set_ylim(0, ymax); ax.set_xlim(lo, hi)
        ax.set_ylabel("kOH [s⁻¹]", fontsize=8); ax.set_xlabel("hour (UTC)", fontsize=8)
        ax.set_title("OH reactivity: measured vs expected from the injected amounts (dilution only, no chemistry)", fontsize=8)
        ax.legend(fontsize=7, loc="upper left")
        return note

    def _injections(self, ax, hours, cols, spec=None):
        qs = [k for k in cols if k.upper().startswith("Q")]
        for i, k in enumerate(qs):
            y = cols[k]
            if np.any(y > 0):
                ax.fill_between(hours, i, i + 0.8 * (y > 0), step="post", color=ACCENT, alpha=0.6)
        labels = [injection_label(k, spec) for k in qs]       # display only: the files keep the Q column names
        ax.set_yticks([i + 0.4 for i in range(len(qs))]); ax.set_yticklabels(labels, fontsize=7); ax.set_ylim(-0.2, max(1, len(qs)))
        if "hv" in cols:
            axb = ax.twinx(); axb.fill_between(hours, 0, cols["hv"], color="#ffd966", alpha=0.25, step="post"); axb.set_ylim(0, 1.05); axb.set_yticks([])
        ax.set_xlabel("hour (UTC)", fontsize=8); ax.set_title("injection pulses (light on in yellow)", fontsize=8)

    def _chamber(self, ax, spec, hours, cols):
        import experiment_plots as ep
        src = ep.chamber_sources(spec, cols, len(hours))
        colors = {"HONO wall source": ACCENT, "HCHO wall source": "#3b6fc4", "dark HONO": "#888"}
        for name, col in colors.items():
            if name in src:
                ax.plot(hours, src[name], color=col, lw=1.3 if name != "dark HONO" else 1, ls=":" if name == "dark HONO" else "-", label=name)
        if "dilution" in src:
            axb = ax.twinx(); axb.plot(hours, src["dilution"], color="#2e9e4f", lw=1, ls="--"); axb.set_ylabel("dilution [%/h]", color="#2e9e4f", fontsize=8)
        ax.set_ylabel("source [ppb/h]", fontsize=8); ax.set_xlabel("hour (UTC)", fontsize=8)
        ax.set_title("what the chamber adds per hour (dilution dashed)", fontsize=8)
        if ax.get_legend_handles_labels()[0]:
            ax.legend(fontsize=7, loc="upper left")


# ─────────────────────────────────────────────────────────────────────────────
class ReviewPanel(QWidget):
    """Always visible: summary, what goes into the files, validation, the inputs plot, the files."""

    def __init__(self, win: "JuiceWindow"):
        super().__init__()
        self.win = win
        lay = QVBoxLayout(self); lay.setContentsMargins(4, 4, 4, 4)
        hb0 = QHBoxLayout(); hb0.addWidget(QLabel("<b>Review</b> — follows your entries; press Refresh if a plot looks stale"))
        hb0.addStretch(); rb = QPushButton("Refresh"); rb.clicked.connect(lambda: self.win.refresh_all(force_plot=True)); hb0.addWidget(rb)
        lay.addLayout(hb0)
        self.tabs = QTabWidget(); lay.addWidget(self.tabs, stretch=1)
        st = QWidget(); sl = QVBoxLayout(st); sl.setSpacing(4)
        self.summary = QLabel(""); self.summary.setWordWrap(True); self.summary.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        sl.addWidget(self.summary)                                 # one line: name, day, time window, conditions
        # 'Needs attention' and 'What goes into the input files' share the space below; drag the
        # handle between them to give either more room (the split is remembered)
        top = QWidget(); tl = QVBoxLayout(top); tl.setContentsMargins(0, 0, 0, 0); tl.setSpacing(4)
        self.checks_title = QLabel("<b>Needs attention</b>"); tl.addWidget(self.checks_title)
        self.checks = QListWidget(); self.checks.setMinimumHeight(40); self.checks.setWordWrap(True)
        self.checks.setToolTip("click an entry: the step it belongs to opens and the field is pointed out"); tl.addWidget(self.checks, stretch=1)
        self.checks.itemClicked.connect(self._guide_clicked)
        self.fix = QLabel(""); self.fix.setWordWrap(True); self.fix.setStyleSheet("color: #9a6b00; font-size: 10px; padding: 2px 4px;"); self.fix.hide()
        tl.addWidget(self.fix)
        bot = QWidget(); bl = QVBoxLayout(bot); bl.setContentsMargins(0, 4, 0, 0); bl.setSpacing(4)
        bl.addWidget(QLabel("<b>What goes into the input files</b>"))
        self.bullets = QListWidget(); self.bullets.setMinimumHeight(40); self.bullets.setWordWrap(True); bl.addWidget(self.bullets, stretch=1)
        self.split = QSplitter(Qt.Orientation.Vertical); self.split.addWidget(top); self.split.addWidget(bot)
        self.split.setChildrenCollapsible(False); self.split.setHandleWidth(7)
        self.split.setStyleSheet("QSplitter::handle:vertical { background: #d8d8d8; margin: 2px 40px; border-radius: 2px; }"
                                 "QSplitter::handle:vertical:hover { background: #6b3e85; }")
        self.split.handle(1).setToolTip("drag to give 'Needs attention' or 'What goes into the input files' more room")
        self.split.handle(1).setCursor(Qt.CursorShape.SplitVCursor)
        self.split.setStretchFactor(0, 1); self.split.setStretchFactor(1, 2); self.split.setSizes([160, 320])
        try:
            from PyQt6.QtCore import QSettings
            self._settings = QSettings("FZJ", "JUICE")
            state = self._settings.value("review/attention_split")
            if state is not None:
                self.split.restoreState(state)
            self.split.splitterMoved.connect(lambda *_: self._settings.setValue("review/attention_split", self.split.saveState()))
        except Exception:
            pass
        sl.addWidget(self.split, stretch=1)
        self.tabs.addTab(st, "Summary")
        self.plot = InputsPlot(win, "all", height=6.0); self.tabs.addTab(self.plot, "Inputs plot")
        ft = QWidget(); fl = QVBoxLayout(ft)
        fl.addWidget(_info("Generated in the project folder when you press 'Prepare input files' on the Run step. "
                           "The .easy file is the same input FACSIMILE uses; you do not need to read it."))
        self.files = QListWidget(); self.files.setWordWrap(True); self.files.itemDoubleClicked.connect(lambda it: _open_path(it.data(Qt.ItemDataRole.UserRole))); fl.addWidget(self.files, stretch=1)
        hb2 = QHBoxLayout()
        fb = QPushButton("Open project folder"); fb.clicked.connect(lambda: _open_path(self.win.out_dir()) if os.path.isdir(self.win.out_dir()) else None); hb2.addWidget(fb)
        ob = QPushButton("Open selected file"); ob.clicked.connect(self._open_selected); hb2.addWidget(ob); hb2.addStretch(); fl.addLayout(hb2)
        self.tabs.addTab(ft, "Files")
        self.tabs.currentChanged.connect(lambda i: self.plot.plot() if i == 1 and not self.plot.plotted else None)

    def refresh(self):
        try:
            spec = self.win.collect_spec()
        except Exception as e:
            self.summary.setText(f"<i>cannot read the pages: {e}</i>"); return []
        self.checks.clear(); self.bullets.clear(); self.plot.plotted = False
        if is_blank(spec):
            self.summary.setText("<i>Nothing entered yet.</i><br>Give the experiment a name, then open a saved experiment, "
                                 "import the SAPHIR wiki page, paste your notes, or pick a template.")
            self.refresh_files(); return []
        known = known_species_for(spec, self.win.include_dir, self.win.out_dir())
        msgs = validate(spec, known)
        c = spec.conditions
        self.summary.setText(f"<b>{spec.name or '(no name)'}</b> — {c.date or '(no date)'} {c.start}–{c.end} UTC, "
                             f"{'measured conditions from the archive' if c.source == 'archive' else 'manual conditions'}"
                             f"{', step ' + str(c.step_s) + ' s' if getattr(c, 'step_s', None) else ''}")
        self.fix.hide()
        n_err = sum(1 for lvl, _m in msgs if lvl == "error")
        n_warn = sum(1 for lvl, _m in msgs if lvl == "warning")
        self.checks_title.setText("<b>Needs attention</b>" + (f" — {n_err} to fix before running, {n_warn} to check" if n_err + n_warn else " — nothing"))
        if not n_err + n_warn:
            self.checks.addItem(QListWidgetItem("✔ No problems found"))
        for lvl, m in sorted(msgs, key=lambda x: {"error": 0, "warning": 1}.get(x[0], 2)):   # notes last
            it = QListWidgetItem(level_icon(lvl) + m); it.setData(Qt.ItemDataRole.UserRole, (lvl, m))
            if lvl == "info":
                it.setForeground(QColor("#666666")); it.setToolTip("a note on this experiment, nothing to fix")
            else:
                it.setToolTip("click: go to the step and see how to fix it")
            self.checks.addItem(it)
        for b in spec_bullets(spec):
            self.bullets.addItem(QListWidgetItem("• " + b))
        self.refresh_files()
        return msgs

    # ── guidance: a clicked warning opens its step and points at the field ──
    def _guide(self, msg: str):
        """(step index, how to fix, action) for a validation message."""
        import re as _re
        w = self.win; c, s0, ch, inj = w.page_cond, w.page_start, w.page_chem, w.page_inj

        def focus(widget):
            def act():
                try:
                    widget.setFocus(); widget.selectAll() if hasattr(widget, "selectAll") else None
                except Exception:
                    pass
            return act

        def inj_row(name):
            def act():
                t = inj.table
                for r in range(t.rowCount()):
                    it = t.item(r, 1)
                    if it is not None and it.text().strip() == name:
                        t.selectRow(r); t.scrollToItem(it); t.setFocus(); return
            return act

        if msg.startswith("Import:"):
            m = _re.search(r"^Import: (\d\d:\d\d): ([^:]+): no MCM species and no surrogate", msg)
            if m:
                return (2, f"'{m.group(2)}' (injected {m.group(1)}) is in the mixture but has no MCM species. Open the mixture on the "
                           f"Injections step ('+ mixture …' or double-click its row) and give it a surrogate — the suggestion in the "
                           f"message is a start — or leave it unmodelled on purpose. The note stays as a record either way.",
                        None)
            return (2, "A remark from the wiki import: check the injection it refers to; the note stays as a record of what the import assumed.", None)
        m = _re.search(r"Injection '([^']+)': product '([^']+)' is not a species", msg)
        if m:
            name, sp = m.group(1), m.group(2)

            def act():
                ch.mode_compounds.setChecked(True); ch.search.setText(sp); ch.search.setFocus()
            return (1, f"'{sp}' is injected but not in the mechanism. Either add it here — the search is prefilled, pick it and press "
                       f"'Add' — or go to Injections and replace it in '{name}' by a surrogate that is in the mechanism "
                       f"(the mixture dialog's formula search suggests one).", act)
        m = _re.search(r"Injection '([^']+)': start time", msg)
        if m:
            return (2, f"Enter the start of '{m.group(1)}' as HH:MM (UTC), e.g. 08:30.", inj_row(m.group(1)))
        m = _re.search(r"Injection '([^']+)': amount", msg)
        if m:
            return (2, f"'{m.group(1)}' uses an expression instead of a number; fine if intended, otherwise type the total ppb.", inj_row(m.group(1)))
        m = _re.search(r"Injection '([^']+)': neither modelled nor timed", msg)
        if m:
            return (2, f"'{m.group(1)}' has neither a product nor a start time, so it does nothing: give it a product (modelled) "
                       f"or a start (recorded in the ENZ), or remove the row.", inj_row(m.group(1)))
        if msg.startswith("No injections"):
            return (2, "Add the injections of the day: 'Add action' for NO, O3, humidification …, '+ mixture …' for a solution, "
                       "or import them from the SAPHIR wiki page on the Experiment step.", None)
        if msg.startswith("Enter the experiment date") or msg.startswith("Time window") or msg.startswith("Date or time"):
            return (0, "Fill in the date (YYYY-MM-DD) and the start and end time (HH:MM UTC) under 'Where and when'.", focus(c.date))
        if msg.startswith("Measured inputs from a file") or msg.startswith("Measured inputs file"):
            return (0, "Choose the table with the measured conditions (Browse …), or press 'Template …' to get an example with the column names, units and time formats.", focus(c.input_path))
        if msg.startswith("Photolysis: measured J values selected but none found"):
            return (0, "No measured J values in the archive / file for this day: switch 'J values' under Advanced to clear sky, or add jNO2, jO1D, … columns to the input file.", focus(c.jsrc))
        if msg.startswith("Inputs:"):
            return (0, "Some conditions are not in the measurements and stay constant; check the manual values or add the columns to the input file.", focus(c.temp))
        if msg.startswith("Held species without measurements"):
            return (0, "A species can only be held to measurements that exist: add the column (NO, NO2, O3) to the input file or remove it from 'held species' under Advanced.", focus(c.constrain))
        if msg.startswith("Choose where the chamber conditions"):
            return (0, "Tick 'measured (archive)' when the day is on the Z: drive, otherwise 'manual' and enter T, p, humidity and flow.", focus(c.src_manual))
        if msg.startswith("No chemistry"):
            return (1, "Search a compound by name or formula and add it to the chosen compounds, choose an existing include set, "
                       "or tick 'Inorganic chemistry only' for a run without VOC.",
                    lambda: (ch.search.setFocus()))
        if msg.startswith("Measured J values"):
            return (0, "There is no archive data for this day: switch the conditions to manual (clear-sky photolysis), "
                       "or correct the date / archive path.", focus(c.src_manual))
        if msg.startswith(("Roof never opens", "Lights never on")):
            return (0, "Nothing to fix for a dark experiment. If the roof was opened (lights switched on), add the open and "
                       "close times to the roof table.", lambda: (c.roof.setFocus()))
        return (None, "", None)

    def _guide_clicked(self, item):
        data = item.data(Qt.ItemDataRole.UserRole)
        if not data:
            self.fix.hide(); return
        lvl, msg = data
        step, hint, act = self._guide(msg)
        if step is None:
            self.fix.setText("No guidance for this message yet."); self.fix.show(); return
        self.win.steps.setCurrentRow(step)
        if act is not None:
            try:
                act()
            except Exception:
                pass
        self.fix.setText(f"<b>How to fix</b> ({STEPS[step]} step): {hint}"); self.fix.show()

    def refresh_files(self):
        self.files.clear()
        out = self.win.out_dir()
        if not os.path.isdir(out):
            return
        for f in sorted(os.listdir(out)):
            p = os.path.join(out, f)
            label = {".enz": "measured conditions and injection switches (ENZ)", ".easy": "model input shared with FACSIMILE (EASY)",
                     ".json": "scenario for the solver" if not f.endswith(".experiment.json") else "experiment description (this GUI)",
                     ".cs": "results: concentrations", ".kv": "results: rate constants and turnover", ".pd": "results: production and loss",
                     ".elog": "expanded input (all includes inlined)"}.get(os.path.splitext(f)[1], "")
            if os.path.isdir(p):
                label = "chemistry include files"
            it = QListWidgetItem(f + (f"   — {label}" if label else "")); it.setData(Qt.ItemDataRole.UserRole, p); self.files.addItem(it)

    def _open_selected(self):
        it = self.files.currentItem()
        if it:
            _open_path(it.data(Qt.ItemDataRole.UserRole))


# ─────────────────────────────────────────────────────────────────────────────
class RunPage(QWidget):
    """Last step: prepare the input files → run scenario_runner → open the results."""

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
        self.kvpd = QCheckBox("Also write the reaction rates and production/loss budgets "
                              "(large files; needed for the process analysis and reactivity views in the viewer)")
        self.kvpd.setChecked(True); g.addWidget(self.kvpd, r, 0, 1, 4); r += 1
        # Solver profile: the validated ways to run (FR_GUI.toml [solver]); the TOML's profile is preselected.
        g.addWidget(QLabel("Solver"), r, 0); self.solver = QComboBox()
        self._profiles = solver_profiles(self.win.toml_path)
        for name, p in self._profiles.items():
            self.solver.addItem(p["label"], name)
        dflt = default_profile(self.win.toml_path)
        if dflt in self._profiles:
            self.solver.setCurrentIndex(self.solver.findData(dflt))
        else:                                   # a settings file from before profiles: its own rtol / restart mode
            sol = solver_from_toml(self.win.toml_path)
            self.solver.insertItem(0, f"From the settings file (rtol {sol['rtol']:g}, "
                                      f"{'warm' if sol['perstep_warm'] else 'cold restart'})", None)
            self.solver.setCurrentIndex(0)
        self.solver.setToolTip("\n\n".join(f"{p['label']}:\n{p['summary']}" for p in self._profiles.values()))
        g.addWidget(self.solver, r, 1, 1, 3); r += 1
        self.solver_note = QLabel(""); self.solver_note.setWordWrap(True)
        self.solver_note.setStyleSheet("color: #555555;"); g.addWidget(self.solver_note, r, 1, 1, 3); r += 1
        self.solver.currentIndexChanged.connect(self._solver_changed); self._solver_changed()
        lay.addLayout(g)
        lay.addWidget(_info(f"The solver profiles are validated against FACSIMILE (settings in {os.path.basename(self.win.toml_path) or 'the TOML'}). "
                            "Fast is the everyday choice; Precise tracks FACSIMILE 5-10x more tightly for validation runs; "
                            "TWOSTEP needs no Jacobian and is for mechanisms far larger than the MCM (GECKO-A). "
                            "Compare runs only within one profile. Benchmark and tolerance experiments belong in the expert tool, not here."))
        hb = QHBoxLayout()
        self.prep_btn = QPushButton("1. Prepare input files"); self.prep_btn.setObjectName("primary"); self.prep_btn.clicked.connect(self.prepare); hb.addWidget(self.prep_btn)
        self.run_btn = QPushButton("2. Run simulation"); self.run_btn.setObjectName("primary"); self.run_btn.setEnabled(False); self.run_btn.clicked.connect(self.run); hb.addWidget(self.run_btn)
        self.stop_btn = QPushButton("Stop"); self.stop_btn.setEnabled(False); self.stop_btn.clicked.connect(self.stop); hb.addWidget(self.stop_btn)
        self.view_btn = QPushButton("3. Open results"); self.view_btn.setObjectName("primary"); self.view_btn.setEnabled(False); self.view_btn.clicked.connect(self.open_results); hb.addWidget(self.view_btn)
        # an imported FACSIMILE experiment remembers the FACSIMILE results it came with
        self.cmp_btn = QPushButton("Compare with FACSIMILE"); self.cmp_btn.setEnabled(False); self.cmp_btn.clicked.connect(self.compare)
        self.cmp_btn.setToolTip(f"The FACSIMILE results this experiment was imported with, against this run: max relative deviation over "
                                f"the run, ±{COMPARE_TOL*100:g} %; concentrations measured against at least {COMPARE_FLOOR:g} molec/cm³ "
                                "(≈ 1e-3 ppt); settings (J values, rate constants) without a floor")
        hb.addWidget(self.cmp_btn)
        self.study_btn = QPushButton("Sweep this experiment …")
        self.study_btn.setToolTip("Saves this experiment and opens the study window on it: pick which of its injections, sources or "
                                  "initial values to vary over a grid (isopleth, sensitivities).")
        self.study_btn.clicked.connect(self._study_from_here); hb.addWidget(self.study_btn)
        hb.addStretch(); lay.addLayout(hb)
        self.files = QLabel(""); self.files.setWordWrap(True); lay.addWidget(self.files)
        self.bar = QProgressBar(); self.bar.setRange(0, 100); self.bar.setTextVisible(True); self.bar.setFormat("not started"); lay.addWidget(self.bar)
        self.log = QTextEdit(); self.log.setReadOnly(True); self.log.setMinimumHeight(260)
        self.log.setStyleSheet("font-family: Consolas, monospace; font-size: 11px; background: #fbfbfb; color: #222222;")
        lay.addWidget(self.log, stretch=1)
        lay.addWidget(self._earlier_group())

    def refresh(self):
        self.folder.setText(self.win.out_dir())
        self.cmp_btn.setVisible(bool(self.win._base.reference_cs))
        self.cmp_btn.setEnabled(bool(self.win._base.reference_cs) and os.path.isfile(self.win._base.reference_cs) and os.path.isfile(self.cs_path()))

    def _solver_changed(self, *_a):
        name = self.solver.currentData()
        self.solver_note.setText(self._profiles[name]["summary"] if name in self._profiles
                                 else "The rtol and restart mode written in the settings file (no profile set there).")
        self.json_path = ""                     # the prepared JSON carries the solver: prepare again
        if hasattr(self, "run_btn"):
            self.run_btn.setEnabled(False)

    def prepare(self) -> str:
        """Write ENZ/EASY/experiment.json and the scenario JSON.  Returns the JSON path ('' on error)."""
        self.refresh()
        spec = self.win.collect_spec()
        known = known_species_for(spec, self.win.include_dir, self.win.out_dir())
        errors = [m for lvl, m in validate(spec, known) if lvl == "error"]
        if errors:
            self.log.append("Cannot prepare:\n  " + "\n  ".join(errors)); self.bar.setValue(0); self.bar.setFormat("fix the errors first"); return ""
        self.bar.setValue(0); self.bar.setFormat("preparing …"); QApplication.processEvents()    # white while working, green when done
        try:
            res = self.win.write_all()
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):           # the EASY-version line and the parser's warnings go to the log
                self.json_path = prepare_scenario(res, self.win.include_dir, self.kvpd.isChecked(), self.win.toml_path,
                                                  extra_dirs=list(spec.include_dirs), easy_version=spec.easy_version or None,
                                                  solver=solver_from_toml(self.win.toml_path, self.solver.currentData()))
        except Exception as e:
            self.log.append(f"Prepare failed: {e}"); self.bar.setValue(0); self.bar.setFormat("prepare failed"); return ""
        out = os.path.dirname(self.json_path)
        for line in buf.getvalue().splitlines():
            if "[scenario_prep]" in line:
                self.log.append(line)
        missing = []
        try:
            with open(os.path.join(out, spec.name + "_juice.elog"), encoding="utf-8", errors="replace") as fh:
                missing = sorted({m.group(1).strip() for m in re.finditer(r"^;#include (.+?) ! NOT FOUND", fh.read(), re.M)})
        except OSError:
            pass
        if missing:
            self.log.append(f"✖ {len(missing)} include file(s) not found: " + ", ".join(missing) +
                            "\n  Add the folder that holds them under 'Include folders' on the Chemistry step.")
        self.files.setText("Written: " + ", ".join(os.path.basename(p) for p in (res["enz"], res["easy"], res["spec_path"], self.json_path)))
        try:                                   # the measurements the results will be compared with, extracted once
            import experiment_plots as _ep
            note: list = []
            obs = _ep.load_observations(spec.conditions.archive_root, spec.conditions.date, cache_dir=out, status=note)
            if obs:
                self.log.append(f"observations {note[0] if note else ''}: {', '.join(sorted(obs))}")
        except Exception as e:
            self.log.append(f"note: observations not extracted ({e})")
        for n in res.get("notes", []):
            self.log.append("note: " + n)
        self.log.append(f"Inputs ready in {out}")
        self.bar.setValue(100); self.bar.setFormat("inputs ready" if not missing else f"inputs ready — {len(missing)} include(s) missing")
        self.run_btn.setEnabled(True)
        self.view_btn.setEnabled(os.path.isfile(self.cs_path()))
        self.win.review.refresh_files()
        self.win.set_step_state(RUN_STEP, "ok")
        self.refresh()
        return self.json_path

    def compare(self):
        """This run against the FACSIMILE results the experiment was imported with (compare_juice_facsimile)."""
        import compare_juice_facsimile as cjf
        ref = self.win._base.reference_cs
        cs = self.cs_path()
        if not (ref and os.path.isfile(ref) and cs and os.path.isfile(cs)):
            return
        stem = os.path.splitext(os.path.basename(cs))[0]
        work = os.path.dirname(cs); log_path = os.path.join(work, f"{stem}_difference_JUICE_Facsimile.log")
        self.bar.setFormat("comparing …"); QApplication.processEvents()
        buf = io.StringIO()
        try:
            with contextlib.redirect_stdout(buf):
                s = cjf.compare(stem, COMPARE_TOL, no_stop=True, floor=COMPARE_FLOOR, juice_dir_override=work,
                                fac_dir_override=os.path.dirname(ref), fac_stem=os.path.splitext(os.path.basename(ref))[0])
        except Exception as e:
            self.log.append(f"Comparison failed: {e}"); self.bar.setFormat("comparison failed"); return
        with open(log_path, "w", encoding="utf-8") as fh:
            fh.write(buf.getvalue())
        self.log.append(f"JUICE vs FACSIMILE ({os.path.basename(ref)}), ±{COMPARE_TOL*100:g} % (concentrations measured against ≥ {COMPARE_FLOOR:g} molec/cm³):")
        for label, g in (s.get("groups") or {}).items():
            worst = ", ".join(f"{n} {r*100:.2g} %" for r, n in g["worst"][:3])
            self.log.append(f"  {label}: {g['n'] - g['n_exceed']}/{g['n']} within  ·  worst {worst}")
        self.log.append(("✔ the two runs agree" if s.get("ok") else "✖ differences above the tolerance") + f"  (full report: {os.path.basename(log_path)})")
        self.bar.setFormat("agree with FACSIMILE" if s.get("ok") else "differences — see the log")

    def cs_path(self) -> str:
        return self.json_path[:-5] + ".cs" if self.json_path else ""

    def _earlier_group(self) -> QGroupBox:
        """Earlier experiments: look at their results, or sweep / batch them without reopening each one."""
        g = QGroupBox("Earlier experiments"); gl = QVBoxLayout(g)
        gl.addWidget(_info("Work with experiments that were run before, without loading them into the steps above."))
        hb = QHBoxLayout()
        b = QPushButton("Open earlier results in the viewer …"); b.clicked.connect(lambda: self.win._open_results_file()); hb.addWidget(b)
        b = QPushButton("Sweep an earlier experiment …"); b.setToolTip("Choose a saved .experiment.json: its injections, sources and initial values can then be varied over a grid.")
        b.clicked.connect(lambda: self._study_earlier("sweep")); hb.addWidget(b)
        b = QPushButton("Run several earlier experiments with one change …"); b.setToolTip("Choose saved experiments; one change (rate, add-on, chamber parameter, scaled injections, MCM version) is applied to all of them.")
        b.clicked.connect(lambda: self._study_earlier("batch")); hb.addWidget(b)
        b = QPushButton("Open an earlier study …"); b.setToolTip("A saved .study.json: its finished runs, results table and plots, without running anything again.")
        b.clicked.connect(lambda: self._study_earlier("open")); hb.addWidget(b)
        hb.addStretch(); gl.addLayout(hb)
        return g

    def _study_earlier(self, mode: str):
        start = self.win.output_root if os.path.isdir(self.win.output_root) else self.win.easy_root
        if mode == "open":
            p, _ = QFileDialog.getOpenFileName(self.win, "Open a study", start, "Study (*.study.json)")
            if p:
                self.win._study(study=p)
            return
        if mode == "sweep":
            p, _ = QFileDialog.getOpenFileName(self.win, "Experiment to sweep", start, "Experiment (*.experiment.json)")
            if p:
                self.win._study(p, mode="sweep")
        else:
            files, _ = QFileDialog.getOpenFileNames(self.win, "Experiments to run with one change", start, "Experiment (*.experiment.json)")
            if files:
                self.win._study(files, mode="batch")

    def _study_from_here(self):
        """Save the experiment description, then open the study window with it as the base experiment."""
        spec = self.win.collect_spec()
        path = os.path.join(self.win.out_dir(), spec.name + ".experiment.json")
        if not os.path.isdir(os.path.join(self.win.out_dir(), "chemistry")) and not spec.inorganic_only:
            if QMessageBox.question(self.win, "Study", "A study needs this experiment's input files prepared once (step 1).\n"
                                    "Prepare them now and open the study?") != QMessageBox.StandardButton.Yes:
                return
            if not self.prepare():
                return
        self.win._save(path); self.win._study(path, mode="sweep")

    def run(self):
        if not self.json_path:
            return
        # no --gui: that flag makes the runner spawn the viewer itself; here "Open results" does it
        args = ["-u", os.path.join(HERE, "scenario_runner.py"), self.json_path]
        if self.win.toml_path:
            args += ["--toml", self.win.toml_path]
        self._proc = QProcess(self)
        env = QProcessEnvironment.systemEnvironment()      # the runner prints '·' and '—': make it write UTF-8
        env.insert("PYTHONIOENCODING", "utf-8"); env.insert("PYTHONUTF8", "1")
        self._proc.setProcessEnvironment(env)
        self._proc.setProcessChannelMode(QProcess.ProcessChannelMode.MergedChannels)
        self._proc.readyReadStandardOutput.connect(self._on_output)
        self._proc.finished.connect(self._on_finished)
        self._t0 = time.monotonic()
        self.bar.setValue(0); self.bar.setFormat("running   0%")
        self.run_btn.setEnabled(False); self.prep_btn.setEnabled(False); self.stop_btn.setEnabled(True); self.view_btn.setEnabled(False)
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
        try:                                   # <name>.run.json: the study window uses the time as its first estimate
            info = {"exit": code, "elapsed": time.monotonic() - self._t0, "finished": time.strftime("%Y-%m-%d %H:%M:%S"),
                    "cs_bytes": os.path.getsize(self.cs_path()) if ok else 0}
            with open(self.json_path[:-5] + ".run.json", "w", encoding="utf-8") as fh:
                json.dump(info, fh, indent=1)
        except Exception:
            pass
        self.bar.setValue(100 if ok else 0)
        self.bar.setFormat(f"done   {mm}m {ss:02d}s" if ok else f"failed (exit {code})   {mm}m {ss:02d}s")
        self.run_btn.setEnabled(True); self.prep_btn.setEnabled(True); self.stop_btn.setEnabled(False); self.view_btn.setEnabled(ok)
        self.win.review.refresh_files()
        self.win.set_step_state(RUN_STEP, "ok" if ok else "error")
        self.refresh()
        if ok and self.cmp_btn.isEnabled():
            self.log.append("The FACSIMILE results this experiment came with are at hand: 'Compare with FACSIMILE' sets the two runs against each other.")

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
        self._enz_cache: Tuple[str, tuple] = ("", ())
        self.setWindowTitle("JUICE — Jülich Utility for Integrating Chemical Equations")
        self.resize(1550, 920)
        self.setWindowIcon(app_icon())

        # ── left: logo + steps ──
        left = QWidget(); ll = QVBoxLayout(left); ll.setContentsMargins(6, 6, 6, 6); left.setFixedWidth(220)
        if os.path.isfile(LOGO):
            lg = QLabel(); lg.setPixmap(QPixmap(LOGO).scaledToWidth(150, Qt.TransformationMode.SmoothTransformation)); lg.setAlignment(Qt.AlignmentFlag.AlignCenter); ll.addWidget(lg)
        t = QLabel("JUICE"); t.setObjectName("title"); t.setAlignment(Qt.AlignmentFlag.AlignCenter); ll.addWidget(t)
        st = QLabel("chamber experiment → model run"); st.setObjectName("subtitle"); st.setAlignment(Qt.AlignmentFlag.AlignCenter); ll.addWidget(st)
        self.steps = QListWidget(); self.steps.setObjectName("steps")
        for i, s in enumerate(STEPS):
            self.steps.addItem(QListWidgetItem(f"{i + 1}   {s}"))
        self._state = ["todo"] * len(STEPS)
        self.steps.currentRowChanged.connect(self._go)
        ll.addWidget(self.steps, stretch=1)

        # ── middle: pages, each with its plot at the bottom ──
        self.plot_cond = InputsPlot(self, "conditions", 2.8)
        self.plot_inj = InputsPlot(self, "injections", 4.4)
        self.plot_chamber = InputsPlot(self, "chamber", 2.8)
        self.stack = QStackedWidget()
        self.stack.addWidget(_scroll(_boxed(self.page_start, "Experiment"), _boxed(self.page_cond, "Where and when"),
                                     _boxed(self.plot_cond, "Conditions the model will see")))
        self.stack.addWidget(_scroll(_boxed(self.page_chem, "Chemistry")))
        self.stack.addWidget(_scroll(_boxed(self.page_inj, "Injections"), _boxed(self.plot_inj, "Injection timeline")))
        self.stack.addWidget(_scroll(_boxed(self.page_chamber, "What the chamber adds on its own"), _boxed(self.plot_chamber, "Chamber contributions")))
        self.page_run = RunPage(self)
        self.stack.addWidget(_scroll(_boxed(self.page_run, "Run the model")))
        nav = QHBoxLayout()
        self.back_btn = QPushButton("◀ Back"); self.back_btn.clicked.connect(lambda: self._step(-1)); nav.addWidget(self.back_btn)
        nav.addStretch()
        self.next_btn = QPushButton("Next ▶"); self.next_btn.setObjectName("primary"); self.next_btn.clicked.connect(lambda: self._step(+1)); nav.addWidget(self.next_btn)
        mid = QWidget(); ml = QVBoxLayout(mid); ml.setContentsMargins(0, 0, 0, 0); ml.addWidget(self.stack, stretch=1); ml.addLayout(nav)

        # ── right: review ──
        self.review = ReviewPanel(self); self.review.setMinimumWidth(360); self.review.setMaximumWidth(560)
        self.page_review.hide()                  # keeps the overrides text for collect_spec(); edited via the expert tool

        split = QSplitter(); split.addWidget(left); split.addWidget(mid); split.addWidget(self.review)
        split.setStretchFactor(0, 0); split.setStretchFactor(1, 5); split.setStretchFactor(2, 2)
        split.setSizes([220, 950, 420]); split.setChildrenCollapsible(False)
        self.setCentralWidget(split)

        # ── menu ──
        m = self.menuBar().addMenu("File")
        m.addAction("New experiment (blank)", QKeySequence("Ctrl+N"), self._new)
        m.addAction("Open experiment (.experiment.json) …", QKeySequence("Ctrl+O"), self._open)
        m.addAction("Import a FACSIMILE experiment (.easy + ENZ) …", lambda: self.import_from_facsimile(self))
        m.addAction("Save experiment", QKeySequence("Ctrl+S"), lambda: self._save())
        m.addAction("Save experiment as …", self._save_as)
        m.addAction("Open experiment folder", self._open_folder)
        m.addSeparator()
        m.addAction("Open results of an earlier run in the viewer …", self._open_results_file)
        m.addAction("Study: several experiments at once, or a sweep …", self._study)
        m.addSeparator()
        m.addAction("Settings …", self._settings)
        m.addSeparator()
        m.addAction("Quit", self.close)
        m.setToolTipsVisible(True)
        for a in m.actions():
            if a.text() == "Save experiment":
                a.setToolTip("Saves what you entered as <name>.experiment.json in the experiment folder.\n"
                             "The model input files are written by 'Prepare input files' in the Run step.")

        self.load_spec(blank_spec())
        self.page_start.name_edit.setPlaceholderText("e.g. isoprene_2024-07-09")
        self.steps.setCurrentRow(0)
        self._wire_live_refresh()

    # ── live refresh: every edit updates review, step marks, missing-field marks and the plot ──
    def _wire_live_refresh(self):
        self._refresh_timer = QTimer(self); self._refresh_timer.setSingleShot(True); self._refresh_timer.setInterval(500)
        self._refresh_timer.timeout.connect(self.refresh_all)
        c, s0, ch, inj, cb = self.page_cond, self.page_start, self.page_chem, self.page_inj, self.page_chamber
        for w in (c.date, c.start, c.end, c.archive, c.temp, c.press, c.h2o, c.flow, c.constrain, s0.name_edit, s0.out_edit,
                  cb.hono_factor, cb.y_inline, cb.dark_factor):
            w.textChanged.connect(self._schedule_refresh)
        for w in (c.src_archive, c.src_manual, c.src_file):
            w.toggled.connect(self._schedule_refresh)
        c.input_path.textChanged.connect(self._schedule_refresh)
        for w in (c.step,):
            w.valueChanged.connect(self._schedule_refresh)
        for w in (c.jsrc, c.temp_mode, c.chamber, ch.base_combo, ch.mech_combo):
            w.currentIndexChanged.connect(self._schedule_refresh)
        for t in (c.roof, inj.table, inj.chain, ch.adjust, self.page_init.table):
            t.itemChanged.connect(self._schedule_refresh); t.model().rowsRemoved.connect(self._schedule_refresh); t.model().rowsInserted.connect(self._schedule_refresh)
        ch.chosen.model().rowsInserted.connect(self._schedule_refresh); ch.chosen.model().rowsRemoved.connect(self._schedule_refresh)
        for w in (ch.mode_compounds, ch.mode_existing, ch.mode_inorganic):
            w.toggled.connect(self._schedule_refresh)
        for w in list(cb.checks.values()):
            w.toggled.connect(self._schedule_refresh)
        ch.addon_combo.changed.connect(self._schedule_refresh)

    def _schedule_refresh(self, *_a):
        self._refresh_timer.start()

    def refresh_all(self, force_plot: bool = False):
        msgs = self.review.refresh()
        self._update_states(msgs or [])
        self._mark_missing()
        if self.steps.currentRow() == 2:                 # kOH column follows the temperature typed on the Experiment step
            self.page_inj.refresh_koh()
        plot = {0: self.plot_cond, 2: self.plot_inj, 3: self.plot_chamber}.get(self.steps.currentRow())
        if plot is not None and not is_blank(self.collect_spec()):
            plot.plot()
        if force_plot:
            self.review.plot.plot()

    def _mark_missing(self):
        """Soft yellow on what the experiment still needs: name, date, time window, conditions source, chemistry."""
        c, s0, ch = self.page_cond, self.page_start, self.page_chem
        theme.mark_missing(s0.name_edit, not s0.name_edit.text().strip())
        theme.mark_missing(c.date, not c.date.text().strip())
        theme.mark_missing(c.start, not c.start.text().strip())
        theme.mark_missing(c.end, not c.end.text().strip())
        no_src = not (c.src_archive.isChecked() or c.src_manual.isChecked() or c.src_file.isChecked())
        theme.mark_missing(c.src_archive, no_src); theme.mark_missing(c.src_manual, no_src); theme.mark_missing(c.src_file, no_src)
        theme.mark_missing(c.input_path, c.src_file.isChecked() and not c.input_path.text().strip())
        no_chem = (ch.mode_compounds.isChecked() and ch.chosen.count() == 0) or (ch.mode_existing.isChecked() and not ch.mech_combo.currentText().strip())
        theme.mark_missing(ch.chosen, no_chem and ch.mode_compounds.isChecked())
        theme.mark_missing(ch.mech_combo, no_chem and ch.mode_existing.isChecked())

    def load_spec(self, spec: ExperimentSpec, keep_name: bool = False):
        """A template, an import or an opened file arrived: fill the pages, then redraw
        the review and the plot of the step on screen."""
        super().load_spec(spec, keep_name)
        if hasattr(self, "review"):
            self.refresh_all()

    # ── inputs cache for the plots (archive reads are slow) ──
    def inputs_for(self, spec: ExperimentSpec):
        key = json.dumps({"c": spec.conditions.__dict__, "i": [(i.name, i.start, i.duration_min, i.column) for i in spec.injections]},
                         default=str, sort_keys=True)
        if self._enz_cache[0] != key:
            self._enz_cache = (key, enz_columns(spec))
        return self._enz_cache[1]

    # ── navigation ──
    def _go(self, row: int):
        if row < 0:
            return
        self.stack.setCurrentIndex(row)
        self.back_btn.setEnabled(row > 0); self.next_btn.setEnabled(row < RUN_STEP)
        if row == RUN_STEP:
            self.page_run.refresh()
        self.refresh_all()

    def _step(self, d: int):
        row = self.steps.currentRow()
        if d > 0:
            if row == 0 and not self.page_start.validatePage():
                return
            if row == 1 and not self.page_chem.validatePage():
                return
        self.steps.setCurrentRow(max(0, min(RUN_STEP, row + d)))

    def set_step_state(self, i: int, state: str):
        self._state[i] = state
        self._paint_states()

    def _update_states(self, msgs):
        if not msgs and is_blank(self.collect_spec()):
            self._state[:RUN_STEP] = ["todo"] * RUN_STEP; self._paint_states(); return
        txt = " ".join(m for lvl, m in msgs if lvl != "info").lower()      # notes (e.g. a dark experiment) mark no step
        self._state[0] = "error" if any(k in txt for k in ("date", "time window", "conditions come from")) else ("warn" if "j values" in txt or "roof" in txt else "ok")
        self._state[1] = "error" if "no chemistry" in txt else "ok"
        self._state[2] = "warn" if "injection" in txt else "ok"
        self._state[3] = "ok"
        self._paint_states()

    def _paint_states(self):
        marks = {"todo": "", "ok": "  ✔", "warn": "  ⚠", "error": "  ✖"}
        for i, s in enumerate(STEPS):
            self.steps.item(i).setText(f"{i + 1}   {s}" + marks.get(self._state[i], ""))

    # ── file menu ──
    def _new(self):
        self.wiki_memory = None                        # a new experiment starts without the last import
        self.load_spec(blank_spec()); self.steps.setCurrentRow(0); self.review.refresh()

    def _open(self, p: str = ""):
        if not p:
            p, _ = QFileDialog.getOpenFileName(self, "Open experiment", self.easy_root, "Experiment (*.experiment.json)")
        if not p:
            return
        try:
            spec = load_spec(p)
        except Exception as e:
            QMessageBox.critical(self, "Open failed", str(e)); return
        self.wiki_memory = None                        # another experiment: forget the last import
        self.load_spec(spec)
        self.page_start.out_edit.setText(os.path.dirname(os.path.dirname(p)) or self.easy_root)
        self.steps.setCurrentRow(0); self.review.refresh()

    def _save(self, path: str = ""):
        """The experiment description only (<name>.experiment.json) — works for an unfinished experiment too.
        The model input files (.easy, .enz, chemistry, scenario JSON) are written by 'Prepare input files'."""
        try:
            spec = self.collect_spec()
            path = path or os.path.join(self.out_dir(), spec.name + ".experiment.json")
            os.makedirs(os.path.dirname(path), exist_ok=True)
            save_spec(spec, path)
            self.statusBar().showMessage(f"saved {path}", 8000); self.review.refresh_files()
        except Exception as e:
            QMessageBox.critical(self, "Save failed", str(e))

    def _save_as(self):
        spec = self.collect_spec()
        start = os.path.join(self.out_dir(), spec.name + ".experiment.json")
        p, _ = QFileDialog.getSaveFileName(self, "Save experiment as", start, "Experiment (*.experiment.json)")
        if p:
            self._save(p if p.endswith(".experiment.json") else os.path.splitext(p)[0] + ".experiment.json")

    def _open_folder(self):
        if os.path.isdir(self.out_dir()):
            _open_path(self.out_dir())
        else:
            QMessageBox.information(self, "Experiment folder", f"{self.out_dir()} does not exist yet — save the experiment first.")

    def _study(self, experiment=None, mode: str = "", study: str = ""):
        """The study window (study_window.py): sweep one experiment, or run several with one change."""
        args = [os.path.join(HERE, "study_window.py")] + (["--toml", self.toml_path] if self.toml_path else [])
        if study:
            args.append(study)
        if mode:
            args += ["--mode", mode]
        for e in ([experiment] if isinstance(experiment, str) else (experiment or [])):
            if e:
                args += ["--experiment", e]
        p = QProcess(self); p.start(sys.executable, args); self._studies = getattr(self, "_studies", []) + [p]

    def _settings(self):
        """Where JUICE reads and writes, and the solver every run uses — from the settings file (TOML)."""
        d = QDialog(self); d.setWindowTitle("Settings"); d.setMinimumWidth(680)
        v = QVBoxLayout(d)
        v.addWidget(_info("Read from the settings file when the window opens. Edit the file and reopen JUICE to change them."))
        f = QFormLayout(); v.addLayout(f)
        paths = load_paths(self.toml_path)

        def row(label, value, tip=""):
            lb = QLabel(value or "—"); lb.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse); lb.setWordWrap(True)
            if tip:
                lb.setToolTip(tip)
            f.addRow(label, lb)
        row("Settings file", self.toml_path or "none (built-in defaults)")
        row("New experiments go to", paths.get("output_root", ""), "each experiment gets its own folder <here>/<name>")
        row("Mechanism include sets", self.include_dir, "the choices under 'Existing mechanism include set' on the Chemistry step")
        try:
            row("Base mechanisms", ", ".join(mreg.load_bases()))
        except Exception as e:
            row("Base mechanisms", f"not readable: {e}")
        sol = solver_from_toml(self.toml_path)
        row("Solver", f"{sol['solver_method']} ({sol['jac_mode']} Jacobian, {sol['ode_backend']}), "
                      f"{sol['interpolation']} inputs, first step {sol['first_step']:g} s")
        row("Default solver profile", f"{sol['solver_profile']}: rtol {sol['rtol']:g}, "
                                      f"{'warm continuation' if sol['perstep_warm'] else 'cold restart per ENZ row'} "
                                      "(the Run step can pick the other profile per run)")
        row("Python", sys.executable)
        bb = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        if self.toml_path:
            ob = bb.addButton("Open settings file", QDialogButtonBox.ButtonRole.ActionRole); ob.clicked.connect(lambda: _open_path(self.toml_path))
        bb.rejected.connect(d.reject); v.addWidget(bb)
        d.exec()

    def _open_results_file(self, path: str = ""):
        """Any result file (.cs), from this or an earlier experiment, in its own viewer window."""
        if not path:
            start = self.out_dir() if os.path.isdir(self.out_dir()) else (self.output_root if os.path.isdir(self.output_root) else self.easy_root)
            path, _ = QFileDialog.getOpenFileName(self, "Open results", start, "JUICE / FACSIMILE results (*.cs);;All files (*)")
        if not path:
            return
        args = [os.path.join(HERE, "results_viewer.py"), path]
        if self.toml_path:
            args += ["--toml", self.toml_path]
        p = QProcess(self); p.start(sys.executable, args)
        self._viewers = getattr(self, "_viewers", []) + [p]


def main(argv=None):
    argv = list(sys.argv if argv is None else argv)
    theme.taskbar_identity()                 # same "FZJ.JUICE" identity as the viewer and the expert tools: one taskbar icon for all JUICE windows
    app = QApplication(argv)
    theme.apply_theme(app)
    w = JuiceWindow(pick_toml(argv))
    w.show()
    theme.force_taskbar_icon(w)              # the taskbar reads the icon from the native window
    exp = next((a for a in argv[1:] if a.endswith(".experiment.json") and os.path.isfile(a)), "")
    if exp:                                  # from the start menu or a file dropped on the icon
        w._open(exp)
    fac = next((a for a in argv[1:] if a.lower().endswith((".easy", ".elog")) and os.path.isfile(a)), "")
    if fac:                                  # a FACSIMILE experiment: the import dialog opens on it
        QTimer.singleShot(0, lambda: w.import_from_facsimile(w, fac))
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
