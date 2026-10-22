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
and can be opened from the review panel, but never need to be.  FR_GUI.py remains
the expert tool for file-level editing.

    python juice_gui.py                 # paths from GG_GUI.toml if present, else FR_GUI.toml
    python juice_gui.py -t FR_GUI.toml  # another TOML

The step pages are the ones of experiment_wizard.py (shared through ExperimentHost).
Solver numerics come from the TOML (solver_from_toml), never from this window.
"""
from __future__ import annotations

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
                             QCheckBox, QProgressBar, QTextEdit, QMessageBox, QFileDialog, QSplitter, QTabWidget, QFrame)
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

from experiment_spec import load_spec, validate, mechanism_species, enz_columns, ExperimentSpec, ACTIONS
from experiment_wizard import ExperimentHost, ReviewPage, CHAMBER_ITEMS, _info
from elog_generator import generate_both
from idl_parser import parse_idl_easy
import injection_calc as ic
import juice_theme as theme

HERE = os.path.dirname(os.path.abspath(__file__))
LOGO = theme.LOGO
app_icon = theme.app_icon
STEPS = ["Experiment", "Chemistry", "Injections", "Chamber", "Run"]
RUN_STEP = len(STEPS) - 1
ACCENT, ACCENT_DARK, ACCENT_SOFT = theme.ACCENT, theme.ACCENT_DARK, theme.ACCENT_SOFT


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


def solver_from_toml(toml_path: str) -> dict:
    """The solver keys scenario_runner reads from the scenario JSON (not from the TOML),
    taken from the TOML's [solver] section so every run uses the validated set.
    Fallbacks are the validated values, never the runner's own (python / linear / 1e-2)."""
    sol: dict = {}
    for p in (toml_path, toml_path[:-5] + ".local.toml" if toml_path.endswith(".toml") else ""):
        if p and os.path.isfile(p):
            try:
                with open(p, "rb") as fh:
                    sol.update(tomllib.load(fh).get("solver", {}))
            except Exception as e:
                print(f"[juice_gui] could not read {p}: {e}")

    def _f(key, default):
        try:
            return float(sol.get(key, default))
        except (TypeError, ValueError):
            return default
    return {"rtol": _f("rtol_default", 1e-4),
            "ode_backend": str(sol.get("ode_backend", "numba")),
            "solver_method": str(sol.get("solver_method", "BDF")),
            "jac_mode": str(sol.get("jac_mode", "analytic")),
            "interpolation": str(sol.get("interpolation", "piecewise")),
            "first_step": _f("first_step", 1e-6)}


def prepare_scenario(res: dict, include_dir: str, export_kv_pd: bool = True, toml_path: str = "") -> str:
    """EASY + ENZ (from write_all) → canonical elog → scenario JSON next to them.
    Returns the JSON path.  The JSON stem = experiment name, so the runner writes <name>.cs."""
    out = os.path.dirname(res["easy"])
    raw, canon = generate_both(res["easy"], search_dirs=[out, os.path.join(out, "chemistry"), include_dir], output_dir=out)
    cfg = parse_idl_easy(canon, enz_path=res["enz"], chemistry_folder=out,
                         include_fallbacks=[os.path.join(out, "chemistry"), include_dir])
    cfg.update({"folder": out, "enz_path": res["enz"], "elog_path": canon, "export_kv_pd": bool(export_kv_pd)})
    cfg.update(solver_from_toml(toml_path))
    name = os.path.splitext(os.path.basename(res["easy"]))[0]
    path = os.path.join(out, name + ".json")
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(cfg, fh, indent=1, default=str)
    return path


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
    else:
        out.append("Chemistry: not chosen yet")
    for rxn, rate, _cmt in spec.reaction_adjustments:
        out.append(f"Reaction adjusted: {rxn}  →  " + ("switched off" if (rate or "0").strip() in ("0", "0.", "0.0") else f"k = {rate}"))
    if c.source == "archive":
        out.append(f"Conditions: measured T, P, humidity, flow and J values from the archive for {c.date}"
                   + (", clear-sky J" if c.j_source == "clear_sky" else ""))
    elif c.source == "manual":
        out.append(f"Conditions: constant T {c.temperature_K:g} K, P {c.pressure_hPa:g} hPa, H2O {c.h2o_percent:g} %, "
                   f"flow {c.flow_m3h:g} m³/h, {'clear-sky' if c.j_source == 'clear_sky' else 'measured'} photolysis")
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
            self._conditions(axes[0], axes[1], hours, cols); self._injections(axes[2], hours, cols)
        elif self.mode == "conditions":
            self._conditions(self.fig.add_subplot(121), self.fig.add_subplot(122), hours, cols)
        elif self.mode == "injections":
            self._injections(self.fig.add_subplot(111), hours, cols)
        else:
            self._chamber(self.fig.add_subplot(111), spec, hours, cols)
        for ax in self.fig.axes:
            ax.tick_params(labelsize=7)
        if notes:
            self.note.setText("; ".join(n[:70] for n in notes[:3]))
        self.canvas.draw(); self.plotted = True

    def _conditions(self, ax1, ax2, hours, cols):
        if "t" in cols:
            ax1.plot(hours, cols["t"], color=ACCENT, lw=1.2); ax1.set_ylabel("T [K]", color=ACCENT, fontsize=8)
        if "hv" in cols:
            axb = ax1.twinx(); axb.fill_between(hours, 0, cols["hv"], color="#ffd966", alpha=0.4, step="post"); axb.set_ylim(0, 1.05); axb.set_yticks([])
        ax1.set_title("temperature, roof open (yellow)", fontsize=8); ax1.set_xlabel("hour (UTC)", fontsize=8)
        jn = [k for k in cols if k.lower() == "jno2"] or [k for k in cols if k.lower().startswith("j")][:1]
        for k in jn:
            ax2.plot(hours, cols[k], lw=1, color="#3b6fc4", label=k)
        ax2.set_ylabel("jNO2 [1/s]", color="#3b6fc4", fontsize=8)
        if "H2O" in cols:
            axc = ax2.twinx(); axc.plot(hours, cols["H2O"], color="#2e9e4f", lw=1, ls="--"); axc.set_ylabel("H2O [%]", color="#2e9e4f", fontsize=8)
        ax2.set_title("photolysis, humidity (dashed)", fontsize=8); ax2.set_xlabel("hour (UTC)", fontsize=8)

    def _injections(self, ax, hours, cols):
        qs = [k for k in cols if k.upper().startswith("Q")]
        for i, k in enumerate(qs):
            y = cols[k]
            if np.any(y > 0):
                ax.fill_between(hours, i, i + 0.8 * (y > 0), step="post", color=ACCENT, alpha=0.6)
        ax.set_yticks([i + 0.4 for i in range(len(qs))]); ax.set_yticklabels(qs, fontsize=7); ax.set_ylim(-0.2, max(1, len(qs)))
        if "hv" in cols:
            axb = ax.twinx(); axb.fill_between(hours, 0, cols["hv"], color="#ffd966", alpha=0.25, step="post"); axb.set_ylim(0, 1.05); axb.set_yticks([])
        ax.set_xlabel("hour (UTC)", fontsize=8); ax.set_title("injection pulses (roof open in yellow)", fontsize=8)

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
        st = QWidget(); sl = QVBoxLayout(st)
        self.summary = QLabel(""); self.summary.setWordWrap(True); self.summary.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        sl.addWidget(self.summary)
        self.checks = QListWidget(); self.checks.setMaximumHeight(120); self.checks.setWordWrap(True); sl.addWidget(self.checks)
        sl.addWidget(QLabel("<b>What goes into the input files</b>"))
        self.bullets = QListWidget(); self.bullets.setWordWrap(True); sl.addWidget(self.bullets, stretch=1)
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
        known = mechanism_species(self.win.include_dir, spec.mechanism_stem) if spec.mechanism_stem else None
        msgs = validate(spec, known)
        self.summary.setText("<b>" + (spec.name or "(no name)") + "</b><br>" + ReviewPage._summary(spec).replace("\n", "<br>"))
        for lvl, m in msgs:
            self.checks.addItem(QListWidgetItem(("✖ " if lvl == "error" else "⚠ ") + m))
        if not msgs:
            self.checks.addItem(QListWidgetItem("✔ No problems found"))
        for b in spec_bullets(spec):
            self.bullets.addItem(QListWidgetItem("• " + b))
        self.refresh_files()
        return msgs

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
        lay.addLayout(g)
        lay.addWidget(_info(f"The solver runs with the validated settings from {os.path.basename(self.win.toml_path) or 'the TOML'} "
                            "(the same for every experiment, so runs stay comparable with FACSIMILE). "
                            "Benchmark and tolerance experiments belong in the expert tool, not here."))
        hb = QHBoxLayout()
        self.prep_btn = QPushButton("1. Prepare input files"); self.prep_btn.setObjectName("primary"); self.prep_btn.clicked.connect(self.prepare); hb.addWidget(self.prep_btn)
        self.run_btn = QPushButton("2. Run simulation"); self.run_btn.setObjectName("primary"); self.run_btn.setEnabled(False); self.run_btn.clicked.connect(self.run); hb.addWidget(self.run_btn)
        self.stop_btn = QPushButton("Stop"); self.stop_btn.setEnabled(False); self.stop_btn.clicked.connect(self.stop); hb.addWidget(self.stop_btn)
        self.view_btn = QPushButton("3. Open results"); self.view_btn.setObjectName("primary"); self.view_btn.setEnabled(False); self.view_btn.clicked.connect(self.open_results); hb.addWidget(self.view_btn)
        hb.addStretch(); lay.addLayout(hb)
        self.files = QLabel(""); self.files.setWordWrap(True); lay.addWidget(self.files)
        self.bar = QProgressBar(); self.bar.setRange(0, 100); self.bar.setTextVisible(True); self.bar.setFormat("not started"); lay.addWidget(self.bar)
        self.log = QTextEdit(); self.log.setReadOnly(True); self.log.setMinimumHeight(260)
        self.log.setStyleSheet("font-family: Consolas, monospace; font-size: 11px; background: #fbfbfb; color: #222222;")
        lay.addWidget(self.log, stretch=1)

    def refresh(self):
        self.folder.setText(self.win.out_dir())

    def prepare(self) -> str:
        """Write ENZ/EASY/experiment.json and the scenario JSON.  Returns the JSON path ('' on error)."""
        self.refresh()
        spec = self.win.collect_spec()
        known = mechanism_species(self.win.include_dir, spec.mechanism_stem) if spec.mechanism_stem else None
        errors = [m for lvl, m in validate(spec, known) if lvl == "error"]
        if errors:
            self.log.append("Cannot prepare:\n  " + "\n  ".join(errors)); self.bar.setFormat("fix the errors first"); return ""
        self.bar.setFormat("preparing …"); QApplication.processEvents()
        try:
            res = self.win.write_all()
            self.json_path = prepare_scenario(res, self.win.include_dir, self.kvpd.isChecked(), self.win.toml_path)
        except Exception as e:
            self.log.append(f"Prepare failed: {e}"); self.bar.setFormat("prepare failed"); return ""
        out = os.path.dirname(self.json_path)
        self.files.setText("Written: " + ", ".join(os.path.basename(p) for p in (res["enz"], res["easy"], res["spec_path"], self.json_path)))
        for n in res.get("notes", []):
            self.log.append("note: " + n)
        self.log.append(f"Inputs ready in {out}")
        self.bar.setFormat("inputs ready"); self.run_btn.setEnabled(True)
        self.view_btn.setEnabled(os.path.isfile(self.cs_path()))
        self.win.review.refresh_files()
        self.win.set_step_state(RUN_STEP, "ok")
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
        self.bar.setValue(100 if ok else self.bar.value())
        self.bar.setFormat(f"done   {mm}m {ss:02d}s" if ok else f"failed (exit {code})   {mm}m {ss:02d}s")
        self.run_btn.setEnabled(True); self.prep_btn.setEnabled(True); self.stop_btn.setEnabled(False); self.view_btn.setEnabled(ok)
        self.win.review.refresh_files()
        self.win.set_step_state(RUN_STEP, "ok" if ok else "error")

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
        self.plot_inj = InputsPlot(self, "injections", 2.6)
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
        m.addAction("New experiment (blank)", self._new)
        m.addAction("Open experiment (.experiment.json) …", self._open)
        m.addAction("Save experiment", self._save)
        m.addSeparator()
        m.addAction("Expert tool (FR_GUI) …", self._expert)
        m.addSeparator()
        m.addAction("Quit", self.close)
        self.statusBar().showMessage(f"settings: {os.path.basename(toml_path) if toml_path else 'defaults'}   |   include files: {os.path.basename(self.include_dir)}")

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
        for w in (c.src_archive, c.src_manual):
            w.toggled.connect(self._schedule_refresh)
        for w in (c.step,):
            w.valueChanged.connect(self._schedule_refresh)
        for w in (c.jsrc, c.temp_mode, c.chamber, ch.base_combo, ch.mech_combo):
            w.currentIndexChanged.connect(self._schedule_refresh)
        for t in (c.roof, inj.table, inj.chain, ch.adjust, self.page_init.table):
            t.itemChanged.connect(self._schedule_refresh); t.model().rowsRemoved.connect(self._schedule_refresh); t.model().rowsInserted.connect(self._schedule_refresh)
        ch.chosen.model().rowsInserted.connect(self._schedule_refresh); ch.chosen.model().rowsRemoved.connect(self._schedule_refresh)
        for w in (ch.mode_compounds, ch.mode_existing):
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
        no_src = not (c.src_archive.isChecked() or c.src_manual.isChecked())
        theme.mark_missing(c.src_archive, no_src); theme.mark_missing(c.src_manual, no_src)
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
        txt = " ".join(m for _, m in msgs).lower()
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
        self.load_spec(blank_spec()); self.steps.setCurrentRow(0); self.review.refresh()

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
        self.steps.setCurrentRow(0); self.review.refresh()

    def _save(self):
        try:
            res = self.write_all()
            self.statusBar().showMessage(f"saved {res['spec_path']}"); self.review.refresh_files()
        except Exception as e:
            QMessageBox.critical(self, "Save failed", str(e))

    def _expert(self):
        p = QProcess(self); args = [os.path.join(HERE, "FR_GUI.py"), "-timeout", "0"]
        if self.toml_path:
            args += ["-t", self.toml_path]
        p.start(sys.executable, args); self._expert_proc = p


def main(argv=None):
    argv = list(sys.argv if argv is None else argv)
    theme.taskbar_identity("FZJ.JUICE.GUI")
    app = QApplication(argv)
    theme.apply_theme(app)
    w = JuiceWindow(pick_toml(argv))
    w.show()
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
