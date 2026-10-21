"""
juice_gui.py — JUICE, step by step.

One window, six steps down the left, the active step in the middle, a review panel on
the right that always shows what the experiment contains, and the solver run as the
last step.  No pop-ups on the main path.

    1 Experiment        name, folder, load (saved / wiki / pasted notes), template, date, conditions, roof
    2 Chemistry         precursors (MCM subset) or an existing mechanism
    3 Injections        actions, amounts, mixtures
    4 Chamber           wall deposition, sources, background reactivity, dilution
    5 Initial values    concentrations at the start
    6 Run               prepare the input files, run the solver, open the results

The only thing a user edits is the experiment description (<name>.experiment.json).
The input files (ENZ, EASY, scenario JSON) are generated from it in the project
folder with the IDL/FACSIMILE naming (QNO, QO3, QH2O, Q_<name>; q<NAME> switches)
and can be opened from the review panel, but never need to be.  FR_GUI.py remains
the expert tool for file-level editing.

    python juice_gui.py                 # paths from GG_GUI.toml if present, else FR_GUI.toml
    python juice_gui.py -t FR_GUI.toml  # another TOML

The step pages are the ones of experiment_wizard.py (shared through ExperimentHost).
"""
from __future__ import annotations

import json
import os
import re
import subprocess
import sys
import time
import tomllib
from typing import Dict, List, Optional

import numpy as np
from PyQt6.QtCore import Qt, QProcess, QProcessEnvironment
from PyQt6.QtGui import QPixmap, QIcon, QFont, QPalette, QColor
from PyQt6.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QGridLayout, QLabel,
                             QListWidget, QListWidgetItem, QStackedWidget, QPushButton, QGroupBox, QScrollArea,
                             QComboBox, QLineEdit, QCheckBox, QProgressBar, QTextEdit, QPlainTextEdit, QMessageBox,
                             QFileDialog, QSplitter, QTabWidget, QFrame)
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

from experiment_spec import load_spec, validate, mechanism_species, enz_columns, ExperimentSpec, ACTIONS
from experiment_wizard import ExperimentHost, ReviewPage, CHAMBER_ITEMS, _info
from elog_generator import generate_both
from idl_parser import parse_idl_easy
import injection_calc as ic

HERE = os.path.dirname(os.path.abspath(__file__))
LOGO = os.path.join(HERE, "JUICE Logo.png")
STEPS = ["Experiment", "Chemistry", "Injections", "Chamber", "Initial values", "Run"]
ACCENT, ACCENT_DARK, ACCENT_SOFT = "#6b3e85", "#4a3e85", "#efe6f7"

QSS = f"""
QWidget {{ background: #ffffff; color: #222222; }}
QMainWindow, QScrollArea, QSplitter {{ background: #ffffff; }}
QMenuBar, QMenu, QStatusBar {{ background: #ffffff; color: #222222; }}
QMenu::item:selected {{ background: {ACCENT_SOFT}; }}
QLabel {{ background: transparent; color: #222222; }}
QLineEdit, QPlainTextEdit, QTextEdit, QComboBox, QSpinBox, QListWidget, QTableWidget {{
    background: #ffffff; color: #222222; border: 1px solid #c8c8c8; border-radius: 3px; selection-background-color: {ACCENT_SOFT}; selection-color: #222222; }}
QLineEdit:disabled, QComboBox:disabled, QSpinBox:disabled {{ color: #999999; background: #f4f4f4; }}
QComboBox QAbstractItemView {{ background: #ffffff; color: #222222; selection-background-color: {ACCENT_SOFT}; }}
QHeaderView::section {{ background: #f0f0f0; color: #222222; border: 1px solid #d8d8d8; padding: 3px; }}
QTableWidget {{ gridline-color: #e0e0e0; alternate-background-color: #fafafa; }}
QTableWidget::item:selected, QListWidget::item:selected {{ background: {ACCENT_SOFT}; color: #222222; }}
QPushButton {{ background: #f3f3f3; color: #222222; border: 1px solid #c8c8c8; border-radius: 4px; padding: 5px 12px; }}
QPushButton:hover {{ background: #e9e9e9; }}
QPushButton:disabled {{ color: #a0a0a0; background: #f7f7f7; }}
QPushButton#primary {{ background: {ACCENT}; color: #ffffff; font-weight: bold; padding: 8px 18px; border: none; }}
QPushButton#primary:hover {{ background: {ACCENT_DARK}; }}
QPushButton#primary:disabled {{ background: #c4c4c4; color: #ffffff; }}
QCheckBox, QRadioButton {{ background: transparent; color: #222222; }}
QGroupBox {{ font-weight: bold; color: {ACCENT_DARK}; border: 1px solid #dddddd; border-radius: 6px; margin-top: 12px; background: #ffffff; }}
QGroupBox::title {{ subcontrol-origin: margin; left: 10px; padding: 0 4px; background: #ffffff; }}
QGroupBox::indicator {{ width: 13px; height: 13px; }}
QListWidget#steps {{ background: #ffffff; border: none; font-size: 14px; }}
QListWidget#steps::item {{ padding: 10px 8px; border-bottom: 1px solid #eeeeee; color: #222222; }}
QListWidget#steps::item:selected {{ background: {ACCENT}; color: #ffffff; border-radius: 4px; }}
QLabel#title {{ font-size: 22px; font-weight: bold; color: {ACCENT_DARK}; }}
QLabel#subtitle {{ color: #666666; }}
QTabWidget::pane {{ border: 1px solid #dddddd; background: #ffffff; }}
QTabBar::tab {{ background: #f0f0f0; color: #222222; padding: 6px 12px; border: 1px solid #dddddd; border-bottom: none; }}
QTabBar::tab:selected {{ background: #ffffff; color: {ACCENT_DARK}; font-weight: bold; }}
QProgressBar {{ background: #f0f0f0; color: #222222; border: 1px solid #c8c8c8; border-radius: 3px; text-align: center; }}
QProgressBar::chunk {{ background-color: #2ecc40; }}
QScrollBar:vertical {{ background: #f4f4f4; width: 12px; }} QScrollBar::handle:vertical {{ background: #cccccc; border-radius: 4px; min-height: 20px; }}
QScrollBar:horizontal {{ background: #f4f4f4; height: 12px; }} QScrollBar::handle:horizontal {{ background: #cccccc; border-radius: 4px; min-width: 20px; }}
QToolTip {{ background: #ffffff; color: #222222; border: 1px solid #c8c8c8; }}
"""

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
        v = sol.get(key, default)
        try:
            return float(v)
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
    Returns the JSON path.  The JSON stem = experiment name, so the runner writes <name>.cs.
    Solver numerics are the validated ones from the TOML (solver_from_toml); the user
    never chooses them here."""
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


def spec_bullets(spec: ExperimentSpec, include_dir: str = "") -> List[str]:
    """Plain-language list of what goes into the input files (no file syntax)."""
    c = spec.conditions
    ch = ic.chamber(c.chamber_name)
    out = []
    if spec.precursors:
        out.append(f"Chemistry: MCM v3.3.1 subset for {', '.join(spec.precursors)} + the complete inorganic core and photolysis")
    elif spec.mechanism_stem:
        out.append(f"Chemistry: mechanism include set '{spec.mechanism_stem}' + inorganic core and photolysis")
    else:
        out.append("Chemistry: not chosen yet")
    if c.source == "archive":
        out.append(f"Conditions: measured T, P, humidity, flow and J values from the archive for {c.date}"
                   + (", clear-sky J" if c.j_source == "clear_sky" else ""))
    else:
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
        out.append(f"Injection {inj.start or '--:--'} {label}: {what}{extra}" + (f" [from {mode}]" if mode and mode != "ppb" else ""))
    on = [t for k, t, _d in CHAMBER_ITEMS if getattr(spec.chamber, k)]
    out.append("Chamber effects: " + (", ".join(t.split(" (")[0] for t in on) or "none"))
    out.append(f"Initial values: {len(spec.initials)} species set, all others at the floor value")
    if spec.overrides:
        out.append(f"Expert overrides: {len(spec.overrides)} statement(s) appended")
    return out


def _boxed(page: QWidget, title: str = "") -> QWidget:
    """A wizard page shown as a plain section: title as group box, page inside."""
    box = QGroupBox(title or page.title())
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


# ─────────────────────────────────────────────────────────────────────────────
class ReviewPanel(QWidget):
    """Always visible: summary, what goes into the files, validation, a plot of the
    inputs (conditions and injections against time), and the generated files."""

    def __init__(self, win: "JuiceWindow"):
        super().__init__()
        self.win = win
        lay = QVBoxLayout(self); lay.setContentsMargins(4, 4, 4, 4)
        self.tabs = QTabWidget(); lay.addWidget(self.tabs, stretch=1)
        # summary tab
        st = QWidget(); sl = QVBoxLayout(st)
        self.summary = QLabel(""); self.summary.setWordWrap(True); self.summary.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        sl.addWidget(self.summary)
        self.checks = QListWidget(); self.checks.setMaximumHeight(120); self.checks.setWordWrap(True); sl.addWidget(self.checks)
        sl.addWidget(QLabel("<b>What goes into the input files</b>"))
        self.bullets = QListWidget(); self.bullets.setWordWrap(True); sl.addWidget(self.bullets, stretch=1)
        self.tabs.addTab(st, "Summary")
        # plot tab
        pt = QWidget(); pl = QVBoxLayout(pt)
        self.fig = Figure(figsize=(4, 6), tight_layout=True); self.canvas = FigureCanvas(self.fig); pl.addWidget(self.canvas, stretch=1)
        hb = QHBoxLayout(); b = QPushButton("Update plot"); b.clicked.connect(self.plot); hb.addWidget(b)
        self.plot_note = _info("Conditions and injections as the model will see them (archive data is read when you update)."); hb.addWidget(self.plot_note, stretch=1)
        pl.addLayout(hb)
        self.tabs.addTab(pt, "Inputs plot")
        # files tab
        ft = QWidget(); fl = QVBoxLayout(ft)
        fl.addWidget(_info("Generated in the project folder when you press 'Prepare input files' on the Run step. "
                           "The .easy file is the same input FACSIMILE uses; you do not need to read it."))
        self.files = QListWidget(); self.files.setWordWrap(True); self.files.itemDoubleClicked.connect(lambda it: _open_path(it.data(Qt.ItemDataRole.UserRole))); fl.addWidget(self.files, stretch=1)
        hb2 = QHBoxLayout()
        fb = QPushButton("Open project folder"); fb.clicked.connect(lambda: _open_path(self.win.out_dir()) if os.path.isdir(self.win.out_dir()) else None); hb2.addWidget(fb)
        ob = QPushButton("Open selected file"); ob.clicked.connect(self._open_selected); hb2.addWidget(ob); hb2.addStretch(); fl.addLayout(hb2)
        self.tabs.addTab(ft, "Files")
        self.tabs.currentChanged.connect(lambda i: self.plot() if i == 1 and not self._plotted else None)
        self._plotted = False

    def refresh(self):
        try:
            spec = self.win.collect_spec()
        except Exception as e:
            self.summary.setText(f"<i>cannot read the pages: {e}</i>"); return
        self.checks.clear(); self.bullets.clear()
        if is_blank(spec):
            self.summary.setText("<i>Nothing entered yet.</i><br>Give the experiment a name, then open a saved experiment, "
                                 "import the SAPHIR wiki page, paste your notes, or pick a template.")
            self._plotted = False; self.refresh_files()
            return []
        known = mechanism_species(self.win.include_dir, spec.mechanism_stem) if spec.mechanism_stem else None
        msgs = validate(spec, known)
        self.summary.setText("<b>" + (spec.name or "(no name)") + "</b><br>" + ReviewPage._summary(spec).replace("\n", "<br>"))
        for lvl, m in msgs:
            self.checks.addItem(QListWidgetItem(("✖ " if lvl == "error" else "⚠ ") + m))
        if not msgs:
            self.checks.addItem(QListWidgetItem("✔ No problems found"))
        for b in spec_bullets(spec, self.win.include_dir):
            self.bullets.addItem(QListWidgetItem("• " + b))
        self._plotted = False
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

    def plot(self):
        self.fig.clear()
        try:
            spec = self.win.collect_spec()
            grid, cols, notes = enz_columns(spec)
        except Exception as e:
            ax = self.fig.add_subplot(111); ax.text(0.5, 0.5, f"cannot build the inputs:\n{e}", ha="center", va="center", wrap=True); self.canvas.draw(); return
        hours = (grid - grid[0]) / 3600.0 + float(spec.conditions.start.split(":")[0]) + float(spec.conditions.start.split(":")[1]) / 60.0
        ax1 = self.fig.add_subplot(311); ax2 = self.fig.add_subplot(312, sharex=ax1); ax3 = self.fig.add_subplot(313, sharex=ax1)
        if "t" in cols:
            ax1.plot(hours, cols["t"], color=ACCENT, lw=1.2, label="T [K]"); ax1.set_ylabel("T [K]", color=ACCENT)
        if "hv" in cols:
            ax1b = ax1.twinx(); ax1b.fill_between(hours, 0, cols["hv"], color="#ffd966", alpha=0.4, step="post", label="roof open"); ax1b.set_ylim(0, 1.05); ax1b.set_yticks([])
        ax1.set_title("Conditions", fontsize=9)
        jn = [k for k in cols if k.lower() == "jno2"] or [k for k in cols if k.lower().startswith("j")][:1]
        for k in jn:
            ax2.plot(hours, cols[k], lw=1, color="#3b6fc4", label=k)
        if "H2O" in cols:
            ax2b = ax2.twinx(); ax2b.plot(hours, cols["H2O"], color="#3b82c4", lw=1, ls="--", label="H2O [%]"); ax2b.set_ylabel("H2O [%]", color="#3b82c4")
        ax2.set_ylabel("J [1/s]"); ax2.legend(fontsize=7, loc="upper left") if jn else None
        qs = [k for k in cols if k.upper().startswith("Q")]
        for i, k in enumerate(qs):
            y = cols[k]
            if np.any(y > 0):
                ax3.fill_between(hours, i, i + 0.8 * (y > 0), step="post", color=ACCENT, alpha=0.5)
        ax3.set_yticks([i + 0.4 for i in range(len(qs))]); ax3.set_yticklabels(qs, fontsize=7); ax3.set_ylim(-0.2, max(1, len(qs)))
        ax3.set_xlabel("hour of day (UTC)"); ax3.set_title("Injections", fontsize=9)
        for ax in (ax1, ax2):
            ax.tick_params(labelsize=7)
        ax3.tick_params(labelsize=7)
        self.plot_note.setText("; ".join(notes[:3]) if notes else "Conditions and injections as the model will see them.")
        self.canvas.draw(); self._plotted = True


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
        self.bar = QProgressBar(); self.bar.setRange(0, 100); self.bar.setTextVisible(True); self.bar.setFormat("not started")
        lay.addWidget(self.bar)
        self.log = QTextEdit(); self.log.setReadOnly(True)
        self.log.setStyleSheet("font-family: Consolas, monospace; font-size: 11px; background: #fbfbfb; color: #222222;")
        lay.addWidget(self.log, stretch=1)

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
        self.resize(1550, 920)
        self.setStyleSheet(QSS)
        if os.path.isfile(LOGO):
            self.setWindowIcon(QIcon(QPixmap(LOGO)))

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

        # ── middle: pages ──
        self.stack = QStackedWidget()
        self.stack.addWidget(_scroll(_boxed(self.page_start, "Experiment"), _boxed(self.page_cond, "Where and when")))
        self.stack.addWidget(_scroll(_boxed(self.page_chem, "Chemistry")))
        self.stack.addWidget(_scroll(_boxed(self.page_inj, "Injections")))
        self.stack.addWidget(_scroll(_boxed(self.page_chamber, "What the chamber adds on its own")))
        self.stack.addWidget(_scroll(_boxed(self.page_init, "Initial values")))
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
        m.addAction("New experiment (SAPHIR defaults)", self._new)
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

    # ── navigation ──
    def _go(self, row: int):
        if row < 0:
            return
        self.stack.setCurrentIndex(row)
        self.back_btn.setEnabled(row > 0); self.next_btn.setEnabled(row < len(STEPS) - 1)
        if row == len(STEPS) - 1:
            self.page_run.refresh()
        msgs = self.review.refresh()
        self._update_states(msgs or [])

    def _step(self, d: int):
        row = self.steps.currentRow()
        if d > 0:
            if row == 0 and not (self.page_start.validatePage() and True):
                return
            if row == 1 and not self.page_chem.validatePage():
                return
        self.steps.setCurrentRow(max(0, min(len(STEPS) - 1, row + d)))

    def set_step_state(self, i: int, state: str):
        self._state[i] = state
        self._paint_states()

    def _update_states(self, msgs):
        if not msgs and is_blank(self.collect_spec()):
            self._state[:5] = ["todo"] * 5; self._paint_states(); return
        txt = " ".join(m for _, m in msgs).lower()
        self._state[0] = "error" if any(k in txt for k in ("date", "time window", "conditions come from")) else ("warn" if "j values" in txt or "roof" in txt else "ok")
        self._state[1] = "error" if "no chemistry" in txt else "ok"
        self._state[2] = "warn" if "injection" in txt else "ok"
        self._state[3] = "ok"; self._state[4] = "ok"
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
    app = QApplication(argv)
    app.setStyle("Fusion")                     # same look on every machine; ignores the OS dark mode
    pal = QPalette()
    for role, col in ((QPalette.ColorRole.Window, "#ffffff"), (QPalette.ColorRole.Base, "#ffffff"), (QPalette.ColorRole.AlternateBase, "#fafafa"),
                      (QPalette.ColorRole.Button, "#f3f3f3"), (QPalette.ColorRole.ToolTipBase, "#ffffff"), (QPalette.ColorRole.Highlight, ACCENT_SOFT)):
        pal.setColor(role, QColor(col))
    for role in (QPalette.ColorRole.WindowText, QPalette.ColorRole.Text, QPalette.ColorRole.ButtonText, QPalette.ColorRole.ToolTipText,
                 QPalette.ColorRole.HighlightedText, QPalette.ColorRole.PlaceholderText):
        pal.setColor(role, QColor("#222222" if role != QPalette.ColorRole.PlaceholderText else "#888888"))
    app.setPalette(pal)
    app.setFont(QFont("Segoe UI", 10))
    w = JuiceWindow(pick_toml(argv))
    w.show()
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
