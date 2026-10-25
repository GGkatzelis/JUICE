"""
overview_tab.py — the landing page of the results viewer.

Built automatically from the run: the experiment description next to the result file
(.experiment.json) says what was injected and what the chamber adds, the result file gives
the model, and the SAPHIR observations of that day come from `<date>.observations.csv`
next to the result (extracted once from the archive by experiment_plots.load_observations).
Panels:

    Injected compounds · NOx and ozone · Radicals · HONO and CO · OH reactivity ·
    Chamber sources and conditions

Compact and interactive without a toolbar per panel: the time axis is shared by all panels
— scroll wheel zooms the time around the cursor (Ctrl+wheel, or a right-button drag up/down,
zooms the value axis of that panel), left-drag pans, double-click resets, moving the mouse reads the values at that time,
clicking a legend entry hides or shows a line.  One row of buttons serves all panels.
"""
from __future__ import annotations

import os
from typing import Dict, List, Optional

import numpy as np
from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QGridLayout, QLabel, QCheckBox, QPushButton,
                             QScrollArea, QFrame, QLineEdit, QFileDialog, QToolButton, QSizePolicy, QComboBox)
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

import experiment_plots as ep

ACCENT = "#6b3e85"
# pages of the Overview: name → builder(self) returning panel dicts (see experiment_plots.landing_panels for the keys)
PAGES = {
    "1  Experiment overview": "page_overview",
    "2  Model vs measurements (1:1)": "page_one_to_one",
    "3  Model − measurement residuals": "page_residuals",
    "4  Products and reservoirs": "page_products",
    "5  Chamber conditions": "page_conditions",
}
MODEL_COLORS = ["#6b3e85", "#3b6fc4", "#2e9e4f", "#d9822b", "#c0392b", "#8e6c3a", "#1aa3a3", "#7f7f7f"]
PANEL_MIN_H = 270
try:
    from juice_theme import PLOT_FONT as F
except Exception:
    F = dict(base=13, title=14, label=13, tick=11, legend=11, small=10, readout=11)


class PanelGroup:
    """The panels of one page share the time axis: a zoom or pan in one moves all."""

    def __init__(self):
        self.panels: List["InteractivePanel"] = []
        self.home = None            # (lo, hi) of the full time range

    def add(self, p: "InteractivePanel"):
        self.panels.append(p)

    def set_xlim(self, lo: float, hi: float, source=None):
        if not np.isfinite(lo) or not np.isfinite(hi) or hi - lo < 1e-3:
            return
        for p in self.panels:
            p.ax.set_xlim(lo, hi)
            p.canvas.draw_idle()

    def reset(self):
        for p in self.panels:
            p.ax.autoscale(enable=True, axis="both")
            if p.ax2 is not None:
                p.ax2.autoscale(enable=True, axis="y")
            p.ax.relim(); p.ax.autoscale_view()
            if p.ax2 is not None:
                p.ax2.relim(); p.ax2.autoscale_view()
            if p.log:
                p.ax.set_ylim(bottom=1e4)
        if self.home is not None:
            for p in self.panels:
                p.ax.set_xlim(*self.home)
        for p in self.panels:
            p.canvas.draw_idle()


class InteractivePanel(QWidget):
    """One compact matplotlib panel: model lines, observation markers, hover readout, legend toggle,
    wheel zoom, drag pan, double-click reset.  No toolbar — the page has one row for all panels."""

    def __init__(self, title: str, group: PanelGroup, log: bool = False, xlabel: bool = False, parent=None):
        super().__init__(parent)
        self.title = title; self.group = group; self.log = log; self.xlabel = xlabel
        lay = QVBoxLayout(self); lay.setContentsMargins(0, 0, 0, 0); lay.setSpacing(0)
        self.fig = Figure(figsize=(5.0, 2.4)); self.canvas = FigureCanvas(self.fig)
        self.canvas.setMinimumHeight(PANEL_MIN_H)
        self.ax = self.fig.add_subplot(111); self.ax2 = None
        if log:
            self.ax.set_yscale("log")
        lay.addWidget(self.canvas, stretch=1)
        self._lines = {}                 # legend handle → artist
        self._vline = None; self._text = None
        self._drag = None                # (x0 data, xlim at press)
        self.canvas.mpl_connect("motion_notify_event", self._motion)
        self.canvas.mpl_connect("axes_leave_event", lambda _e: self._clear())
        self.canvas.mpl_connect("pick_event", self._pick)
        self.canvas.mpl_connect("scroll_event", self._scroll)
        self.canvas.mpl_connect("button_press_event", self._press)
        self.canvas.mpl_connect("button_release_event", self._release)
        try:                                   # drag the legend; the panel's own pan is cancelled while it moves
            from plot_interaction import attach_legends
            attach_legends(self.canvas, on_grab=lambda: setattr(self, "_drag", None))
        except Exception:
            pass
        group.add(self)

    # ── layout ──
    def finish(self, ylabel: str, ylabel2: str = ""):
        ax = self.ax
        ax.set_title(self.title, loc="left", fontsize=F["title"], color=ACCENT, fontweight="bold", pad=4)   # above the axes: never under a legend
        ax.set_ylabel(ylabel, fontsize=F["label"], labelpad=3)
        if self.xlabel:
            ax.set_xlabel("hour of day (UTC)", fontsize=F["label"], labelpad=2)
        ax.tick_params(labelsize=F["tick"], pad=2, length=3); ax.grid(True, alpha=0.25)
        if self.ax2 is not None:
            self.ax2.set_ylabel(ylabel2, fontsize=F["label"], labelpad=3); self.ax2.tick_params(labelsize=F["tick"], pad=2, length=3)
        handles, labels = ax.get_legend_handles_labels()
        if self.ax2 is not None:
            h2, l2 = self.ax2.get_legend_handles_labels(); handles += h2; labels += l2
        if handles:
            leg = ax.legend(handles, labels, fontsize=F["legend"], loc="best", framealpha=0.85, borderpad=0.3, ncol=2 if len(handles) > 4 else 1,
                            labelspacing=0.2, handlelength=1.5, borderaxespad=0.2)
            for lh, art in zip(leg.get_lines(), handles):
                lh.set_picker(True); lh.set_pickradius(6); self._lines[lh] = art
        self.fig.subplots_adjust(left=0.13, right=0.985 if self.ax2 is None else 0.88, top=0.86, bottom=0.2 if self.xlabel else 0.1)
        self.canvas.draw_idle()

    # ── interaction ──
    def _pick(self, event):
        art = self._lines.get(event.artist)
        if art is None:
            return
        art.set_visible(not art.get_visible()); event.artist.set_alpha(1.0 if art.get_visible() else 0.25)
        self.canvas.draw_idle()

    @staticmethod
    def _modifiers():
        """Keyboard modifiers straight from Qt: matplotlib's event.key is only filled after a key press
        on a focused canvas, so Ctrl+wheel would otherwise never be recognised."""
        from PyQt6.QtWidgets import QApplication
        from PyQt6.QtCore import Qt as _Qt
        m = QApplication.keyboardModifiers()
        return bool(m & _Qt.KeyboardModifier.ControlModifier), bool(m & _Qt.KeyboardModifier.ShiftModifier)

    def _axes(self):
        return [a for a in (self.ax, self.ax2) if a is not None]

    @staticmethod
    def _zoom_axis(ax, y: float, f: float, lim=None):
        """Zoom one value axis by factor f around y in its own units (log-aware)."""
        lo, hi = lim if lim is not None else ax.get_ylim()
        if ax.get_yscale() == "log":
            import math
            if y <= 0 or lo <= 0:
                return
            ly, llo, lhi = math.log10(y), math.log10(lo), math.log10(hi)
            ax.set_ylim(10 ** (ly - (ly - llo) * f), 10 ** (ly + (lhi - ly) * f))
        else:
            ax.set_ylim(y - (y - lo) * f, y + (hi - y) * f)

    def _zoom_y(self, pixel_y: float, f: float):
        """Zoom both value axes of this panel (left and, if present, right) by factor f around the
        cursor height — each axis around that height in its own units, so they stay aligned."""
        for ax in self._axes():
            y = ax.transData.inverted().transform((0.0, pixel_y))[1]
            self._zoom_axis(ax, y, f)
        self.canvas.draw_idle()

    def _scroll(self, event):
        if event.inaxes is None or event.xdata is None:
            return
        f = 0.8 if event.button == "up" else 1.25
        ctrl, shift = self._modifiers()
        if ctrl or shift:                          # value axes of the panel under the cursor only
            self._zoom_y(event.y, f); return
        lo, hi = self.ax.get_xlim(); x = event.xdata
        self.group.set_xlim(x - (x - lo) * f, x + (hi - x) * f, source=self)

    def _press(self, event):
        if event.inaxes is None:
            return
        if event.button == 1:
            if event.dblclick:
                self.group.reset(); return
            self._drag = ("x", event.x, self.ax.get_xlim(), None)
        elif event.button == 3:                     # right button: drag up/down = zoom the values (both axes)
            self._drag = ("y", event.y, [a.get_ylim() for a in self._axes()], None)

    def _release(self, _event):
        self._drag = None

    def _motion(self, event):
        if self._drag is not None and event.x is not None:
            kind, p0, lim, ax = self._drag
            bbox = self.ax.get_window_extent()
            if kind == "x" and bbox.width > 0:
                lo, hi = lim
                dx = (event.x - p0) / bbox.width * (hi - lo)
                self.group.set_xlim(lo - dx, hi - dx, source=self)
            elif kind == "y" and bbox.height > 0:
                import math
                f = math.exp(-(event.y - p0) / bbox.height * 2.0)    # drag up = zoom in, down = zoom out
                for a, (lo, hi) in zip(self._axes(), lim):           # around the middle of each axis' start view
                    if a.get_yscale() == "log" and lo > 0:
                        c = 10 ** (0.5 * (math.log10(lo) + math.log10(hi)))
                    else:
                        c = 0.5 * (lo + hi)
                    self._zoom_axis(a, c, f, (lo, hi))
                self.canvas.draw_idle()
            return
        self._hover(event)

    def _hover(self, event):
        if event.inaxes is None or event.xdata is None:
            return
        x = event.xdata
        rows = []
        for ax in (self.ax, self.ax2):
            if ax is None:
                continue
            for ln in ax.get_lines():
                if not ln.get_visible() or ln.get_label().startswith("_"):
                    continue
                xd, yd = ln.get_xdata(), ln.get_ydata()
                if len(xd) == 0:
                    continue
                i = int(np.argmin(np.abs(np.asarray(xd, dtype=float) - x)))
                if abs(float(xd[i]) - x) > 0.5:
                    continue
                rows.append(f"{ln.get_label()}: {float(yd[i]):.3g}")
        self._clear(draw=False)
        self._vline = self.ax.axvline(x, color="#888", lw=0.8, ls=":")
        h, m = int(x), int(round((x - int(x)) * 60))
        self._text = self.ax.annotate(f"{h:02d}:{m:02d}\n" + "\n".join(rows), xy=(x, 0.02), xycoords=("data", "axes fraction"),
                                      va="bottom", ha="left", fontsize=F["readout"], xytext=(4, 0), textcoords="offset points",
                                      bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#aaa", alpha=0.9))
        self.canvas.draw_idle()

    def _clear(self, draw=True):
        for a in (self._vline, self._text):
            if a is not None:
                try:
                    a.remove()
                except Exception:
                    pass
        self._vline = self._text = None
        if draw:
            self.canvas.draw_idle()


def _small_button(text: str, tip: str, slot) -> QToolButton:
    b = QToolButton(); b.setText(text); b.setToolTip(tip); b.clicked.connect(slot); b.setAutoRaise(True)
    return b


class OverviewTab(QWidget):
    """Landing page: panels from experiment_plots.landing_panels with observations overlaid."""

    def __init__(self, viewer, parent=None):
        super().__init__(parent)
        self.v = viewer
        self.spec = ep.spec_for_run(viewer.csv_path) if getattr(viewer, "csv_path", "") else None
        self.group = PanelGroup()
        lay = QVBoxLayout(self); lay.setContentsMargins(4, 2, 4, 2); lay.setSpacing(2)
        self.head = QLabel(""); self.head.setWordWrap(True)
        self.head.setStyleSheet("color: #6a6a6a; font-size: 9pt; font-style: italic;")     # experiment notes: small, grey, italic
        lay.addWidget(self.head)
        bar = QHBoxLayout(); bar.setSpacing(6)
        bar.addWidget(QLabel("Page"))
        self.page = QComboBox(); self.page.addItems(list(PAGES)); self.page.setToolTip("which set of panels to show; every page can be added to the figure list")
        self.page.currentIndexChanged.connect(self.rebuild); bar.addWidget(self.page)
        self.obs_box = QCheckBox("SAPHIR observations"); self.obs_box.setChecked(True); self.obs_box.toggled.connect(self.rebuild)
        bar.addWidget(self.obs_box)
        self.archive = QLineEdit(self.spec.conditions.archive_root if self.spec else r"Z:\IEK8-SAPHIR\data")
        self.archive.setFixedWidth(150); self.archive.setToolTip("archive root; only read when no observation file exists next to the result")
        bar.addWidget(QLabel("archive")); bar.addWidget(self.archive)
        bar.addWidget(_small_button("Re-extract", "read the archive again and rewrite the observation file next to the result", self._reextract))
        bar.addWidget(_small_button("Refresh", "rebuild the panels", self.rebuild))
        bar.addWidget(_small_button("Reset view", "full time range in all panels (or double-click a panel)", lambda: self.group.reset()))
        bar.addWidget(_small_button("Save image", "save the page as a PNG", self._save_png))
        bar.addStretch()
        hint = QLabel("wheel: zoom time (all panels) · Ctrl+wheel or right-drag up/down: zoom values · left-drag: pan · double-click: reset · click legend: hide")
        hint.setStyleSheet("color: #888; font-size: 10px;"); hint.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Preferred); bar.addWidget(hint, stretch=1)
        lay.addLayout(bar)
        self.status = QLabel(""); self.status.setStyleSheet("color: #666; font-size: 10px;")
        self.status.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Preferred); lay.addWidget(self.status)
        self.area = QScrollArea(); self.area.setWidgetResizable(True); self.area.setFrameShape(QFrame.Shape.NoFrame)
        self.inner = QWidget(); self.grid = QGridLayout(self.inner); self.grid.setSpacing(2); self.grid.setContentsMargins(0, 0, 0, 0)
        self.area.setWidget(self.inner)
        lay.addWidget(self.area, stretch=1)
        self.rebuild()

    # ── data helpers ──
    def _date(self) -> str:
        if self.spec and self.spec.conditions.date:
            return self.spec.conditions.date
        try:
            import datetime as _dt
            return _dt.datetime.fromtimestamp(self.v._midnight_unix, _dt.timezone.utc).strftime("%Y-%m-%d")
        except Exception:
            return ""

    def _run_dir(self) -> str:
        return os.path.dirname(os.path.abspath(self.v.csv_path)) if getattr(self.v, "csv_path", "") else ""

    def _model(self, species: str, kind: str):
        df = self.v.df
        if species not in df.columns:
            return None, None
        y = df[species].values.astype(float)
        if kind == "ppb" and "M" in df.columns:
            y = y / df["M"].values.astype(float) * 1e9
        return df.index.values.astype(float), y

    def _koh_model(self):
        """kOH from the .pd file next to the result (Σ k·[X]) if present."""
        if getattr(self, "_koh_cache", None) is not None:
            return self._koh_cache
        self._koh_cache = (None, None)
        try:
            p = os.path.splitext(self.v.csv_path)[0] + ".pd"
            if os.path.isfile(p):
                t, cols = ep.read_columns(p, ["kOH"])           # only that column, not the whole file
                if "kOH" in cols:
                    self._koh_cache = ((t + ep.ENZ_EPOCH_UNIX - self.v._midnight_unix) / 3600.0, cols["kOH"])
        except Exception:
            pass
        return self._koh_cache

    # ── build ──
    def _reextract(self):
        date = self._date(); d = self._run_dir()
        if date and d:
            f = ep.observation_file(d, date)
            if os.path.isfile(f):
                try:
                    os.remove(f)
                except OSError:
                    pass
        self.rebuild()

    def rebuild(self):
        while self.grid.count():
            w = self.grid.takeAt(0).widget()
            if w is not None:
                w.deleteLater()
        self.group = PanelGroup()
        date = self._date()
        panels = getattr(self, PAGES[self.page.currentText()])()
        obs: Dict = {}; note: List[str] = []
        if (self.obs_box.isChecked() or any(p.get("needs_obs") for p in panels)) and date:
            keys = sorted({k for p in panels for k in p["obs"]})
            obs = ep.load_observations(self.archive.text().strip(), date, keys, self.v._midnight_unix,
                                       cache_dir=self._run_dir(), status=note)
        if any(p.get("needs_obs") for p in panels):
            panels = [p for p in panels if not p.get("needs_obs") or any(k in obs for k in p["obs"])] or \
                     [dict(title="no measurements for this day", species=[], kind="ppb", obs=[], log=False, empty=True)]
        self.head.setText(self._headline(date))
        self.status.setText(("observations " + (note[0] if note else "") + ": " + (", ".join(sorted(obs)) or "none of the mapped series present"))
                            if self.obs_box.isChecked() else "model only")
        n = len(panels); grid_cols = next((p["cols"] for p in panels if p.get("cols")), 3 if n >= 5 else 2)
        rows = (n + grid_cols - 1) // grid_cols
        self.last_grid = (rows, grid_cols)               # QGridLayout never forgets old columns; this is the real arrangement
        pairs = any(p.get("scatter11") for p in panels) and grid_cols == 2
        for c in range(max(3, grid_cols)):               # time series wide, the square correlation narrow
            self.grid.setColumnStretch(c, (5 if c == 0 else 2) if pairs else 1)
        spans = self._hv_spans()
        for i, p in enumerate(panels):
            if p.get("scatter11"):
                self.grid.addWidget(self._scatter_panel(p, obs), i // grid_cols, i % grid_cols); continue
            panel = InteractivePanel(p["title"], self.group, log=p["log"], xlabel=(i // grid_cols == rows - 1))
            for lo, hi in spans:                               # roof open: the same shading in every panel
                panel.ax.axvspan(lo, hi, color="#fffacd", alpha=0.35, zorder=0, linewidth=0)
            self._fill(panel, p, obs)
            self.grid.addWidget(panel, i // grid_cols, i % grid_cols)
        self._set_home()

    def _hv_spans(self):
        """Roof-open periods (hours) from the model's HV column, else from the experiment description."""
        df = self.v.df
        col = next((c for c in df.columns if c.upper() in ("SAPHIR_HV", "HV")), None)
        spans, start = [], None
        if col is not None:
            t = df.index.values.astype(float); hv = df[col].values.astype(float) > 0.5
            for ti, on in zip(t, hv):
                if on and start is None:
                    start = ti
                elif not on and start is not None:
                    spans.append((start, ti)); start = None
            if start is not None:
                spans.append((start, float(t[-1])))
        elif self.spec is not None:
            for o, c in self.spec.conditions.roof_open:
                try:
                    h = lambda x: int(x.split(":")[0]) + int(x.split(":")[1]) / 60.0
                    spans.append((h(o), h(c)))
                except Exception:
                    pass
        return spans

    def _set_home(self):
        """Full time range = the model run; shared by all panels."""
        try:
            t = self.v.df.index.values.astype(float)
            lo, hi = float(np.nanmin(t)), float(np.nanmax(t))
            if hi > lo:
                self.group.home = (lo, hi); self.group.set_xlim(lo, hi)
        except Exception:
            pass

    def _save_png(self):
        start = os.path.join(self._run_dir(), "overview.png") if self._run_dir() else "overview.png"
        path, _ = QFileDialog.getSaveFileName(self, "Save the overview", start, "PNG image (*.png)")
        if path:
            self.inner.grab().save(path)
            self.status.setText(f"saved {path}")

    # ── pages ──
    def page_overview(self):
        return ep.landing_panels(self.spec, list(self.v.df.columns))

    def _measured_species(self):
        cols = set(self.v.df.columns)
        return [k for k in ("NO", "NO2", "O3", "CO", "HONO", "OH", "HO2", "RO2") if k in cols]

    def page_one_to_one(self):
        """One row per measured species: the time series (model line, measurement dots) on the left,
        the model-vs-measurement correlation on the right."""
        out = []
        for k in self._measured_species():
            kind = ep.OBS_MAP[k][5]
            out.append(dict(title=k, species=[k], kind=kind, obs=[k], log=False, needs_obs=True, cols=2))
            out.append(dict(title=f"{k}: model vs {ep.OBS_MAP[k][4]}", species=[k], kind=kind, obs=[k], log=False, scatter11=True, needs_obs=True, cols=2))
        if os.path.isfile(os.path.splitext(self.v.csv_path)[0] + ".pd"):
            out.append(dict(title="OH reactivity", species=[], kind="s-1", obs=["kOH"], log=False, koh=True, needs_obs=True, cols=2))
            out.append(dict(title="kOH: model vs LP-LIF", species=[], kind="s-1", obs=["kOH"], log=False, scatter11=True, koh=True, needs_obs=True, cols=2))
        return out

    def page_residuals(self):
        out = [dict(title=f"{k}: model − measurement", species=[k], kind=ep.OBS_MAP[k][5], obs=[k], log=False, residual=True, needs_obs=True)
               for k in self._measured_species()]
        if os.path.isfile(os.path.splitext(self.v.csv_path)[0] + ".pd"):
            out.append(dict(title="kOH: model − measurement", species=[], kind="s-1", obs=["kOH"], log=False, residual=True, koh=True, needs_obs=True))
        return out

    def page_products(self):
        df = self.v.df; cols = set(df.columns)
        inj = set(ep.injected_species(self.spec)) if self.spec else set()
        inorg = {"NO", "NO2", "O3", "CO", "CO2", "HONO", "HNO3", "H2O2", "OH", "HO2", "RO2", "H2", "H2O", "SO2", "N2O5", "NO3", "O", "O1D", "CH4", "M"}
        try:
            from mcm_species_db import MCM_DB
            def closed_shell(n):           # no peroxy (…O2) or alkoxy (…O after a non-H letter/digit) radicals; aldehydes (…CHO) stay
                u = n.upper()
                return not u.endswith("O2") and not (u.endswith("O") and not u.endswith("HO"))
            org = [c for c in cols if c.upper() not in inorg and c.upper() not in inj and c in MCM_DB and MCM_DB[c].n_C > 0
                   and closed_shell(c) and not c.upper().startswith("SAPHIR_")]
        except Exception:
            org = [c for c in cols if c.upper() not in inorg and c.upper() not in inj and c.isupper()]
        M = df["M"].values.astype(float) if "M" in cols else None
        def peak(c):
            v = df[c].values.astype(float)
            return float(np.nanmax(v / M * 1e9)) if M is not None else float(np.nanmax(v))
        top = sorted(org, key=peak, reverse=True)[:8]
        have = lambda names: [s for s in names if s in cols]
        return [dict(title="Main products (largest peak)", species=top, kind="ppb", obs=[], log=False),
                dict(title="Carbonyls", species=have(["HCHO", "CH3CHO", "CH3COCH3", "MEK", "GLYOX", "MGLYOX", "MACR", "MVK"]), kind="ppb", obs=[], log=False),
                dict(title="Peroxides and acids", species=have(["H2O2", "CH3OOH", "HCOOH", "CH3CO2H", "CH3CO3H"]), kind="ppb", obs=[], log=False),
                dict(title="Nitrogen reservoirs", species=have(["PAN", "HNO3", "HONO", "N2O5", "HNO4"]), kind="ppb", obs=["HONO"], log=False)]

    def page_conditions(self):
        cols = set(self.v.df.columns); have = lambda names: [s for s in names if s in cols]
        return [dict(title="Temperature and humidity", species=have(["T"]) + have(["RH"]), kind="native", obs=[], log=False, twin="RH", units=("K", "RH [%]")),
                dict(title="Photolysis: jNO2 and jO1D", species=have(["JNO2", "J4"])[:1] + have(["JO1D", "J1"])[:1], kind="native", obs=["jNO2"], log=False,
                     twin=have(["JO1D", "J1"])[:1][0] if have(["JO1D", "J1"]) else None, units=("jNO2 [1/s]", "jO1D [1/s]")),
                dict(title="Chamber sources and conditions", species=[], kind="sources", obs=["jNO2"], log=False, sources=True),
                dict(title="NOx and ozone", species=have(["NO", "NO2", "O3"]), kind="ppb", obs=["NO", "NO2", "O3"], log=False,
                     twin="O3", units=("NO, NO2 [ppb]", "O3 [ppb]"))]

    def _model_at(self, species: str, kind: str, hours: np.ndarray):
        t, y = self._model(species, kind)
        if t is None:
            return None
        ok = np.isfinite(y)
        return np.interp(hours, t[ok], y[ok], left=np.nan, right=np.nan)

    def _scatter_panel(self, p: dict, obs: Dict):
        from timeseries_tab import CorrelationPanel
        w = CorrelationPanel()
        key = p["obs"][0]
        if key not in obs:
            w.ax.text(0.5, 0.5, "no measurement", ha="center", va="center", color="#888", transform=w.ax.transAxes); w.canvas.draw_idle(); return w
        h, v, lab, kind = obs[key]
        if p.get("koh"):
            t, k = self._koh_model()
            if t is None:
                w.ax.text(0.5, 0.5, "no kOH in the .pd file", ha="center", va="center", color="#888", transform=w.ax.transAxes); w.canvas.draw_idle(); return w
            ok = np.isfinite(k); ym = np.interp(h, t[ok], k[ok], left=np.nan, right=np.nan)
        else:
            ym = self._model_at(key, kind, h)
        unit = {"ppb": "ppb", "cm-3": "cm⁻³", "s-1": "1/s"}.get(kind, "")
        w.plot(f"{lab} [{unit}]", np.asarray(v, float), [(f"{key} model [{unit}]", ym)], np.asarray(h, float), p["title"], one_to_one=True)
        return w

    def _headline(self, date: str) -> str:
        if self.spec is None:
            return f"<b>{os.path.basename(self.v.csv_path or '')}</b> — no experiment description next to the file; panels from the model species only."
        s = self.spec; c = s.conditions
        items = [f"{i.start or '--:--'} {i.name} {i.total_ppb} ppb" for i in s.injections if i.products]
        inj = "; ".join(items[:6]) + (f"; … {len(items) - 6} more" if len(items) > 6 else "")
        chem = ", ".join(s.precursors) if s.precursors else s.mechanism_stem
        return (f"<b>{s.name}</b> — {c.date} {c.start}–{c.end} UTC, {'measured' if c.source == 'archive' else 'manual'} conditions, "
                f"chemistry {s.base_mechanism} ({chem}){' + ' + ', '.join(s.addons) if s.addons else ''}; "
                f"injections: {inj or 'none'}")

    def _fill(self, panel: InteractivePanel, p: dict, obs: Dict):
        ax = panel.ax
        if p.get("empty"):
            ax.text(0.5, 0.5, p["title"], ha="center", va="center", color="#888", transform=ax.transAxes); panel.finish(""); return
        if p.get("sources"):
            self._fill_sources(panel, obs); return
        if p.get("residual"):
            key = p["obs"][0]
            if key in obs:
                h, v, lab, kind = obs[key]
                if p.get("koh"):
                    t, k = self._koh_model(); ym = None if t is None else np.interp(h, t[np.isfinite(k)], k[np.isfinite(k)], left=np.nan, right=np.nan)
                else:
                    ym = self._model_at(key, kind, h)
                if ym is not None:
                    d = ym - np.asarray(v, float)
                    ax.axhline(0, color="#444", lw=0.8, ls=":")
                    ax.plot(h, d, ".", ms=2.5, color=MODEL_COLORS[0], alpha=0.6, label=f"{key} model − {lab}")
                    ok = np.isfinite(d)
                    if ok.sum() > 2:
                        ax.text(0.99, 0.03, f"mean {np.nanmean(d):.3g}, σ {np.nanstd(d):.3g}", transform=ax.transAxes, ha="right", va="bottom", fontsize=7, color="#444")
            unit = {"ppb": "ppb", "cm-3": "molecules cm⁻³", "s-1": "1/s"}.get(p["kind"], "")
            panel.finish(f"Δ [{unit}]"); return
        if p.get("koh"):
            t, k = self._koh_model()
            if t is not None:
                ax.plot(t, k, color=MODEL_COLORS[0], lw=1.3, label="kOH model (Σ k·[X])")
            if "kOH" in obs:
                h, v, lab, _k = obs["kOH"]; ax.plot(h, v, ".", ms=3, color="#444", alpha=0.6, label=lab)
            panel.finish("kOH [1/s]"); return
        twin = p.get("twin")
        for i, sp in enumerate(p["species"]):
            t, y = self._model(sp, p["kind"])
            if t is None:
                continue
            target = ax
            if twin and sp == twin:
                panel.ax2 = panel.ax2 or ax.twinx(); target = panel.ax2
            target.plot(t, y, color=MODEL_COLORS[i % len(MODEL_COLORS)], lw=1.3, label=f"{sp} model")
        for i, key in enumerate(p["obs"]):
            if key not in obs:
                continue
            h, v, lab, kind = obs[key]
            target = panel.ax2 if (twin and key == twin and panel.ax2 is not None) else ax
            if p["log"]:
                v = np.where(v > 0, v, np.nan)
            target.plot(h, v, ".", ms=2.5, color=MODEL_COLORS[i % len(MODEL_COLORS)], alpha=0.5, label=lab)
        unit = {"ppb": "ppb", "cm-3": "molecules cm⁻³", "s-1": "1/s"}.get(p["kind"], "")
        if p["log"]:
            ax.set_ylim(bottom=1e4)                       # hide the floor value of the flushed chamber
        if p.get("units"):
            panel.finish(p["units"][0], p["units"][1] if twin else ""); return
        panel.finish((f"HO2, RO2 [{unit}]" if twin == "OH" else unit), f"{twin} [{unit}]" if twin else "")

    def _fill_sources(self, panel: InteractivePanel, obs: Dict):
        ax = panel.ax
        if self.spec is None:
            ax.text(0.5, 0.5, "no experiment description", ha="center", va="center", transform=ax.transAxes); panel.finish(""); return
        try:
            enz = ep.enz_next_to(self.v.csv_path) if getattr(self.v, "csv_path", "") else ""
            if enz:                                   # the ENZ the run used: no archive access needed
                grid, cols = ep.read_columns(enz, ["sz", "hv", "dilute", "t", "H2O", "jNO2"])
            else:
                from experiment_spec import enz_columns
                grid, cols, _notes = enz_columns(self.spec)
        except Exception as e:
            ax.text(0.5, 0.5, f"cannot rebuild the inputs: {e}", ha="center", va="center", fontsize=8, wrap=True, transform=ax.transAxes)
            panel.finish(""); return
        hours = (grid + ep.ENZ_EPOCH_UNIX - self.v._midnight_unix) / 3600.0
        src = ep.chamber_sources(self.spec, cols, len(hours))
        colors = iter(MODEL_COLORS)
        for name in ("HONO wall source", "HCHO wall source", "dark HONO"):
            if name in src:
                ax.plot(hours, src[name], color=next(colors), lw=1.2, label=name + " [ppb/h]")
        if "dilution" in src:
            panel.ax2 = ax.twinx(); panel.ax2.plot(hours, src["dilution"], color="#2e9e4f", lw=1, ls="--", label="dilution [%/h]")
        if "jNO2" in obs:
            h, v, lab, _k = obs["jNO2"]
            axj = panel.ax2 or ax.twinx(); panel.ax2 = axj
            axj.plot(h, v * 1e3, ".", ms=2, color="#3b6fc4", alpha=0.4, label="jNO2 measured [1e-3/s]")
        panel.finish("source [ppb/h]", "dilution [%/h] / jNO2")
