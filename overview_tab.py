"""
overview_tab.py — the landing page of the results viewer.

Built automatically from the run: the experiment description next to the result file
(.experiment.json) says what was injected and what the chamber adds, the result file gives
the model, and the SAPHIR archive gives the observations for that day.  Panels:

    Injected compounds · NOx and ozone · Radicals · HONO and CO · OH reactivity ·
    Chamber sources and conditions

Every panel is interactive: move the mouse to read the values at that time (model and
observations), click a legend entry to hide or show a line, use the toolbar to zoom.
"""
from __future__ import annotations

import os
from typing import Dict, List, Optional

import numpy as np
from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QGridLayout, QLabel, QCheckBox, QPushButton,
                             QScrollArea, QFrame, QLineEdit)
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure

import experiment_plots as ep

ACCENT = "#6b3e85"
MODEL_COLORS = ["#6b3e85", "#3b6fc4", "#2e9e4f", "#d9822b", "#c0392b", "#8e6c3a", "#1aa3a3", "#7f7f7f"]


class InteractivePanel(QWidget):
    """One matplotlib panel: model lines, observation markers, hover readout, legend toggle."""

    def __init__(self, title: str, log: bool = False, parent=None):
        super().__init__(parent)
        self.title = title
        lay = QVBoxLayout(self); lay.setContentsMargins(2, 2, 2, 2)
        self.fig = Figure(figsize=(5.2, 3.0), tight_layout=True); self.canvas = FigureCanvas(self.fig)
        self.canvas.setMinimumHeight(240)
        self.ax = self.fig.add_subplot(111); self.ax2 = None
        if log:
            self.ax.set_yscale("log")
        self.toolbar = NavigationToolbar(self.canvas, self); self.toolbar.setIconSize(self.toolbar.iconSize() * 0.7)
        lay.addWidget(self.toolbar); lay.addWidget(self.canvas, stretch=1)
        self._lines = {}                 # legend handle → (artist, axis)
        self._vline = None; self._text = None
        self.canvas.mpl_connect("motion_notify_event", self._hover)
        self.canvas.mpl_connect("axes_leave_event", lambda _e: self._clear())
        self.canvas.mpl_connect("pick_event", self._pick)

    def finish(self, ylabel: str, ylabel2: str = ""):
        self.ax.set_title(self.title, fontsize=9, loc="left", color=ACCENT, fontweight="bold")
        self.ax.set_ylabel(ylabel, fontsize=8); self.ax.set_xlabel("hour of day (UTC)", fontsize=8)
        self.ax.tick_params(labelsize=7); self.ax.grid(True, alpha=0.25)
        if self.ax2 is not None:
            self.ax2.set_ylabel(ylabel2, fontsize=8); self.ax2.tick_params(labelsize=7)
        handles, labels = self.ax.get_legend_handles_labels()
        if self.ax2 is not None:
            h2, l2 = self.ax2.get_legend_handles_labels(); handles += h2; labels += l2
        if handles:
            leg = self.ax.legend(handles, labels, fontsize=7, loc="best", framealpha=0.85)
            for lh, art in zip(leg.get_lines(), handles):
                lh.set_picker(True); lh.set_pickradius(6); self._lines[lh] = art
        self.canvas.draw_idle()

    def _pick(self, event):
        art = self._lines.get(event.artist)
        if art is None:
            return
        art.set_visible(not art.get_visible()); event.artist.set_alpha(1.0 if art.get_visible() else 0.25)
        self.canvas.draw_idle()

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
        self._text = self.ax.annotate(f"{h:02d}:{m:02d}\n" + "\n".join(rows), xy=(x, 0.98), xycoords=("data", "axes fraction"),
                                      va="top", ha="left", fontsize=7, bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#aaa", alpha=0.9))
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


class OverviewTab(QWidget):
    """Landing page: panels from experiment_plots.landing_panels with observations overlaid."""

    def __init__(self, viewer, parent=None):
        super().__init__(parent)
        self.v = viewer
        self.spec = ep.spec_for_run(viewer.csv_path) if getattr(viewer, "csv_path", "") else None
        lay = QVBoxLayout(self); lay.setContentsMargins(4, 4, 4, 4)
        top = QHBoxLayout()
        self.head = QLabel(""); self.head.setWordWrap(True); top.addWidget(self.head, stretch=1)
        self.obs_box = QCheckBox("Overlay SAPHIR observations"); self.obs_box.setChecked(True); self.obs_box.toggled.connect(self.rebuild)
        top.addWidget(self.obs_box)
        top.addWidget(QLabel("archive")); self.archive = QLineEdit(self.spec.conditions.archive_root if self.spec else r"Z:\IEK8-SAPHIR\data")
        self.archive.setFixedWidth(180); top.addWidget(self.archive)
        rb = QPushButton("Refresh"); rb.clicked.connect(self.rebuild); top.addWidget(rb)
        lay.addLayout(top)
        self.status = QLabel(""); self.status.setStyleSheet("color: #666;"); lay.addWidget(self.status)
        self.area = QScrollArea(); self.area.setWidgetResizable(True); self.area.setFrameShape(QFrame.Shape.NoFrame)
        self.inner = QWidget(); self.grid = QGridLayout(self.inner); self.area.setWidget(self.inner)
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
        try:
            stem = os.path.splitext(self.v.csv_path)[0]
            p = stem + ".pd"
            if not os.path.isfile(p):
                return None, None
            from column_file_reader import load_column_file
            pdf = load_column_file(p)
            if "kOH" not in pdf.columns:
                return None, None
            t = (pdf.index.values.astype(float) + ep.ENZ_EPOCH_UNIX - self.v._midnight_unix) / 3600.0
            return t, pdf["kOH"].values.astype(float)
        except Exception:
            return None, None

    # ── build ──
    def rebuild(self):
        while self.grid.count():
            w = self.grid.takeAt(0).widget()
            if w is not None:
                w.deleteLater()
        date = self._date()
        panels = ep.landing_panels(self.spec, list(self.v.df.columns))
        obs: Dict = {}
        day = ep.archive_day(self.archive.text().strip(), date) if date else ""
        if self.obs_box.isChecked() and day:
            keys = sorted({k for p in panels for k in p["obs"]})
            obs = ep.load_observations(self.archive.text().strip(), date, keys, self.v._midnight_unix)
        self.head.setText(self._headline(date))
        self.status.setText(("observations from " + day + ": " + (", ".join(sorted(obs)) or "none of the mapped files present")) if day
                            else f"no archive folder for {date or 'this day'} — model only")
        grid_cols = 2
        for i, p in enumerate(panels):
            panel = InteractivePanel(p["title"], log=p["log"])
            self._fill(panel, p, obs)
            self.grid.addWidget(panel, i // grid_cols, i % grid_cols)

    def _headline(self, date: str) -> str:
        if self.spec is None:
            return f"<b>{os.path.basename(self.v.csv_path or '')}</b> — no experiment description next to the file; panels from the model species only."
        s = self.spec; c = s.conditions
        inj = "; ".join(f"{i.start or '--:--'} {i.name} {i.total_ppb} ppb" for i in s.injections if i.products)
        chem = ", ".join(s.precursors) if s.precursors else s.mechanism_stem
        return (f"<b>{s.name}</b> — {c.date} {c.start}–{c.end} UTC, {'measured' if c.source == 'archive' else 'manual'} conditions, "
                f"chemistry {s.base_mechanism} ({chem}){' + ' + ', '.join(s.addons) if s.addons else ''}<br>"
                f"injections: {inj or 'none'}")

    def _fill(self, panel: InteractivePanel, p: dict, obs: Dict):
        ax = panel.ax
        if p.get("sources"):
            self._fill_sources(panel, obs); return
        if p.get("koh"):
            t, k = self._koh_model()
            if t is not None:
                ax.plot(t, k, color=MODEL_COLORS[0], lw=1.4, label="kOH model (Σ k·[X])")
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
            target.plot(t, y, color=MODEL_COLORS[i % len(MODEL_COLORS)], lw=1.4, label=f"{sp} model")
        for i, key in enumerate(p["obs"]):
            if key not in obs:
                continue
            h, v, lab, kind = obs[key]
            target = panel.ax2 if (twin and key == twin and panel.ax2 is not None) else ax
            if p["log"]:
                v = np.where(v > 0, v, np.nan)
            target.plot(h, v, ".", ms=3, color=MODEL_COLORS[i % len(MODEL_COLORS)], alpha=0.5, label=lab)
        unit = {"ppb": "ppb", "cm-3": "molecules cm⁻³", "s-1": "1/s"}.get(p["kind"], "")
        if p["log"]:
            ax.set_ylim(bottom=1e4)                       # hide the floor value of the flushed chamber
        panel.finish(unit, f"{twin} [ppb]" if twin else "")

    def _fill_sources(self, panel: InteractivePanel, obs: Dict):
        ax = panel.ax
        if self.spec is None:
            ax.text(0.5, 0.5, "no experiment description", ha="center", va="center", transform=ax.transAxes); panel.finish(""); return
        try:
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
                ax.plot(hours, src[name], color=next(colors), lw=1.3, label=name + " [ppb/h]")
        if "hv" in cols:
            ax.fill_between(hours, 0, np.asarray(cols["hv"], float) * (ax.get_ylim()[1] or 1), color="#ffd966", alpha=0.25, step="post", label="roof open")
        if "dilution" in src:
            panel.ax2 = ax.twinx(); panel.ax2.plot(hours, src["dilution"], color="#2e9e4f", lw=1, ls="--", label="dilution [%/h]")
        if "jNO2" in obs:
            h, v, lab, _k = obs["jNO2"]
            axj = panel.ax2 or ax.twinx(); panel.ax2 = axj
            axj.plot(h, v * 1e3, ".", ms=2, color="#3b6fc4", alpha=0.4, label="jNO2 measured [1e-3/s]")
        panel.finish("source [ppb/h]", "dilution [%/h] / jNO2")
