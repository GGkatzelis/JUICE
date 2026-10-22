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
— scroll wheel zooms the time around the cursor (Ctrl+wheel zooms the value axis of that
panel), left-drag pans, double-click resets, moving the mouse reads the values at that time,
clicking a legend entry hides or shows a line.  One row of buttons serves all panels.
"""
from __future__ import annotations

import os
from typing import Dict, List, Optional

import numpy as np
from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QGridLayout, QLabel, QCheckBox, QPushButton,
                             QScrollArea, QFrame, QLineEdit, QFileDialog, QToolButton)
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

import experiment_plots as ep

ACCENT = "#6b3e85"
MODEL_COLORS = ["#6b3e85", "#3b6fc4", "#2e9e4f", "#d9822b", "#c0392b", "#8e6c3a", "#1aa3a3", "#7f7f7f"]
PANEL_MIN_H = 190


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
        group.add(self)

    # ── layout ──
    def finish(self, ylabel: str, ylabel2: str = ""):
        ax = self.ax
        ax.text(0.01, 0.97, self.title, transform=ax.transAxes, fontsize=8.5, color=ACCENT, fontweight="bold", va="top", ha="left")
        ax.set_ylabel(ylabel, fontsize=7.5, labelpad=2)
        if self.xlabel:
            ax.set_xlabel("hour of day (UTC)", fontsize=7.5, labelpad=1)
        ax.tick_params(labelsize=7, pad=2, length=2); ax.grid(True, alpha=0.25)
        if self.ax2 is not None:
            self.ax2.set_ylabel(ylabel2, fontsize=7.5, labelpad=2); self.ax2.tick_params(labelsize=7, pad=2, length=2)
        handles, labels = ax.get_legend_handles_labels()
        if self.ax2 is not None:
            h2, l2 = self.ax2.get_legend_handles_labels(); handles += h2; labels += l2
        if handles:
            leg = ax.legend(handles, labels, fontsize=6.5, loc="best", framealpha=0.8, borderpad=0.3, ncol=2 if len(handles) > 4 else 1,
                            labelspacing=0.2, handlelength=1.5, borderaxespad=0.2)
            for lh, art in zip(leg.get_lines(), handles):
                lh.set_picker(True); lh.set_pickradius(6); self._lines[lh] = art
        self.fig.subplots_adjust(left=0.11, right=0.99 if self.ax2 is None else 0.9, top=0.97, bottom=0.17 if self.xlabel else 0.1)
        self.canvas.draw_idle()

    # ── interaction ──
    def _pick(self, event):
        art = self._lines.get(event.artist)
        if art is None:
            return
        art.set_visible(not art.get_visible()); event.artist.set_alpha(1.0 if art.get_visible() else 0.25)
        self.canvas.draw_idle()

    def _scroll(self, event):
        if event.inaxes is None or event.xdata is None:
            return
        f = 0.8 if event.button == "up" else 1.25
        ctrl = bool(event.key and "control" in event.key)
        if ctrl:                                   # value axis of this panel only
            ax = event.inaxes
            lo, hi = ax.get_ylim(); y = event.ydata
            if ax.get_yscale() == "log" and y > 0 and lo > 0:
                import math
                ly, llo, lhi = math.log10(y), math.log10(lo), math.log10(hi)
                ax.set_ylim(10 ** (ly - (ly - llo) * f), 10 ** (ly + (lhi - ly) * f))
            else:
                ax.set_ylim(y - (y - lo) * f, y + (hi - y) * f)
            self.canvas.draw_idle(); return
        lo, hi = self.ax.get_xlim(); x = event.xdata
        self.group.set_xlim(x - (x - lo) * f, x + (hi - x) * f, source=self)

    def _press(self, event):
        if event.inaxes is None or event.button != 1:
            return
        if event.dblclick:
            self.group.reset(); return
        self._drag = (event.x, self.ax.get_xlim())

    def _release(self, _event):
        self._drag = None

    def _motion(self, event):
        if self._drag is not None and event.x is not None:
            x0, (lo, hi) = self._drag
            bbox = self.ax.get_window_extent()
            if bbox.width > 0:
                dx = (event.x - x0) / bbox.width * (hi - lo)
                self.group.set_xlim(lo - dx, hi - dx, source=self)
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
                                      va="bottom", ha="left", fontsize=6.5, xytext=(4, 0), textcoords="offset points",
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
        self.head = QLabel(""); self.head.setWordWrap(True); lay.addWidget(self.head)
        bar = QHBoxLayout(); bar.setSpacing(6)
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
        hint = QLabel("wheel: zoom time (all panels) · Ctrl+wheel: zoom values · drag: pan · double-click: reset · hover: read values · click legend: hide")
        hint.setStyleSheet("color: #888; font-size: 10px;"); bar.addWidget(hint)
        lay.addLayout(bar)
        self.status = QLabel(""); self.status.setStyleSheet("color: #666; font-size: 10px;"); lay.addWidget(self.status)
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
        panels = ep.landing_panels(self.spec, list(self.v.df.columns))
        obs: Dict = {}; note: List[str] = []
        if self.obs_box.isChecked() and date:
            keys = sorted({k for p in panels for k in p["obs"]})
            obs = ep.load_observations(self.archive.text().strip(), date, keys, self.v._midnight_unix,
                                       cache_dir=self._run_dir(), status=note)
        self.head.setText(self._headline(date))
        self.status.setText(("observations " + (note[0] if note else "") + ": " + (", ".join(sorted(obs)) or "none of the mapped series present"))
                            if self.obs_box.isChecked() else "model only")
        n = len(panels); grid_cols = 3 if n >= 5 else 2
        rows = (n + grid_cols - 1) // grid_cols
        for i, p in enumerate(panels):
            panel = InteractivePanel(p["title"], self.group, log=p["log"], xlabel=(i // grid_cols == rows - 1))
            self._fill(panel, p, obs)
            self.grid.addWidget(panel, i // grid_cols, i % grid_cols)
        self._set_home()

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

    def _headline(self, date: str) -> str:
        if self.spec is None:
            return f"<b>{os.path.basename(self.v.csv_path or '')}</b> — no experiment description next to the file; panels from the model species only."
        s = self.spec; c = s.conditions
        inj = "; ".join(f"{i.start or '--:--'} {i.name} {i.total_ppb} ppb" for i in s.injections if i.products)
        chem = ", ".join(s.precursors) if s.precursors else s.mechanism_stem
        return (f"<b>{s.name}</b> — {c.date} {c.start}–{c.end} UTC, {'measured' if c.source == 'archive' else 'manual'} conditions, "
                f"chemistry {s.base_mechanism} ({chem}){' + ' + ', '.join(s.addons) if s.addons else ''}; "
                f"injections: {inj or 'none'}")

    def _fill(self, panel: InteractivePanel, p: dict, obs: Dict):
        ax = panel.ax
        if p.get("sources"):
            self._fill_sources(panel, obs); return
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
        panel.finish(unit, f"{twin} [ppb]" if twin else "")

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
        if "hv" in cols:
            ax.fill_between(hours, 0, np.asarray(cols["hv"], float) * (ax.get_ylim()[1] or 1), color="#ffd966", alpha=0.25, step="post", label="roof open")
        if "dilution" in src:
            panel.ax2 = ax.twinx(); panel.ax2.plot(hours, src["dilution"], color="#2e9e4f", lw=1, ls="--", label="dilution [%/h]")
        if "jNO2" in obs:
            h, v, lab, _k = obs["jNO2"]
            axj = panel.ax2 or ax.twinx(); panel.ax2 = axj
            axj.plot(h, v * 1e3, ".", ms=2, color="#3b6fc4", alpha=0.4, label="jNO2 measured [1e-3/s]")
        panel.finish("source [ppb/h]", "dilution [%/h] / jNO2")
