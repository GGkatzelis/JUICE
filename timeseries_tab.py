"""
timeseries_tab.py — the Timeseries tab of the results viewer as a grid of linked panels.

One to six panels (1, 2 or 3 stacked, 2×2, 3×2).  Each panel is either

  • a time series: any model species (ppb or native units, linear or log) with the matching
    SAPHIR observation of that day drawn as dots in the same colour (from the day's
    observation file next to the result, see experiment_plots.load_observations), or
  • a correlation: X species against one or more Y species, coloured by hour of day when
    there is one Y, with a least-squares line (slope, intercept, r²) per pair.

Species come from the tree on the left (select, then "Add to panel", or double-click).
Interaction is the same as on the Overview: wheel zooms the time in all panels at once,
Ctrl+wheel or a right-button drag zooms the values, left-drag pans, double-click resets,
moving the mouse reads the values, clicking a legend entry hides a line.
"""
from __future__ import annotations

import os
from typing import Dict, List, Optional, Tuple

import numpy as np
from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QGridLayout, QLabel, QComboBox, QSpinBox, QCheckBox,
                             QPushButton, QToolButton, QListWidget, QListWidgetItem, QScrollArea, QFrame, QFileDialog)
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

import experiment_plots as ep
from overview_tab import InteractivePanel, PanelGroup, ACCENT

PALETTE = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
           "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5", "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5"]
LAYOUTS = {1: (1, 1), 2: (2, 1), 3: (3, 1), 4: (2, 2), 6: (2, 3)}      # panels → (rows, cols)
OBS_SPECIES = {"NO", "NO2", "O3", "CO", "HONO", "OH", "HO2", "RO2"}      # keys of OBS_MAP that are model species


class CorrelationPanel(QWidget):
    """X species against Y species: scatter coloured by hour (one Y) or by species (several),
    least-squares line per pair, hover readout of the nearest point, wheel zoom, double-click reset."""

    def __init__(self, parent=None):
        super().__init__(parent)
        lay = QVBoxLayout(self); lay.setContentsMargins(0, 0, 0, 0)
        self.fig = Figure(figsize=(5.0, 2.6)); self.canvas = FigureCanvas(self.fig); self.canvas.setMinimumHeight(190)
        self.ax = self.fig.add_subplot(111); lay.addWidget(self.canvas)
        self._pts: List[Tuple[str, np.ndarray, np.ndarray, np.ndarray]] = []    # (label, x, y, hours)
        self._text = None; self._home = None
        self.canvas.mpl_connect("motion_notify_event", self._hover)
        self.canvas.mpl_connect("axes_leave_event", lambda _e: self._clear())
        self.canvas.mpl_connect("scroll_event", self._scroll)
        self.canvas.mpl_connect("button_press_event", self._press)

    def plot(self, x_label: str, x: np.ndarray, ys: List[Tuple[str, np.ndarray]], hours: np.ndarray, title: str):
        ax = self.ax; ax.clear(); self._pts = []
        if len(ys) == 1:
            lab, y = ys[0]
            ok = np.isfinite(x) & np.isfinite(y)
            sc = ax.scatter(x[ok], y[ok], c=hours[ok], cmap="viridis", s=9, alpha=0.8, label=lab)
            cb = self.fig.colorbar(sc, ax=ax, pad=0.01); cb.set_label("hour (UTC)", fontsize=7); cb.ax.tick_params(labelsize=6.5)
            self._pts.append((lab, x[ok], y[ok], hours[ok])); self._fit(ax, x[ok], y[ok], lab, "#c0392b")
            ax.set_ylabel(lab, fontsize=7.5, labelpad=2)
        else:
            for i, (lab, y) in enumerate(ys):
                ok = np.isfinite(x) & np.isfinite(y); c = PALETTE[i % len(PALETTE)]
                ax.scatter(x[ok], y[ok], color=c, s=8, alpha=0.7, label=lab)
                self._pts.append((lab, x[ok], y[ok], hours[ok])); self._fit(ax, x[ok], y[ok], lab, c)
            ax.set_ylabel("y", fontsize=7.5, labelpad=2)
        ax.set_xlabel(x_label, fontsize=7.5, labelpad=1)
        ax.text(0.01, 0.97, title, transform=ax.transAxes, fontsize=8.5, color=ACCENT, fontweight="bold", va="top", ha="left")
        ax.tick_params(labelsize=7, pad=2, length=2); ax.grid(True, alpha=0.25)
        if ax.get_legend_handles_labels()[0]:
            ax.legend(fontsize=6.5, loc="lower right", framealpha=0.8, borderpad=0.3, labelspacing=0.2)
        self.fig.subplots_adjust(left=0.12, right=0.98 if len(ys) != 1 else 0.9, top=0.97, bottom=0.17)
        self._home = (ax.get_xlim(), ax.get_ylim()); self.canvas.draw_idle()

    @staticmethod
    def _fit(ax, x, y, lab, color):
        if len(x) < 3 or np.ptp(x) == 0:
            return
        a, b = np.polyfit(x, y, 1)
        r = np.corrcoef(x, y)[0, 1]
        xx = np.array([np.nanmin(x), np.nanmax(x)])
        ax.plot(xx, a * xx + b, color=color, lw=1, ls="--", label=f"{lab}: slope {a:.3g}, r² {r * r:.2f}")

    # ── interaction ──
    def _hover(self, event):
        if event.inaxes is not self.ax or event.xdata is None or not self._pts:
            return
        best = None
        xr = np.ptp(self.ax.get_xlim()) or 1.0; yr = np.ptp(self.ax.get_ylim()) or 1.0
        for lab, x, y, h in self._pts:
            if len(x) == 0:
                continue
            d = ((x - event.xdata) / xr) ** 2 + ((y - event.ydata) / yr) ** 2
            i = int(np.argmin(d))
            if best is None or d[i] < best[0]:
                best = (d[i], lab, x[i], y[i], h[i])
        self._clear(draw=False)
        if best is not None and best[0] < 1e-3:
            _d, lab, xv, yv, hv = best
            hh, mm = int(hv), int(round((hv - int(hv)) * 60))
            self._text = self.ax.annotate(f"{hh:02d}:{mm:02d}  {lab}\nx {xv:.3g}, y {yv:.3g}", xy=(xv, yv), xytext=(6, 6), textcoords="offset points",
                                          fontsize=6.5, bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#aaa", alpha=0.9))
        self.canvas.draw_idle()

    def _clear(self, draw=True):
        if self._text is not None:
            try:
                self._text.remove()
            except Exception:
                pass
            self._text = None
        if draw:
            self.canvas.draw_idle()

    def _scroll(self, event):
        if event.inaxes is not self.ax or event.xdata is None:
            return
        f = 0.8 if event.button == "up" else 1.25
        for get, set_, v in ((self.ax.get_xlim, self.ax.set_xlim, event.xdata), (self.ax.get_ylim, self.ax.set_ylim, event.ydata)):
            lo, hi = get(); set_(v - (v - lo) * f, v + (hi - v) * f)
        self.canvas.draw_idle()

    def _press(self, event):
        if event.dblclick and self._home is not None:
            self.ax.set_xlim(*self._home[0]); self.ax.set_ylim(*self._home[1]); self.canvas.draw_idle()


class TimeseriesTab(QWidget):
    """Grid of linked time-series / correlation panels fed from the species tree."""

    def __init__(self, viewer, parent=None):
        super().__init__(parent)
        self.v = viewer
        self.n = 2
        self.specs: List[dict] = [self._blank() for _ in range(6)]
        self.group = PanelGroup()
        self._obs: Optional[Dict] = None
        self._obs_note = ""
        lay = QVBoxLayout(self); lay.setContentsMargins(4, 2, 4, 2); lay.setSpacing(2)

        bar = QHBoxLayout(); bar.setSpacing(6)
        bar.addWidget(QLabel("Panels"))
        self.layout_combo = QComboBox(); self.layout_combo.addItems([str(k) for k in LAYOUTS]); self.layout_combo.setCurrentText("2")
        self.layout_combo.currentTextChanged.connect(self._layout_changed); bar.addWidget(self.layout_combo)
        bar.addWidget(QLabel("target"))
        self.target = QSpinBox(); self.target.setRange(1, 2); self.target.setToolTip("the panel that 'Add' and 'Remove' act on")
        self.target.valueChanged.connect(self._refresh_contents); bar.addWidget(self.target)
        bar.addWidget(QLabel("unit"))
        self.unit = QComboBox(); self.unit.addItems(["ppb", "native"]); bar.addWidget(self.unit)
        self.log = QCheckBox("log"); self.log.toggled.connect(self._log_toggled); bar.addWidget(self.log)
        add = QPushButton("Add to panel"); add.setToolTip("add the species selected in the tree to the target panel"); add.clicked.connect(self.add_selected)
        bar.addWidget(add)
        corr = QPushButton("Correlation"); corr.setToolTip("target panel = first selected species (X) against the other selected species (Y)")
        corr.clicked.connect(self.correlation_selected); bar.addWidget(corr)
        bar.addWidget(_tool("Remove", "remove the entries selected in the contents list from the target panel", self.remove_selected))
        bar.addWidget(_tool("Clear panel", "empty the target panel", self.clear_panel))
        bar.addWidget(_tool("Clear all", "empty every panel", self.clear_all))
        self.obs_box = QCheckBox("observations"); self.obs_box.setChecked(True); self.obs_box.setToolTip("SAPHIR observations of the day as dots, same colour as the model line")
        self.obs_box.toggled.connect(self.rebuild); bar.addWidget(self.obs_box)
        bar.addWidget(_tool("Reset view", "full range in all panels (or double-click)", lambda: self.group.reset()))
        bar.addWidget(_tool("Save image", "save the panels as a PNG", self._save_png))
        bar.addStretch()
        hint = QLabel("wheel: zoom time · Ctrl+wheel / right-drag: zoom values · left-drag: pan · double-click: reset · click legend: hide")
        hint.setStyleSheet("color: #888; font-size: 10px;"); bar.addWidget(hint)
        lay.addLayout(bar)

        row2 = QHBoxLayout(); row2.setSpacing(6)
        self.contents_label = QLabel("panel 1:"); row2.addWidget(self.contents_label)
        self.contents = QListWidget(); self.contents.setFlow(QListWidget.Flow.LeftToRight); self.contents.setFixedHeight(26)
        self.contents.setSelectionMode(QListWidget.SelectionMode.ExtendedSelection); self.contents.itemDoubleClicked.connect(lambda _i: self.remove_selected())
        self.contents.setToolTip("entries of the target panel — double-click one to remove it"); row2.addWidget(self.contents, stretch=1)
        self.status = QLabel(""); self.status.setStyleSheet("color: #666; font-size: 10px;"); row2.addWidget(self.status)
        lay.addLayout(row2)

        self.area = QScrollArea(); self.area.setWidgetResizable(True); self.area.setFrameShape(QFrame.Shape.NoFrame)
        self.inner = QWidget(); self.grid = QGridLayout(self.inner); self.grid.setSpacing(2); self.grid.setContentsMargins(0, 0, 0, 0)
        self.area.setWidget(self.inner); lay.addWidget(self.area, stretch=1)

    # ── model of the panels ──
    @staticmethod
    def _blank() -> dict:
        return dict(kind="time", entries=[], log=False, x=None)      # entries: [(name, unit)]; x: name for correlations

    def _target(self) -> dict:
        return self.specs[self.target.value() - 1]

    def _layout_changed(self, text: str):
        self.n = int(text); self.target.setRange(1, self.n)
        self.rebuild()

    def _log_toggled(self, on: bool):
        self._target()["log"] = on; self.rebuild()

    # ── public API (tree, defaults, other tabs) ──
    def add_species(self, names: List[str], unit: Optional[str] = None, panel: Optional[int] = None, log: Optional[bool] = None):
        """Add model species to a panel (1-based; default the target panel)."""
        idx = (panel or self.target.value()) - 1
        spec = self.specs[idx]
        if spec["kind"] != "time":
            spec.update(self._blank())
        u = unit or self.unit.currentText()
        for n in names:
            if n in self.v.df.columns and (n, u) not in spec["entries"]:
                spec["entries"].append((n, u))
        if log is not None:
            spec["log"] = log
        self.rebuild()

    def set_correlation(self, x: str, ys: List[str], panel: Optional[int] = None):
        idx = (panel or self.target.value()) - 1
        u = self.unit.currentText()
        self.specs[idx] = dict(kind="corr", entries=[(y, u) for y in ys if y in self.v.df.columns], log=False, x=x if x in self.v.df.columns else None)
        self.rebuild()

    def _selected_names(self) -> List[str]:
        tree = getattr(self.v, "var_tree", None)
        if tree is None:
            return []
        out = []
        for it in tree.selectedItems():
            if it.parent() is None or it.childCount() != 0:
                continue
            key = it.data(0, Qt.ItemDataRole.UserRole)
            name = it.text(0)
            if key and key in getattr(self.v, "_inst_data", {}):
                m = self.v._inst_data[key].get("mcm_match") or ""
                if m and m in self.v.df.columns:
                    name = m
                else:
                    continue
            if name in self.v.df.columns and name not in out:
                out.append(name)
        return out

    def add_selected(self):
        names = self._selected_names()
        if names:
            self.add_species(names)
            try:
                self.v.var_tree.clearSelection()
            except Exception:
                pass

    def correlation_selected(self):
        names = self._selected_names()
        if len(names) >= 2:
            self.set_correlation(names[0], names[1:])
        else:
            self.status.setText("select the X species first, then at least one Y species")

    def remove_selected(self):
        spec = self._target(); rows = sorted({self.contents.row(i) for i in self.contents.selectedItems()}, reverse=True)
        for r in rows:
            if spec["kind"] == "corr" and r == 0:
                spec["x"] = None
            else:
                k = r - (1 if spec["kind"] == "corr" else 0)
                if 0 <= k < len(spec["entries"]):
                    spec["entries"].pop(k)
        self.rebuild()

    def clear_panel(self):
        self.specs[self.target.value() - 1] = self._blank(); self.rebuild()

    def clear_all(self):
        self.specs = [self._blank() for _ in range(6)]; self.rebuild()

    def apply_defaults(self, left: List[str], right: List[str]):
        """Quick-plot defaults of the TOML: panel 1 = left list in ppb, panel 2 = right list native, log."""
        self.specs = [self._blank() for _ in range(6)]
        self.specs[0]["entries"] = [(s, "ppb") for s in left if s in self.v.df.columns]
        self.specs[1]["entries"] = [(s, "native") for s in right if s in self.v.df.columns]
        self.specs[1]["log"] = bool(self.specs[1]["entries"])
        self.layout_combo.setCurrentText("2" if self.specs[1]["entries"] else "1")
        self.rebuild()

    # ── data ──
    def _series(self, name: str, unit: str):
        df = self.v.df
        y = df[name].values.astype(float)
        if unit == "ppb" and "M" in df.columns and not name.upper().startswith("SAPHIR_"):
            y = y / df["M"].values.astype(float) * 1e9
        return df.index.values.astype(float), y

    def _observations(self) -> Dict:
        if self._obs is not None:
            return self._obs
        self._obs = {}
        try:
            spec = ep.spec_for_run(self.v.csv_path) if getattr(self.v, "csv_path", "") else None
            root = spec.conditions.archive_root if spec else r"Z:\IEK8-SAPHIR\data"
            if spec and spec.conditions.date:
                date = spec.conditions.date
            else:
                import datetime as _dt
                date = _dt.datetime.fromtimestamp(self.v._midnight_unix, _dt.timezone.utc).strftime("%Y-%m-%d")
            note: List[str] = []
            self._obs = ep.load_observations(root, date, sorted(OBS_SPECIES), self.v._midnight_unix,
                                             cache_dir=os.path.dirname(os.path.abspath(self.v.csv_path)) if getattr(self.v, "csv_path", "") else "",
                                             status=note)
            self._obs_note = note[0] if note else ""
        except Exception as e:
            self._obs_note = f"observations unavailable ({e})"
        return self._obs

    # ── build ──
    def _refresh_contents(self):
        spec = self._target(); self.contents.clear()
        self.contents_label.setText(f"panel {self.target.value()}:")
        if spec["kind"] == "corr":
            self.contents.addItem(QListWidgetItem(f"X = {spec['x'] or '?'}"))
            for n, u in spec["entries"]:
                self.contents.addItem(QListWidgetItem(f"Y = {n} [{u}]"))
        else:
            for n, u in spec["entries"]:
                self.contents.addItem(QListWidgetItem(f"{n} [{u}]"))
        self.log.blockSignals(True); self.log.setChecked(spec["log"]); self.log.blockSignals(False)

    def rebuild(self):
        while self.grid.count():
            w = self.grid.takeAt(0).widget()
            if w is not None:
                w.deleteLater()
        self.group = PanelGroup()
        rows, cols = LAYOUTS[self.n]
        obs = self._observations() if self.obs_box.isChecked() else {}
        hv = self._hv_spans()
        for i in range(self.n):
            spec = self.specs[i]
            r, c = i // cols, i % cols
            if spec["kind"] == "corr":
                w = self._build_corr(i, spec)
            else:
                w = InteractivePanel(self._title(i, spec), self.group, log=spec["log"], xlabel=(r == rows - 1))
                self._fill_time(w, spec, obs, hv)
            self.grid.addWidget(w, r, c)
        try:
            t = self.v.df.index.values.astype(float); lo, hi = float(np.nanmin(t)), float(np.nanmax(t))
            if hi > lo:
                self.group.home = (lo, hi); self.group.set_xlim(lo, hi)
        except Exception:
            pass
        self._refresh_contents()
        self.status.setText(("observations " + self._obs_note) if self.obs_box.isChecked() and self._obs_note else "")

    def _title(self, i: int, spec: dict) -> str:
        names = [n for n, _u in spec["entries"]]
        s = ", ".join(names[:6]) + (" …" if len(names) > 6 else "")
        return f"{i + 1}  {s}" if s else f"{i + 1}  (empty — select species in the tree and press 'Add to panel')"

    def _hv_spans(self) -> List[Tuple[float, float]]:
        df = self.v.df
        col = next((c for c in df.columns if c.upper() in ("SAPHIR_HV", "HV")), None)
        if col is None:
            return []
        t = df.index.values.astype(float); hv = df[col].values.astype(float) > 0.5
        spans, start = [], None
        for ti, on in zip(t, hv):
            if on and start is None:
                start = ti
            elif not on and start is not None:
                spans.append((start, ti)); start = None
        if start is not None:
            spans.append((start, float(t[-1])))
        return spans

    def _fill_time(self, panel: InteractivePanel, spec: dict, obs: Dict, hv_spans):
        ax = panel.ax
        units = {u for _n, u in spec["entries"]}
        for lo, hi in hv_spans:
            ax.axvspan(lo, hi, color="#fffacd", alpha=0.35, zorder=0, linewidth=0)
        for i, (name, unit) in enumerate(spec["entries"]):
            t, y = self._series(name, unit); color = PALETTE[i % len(PALETTE)]
            if spec["log"]:
                y = np.where(y > 0, y, np.nan)
            ax.plot(t, y, color=color, lw=1.3, label=f"{name} [{'ppb' if unit == 'ppb' else 'native'}]")
            key = name.upper()
            if key in obs:
                h, v, lab, kind = obs[key]
                if (unit == "ppb") == (kind == "ppb"):               # only when the units agree
                    if spec["log"]:
                        v = np.where(v > 0, v, np.nan)
                    ax.plot(h, v, ".", ms=2.5, color=color, alpha=0.5, label=lab)
        ylabel = "ppb" if units == {"ppb"} else ("native units" if units == {"native"} else "mixed units")
        if spec["log"] and spec["entries"]:                       # hide the floor of the flushed chamber (1e-7 …)
            top = max((np.nanmax(self._series(n, u)[1]) for n, u in spec["entries"]), default=0.0)
            if top > 0:
                ax.set_ylim(bottom=top * 1e-5)
        panel.finish(ylabel)

    def _build_corr(self, i: int, spec: dict) -> QWidget:
        w = CorrelationPanel()
        x = spec.get("x")
        if not x or not spec["entries"]:
            w.ax.text(0.5, 0.5, "correlation: choose X then Y species in the tree and press 'Correlation'", ha="center", va="center", fontsize=8, color="#888", transform=w.ax.transAxes)
            w.canvas.draw_idle(); return w
        u = spec["entries"][0][1]
        t, xv = self._series(x, u)
        ys = [(f"{n} [{un}]", self._series(n, un)[1]) for n, un in spec["entries"]]
        w.plot(f"{x} [{u}]", xv, ys, t, f"{i + 1}  {x} vs " + ", ".join(n for n, _u in spec["entries"]))
        return w

    def _save_png(self):
        d = os.path.dirname(os.path.abspath(self.v.csv_path)) if getattr(self.v, "csv_path", "") else ""
        path, _ = QFileDialog.getSaveFileName(self, "Save the panels", os.path.join(d, "timeseries.png"), "PNG image (*.png)")
        if path:
            self.inner.grab().save(path); self.status.setText(f"saved {path}")


def _tool(text: str, tip: str, slot) -> QToolButton:
    b = QToolButton(); b.setText(text); b.setToolTip(tip); b.clicked.connect(slot); b.setAutoRaise(True)
    return b
