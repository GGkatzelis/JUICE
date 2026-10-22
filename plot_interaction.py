"""
plot_interaction.py — the Overview's mouse interaction for any matplotlib canvas in the viewer.

    from plot_interaction import attach
    attach(canvas)                 # wheel: zoom time/x · Ctrl+wheel or right-drag: zoom values ·
                                   # left-drag: pan · double-click: reset · hover: read the lines
    attach(canvas, hover=False)    # for canvases that already have their own hover/pick handlers

Works on whatever axes are under the cursor, including twin axes drawn at the same place, and
survives fig.clear()/redraws (the 'home' view is remembered per axes object on first use).  It
stays out of the way while a NavigationToolbar pan/zoom mode is active.
"""
from __future__ import annotations

import math
from typing import Dict, Tuple

import numpy as np


class _Interaction:
    def __init__(self, canvas, hover: bool):
        self.canvas = canvas
        self.hover_on = hover
        self._home: Dict[int, Tuple[tuple, tuple]] = {}
        self._drag = None
        self._marks = []
        canvas.mpl_connect("scroll_event", self._scroll)
        canvas.mpl_connect("button_press_event", self._press)
        canvas.mpl_connect("button_release_event", self._release)
        canvas.mpl_connect("motion_notify_event", self._motion)
        canvas.mpl_connect("axes_leave_event", lambda _e: self._clear())

    # ── helpers ──
    def _toolbar_busy(self) -> bool:
        tb = getattr(self.canvas, "toolbar", None)
        return bool(tb is not None and getattr(tb, "mode", ""))

    def _stack(self, ax):
        """All axes drawn at the same place as `ax` (twinx/twiny partners), ax first."""
        fig = ax.figure
        bb = ax.get_position().bounds
        others = [a for a in fig.axes if a is not ax and np.allclose(a.get_position().bounds, bb)]
        return [ax] + others

    def _remember(self, ax):
        for a in self._stack(ax):
            self._home.setdefault(id(a), (a.get_xlim(), a.get_ylim()))

    @staticmethod
    def _mods():
        try:
            from PyQt6.QtWidgets import QApplication
            from PyQt6.QtCore import Qt
            m = QApplication.keyboardModifiers()
            return bool(m & Qt.KeyboardModifier.ControlModifier), bool(m & Qt.KeyboardModifier.ShiftModifier)
        except Exception:
            return False, False

    @staticmethod
    def _zoom_axis(ax, which: str, v: float, f: float):
        get, set_ = (ax.get_xlim, ax.set_xlim) if which == "x" else (ax.get_ylim, ax.set_ylim)
        lo, hi = get()
        scale = ax.get_xscale() if which == "x" else ax.get_yscale()
        if scale == "log":
            if v <= 0 or lo <= 0:
                return
            lv, llo, lhi = math.log10(v), math.log10(lo), math.log10(hi)
            set_(10 ** (lv - (lv - llo) * f), 10 ** (lv + (lhi - lv) * f))
        else:
            set_(v - (v - lo) * f, v + (hi - v) * f)

    # ── events ──
    def _scroll(self, event):
        ax = event.inaxes
        if ax is None or event.xdata is None or self._toolbar_busy():
            return
        self._remember(ax)
        f = 0.8 if event.button == "up" else 1.25
        ctrl, shift = self._mods()
        if ctrl or shift:
            for a in self._stack(ax):                       # values of every axis drawn here, around the cursor height
                y = a.transData.inverted().transform((0.0, event.y))[1]
                self._zoom_axis(a, "y", y, f)
        else:
            self._zoom_axis(ax, "x", event.xdata, f)        # twins share x automatically
        self.canvas.draw_idle()

    def _press(self, event):
        ax = event.inaxes
        if ax is None or self._toolbar_busy():
            return
        self._remember(ax)
        if event.button == 1:
            if event.dblclick:
                for a in self._stack(ax):
                    h = self._home.get(id(a))
                    if h:
                        a.set_xlim(*h[0]); a.set_ylim(*h[1])
                self.canvas.draw_idle(); return
            pan_y = not ax.lines                              # scatter-only panels pan in both directions
            self._drag = ("pan", ax, event.x, event.y, ax.get_xlim(), [(a, a.get_ylim()) for a in self._stack(ax)], pan_y)
        elif event.button == 3:
            self._drag = ("zoomy", ax, event.x, event.y, ax.get_xlim(), [(a, a.get_ylim()) for a in self._stack(ax)], False)

    def _release(self, _event):
        self._drag = None

    def _motion(self, event):
        if self._drag is not None and event.x is not None:
            kind, ax, x0, y0, xlim, ylims, pan_y = self._drag
            bb = ax.get_window_extent()
            if bb.width <= 0 or bb.height <= 0:
                return
            if kind == "pan":
                lo, hi = xlim
                if ax.get_xscale() == "log" and lo > 0:
                    d = (event.x - x0) / bb.width * (math.log10(hi) - math.log10(lo))
                    ax.set_xlim(10 ** (math.log10(lo) - d), 10 ** (math.log10(hi) - d))
                else:
                    d = (event.x - x0) / bb.width * (hi - lo); ax.set_xlim(lo - d, hi - d)
                if pan_y:
                    for a, (ylo, yhi) in ylims:
                        if a.get_yscale() == "log" and ylo > 0:
                            d = (event.y - y0) / bb.height * (math.log10(yhi) - math.log10(ylo))
                            a.set_ylim(10 ** (math.log10(ylo) - d), 10 ** (math.log10(yhi) - d))
                        else:
                            d = (event.y - y0) / bb.height * (yhi - ylo); a.set_ylim(ylo - d, yhi - d)
            else:                                             # right-drag: up = zoom in, down = zoom out
                f = math.exp(-(event.y - y0) / bb.height * 2.0)
                for a, (ylo, yhi) in ylims:
                    if a.get_yscale() == "log" and ylo > 0:
                        c = 10 ** (0.5 * (math.log10(ylo) + math.log10(yhi)))
                        lv, llo, lhi = math.log10(c), math.log10(ylo), math.log10(yhi)
                        a.set_ylim(10 ** (lv - (lv - llo) * f), 10 ** (lv + (lhi - lv) * f))
                    else:
                        c = 0.5 * (ylo + yhi); a.set_ylim(c - (c - ylo) * f, c + (yhi - c) * f)
            self.canvas.draw_idle(); return
        if self.hover_on:
            self._hover(event)

    def _hover(self, event):
        ax = event.inaxes
        if ax is None or event.xdata is None or self._toolbar_busy():
            return
        x = event.xdata; rows = []
        for a in self._stack(ax):
            for ln in a.get_lines():
                lab = ln.get_label()
                if not ln.get_visible() or lab.startswith("_"):
                    continue
                xd = np.asarray(ln.get_xdata(), dtype=float)
                if xd.size == 0:
                    continue
                i = int(np.argmin(np.abs(xd - x)))
                if abs(float(xd[i]) - x) > 0.05 * (np.ptp(a.get_xlim()) or 1.0):
                    continue
                rows.append(f"{lab}: {float(np.asarray(ln.get_ydata(), dtype=float)[i]):.3g}")
        self._clear(draw=False)
        if not rows:
            self.canvas.draw_idle(); return
        head = f"{int(x):02d}:{int(round((x - int(x)) * 60)) % 60:02d}" if 0 <= x < 48 else f"x = {x:.3g}"
        self._marks = [ax.axvline(x, color="#888", lw=0.8, ls=":"),
                       ax.annotate(head + "\n" + "\n".join(rows[:14]), xy=(x, 0.02), xycoords=("data", "axes fraction"), va="bottom", ha="left",
                                   fontsize=6.5, xytext=(4, 0), textcoords="offset points",
                                   bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#aaa", alpha=0.9))]
        self.canvas.draw_idle()

    def _clear(self, draw=True):
        for m in self._marks:
            try:
                m.remove()
            except Exception:
                pass
        self._marks = []
        if draw:
            self.canvas.draw_idle()


def attach(canvas, hover: bool = True):
    """Give a FigureCanvas the Overview interaction.  Keeps the handler alive on the canvas."""
    if getattr(canvas, "_juice_interaction", None) is None:
        canvas._juice_interaction = _Interaction(canvas, hover)
    return canvas._juice_interaction
