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


class _Dividers:
    """Draggable gaps between side-by-side (or stacked) axes of one figure, like a Qt splitter:
    hover the gap → split cursor, drag → the boundary moves and both neighbours resize.  The chosen
    boundaries survive redraws: they are keyed by gap order and re-applied on every draw."""
    MIN_W = 0.08          # axes narrower than this (figure fraction) are colourbars etc.: not neighbours

    def __init__(self, canvas):
        self.canvas = canvas
        self._drag = None
        self._cursor_set = False
        self.overrides: Dict[str, float] = {}
        canvas.mpl_connect("motion_notify_event", self._motion)
        canvas.mpl_connect("button_press_event", self._press)
        canvas.mpl_connect("button_release_event", self._release)
        canvas.mpl_connect("draw_event", self._on_draw)

    # ── geometry ──
    def _panels(self):
        fig = self.canvas.figure
        out = []
        for ax in fig.axes:
            b = ax.get_position()
            if b.width >= self.MIN_W and b.height >= self.MIN_W and ax.get_visible():
                out.append(ax)
        return out

    def _gaps(self):
        """[(kind, lo, hi, left_group, right_group)] in figure fraction; kind 'h' = vertical divider."""
        axes = self._panels()
        pos = {id(a): a.get_position() for a in axes}
        seen = set(); gaps = []
        for a in axes:
            pa = pos[id(a)]
            for b in axes:
                if a is b:
                    continue
                pb = pos[id(b)]
                # horizontal neighbours: b right of a, vertical overlap
                if pb.x0 >= pa.x1 - 1e-3 and min(pa.y1, pb.y1) - max(pa.y0, pb.y0) > 0.3 * min(pa.height, pb.height):
                    lo, hi = pa.x1, pb.x0
                    if 0 <= hi - lo < 0.25:
                        key = ("h", round(lo, 3), round(hi, 3))
                        if key not in seen:
                            seen.add(key)
                            left = [c for c in axes if abs(pos[id(c)].x1 - lo) < 0.01]
                            right = [c for c in axes if abs(pos[id(c)].x0 - hi) < 0.01]
                            gaps.append(("h", lo, hi, left, right))
                # vertical neighbours: b below a
                if pb.y1 <= pa.y0 + 1e-3 and min(pa.x1, pb.x1) - max(pa.x0, pb.x0) > 0.3 * min(pa.width, pb.width):
                    lo, hi = pb.y1, pa.y0
                    if 0 <= hi - lo < 0.25:
                        key = ("v", round(lo, 3), round(hi, 3))
                        if key not in seen:
                            seen.add(key)
                            below = [c for c in axes if abs(pos[id(c)].y1 - lo) < 0.01]
                            above = [c for c in axes if abs(pos[id(c)].y0 - hi) < 0.01]
                            gaps.append(("v", lo, hi, below, above))
        gaps.sort(key=lambda g: (g[0], g[1]))
        return gaps

    def _fig_xy(self, event):
        w, h = self.canvas.figure.bbox.width, self.canvas.figure.bbox.height
        return event.x / w, event.y / h

    def _hit(self, event):
        if event.x is None or event.y is None:
            return None
        fx, fy = self._fig_xy(event)
        for i, g in enumerate(self._gaps()):
            kind, lo, hi = g[0], g[1], g[2]
            pad = 0.006
            if kind == "h" and lo - pad <= fx <= hi + pad and g[3] and g[4]:
                ys = [a.get_position() for a in g[3] + g[4]]
                if min(p.y0 for p in ys) <= fy <= max(p.y1 for p in ys):
                    return i, g
            if kind == "v" and lo - pad <= fy <= hi + pad and g[3] and g[4]:
                xs = [a.get_position() for a in g[3] + g[4]]
                if min(p.x0 for p in xs) <= fx <= max(p.x1 for p in xs):
                    return i, g
        return None

    # ── apply ──
    @staticmethod
    def _move(g, new_mid):
        kind, lo, hi, first, second = g
        gap = hi - lo
        a, b = new_mid - gap / 2, new_mid + gap / 2
        for ax in first:
            p = ax.get_position()
            if kind == "h":
                if a - p.x0 > 0.05:
                    ax.set_position([p.x0, p.y0, a - p.x0, p.height])
            else:
                if a - p.y0 > 0.05:
                    ax.set_position([p.x0, p.y0, p.width, a - p.y0])
        for ax in second:
            p = ax.get_position()
            if kind == "h":
                if p.x1 - b > 0.05:
                    ax.set_position([b, p.y0, p.x1 - b, p.height])
            else:
                if p.y1 - b > 0.05:
                    ax.set_position([p.x0, b, p.width, p.y1 - b])

    def _on_draw(self, _event):
        if not self.overrides or self._drag is not None:
            return
        gaps = self._gaps()
        changed = False
        for i, g in enumerate(gaps):
            key = f"{g[0]}{i}"
            if key in self.overrides and abs((g[1] + g[2]) / 2 - self.overrides[key]) > 1e-3:
                self._move(g, self.overrides[key]); changed = True
        if changed:
            self.canvas.draw_idle()

    # ── events ──
    def _busy(self):
        tb = getattr(self.canvas, "toolbar", None)
        return bool(tb is not None and getattr(tb, "mode", ""))

    def _set_cursor(self, kind):
        try:
            from PyQt6.QtCore import Qt
            from PyQt6.QtGui import QCursor
            if kind is None:
                if self._cursor_set:
                    self.canvas.unsetCursor(); self._cursor_set = False
            else:
                self.canvas.setCursor(QCursor(Qt.CursorShape.SplitHCursor if kind == "h" else Qt.CursorShape.SplitVCursor)); self._cursor_set = True
        except Exception:
            pass

    def _press(self, event):
        if event.button != 1 or self._busy():
            return
        hit = self._hit(event)
        if hit is not None:
            self._drag = hit
            inter = getattr(self.canvas, "_juice_interaction", None)
            if inter is not None:
                inter._drag = None          # the pan/zoom handler must not act on this drag

    def _release(self, _event):
        if self._drag is not None:
            i, g = self._drag
            self.overrides[f"{g[0]}{i}"] = (g[1] + g[2]) / 2
            self._drag = None

    def _motion(self, event):
        if self._drag is not None:
            if event.x is None:
                return
            i, g = self._drag
            fx, fy = self._fig_xy(event)
            new_mid = fx if g[0] == "h" else fy
            self._move(g, new_mid)
            lo, hi = g[1], g[2]; gap = hi - lo
            self._drag = (i, (g[0], new_mid - gap / 2, new_mid + gap / 2, g[3], g[4]))
            self.canvas.draw_idle(); return
        if self._busy():
            return
        hit = self._hit(event)
        self._set_cursor(hit[1][0] if hit else None)


def attach_dividers(canvas):
    """Draggable dividers between the sub-plots of a figure (see _Dividers)."""
    if getattr(canvas, "_juice_dividers", None) is None:
        canvas._juice_dividers = _Dividers(canvas)
    return canvas._juice_dividers


def attach(canvas, hover: bool = True):
    """Give a FigureCanvas the Overview interaction.  Keeps the handler alive on the canvas."""
    if getattr(canvas, "_juice_interaction", None) is None:
        canvas._juice_interaction = _Interaction(canvas, hover)
    return canvas._juice_interaction
