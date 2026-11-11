"""
plot_interaction.py — the Overview's mouse interaction for any matplotlib canvas in the viewer.

    from plot_interaction import attach
    attach(canvas)                 # wheel: zoom time/x · Ctrl+wheel or right-drag: zoom values ·
                                   # left-drag: pan · double-click: reset · hover: read the lines
    attach(canvas, hover=False)    # for canvases that already have their own hover/pick handlers

    right-click (without dragging): menu with zoom in/out at that point, reset, and linear/log
    for each axis (set_log, popup — also used by the Overview and Timeseries panels)

Works on whatever axes are under the cursor, including twin axes drawn at the same place, and
survives fig.clear()/redraws (the 'home' view is remembered per axes object on first use).  It
stays out of the way while a NavigationToolbar pan/zoom mode is active.
"""
from __future__ import annotations

import math
from typing import Dict, List, Tuple

import numpy as np

CLICK_PX = 4          # a right-button press and release closer than this is a click (menu), not a drag (zoom)


def zoom_axis(ax, which: str, v: float, f: float):
    """Zoom one axis of `ax` by factor f around v in data units (log-aware); f < 1 zooms in."""
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


def is_log(ax, which: str) -> bool:
    return (ax.get_xscale() if which == "x" else ax.get_yscale()) == "log"


def _data_range(ax, which: str, positive: bool):
    """Smallest and largest value on one axis of the data drawn in `ax` (lines and scatter points in
    data coordinates; axvline/axvspan & co. are left out).  For the values only what lies in the current
    x view counts, so a zoomed-in time window is fitted, not the whole day."""
    xlo, xhi = sorted(ax.get_xlim())
    found = []
    pairs = [(ln.get_xdata(), ln.get_ydata()) for ln in ax.get_lines() if ln.get_visible() and ln.get_transform() == ax.transData]
    for c in ax.collections:
        try:
            if c.get_visible() and c.get_offset_transform() == ax.transData:
                off = np.asarray(c.get_offsets(), dtype=float)
                if off.ndim == 2 and off.shape[1] == 2:
                    pairs.append((off[:, 0], off[:, 1]))
        except Exception:
            pass
    for xd, yd in pairs:
        try:
            x = np.asarray(xd, dtype=float).ravel(); y = np.asarray(yd, dtype=float).ravel()
        except (TypeError, ValueError):
            continue
        if x.size != y.size:
            continue
        v = y[(x >= xlo) & (x <= xhi)] if which == "y" else x
        v = v[np.isfinite(v)]
        if positive:
            v = v[v > 0]
        if v.size:
            found.append((float(v.min()), float(v.max())))
    if not found:
        return None
    return min(f[0] for f in found), max(f[1] for f in found)


def set_log(ax, which: str, on: bool) -> None:
    """Switch one axis between linear and log and fit its range to the data on the new scale.  On log,
    values <= 0 are left out and the range ends at 1e-5 of the largest value (else the floor of a
    flushed chamber, 1e-7 ppb and less, squeezes the real data into the top of the plot)."""
    if is_log(ax, which) == on:
        return
    scale = dict(value="log", nonpositive="mask") if on else dict(value="linear")
    (ax.set_xscale if which == "x" else ax.set_yscale)(**scale)
    set_ = ax.set_xlim if which == "x" else ax.set_ylim
    rng = _data_range(ax, which, positive=on)
    if rng is None:                                        # bars, images …: let matplotlib decide
        ax.relim(); ax.autoscale(enable=True, axis=which); return
    lo, hi = rng
    if on:
        lo = max(lo, hi * 1e-5)
        pad = 10 ** (0.04 * max(math.log10(hi / lo), 0.5))
        set_(lo / pad, hi * pad)
    else:
        pad = 0.05 * (hi - lo) or 0.05 * abs(hi) or 1.0
        set_(lo - pad, hi + pad)


def popup(canvas, entries) -> None:
    """Context menu at the mouse.  entries: (text, callback), (text, callback, checked) for a tick
    box, or None for a separator."""
    from PyQt6.QtGui import QCursor
    from PyQt6.QtWidgets import QMenu
    menu = QMenu(canvas)
    for e in entries:
        if e is None:
            menu.addSeparator(); continue
        act = menu.addAction(e[0])
        if len(e) > 2:
            act.setCheckable(True); act.setChecked(bool(e[2]))
        act.triggered.connect(lambda _checked=False, f=e[1]: f())
    menu.exec(QCursor.pos())


def is_click(drag_xy, event) -> bool:
    """True when the button went up (about) where it went down."""
    return event.x is not None and abs(event.x - drag_xy[0]) + abs(event.y - drag_xy[1]) < CLICK_PX


class _Interaction:
    def __init__(self, canvas, hover: bool):
        self.canvas = canvas
        self.hover_on = hover
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
        for a in self._stack(ax):                  # on the axes object: ids are reused after fig.clear()
            if getattr(a, "_juice_home", None) is None:
                a._juice_home = (a.get_xlim(), a.get_ylim())

    @staticmethod
    def _mods():
        try:
            from PyQt6.QtWidgets import QApplication
            from PyQt6.QtCore import Qt
            m = QApplication.keyboardModifiers()
            return bool(m & Qt.KeyboardModifier.ControlModifier), bool(m & Qt.KeyboardModifier.ShiftModifier)
        except Exception:
            return False, False

    _zoom_axis = staticmethod(zoom_axis)

    def _reset(self, ax):
        for a in self._stack(ax):
            h = getattr(a, "_juice_home", None)
            if h:
                a.set_xlim(*h[0]); a.set_ylim(*h[1])
        self.canvas.draw_idle()

    def _menu(self, ax, event):
        """Right-click menu: zoom at the clicked point, reset, linear/log for each axis drawn here."""
        stack = self._stack(ax)
        px, py = event.x, event.y

        def zoom(which, f):
            for a in (stack[:1] if which == "x" else stack):          # twins share x
                xv, yv = a.transData.inverted().transform((px, py))
                self._zoom_axis(a, which, xv if which == "x" else yv, f)
            self.canvas.draw_idle()

        def log(a, which, on):
            set_log(a, which, on)
            for b in stack:                                          # the new view is the one double-click returns to
                b._juice_home = (b.get_xlim(), b.get_ylim())
            self.canvas.draw_idle()

        entries = [("Zoom in here", lambda: zoom("x", 0.5)), ("Zoom out", lambda: zoom("x", 2.0)),
                   ("Zoom in values here", lambda: zoom("y", 0.5)), ("Zoom out values", lambda: zoom("y", 2.0)),
                   ("Reset view   (double-click)", lambda: self._reset(ax)), None]
        right = lambda a: a.yaxis.get_ticks_position() == "right"    # the clicked axes may be the twin
        for a in sorted(stack, key=right):
            name = "Log scale — values" + ("" if len(stack) == 1 else (" (right axis)" if right(a) else " (left axis)"))
            entries.append((name, lambda a=a: log(a, "y", not is_log(a, "y")), is_log(a, "y")))
        entries.append(("Log scale — horizontal axis", lambda: log(ax, "x", not is_log(ax, "x")), is_log(ax, "x")))
        popup(self.canvas, entries)

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
                self._reset(ax); return
            pan_y = not ax.lines                              # scatter-only panels pan in both directions
            self._drag = ("pan", ax, event.x, event.y, ax.get_xlim(), [(a, a.get_ylim()) for a in self._stack(ax)], pan_y)
        elif event.button == 3:
            self._drag = ("zoomy", ax, event.x, event.y, ax.get_xlim(), [(a, a.get_ylim()) for a in self._stack(ax)], False)

    def _release(self, event):
        drag, self._drag = self._drag, None
        if drag is not None and drag[0] == "zoomy" and is_click(drag[2:4], event):
            for a, lim in drag[5]:                            # undo the sub-pixel zoom of the 'drag'
                a.set_ylim(*lim)
            self._menu(drag[1], event)

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
    """Resize the sub-plots of a figure with the mouse, like Qt splitters:

      • the middle of the gap between two neighbours → split cursor, drag moves the boundary and
        both neighbours resize (the gap keeps its width);
      • the edge of one plot (also an outer edge) → resize cursor, drag moves only that edge, the
        other plot stays — this widens or narrows the gap, e.g. for long tick labels.

    Choices survive redraws: they are recorded as operations (keyed by plot rank and edge, or by
    gap order) and re-applied on every draw.  Colourbars (narrow axes) are left alone."""
    MIN_W = 0.08          # axes narrower than this (figure fraction) are colourbars etc.
    PAD = 0.007           # hit distance to an edge, figure fraction
    MIN_SIZE = 0.05

    def __init__(self, canvas):
        self.canvas = canvas
        self._drag = None
        self._cursor_set = False
        self.ops: List[tuple] = []          # ("edge", rank, side, value) | ("gap", kind, index, mid)
        canvas.mpl_connect("motion_notify_event", self._motion)
        canvas.mpl_connect("button_press_event", self._press)
        canvas.mpl_connect("button_release_event", self._release)
        canvas.mpl_connect("draw_event", self._on_draw)

    # ── geometry ──
    def _panels(self):
        fig = self.canvas.figure
        axes = [ax for ax in fig.axes if ax.get_visible() and ax.get_label() != "<colorbar>"
                and ax.get_position().width >= self.MIN_SIZE and ax.get_position().height >= self.MIN_SIZE]
        return sorted(axes, key=lambda a: (round(a.get_position().x0, 3), -round(a.get_position().y0, 3)))

    def _gaps(self, axes):
        pos = {id(a): a.get_position() for a in axes}
        seen = set(); gaps = []
        for a in axes:
            pa = pos[id(a)]
            for b in axes:
                if a is b:
                    continue
                pb = pos[id(b)]
                if pb.x0 >= pa.x1 - 1e-3 and min(pa.y1, pb.y1) - max(pa.y0, pb.y0) > 0.3 * min(pa.height, pb.height):
                    lo, hi = pa.x1, pb.x0
                    if 0 <= hi - lo < 0.6:
                        key = ("h", round(lo, 3), round(hi, 3))
                        if key not in seen:
                            seen.add(key)
                            gaps.append(("h", lo, hi, [c for c in axes if abs(pos[id(c)].x1 - lo) < 0.01], [c for c in axes if abs(pos[id(c)].x0 - hi) < 0.01]))
                if pb.y1 <= pa.y0 + 1e-3 and min(pa.x1, pb.x1) - max(pa.x0, pb.x0) > 0.3 * min(pa.width, pb.width):
                    lo, hi = pb.y1, pa.y0
                    if 0 <= hi - lo < 0.6:
                        key = ("v", round(lo, 3), round(hi, 3))
                        if key not in seen:
                            seen.add(key)
                            gaps.append(("v", lo, hi, [c for c in axes if abs(pos[id(c)].y1 - lo) < 0.01], [c for c in axes if abs(pos[id(c)].y0 - hi) < 0.01]))
        gaps.sort(key=lambda g: (g[0], g[1]))
        return gaps

    def _fig_xy(self, event):
        w, h = self.canvas.figure.bbox.width, self.canvas.figure.bbox.height
        return event.x / w, event.y / h

    def _hit(self, event):
        """('edge', rank, side) | ('gap', index, gap) | None"""
        if event.x is None or event.y is None:
            return None
        fx, fy = self._fig_xy(event)
        axes = self._panels()
        for rank, ax in enumerate(axes):                       # edges first: they take priority over the gap middle
            p = ax.get_position()
            if p.y0 - self.PAD <= fy <= p.y1 + self.PAD:
                if abs(fx - p.x0) <= self.PAD:
                    return ("edge", rank, "left")
                if abs(fx - p.x1) <= self.PAD:
                    return ("edge", rank, "right")
            if p.x0 - self.PAD <= fx <= p.x1 + self.PAD:
                if abs(fy - p.y0) <= self.PAD:
                    return ("edge", rank, "bottom")
                if abs(fy - p.y1) <= self.PAD:
                    return ("edge", rank, "top")
        for i, g in enumerate(self._gaps(axes)):
            kind, lo, hi = g[0], g[1], g[2]
            if kind == "h" and lo <= fx <= hi and g[3] and g[4]:
                ys = [a.get_position() for a in g[3] + g[4]]
                if min(q.y0 for q in ys) <= fy <= max(q.y1 for q in ys):
                    return ("gap", i, g)
            if kind == "v" and lo <= fy <= hi and g[3] and g[4]:
                xs = [a.get_position() for a in g[3] + g[4]]
                if min(q.x0 for q in xs) <= fx <= max(q.x1 for q in xs):
                    return ("gap", i, g)
        return None

    # ── apply ──
    def _set_edge(self, ax, side, value):
        p = ax.get_position()
        if side == "left" and p.x1 - value > self.MIN_SIZE:
            ax.set_position([value, p.y0, p.x1 - value, p.height])
        elif side == "right" and value - p.x0 > self.MIN_SIZE:
            ax.set_position([p.x0, p.y0, value - p.x0, p.height])
        elif side == "bottom" and p.y1 - value > self.MIN_SIZE:
            ax.set_position([p.x0, value, p.width, p.y1 - value])
        elif side == "top" and value - p.y0 > self.MIN_SIZE:
            ax.set_position([p.x0, p.y0, p.width, value - p.y0])

    def _move_gap(self, g, new_mid):
        kind, lo, hi, first, second = g
        gap = hi - lo; a, b = new_mid - gap / 2, new_mid + gap / 2
        for ax in first:
            self._set_edge(ax, "right" if kind == "h" else "top", a)
        for ax in second:
            self._set_edge(ax, "left" if kind == "h" else "bottom", b)

    def _apply_ops(self):
        axes = self._panels()
        for op in self.ops:
            if op[0] == "edge":
                _k, rank, side, value = op
                if rank < len(axes):
                    self._set_edge(axes[rank], side, value)
            else:
                _k, kind, index, mid = op
                gaps = [g for g in self._gaps(axes) if g[0] == kind]
                if index < len(gaps):
                    self._move_gap(gaps[index], mid)

    def _on_draw(self, _event):
        if not self.ops or self._drag is not None:
            return
        before = [tuple(round(v, 4) for v in a.get_position().bounds) for a in self._panels()]
        self._apply_ops()
        after = [tuple(round(v, 4) for v in a.get_position().bounds) for a in self._panels()]
        if before != after:
            self.canvas.draw_idle()

    # ── events ──
    def _busy(self):
        tb = getattr(self.canvas, "toolbar", None)
        return bool(tb is not None and getattr(tb, "mode", ""))

    def _set_cursor(self, shape):
        try:
            from PyQt6.QtCore import Qt
            from PyQt6.QtGui import QCursor
            if shape is None:
                if self._cursor_set:
                    self.canvas.unsetCursor(); self._cursor_set = False
            else:
                self.canvas.setCursor(QCursor(shape)); self._cursor_set = True
        except Exception:
            pass

    def _cursor_for(self, hit):
        from PyQt6.QtCore import Qt
        if hit is None:
            return None
        if hit[0] == "gap":
            return Qt.CursorShape.SplitHCursor if hit[2][0] == "h" else Qt.CursorShape.SplitVCursor
        return Qt.CursorShape.SizeHorCursor if hit[2] in ("left", "right") else Qt.CursorShape.SizeVerCursor

    def _press(self, event):
        if event.button != 1 or self._busy():
            return
        lg = getattr(self.canvas, "_juice_legends", None)
        if lg is not None and lg.hit(event) is not None:
            return                            # a legend is under the cursor: it is moved, not the plot edge
        hit = self._hit(event)
        if hit is not None:
            self._drag = hit
            inter = getattr(self.canvas, "_juice_interaction", None)
            if inter is not None:
                inter._drag = None            # the pan handler must not act on this drag

    def _release(self, event):
        if self._drag is None:
            return
        hit = self._drag; self._drag = None
        if hit[0] == "edge":
            _k, rank, side = hit
            axes = self._panels()
            if rank < len(axes):
                p = axes[rank].get_position()
                value = {"left": p.x0, "right": p.x1, "bottom": p.y0, "top": p.y1}[side]
                self.ops = [o for o in self.ops if not (o[0] == "edge" and o[1] == rank and o[2] == side)] + [("edge", rank, side, value)]
        else:
            _k, index, g = hit
            self.ops = [o for o in self.ops if not (o[0] == "gap" and o[1] == g[0] and o[2] == index)] + [("gap", g[0], index, (g[1] + g[2]) / 2)]

    def _motion(self, event):
        if self._drag is not None:
            if event.x is None:
                return
            fx, fy = self._fig_xy(event)
            hit = self._drag
            if hit[0] == "edge":
                _k, rank, side = hit
                axes = self._panels()
                if rank < len(axes):
                    self._set_edge(axes[rank], side, fx if side in ("left", "right") else fy)
            else:
                _k, index, g = hit
                new_mid = fx if g[0] == "h" else fy
                self._move_gap(g, new_mid)
                gap = g[2] - g[1]
                self._drag = ("gap", index, (g[0], new_mid - gap / 2, new_mid + gap / 2, g[3], g[4]))
            self.canvas.draw_idle(); return
        if self._busy():
            return
        self._set_cursor(self._cursor_for(self._hit(event)))


class _Legends:
    """Drag a legend with the left mouse button; double-click it to put it back where the plot placed it.

    The position is kept in the coordinates of the legend's axes (or of the figure for figure legends),
    remembered per legend (axes rank in the figure, or figure-legend index) and re-applied on every
    draw, so a legend stays where it was put when the plot is redrawn (new metric, new data …).
    While a legend is under the cursor the other drags of the canvas (pan, plot edges) are suppressed;
    a click without moving still reaches pick handlers (legend entries that hide a line)."""

    def __init__(self, canvas, on_grab=None):
        self.canvas = canvas
        self.on_grab = on_grab                # e.g. cancels a panel's own pan
        self._drag = None
        self._cursor = False
        self.ops: Dict[tuple, Tuple[float, float]] = {}       # legend key → lower-left corner in parent fractions
        self._orig: Dict[tuple, tuple] = {}                    # legend key → (loc, bbox_to_anchor) before the first move
        canvas.mpl_connect("button_press_event", self._press)
        canvas.mpl_connect("motion_notify_event", self._motion)
        canvas.mpl_connect("button_release_event", self._release)
        canvas.mpl_connect("draw_event", self._on_draw)

    def _legends(self):
        fig = self.canvas.figure
        out = []
        for i, ax in enumerate(fig.axes):
            lg = ax.get_legend()
            if lg is not None and lg.get_visible():
                out.append((("ax", i), lg, ax))
        for j, lg in enumerate(getattr(fig, "legends", [])):
            if lg.get_visible():
                out.append((("fig", j), lg, None))
        return out

    def hit(self, event):
        if getattr(event, "x", None) is None:
            return None
        for key, lg, ax in reversed(self._legends()):
            try:
                if lg.get_window_extent().contains(event.x, event.y):
                    return key, lg, ax
            except Exception:
                continue
        return None

    def _parent(self, ax):
        return ax.transAxes if ax is not None else self.canvas.figure.transFigure

    def _anchor(self, lg, ax):
        """The legend's bbox_to_anchor in fractions of its parent (None = the parent itself), so a reset still
        fits after the plot was rebuilt with new axes."""
        try:
            bb = lg.get_bbox_to_anchor()
            par = (ax.bbox if ax is not None else self.canvas.figure.bbox)
            if np.allclose(bb.bounds, par.bounds, atol=0.5):
                return None
            (x0, y0), (x1, y1) = self._parent(ax).inverted().transform([(bb.x0, bb.y0), (bb.x1, bb.y1)])
            return (x0, y0, x1 - x0, y1 - y0)
        except Exception:
            return None

    def _place(self, lg, ax, xy):
        lg.set_bbox_to_anchor(None)           # the parent itself: loc (x, y) = lower-left corner in its fractions
        try:
            lg.set_loc(tuple(xy))
        except Exception:
            lg._loc = tuple(xy)

    def _press(self, event):
        if event.button != 1:
            return
        tb = getattr(self.canvas, "toolbar", None)
        if tb is not None and getattr(tb, "mode", ""):
            return
        h = self.hit(event)
        if h is None:
            return
        key, lg, ax = h
        if event.dblclick:                    # back to where the plot put it
            self.ops.pop(key, None)
            if key in self._orig:
                loc, bounds = self._orig.pop(key)
                try:
                    lg.set_bbox_to_anchor(bounds, transform=self._parent(ax) if bounds is not None else None)
                    lg.set_loc(loc)
                except Exception:
                    lg._loc = loc
            self.canvas.draw_idle()
        else:
            bb = lg.get_window_extent()
            self._drag = (key, lg, ax, event.x, event.y, bb.x0, bb.y0, False)
        for other in ("_juice_interaction", "_juice_dividers"):
            o = getattr(self.canvas, other, None)
            if o is not None:
                o._drag = None
        if self.on_grab is not None:
            self.on_grab()

    def _motion(self, event):
        if self._drag is not None:
            if event.x is None:
                return
            key, lg, ax, x0, y0, bx, by, moved = self._drag
            if key not in self._orig:
                self._orig[key] = (getattr(lg, "_loc", 0), self._anchor(lg, ax))
            xy = self._parent(ax).inverted().transform((bx + event.x - x0, by + event.y - y0))
            self._place(lg, ax, xy)
            self._drag = (key, lg, ax, x0, y0, bx, by, True)
            self.canvas.draw_idle()
            return
        h = self.hit(event)
        try:
            from PyQt6.QtCore import Qt
            from PyQt6.QtGui import QCursor
            if h is not None:
                self.canvas.setCursor(QCursor(Qt.CursorShape.OpenHandCursor)); self._cursor = True
            elif self._cursor:
                div = getattr(self.canvas, "_juice_dividers", None)
                if div is None or not div._cursor_set:
                    self.canvas.unsetCursor()
                self._cursor = False
        except Exception:
            pass

    def _release(self, _event):
        if self._drag is None:
            return
        key, lg, ax, _x0, _y0, _bx, _by, moved = self._drag
        self._drag = None
        if moved:
            loc = getattr(lg, "_loc", None)
            if isinstance(loc, tuple):
                self.ops[key] = loc

    def _on_draw(self, _event):
        if not self.ops or self._drag is not None:
            return
        changed = False
        for key, lg, ax in self._legends():
            xy = self.ops.get(key)
            if xy is not None and getattr(lg, "_loc", None) != xy:
                if key not in self._orig:
                    self._orig[key] = (getattr(lg, "_loc", 0), self._anchor(lg, ax))
                self._place(lg, ax, xy); changed = True
        if changed:
            self.canvas.draw_idle()


def attach_legends(canvas, on_grab=None):
    """Draggable legends on a canvas (see _Legends)."""
    if getattr(canvas, "_juice_legends", None) is None:
        canvas._juice_legends = _Legends(canvas, on_grab)
    return canvas._juice_legends


def attach_dividers(canvas):
    """Draggable dividers between the sub-plots of a figure (see _Dividers)."""
    if getattr(canvas, "_juice_dividers", None) is None:
        canvas._juice_dividers = _Dividers(canvas)
    return canvas._juice_dividers


def attach(canvas, hover: bool = True):
    """Give a FigureCanvas the Overview interaction.  Keeps the handler alive on the canvas."""
    if getattr(canvas, "_juice_interaction", None) is None:
        canvas._juice_interaction = _Interaction(canvas, hover)
        attach_legends(canvas)
    return canvas._juice_interaction
