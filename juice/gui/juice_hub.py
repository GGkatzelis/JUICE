"""
juice_hub.py — the JUICE start window: the logo in the middle, one circle for each of the two worlds.

    Chambers  →  the chamber menu (new / saved experiment, study, FACSIMILE import, results)
    Field     →  the field menu   (new / saved field case, examples, results)

The hub stays open: a circle opens its menu as a window of its own, and every choice made there
starts its own process (juice_start.launch).  Clicking a circle whose menu is already open brings
that menu to the front.  Everything is drawn with QPainter in the colours of the logo, so there is
no image besides the logo itself.
"""
from __future__ import annotations

import math
import os

from PyQt6.QtCore import QPointF, QRectF, Qt
from PyQt6.QtGui import (QBrush, QColor, QFont, QLinearGradient, QPainter, QPainterPath, QPen, QPixmap,
                         QPolygonF, QRadialGradient)
from PyQt6.QtWidgets import QPushButton, QWidget

# the colours of the logo
NAVY, NAVY_DEEP, NAVY_LIGHT = QColor("#17486f"), QColor("#0b2236"), QColor("#2a6da0")
ORANGE, SKY, WHITE, DIM = QColor("#f59a23"), QColor("#8fc1e6"), QColor("#ffffff"), QColor("#9db4c8")

WIN_W, WIN_H = 880, 620
CX, CY = WIN_W // 2, 286
LOGO_R, ORBIT_R, SAT_R = 120, 286, 88

WORLDS = [
    dict(key="chamber", title="Chambers", angle=270, glyph="chamber", menu="StartWindow",
         text="Simulation chambers: experiments,\nstudies, FACSIMILE imports, results"),
    dict(key="field", title="Field", angle=90, glyph="field", menu="FieldMenu",
         text="Air outside: a site, a plume or\nan air parcel along a trajectory"),
]


def _xy(angle_deg: float, radius: float):
    """Clock angle (0 = top, clockwise) → window coordinates."""
    a = math.radians(angle_deg - 90)
    return CX + radius * math.cos(a), CY + radius * math.sin(a)


def _font(size: int, bold: bool = False) -> QFont:
    f = QFont("Segoe UI", size)
    f.setBold(bold)
    return f


class Satellite(QPushButton):
    """One round button of the hub: a small drawing and the name of the world it opens."""

    def __init__(self, world: dict, parent=None):
        super().__init__(parent)
        self.world = world
        self._hover = False
        d = 2 * (SAT_R + 22)                           # room for the glow around the circle
        self.setFixedSize(d, d)
        self.setCursor(Qt.CursorShape.PointingHandCursor)
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self.setAccessibleName(world["title"])
        self.setToolTip(world["text"].replace("\n", " "))
        self.setFlat(True)
        self.setStyleSheet("background: transparent; border: none;")

    def title(self) -> str:
        return self.world["title"]

    def enterEvent(self, e):
        self._hover = True; self.update()

    def leaveEvent(self, e):
        self._hover = False; self.update()

    def hitButton(self, pos) -> bool:                  # only the circle is clickable, not the square around it
        c = self.width() / 2
        return (pos.x() - c) ** 2 + (pos.y() - c) ** 2 <= SAT_R ** 2

    def paintEvent(self, _):
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing)
        c = self.width() / 2
        mid = QPointF(c, c)
        lit = self._hover or self.hasFocus()

        if lit:
            glow = QRadialGradient(mid, SAT_R + 22)
            glow.setColorAt(0.74, QColor(245, 154, 35, 110))
            glow.setColorAt(1.0, QColor(245, 154, 35, 0))
            p.setPen(Qt.PenStyle.NoPen); p.setBrush(QBrush(glow))
            p.drawEllipse(mid, SAT_R + 22, SAT_R + 22)

        fill = QRadialGradient(QPointF(c - SAT_R * 0.3, c - SAT_R * 0.35), SAT_R * 1.5)
        fill.setColorAt(0.0, NAVY_LIGHT if lit else QColor("#215c8c"))
        fill.setColorAt(1.0, NAVY)
        p.setBrush(QBrush(fill))
        p.setPen(QPen(ORANGE if lit else SKY, 2.5 if lit else 1.5))
        p.drawEllipse(mid, SAT_R, SAT_R)

        p.save(); p.translate(c, c - 24)
        (_draw_chamber if self.world["glyph"] == "chamber" else _draw_field)(p)
        p.restore()

        f = _font(15, True); p.setFont(f); p.setPen(WHITE)
        p.drawText(QRectF(0, c + 12, self.width(), 30), Qt.AlignmentFlag.AlignHCenter | Qt.AlignmentFlag.AlignTop, self.title())


def _draw_chamber(p: QPainter):
    """A chamber seen from the side: a box on feet, light from above, a few molecules inside."""
    pen = QPen(WHITE, 2.2); pen.setCapStyle(Qt.PenCapStyle.RoundCap); pen.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
    p.setPen(pen); p.setBrush(Qt.BrushStyle.NoBrush)
    p.drawRoundedRect(QRectF(-30, -16, 60, 38), 8, 8)
    p.drawLine(QPointF(-18, 22), QPointF(-18, 28)); p.drawLine(QPointF(18, 22), QPointF(18, 28))
    p.setPen(QPen(ORANGE, 2.2, cap=Qt.PenCapStyle.RoundCap))
    for x in (-16, 0, 16):                             # light
        p.drawLine(QPointF(x, -30), QPointF(x, -23))
    p.setPen(QPen(SKY, 1.4))
    pts = [QPointF(-14, 8), QPointF(-2, -3), QPointF(12, 6)]
    p.drawPolyline(QPolygonF(pts))
    p.setPen(Qt.PenStyle.NoPen)
    for q, col, r in zip(pts, (WHITE, ORANGE, WHITE), (3.6, 4.4, 3.6)):
        p.setBrush(col); p.drawEllipse(q, r, r)


def _draw_field(p: QPainter):
    """Outside: the sun over two hills, wind above them."""
    p.setPen(Qt.PenStyle.NoPen); p.setBrush(ORANGE)
    p.drawEllipse(QPointF(17, -17), 9, 9)
    pen = QPen(WHITE, 2.2); pen.setCapStyle(Qt.PenCapStyle.RoundCap); pen.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
    p.setPen(pen); p.setBrush(Qt.BrushStyle.NoBrush)
    hills = QPainterPath(QPointF(-34, 24))
    hills.cubicTo(QPointF(-24, -2), QPointF(-10, -2), QPointF(0, 14))
    hills.cubicTo(QPointF(8, 2), QPointF(22, 2), QPointF(34, 24))
    p.drawPath(hills)
    p.setPen(QPen(SKY, 1.8, cap=Qt.PenCapStyle.RoundCap))
    for x0, x1, y in ((-30, -8, -22), (-24, 0, -14)):  # wind
        p.drawLine(QPointF(x0, y), QPointF(x1, y))


def _hexagon(cx: float, cy: float, r: float, turn: float = 0.0) -> QPolygonF:
    return QPolygonF([QPointF(cx + r * math.cos(math.radians(60 * k + turn)), cy + r * math.sin(math.radians(60 * k + turn)))
                      for k in range(6)])


# faint ring structures in the background, as in the logo: (x, y, radius, turn in degrees)
_RINGS = [(92, 96, 34, 30), (151, 130, 34, 30), (770, 84, 30, 0), (806, 520, 36, 30), (60, 470, 28, 0), (108, 498, 28, 0),
          (640, 590, 24, 30), (300, 40, 22, 0), (560, 36, 26, 30)]


class Hub(QWidget):
    """The start window.  `toml` is handed on to every window opened from here."""

    def __init__(self, toml: str = ""):
        super().__init__()
        from juice.gui import juice_theme as theme
        self.toml = toml
        self._menus = {}                               # world key → its open menu window
        self.setWindowTitle("JUICE"); self.setWindowIcon(theme.app_icon())
        self.setFixedSize(WIN_W, WIN_H)
        self._logo = self._logo_pixmap(theme.LOGO)
        self.satellites = []
        for w in WORLDS:
            b = Satellite(w, self)
            x, y = _xy(w["angle"], ORBIT_R)
            b.move(int(x - b.width() / 2), int(y - b.height() / 2))
            b.clicked.connect(lambda _checked=False, world=w: self.open_world(world["key"]))
            self.satellites.append(b)

    def titles(self):
        return [b.title() for b in self.satellites]

    def _logo_pixmap(self, path: str) -> QPixmap:
        if not os.path.isfile(path):
            return QPixmap()
        dpr = self.devicePixelRatioF() or 1.0
        side = int(2 * LOGO_R * dpr)
        pix = QPixmap(path).scaled(side, side, Qt.AspectRatioMode.KeepAspectRatio, Qt.TransformationMode.SmoothTransformation)
        pix.setDevicePixelRatio(dpr)
        return pix

    # ── painting ──
    def paintEvent(self, _):
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing)
        bg = QRadialGradient(QPointF(CX, CY), WIN_W * 0.62)
        bg.setColorAt(0.0, QColor("#123a5a")); bg.setColorAt(1.0, NAVY_DEEP)
        p.fillRect(self.rect(), QBrush(bg))

        p.setBrush(Qt.BrushStyle.NoBrush); p.setPen(QPen(QColor(143, 193, 230, 34), 1.3))
        for x, y, r, turn in _RINGS:
            p.drawPolygon(_hexagon(x, y, r, turn))

        p.setPen(QPen(QColor(143, 193, 230, 70), 1, Qt.PenStyle.DashLine))
        p.drawEllipse(QPointF(CX, CY), ORBIT_R, ORBIT_R)
        for w in WORLDS:                               # spokes from the logo to the circles
            a, b = QPointF(*_xy(w["angle"], LOGO_R + 10)), QPointF(*_xy(w["angle"], ORBIT_R - SAT_R - 4))
            g = QLinearGradient(a, b)
            g.setColorAt(0.0, QColor(245, 154, 35, 200)); g.setColorAt(1.0, QColor(143, 193, 230, 200))
            p.setPen(QPen(QBrush(g), 2)); p.drawLine(a, b)

        halo = QRadialGradient(QPointF(CX, CY), LOGO_R + 36)
        halo.setColorAt(0.70, QColor(245, 154, 35, 70)); halo.setColorAt(1.0, QColor(245, 154, 35, 0))
        p.setPen(Qt.PenStyle.NoPen); p.setBrush(QBrush(halo))
        p.drawEllipse(QPointF(CX, CY), LOGO_R + 36, LOGO_R + 36)
        if not self._logo.isNull():
            clip = QPainterPath(); clip.addEllipse(QPointF(CX, CY), LOGO_R - 1, LOGO_R - 1)
            p.save(); p.setClipPath(clip)
            p.drawPixmap(QPointF(CX - LOGO_R, CY - LOGO_R), self._logo)
            p.restore()
        else:
            p.setBrush(NAVY); p.drawEllipse(QPointF(CX, CY), LOGO_R, LOGO_R)
            p.setPen(WHITE); p.setFont(_font(26, True))
            p.drawText(QRectF(CX - LOGO_R, CY - LOGO_R, 2 * LOGO_R, 2 * LOGO_R), Qt.AlignmentFlag.AlignCenter, "JUICE")
        p.setBrush(Qt.BrushStyle.NoBrush); p.setPen(QPen(WHITE, 2.5))
        p.drawEllipse(QPointF(CX, CY), LOGO_R, LOGO_R)

        p.setFont(_font(10)); p.setPen(DIM)                # what each circle leads to, under it
        for w in WORLDS:
            x, y = _xy(w["angle"], ORBIT_R)
            p.drawText(QRectF(x - 150, y + SAT_R + 16, 300, 44), Qt.AlignmentFlag.AlignHCenter | Qt.AlignmentFlag.AlignTop, w["text"])

        p.setPen(WHITE); p.setFont(_font(13, True))
        p.drawText(QRectF(0, WIN_H - 92, WIN_W, 26), Qt.AlignmentFlag.AlignCenter, "Jülich Utility for Integrating Chemical Equations")
        p.setPen(DIM); p.setFont(_font(9))
        p.drawText(QRectF(0, WIN_H - 62, WIN_W, 20), Qt.AlignmentFlag.AlignCenter, self._footer())

    @staticmethod
    def _footer() -> str:
        try:
            from juice.version import __version__ as v
        except Exception:
            v = ""
        return (f"version {v}   ·   " if v else "") + "this window stays open: every choice starts in a window of its own"

    # ── actions ──
    def open_world(self, key: str):
        """The menu of a world; an open one is brought to the front."""
        from juice.gui import juice_start as js
        from juice.gui import juice_theme as theme
        old = self._menus.get(key)
        if old is not None and old.isVisible():
            old.raise_(); old.activateWindow()
            return old
        cls = getattr(js, next(w["menu"] for w in WORLDS if w["key"] == key))
        menu = cls(self.toml, on_back=lambda: self._back(key))
        self._menus[key] = menu
        menu.show(); theme.force_taskbar_icon(menu)
        return menu

    def open_menu(self, cls):
        """The menu given by its class (juice_start.StartWindow or FieldMenu)."""
        return self.open_world(next(w["key"] for w in WORLDS if w["menu"] == cls.__name__))

    def _back(self, key: str):
        menu = self._menus.pop(key, None)
        if menu is not None:
            menu.close()
        self.raise_(); self.activateWindow()
