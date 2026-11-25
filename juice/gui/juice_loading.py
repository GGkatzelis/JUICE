"""
juice_loading.py — the small "JUICE is opening …" window shown while another JUICE window starts.

Every JUICE window runs in its own process and needs some seconds before it is on screen (Python, Qt,
matplotlib, and for a study or a result the files themselves).  The window that started it shows this
card at once and closes it when the new window is there:

    start(module, args, cwd, what)   start `python -m module args` and show the card for it

How the card knows: the new process gets the name of a marker file in the environment variable
JUICE_READY_FILE and removes the file once its window is shown (juice_theme.signal_ready, called by
every window right after show()).  The card also closes when that process has ended without a window,
or after TIMEOUT_S.
"""
from __future__ import annotations

import os
import sys
import tempfile
import time
import uuid

from PyQt6.QtCore import QPointF, QProcess, QProcessEnvironment, QRectF, Qt, QTimer
from PyQt6.QtGui import QBrush, QColor, QFont, QGuiApplication, QLinearGradient, QPainter, QPainterPath, QPen, QPixmap
from PyQt6.QtWidgets import QWidget

from juice.gui.juice_theme import READY_ENV

TIMEOUT_S = 120.0
W, H, LOGO = 440, 150, 84
NAVY, NAVY_DEEP, ORANGE, WHITE, DIM = QColor("#17486f"), QColor("#0b2236"), QColor("#f59a23"), QColor("#ffffff"), QColor("#9db4c8")

_open = []                      # the cards on screen (a window without a reference would be collected)


def _alive(pid: int) -> bool:
    """Whether the started process still runs (unknown counts as running)."""
    if not pid:
        return True
    try:
        if sys.platform.startswith("win"):
            import ctypes
            k32 = ctypes.windll.kernel32
            h = k32.OpenProcess(0x00100000, False, int(pid))          # SYNCHRONIZE
            if not h:
                return False
            try:
                return k32.WaitForSingleObject(h, 0) != 0             # 0 = ended
            finally:
                k32.CloseHandle(h)
        os.kill(int(pid), 0)
        return True
    except Exception:
        return sys.platform.startswith("win")


class LoadingWindow(QWidget):
    """The card: logo, what is opening, and a bar that keeps moving until the new window is there."""

    def __init__(self, what: str, marker: str = "", pid: int = 0):
        # an ordinary window (not a tool window): it keeps the application alive while it is the only one on screen
        super().__init__(None, Qt.WindowType.FramelessWindowHint | Qt.WindowType.WindowStaysOnTopHint)
        from juice.gui import juice_theme as theme
        self.what, self.marker, self.pid = what, marker, pid
        self._t0, self._phase = time.monotonic(), 0.0
        self.setAttribute(Qt.WidgetAttribute.WA_TranslucentBackground)
        self.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose)
        self.setFixedSize(W, H)
        self.setWindowTitle("JUICE")
        dpr = self.devicePixelRatioF() or 1.0
        self._logo = QPixmap()
        if os.path.isfile(theme.LOGO):
            self._logo = QPixmap(theme.LOGO).scaled(int(LOGO * dpr), int(LOGO * dpr), Qt.AspectRatioMode.KeepAspectRatio,
                                                    Qt.TransformationMode.SmoothTransformation)
            self._logo.setDevicePixelRatio(dpr)
        scr = QGuiApplication.screenAt(self.cursor().pos()) or QGuiApplication.primaryScreen()
        if scr is not None:
            g = scr.availableGeometry()
            self.move(g.center().x() - W // 2, g.center().y() - H // 2)
        self._timer = QTimer(self); self._timer.timeout.connect(self._tick); self._timer.start(30)

    def _tick(self):
        self._phase = (self._phase + 0.012) % 1.0
        waited = time.monotonic() - self._t0
        there = bool(self.marker) and not os.path.exists(self.marker)
        if there or waited > TIMEOUT_S or (waited > 1.0 and not _alive(self.pid)):
            self.close(); return
        self.update()

    def closeEvent(self, ev):
        self._timer.stop()
        try:
            if self.marker and os.path.exists(self.marker):
                os.remove(self.marker)
        except OSError:
            pass
        if self in _open:
            _open.remove(self)
        super().closeEvent(ev)

    def paintEvent(self, _):
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing)
        card = QRectF(1, 1, W - 2, H - 2)
        bg = QLinearGradient(0, 0, W, H); bg.setColorAt(0.0, QColor("#123a5a")); bg.setColorAt(1.0, NAVY_DEEP)
        p.setBrush(QBrush(bg)); p.setPen(QPen(QColor(143, 193, 230, 120), 1.2))
        p.drawRoundedRect(card, 14, 14)

        cx, cy = 26 + LOGO / 2, H / 2 - 8
        if not self._logo.isNull():
            clip = QPainterPath(); clip.addEllipse(QPointF(cx, cy), LOGO / 2 - 1, LOGO / 2 - 1)
            p.save(); p.setClipPath(clip); p.drawPixmap(QPointF(cx - LOGO / 2, cy - LOGO / 2), self._logo); p.restore()
        p.setBrush(Qt.BrushStyle.NoBrush); p.setPen(QPen(WHITE, 2)); p.drawEllipse(QPointF(cx, cy), LOGO / 2, LOGO / 2)

        x = 26 + LOGO + 22
        f = QFont("Segoe UI", 15); f.setBold(True); p.setFont(f); p.setPen(WHITE)
        p.drawText(QRectF(x, 30, W - x - 20, 28), Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter, "JUICE")
        p.setFont(QFont("Segoe UI", 10)); p.setPen(DIM)
        p.drawText(QRectF(x, 60, W - x - 20, 22), Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter, f"Opening {self.what} …")

        track = QRectF(26, H - 30, W - 52, 5)                      # a bar that keeps moving: no known end
        p.setPen(Qt.PenStyle.NoPen); p.setBrush(QColor(143, 193, 230, 50)); p.drawRoundedRect(track, 2.5, 2.5)
        seg = track.width() * 0.28
        left = track.left() - seg + (track.width() + seg) * self._phase
        p.save(); p.setClipRect(track)
        p.setBrush(ORANGE); p.drawRoundedRect(QRectF(left, track.top(), seg, track.height()), 2.5, 2.5)
        p.restore()


def start(module: str, args: list, cwd: str, what: str) -> bool:
    """Start `python -m module args` as an independent process; with a running application, show the card until its
    window is on screen."""
    from PyQt6.QtWidgets import QApplication
    marker = os.path.join(tempfile.gettempdir(), f"juice_opening_{uuid.uuid4().hex}.tmp")
    try:
        open(marker, "w").close()
    except OSError:
        marker = ""
    proc = QProcess()
    proc.setProgram(sys.executable); proc.setArguments(["-m", module] + list(args)); proc.setWorkingDirectory(cwd)
    env = QProcessEnvironment.systemEnvironment()
    if marker:
        env.insert(READY_ENV, marker)
    else:
        env.remove(READY_ENV)
    proc.setProcessEnvironment(env)
    ok, pid = proc.startDetached()
    if not ok:
        if marker:
            os.remove(marker)
        return False
    if QApplication.instance() is not None and marker:
        card = LoadingWindow(what, marker, int(pid or 0))
        _open.append(card); card.show()
    return True
