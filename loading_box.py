"""
loading_box.py — a small window that says what JUICE is doing while a heavy window is built.

Imports only PyQt so it can appear before pandas / matplotlib / the viewer are loaded:

    from loading_box import LoadingBox
    box = LoadingBox.start("Opening results")          # creates the QApplication if needed
    box.step("Reading saphir_20240709.cs (25 MB)", 10)  # message + percent (None = busy bar)
    ...
    box.done()

Every step repaints immediately (processEvents), so the text is visible although the work
runs in the main thread.
"""
from __future__ import annotations

import os
import sys
import time
from typing import Optional

from PyQt6.QtCore import Qt
from PyQt6.QtGui import QPixmap
from PyQt6.QtWidgets import QApplication, QWidget, QVBoxLayout, QHBoxLayout, QLabel, QProgressBar, QFrame

ACCENT = "#6b3e85"
_APP = None                      # the QApplication this module created, if it had to
BOX_W, BOX_H = 480, 150          # fixed: the box never grows or shrinks with its text


class LoadingBox(QWidget):
    def __init__(self, title: str):
        super().__init__(None, Qt.WindowType.SplashScreen | Qt.WindowType.WindowStaysOnTopHint | Qt.WindowType.FramelessWindowHint)
        self.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, False)
        self.setFixedSize(BOX_W, BOX_H)
        outer = QVBoxLayout(self); outer.setContentsMargins(0, 0, 0, 0)
        frame = QFrame(); frame.setObjectName("loadingFrame"); outer.addWidget(frame)
        frame.setStyleSheet(f"#loadingFrame {{ background: #ffffff; border: 1px solid #c8c8c8; border-radius: 4px; }} "
                            f"QLabel {{ color: #222; background: transparent; border: none; }} "
                            f"QProgressBar {{ border: 1px solid #d0d0d0; border-radius: 3px; background: #f3f3f3; text-align: center; font-size: 9px; }} "
                            f"QProgressBar::chunk {{ background: {ACCENT}; }}")
        lay = QVBoxLayout(frame); lay.setContentsMargins(18, 12, 18, 12); lay.setSpacing(6)
        head = QHBoxLayout()
        logo = QLabel()
        pm = QPixmap(os.path.join(os.path.dirname(os.path.abspath(__file__)), "JUICE Logo.png"))
        if not pm.isNull():
            logo.setPixmap(pm.scaledToHeight(40, Qt.TransformationMode.SmoothTransformation))
        head.addWidget(logo)
        t = QLabel(f"<b style='color:{ACCENT}; font-size: 13pt;'>JUICE</b> &nbsp; {title}"); head.addWidget(t, stretch=1)
        lay.addLayout(head)
        self.msg = QLabel("starting …"); self.msg.setFixedHeight(20); lay.addWidget(self.msg)
        self.bar = QProgressBar(); self.bar.setRange(0, 0); self.bar.setTextVisible(False); self.bar.setFixedHeight(10); lay.addWidget(self.bar)
        self.why = QLabel(""); self.why.setStyleSheet("color: #777; font-size: 9pt;"); self.why.setWordWrap(True)
        self.why.setFixedHeight(34); self.why.setAlignment(Qt.AlignmentFlag.AlignTop | Qt.AlignmentFlag.AlignLeft); lay.addWidget(self.why)
        lay.addStretch()
        self._t0 = time.time()
        scr = QApplication.primaryScreen()
        if scr is not None:
            c = scr.availableGeometry().center()
            self.move(c.x() - BOX_W // 2, c.y() - BOX_H // 2)

    @classmethod
    def start(cls, title: str) -> "LoadingBox":
        global _APP
        app = QApplication.instance()
        if app is None:
            app = _APP = QApplication(sys.argv)      # keep a reference: a collected QApplication takes the box with it
        try:
            from juice_theme import app_icon, taskbar_identity
            taskbar_identity(); app.setWindowIcon(app_icon())
        except Exception:
            pass
        box = cls(title); box.show(); app.processEvents()
        return box

    def step(self, message: str, percent: Optional[int] = None, why: str = ""):
        """Say what happens now (and why it takes a moment).  percent None = busy bar."""
        elapsed = time.time() - self._t0
        fm = self.msg.fontMetrics()
        text = fm.elidedText(message, Qt.TextElideMode.ElideRight, BOX_W - 36 - fm.horizontalAdvance("  (000 s)"))
        self.msg.setText(f"{text}  <span style='color:#999'>({elapsed:.0f} s)</span>")
        self.why.setText(why)
        if percent is None:
            self.bar.setRange(0, 0)
        else:
            self.bar.setRange(0, 100); self.bar.setValue(max(0, min(100, int(percent))))
        self.repaint()
        QApplication.processEvents()

    def done(self):
        self.step("ready", 100)
        self.close()
