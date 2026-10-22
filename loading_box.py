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
from PyQt6.QtWidgets import QApplication, QWidget, QVBoxLayout, QHBoxLayout, QLabel, QProgressBar

ACCENT = "#6b3e85"


class LoadingBox(QWidget):
    def __init__(self, title: str):
        super().__init__(None, Qt.WindowType.SplashScreen | Qt.WindowType.WindowStaysOnTopHint | Qt.WindowType.FramelessWindowHint)
        self.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, False)
        self.setStyleSheet(f"QWidget {{ background: #ffffff; }} QLabel {{ color: #222; }} "
                           f"QProgressBar {{ border: 1px solid #d0d0d0; border-radius: 3px; background: #f3f3f3; height: 10px; text-align: center; font-size: 9px; }} "
                           f"QProgressBar::chunk {{ background: {ACCENT}; }}")
        self.setFixedWidth(460)
        lay = QVBoxLayout(self); lay.setContentsMargins(18, 14, 18, 14); lay.setSpacing(8)
        head = QHBoxLayout()
        logo = QLabel()
        pm = QPixmap(os.path.join(os.path.dirname(os.path.abspath(__file__)), "JUICE Logo.png"))
        if not pm.isNull():
            logo.setPixmap(pm.scaledToHeight(40, Qt.TransformationMode.SmoothTransformation))
        head.addWidget(logo)
        t = QLabel(f"<b style='color:{ACCENT}; font-size: 13pt;'>JUICE</b> &nbsp; {title}"); head.addWidget(t, stretch=1)
        lay.addLayout(head)
        self.msg = QLabel("starting …"); self.msg.setWordWrap(True); lay.addWidget(self.msg)
        self.bar = QProgressBar(); self.bar.setRange(0, 0); self.bar.setTextVisible(False); lay.addWidget(self.bar)
        self.why = QLabel(""); self.why.setStyleSheet("color: #777; font-size: 9pt;"); self.why.setWordWrap(True); lay.addWidget(self.why)
        self._t0 = time.time()
        self._frame = self.frameGeometry()
        scr = QApplication.primaryScreen()
        if scr is not None:
            c = scr.availableGeometry().center()
            self.move(c.x() - 230, c.y() - 60)

    @classmethod
    def start(cls, title: str) -> "LoadingBox":
        app = QApplication.instance() or QApplication(sys.argv)
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
        self.msg.setText(f"{message}  <span style='color:#999'>({elapsed:.0f} s)</span>")
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
