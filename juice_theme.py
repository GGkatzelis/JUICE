"""
juice_theme.py — one look for every JUICE window (step window, results viewer).

    from juice_theme import apply_theme, app_icon
    app = QApplication(argv); apply_theme(app)     # Fusion style, white palette, QSS, font, icon

Forced light theme: the operating system's dark mode never reaches a widget, every text is
dark on white, the accent is the JUICE purple.  The taskbar identity (AppUserModelID) is
set by `taskbar_identity()` before the QApplication exists.
"""
from __future__ import annotations

import os
import sys

from PyQt6.QtCore import QObject, QEvent
from PyQt6.QtGui import QColor, QFont, QIcon, QPalette, QPixmap
from PyQt6.QtWidgets import QComboBox

HERE = os.path.dirname(os.path.abspath(__file__))
LOGO = os.path.join(HERE, "JUICE Logo.png")
ICON = os.path.join(HERE, "JUICE.ico")          # multi-size icon: Windows takes the taskbar image from here
ACCENT, ACCENT_DARK, ACCENT_SOFT = "#6b3e85", "#4a3e85", "#efe6f7"

QSS = f"""
QWidget {{ background: #ffffff; color: #222222; }}
QMainWindow, QScrollArea, QSplitter, QDockWidget {{ background: #ffffff; }}
QMenuBar, QMenu, QStatusBar, QToolBar {{ background: #ffffff; color: #222222; }}
QMenu::item:selected {{ background: {ACCENT_SOFT}; }}
QLabel {{ background: transparent; color: #222222; }}
QLineEdit, QPlainTextEdit, QTextEdit, QComboBox, QSpinBox, QDoubleSpinBox, QListWidget, QTableWidget, QTreeWidget, QTableView, QListView, QTreeView {{
    background: #ffffff; color: #222222; border: 1px solid #c8c8c8; border-radius: 3px; selection-background-color: {ACCENT_SOFT}; selection-color: #222222; }}
QLineEdit:disabled, QComboBox:disabled, QSpinBox:disabled, QDoubleSpinBox:disabled {{ color: #999999; background: #f4f4f4; }}
QComboBox QAbstractItemView {{ background: #ffffff; color: #222222; selection-background-color: {ACCENT_SOFT}; }}
QHeaderView::section {{ background: #f0f0f0; color: #222222; border: 1px solid #d8d8d8; padding: 3px; }}
QTableWidget, QTableView {{ gridline-color: #e0e0e0; alternate-background-color: #fafafa; }}
QTableWidget::item:selected, QListWidget::item:selected, QTreeWidget::item:selected, QTableView::item:selected {{ background: {ACCENT_SOFT}; color: #222222; }}
QPushButton {{ background: #f3f3f3; color: #222222; border: 1px solid #c8c8c8; border-radius: 4px; padding: 5px 12px; }}
QPushButton:hover {{ background: #e9e9e9; }}
QPushButton:checked {{ background: {ACCENT_SOFT}; border-color: {ACCENT}; }}
QPushButton:disabled {{ color: #a0a0a0; background: #f7f7f7; }}
QPushButton#primary {{ background: {ACCENT}; color: #ffffff; font-weight: bold; padding: 8px 18px; border: none; }}
QPushButton#primary:hover {{ background: {ACCENT_DARK}; }}
QPushButton#primary:disabled {{ background: #c4c4c4; color: #ffffff; }}
QCheckBox, QRadioButton {{ background: transparent; color: #222222; }}
QGroupBox {{ font-weight: bold; color: {ACCENT_DARK}; border: 1px solid #dddddd; border-radius: 6px; margin-top: 12px; background: #ffffff; }}
QGroupBox::title {{ subcontrol-origin: margin; left: 10px; padding: 0 4px; background: #ffffff; }}
QGroupBox::indicator {{ width: 13px; height: 13px; }}
QListWidget#steps {{ background: #ffffff; border: none; font-size: 14px; }}
QListWidget#steps::item {{ padding: 10px 8px; border-bottom: 1px solid #eeeeee; color: #222222; }}
QListWidget#steps::item:selected {{ background: {ACCENT}; color: #ffffff; border-radius: 4px; }}
QLabel#title {{ font-size: 22px; font-weight: bold; color: {ACCENT_DARK}; }}
QLabel#subtitle {{ color: #666666; }}
QTabWidget::pane {{ border: 1px solid #dddddd; background: #ffffff; }}
QTabBar::tab {{ background: #f0f0f0; color: #222222; padding: 6px 12px; border: 1px solid #dddddd; border-bottom: none; }}
QTabBar::tab:selected {{ background: #ffffff; color: {ACCENT_DARK}; font-weight: bold; }}
QProgressBar {{ background: #f0f0f0; color: #222222; border: 1px solid #c8c8c8; border-radius: 3px; text-align: center; }}
QProgressBar::chunk {{ background-color: #2ecc40; }}
QScrollBar:vertical {{ background: #f4f4f4; width: 12px; }} QScrollBar::handle:vertical {{ background: #cccccc; border-radius: 4px; min-height: 20px; }}
QScrollBar:horizontal {{ background: #f4f4f4; height: 12px; }} QScrollBar::handle:horizontal {{ background: #cccccc; border-radius: 4px; min-width: 20px; }}
QSplitter::handle {{ background: #e6e6e6; }}
QToolTip {{ background: #ffffff; color: #222222; border: 1px solid #c8c8c8; }}
QLineEdit[missing="true"], QListWidget[missing="true"], QComboBox[missing="true"] {{ background: #fff4cc; border: 1px solid #e0b84a; }}
QRadioButton[missing="true"], QLabel[missing="true"] {{ color: #9a6b00; font-weight: bold; }}
"""


def mark_missing(widget, missing: bool) -> None:
    """Soft yellow on a field that still needs input; plain again once filled."""
    if widget is None or widget.property("missing") == missing:
        return
    widget.setProperty("missing", missing)
    widget.style().unpolish(widget); widget.style().polish(widget); widget.update()


def taskbar_identity(name: str = "FZJ.JUICE") -> None:
    """Own taskbar entry and icon instead of python.exe's (call before QApplication)."""
    if sys.platform.startswith("win"):
        try:
            import ctypes
            ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(name)
        except Exception:
            pass


def app_icon() -> QIcon:
    ico = QIcon(ICON) if os.path.isfile(ICON) else QIcon()
    if ico.isNull() and os.path.isfile(LOGO):
        ico = QIcon(QPixmap(LOGO))
    return ico


class _ComboFit(QObject):
    """App-wide: every QComboBox is wide enough for its entries and its popup list is at least as
    wide as the longest entry — no more 'p…' in a narrow control row (Georgios 2026-10-22)."""
    MAX_CLOSED = 340            # px: wider entries only widen the popup, not the closed box

    def eventFilter(self, obj, ev):
        try:
            if isinstance(obj, QComboBox):
                t = ev.type()
                if t == QEvent.Type.Show and not obj.property("juice_fit"):
                    obj.setProperty("juice_fit", True)
                    if self._longest(obj) <= self.MAX_CLOSED and obj.minimumWidth() == 0 and obj.maximumWidth() >= 16777215:
                        obj.setSizeAdjustPolicy(QComboBox.SizeAdjustPolicy.AdjustToContents)
                elif t in (QEvent.Type.MouseButtonPress, QEvent.Type.KeyPress):
                    w = self._longest(obj) + 48
                    if obj.view() is not None and w > obj.view().minimumWidth():
                        obj.view().setMinimumWidth(w)
        except Exception:
            pass
        return False

    @staticmethod
    def _longest(combo) -> int:
        fm = combo.fontMetrics()
        return max((fm.horizontalAdvance(combo.itemText(i)) for i in range(min(combo.count(), 2000))), default=0)


def apply_theme(app) -> None:
    """Fusion style + white palette + QSS + font + icon on a QApplication."""
    app.setStyle("Fusion")
    if getattr(app, "_juice_combo_fit", None) is None:
        app._juice_combo_fit = _ComboFit(app); app.installEventFilter(app._juice_combo_fit)
    pal = QPalette()
    for role, col in ((QPalette.ColorRole.Window, "#ffffff"), (QPalette.ColorRole.Base, "#ffffff"), (QPalette.ColorRole.AlternateBase, "#fafafa"),
                      (QPalette.ColorRole.Button, "#f3f3f3"), (QPalette.ColorRole.ToolTipBase, "#ffffff"), (QPalette.ColorRole.Highlight, ACCENT_SOFT)):
        pal.setColor(role, QColor(col))
    for role in (QPalette.ColorRole.WindowText, QPalette.ColorRole.Text, QPalette.ColorRole.ButtonText, QPalette.ColorRole.ToolTipText,
                 QPalette.ColorRole.HighlightedText, QPalette.ColorRole.PlaceholderText):
        pal.setColor(role, QColor("#222222" if role != QPalette.ColorRole.PlaceholderText else "#888888"))
    app.setPalette(pal)
    app.setStyleSheet(QSS)
    app.setFont(QFont("Segoe UI", 10))
    app.setWindowIcon(app_icon())
