"""
report_dock.py — the figure list in the results viewer.

Two buttons sit in the corner of the tab bar: "Add to figure list" snapshots the figures of the
tab you are looking at (every panel of a grid, both budgets, the four Chemical Space panels, …)
into the list, "Figure list (n)" shows this dock.  In the dock you reorder, rename, caption and
remove entries, add a screen image of the whole current page, and "Create document…" asks for a
folder and writes the numbered SVG/PNG set, an index and a PowerPoint there (report_basket).
"""
from __future__ import annotations

import os
import tempfile
from typing import List, Tuple

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (QDockWidget, QWidget, QVBoxLayout, QHBoxLayout, QListWidget, QListWidgetItem, QPushButton,
                             QToolButton, QLabel, QFileDialog, QInputDialog, QMessageBox, QLineEdit, QDialog, QDialogButtonBox,
                             QFormLayout)

from juice.gui.report_basket import FigureList


class FigureListDock(QDockWidget):
    def __init__(self, viewer, figure_list: FigureList, parent=None):
        super().__init__("Figure list", parent)
        self.v = viewer; self.fl = figure_list
        self.setAllowedAreas(Qt.DockWidgetArea.RightDockWidgetArea | Qt.DockWidgetArea.LeftDockWidgetArea | Qt.DockWidgetArea.BottomDockWidgetArea)
        w = QWidget(); lay = QVBoxLayout(w); lay.setContentsMargins(6, 6, 6, 6); lay.setSpacing(4)
        self.info = QLabel(""); self.info.setWordWrap(True); self.info.setStyleSheet("color: #666; font-size: 10px;"); lay.addWidget(self.info)
        self.list = QListWidget(); self.list.setSelectionMode(QListWidget.SelectionMode.ExtendedSelection)
        self.list.itemDoubleClicked.connect(lambda _i: self._edit()); lay.addWidget(self.list, stretch=1)
        row = QHBoxLayout()
        for text, tip, slot in (("▲", "move up", lambda: self._move(-1)), ("▼", "move down", lambda: self._move(1)),
                                ("Edit…", "title and caption", self._edit), ("Remove", "remove the selected entries (files too)", self._remove),
                                ("Page image", "add a screen image of the whole current page (PNG only)", self._add_page_image)):
            b = QToolButton(); b.setText(text); b.setToolTip(tip); b.setAutoRaise(True); b.clicked.connect(slot); row.addWidget(b)
        row.addStretch(); lay.addLayout(row)
        row2 = QHBoxLayout()
        self.build_btn = QPushButton("Create document…"); self.build_btn.setObjectName("primary"); self.build_btn.setToolTip("choose a folder: numbered SVG + PNG figures, an index and a PowerPoint are written there")
        self.build_btn.clicked.connect(self._build); row2.addWidget(self.build_btn)
        clr = QToolButton(); clr.setText("Clear"); clr.setAutoRaise(True); clr.setToolTip("empty the list"); clr.clicked.connect(self._clear); row2.addWidget(clr)
        row2.addStretch(); lay.addLayout(row2)
        self.setWidget(w); self.setMinimumWidth(300)
        self.refresh()

    # ── list ──
    def refresh(self):
        self.list.clear()
        for n, it in enumerate(self.fl.items, 1):
            li = QListWidgetItem(f"{n:02d}  {it.title}" + ("" if it.svg else "  (image)"))
            li.setData(Qt.ItemDataRole.UserRole, it.key); li.setToolTip(f"{it.caption}\nfrom {it.source}, added {it.added}")
            self.list.addItem(li)
        self.info.setText(f"{len(self.fl.items)} figure(s), kept in {self.fl.store_dir}")
        self.build_btn.setEnabled(bool(self.fl.items))
        cb = getattr(self.v, "_figlist_btn", None)
        if cb is not None:
            cb.setText(f"Figure list ({len(self.fl.items)})")

    def _keys(self) -> List[str]:
        return [i.data(Qt.ItemDataRole.UserRole) for i in self.list.selectedItems()]

    def _move(self, d: int):
        keys = self._keys()
        for k in (keys if d < 0 else reversed(keys)):
            self.fl.move(k, d)
        self.refresh()
        for i in range(self.list.count()):
            if self.list.item(i).data(Qt.ItemDataRole.UserRole) in keys:
                self.list.item(i).setSelected(True)

    def _remove(self):
        for k in self._keys():
            self.fl.remove(k)
        self.refresh()

    def _clear(self):
        if self.fl.items and QMessageBox.question(self, "Figure list", "Remove every figure from the list?") == QMessageBox.StandardButton.Yes:
            self.fl.clear(); self.refresh()

    def _edit(self):
        keys = self._keys()
        if not keys:
            return
        it = self.fl.get(keys[0])
        if it is None:
            return
        dlg = QDialog(self); dlg.setWindowTitle("Figure"); form = QFormLayout(dlg)
        t = QLineEdit(it.title); c = QLineEdit(it.caption); c.setMinimumWidth(420)
        form.addRow("Title", t); form.addRow("Caption", c)
        bb = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        bb.accepted.connect(dlg.accept); bb.rejected.connect(dlg.reject); form.addRow(bb)
        if dlg.exec() == QDialog.DialogCode.Accepted:
            self.fl.set_caption(it.key, c.text().strip() or it.caption, t.text().strip() or it.title); self.refresh()

    # ── adding ──
    def add_figures(self, figs: List[Tuple[str, object]], source: str):
        n = 0
        for title, fig in figs:
            try:
                self.fl.add(fig, title, source=source); n += 1
            except Exception as e:
                print(f"[figure list] could not add {title}: {e}")
        self.refresh()
        return n

    def _add_page_image(self):
        tab = self.v.right_tabs.currentWidget(); name = self.v.right_tabs.tabText(self.v.right_tabs.currentIndex())
        if tab is None:
            return
        tmp = os.path.join(tempfile.gettempdir(), "juice_page.png")
        tab.grab().save(tmp)
        self.fl.add_image(tmp, f"{name} page", source=name); self.refresh()

    # ── building ──
    def _build(self):
        if not self.fl.items:
            return
        start = getattr(self.v, "_last_report_dir", "") or os.path.dirname(self.fl.store_dir)
        folder = QFileDialog.getExistingDirectory(self, "Folder for the figures and the presentation", start)
        if not folder:
            return
        title, ok = QInputDialog.getText(self, "Document title", "Title for the title slide and the index:", text=self.fl.run_name)
        if not ok:
            return
        try:
            res = self.fl.build(folder, title=title.strip() or self.fl.run_name)
        except Exception as e:
            QMessageBox.critical(self, "Create document", str(e)); return
        self.v._last_report_dir = folder
        msg = (f"{len(res['figures'])} files in {os.path.join(folder, 'figures')}\n"
               f"index: {os.path.basename(res['index'])}\n" + (f"presentation: {os.path.basename(res['pptx'])}" if res['pptx'] else res['warning']))
        box = QMessageBox(self); box.setWindowTitle("Document created"); box.setText(msg)
        open_btn = box.addButton("Open folder", QMessageBox.ButtonRole.ActionRole); box.addButton(QMessageBox.StandardButton.Ok)
        box.exec()
        if box.clickedButton() is open_btn:
            try:
                os.startfile(folder)
            except Exception:
                pass
