"""'Compare methods' of the partitioning panel: where the experiment's species sit in volatility and oxidation state,
by every volatility method (one cloud each, the chosen one bold) and, for one method, by structural group.
The numbers come from juice/chemistry/volatility_map.py; this window only shows them."""
from __future__ import annotations

import os
from typing import Callable, Optional

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from PyQt6.QtCore import QThread, Qt, pyqtSignal
from PyQt6.QtWidgets import (QApplication, QComboBox, QDialog, QFileDialog, QHBoxLayout, QLabel, QProgressBar,
                             QPushButton, QTabWidget, QVBoxLayout, QWidget)

from juice.chemistry import volatility as vol
from juice.chemistry import volatility_map as vmap

WEIGHTS = [("√(peak concentration)", 0.5), ("each species once", 0.0), ("peak concentration", 1.0)]


class _Builder(QThread):
    """volatility_map.build in the background; step(text, fraction) while it runs, done(map or None, error) at the end."""
    step = pyqtSignal(str, float)
    done = pyqtSignal(object, str)

    def __init__(self, spec, run_file: str):
        super().__init__()
        self.spec, self.run_file = spec, run_file

    def run(self):
        try:
            vm = vmap.build(self.spec, self.run_file, progress=lambda t, f: self.step.emit(t, f))
            self.done.emit(vm, "")
        except Exception as e:
            self.done.emit(None, f"{type(e).__name__}: {e}")


class VolatilityWindow(QDialog):
    """`spec`: the experiment as set up so far; `chosen`: the method id of the panel; `run_file`: a .cs of this
    experiment for the weights ('' = none yet)."""

    def __init__(self, spec, chosen: str, run_file: str = "", parent=None,
                 on_choose: Optional[Callable[[str], None]] = None):
        super().__init__(parent)
        self.setWindowTitle(f"Volatility methods compared: {spec.name}")
        self.resize(1150, 820)
        self.spec, self.chosen, self.on_choose = spec, chosen, on_choose
        lay = QVBoxLayout(self)
        top = QHBoxLayout()
        top.addWidget(QLabel("Weight"))
        self.weight = QComboBox()
        for label, p in WEIGHTS:
            self.weight.addItem(label, p)
        self.weight.setToolTip("How much each species counts in a cloud, from its peak concentration in the run:\n"
                               "the square root keeps minor products visible without letting them dominate.")
        top.addWidget(self.weight)
        self.run_label = QLabel(); self.run_label.setStyleSheet("color: #666;")
        top.addWidget(self.run_label, stretch=1)
        b = QPushButton("Run results…"); b.setToolTip("A .cs of this experiment, for the weights"); b.clicked.connect(self._browse)
        top.addWidget(b)
        lay.addLayout(top)
        self.info = QLabel(); self.info.setWordWrap(True); self.info.setStyleSheet("color: #666;")
        lay.addWidget(self.info)

        self.tabs = QTabWidget(); lay.addWidget(self.tabs, stretch=1)
        self.fig_m, self.can_m = self._tab("All methods")
        w = self.tabs.widget(0)
        self.fig_g, self.can_g = self._tab("By group")
        grow = QHBoxLayout(); grow.addWidget(QLabel("Method"))
        self.group_method = QComboBox()
        for m in vol.METHODS.values():
            self.group_method.addItem(m.label, m.id)
        grow.addWidget(self.group_method, stretch=1)
        self.use_btn = QPushButton("Use this method"); self.use_btn.clicked.connect(self._use)
        self.use_btn.setToolTip("Set the partitioning panel to this method")
        grow.addWidget(self.use_btn)
        self.tabs.widget(1).layout().insertLayout(0, grow)

        # the species and their C* are worked out in a thread: the window shows what it is doing meanwhile
        self.busy = QWidget(); bl = QVBoxLayout(self.busy)
        self.busy_label = QLabel("Starting…"); self.bar = QProgressBar(); self.bar.setRange(0, 100)
        bl.addStretch(); bl.addWidget(self.busy_label); bl.addWidget(self.bar); bl.addStretch()
        lay.insertWidget(lay.indexOf(self.tabs), self.busy)
        self.tabs.setVisible(False); self.weight.setEnabled(False); self.group_method.setEnabled(False)
        self.vm, self.error = None, ""
        i = self.group_method.findData(chosen)
        self.group_method.setCurrentIndex(i if i >= 0 else 0)
        self._worker = _Builder(spec, run_file if run_file and os.path.isfile(run_file) else "")
        self._worker.step.connect(self._step)
        self._worker.done.connect(self._built)
        self._worker.start()

    def _step(self, text: str, frac: float):
        self.busy_label.setText(text + "…")
        self.bar.setValue(int(100 * frac))

    def _built(self, vm, error: str):
        self.vm, self.error = vm, error
        self._step("Drawing the clouds", 0.95)
        QApplication.processEvents()
        self.busy.setVisible(False); self.tabs.setVisible(True); self.group_method.setEnabled(True)
        if self.vm is not None and self.vm.peak is None:
            self.weight.setCurrentIndex(1)
        self.weight.currentIndexChanged.connect(lambda *_: self.redraw())
        self.group_method.currentIndexChanged.connect(lambda *_: self.redraw(methods=False))
        self.redraw()

    def closeEvent(self, e):
        if self._worker.isRunning():
            self._worker.wait(30000)                  # it cannot be interrupted; it ends within seconds
        super().closeEvent(e)

    def _tab(self, title):
        w = QWidget(); v = QVBoxLayout(w); v.setContentsMargins(0, 0, 0, 0)
        fig = Figure(figsize=(11, 7)); can = FigureCanvas(fig)
        v.addWidget(NavigationToolbar(can, self)); v.addWidget(can, stretch=1)
        self.tabs.addTab(w, title)
        return fig, can

    def _browse(self):
        p, _ = QFileDialog.getOpenFileName(self, "Results of this experiment", "", "JUICE results (*.cs);;All files (*)")
        if p and self.vm is not None:
            try:
                vmap.load_run(self.vm, p)
            except Exception as e:
                self.info.setText(f"Not read: {e}"); return
            self.weight.setCurrentIndex(0)
            self.redraw()

    def _use(self):
        mid = self.group_method.currentData()
        if self.on_choose and mid in (self.vm.log_cstar if self.vm else {}):
            self.on_choose(mid)
            self.chosen = mid
            self.redraw()

    def redraw(self, methods: bool = True):
        if self.vm is None:
            self.info.setText(f"Not drawn: {self.error}")
            return
        p = float(self.weight.currentData()) if self.vm.peak is not None else 0.0
        self.weight.setEnabled(self.vm.peak is not None)
        self.run_label.setText(f"run: {os.path.basename(self.vm.run_file)}" if self.vm.run_file else
                               "no run of this experiment yet: each species counted once")
        miss = ", ".join(vol.METHODS[m].label for m in self.vm.unavailable)
        self.info.setText(f"{len(self.vm.names)} closed-shell organic species of the experiment's mechanism, at 298 K. "
                          "Each cloud is the smoothed region holding 80 % of them; the shaded bands are the volatility "
                          "classes" + (", the dashed line the run's largest organic aerosol mass" if self.vm.coa_max else "")
                          + "." + (f" Not installed: {miss}." if miss else ""))
        QApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
        try:
            if methods:
                vmap.draw_methods(self.fig_m, self.vm, self.chosen, p); self.can_m.draw()
            mid = self.group_method.currentData()
            vmap.draw_groups(self.fig_g, self.vm, mid, p); self.can_g.draw()
            self.use_btn.setEnabled(mid != self.chosen and mid in self.vm.log_cstar)
        finally:
            QApplication.restoreOverrideCursor()
