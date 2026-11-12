"""Gas <-> particle <-> wall partitioning on the Chemistry step: the particles of the run (none, one size, or a
size-distribution file), the scheme as a drawing, the rate constants the run uses, the C* (volatility) method and,
for a GECKO-A base, the condensable species with their C* and C_w.

Only a GECKO-A base has the scheme today (its <id>.SOA.include, written by make_soa_condensation.emit_soa), so
the panel is live for GECKO-A and greyed out, marked pending, for every other base. The values shown are those of
the run: read from the include it loads, with the condensation sink of the chosen particles and the wall terms
scaled to the chamber (experiment_spec.partitioning_lines writes the same into the EASY file).

Per condensable species X (Krechmer, Day & Jimenez 2020; Krechmer et al. 2016 for C_w):
  X -> SOA_X        k_cond = CS                      SOA_X -> X     k_evap = CS * C*/C_OA
  X -> WALL_X       k_gw                             WALL_X -> X    k_wg = k_gw * C*/C_w,  C_w = s * 5.33 * max(1, C*^0.6)
  SOA_X -> WALL_PL_X  k_pw (particles lost to the wall, irreversible)      SOA_X + DIL -> (dilution)
"""
from __future__ import annotations

import math
import os
from typing import Dict

import numpy as np
from PyQt6.QtCore import QPointF, QRectF, Qt
from PyQt6.QtGui import QColor, QFont, QPainter, QPalette, QPen, QPolygonF, QStandardItemModel
from PyQt6.QtWidgets import (QAbstractItemView, QButtonGroup, QCheckBox, QComboBox, QDoubleSpinBox, QFileDialog,
                             QGridLayout, QGroupBox, QHBoxLayout, QHeaderView, QLabel, QLineEdit, QPushButton,
                             QRadioButton, QSizePolicy, QTableWidget, QTableWidgetItem, QVBoxLayout, QWidget)

from juice.core.particles import read_soa_include, soa_defaults as _defaults

# volatility methods: (label, available). Only GECKO-A's own vapour pressures are wired to a run.
CSTAR_METHODS = [
    ("Nannoolal et al. (2008), as GECKO-A computes it (pvap.nan.dat)", True),
    ("SIMPOL.1 (Pankow and Asher, 2008), group contributions from the SMILES (pending)", False),
    ("EVAPORATION (Compernolle et al., 2011) (pending)", False),
    ("Measured or literature C* from a table (pending)", False),
]
PARTICLE_MODES = ("none", "mono", "file")


def _tau(k: float) -> str:
    if not k:
        return "never"
    s = 1.0 / k
    return f"τ {s:.0f} s" if s < 120 else (f"τ {s / 60:.0f} min" if s < 7200 else f"τ {s / 3600:.1f} h")


class PartitionDiagram(QWidget):
    """Gas in the middle, particles to the left, the wall to the right; particles lost to the wall below."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.values: Dict = _defaults()
        self.setMinimumHeight(250)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)

    def set_values(self, values: Dict):
        self.values = values
        self.update()

    def paintEvent(self, _e):
        p = QPainter(self); p.setRenderHint(QPainter.RenderHint.Antialiasing)
        pal = self.palette()
        grp = QPalette.ColorGroup.Active if self.isEnabled() else QPalette.ColorGroup.Disabled
        fg = pal.color(grp, QPalette.ColorRole.WindowText)
        base = pal.color(grp, QPalette.ColorRole.Base)
        accent = QColor("#6a4c93") if self.isEnabled() else fg
        v = self.values
        w = self.width()
        bw, bh = min(170, w * 0.2), 64
        y = 20
        boxes = {"particle": QRectF(w * 0.03, y, bw, bh), "gas": QRectF((w - bw) / 2, y, bw, bh),
                 "wall": QRectF(w * 0.97 - bw, y, bw, bh)}
        lost = QRectF(boxes["particle"].left(), y + bh + 110, bw, 44)
        dil = QRectF((w - bw) / 2, y + bh + 110, bw, 44)
        titles = {"particle": ("Particles", "SOA_X"), "gas": ("Gas", "X"), "wall": ("Chamber wall", "WALL_X")}
        f = QFont(self.font()); fb = QFont(f); fb.setBold(True); fs = QFont(f); fs.setPointSizeF(max(7.0, f.pointSizeF() - 1))
        for key, r in boxes.items():
            p.setPen(QPen(accent, 1.6)); p.setBrush(base); p.drawRoundedRect(r, 8, 8)
            p.setPen(fg); p.setFont(fb); p.drawText(r.adjusted(0, 8, 0, -bh / 2), Qt.AlignmentFlag.AlignHCenter, titles[key][0])
            p.setFont(f); p.drawText(r.adjusted(0, bh / 2, 0, 0), Qt.AlignmentFlag.AlignHCenter, titles[key][1])
        for r, t1, t2 in ((lost, "Lost to the wall", "WALL_PL_X"), (dil, "Diluted out", "")):
            p.setPen(QPen(fg, 1, Qt.PenStyle.DashLine)); p.setBrush(Qt.BrushStyle.NoBrush); p.drawRoundedRect(r, 6, 6)
            p.setPen(fg); p.setFont(f); p.drawText(r.adjusted(0, 4, 0, 0), Qt.AlignmentFlag.AlignHCenter, t1)
            if t2:
                p.setFont(fs); p.drawText(r.adjusted(0, 22, 0, 0), Qt.AlignmentFlag.AlignHCenter, t2)

        def arrow(a: QPointF, b: QPointF, label: str, above: bool, dashed: bool = False):
            p.setPen(QPen(accent, 1.4, Qt.PenStyle.DashLine if dashed else Qt.PenStyle.SolidLine))
            p.drawLine(a, b)
            ang = math.atan2(b.y() - a.y(), b.x() - a.x())
            head = QPolygonF([b, QPointF(b.x() - 9 * math.cos(ang - 0.4), b.y() - 9 * math.sin(ang - 0.4)),
                              QPointF(b.x() - 9 * math.cos(ang + 0.4), b.y() - 9 * math.sin(ang + 0.4))])
            p.setBrush(accent); p.setPen(Qt.PenStyle.NoPen); p.drawPolygon(head); p.setBrush(Qt.BrushStyle.NoBrush)
            p.setPen(fg); p.setFont(fs)
            mid = QPointF((a.x() + b.x()) / 2, (a.y() + b.y()) / 2)
            lines = label.split("\n")
            tw = max(p.fontMetrics().horizontalAdvance(t) for t in lines) + 6
            th = p.fontMetrics().height() * len(lines)
            if abs(b.x() - a.x()) > abs(b.y() - a.y()):       # horizontal: text above or below the arrow
                rect = QRectF(mid.x() - tw / 2, mid.y() - th - 3 if above else mid.y() + 3, tw, th)
            else:                                            # vertical: text beside it
                rect = QRectF(mid.x() + 6, mid.y() - th / 2, tw, th)
            p.drawText(rect, Qt.AlignmentFlag.AlignCenter, label)

        g, pa, wa = boxes["gas"], boxes["particle"], boxes["wall"]
        cs, kgw, kpw = v.get("KCOND", 0), v.get("KFWALL", 0), v.get("KPWALL", 0)
        arrow(QPointF(g.left(), g.top() + 20), QPointF(pa.right(), pa.top() + 20),
              v.get("CS_TEXT") or f"k_cond = CS = {cs:.2g} s⁻¹  ({_tau(cs)})", True)
        arrow(QPointF(pa.right(), pa.bottom() - 20), QPointF(g.left(), g.bottom() - 20),
              f"k_evap = CS · C* / C_OA\n(C_OA from {v.get('COA0', 0.5):g} µg m⁻³ + the SOA formed)", False)
        arrow(QPointF(g.right(), g.top() + 20), QPointF(wa.left(), wa.top() + 20),
              f"k_gw = {kgw:.2g} s⁻¹  ({_tau(kgw)})", True)
        arrow(QPointF(wa.left(), wa.bottom() - 20), QPointF(g.right(), g.bottom() - 20),
              f"k_wg = k_gw · C* / C_w\nC_w = {v.get('CW_SCALE', 1):.3g} · {v.get('CW_PRE', 5.33):.3g} · max(1, C*^{v.get('CW_N', 0.6):g}) µg m⁻³", False)
        arrow(QPointF(pa.center().x(), pa.bottom()), QPointF(lost.center().x(), lost.top()),
              f"k_pw = {kpw:.2g} s⁻¹\n({_tau(kpw)}, irreversible)", True)
        arrow(QPointF(g.center().x(), g.bottom()), QPointF(dil.center().x(), dil.top()), "k_dil = DILUTE\n(chamber flow)", True, dashed=True)
        p.end()


class PartitioningBox(QGroupBox):
    """The panel. `enable` is the run switch (ChemistryPage.gecko_part, ExperimentSpec.gecko_partitioning);
    load(spec) / store(spec) carry the particle choice; `context` gives (chamber profile, T in K, date)."""

    def __init__(self, parent=None):
        super().__init__("Gas–particle–wall partitioning", parent)
        self._b = None
        self.context = lambda: ("", 298.15, "")
        lay = QVBoxLayout(self)
        top = QHBoxLayout()
        self.enable = QCheckBox("Partition the condensable species between gas, particles and the chamber wall")
        self.enable.setChecked(True); top.addWidget(self.enable, stretch=1)
        self.toggle = QPushButton("ⓘ Scheme and values"); self.toggle.setCheckable(True); self.toggle.setFlat(True)
        self.toggle.setStyleSheet("QPushButton { color: #6a4c93; border: none; padding: 0 4px; } QPushButton:checked { font-weight: bold; }")
        self.toggle.setMinimumWidth(self.toggle.fontMetrics().horizontalAdvance(self.toggle.text()) + 24)
        top.addWidget(self.toggle); lay.addLayout(top)
        self.status = QLabel(); self.status.setWordWrap(True); self.status.setStyleSheet("color: #888;"); lay.addWidget(self.status)

        # the particles of the run: they set the condensation sink
        self.pbox = QWidget(); pg = QGridLayout(self.pbox); pg.setContentsMargins(0, 0, 0, 0)
        pg.addWidget(QLabel("Particles"), 0, 0)
        self.p_none = QRadioButton("None: the gas exchanges with the walls only")
        self.p_mono = QRadioButton("One size"); self.p_file = QRadioButton("Size distribution file (SMPS)")
        self.p_group = QButtonGroup(self)
        for i, b in enumerate((self.p_none, self.p_mono, self.p_file)):
            self.p_group.addButton(b, i)
        self.p_none.setChecked(True)
        pg.addWidget(self.p_none, 0, 1, 1, 2)
        pg.addWidget(self.p_mono, 1, 1)
        self.p_d = QDoubleSpinBox(); self.p_d.setRange(1, 10000); self.p_d.setDecimals(0); self.p_d.setValue(100); self.p_d.setSuffix(" nm")
        self.p_n = QLineEdit("1e4"); self.p_n.setMaximumWidth(110); self.p_n.setToolTip("number concentration in cm⁻³, e.g. 1e4")
        mh = QHBoxLayout(); mh.addWidget(self.p_d); mh.addWidget(self.p_n); mh.addWidget(QLabel("cm⁻³")); mh.addStretch()
        pg.addLayout(mh, 1, 2)
        pg.addWidget(self.p_file, 2, 1)
        self.p_path = QLineEdit(); self.p_path.setPlaceholderText("time column, then one column per size bin with its diameter (nm) as header")
        self.p_path.setToolTip("CSV, tab, semicolon or Excel; lines starting with # are skipped.\n"
                               "First column: the time (YYYY-MM-DD HH:MM, HH:MM, hours of the day, seconds since 2000 or unix seconds).\n"
                               "Other columns: one size bin each, the header its mid diameter in nm ('100', '100.5 nm', 'Dp_100.5';\n"
                               "read as µm when all are below 5). Values: dN/dlogDp or the number per bin, in cm⁻³.")
        fh = QHBoxLayout(); fh.addWidget(self.p_path, stretch=1)
        b = QPushButton("Browse…"); b.clicked.connect(self._browse); fh.addWidget(b)
        self.p_values = QComboBox(); self.p_values.addItem("dN/dlogDp", "dNdlogDp"); self.p_values.addItem("dN per bin", "dN")
        fh.addWidget(self.p_values); pg.addLayout(fh, 2, 2)
        self.p_info = QLabel(); self.p_info.setWordWrap(True); self.p_info.setStyleSheet("color: #888;")
        pg.addWidget(self.p_info, 3, 1, 1, 2); pg.setColumnStretch(2, 1)
        lay.addWidget(self.pbox)
        self.p_group.idToggled.connect(lambda *_: self.refresh())
        self.p_d.valueChanged.connect(lambda *_: self.refresh())
        self.p_values.currentIndexChanged.connect(lambda *_: self.refresh())
        self.p_n.editingFinished.connect(self.refresh); self.p_path.editingFinished.connect(self.refresh)

        self.details = QWidget(); dl = QVBoxLayout(self.details); dl.setContentsMargins(0, 0, 0, 0)
        self.diagram = PartitionDiagram(); dl.addWidget(self.diagram)
        mrow = QHBoxLayout(); mrow.addWidget(QLabel("C* (volatility) from"))
        self.method = QComboBox(); model: QStandardItemModel = self.method.model()
        for label, ok in CSTAR_METHODS:
            self.method.addItem(label)
            model.item(self.method.count() - 1).setEnabled(ok)
        self.method.setToolTip("C* = p_sat · M / (R T) in µg m⁻³. Only the GECKO-A vapour pressures are wired to a run so far;\n"
                               "the other methods are listed for the discussion and not selectable yet.")
        mrow.addWidget(self.method, stretch=1); dl.addLayout(mrow)
        self.params = QTableWidget(0, 4); self.params.setHorizontalHeaderLabels(["Symbol", "Meaning", "Value in the run", "Source"])
        self._table_style(self.params); dl.addWidget(self.params)
        self.species_label = QLabel(); self.species_label.setStyleSheet("color: #888;"); dl.addWidget(self.species_label)
        self.species = QTableWidget(0, 6)
        self.species.setHorizontalHeaderLabels(["Species", "log10 C*", "C* (µg m⁻³)", "C_w (µg m⁻³)",
                                                "in particles at C_OA start", "k_evap at C_OA start (s⁻¹)"])
        self._table_style(self.species); self.species.horizontalHeader().setStretchLastSection(False)
        self.species.setMaximumHeight(170); dl.addWidget(self.species)
        self.species.horizontalHeaderItem(4).setToolTip("Gas–particle equilibrium without the wall: 1 / (1 + C*/C_OA)")
        lay.addWidget(self.details)
        self.details.setVisible(False)
        self.toggle.toggled.connect(self.details.setVisible)
        self.refresh()

    @staticmethod
    def _table_style(t: QTableWidget):
        t.verticalHeader().setVisible(False); t.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        t.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
        t.horizontalHeader().setStretchLastSection(True)

    def showEvent(self, e):                           # the chamber or temperature may have changed on an earlier step
        super().showEvent(e)
        self.refresh()

    # ── the particle choice ──
    def _browse(self):
        p, _ = QFileDialog.getOpenFileName(self, "Size distribution", os.path.dirname(self.p_path.text()),
                                           "Tables (*.csv *.txt *.tsv *.dat *.xlsx *.xls);;All files (*)")
        if p:
            self.p_path.setText(p); self.refresh()

    def particle_mode(self) -> str:
        return PARTICLE_MODES[max(0, self.p_group.checkedId())]

    def load(self, spec):
        {"mono": self.p_mono, "file": self.p_file}.get(spec.particles, self.p_none).setChecked(True)
        self.p_d.setValue(float(spec.particle_d_nm)); self.p_n.setText(f"{float(spec.particle_n_cm3):g}")
        self.p_path.setText(spec.particle_file or "")
        self.p_values.setCurrentIndex(max(0, self.p_values.findData(spec.particle_file_values)))
        self.refresh()

    def store(self, spec):
        spec.particles = self.particle_mode()
        spec.particle_d_nm = float(self.p_d.value())
        try:
            spec.particle_n_cm3 = float(self.p_n.text())
        except ValueError:
            pass
        spec.particle_file = self.p_path.text().strip()
        spec.particle_file_values = self.p_values.currentData()

    def _particles(self, T: float, date: str):
        """(CS for the drawing and tables, the diagram's CS label or '', one line on what was read)."""
        from juice.core import particles as pt
        mode = self.particle_mode()
        if mode == "mono":
            try:
                n = float(self.p_n.text())
            except ValueError:
                return 0.0, "", "the number concentration is not a number"
            cs = pt.condensation_sink(self.p_d.value(), n, T)
            return cs, "", f"CS = {cs:.3g} s⁻¹ ({_tau(cs)}) for {self.p_d.value():g} nm at {n:g} cm⁻³ and {T:g} K (Fuchs–Sutugin)"
        if mode == "file":
            path = self.p_path.text().strip()
            if not path:
                return 0.0, "k_cond = CS from the size distribution", "choose the size-distribution file"
            try:
                sd = pt.read_size_distribution(path, date, self.p_values.currentData())
                cs = sd.cs(T); lo, hi = float(np.nanmin(cs)), float(np.nanmax(cs))
            except Exception as e:
                return 0.0, "k_cond = CS from the size distribution", f"not read: {e}"
            return (float(np.nanmedian(cs)), f"k_cond = CS from the size distribution\n({lo:.2g}–{hi:.2g} s⁻¹)",
                    f"{sd.summary()}; CS {lo:.3g}–{hi:.3g} s⁻¹ at {T:g} K (ENZ column CS)" + (f"; {'; '.join(sd.notes)}" if sd.notes else ""))
        return 0.0, "k_cond = CS = 0\n(no particles)", "CS = 0: nothing condenses; the gas exchanges with the walls only"

    # ── everything shown ──
    def set_base(self, b):
        """Show the scheme of base mechanism `b` (a mechanism_registry.BaseMechanism, or None)."""
        self._b = b
        self.refresh()

    def refresh(self):
        from juice.core import chambers, particles as pt
        b = self._b
        live = bool(b is not None and getattr(b, "is_gecko", False) and b.path("soa"))
        v = read_soa_include(b.path("soa") if live else "")
        src = os.path.basename(v["path"]) if live else "make_soa_condensation.py defaults (proposal)"
        try:
            chamber, T, date = self.context()
        except Exception:
            chamber, T, date = "", 298.15, ""
        T = T or 298.15
        cs, cs_label, p_line = self._particles(T, date)
        v["KCOND"], v["CS_TEXT"] = cs, cs_label
        self.p_info.setText(p_line)
        try:
            prof = chambers.chamber(chamber) if chamber else None
        except Exception:
            prof = None
        scale = pt.wall_scale(prof.surface_m2, prof.volume_m3) if prof is not None else None
        if scale is not None:
            v["KFWALL"] *= scale; v["CW_SCALE"] *= scale
            walls = (f"{chamber}, A/V {prof.surface_m2 / prof.volume_m3:.3g} m⁻¹ = {scale:.3g} × SAPHIR's"
                     if abs(scale - 1) > 1e-9 else f"{chamber}, the chamber the values were set for")
        else:
            walls = (f"{chamber}: wall area not known (chamber profile), so SAPHIR's" if chamber
                     else "SAPHIR's (no chamber chosen yet)")
        mode = self.particle_mode()

        self.setTitle("Gas–particle–wall partitioning" + ("" if live else " (pending)"))
        self.enable.setVisible(live)                  # its tick is ExperimentSpec.gecko_partitioning: meaningless elsewhere
        self.species.setVisible(live)
        for w in (self.diagram, self.method, self.params, self.species, self.pbox):
            w.setEnabled(live)
        self.diagram.set_values(v)
        name = getattr(b, "id", "this mechanism")
        self.status.setText(
            f"Active for {name}: {len(v['species'])} condensable species, values from {src}; walls: {walls}." if live else
            f"Pending for {name}: only GECKO-A mechanisms have the scheme so far (it needs C* for every product). "
            "Shown greyed out with the values a GECKO-A run uses, as the starting point for deciding how other "
            "mechanisms get it.")
        rows = [("CS = k_cond", "condensation sink onto the particles",
                 "from the size distribution, at every time" if mode == "file" else f"{cs:.3g} s⁻¹ ({_tau(cs)})",
                 {"none": "no particles", "mono": "one size and number, Fuchs–Sutugin",
                  "file": "size-distribution file, Fuchs–Sutugin"}[mode]),
                ("k_gw", "gas → wall deposition", f"{v['KFWALL']:.3g} s⁻¹ ({_tau(v['KFWALL'])})",
                 f"1e-3 s⁻¹ for SAPHIR (Krechmer et al., 2020), scaled with A/V: {walls}"),
                ("k_wg", "wall → gas evaporation", "k_gw · C*/C_w, per species", "Krechmer et al. (2020)"),
                ("C_w", "wall absorbing capacity", f"{v['CW_SCALE']:.3g} · {v['CW_PRE']:.3g} · max(1, C*^{v['CW_N']:g}) µg m⁻³",
                 "Krechmer et al. (2016) for their 8 m³ chamber; SAPHIR's A/V, then scaled with A/V"),
                ("k_evap", "particle → gas evaporation", "CS · C*/C_OA, per species", "absorptive partitioning"),
                ("C_OA", "organic aerosol mass", f"{v['COA0']:g} µg m⁻³ at the start + the SOA formed", "seed organic floor"),
                ("k_pw", "particle loss to the wall (irreversible)", f"{v['KPWALL']:.3g} s⁻¹ ({_tau(v['KPWALL'])})",
                 "not yet per chamber; SAPHIR seed sulfate gives about 1.4e-5 s⁻¹ (τ ≈ 20 h)"),
                ("k_dil", "dilution by the replenishment flow", "DILUTE of the chamber profile", "chamber profile")]
        self.params.setRowCount(len(rows))
        for r, row in enumerate(rows):
            for c, txt in enumerate(row):
                self.params.setItem(r, c, QTableWidgetItem(txt))
        self.params.setFixedHeight(self.params.horizontalHeader().height() + sum(self.params.rowHeight(i) for i in range(len(rows))) + 4)
        sp = sorted(v["species"], key=lambda s: s[1])
        self.species_label.setText(f"Condensable species of {name}, lowest C* first ({len(sp)})" if live else
                                   "Condensable species: none (no partitioning include for this mechanism)")
        self.species.setRowCount(len(sp))
        coa = v["COA0"] or 0.5
        for r, (n, cstar, cw) in enumerate(sp):
            fp = 1.0 / (1.0 + cstar / coa)
            cw_run = cw * (scale or 1.0) if cw == cw else cw
            for c, txt in enumerate((n, f"{math.log10(cstar):.2f}" if cstar > 0 else "–", f"{cstar:.3g}",
                                     f"{cw_run:.3g}" if cw_run == cw_run else "–", f"{100 * fp:.3g} %", f"{cs * cstar / coa:.3g}")):
                it = QTableWidgetItem(txt)
                if c:
                    it.setTextAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
                self.species.setItem(r, c, it)
