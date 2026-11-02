"""
chamber_editor.py — edit a chamber profile or add a new one (Qt dialog around chambers.py)

Opened from "Where and when" (Edit … / New …).  Every field of a ChamberProfile is a widget; Save
writes chambers/<id>.chamber.json through chambers.save_profile() and the caller refreshes its list.
SAPHIR is the reference profile of the FACSIMILE scenarios: it opens read-only, with "Save as new
chamber …" as the way to derive another chamber from it.

    dlg = ChamberDialog(profile, parent)                 # edit (SAPHIR: read-only)
    dlg = ChamberDialog(profile, parent, as_new=True)    # a new chamber, pre-filled with a copy of `profile`
    if dlg.exec(): dlg.saved_id                          # id of the profile written to chambers/
"""
from __future__ import annotations

import copy
from typing import Optional

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (QDialog, QDialogButtonBox, QVBoxLayout, QHBoxLayout, QGridLayout, QGroupBox, QLabel,
                             QLineEdit, QComboBox, QCheckBox, QPlainTextEdit, QTableWidget, QTableWidgetItem,
                             QHeaderView, QAbstractItemView, QPushButton, QMessageBox, QWidget, QScrollArea,
                             QInputDialog)

import chambers as C


def _info(text: str) -> QLabel:
    lab = QLabel(text); lab.setWordWrap(True); lab.setStyleSheet("color: #888;")
    return lab


def _float(edit: QLineEdit, default: float = 0.0) -> float:
    t = edit.text().strip().replace(",", ".")
    return float(t) if t else default


class ChamberDialog(QDialog):
    def __init__(self, profile: C.ChamberProfile, parent=None, as_new: bool = False):
        super().__init__(parent)
        self.p = copy.deepcopy(profile)
        self.as_new = as_new
        self.locked = (not as_new) and profile.id == C.DEFAULT_CHAMBER   # the SAPHIR profile stays as shipped
        self.saved_id: str = ""
        self.setWindowTitle("New chamber" if as_new else f"Chamber — {profile.name}")
        self.resize(820, 760)

        outer = QVBoxLayout(self)
        scroll = QScrollArea(); scroll.setWidgetResizable(True); scroll.setFrameShape(QScrollArea.Shape.NoFrame)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        body = QWidget(); lay = QVBoxLayout(body); scroll.setWidget(body); outer.addWidget(scroll, stretch=1)
        if self.locked:
            lay.addWidget(_info("SAPHIR is the reference profile of the FACSIMILE scenarios and is not edited here. "
                                "Use 'Save as new chamber …' to start another chamber from its values."))
        elif as_new:
            lay.addWidget(_info("Pre-filled with a copy of the chamber you had selected. Give it an id (the file name) "
                                "and a name, then set the volume, the light source and the wall losses."))

        # ── identity ──
        box = QGroupBox("Chamber"); g = QGridLayout(box); r = 0
        g.addWidget(QLabel("Id (file name)"), r, 0); self.id = QLineEdit("" if as_new else self.p.id); g.addWidget(self.id, r, 1)
        self.id.setPlaceholderText("letters, digits, _ - .   →  chambers/<id>.chamber.json")
        self.id.setEnabled(as_new)                     # renaming = a new file; keep it simple
        g.addWidget(QLabel("Name"), r, 2); self.name = QLineEdit("" if as_new else self.p.name); g.addWidget(self.name, r, 3); r += 1
        g.addWidget(QLabel("Institution"), r, 0); self.inst = QLineEdit(self.p.institution); g.addWidget(self.inst, r, 1)
        g.addWidget(QLabel("Location"), r, 2); self.loc = QLineEdit(self.p.location); g.addWidget(self.loc, r, 3); r += 1
        g.addWidget(QLabel("Note"), r, 0); self.note = QPlainTextEdit(self.p.note); self.note.setMaximumHeight(60); g.addWidget(self.note, r, 1, 1, 3); r += 1
        g.addWidget(QLabel("References"), r, 0); self.refs = QPlainTextEdit("\n".join(self.p.references)); self.refs.setMaximumHeight(60)
        self.refs.setPlaceholderText("one per line"); g.addWidget(self.refs, r, 1, 1, 3); r += 1
        self.verified = QCheckBox("values checked against the chamber's publications or with its operators"); self.verified.setChecked(bool(self.p.verified))
        g.addWidget(self.verified, r, 1, 1, 3)
        lay.addWidget(box)

        # ── geometry and dilution ──
        box = QGroupBox("Volume, walls and dilution"); g = QGridLayout(box); r = 0
        g.addWidget(QLabel("Volume [m³]"), r, 0); self.vol = QLineEdit(f"{self.p.volume_m3:g}"); g.addWidget(self.vol, r, 1)
        g.addWidget(QLabel("Wall surface [m²]"), r, 2); self.surf = QLineEdit(f"{self.p.surface_m2:g}" if self.p.surface_m2 else ""); self.surf.setPlaceholderText("optional"); g.addWidget(self.surf, r, 3); r += 1
        g.addWidget(QLabel("Wall material"), r, 0); self.mat = QLineEdit(self.p.material); self.mat.setPlaceholderText("FEP film, stainless steel, …"); g.addWidget(self.mat, r, 1)
        g.addWidget(QLabel("Leak / extra inflow [m³/h]"), r, 2); self.leak = QLineEdit(f"{self.p.leak_flow_m3h:g}"); g.addWidget(self.leak, r, 3); r += 1
        g.addWidget(QLabel("DILUTE expression"), r, 0); self.dil = QLineEdit(self.p.dilute_expr); g.addWidget(self.dil, r, 1, 1, 3); r += 1
        self.dil.setPlaceholderText("empty = built from volume and leak flow")
        self.dil_preview = _info(""); g.addWidget(self.dil_preview, r, 1, 1, 3); r += 1
        g.addWidget(_info("Dilution per second = (FLOW column + leak) / (volume · 3600). Only type an expression when the "
                          "chamber's dilution is not of that form (the FACSIMILE/EASY syntax of the mechanism applies)."), r, 0, 1, 4)
        for w in (self.vol, self.leak, self.dil):
            w.textChanged.connect(self._dilute_preview)
        lay.addWidget(box)

        # ── light ──
        box = QGroupBox("Light source"); g = QGridLayout(box); r = 0
        g.addWidget(QLabel("Type"), r, 0); self.ltype = QComboBox()
        self.ltype.addItems(["sun — outdoor chamber", "lamps — indoor chamber with UV lamps", "dark — no light source"])
        self.ltype.setSizeAdjustPolicy(QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
        self.ltype.setCurrentIndex(C.LIGHT_TYPES.index(self.p.light_type)); g.addWidget(self.ltype, r, 1, 1, 3); r += 1
        lat, lon = self.p.lat_lon
        self.lat_lab = QLabel("Latitude [°N]"); g.addWidget(self.lat_lab, r, 0); self.lat = QLineEdit(f"{lat:g}"); g.addWidget(self.lat, r, 1)
        self.lon_lab = QLabel("Longitude [°E]"); g.addWidget(self.lon_lab, r, 2); self.lon = QLineEdit(f"{lon:g}"); g.addWidget(self.lon, r, 3); r += 1
        self.lamps_lab = QLabel("Lamps"); g.addWidget(self.lamps_lab, r, 0); self.lamps = QLineEdit(str(self.p.light.get("lamps", "") or ""))
        self.lamps.setPlaceholderText("e.g. 40 × Philips TL-K 40W/05 blacklights (350 nm), Xe arc, …"); g.addWidget(self.lamps, r, 1, 1, 3); r += 1
        self.jno2_lab = QLabel("jNO2 at full power [1/s]"); g.addWidget(self.jno2_lab, r, 0)
        j = self.p.light.get("jno2_full_s1", 0) or 0
        self.jno2 = QLineEdit(f"{float(j):g}" if j else ""); self.jno2.setPlaceholderText("actinometry or lamp spectrum, e.g. 2e-3"); g.addWidget(self.jno2, r, 1); r += 1
        self.light_note = _info(""); g.addWidget(self.light_note, r, 0, 1, 4)
        self.ltype.currentIndexChanged.connect(self._light_changed)
        lay.addWidget(box)

        # ── archive ──
        box = QGroupBox("Data archive JUICE can read"); g = QGridLayout(box)
        g.addWidget(QLabel("Kind"), 0, 0); self.akind = QComboBox(); self.akind.addItems(["none — conditions from a file or typed", "saphir_nc — SAPHIR netCDF day folders"])
        self.akind.setSizeAdjustPolicy(QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
        self.akind.setCurrentIndex(1 if (self.p.archive or {}).get("kind") == "saphir_nc" else 0); g.addWidget(self.akind, 0, 1)
        g.addWidget(QLabel("Root folder"), 1, 0); self.aroot = QLineEdit(self.p.archive_root); g.addWidget(self.aroot, 1, 1)
        self.akind.currentIndexChanged.connect(lambda i: self.aroot.setEnabled(i > 0)); self.aroot.setEnabled(self.akind.currentIndex() > 0)
        lay.addWidget(box)

        # ── default effects ──
        box = QGroupBox("Chamber effects switched on by default when this chamber is chosen"); v = QVBoxLayout(box)
        titles = {"walldep": "Wall deposition of HNO3, HONO, H2O2, HO2, O3, N2O5 (SAPHIR rates)",
                  "hono": "Photolytic HONO wall source (SAPHIR parametrisation)",
                  "hcho": "Photolytic HCHO wall source (SAPHIR parametrisation)",
                  "bgreac": "Background OH reactivity (CO-equivalent Y)",
                  "dilution": "Dilution by the replenishment flow",
                  "dark_hono": "Dark HONO formation on the walls (SAPHIR parametrisation)",
                  "light_source": "NO2 source under illumination (SNO2)"}
        self.effects = {}
        for k in C.EFFECT_KEYS:
            cb = QCheckBox(titles[k]); cb.setChecked(self.p.effect(k)); self.effects[k] = cb; v.addWidget(cb)
        v.addWidget(_info("The SAPHIR wall sources are fitted to SAPHIR's FEP film; for another chamber leave them off unless you "
                          "have checked that the parametrisation applies. The experiment can still switch each one."))
        lay.addWidget(box)

        # ── wall losses ──
        box = QGroupBox("First-order wall loss per species [1/s]"); v = QVBoxLayout(box)
        h = QHBoxLayout(); h.addWidget(_info("Measured decay of the species in the clean, dark chamber (e.g. O3 2.1e-6, NO2 9.4e-7, HNO3 8e-5)."), stretch=1)
        b = QPushButton("+"); b.clicked.connect(lambda: self._add_wall("", "")); h.addWidget(b)
        b = QPushButton("−"); b.clicked.connect(self._del_wall); h.addWidget(b); v.addLayout(h)
        self.wall = QTableWidget(0, 2); self.wall.setHorizontalHeaderLabels(["species (MCM name)", "k [1/s]"])
        self.wall.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.wall.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows); self.wall.setMaximumHeight(140)
        for sp, k in (self.p.wall_loss_s1 or {}).items():
            self._add_wall(sp, f"{float(k):g}")
        v.addWidget(self.wall)
        lay.addWidget(box)

        # ── injection standards ──
        box = QGroupBox("Standards of the injection calculator"); g = QGridLayout(box)
        g.addWidget(QLabel("Mass-flow controller standard T [K]"), 0, 0); self.mfc_t = QLineEdit(f"{self.p.mfc_std_T_K:g}"); g.addWidget(self.mfc_t, 0, 1)
        g.addWidget(QLabel("standard p [hPa]"), 0, 2); self.mfc_p = QLineEdit(f"{self.p.mfc_std_P_hPa:g}"); g.addWidget(self.mfc_p, 0, 3)
        g.addWidget(QLabel("Gas-syringe filling T [K]"), 1, 0); self.syr_t = QLineEdit(f"{self.p.syringe_T_K:g}"); g.addWidget(self.syr_t, 1, 1)
        lay.addWidget(box)
        lay.addStretch()

        # ── buttons ──
        bb = QDialogButtonBox()
        self.save_btn = bb.addButton("Save", QDialogButtonBox.ButtonRole.AcceptRole)
        self.save_as_btn = bb.addButton("Save as new chamber …", QDialogButtonBox.ButtonRole.ActionRole)
        bb.addButton(QDialogButtonBox.StandardButton.Cancel)
        self.save_btn.clicked.connect(self._save); self.save_as_btn.clicked.connect(self._save_as); bb.rejected.connect(self.reject)
        if self.locked:
            self.save_btn.setEnabled(False); self.save_btn.setToolTip("SAPHIR is the reference profile; derive a new chamber instead")
            for w in body.findChildren((QLineEdit, QPlainTextEdit, QCheckBox, QComboBox, QTableWidget, QPushButton)):
                w.setEnabled(False)
        if as_new:
            self.save_as_btn.hide()
        outer.addWidget(bb)
        self._dilute_preview(); self._light_changed()

    # ── helpers ──
    def _add_wall(self, sp: str, k: str):
        r = self.wall.rowCount(); self.wall.insertRow(r)
        self.wall.setItem(r, 0, QTableWidgetItem(sp)); self.wall.setItem(r, 1, QTableWidgetItem(k))

    def _del_wall(self):
        if self.wall.currentRow() >= 0:
            self.wall.removeRow(self.wall.currentRow())

    def _dilute_preview(self):
        try:
            tmp = C.ChamberProfile(name="x", volume_m3=_float(self.vol), leak_flow_m3h=_float(self.leak), dilute_expr=self.dil.text().strip())
            self.dil_preview.setText("DILUTE = " + C.render_dilute(tmp))
        except ValueError:
            self.dil_preview.setText("DILUTE = ? (volume or leak flow is not a number)")

    def _light_changed(self):
        t = C.LIGHT_TYPES[self.ltype.currentIndex()]
        sun = t == "sun"; lamps = t == "lamps"
        for w in (self.lat_lab, self.lat, self.lon_lab, self.lon):
            w.setVisible(sun)
        for w in (self.lamps_lab, self.lamps):
            w.setVisible(lamps or sun)
        self.lamps_lab.setText("Lamps" if lamps else "Auxiliary lamps (if any)")
        for w in (self.jno2_lab, self.jno2):
            w.setVisible(lamps)
        self.light_note.setText({
            "sun": "Photolysis frequencies come from measured J values when the experiment has them, otherwise from the "
                   "clear-sky parametrisation of the mechanism at this site.",
            "lamps": "jNO2 at full power scales the lamp photolysis; the ratios of the other J values to jNO2 depend on the "
                     "lamp spectrum (lamp J table: next step; until then the sunlight spectrum scaled by jNO2 is used).",
            "dark": "No photolysis: HV stays 0 and all J values are 0 unless the experiment gives measured ones."}[t])

    # ── collect / save ──
    def profile(self) -> C.ChamberProfile:
        p = copy.deepcopy(self.p)
        p.id = C.safe_id(self.id.text()) if self.id.text().strip() else p.id
        p.name = self.name.text().strip() or p.name
        p.institution, p.location = self.inst.text().strip(), self.loc.text().strip()
        p.note = self.note.toPlainText().strip()
        p.references = [l.strip() for l in self.refs.toPlainText().splitlines() if l.strip()]
        p.verified = self.verified.isChecked()
        p.volume_m3 = _float(self.vol); p.surface_m2 = _float(self.surf); p.material = self.mat.text().strip()
        p.leak_flow_m3h = _float(self.leak); p.dilute_expr = self.dil.text().strip()
        t = C.LIGHT_TYPES[self.ltype.currentIndex()]
        light = {"type": t}
        if t == "sun":
            light.update(lat=_float(self.lat, C.SAPHIR_LAT), lon=_float(self.lon, C.SAPHIR_LON))
        light["lamps"] = self.lamps.text().strip()           # kept for every type: a sun chamber may have auxiliary lamps
        light["jno2_full_s1"] = _float(self.jno2)
        p.light = light
        p.archive = {"kind": "saphir_nc", "root": self.aroot.text().strip()} if self.akind.currentIndex() > 0 and self.aroot.text().strip() else {}
        p.effects = {k: cb.isChecked() for k, cb in self.effects.items()}
        wl = {}
        for r in range(self.wall.rowCount()):
            sp = self.wall.item(r, 0).text().strip() if self.wall.item(r, 0) else ""
            kv = self.wall.item(r, 1).text().strip() if self.wall.item(r, 1) else ""
            if sp and kv:
                wl[sp] = float(kv)
        p.wall_loss_s1 = wl
        p.mfc_std_T_K = _float(self.mfc_t, 273.15); p.mfc_std_P_hPa = _float(self.mfc_p, 1013.25); p.syringe_T_K = _float(self.syr_t, 293.15)
        return p

    def _checked_profile(self) -> Optional[C.ChamberProfile]:
        try:
            p = self.profile()
        except ValueError as e:
            QMessageBox.warning(self, "Chamber", f"A number could not be read: {e}"); return None
        if not p.name.strip():
            QMessageBox.warning(self, "Chamber", "Give the chamber a name."); return None
        if p.volume_m3 <= 0:
            QMessageBox.warning(self, "Chamber", "The volume must be a positive number of m³."); return None
        if p.light_type == "lamps" and not p.light.get("jno2_full_s1"):
            if QMessageBox.question(self, "Chamber", "No jNO2 at full lamp power given — the lamp photolysis cannot be scaled.\nSave anyway?") != QMessageBox.StandardButton.Yes:
                return None
        return p

    def _write(self, p: C.ChamberProfile) -> bool:
        if p.id == C.DEFAULT_CHAMBER:              # never from here: the self-test checks it against the built-in
            QMessageBox.warning(self, "Chamber", "The id SAPHIR is reserved for the reference profile; choose another id."); return False
        try:
            path = C.save_profile(p)
        except OSError as e:
            QMessageBox.critical(self, "Chamber", f"The profile could not be written:\n{e}"); return False
        self.saved_id = p.id
        self.setToolTip(path)
        return True

    def _save(self):
        p = self._checked_profile()
        if p is None:
            return
        if self.as_new:
            if not self.id.text().strip():
                QMessageBox.warning(self, "Chamber", "Give the new chamber an id (the file name)."); return
            if p.id in C.CHAMBERS:
                if QMessageBox.question(self, "Chamber", f"chambers/{p.id}{C.SUFFIX} exists. Replace it?") != QMessageBox.StandardButton.Yes:
                    return
        if self._write(p):
            self.accept()

    def _save_as(self):
        p = self._checked_profile()
        if p is None:
            return
        new_id, ok = QInputDialog.getText(self, "Save as new chamber", "Id of the new chamber (file name, letters/digits/_-.):",
                                          text=C.safe_id(p.name) if p.id == C.DEFAULT_CHAMBER else p.id + "_copy")
        if not ok or not new_id.strip():
            return
        p.id = C.safe_id(new_id)
        if p.id in C.CHAMBERS and QMessageBox.question(self, "Chamber", f"chambers/{p.id}{C.SUFFIX} exists. Replace it?") != QMessageBox.StandardButton.Yes:
            return
        if self.locked:
            p.verified = False
            p.name = p.name if p.name != "SAPHIR" else new_id
            p.note = (p.note + " " if p.note else "") + "Derived from the SAPHIR profile — edit to your chamber's values."
        if self._write(p):
            self.accept()
