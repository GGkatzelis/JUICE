"""
experiment_wizard.py — guided scenario builder (Step 3 of the JUICE workflow)
============================================================================

A QWizard that asks experiment questions and fills an `ExperimentSpec`
(experiment_spec.py).  On Finish the backend writes <name>.enz, <name>.easy,
the mechanism subset (if compounds were chosen) and <name>.experiment.json.
Its pages are the steps of juice_gui.py (shared through ExperimentHost).  The user never edits a
file; the review page shows the generated EASY and keeps an "expert overrides"
box and a button to the block-level editor (easy_builder.EasyBuilderDialog).

Pages:  Start (template, name, folder) → Where and when → Chemistry →
        Injections → Chamber effects → Initial conditions → Review
"""

from __future__ import annotations

import os
import re
from typing import Dict, List, Optional, Tuple

from PyQt6.QtCore import Qt, pyqtSignal, QThread, QStringListModel
from PyQt6.QtGui import QStandardItemModel, QStandardItem, QColor
from PyQt6.QtWidgets import (QSplitter, QWizard, QWizardPage, QVBoxLayout, QHBoxLayout, QGridLayout, QLabel,
                             QLineEdit, QPushButton, QCheckBox, QComboBox, QListWidget, QListWidgetItem,
                             QTableWidget, QTableWidgetItem, QPlainTextEdit, QFileDialog, QMessageBox,
                             QGroupBox, QHeaderView, QRadioButton, QButtonGroup, QWidget, QSpinBox,
                             QAbstractItemView, QDialog, QDialogButtonBox, QTreeWidget, QTreeWidgetItem, QTabWidget,
                             QStyledItemDelegate, QCompleter)

from experiment_spec import (ExperimentSpec, Conditions, InjectionEvent, ChamberEffects,
                             list_templates, load_template, load_spec, save_spec, validate,
                             mechanism_species, render_easy, easy_form, write_experiment,
                             H2O_IMPURITY_DEFAULTS, STANDARD_COLUMNS, ACTIONS, ACTION_ORDER, level_icon)
from easy_builder import scan_mechanisms, EasyBuilderDialog, DILUTE_DEF
import saphir_wiki
import injection_calc as ic
import mechanism_registry as mreg

try:
    from mcm_species_db import MCM_DB, reactions_of, search_species
except Exception:                               # database optional
    MCM_DB = {}


def _info(text: str) -> QLabel:
    lab = QLabel(text)
    lab.setWordWrap(True)
    lab.setStyleSheet("color: #888;")
    return lab


def _more_info(summary: str, bullets: List[str], button: str = "ⓘ How it works") -> QWidget:
    """One grey line + a small toggle that unfolds the details as bullet points (folded by default)."""
    w = QWidget(); v = QVBoxLayout(w); v.setContentsMargins(0, 0, 0, 0); v.setSpacing(2)
    h = QHBoxLayout(); h.setContentsMargins(0, 0, 0, 0)
    h.addWidget(_info(summary), stretch=1)
    b = QPushButton(button); b.setCheckable(True); b.setFlat(True); b.setCursor(Qt.CursorShape.PointingHandCursor)
    b.setStyleSheet("QPushButton { color: #6a4c93; border: none; padding: 0 4px; } QPushButton:checked { font-weight: bold; }")
    h.addWidget(b, alignment=Qt.AlignmentFlag.AlignTop); v.addLayout(h)
    details = QLabel("<ul style='margin-left:-20px'>" + "".join(f"<li>{t}</li>" for t in bullets) + "</ul>")
    details.setWordWrap(True); details.setStyleSheet("color: #666;"); details.setVisible(False)
    v.addWidget(details)
    b.toggled.connect(details.setVisible)
    w.details, w.toggle = details, b
    return w


WIKI_LINE_HELP = [
    "<b>Time</b> at the start of the line: <code>HH:MM</code> or <code>HH:MM:SS</code>.",
    "<b>Code</b>: one letter on its own after the time. Upper case starts, lower case stops: "
    "<code>R</code>/<code>r</code> roof open/close, <code>H</code>/<code>h</code> humidification, "
    "<code>S</code>/<code>s</code> seed aerosol, <code>X</code> VOC or mixture, <code>O</code> ozone, "
    "<code>M</code> NO, <code>D</code> NO2, <code>c</code> CO, <code>C</code> CO2 tracer. "
    "A stop written as <code>H … stop</code> is also understood.",
    "<b>The rest is free text.</b> Initials (F.R.) and words like CHANGE: are ignored. "
    "Read from it: the compound or mixture name, amounts (<code>23 ppb</code>, <code>ppm</code>, <code>180 µl</code>, "
    "<code>ml</code>, <code>sccm</code>, <code>FC 0.5l 62.5%</code>), bottle concentration <code>(1%)</code>, "
    "durations (<code>50 sec</code>, <code>5 min</code>) and <code>ppbC</code>.",
    "<b>No code?</b> Logs written as <code>08:01 A.N. CHANGE: Roof open</code> are read from the words: roof open / closed, "
    "humidification on / off, seed, O3 / NO / NO2 / CO / CO2 additions, and any other injection → X. "
    "<code>ZUSTAND:</code> / <code>STATUS:</code> lines and flushing, valves or fans stay notes. Check the Code column before importing.",
    "Lines without a time are skipped. Lines that do not describe an action are kept as notes and not modelled.",
]


# ── MCM species completer: type a name, a common name or a formula; pick the MCM species ──
_MCM_ENTRIES: List[str] = []
_MCM_ENTRY_TO_NAME: Dict[str, str] = {}


def mcm_completer_entries() -> List[str]:
    """'C5H8   C5H8   68.1 g/mol   C=CC(=C)C' for every MCM species, plus
    'isoprene → C5H8 …' for the common names of injection_calc.COMPOUNDS (built once)."""
    if _MCM_ENTRIES:
        return _MCM_ENTRIES
    for key, (mcm, mw, _dens) in sorted(ic.COMPOUNDS.items(), key=lambda kv: (len(MCM_DB[kv[1][0]].formula or "") if kv[1][0] in MCM_DB else 99, kv[0])):   # common names first, shortest formula first
        if mcm and mcm in MCM_DB:
            e = MCM_DB[mcm]
            txt = f"{key} → {mcm}   {e.formula or ''}   {mw:.1f} g/mol"
            _MCM_ENTRIES.append(txt); _MCM_ENTRY_TO_NAME[txt] = mcm
    # then the MCM species, shortest formula first: typing a formula shows the exact
    # match (C5H8) before the species that merely contain it (C5H8O5 …)
    for name, e in sorted(MCM_DB.items(), key=lambda kv: (len(kv[1].formula or ""), kv[0])):
        if getattr(e, "is_peroxy", False):
            continue                                   # radicals are never injected
        txt = f"{name}   {e.formula or ''}   {e.mass or 0:.1f} g/mol   {e.smiles or ''}"
        _MCM_ENTRIES.append(txt); _MCM_ENTRY_TO_NAME[txt] = name
    return _MCM_ENTRIES


class MCMSpeciesDelegate(QStyledItemDelegate):
    """Cell editor with an MCM completer: matches anywhere in name, common name, formula
    or SMILES; the dropdown shows formula, molar mass and SMILES; choosing an entry writes
    the MCM name into the cell.  Free text (mixtures, surrogates, 0.9*A+0.1*B) still works."""

    def createEditor(self, parent, option, index):
        ed = QLineEdit(parent)
        comp = QCompleter(QStringListModel(mcm_completer_entries(), ed), ed)
        comp.setFilterMode(Qt.MatchFlag.MatchContains)
        comp.setCaseSensitivity(Qt.CaseSensitivity.CaseInsensitive)
        comp.setCompletionMode(QCompleter.CompletionMode.PopupCompletion)
        comp.setMaxVisibleItems(12)
        comp.popup().setMinimumWidth(560)
        comp.activated.connect(lambda txt: ed.setText(_MCM_ENTRY_TO_NAME.get(txt, txt)))
        ed.setCompleter(comp)
        ed.setPlaceholderText("MCM name, common name or formula")
        return ed

    def setModelData(self, editor, model, index):
        txt = editor.text().strip()
        model.setData(index, _MCM_ENTRY_TO_NAME.get(txt, txt))


class MultiSelectCombo(QComboBox):
    """A dropdown whose entries are tick boxes: several can be on at once, the closed
    combo shows what is on, each entry has a tooltip with the details."""
    changed = pyqtSignal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self._model = QStandardItemModel(self); self.setModel(self._model)
        self._model.itemChanged.connect(self._on_item_changed)
        self.setEditable(True); self.lineEdit().setReadOnly(True)
        self.lineEdit().setPlaceholderText("none")
        self._ids: List[str] = []
        self.view().pressed.connect(self._toggle_row)
        self._updating = False

    def set_items(self, items, empty_text: str = "none"):
        """items: [(id, label, tooltip, checked)]"""
        self._updating = True
        self._model.clear(); self._ids = []
        for i, label, tip, checked in items:
            it = QStandardItem(label); it.setFlags(Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsUserCheckable)
            it.setCheckState(Qt.CheckState.Checked if checked else Qt.CheckState.Unchecked); it.setToolTip(tip)
            self._model.appendRow(it); self._ids.append(i)
        self._updating = False
        self.lineEdit().setPlaceholderText(empty_text if not items else "none")
        self._refresh_text()

    def _toggle_row(self, index):
        it = self._model.itemFromIndex(index)
        if it is not None:
            it.setCheckState(Qt.CheckState.Unchecked if it.checkState() == Qt.CheckState.Checked else Qt.CheckState.Checked)

    def _on_item_changed(self, _it):
        if not self._updating:
            self._refresh_text(); self.changed.emit()

    def checked_ids(self) -> List[str]:
        return [self._ids[r] for r in range(self._model.rowCount()) if self._model.item(r).checkState() == Qt.CheckState.Checked]

    def set_checked(self, ids: List[str]):
        self._updating = True
        for r in range(self._model.rowCount()):
            self._model.item(r).setCheckState(Qt.CheckState.Checked if self._ids[r] in ids else Qt.CheckState.Unchecked)
        self._updating = False
        self._refresh_text()

    def _refresh_text(self):
        ids = self.checked_ids()
        self.lineEdit().setText(", ".join(ids) if ids else "")
        self.setToolTip("\n".join(self._model.item(r).text() for r in range(self._model.rowCount())
                                   if self._model.item(r).checkState() == Qt.CheckState.Checked) or "no add-on chosen")

    def hidePopup(self):                       # keep the list open while ticking several entries
        if self.view().underMouse():
            return
        super().hidePopup()


def _table(headers: List[str]) -> QTableWidget:
    t = QTableWidget(0, len(headers))
    t.setHorizontalHeaderLabels(headers)
    t.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
    t.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
    return t


def _rows(t: QTableWidget) -> List[List[str]]:
    out = []
    for r in range(t.rowCount()):
        vals = [(t.item(r, c).text().strip() if t.item(r, c) else "") for c in range(t.columnCount())]
        if any(vals):
            out.append(vals)
    return out


def _set_rows(t: QTableWidget, rows: List[List[str]]):
    t.setRowCount(0)
    for vals in rows:
        r = t.rowCount(); t.insertRow(r)
        for c, v in enumerate(vals):
            t.setItem(r, c, QTableWidgetItem(str(v)))


def _add_row(t: QTableWidget, vals: Optional[List[str]] = None):
    r = t.rowCount(); t.insertRow(r)
    for c, v in enumerate(vals or []):
        t.setItem(r, c, QTableWidgetItem(str(v)))


def _del_row(t: QTableWidget):
    r = t.currentRow()
    if r >= 0:
        t.removeRow(r)


# ─────────────────────────────────────────────────────────────────────────────
class StartPage(QWizardPage):
    def __init__(self, wiz: "ExperimentHost"):
        super().__init__()
        self.wiz = wiz
        self.setTitle("Experiment")
        self.setSubTitle("Name the experiment, choose where its files go, then load what you already have "
                         "(a saved experiment, the SAPHIR wiki page, or pasted notes). A template is optional.")
        lay = QVBoxLayout(self)
        g = QGridLayout()
        g.addWidget(QLabel("Experiment name"), 0, 0)
        self.name_edit = QLineEdit("my_experiment"); g.addWidget(self.name_edit, 0, 1, 1, 2)
        g.addWidget(QLabel("Project folder"), 1, 0)
        self.out_edit = QLineEdit(wiz.easy_root); g.addWidget(self.out_edit, 1, 1)
        b = QPushButton("…"); b.clicked.connect(self._browse); g.addWidget(b, 1, 2)
        lay.addLayout(g)
        lay.addWidget(_info("A sub-folder with the experiment name is created inside the project folder for the "
                            "generated input files and the experiment description (.experiment.json)."))
        lay.addWidget(QLabel("<b>Start from</b>"))
        hb = QHBoxLayout()
        ob = QPushButton("Open a saved experiment …"); ob.clicked.connect(self._open_saved); hb.addWidget(ob)
        wb = QPushButton("Import experiment log …"); wb.clicked.connect(lambda: self.wiz.import_from_wiki(self)); hb.addWidget(wb)
        wb.setToolTip("Fetch the day's SAPHIR wiki page, load a saved page or text file, or paste your own notes "
                      "('HH:MM <code> …' lines) — all end in the same action-log check")
        hb.addStretch(); lay.addLayout(hb)
        lay.addWidget(_info("Wiki import and pasted notes read the action log (roof, humidification, injections), "
                            "switch to the archive for the measured conditions when the day folder exists, and leave "
                            "everything editable in the following steps."))
        hb2 = QHBoxLayout(); hb2.addWidget(QLabel("Template (optional)"))
        self.tpl_combo = QComboBox()
        self.templates = list_templates()
        self.tpl_combo.addItem("Empty — SAPHIR defaults", "")
        self._tpl_desc: Dict[str, str] = {}
        for name, path in self.templates.items():
            try:
                title, desc, _ = load_template(path)
            except Exception as e:
                title, desc = name, f"(unreadable: {e})"
            self.tpl_combo.addItem(title, path); self._tpl_desc[path] = desc
        self.tpl_combo.currentIndexChanged.connect(self._template_chosen)
        hb2.addWidget(self.tpl_combo, stretch=1); lay.addLayout(hb2)
        self.tpl_info = _info("Templates pre-fill the typical SAPHIR day (humidification, VOC, ozone, roof times). "
                              "Choosing one replaces the current entries.")
        lay.addWidget(self.tpl_info)
        lay.addStretch()

    def _browse(self):
        d = QFileDialog.getExistingDirectory(self, "Project folder", self.out_edit.text() or self.wiz.easy_root)
        if d:
            self.out_edit.setText(d)

    def _template_chosen(self, _idx: int = 0):
        self.wiz.wiki_memory = None                    # another experiment: forget the last import
        path = self.tpl_combo.currentData()
        if path:
            try:
                _t, desc, spec = load_template(path)
            except Exception as e:
                QMessageBox.warning(self, "Template", str(e)); return
            self.tpl_info.setText(desc)
        else:
            spec = ExperimentSpec()
            self.tpl_info.setText("SAPHIR defaults, nothing pre-filled.")
        spec.name = self.name_edit.text().strip() or spec.name
        self.wiz.load_spec(spec, keep_name=True)

    def _open_saved(self):
        p, _ = QFileDialog.getOpenFileName(self, "Open experiment", self.wiz.easy_root,
                                           "Experiment (*.experiment.json);;All files (*)")
        if not p:
            return
        try:
            spec = load_spec(p)
        except Exception as e:
            QMessageBox.critical(self, "Open failed", str(e)); return
        self.wiz.load_spec(spec)
        self.name_edit.setText(spec.name)
        self.out_edit.setText(os.path.dirname(os.path.dirname(p)) or self.wiz.easy_root)

    def validatePage(self) -> bool:
        if not re.fullmatch(r"[A-Za-z0-9_.\-]+", self.name_edit.text().strip()):
            QMessageBox.warning(self, "Name", "Use letters, digits, '_' or '-' for the experiment name."); return False
        return True


# ─────────────────────────────────────────────────────────────────────────────
class ConditionsPage(QWizardPage):
    def __init__(self, wiz):
        super().__init__()
        self.wiz = wiz
        self.setTitle("Where and when")
        self.setSubTitle("Date and time window, where the chamber conditions come from, and the roof schedule.")
        lay = QVBoxLayout(self)
        g = QGridLayout(); r = 0
        g.addWidget(QLabel("Date (UTC)"), r, 0); self.date = QLineEdit("2024-07-09"); self.date.textChanged.connect(self._archive_status); g.addWidget(self.date, r, 1)
        g.addWidget(QLabel("from"), r, 2); self.start = QLineEdit("04:00"); g.addWidget(self.start, r, 3)
        g.addWidget(QLabel("to"), r, 4); self.end = QLineEdit("21:59"); g.addWidget(self.end, r, 5)
        g.addWidget(QLabel("step [s]"), r, 6); self.step = QSpinBox(); self.step.setRange(10, 3600); self.step.setValue(60); g.addWidget(self.step, r, 7); r += 1
        g.addWidget(QLabel("Chamber"), r, 0)
        self.chamber = QComboBox(); self.chamber.addItems(list(ic.CHAMBERS)); g.addWidget(self.chamber, r, 1)
        self.chamber.setToolTip("chambers/<id>.chamber.json — volume, light source, data archive, default chamber effects")
        hc = QHBoxLayout(); hc.setContentsMargins(0, 0, 0, 0)
        be = QPushButton("Edit …"); be.setToolTip("change the values of this chamber profile (SAPHIR: read-only, derive a copy)"); be.clicked.connect(self._edit_chamber); hc.addWidget(be)
        bn = QPushButton("New …"); bn.setToolTip("add your own chamber, pre-filled with a copy of the selected one"); bn.clicked.connect(self._new_chamber); hc.addWidget(bn)
        g.addLayout(hc, r, 2)
        self.chamber_info = _info(""); g.addWidget(self.chamber_info, r, 3, 1, 5)
        self._loading = False                     # set_from() in progress: a chamber change then does not reset the chamber effects
        self.chamber.currentTextChanged.connect(self._chamber_changed); r += 1
        lay.addLayout(g)

        lay.addWidget(QLabel("<b>Chamber conditions</b> — temperature, pressure, humidity, flow and photolysis frequencies "
                             "drive the rate constants and the air density; the model needs them for every minute."))
        self.src_archive = QRadioButton("Measured — from the SAPHIR archive (T, P, humidity, flow, J values, NO, NO2, O3)")
        self.src_file = QRadioButton("Measured — from a file you provide (CSV or Excel: time + T, p, humidity, flow, roof, J values, NO, NO2, O3)")
        self.src_manual = QRadioButton("Manual — constant values, clear-sky photolysis (no measurements available)")
        for rb in (self.src_archive, self.src_file, self.src_manual):
            rb.toggled.connect(self._source_changed)
        lay.addWidget(self.src_archive)
        ha = QHBoxLayout(); ha.addSpacing(24); ha.addWidget(QLabel("Archive root"))
        self.archive = QLineEdit(r"Z:\IEK8-SAPHIR\data"); self.archive.textChanged.connect(self._archive_status); ha.addWidget(self.archive, stretch=1)
        self.archive_status = QLabel(""); ha.addWidget(self.archive_status); lay.addLayout(ha)
        lay.addWidget(self.src_file)
        hf = QHBoxLayout(); hf.addSpacing(24); hf.addWidget(QLabel("Input file"))
        self.input_path = QLineEdit(); self.input_path.setPlaceholderText("measured conditions table — see 'Template…' for the columns and units")
        self.input_path.textChanged.connect(self._file_status); hf.addWidget(self.input_path, stretch=1)
        bf = QPushButton("Browse …"); bf.clicked.connect(self._browse_input); hf.addWidget(bf)
        tf = QPushButton("Template …"); tf.setToolTip("write an example CSV with the recognised columns, units and time formats"); tf.clicked.connect(self._write_template); hf.addWidget(tf)
        lay.addLayout(hf)
        self.file_status = _info(""); self.file_status.setContentsMargins(24, 0, 0, 0); lay.addWidget(self.file_status)
        lay.addWidget(self.src_manual)

        self.manual_box = QGroupBox("Manual values")
        mg = QGridLayout(self.manual_box)
        mg.addWidget(QLabel("Temperature [K]"), 0, 0); self.temp = QLineEdit("298"); mg.addWidget(self.temp, 0, 1)
        mg.addWidget(QLabel("Pressure [hPa]"), 0, 2); self.press = QLineEdit("1006.5"); mg.addWidget(self.press, 0, 3)
        mg.addWidget(QLabel("Water vapour [vol-%]"), 1, 0); self.h2o = QLineEdit("1.0"); mg.addWidget(self.h2o, 1, 1)
        mg.addWidget(QLabel("Replenishment flow [m³/h]"), 1, 2); self.flow = QLineEdit("6.0"); mg.addWidget(self.flow, 1, 3)
        lay.addWidget(self.manual_box)

        hr = QHBoxLayout(); self.roof_label = QLabel("<b>Roof open / close</b> (HH:MM, UTC)"); hr.addWidget(self.roof_label)
        a = QPushButton("+"); a.clicked.connect(lambda: _add_row(self.roof, ["09:30", "18:30"])); hr.addWidget(a)
        d = QPushButton("−"); d.clicked.connect(lambda: _del_row(self.roof)); hr.addWidget(d); hr.addStretch(); lay.addLayout(hr)
        self.roof = _table(["open", "close"]); self.roof.setMaximumHeight(100); lay.addWidget(self.roof)

        self.adv = QGroupBox("Advanced"); self.adv.setCheckable(True); self.adv.setChecked(False)
        ag = QGridLayout(self.adv)
        ag.addWidget(QLabel("Photolysis frequencies"), 0, 0)
        self.jsrc = QComboBox(); self.jsrc.addItems(["measured (archive J values)", "clear sky from solar zenith angle"]); ag.addWidget(self.jsrc, 0, 1)
        ag.addWidget(QLabel("Temperature in the model"), 1, 0)
        self.temp_mode = QComboBox(); self.temp_mode.addItems(["follow the measured / tabulated temperature", "hold constant at the manual value"]); ag.addWidget(self.temp_mode, 1, 1)
        ag.addWidget(QLabel("Hold species to measured values"), 2, 0)
        self.constrain = QLineEdit(""); self.constrain.setPlaceholderText("e.g. NO2, O3  (archive only; the model then does not predict them)"); ag.addWidget(self.constrain, 2, 1)
        self.adv_body = [self.jsrc, self.temp_mode, self.constrain]
        self.adv.toggled.connect(lambda on: [w.setVisible(on) for w in self.adv_body + [ag.itemAt(i).widget() for i in range(ag.count()) if isinstance(ag.itemAt(i).widget(), QLabel)]])
        lay.addWidget(self.adv)
        self.adv.toggled.emit(False)
        self.manual_box.setVisible(False)
        self._archive_status()

    def _edit_chamber(self):
        self._chamber_dialog(as_new=False)

    def _new_chamber(self):
        self._chamber_dialog(as_new=True)

    def _chamber_dialog(self, as_new: bool):
        from chamber_editor import ChamberDialog
        from chambers import render_dilute
        before = ic.chamber(self.chamber.currentText())
        dlg = ChamberDialog(before, self, as_new=as_new)
        if dlg.exec() != QDialog.DialogCode.Accepted or not dlg.saved_id:
            return
        after = ic.chamber(dlg.saved_id)
        same_defaults = (dlg.saved_id == before.id and after.effects == before.effects
                         and render_dilute(after) == render_dilute(before))
        self.refresh_chambers(select=dlg.saved_id, apply_defaults=not same_defaults)

    def refresh_chambers(self, select: str = "", apply_defaults: bool = True):
        """Re-fill the chamber list from chambers.CHAMBERS (after a profile was saved) and select
        `select`; apply_defaults=False keeps the experiment's chamber effects (only the profile's
        texts changed)."""
        self.chamber.blockSignals(True)
        self.chamber.clear(); self.chamber.addItems(list(ic.CHAMBERS))
        self.chamber.setCurrentText(select if select in ic.CHAMBERS else ic.DEFAULT_CHAMBER)
        self.chamber.blockSignals(False)
        was = self._loading
        self._loading = not apply_defaults
        try:
            self._chamber_changed(self.chamber.currentText())
        finally:
            self._loading = was

    def _chamber_changed(self, name: str):
        p = ic.chamber(name)
        flag = "" if p.verified else " Numbers from the public record, not verified — check them (Edit …)."
        self.chamber_info.setText(f"{p.name}: {p.summary()}. " + (p.note.split(" — ")[0] if p.note else "") + flag)
        self.src_archive.setText(f"Measured — from the {p.name} data archive (T, P, humidity, flow, J values, NO, NO2, O3)" if p.has_archive
                                 else f"Measured — from a data archive (none that JUICE can read for {p.name})")
        self.src_archive.setEnabled(p.has_archive)
        lt = p.light_type
        self.roof_label.setText({"sun": "<b>Roof open / close</b> (HH:MM, UTC)",
                                 "lamps": "<b>Lights on / off</b> (HH:MM, UTC; a 'lights' column 0–1 in the input file wins)",
                                 "dark": "<b>Light</b> — dark chamber profile: the model runs without light"}[lt])
        self.roof.setEnabled(lt != "dark")
        self.src_manual.setText("Manual — constant values, clear-sky photolysis (no measurements available)" if lt == "sun" else
                                "Manual — constant values; photolysis from the chamber profile's lamps (jNO2 at full power × lights)" if lt == "lamps" else
                                "Manual — constant values (dark chamber)")
        self.jsrc.setEnabled(lt == "sun")
        self.jsrc.setToolTip("" if lt == "sun" else "lamp chambers: jNO2 measured or from the profile, the other J values from the lamp spectrum (lamp_photolysis)")
        if p.has_archive and (self._loading is False or not self.archive.text().strip()):
            self.archive.setText(p.archive_root)
        if not p.has_archive and self.src_archive.isChecked():
            self.src_file.setChecked(True)
        if not self._loading and hasattr(self.wiz, "chamber_selected"):
            self.wiz.chamber_selected(p.id)       # the user picked another chamber: its default effects and DILUTE

    def _source_changed(self, _on: bool = False):
        archive_on, file_on = self.src_archive.isChecked(), self.src_file.isChecked()
        self.manual_box.setVisible(not (archive_on or file_on))
        self.archive.setEnabled(archive_on)
        self.input_path.setEnabled(file_on)
        self.jsrc.setCurrentIndex(0 if (archive_on or file_on) else 1)
        if file_on:
            self._file_status()

    def _browse_input(self):
        p, _ = QFileDialog.getOpenFileName(self, "Measured inputs table", self.input_path.text() or "", "Tables (*.csv *.txt *.xlsx *.xls);;All files (*)")
        if p:
            self.input_path.setText(p); self.src_file.setChecked(True)

    def _write_template(self):
        import measured_inputs as mi
        p, _ = QFileDialog.getSaveFileName(self, "Write the measured-inputs template", "measured_inputs_template.csv", "CSV (*.csv)")
        if p:
            mi.write_template(p)
            self.file_status.setText(f"template written: {p} — fill it with your data and choose it as the input file")

    def _file_status(self):
        p = self.input_path.text().strip()
        if not p:
            self.file_status.setText("no file chosen" if self.src_file.isChecked() else ""); return
        if not os.path.isfile(p):
            self.file_status.setText("✖ file not found"); return
        try:
            import measured_inputs as mi
            inp = mi.cached_inputs(p, self.date.text().strip(), float(self.press.text() or 1006.5))
            txt = "✔ " + inp.summary() + ("; " + "; ".join(inp.notes) if inp.notes else "")
            self.file_status.setText(txt)
        except Exception as e:
            self.file_status.setText(f"✖ cannot read: {e}")

    def _archive_status(self):
        if not self.date.text().strip():
            self.archive_status.setText("enter a date"); self.archive_status.setStyleSheet("color: #888;"); return
        day = saphir_wiki.archive_day_folder(self.archive.text().strip(), self.date.text().strip())
        if day:
            self.archive_status.setText("✔ data found for this date"); self.archive_status.setStyleSheet("color: #2e9e4f;")
        else:
            self.archive_status.setText("✖ no archive folder for this date"); self.archive_status.setStyleSheet("color: #c0392b;")

    def set_from(self, c: Conditions):
        self._loading = True
        try:
            self._set_from(c)
        finally:
            self._loading = False

    def _set_from(self, c: Conditions):
        if self.chamber.count() != len(ic.CHAMBERS):              # a profile was added meanwhile
            self.chamber.blockSignals(True); self.chamber.clear(); self.chamber.addItems(list(ic.CHAMBERS)); self.chamber.blockSignals(False)
        self.chamber.setCurrentText(c.chamber_name if c.chamber_name in ic.CHAMBERS else ic.DEFAULT_CHAMBER)
        self._chamber_changed(self.chamber.currentText())
        self.date.setText(c.date); self.start.setText(c.start); self.end.setText(c.end); self.step.setValue(int(c.step_s))
        self.archive.setText(c.archive_root); self.temp.setText(f"{c.temperature_K:g}")
        self.press.setText(f"{c.pressure_hPa:g}"); self.h2o.setText(f"{c.h2o_percent:g}"); self.flow.setText(f"{c.flow_m3h:g}")
        self.input_path.setText(c.input_file)
        if c.source in ("archive", "manual", "file"):
            {"archive": self.src_archive, "manual": self.src_manual, "file": self.src_file}[c.source].setChecked(True)
        else:                                          # blank experiment: nothing chosen yet
            for rb in (self.src_archive, self.src_file, self.src_manual):
                rb.setAutoExclusive(False); rb.setChecked(False); rb.setAutoExclusive(True)
            self.manual_box.setVisible(False)
        self.jsrc.setCurrentIndex(0 if c.j_source == "measured" else 1)
        self.temp_mode.setCurrentIndex(0 if c.temperature_mode == "enz" else 1)
        _set_rows(self.roof, [[o, cl] for o, cl in c.roof_open])
        self.constrain.setText(", ".join(c.constrain_species))
        self._archive_status()

    def get(self) -> Conditions:
        c = Conditions()
        c.date, c.start, c.end, c.step_s = self.date.text().strip(), self.start.text().strip(), self.end.text().strip(), int(self.step.value())
        c.source = "archive" if self.src_archive.isChecked() else ("file" if self.src_file.isChecked() else ("manual" if self.src_manual.isChecked() else ""))
        c.archive_root = self.archive.text().strip()
        c.input_file = self.input_path.text().strip()
        c.temperature_K = float(self.temp.text() or 298); c.temperature_mode = "enz" if self.temp_mode.currentIndex() == 0 else "const"
        c.pressure_hPa = float(self.press.text() or 1006.5); c.h2o_percent = float(self.h2o.text() or 0); c.flow_m3h = float(self.flow.text() or 0)
        c.j_source = "measured" if self.jsrc.currentIndex() == 0 else "clear_sky"
        c.roof_open = [(o, cl) for o, cl in _rows(self.roof)]
        c.constrain_species = [x.strip() for x in self.constrain.text().split(",") if x.strip()]
        c.chamber_name = self.chamber.currentText()
        return c


# ─────────────────────────────────────────────────────────────────────────────
class ChemistryPage(QWizardPage):
    def __init__(self, wiz):
        super().__init__()
        self.wiz = wiz
        self.setTitle("Chemistry")
        self.setSubTitle("Which compounds take part? The matching part of the MCM is extracted in the background, "
                         "or pick an existing mechanism include set.")
        lay = QVBoxLayout(self)
        # which chemistry: the MCM download (unmodified) and the add-ons loaded on top of it
        mg = QGridLayout()
        mg.addWidget(QLabel("Base mechanism"), 0, 0)
        self.base_combo = QComboBox()
        self._bases = mreg.load_bases()
        for b in self._bases.values():
            self.base_combo.addItem(f"{b.name}  [{b.status}]", b.id)
        self.base_combo.currentIndexChanged.connect(self._base_changed); mg.addWidget(self.base_combo, 0, 1)
        self.base_info = _info(""); mg.addWidget(self.base_info, 0, 2)
        mg.addWidget(QLabel("Add-ons"), 1, 0)
        self.addon_combo = MultiSelectCombo(); self.addon_combo.changed.connect(self._addons_changed)
        mg.addWidget(self.addon_combo, 1, 1)
        self.addon_info = _info("none chosen"); mg.addWidget(self.addon_info, 1, 2)
        lay.addLayout(mg)
        self._base_changed()
        self.mode_compounds = QRadioButton("Compounds (MCM subset is generated)")
        self.mode_existing = QRadioButton("Existing mechanism include set")
        self.mode_inorganic = QRadioButton("Inorganic chemistry only (no VOC)")
        self.mode_inorganic.setToolTip("NOx, O3, HOx, CO, SO2 chemistry and photolysis plus the chamber's own HONO and HCHO sources —\n"
                                       "for tests such as dilution, NO/O3 titration, CO or zero-air experiments.")
        self.mode_compounds.setChecked(True)
        hb = QHBoxLayout(); hb.addWidget(self.mode_compounds); hb.addWidget(self.mode_existing); hb.addWidget(self.mode_inorganic); hb.addStretch(); lay.addLayout(hb)
        g = QGridLayout()
        g.addWidget(QLabel("Search MCM species (name or formula)"), 0, 0)
        self.search = QLineEdit(); self.search.textChanged.connect(self._filter); g.addWidget(self.search, 0, 1)
        self.candidates = QListWidget(); g.addWidget(self.candidates, 1, 0, 1, 2)
        self.candidates.itemDoubleClicked.connect(lambda _it: self._add_selected())
        bb = QHBoxLayout(); a = QPushButton("Add →"); a.clicked.connect(self._add_selected); bb.addWidget(a)
        d = QPushButton("← Remove"); d.clicked.connect(self._remove_selected); bb.addWidget(d); bb.addStretch()
        g.addLayout(bb, 2, 0, 1, 2)
        g.addWidget(QLabel("Chosen compounds"), 0, 2)
        self.chosen = QListWidget(); g.addWidget(self.chosen, 1, 2, 2, 1)
        lay.addLayout(g)
        hb2 = QHBoxLayout(); hb2.addWidget(QLabel("Existing mechanism"))
        self.mech_combo = QComboBox(); self.mech_combo.setEditable(True)
        self.mech_combo.addItems([""] + list(scan_mechanisms(wiz.include_dir))); hb2.addWidget(self.mech_combo, stretch=1)
        lay.addLayout(hb2)
        lay.addWidget(_info("Compounds not in the MCM (e.g. emission mixtures) are represented by MCM surrogates; "
                            "the CHANEL template shows how a mixture is written as mole fractions."))

        # first-generation products of the chosen compounds (from the full MCM .fac): one block per
        # compound, one column per reaction partner, products listed downwards, rate on hover
        tb = QGroupBox("First-generation chemistry of the chosen compounds (untick to hide)"); tb.setCheckable(True); tb.setChecked(True)
        tl = QVBoxLayout(tb)
        self.tree_tabs = QTabWidget(); self.tree_tabs.setDocumentMode(True); self.tree_tabs.setMaximumHeight(230)
        tl.addWidget(self.tree_tabs)
        hb3 = QHBoxLayout(); b = QPushButton("Update"); b.clicked.connect(self.update_tree); hb3.addWidget(b)
        hb3.addWidget(_info("One tab per compound: reaction partners side by side, products below each. "
                            "Hover a product set for the MCM rate; double-click it to adjust the rate below."), stretch=1)
        tl.addLayout(hb3)
        self._tree_widgets = [self.tree_tabs] + [hb3.itemAt(i).widget() for i in range(hb3.count()) if hb3.itemAt(i).widget()]
        tb.toggled.connect(self._toggle_tree)
        lay.addWidget(tb)

        # reaction adjustments = EASY k[...] overrides
        ab = QGroupBox("Reaction adjustments (written exactly as EASY does it)"); al = QVBoxLayout(ab)
        self.adjust = _table(["Reaction (as in the mechanism)", "New rate expression — 0 removes the reaction", "Comment"])
        self.adjust.setMinimumHeight(110); al.addWidget(self.adjust)
        hb4 = QHBoxLayout()
        self.rxn_search = QLineEdit(); self.rxn_search.setPlaceholderText("species name, e.g. CISOPA — lists its MCM reactions")
        hb4.addWidget(self.rxn_search, stretch=1)
        fb = QPushButton("Find reactions"); fb.clicked.connect(self._find_reactions); hb4.addWidget(fb)
        ab_ = QPushButton("+ row"); ab_.clicked.connect(lambda: _add_row(self.adjust, ["", "", ""])); hb4.addWidget(ab_)
        db_ = QPushButton("− row"); db_.clicked.connect(lambda: _del_row(self.adjust)); hb4.addWidget(db_)
        al.addLayout(hb4)
        self.found = QListWidget(); self.found.setMaximumHeight(90); self.found.itemDoubleClicked.connect(self._adjust_from_found); al.addWidget(self.found)
        al.addWidget(_info("Each row becomes  k[reaction] = CONST(expression)  after the mechanism includes; the last definition "
                           "wins, so the MCM rate is replaced. A rate of 0 switches the reaction off (FACSIMILE evaluates it as zero, "
                           "JUICE drops the reaction). Write rates in EASY notation, e.g. 2.70D-11*EXP(390/T)*0.5."))
        lay.addWidget(ab)
        self._filter("")

    # ── base / add-ons ──
    def base_id(self) -> str:
        return self.base_combo.currentData() or mreg.DEFAULT_BASE

    def addon_ids(self) -> List[str]:
        return self.addon_combo.checked_ids()

    def _addons_changed(self):
        ids = self.addon_ids()
        self.addon_info.setText("none chosen — the base mechanism as published" if not ids else
                                "on top of the base: " + ", ".join(ids))

    def _base_changed(self, _i: int = 0):
        b = self._bases.get(self.base_id())
        if b is not None:
            probs = b.verify()
            self.base_info.setText((f"v{b.version}, downloaded {b.source.get('downloaded', '?')[:60]}" if not probs
                                    else "⚠ " + "; ".join(probs)[:120]))
        keep = self.addon_ids()
        addons = mreg.addons_for(self.base_id())
        self.addon_combo.set_items([(a.id, f"{a.name}  v{a.version} [{a.status}]",
                                     f"{a.description}\n\nfiles: {', '.join(a.files)}\nadds: {', '.join(a.adds_species) or '-'}\n"
                                     f"needs from base: {', '.join(a.needs_base_species) or '-'}\nstatus: {a.status}, version {a.version}",
                                     a.id in keep) for a in addons],
                                    empty_text="none installed for this base (mechanisms/addons/<id>/manifest.json)")
        self._addons_changed()

    # ── first-generation tree / adjustments ──
    def _fac(self) -> str:
        from experiment_spec import base_fac
        return base_fac(self.base_id())

    def _toggle_tree(self, on: bool):
        for w in self._tree_widgets:
            w.setVisible(on)
        if on and self.tree_tabs.count() == 0:
            self.update_tree()

    def update_tree(self):
        self.tree_tabs.clear()
        names = [self.chosen.item(i).text() for i in range(self.chosen.count())]
        for n in names:
            try:
                rxns = reactions_of(n, self._fac())
            except Exception as e:
                self.tree_tabs.addTab(_info(f"cannot read the MCM file: {e}"), n); continue
            if not rxns:
                self.tree_tabs.addTab(_info("no reactions found (not an MCM species?)"), n); continue
            order = ["OH", "O3", "NO3", "hv", ""]
            partners = sorted({r["partner"] for r in rxns}, key=lambda p: (order.index(p) if p in order else 9, p))
            cols = {p: [r for r in rxns if r["partner"] == p] for p in partners}
            t = QTableWidget(max(len(v) for v in cols.values()), len(partners))
            t.setHorizontalHeaderLabels([f"+ {p}  ({len(cols[p])})" if p else f"alone ({len(cols[p])})" for p in partners])
            t.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch); t.verticalHeader().setVisible(False)
            t.verticalHeader().setDefaultSectionSize(19); t.setStyleSheet("font-size: 10px;")
            t.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers); t.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
            for c, p in enumerate(partners):
                for r_, rx in enumerate(cols[p]):
                    prods = " + ".join(rx["products"])
                    it = QTableWidgetItem("→ " + (prods if len(prods) <= 34 else prods[:31] + "…"))
                    it.setToolTip(rx["easy"] + "\nk = " + rx["easy_rate"])
                    it.setData(Qt.ItemDataRole.UserRole, (rx["easy"], rx["easy_rate"])); t.setItem(r_, c, it)
            t.itemDoubleClicked.connect(self._adjust_from_tree)
            self.tree_tabs.addTab(t, f"{n} ({len(rxns)})")

    def _adjust_from_tree(self, item, _col=None):
        data = item.data(Qt.ItemDataRole.UserRole)
        if data:
            _add_row(self.adjust, [data[0], data[1], "from the MCM; edit the expression or set 0"])

    def _find_reactions(self):
        self.found.clear()
        q = self.rxn_search.text().strip()
        if not q:
            return
        try:
            rxns = reactions_of(q, self._fac())
        except Exception as e:
            self.found.addItem(f"cannot read the MCM file: {e}"); return
        for r in rxns:
            it = QListWidgetItem(f"{r['easy']}      k = {r['easy_rate']}"); it.setData(Qt.ItemDataRole.UserRole, (r["easy"], r["easy_rate"])); self.found.addItem(it)
        if not rxns:
            self.found.addItem(f"no MCM reactions with {q.upper()} as a reactant")

    def _adjust_from_found(self, item):
        data = item.data(Qt.ItemDataRole.UserRole)
        if data:
            _add_row(self.adjust, [data[0], data[1], "from the MCM; edit the expression or set 0"])

    def _filter(self, text: str):
        self.candidates.clear()
        q = text.strip().upper()
        n = 0
        for name, e in MCM_DB.items():
            formula = getattr(e, "formula", "") or ""
            if not q or q in name or q in formula.upper():
                self.candidates.addItem(f"{name}   {formula}   {getattr(e, 'mass', '') or ''}")
                n += 1
                if n >= 300:
                    break

    def _add_selected(self):
        for it in self.candidates.selectedItems():
            name = it.text().split()[0]
            if not any(self.chosen.item(i).text() == name for i in range(self.chosen.count())):
                self.chosen.addItem(name)
        self.update_tree()

    def _remove_selected(self):
        for it in self.chosen.selectedItems():
            self.chosen.takeItem(self.chosen.row(it))
        self.update_tree()

    def set_from(self, spec: ExperimentSpec):
        i = self.base_combo.findData(spec.base_mechanism)
        self.base_combo.setCurrentIndex(i if i >= 0 else 0)
        self._base_changed()
        self.addon_combo.set_checked(list(spec.addons)); self._addons_changed()
        self.chosen.clear()
        for p in spec.precursors:
            self.chosen.addItem(p)
        self.mech_combo.setCurrentText(spec.mechanism_stem)
        (self.mode_existing if spec.mechanism_stem else self.mode_inorganic if spec.inorganic_only else self.mode_compounds).setChecked(True)
        _set_rows(self.adjust, [[a, b, c] for a, b, c in spec.reaction_adjustments])
        self.update_tree()

    def apply(self, spec: ExperimentSpec):
        spec.inorganic_only = self.mode_inorganic.isChecked()
        if self.mode_existing.isChecked():
            spec.mechanism_stem = self.mech_combo.currentText().strip(); spec.precursors = []
        elif spec.inorganic_only:
            spec.mechanism_stem = ""; spec.precursors = []
        else:
            spec.mechanism_stem = ""; spec.precursors = [self.chosen.item(i).text() for i in range(self.chosen.count())]
        spec.reaction_adjustments = [(a, b, c) for a, b, c in _rows(self.adjust) if a.strip()]
        spec.base_mechanism = self.base_id()
        spec.addons = self.addon_ids()

    def validatePage(self) -> bool:
        if self.mode_existing.isChecked() and not self.mech_combo.currentText().strip():
            QMessageBox.warning(self, "Chemistry", "Choose a mechanism include set."); return False
        if self.mode_compounds.isChecked() and self.chosen.count() == 0:
            QMessageBox.warning(self, "Chemistry", "Add at least one compound."); return False
        return True


# ─────────────────────────────────────────────────────────────────────────────
INJ_HEADERS = ["Action", "Name", "Product(s) — empty = recorded only", "Total [ppb]", "kOH [s⁻¹]", "Start HH:MM", "Duration [min]",
               "ENZ column", "Switch", "Comment"]
KOH_COL = 4          # computed, read-only: OH reactivity of the row's amount (oh_reactivity.injection_koh)
_ACTION_LABELS = [ACTIONS[k]["label"] for k in ACTION_ORDER]


def _action_combo(key: str) -> QComboBox:
    cb = QComboBox(); cb.addItems(_ACTION_LABELS)
    cb.setCurrentIndex(ACTION_ORDER.index(key) if key in ACTION_ORDER else ACTION_ORDER.index("other"))
    return cb


def _inj_cols() -> List[int]:
    """Editable columns after Action, in the order of the row values (the kOH column is skipped)."""
    return [c for c in range(1, len(INJ_HEADERS)) if c != KOH_COL]


def _koh_item(text: str = "", tip: str = "") -> QTableWidgetItem:
    it = QTableWidgetItem(text); it.setFlags(it.flags() & ~Qt.ItemFlag.ItemIsEditable)
    it.setForeground(QColor("#555555")); it.setTextAlignment(int(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter))
    if tip:
        it.setToolTip(tip)
    return it


def _inj_add_row(t: QTableWidget, action: str, vals: List[str]):
    """Append an injection row: combo box in the Action column, text in the others, kOH computed."""
    r = t.rowCount(); t.insertRow(r)
    t.setCellWidget(r, 0, _action_combo(action))
    t.blockSignals(True)
    for c, v in zip(_inj_cols(), vals):
        t.setItem(r, c, QTableWidgetItem(str(v)))
    t.setItem(r, KOH_COL, _koh_item())
    t.blockSignals(False)
    if t.rowCount():
        t.itemChanged.emit(t.item(r, 3) or t.item(r, 1))


def _inj_rows(t: QTableWidget) -> List[List[str]]:
    """[action_key, name, products, ppb, start, dur, column, switch, comment] per non-empty row."""
    out = []
    for r in range(t.rowCount()):
        cb = t.cellWidget(r, 0)
        key = ACTION_ORDER[cb.currentIndex()] if isinstance(cb, QComboBox) else "other"
        vals = [(t.item(r, c).text().strip() if t.item(r, c) else "") for c in _inj_cols()]
        if any(vals):
            out.append([key] + vals)
    return out


class InjectionsPage(QWizardPage):
    def __init__(self, wiz):
        super().__init__()
        self.wiz = wiz
        self.setTitle("Injections")
        self.setSubTitle("What was injected, how much in total, when and over how long. One row per product channel; "
                         "rows with the same name belong to the same injection (e.g. impurities of a humidification).")
        lay = QVBoxLayout(self)
        self._delivery: Dict[str, Dict[str, str]] = {}      # name → how the amount was derived (kept with the spec)
        self.table = _table(INJ_HEADERS); lay.addWidget(self.table)
        self.table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeMode.ResizeToContents)
        self.table.verticalHeader().setDefaultSectionSize(26)
        self.table.setMinimumHeight(26 * 10 + 30)          # ten rows visible without scrolling
        self._species_delegate = MCMSpeciesDelegate(self.table)
        self.table.setItemDelegateForColumn(2, self._species_delegate)   # product: MCM name, common name or formula
        self.table.horizontalHeaderItem(KOH_COL).setToolTip(
            "OH reactivity of this row's amount just after the injection: Σ k_OH(T) × concentration, fully mixed, before any reaction.\n"
            "k_OH from the base MCM's own rate expressions at the experiment temperature. Computed, not an input.")
        self.table.itemChanged.connect(self._koh_changed)
        hb = QHBoxLayout()
        hb.addWidget(QLabel("Add action:"))
        self.action_pick = QComboBox(); self.action_pick.addItems(_ACTION_LABELS); hb.addWidget(self.action_pick)
        for label, fn in [("+ add", self._add_action), ("+ compounds from Chemistry", self._add_compounds),
                          ("+ mixture …", self._mixture_dialog), ("Edit mixture …", self._edit_mixture),
                          ("Amount calculator …", self._amount_dialog),
                          ("Import experiment log …", lambda: self.wiz.import_from_wiki(self)),
                          ("− row", lambda: _del_row(self.table))]:
            b = QPushButton(label); b.clicked.connect(fn); hb.addWidget(b)
        hb.addStretch(); lay.addLayout(hb)
        self.koh_total = QLabel(""); self.koh_total.setStyleSheet("color:#555"); self.koh_total.setWordWrap(True); lay.addWidget(self.koh_total)
        self.help = _more_info("Amount = total ppb injected, spread over the duration. Product = MCM species, or mole fractions (0.9*A+0.1*B).", [
            "<b>Product cell</b>: type a name (isoprene), an MCM name or a formula (C5H8) and pick from the list "
            "(formula, molar mass and SMILES are shown).",
            "<b>Empty product</b> = the injection is recorded in the timeline but not modelled.",
            "<b>ENZ column / Switch</b>: leave empty for the standard names (QNO, QO3, QH2O, Q_&lt;name&gt;).",
            "<b>Same name, empty start</b> = an extra product channel of that injection. The humidification rows are the "
            "impurities the water brings in (HONO 0.002 ppb, CO 0.5 ppb; NO/NO2 0 = placeholders, fill in if measured).",
            "<b>XNOy</b> is the NOy bookkeeping tracer, not a chemical species.",
            "<b>Amount calculator …</b> turns µl by syringe, ml of bottle gas or MFC flow × time into ppb for the selected row and keeps the derivation.",
            "<b>+ mixture …</b> expands a multi-compound solution into MCM species with mole fractions.",
            "<b>Import experiment log …</b> reads the SAPHIR wiki action log (fetch, saved page or pasted text)."])
        lay.addWidget(self.help)
        self.chain_box = QGroupBox("Delayed release (advanced)"); self.chain_box.setCheckable(True); self.chain_box.setChecked(False)
        cl = QVBoxLayout(self.chain_box)
        cl.addWidget(_info("Normally the injected amount appears in the chamber during the injection pulse. Some injections "
                           "evaporate or mix more slowly: the IDL scenarios then let the pulse fill a holding species that "
                           "converts to the real compound with a first-order rate. Example from the ethanol scenario: the pulse "
                           "makes 0.9*VOCin1 + 0.1*VOCin2; VOCin1 --> VOCin_all at 1./60 (1 min), VOCin2 --> VOCin_all at "
                           "1./1800 (30 min), VOCin_all --> C2H5OH at 1./60. Put the holding species in the product of the "
                           "injection row and list the conversion steps here. Leave empty for an instantaneous injection."))
        self.chain = _table(["Reactants", "Products", "Rate expression [1/s]"]); self.chain.setMaximumHeight(110); cl.addWidget(self.chain)
        hb2 = QHBoxLayout(); a = QPushButton("+"); a.clicked.connect(lambda: _add_row(self.chain)); d = QPushButton("−"); d.clicked.connect(lambda: _del_row(self.chain)); hb2.addWidget(a); hb2.addWidget(d); hb2.addStretch(); cl.addLayout(hb2)
        self._chain_widgets = [cl.itemAt(i).widget() for i in range(cl.count()) if cl.itemAt(i).widget()]
        self.chain_box.toggled.connect(lambda on: [w.setVisible(on) for w in self._chain_widgets])
        self.chain_box.toggled.emit(False)
        lay.addWidget(self.chain_box)

    # ── OH reactivity per row ──
    def _koh_conditions(self) -> Tuple[float, float]:
        try:
            c = self.wiz.page_cond.get()
            return float(c.temperature_K or 298.0), float(c.pressure_hPa or 1013.25)
        except Exception:
            return 298.0, 1013.25

    def _koh_changed(self, item: QTableWidgetItem):
        if item is None or item.column() not in (1, 2, 3):
            return
        self._koh_row(item.row()); self._koh_sum()

    def _koh_row(self, r: int, T: float = None, p: float = None):
        import oh_reactivity as R
        if T is None:
            T, p = self._koh_conditions()
        prod = self.table.item(r, 2).text().strip() if self.table.item(r, 2) else ""
        ppb = self.table.item(r, 3).text().strip() if self.table.item(r, 3) else ""
        text, tip = "", ""
        if prod and ppb:
            try:
                val, missing = R.injection_koh(prod, ppb, T, p)
            except Exception as e:
                val, missing = None, [str(e)]
            if val is not None:
                text = "0" if val == 0 else (f"{val:.3g}" if val >= 0.01 else f"{val:.1e}")
                tip = f"{ppb} ppb of {prod} at {T:.0f} K, {p:.0f} hPa"
                if missing:
                    text += " *"; tip += "\n* no OH reaction in the base MCM for: " + ", ".join(missing) + " (counted as 0)"
        self.table.blockSignals(True)
        self.table.setItem(r, KOH_COL, _koh_item(text, tip))
        self.table.blockSignals(False)

    def _koh_sum(self):
        total, n = 0.0, 0
        for r in range(self.table.rowCount()):
            it = self.table.item(r, KOH_COL)
            try:
                total += float(it.text().replace("*", "").strip()); n += 1
            except (AttributeError, ValueError):
                continue
        T, p = self._koh_conditions()
        self.koh_total.setText(f"OH reactivity of all injections together: {total:.3g} s⁻¹ — the sum of each injection just after it is made, "
                               f"fully mixed and before any reaction or dilution, at {T:.0f} K (rows with an amount and an MCM product: {n})."
                               if n else "")

    def refresh_koh(self):
        T, p = self._koh_conditions()
        for r in range(self.table.rowCount()):
            self._koh_row(r, T, p)
        self._koh_sum()

    def _add_compounds(self):
        page = self.wiz.page_chem
        names = [page.chosen.item(i).text() for i in range(page.chosen.count())]
        if not names and page.mech_combo.currentText().strip():
            QMessageBox.information(self, "Compounds", "Add rows by hand: the compounds are defined by the existing mechanism."); return
        for n in names:
            _inj_add_row(self.table, "voc", [n, n, "10", "08:00", "2", "", "", ""])

    def _add_action(self):
        """Append the rows of the chosen action with its SAPHIR defaults (start time left for the user)."""
        key = ACTION_ORDER[self.action_pick.currentIndex()]
        d = ACTIONS[key]
        existing = {r[1] for r in _inj_rows(self.table)}
        name, k = d["name"], 1
        while name and name in existing:               # second NO injection -> NO_2 with its own column/switch
            k += 1; name = f"{d['name']}_{k}"
        suffix = f"_{k}" if k > 1 else ""
        col = d["column"] + suffix if d["column"] else ""
        sw = d["switch"] + suffix if d["switch"] else ""
        _inj_add_row(self.table, key, [name, d["products"], d["total_ppb"], "", f"{d['duration_min']:g}", col, sw, d["comment"]])
        for prod, ppb in d.get("extra_channels", []):
            _inj_add_row(self.table, key, [name, prod, ppb, "", "", col, sw, f"↳ extra channel of '{name}'"])
        self.table.scrollToBottom()

    def _selected_name(self) -> str:
        r = self.table.currentRow()
        it = self.table.item(r, 1) if r >= 0 else None
        return it.text().strip() if it else ""

    def _amount_dialog(self):
        r = self.table.currentRow()
        if r < 0 or not self._selected_name():
            QMessageBox.information(self, "Amount calculator", "Select the injection row first."); return
        name = self._selected_name()
        prod = self.table.item(r, 2).text().strip() if self.table.item(r, 2) else ""
        cond = self.wiz.page_cond.get()
        dlg = AmountDialog(self, compound=prod or name, delivery=self._delivery.get(name, {}), cond=cond)
        if dlg.exec() != QDialog.DialogCode.Accepted or dlg.result_ppb is None:
            return
        self.table.setItem(r, 3, QTableWidgetItem(f"{dlg.result_ppb:.4g}"))
        old = self.table.item(r, 9).text() if self.table.item(r, 9) else ""
        old = re.sub(r"\s*\|\s*calc:.*$", "", old)
        self.table.setItem(r, 9, QTableWidgetItem((old + " | " if old else "") + "calc: " + dlg.result_expl))
        self._delivery[name] = dlg.result_delivery

    def _edit_mixture(self):
        """Selected mixture row → its recipe (surrogates editable) → row and precursors replaced."""
        r = self.table.currentRow()
        name = self._selected_name()
        d = self._delivery.get(name, {}) if name else {}
        if r < 0 or d.get("mode") != "mixture":
            QMessageBox.information(self, "Edit mixture", "Select a mixture row (one that came from a recipe) first."); return
        m = ic.load_mixtures().get(d.get("mixture", ""))
        if m is None:
            QMessageBox.information(self, "Edit mixture", f"Recipe '{d.get('mixture', '')}' not found in {ic.MIXTURE_DIR}"); return
        cond = self.wiz.page_cond.get()
        dlg = MixtureDialog(self, cond=cond)
        dlg.preload(m, d.get("injected_ul", ""), self.table.item(r, 5).text() if self.table.item(r, 5) else "",
                    self.table.item(r, 6).text() if self.table.item(r, 6) else "")
        if dlg.exec() != QDialog.DialogCode.Accepted or dlg.result_event is None:
            return
        ev = dlg.result_event
        self.table.setItem(r, 2, QTableWidgetItem(ev.products)); self.table.setItem(r, 3, QTableWidgetItem(ev.total_ppb))
        self.table.setItem(r, 9, QTableWidgetItem(ev.comment))
        self._delivery[name] = ev.delivery
        chem = self.wiz.page_chem
        have = {chem.chosen.item(i).text() for i in range(chem.chosen.count())}
        for sp in dlg.result_precursors:
            if sp not in have:
                chem.chosen.addItem(sp)

    def _mixture_dialog(self):
        cond = self.wiz.page_cond.get()
        dlg = MixtureDialog(self, cond=cond)
        if dlg.exec() != QDialog.DialogCode.Accepted or dlg.result_event is None:
            return
        ev, precursors = dlg.result_event, dlg.result_precursors
        existing = {r_[1] for r_ in _inj_rows(self.table)}
        if ev.name in existing:
            k = 2
            while f"{ev.name}_{k}" in existing:
                k += 1
            ev.name, ev.column = f"{ev.name}_{k}", f"{ev.column}_{k}"
        _inj_add_row(self.table, "mixture", [ev.name, ev.products, ev.total_ppb, ev.start, f"{ev.duration_min:g}",
                                              ev.column, ev.switch, ev.comment])
        self._delivery[ev.name] = ev.delivery
        chem = self.wiz.page_chem
        chem.mode_compounds.setChecked(True)
        have = {chem.chosen.item(i).text() for i in range(chem.chosen.count())}
        for sp in precursors:
            if sp not in have:
                chem.chosen.addItem(sp)
        self.table.scrollToBottom()

    def set_from(self, spec: ExperimentSpec):
        self.table.setRowCount(0)
        self._delivery = {inj.name: dict(inj.delivery) for inj in spec.injections}
        for inj in spec.injections:
            key = inj.action_key()
            _inj_add_row(self.table, key, [inj.name, inj.products, inj.total_ppb, inj.start, f"{inj.duration_min:g}",
                                           inj.column, inj.switch, inj.comment])
            for prod, ppb in inj.extra_channels:
                _inj_add_row(self.table, key, [inj.name, prod, ppb, "", "", inj.column, inj.switch, f"↳ extra channel of '{inj.name}'"])
        _set_rows(self.chain, [[a, b, c] for inj in spec.injections for a, b, c in inj.chain])
        self.chain_box.setChecked(self.chain.rowCount() > 0)
        self.refresh_koh()

    def get(self) -> List[InjectionEvent]:
        events: Dict[str, InjectionEvent] = {}
        order: List[str] = []
        for action, name, prod, ppb, start, dur, col, sw, comment in _inj_rows(self.table):
            if not name:
                continue
            ev = events.get(name)
            if ev is None:
                ev = InjectionEvent(name, prod, ppb or "0", start, float(dur or 1), col, sw, comment=comment, action=action,
                                    delivery=dict(self._delivery.get(name, {})))
                events[name] = ev; order.append(name)
            else:
                if not ev.start and start:
                    ev.start, ev.duration_min = start, float(dur or 1)
                if prod or ppb:
                    if not ev.products and not ev.extra_channels and prod:
                        ev.products, ev.total_ppb = prod, ppb or "0"
                    else:
                        ev.extra_channels.append((prod, ppb or "0"))
                ev.column = ev.column or col; ev.switch = ev.switch or sw
        evs = [events[n] for n in order]
        chain = [(a, b, c) for a, b, c in _rows(self.chain)]
        if chain and evs:
            evs[-1].chain = chain          # chains are attached to the last modelled injection
            for ev in reversed(evs):
                if ev.products:
                    evs[-1].chain = []; ev.chain = chain; break
        return evs


# ─────────────────────────────────────────────────────────────────────────────
CHAMBER_ITEMS = [
    ("walldep", "Wall deposition of HNO3, HONO, H2O2, HO2, O3, N2O5 (saphir.walldep.incl)",
     "First-order loss to the chamber walls with the standard SAPHIR deposition rates."),
    ("hono", "Photolytic HONO source from the walls (Rohrer et al. 2005)",
     "The Teflon film releases HONO when illuminated; depends on humidity, temperature and jNO2. Main OH source in SAPHIR."),
    ("hcho", "Photolytic HCHO source from the walls",
     "Small light-dependent formaldehyde release."),
    ("bgreac", "Background OH reactivity",
     "Unmeasured trace organics consuming OH, expressed as a CO-equivalent amount Y (ppb)."),
    ("dilution", "Dilution by the replenishment flow",
     "The chamber is kept at slight overpressure; all species are diluted with the flow (DILUTE from FLOW, the leak flow and the chamber volume of the profile)."),
    ("dark_hono", "Dark HONO formation from water on the walls",
     "Slow heterogeneous HONO production without light."),
    ("light_source", "NO2 source under illumination (SNO2)",
     "Light-driven NO2 release from the walls, used in some low-NOx experiments."),
]


def _dark_factor_of(rate: str) -> str:
    """'1e-012' → '1', '1e-012/10' → '0.1', '1e-012*0' → '0'; anything else → '1'."""
    r = (rate or "").replace(" ", "")
    m = re.fullmatch(r"1e-0?12(?:/([0-9.]+)|\*([0-9.]+))?", r)
    if not m:
        return "1"
    if m.group(1):
        return f"{1.0 / float(m.group(1)):g}"
    if m.group(2):
        return f"{float(m.group(2)):g}"
    return "1"


class ChamberPage(QWizardPage):
    def __init__(self, wiz):
        super().__init__()
        self.wiz = wiz
        self.setTitle("What the chamber adds on its own")
        self.setSubTitle("SAPHIR chamber effects with the values of the IDL scenarios. Two numbers are tuned per experiment: "
                         "the strength of the HONO wall source and the background OH reactivity. Everything else is a switch; "
                         "the underlying formulas are only in the expert tool.")
        lay = QVBoxLayout(self)
        self.checks: Dict[str, QCheckBox] = {}
        self.hono_factor = QLineEdit("2.6*0.7"); self.hono_factor.setFixedWidth(90)
        self.y_inline = QLineEdit("250"); self.y_inline.setFixedWidth(70)
        self.dark_factor = QLineEdit("1"); self.dark_factor.setFixedWidth(60)
        for key, title, desc in CHAMBER_ITEMS:
            cb = QCheckBox(title); cb.setChecked(key != "light_source"); self.checks[key] = cb
            row = QHBoxLayout(); row.addWidget(cb)
            if key == "hono":
                row.addSpacing(12); row.addWidget(QLabel("strength factor")); row.addWidget(self.hono_factor)
                desc += " Fitted per experiment to the measured HONO or OH on a clean day (2.0 to 2.6*1.1*0.7 in the reference scenes)."
            elif key == "bgreac":
                row.addSpacing(12); row.addWidget(self.y_inline); row.addWidget(QLabel("ppb CO-equivalent"))
                desc += " Set from the measured OH reactivity of the clean chamber (50 to 250 ppb in the reference scenes)."
            elif key == "dark_hono":
                row.addSpacing(12); row.addWidget(QLabel("×")); row.addWidget(self.dark_factor)
                desc += " Factor on the standard rate (1 = IDL default, 0.1 as in scene_09)."
            row.addStretch(); lay.addLayout(row); lay.addWidget(_info("      " + desc))
        self.wall_loss_cb = QCheckBox("First-order wall losses of the chamber profile"); lay.addWidget(self.wall_loss_cb)
        self.wall_loss_info = _info("      none in this profile"); lay.addWidget(self.wall_loss_info)
        self._wall_loss_s1: Dict[str, float] = {}
        self.params = QGroupBox("Parameters (advanced)"); self.params.setCheckable(True); self.params.setChecked(False)
        g = QGridLayout(self.params)
        self.hono_const = QLineEdit("2.6*0.7"); self.hono_lhs = QLineEdit("SHONO+HV"); self.hono_rhs = QLineEdit("HONO+XNOy")
        g.addWidget(QLabel("HONO source: SHONO = CONST("), 0, 0); g.addWidget(self.hono_const, 0, 1); g.addWidget(QLabel(") reaction"), 0, 2); g.addWidget(self.hono_lhs, 0, 3); g.addWidget(QLabel("-->"), 0, 4); g.addWidget(self.hono_rhs, 0, 5)
        self.hcho_const = QLineEdit("1"); self.hcho_lhs = QLineEdit("qSHCHO+HV"); self.hcho_rhs = QLineEdit("HCHO"); self.hcho_cname = QLineEdit("SHCHO")
        g.addWidget(QLabel("HCHO source:"), 1, 0); g.addWidget(self.hcho_cname, 1, 1); g.addWidget(self.hcho_const, 1, 2); g.addWidget(self.hcho_lhs, 1, 3); g.addWidget(QLabel("-->"), 1, 4); g.addWidget(self.hcho_rhs, 1, 5)
        self.y_ppb = QLineEdit("250"); self.y_name = QLineEdit("qY"); self.y_prod = QLineEdit("CH3O2")
        g.addWidget(QLabel("Background reactivity Y [ppb CO-eq]"), 2, 0); g.addWidget(self.y_ppb, 2, 1); g.addWidget(QLabel("name / product"), 2, 2); g.addWidget(self.y_name, 2, 3); g.addWidget(self.y_prod, 2, 5)
        self.dilute = QLineEdit(DILUTE_DEF); self.tracer = QLineEdit("XNOy, NA")
        g.addWidget(QLabel("DILUTE ="), 3, 0); g.addWidget(self.dilute, 3, 1, 1, 3); g.addWidget(QLabel("tracer dilution"), 3, 4); g.addWidget(self.tracer, 3, 5)
        self.dark_rate = QLineEdit("1e-012/10"); self.dark_rhs = QLineEdit("HONO+XNOy")
        g.addWidget(QLabel("Dark HONO: H2O --> "), 4, 0); g.addWidget(self.dark_rhs, 4, 1); g.addWidget(QLabel("k = CONST("), 4, 2); g.addWidget(self.dark_rate, 4, 3)
        self.light_rate = QLineEdit("jno2/5e-3/3600*HV_in/1.6*0.55"); self.light_rhs = QLineEdit("NO2+XNOy"); self.light_init = QLineEdit("M*1e-9*3.5*0")
        g.addWidget(QLabel("SNO2 --> "), 5, 0); g.addWidget(self.light_rhs, 5, 1); g.addWidget(QLabel("k = CONST("), 5, 2); g.addWidget(self.light_rate, 5, 3); g.addWidget(QLabel("INITIAL[SNO2]"), 5, 4); g.addWidget(self.light_init, 5, 5)
        self.params.hide()                       # formulas stay in the spec (template / expert tool), not on this page

    def embed_initials(self, page: "InitialsPage"):
        """Collapsed by default: one line of explanation; tick the box to open the table
        (it opens by itself when a loaded experiment has initial values)."""
        box = QGroupBox("Initial values (tick to set)"); box.setCheckable(True); box.setChecked(False)
        bl = QVBoxLayout(box)
        bl.addWidget(_info(page.subTitle())); bl.addWidget(page)
        box.toggled.connect(page.setVisible); page.setVisible(False)
        self.initials_box = box
        self.layout().addWidget(box)

    def set_from(self, ch: ChamberEffects):
        for key, cb in self.checks.items():
            cb.setChecked(bool(getattr(ch, key)))
        self.hono_const.setText(ch.hono_const); self.hono_lhs.setText(ch.hono_lhs); self.hono_rhs.setText(ch.hono_rhs)
        self.hcho_const.setText(ch.hcho_const); self.hcho_lhs.setText(ch.hcho_lhs); self.hcho_rhs.setText(ch.hcho_rhs); self.hcho_cname.setText(ch.hcho_const_name)
        self.y_ppb.setText(ch.y_expr or ch.y_ppb); self.y_name.setText(ch.y_name); self.y_prod.setText(ch.y_product)
        self.dilute.setText(ch.dilute_expr); self.tracer.setText(", ".join(ch.tracer_dilution))
        self.dark_rate.setText(ch.dark_hono_rate); self.dark_rhs.setText(ch.dark_hono_rhs)
        self.light_rate.setText(ch.light_source_rate); self.light_rhs.setText(ch.light_source_rhs); self.light_init.setText(ch.light_source_initial)
        self.hono_factor.setText(ch.hono_const); self.y_inline.setText(ch.y_expr or ch.y_ppb)
        self.dark_factor.setText(_dark_factor_of(ch.dark_hono_rate))
        self._wall_loss_s1 = {str(k): float(v) for k, v in (ch.wall_loss_s1 or {}).items()}
        self.wall_loss_cb.setChecked(bool(ch.wall_loss) and bool(self._wall_loss_s1)); self.wall_loss_cb.setEnabled(bool(self._wall_loss_s1))
        self.wall_loss_info.setText("      " + (", ".join(f"{sp} {k:g} 1/s" for sp, k in self._wall_loss_s1.items()) + " (edit the chamber profile to change them)"
                                               if self._wall_loss_s1 else "none in this profile (Edit … the chamber on 'Where and when' to add some)"))

    def get(self, base: ChamberEffects) -> ChamberEffects:
        ch = ChamberEffects(**{k: getattr(base, k) for k in ChamberEffects.__dataclass_fields__})
        for key, cb in self.checks.items():
            setattr(ch, key, cb.isChecked())
        ch.wall_loss, ch.wall_loss_s1 = self.wall_loss_cb.isChecked(), dict(self._wall_loss_s1)
        ch.hono_const, ch.hono_lhs, ch.hono_rhs = self.hono_const.text().strip(), self.hono_lhs.text().strip(), self.hono_rhs.text().strip()
        ch.hcho_const, ch.hcho_lhs, ch.hcho_rhs, ch.hcho_const_name = self.hcho_const.text().strip(), self.hcho_lhs.text().strip(), self.hcho_rhs.text().strip(), self.hcho_cname.text().strip()
        yv = self.y_ppb.text().strip()
        if re.fullmatch(r"[0-9.]+", yv):
            ch.y_ppb, ch.y_expr = yv, ""
        else:
            ch.y_expr = yv
        ch.y_name, ch.y_product = self.y_name.text().strip(), self.y_prod.text().strip()
        ch.dilute_expr = self.dilute.text().strip(); ch.tracer_dilution = [t.strip() for t in self.tracer.text().split(",") if t.strip()]
        ch.dark_hono_rate, ch.dark_hono_rhs = self.dark_rate.text().strip(), self.dark_rhs.text().strip()
        ch.light_source_rate, ch.light_source_rhs, ch.light_source_initial = self.light_rate.text().strip(), self.light_rhs.text().strip(), self.light_init.text().strip()
        # the two tuned numbers and the dark-HONO factor come from the inline fields
        ch.hono_const = self.hono_factor.text().strip() or ch.hono_const
        yv = self.y_inline.text().strip()
        if yv:
            if re.fullmatch(r"[0-9.]+", yv):
                ch.y_ppb, ch.y_expr = yv, ""
            else:
                ch.y_expr = yv
        f = self.dark_factor.text().strip() or "1"
        ch.dark_hono_rate = "1e-012" if f in ("1", "1.", "1.0") else f"1e-012*{f}"
        return ch


# ─────────────────────────────────────────────────────────────────────────────
class InitialsPage(QWizardPage):
    def __init__(self, wiz):
        super().__init__()
        self.wiz = wiz
        self.setTitle("Initial values")
        self.setSubTitle("Concentrations at the start of the time window; usually empty for a flushed chamber. "
                         "Everything not listed starts at the floor value.")
        lay = QVBoxLayout(self)
        hb = QHBoxLayout(); hb.addWidget(QLabel("Floor for all species, INITIAL[*] = CONST(")); self.all_edit = QLineEdit("1."); hb.addWidget(self.all_edit); hb.addWidget(QLabel(") molecules/cm³")); hb.addStretch(); lay.addLayout(hb)
        self.table = _table(["Species", "Expression inside CONST( ), e.g. M*1e-9*14.0 for 14 ppb"]); lay.addWidget(self.table)
        self._species_delegate = MCMSpeciesDelegate(self.table); self.table.setItemDelegateForColumn(0, self._species_delegate)
        hb2 = QHBoxLayout(); a = QPushButton("+ row"); a.clicked.connect(lambda: _add_row(self.table, ["", "M*1e-9*"])); d = QPushButton("− row"); d.clicked.connect(lambda: _del_row(self.table)); hb2.addWidget(a); hb2.addWidget(d); hb2.addStretch(); lay.addLayout(hb2)
        lay.addWidget(_info("Typical: CO 14 ppb (M*1e-9*14.0); with archive data the measured NO2/O3 at the start can be used."))

    def set_from(self, spec: ExperimentSpec):
        self.all_edit.setText(spec.initial_all)
        _set_rows(self.table, [[sp, ex] for sp, ex in spec.initials])
        box = getattr(self.wiz.page_chamber, "initials_box", None) if hasattr(self.wiz, "page_chamber") else None
        if box is not None:
            box.setChecked(bool(spec.initials))          # open only when there is something to show

    def apply(self, spec: ExperimentSpec):
        spec.initial_all = self.all_edit.text().strip() or "1."
        spec.initials = [(sp, ex) for sp, ex in _rows(self.table) if sp]


# ─────────────────────────────────────────────────────────────────────────────
class ReviewPage(QWizardPage):
    def __init__(self, wiz):
        super().__init__()
        self.wiz = wiz
        self.setTitle("Review")
        self.setSubTitle("Check the summary and the generated EASY file. Finish writes all files and opens the scenario.")
        lay = QVBoxLayout(self)
        lay.addWidget(QLabel("Summary"))
        self.summary = QPlainTextEdit(); self.summary.setReadOnly(True); self.summary.setMaximumHeight(170); lay.addWidget(self.summary)
        lay.addWidget(QLabel("Validation — ✖ must be fixed before Finish, ⚠ is a hint, ✔ all good"))
        self.checks = QListWidget(); self.checks.setMaximumHeight(110); lay.addWidget(self.checks)
        hb = QHBoxLayout(); hb.addWidget(QLabel("Generated EASY file (read-only preview)")); hb.addStretch()
        eb = QPushButton("Open expert editor …"); eb.clicked.connect(self._expert); hb.addWidget(eb); lay.addLayout(hb)
        self.preview = QPlainTextEdit(); self.preview.setReadOnly(True); lay.addWidget(self.preview, stretch=1)
        lay.addWidget(QLabel("Expert overrides — EASY statements appended verbatim (last definition wins)"))
        self.overrides = QPlainTextEdit(); self.overrides.setMaximumHeight(80); lay.addWidget(self.overrides)
        self.overrides.textChanged.connect(self.refresh)

    def initializePage(self):
        self.refresh()

    def refresh(self):
        try:
            spec = self.wiz.collect_spec()
            self.preview.setPlainText(render_easy(spec))
            known = mechanism_species(self.wiz.include_dir, spec.mechanism_stem) if spec.mechanism_stem else None
            msgs = validate(spec, known)
            self.checks.clear()
            self._has_error = False
            for level, m in msgs:
                it = QListWidgetItem(level_icon(level) + m)
                self.checks.addItem(it)
                self._has_error |= level == "error"
            if not any(lvl != "info" for lvl, _m in msgs):
                self.checks.addItem(QListWidgetItem("✔ No problems found"))
            self.summary.setPlainText(self._summary(spec))
        except Exception as e:
            self.preview.setPlainText(f"Could not build the scenario: {e}")
            self._has_error = True

    @staticmethod
    def _summary(s: ExperimentSpec) -> str:
        c = s.conditions
        lines = [f"{s.name}: {c.date} {c.start}–{c.end} UTC, step {c.step_s} s, conditions from {c.source}",
                 f"Chemistry: {'MCM subset for ' + ', '.join(s.precursors) if s.precursors else s.mechanism_stem}",
                 f"Roof open: {', '.join(f'{o}–{cl}' for o, cl in c.roof_open) or 'never'};  photolysis: {c.j_source}",
                 f"T {c.temperature_K:g} K ({c.temperature_mode}), P {c.pressure_hPa:g} hPa, H2O {c.h2o_percent:g} %, flow {c.flow_m3h:g} m³/h"]
        for inj in s.injections:
            what = f"{inj.total_ppb} ppb → {inj.products}" if inj.products else "recorded only"
            lines.append(f"  injection {inj.name} at {inj.start or '--:--'} ({inj.duration_min:g} min): {what}"
                         + (f"; +{len(inj.extra_channels)} channel(s)" if inj.extra_channels else ""))
        on = [k for k, _t, _d in CHAMBER_ITEMS if getattr(s.chamber, k)]
        lines.append("Chamber effects: " + ", ".join(on))
        return "\n".join(lines)

    def _expert(self):
        try:
            spec = self.wiz.collect_spec()
            dlg = EasyBuilderDialog(parent=self, include_dir=self.wiz.include_dir, easy_root=self.wiz.easy_root)
            dlg.set_form(easy_form(spec))
            dlg.out_edit.setText(os.path.join(self.wiz.out_dir(), spec.name + ".easy"))
            dlg.exec()
        except Exception as e:
            QMessageBox.critical(self, "Expert editor", str(e))

    def validatePage(self) -> bool:
        if getattr(self, "_has_error", False):
            QMessageBox.warning(self, "Review", "Fix the errors (✖) before finishing."); return False
        return True


# ─────────────────────────────────────────────────────────────────────────────
class AmountDialog(QDialog):
    """Injected amount → ppb with the derivation shown.  Modes: stated ppb, liquid by
    syringe (µl, molar mass, density), gas from a bottle by syringe (ml, bottle
    concentration), MFC from a bottle (flow, time, bottle concentration)."""
    MODES = [("ppb", "Mixing ratio stated (ppb)"), ("liquid", "Liquid by syringe (µl)"),
             ("gas_syringe", "Gas from bottle by syringe (ml)"), ("mfc", "Gas from bottle via MFC (flow × time)")]

    def __init__(self, parent=None, compound: str = "", delivery: Optional[Dict[str, str]] = None, cond: Optional[Conditions] = None):
        super().__init__(parent)
        self.setWindowTitle("Amount calculator")
        self.resize(760, 420)
        self.cond = cond or Conditions()
        self.ch = ic.chamber(self.cond.chamber_name)
        self.result_ppb: Optional[float] = None
        self.result_delivery: Dict[str, str] = {}
        self.result_expl = ""
        d = dict(delivery or {})
        lay = QVBoxLayout(self)
        g = QGridLayout(); r = 0
        g.addWidget(QLabel("Compound"), r, 0); self.compound = QLineEdit(compound); g.addWidget(self.compound, r, 1, 1, 2)
        lb = QPushButton("Look up molar mass / density"); lb.clicked.connect(self._lookup); g.addWidget(lb, r, 3); r += 1
        g.addWidget(QLabel("Delivery"), r, 0); self.mode = QComboBox(); self.mode.addItems([m[1] for m in self.MODES]); g.addWidget(self.mode, r, 1, 1, 3); r += 1
        self.fields: Dict[str, QLineEdit] = {}

        def row(key, label, hint=""):
            nonlocal r
            g.addWidget(QLabel(label), r, 0); e = QLineEdit(d.get(key, "")); e.setPlaceholderText(hint); g.addWidget(e, r, 1)
            self.fields[key] = e; e.textChanged.connect(self._update); r += 1
        row("value", "Mixing ratio [ppb]")
        row("vol_ul", "Volume [µl]")
        row("mw", "Molar mass [g/mol]")
        row("density", "Density [g/ml]")
        row("vol_ml", "Gas volume [ml]")
        row("flow_sccm", "MFC flow [sccm]", "slpm × 1000")
        row("duration_s", "Duration [s]")
        row("bottle_value", "Bottle concentration", "1 for pure gas")
        g.addWidget(QLabel("Bottle unit"), r, 0); self.bottle_unit = QComboBox(); self.bottle_unit.addItems(["ppm", "ppb", "%", "pure"])
        self.bottle_unit.setCurrentText(d.get("bottle_unit", "ppm")); self.bottle_unit.currentTextChanged.connect(self._update); g.addWidget(self.bottle_unit, r, 1); r += 1
        lay.addLayout(g)
        lay.addWidget(_info(f"Chamber {self.ch.name}: {self.ch.volume_m3:g} m³ at {self.cond.temperature_K:g} K, {self.cond.pressure_hPa:g} hPa "
                            f"(from 'Where and when'); MFC standard {self.ch.mfc_std_T_K:g} K / {self.ch.mfc_std_P_hPa:g} hPa; "
                            f"gas syringe at {self.ch.syringe_T_K:g} K."))
        self.result = QPlainTextEdit(); self.result.setReadOnly(True); self.result.setMaximumHeight(90); lay.addWidget(self.result)
        bb = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        bb.button(QDialogButtonBox.StandardButton.Ok).setText("Use this ppb")
        bb.accepted.connect(self._accept); bb.rejected.connect(self.reject); lay.addWidget(bb)
        keys = [m[0] for m in self.MODES]
        self.mode.setCurrentIndex(keys.index(d["mode"]) if d.get("mode") in keys else 0)
        self.mode.currentIndexChanged.connect(self._update)
        if compound and not (d.get("mw") and d.get("density")):
            self._lookup(silent=True)
        self._update()

    def _lookup(self, silent: bool = False):
        rc = ic.resolve_compound(self.compound.text())
        if rc.mw:
            self.fields["mw"].setText(f"{rc.mw:g}")
        if rc.density:
            self.fields["density"].setText(f"{rc.density:g}")
        if not silent and not rc.mw:
            QMessageBox.information(self, "Look up", rc.note)
        self._update()

    def _delivery(self) -> Dict[str, str]:
        mode = self.MODES[self.mode.currentIndex()][0]
        d = {"mode": mode}
        for k, e in self.fields.items():
            if e.text().strip():
                d[k] = e.text().strip()
        d["bottle_unit"] = self.bottle_unit.currentText()
        return d

    def _update(self):
        mode = self.MODES[self.mode.currentIndex()][0]
        need = {"ppb": ["value"], "liquid": ["vol_ul", "mw", "density"], "gas_syringe": ["vol_ml", "bottle_value"],
                "mfc": ["flow_sccm", "duration_s", "bottle_value"]}[mode]
        for k, e in self.fields.items():
            e.setEnabled(k in need)
        self.bottle_unit.setEnabled(mode in ("gas_syringe", "mfc"))
        ppb, expl = ic.resolve_delivery(self._delivery(), self.cond.temperature_K, self.cond.pressure_hPa, self.ch)
        self.result_ppb, self.result_expl = ppb, expl
        self.result.setPlainText((f"= {ppb:.4g} ppb\n" if ppb is not None else "— \n") + expl)

    def _accept(self):
        if self.result_ppb is None:
            QMessageBox.warning(self, "Amount calculator", self.result_expl); return
        self.result_delivery = self._delivery()
        self.accept()


MIX_HEADERS = ["Component (as on the bottle)", "Amount", "Unit", "MCM species", "Surrogate (if no MCM)", "Molar mass", "Density", "Status"]


class MixtureDialog(QDialog):
    """A solution of several compounds injected as one pulse.  Recipe in, one injection
    event out: products = mole fractions of the MCM species (IDL VCP2 convention),
    total ppb = modelled sum.  Recipes are saved in mixtures/*.mixture.json."""

    def __init__(self, parent=None, cond: Optional[Conditions] = None):
        super().__init__(parent)
        self.setWindowTitle("Mixture injection")
        self.resize(1150, 760)
        self.cond = cond or Conditions()
        self.ch = ic.chamber(self.cond.chamber_name)
        self.result_event: Optional[InjectionEvent] = None
        self.result_precursors: List[str] = []
        lay = QVBoxLayout(self)
        g = QGridLayout()
        g.addWidget(QLabel("Mixture name"), 0, 0); self.name = QLineEdit("VCP"); g.addWidget(self.name, 0, 1)
        g.addWidget(QLabel("Aliases (words in the wiki line)"), 0, 2); self.aliases = QLineEdit("VCP solution"); g.addWidget(self.aliases, 0, 3)
        self.saved = QComboBox(); self.saved.addItem("— saved mixtures —")
        self._saved = ic.load_mixtures()
        for n in self._saved:
            self.saved.addItem(n)
        self.saved.currentTextChanged.connect(self._load_saved); g.addWidget(self.saved, 0, 4)
        lb = QPushButton("Load file …"); lb.clicked.connect(self._load_file); g.addWidget(lb, 0, 5)
        sb = QPushButton("Save recipe"); sb.clicked.connect(self._save); g.addWidget(sb, 0, 6)
        lay.addLayout(g)
        lay.addWidget(QLabel("Recipe — paste lines like '* Limonene 5.61 μl' or 'ethanol, 59.19, ul' (units: ul, ml, vol%, mg, g, wt%, mol%) and press Parse"))
        self.paste = QPlainTextEdit(); self.paste.setMaximumHeight(90); lay.addWidget(self.paste)
        hb = QHBoxLayout()
        pb = QPushButton("Parse recipe"); pb.clicked.connect(self._parse); hb.addWidget(pb)
        rb = QPushButton("Resolve in MCM"); rb.clicked.connect(self._resolve); hb.addWidget(rb)
        ab = QPushButton("+ row"); ab.clicked.connect(lambda: _add_row(self.table, ["", "", "ul", "", "", "", "", ""])); hb.addWidget(ab)
        db = QPushButton("− row"); db.clicked.connect(lambda: _del_row(self.table)); hb.addWidget(db)
        hb.addStretch(); lay.addLayout(hb)
        self.table = _table(MIX_HEADERS); lay.addWidget(self.table, stretch=2)
        self.table.itemChanged.connect(lambda _it: self._preview())
        self.table.itemSelectionChanged.connect(self._row_selected)
        # surrogate finder: formula or name → MCM species (same formula first, then same C/O count)
        sb = QGroupBox("Find a surrogate in the MCM for the selected component"); sl = QHBoxLayout(sb)
        self.sur_query = QLineEdit(); self.sur_query.setPlaceholderText("formula (e.g. C6H14O3) or part of an MCM name")
        self.sur_query.returnPressed.connect(self._find_surrogate); sl.addWidget(self.sur_query, stretch=1)
        fb2 = QPushButton("Search"); fb2.clicked.connect(self._find_surrogate); sl.addWidget(fb2)
        self.sur_list = QListWidget(); self.sur_list.setMinimumWidth(420); self.sur_list.setMaximumHeight(110)
        self.sur_list.itemDoubleClicked.connect(self._use_surrogate); sl.addWidget(self.sur_list, stretch=2)
        ub = QPushButton("Use as surrogate"); ub.clicked.connect(lambda: self._use_surrogate(self.sur_list.currentItem())); sl.addWidget(ub)
        lay.addWidget(sb)
        lay.addWidget(_info("Every MCM species has a formula and a SMILES; the list shows exact formula matches first, then species "
                            "with the same carbon and oxygen count, then name matches; radicals last. Double-click to take one."))
        g2 = QGridLayout()
        g2.addWidget(QLabel("Injected volume [µl]"), 0, 0); self.injected = QLineEdit(""); self.injected.setPlaceholderText("empty = recipe amounts as is")
        self.injected.textChanged.connect(self._preview); g2.addWidget(self.injected, 0, 1)
        g2.addWidget(QLabel("Start HH:MM"), 0, 2); self.start = QLineEdit("08:00"); g2.addWidget(self.start, 0, 3)
        g2.addWidget(QLabel("Duration [min]"), 0, 4); self.dur = QLineEdit("2"); g2.addWidget(self.dur, 0, 5)
        lay.addLayout(g2)
        lay.addWidget(QLabel("Result — alerts (⚠ needs a decision, ✖ blocks)"))
        self.preview = QPlainTextEdit(); self.preview.setReadOnly(True); lay.addWidget(self.preview, stretch=1)
        self.save_box = QCheckBox("Save the recipe (with the surrogates) to the mixtures folder for reuse by the wiki import")
        self.save_box.setChecked(True); lay.addWidget(self.save_box)
        bb = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        bb.button(QDialogButtonBox.StandardButton.Ok).setText("Add injection")
        bb.accepted.connect(self._accept); bb.rejected.connect(self.reject); lay.addWidget(bb)

    def preload(self, m: ic.Mixture, injected_ul: str = "", start: str = "", dur: str = ""):
        """Show an existing recipe (e.g. the one a wiki import used) with its surrogates editable."""
        self.name.setText(m.name); self.aliases.setText(", ".join(m.aliases))
        ic.resolve_mixture(m)
        for c in m.components:                       # fill suggested surrogates as notes, never silently as choices
            if not c.mcm and not c.surrogate:
                rc = ic.resolve_compound(c.name)
                c.note = rc.note
        self._set_components(m.components)
        if injected_ul:
            self.injected.setText(injected_ul)
        if start:
            self.start.setText(start)
        if dur:
            self.dur.setText(dur)

    # ── recipe table ↔ Mixture ──
    def _set_components(self, comps: List[ic.MixtureComponent]):
        self.table.blockSignals(True)
        _set_rows(self.table, [[c.name, f"{c.amount:g}", c.unit, c.mcm, c.surrogate,
                                f"{c.mw:g}" if c.mw else "", f"{c.density:g}" if c.density else "", c.note] for c in comps])
        self.table.blockSignals(False)
        self._preview()

    def _mixture(self) -> ic.Mixture:
        comps = []
        for name, amt, unit, mcm, sur, mw, dens, _note in _rows(self.table):
            if not name:
                continue
            try:
                comps.append(ic.MixtureComponent(name, float(amt or 0), (unit or "ul").lower(), mcm.upper(), sur.upper(),
                                                 float(mw) if mw else None, float(dens) if dens else None))
            except ValueError:
                comps.append(ic.MixtureComponent(name, 0.0, unit or "ul", mcm.upper(), sur.upper(), None, None, "bad number"))
        return ic.Mixture(self.name.text().strip() or "mixture", aliases=[a.strip() for a in self.aliases.text().split(",") if a.strip()],
                          components=comps)

    def _row_selected(self):
        r = self.table.currentRow()
        if r < 0 or self.table.item(r, 3) is None or self.table.item(r, 3).text().strip():
            return                                  # has an MCM species already
        rc = ic.resolve_compound(self.table.item(r, 0).text())
        if rc.formula:
            self.sur_query.setText(rc.formula); self._find_surrogate()

    def _find_surrogate(self):
        self.sur_list.clear()
        q = self.sur_query.text().strip()
        if not q:
            return
        try:
            hits = search_species(q)
        except Exception as e:
            self.sur_list.addItem(f"search failed: {e}"); return
        for sp, why in hits:
            it = QListWidgetItem(f"{sp.name:14s} {sp.formula or '':12s} {sp.mass or 0:7.1f} g/mol   {why}   {sp.smiles or ''}")
            it.setData(Qt.ItemDataRole.UserRole, sp.name); self.sur_list.addItem(it)
        if not hits:
            self.sur_list.addItem(f"nothing in the MCM for '{q}' — try the formula with one oxygen more or less, or a shorter carbon chain")

    def _use_surrogate(self, item):
        r = self.table.currentRow()
        if item is None or r < 0:
            return
        name = item.data(Qt.ItemDataRole.UserRole)
        if not name:
            return
        self.table.setItem(r, 4, QTableWidgetItem(name))

    def _parse(self):
        self._set_components(ic.parse_mixture_text(self.paste.toPlainText()))

    def _resolve(self):
        m = self._mixture(); ic.resolve_mixture(m); self._set_components(m.components)

    def _load_saved(self, name: str):
        m = self._saved.get(name)
        if m:
            self.name.setText(m.name); self.aliases.setText(", ".join(m.aliases)); self._set_components(m.components)

    def _load_file(self):
        p, _ = QFileDialog.getOpenFileName(self, "Mixture file", ic.MIXTURE_DIR, "Mixture (*.mixture.json);;All files (*)")
        if p:
            try:
                m = ic.load_mixture(p)
            except Exception as e:
                QMessageBox.critical(self, "Load failed", str(e)); return
            self.name.setText(m.name); self.aliases.setText(", ".join(m.aliases)); self._set_components(m.components)

    def _save(self):
        m = self._mixture()
        try:
            path = ic.save_mixture(m)
        except Exception as e:
            QMessageBox.critical(self, "Save failed", str(e)); return
        QMessageBox.information(self, "Mixture saved", f"{path}\nThe wiki import will now recognise '{m.name}' and its aliases.")

    def _compute(self):
        m = self._mixture()
        inj = self.injected.text().strip()
        ev, precs, alerts = ic.mixture_event(m, float(inj) if inj else None, self.start.text().strip(),
                                             float(self.dur.text() or 1), self.cond.temperature_K, self.cond.pressure_hPa, self.ch)
        return ev, precs, alerts

    def _preview(self):
        try:
            ev, precs, alerts = self._compute()
        except Exception as e:
            self.preview.setPlainText(f"cannot compute: {e}"); return
        lines = [ev.delivery.get("calc", ""), "", f"products: {ev.products or '(none)'}", f"total modelled: {ev.total_ppb} ppb",
                 "precursors for the MCM subset: " + (", ".join(precs) or "none")]
        for lvl, msg in alerts:
            lines.append(("✖ " if lvl == "error" else "⚠ ") + msg)
        self.preview.setPlainText("\n".join(lines))

    def _accept(self):
        try:
            ev, precs, alerts = self._compute()
        except Exception as e:
            QMessageBox.critical(self, "Mixture", str(e)); return
        if any(l == "error" for l, _ in alerts) or not ev.products:
            QMessageBox.warning(self, "Mixture", "Fix the ✖ alerts first (molar mass/density, MCM names, surrogates)."); return
        if self.save_box.isChecked():
            try:
                ic.save_mixture(self._mixture())
            except Exception as e:
                QMessageBox.warning(self, "Mixture", f"Recipe not saved: {e}")
        self.result_event, self.result_precursors = ev, precs
        self.accept()


# ─────────────────────────────────────────────────────────────────────────────
class _WikiFetchThread(QThread):
    done = pyqtSignal(str)
    failed = pyqtSignal(str)

    def __init__(self, date: str, user: str, password: str):
        super().__init__(); self.date, self.user, self.password = date, user, password

    def run(self):
        try:
            self.done.emit(saphir_wiki.fetch_wiki_text(self.date, self.user, self.password))
        except Exception as e:
            self.failed.emit(str(e))


def _wrap(text: str) -> QLabel:
    lb = QLabel(text); lb.setWordWrap(True); return lb


class WikiImportDialog(QDialog):
    """Fetch / load / paste a saphirwiki day page, show the parsed action log for confirmation,
    and hand back a WikiImport (roof, injections with actions, precursors, checks)."""

    def __init__(self, parent=None, date: str = "", archive_root: str = "", T_K: float = 298.0, P_hPa: float = 1013.0,
                 chamber_name: str = "SAPHIR"):
        super().__init__(parent)
        self.T_K, self.P_hPa, self.chamber_name = T_K, P_hPa, chamber_name
        self.setWindowTitle("Import experiment log — SAPHIR wiki, saved page or pasted notes")
        self.resize(1100, 720)
        self.result_import: Optional[saphir_wiki.WikiImport] = None
        self.archive_root = archive_root
        self._events: List[saphir_wiki.WikiEvent] = []
        lay = QVBoxLayout(self)
        g = QGridLayout()
        g.addWidget(QLabel("Date (YYYY-MM-DD)"), 0, 0); self.date = QLineEdit(date); g.addWidget(self.date, 0, 1)
        env_user, env_pw = saphir_wiki.wiki_credentials()
        g.addWidget(QLabel("Wiki user"), 0, 2); self.user = QLineEdit(); g.addWidget(self.user, 0, 3)
        g.addWidget(QLabel("Password"), 0, 4); self.pw = QLineEdit(); self.pw.setEchoMode(QLineEdit.EchoMode.Password); g.addWidget(self.pw, 0, 5)
        if env_user:
            self.user.setPlaceholderText(f"{env_user} (from {saphir_wiki.ENV_USER})")
            self.pw.setPlaceholderText("from " + saphir_wiki.ENV_PASS if env_pw else "not set: type it or set " + saphir_wiki.ENV_PASS)
        else:
            self.user.setPlaceholderText(f"or set {saphir_wiki.ENV_USER}")
            self.pw.setPlaceholderText(f"or set {saphir_wiki.ENV_PASS}")
        self.fetch_btn = QPushButton("Fetch from wiki"); self.fetch_btn.clicked.connect(self._fetch); g.addWidget(self.fetch_btn, 0, 6)
        lb = QPushButton("Load saved page / text file …"); lb.clicked.connect(self._load_file); g.addWidget(lb, 1, 6)
        g.addWidget(_info("Needs the FZJ network and a wiki login (the page answers 'not allowed to view' otherwise). "
                          "Without a login: open the page in the browser, save it (or copy the action log) and load/paste it here."),
                    1, 0, 1, 6)
        lay.addLayout(g)
        self.status = QLabel(""); lay.addWidget(self.status)
        # left: the original page text (fetched, loaded or pasted) — right: what was made of it
        split = QSplitter(Qt.Orientation.Horizontal); lay.addWidget(split, stretch=1)
        left = QWidget(); ll = QVBoxLayout(left); ll.setContentsMargins(0, 0, 0, 0)
        ll.addWidget(_wrap("<b>Original wiki text</b> — as fetched / loaded, or paste an action log here; select a row on the right to find its line"))
        self.paste = QPlainTextEdit()
        self.paste.setPlaceholderText("07:30 H F.R. CHANGE: humidification start\n08:15 h F.R. CHANGE: humidification stop\n"
                                      "08:40 X M.S. CHANGE: 20 µl isoprene injected\n09:30 R F.R. CHANGE: open roof")
        ll.addWidget(self.paste, stretch=1)
        ll.addWidget(_more_info("Each line: time, one-letter code, free text — e.g. 07:30 H F.R. CHANGE: humidification start",
                                WIKI_LINE_HELP, "ⓘ Line format"))
        pb = QPushButton("Parse this text"); pb.clicked.connect(self._parse_paste)
        hb = QHBoxLayout(); hb.addWidget(pb); hb.addStretch(); ll.addLayout(hb)
        split.addWidget(left)
        right = QWidget(); rl = QVBoxLayout(right); rl.setContentsMargins(0, 0, 0, 0)
        rl.addWidget(_wrap("<b>Action log</b> — what the parser read; untick lines that should not be used"))
        self.table = _table(["Use", "Time", "Code", "Action", "Wiki text"])
        for c in (0, 1, 2):
            self.table.horizontalHeader().setSectionResizeMode(c, QHeaderView.ResizeMode.ResizeToContents)
        self.table.itemChanged.connect(lambda _it: self._preview())
        self.table.currentCellChanged.connect(lambda r, _c, _pr, _pc: self._show_source_line(r))
        rl.addWidget(self.table, stretch=2)
        rl.addWidget(_wrap("<b>What will be imported</b> — ⚠ marks things the wiki did not say clearly"))
        self.preview = QPlainTextEdit(); self.preview.setReadOnly(True); rl.addWidget(self.preview, stretch=1)
        split.addWidget(right); split.setSizes([420, 680])
        hbr = QHBoxLayout()
        self.replace_box = QCheckBox("Replace the current injection rows (untick to append)"); self.replace_box.setChecked(True)
        hbr.addWidget(self.replace_box); hbr.addStretch()
        self.sur_btn = QPushButton("Set surrogates for the mixture …"); self.sur_btn.setEnabled(False); self.sur_btn.clicked.connect(self._surrogates)
        hbr.addWidget(self.sur_btn); lay.addLayout(hbr)
        bb = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        bb.button(QDialogButtonBox.StandardButton.Ok).setText("Use in wizard")
        bb.accepted.connect(self._accept); bb.rejected.connect(self.reject); lay.addWidget(bb)

    # ── sources ──
    def _fetch(self):
        d = self.date.text().strip()
        if not re.fullmatch(r"\d{4}-\d{2}-\d{2}", d):
            QMessageBox.warning(self, "Date", "Enter the date as YYYY-MM-DD."); return
        self.fetch_btn.setEnabled(False); self.status.setStyleSheet(""); self.status.setText("Fetching …")
        self._thr = _WikiFetchThread(d, self.user.text().strip(), self.pw.text())
        self._thr.done.connect(self._fetched); self._thr.failed.connect(self._fetch_failed); self._thr.start()

    def _fetched(self, text: str):
        self.fetch_btn.setEnabled(True)
        self.paste.setPlainText(text)                      # keep the original next to the parsed log
        self._set_events(saphir_wiki.parse_action_log(text), f"wiki {self.date.text().strip()}")

    def _fetch_failed(self, msg: str):
        self.fetch_btn.setEnabled(True); self.status.setText("Fetch failed: " + msg[:160])
        self.status.setStyleSheet("color: #e05252;")

    def _load_file(self):
        p, _ = QFileDialog.getOpenFileName(self, "Saved wiki page", "", "Wiki page (*.html *.htm *.txt *.md);;All files (*)")
        if not p:
            return
        try:
            text = saphir_wiki.load_wiki_file(p)
        except Exception as e:
            QMessageBox.critical(self, "Load failed", str(e)); return
        m = re.search(r"(\d{4}-\d{2}-\d{2})", os.path.basename(p) + " " + text[:2000])
        if m and not self.date.text().strip():
            self.date.setText(m.group(1))
        self.paste.setPlainText(text)
        self._set_events(saphir_wiki.parse_action_log(text), os.path.basename(p))

    def _parse_paste(self):
        self._set_events(saphir_wiki.parse_action_log(self.paste.toPlainText()), "pasted text")

    def memory(self) -> dict:
        """What the dialog holds, to reopen it as it was within the same experiment (no password)."""
        return {"date": self.date.text().strip(), "user": self.user.text().strip(), "text": self.paste.toPlainText(),
                "source": getattr(self, "_source", "pasted text"),
                "use": [self.table.item(r, 0).checkState() == Qt.CheckState.Checked for r in range(self.table.rowCount())],
                "replace": self.replace_box.isChecked()}

    def restore(self, mem: dict):
        """Reopen with the text, ticks and options of the last import of this experiment."""
        if not mem:
            return
        if mem.get("date"):
            self.date.setText(mem["date"])
        if mem.get("user"):
            self.user.setText(mem["user"])
        self.replace_box.setChecked(bool(mem.get("replace", True)))
        text = mem.get("text") or ""
        if not text.strip():
            return
        self.paste.setPlainText(text)
        self._set_events(saphir_wiki.parse_action_log(text), mem.get("source", "pasted text"))
        use = mem.get("use") or []
        if len(use) == self.table.rowCount():
            self.table.blockSignals(True)
            for r, on in enumerate(use):
                self.table.item(r, 0).setCheckState(Qt.CheckState.Checked if on else Qt.CheckState.Unchecked)
            self.table.blockSignals(False)
            self._preview()
        self.status.setText(self.status.text() + "  — restored from the last import of this experiment")

    def _set_events(self, events, source: str):
        self._source = source
        self._events = events
        self.status.setStyleSheet("color: #4caf50;" if events else "color: #e05252;")
        guessed = sum(1 for e in events if e.description and not e.line.split(None, 2)[1:2] == [e.code])
        self.status.setText((f"{len(events)} action-log lines from {source}" + (" — codes read from the words where the log has none; check the 'Code' column" if guessed else ""))
                            if events else f"No time-stamped lines ('HH:MM …') found in {source}")
        self.table.blockSignals(True); self.table.setRowCount(0)
        for ev in events:
            r = self.table.rowCount(); self.table.insertRow(r)
            use = QTableWidgetItem(""); use.setFlags(Qt.ItemFlag.ItemIsUserCheckable | Qt.ItemFlag.ItemIsEnabled)
            use.setCheckState(Qt.CheckState.Checked if ev.action != "note" else Qt.CheckState.Unchecked)
            self.table.setItem(r, 0, use)
            for c, v in enumerate([ev.time, ev.code, ev.action, ev.description], start=1):
                it = QTableWidgetItem(v); it.setFlags(it.flags() & ~Qt.ItemFlag.ItemIsEditable); self.table.setItem(r, c, it)
        self.table.blockSignals(False)
        self._preview()

    def _show_source_line(self, row: int):
        """Highlight and scroll to the original line of the selected action-log row."""
        from PyQt6.QtGui import QTextCursor, QTextCharFormat, QColor
        from PyQt6.QtWidgets import QTextEdit as _QTE            # ExtraSelection lives on QTextEdit in PyQt6
        try:
            self._show_source_line_impl(row, QTextCursor, QTextCharFormat, QColor, _QTE)
        except Exception as e:                                   # never let a slot raise: PyQt6 would abort the process
            print(f"[wiki] highlight failed: {e}")

    def _show_source_line_impl(self, row, QTextCursor, QTextCharFormat, QColor, _QTE):
        if row < 0 or row >= len(self._events):
            return
        line = (self._events[row].line or "").strip()
        doc = self.paste.document(); text = self.paste.toPlainText()
        pos = text.find(line) if line else -1
        if pos < 0:
            needle = f"{self._events[row].time} {self._events[row].code}"
            pos = text.find(needle)
        if pos < 0:
            self.paste.setExtraSelections([]); return
        cur = QTextCursor(doc); cur.setPosition(pos)
        cur.movePosition(QTextCursor.MoveOperation.EndOfLine, QTextCursor.MoveMode.KeepAnchor)
        sel = _QTE.ExtraSelection(); fmt = QTextCharFormat(); fmt.setBackground(QColor("#fff2a8")); sel.format = fmt; sel.cursor = cur
        self.paste.setExtraSelections([sel])
        show = QTextCursor(doc); show.setPosition(pos); self.paste.setTextCursor(show); self.paste.centerCursor()

    def _current_import(self) -> saphir_wiki.WikiImport:
        for r, ev in enumerate(self._events):
            it = self.table.item(r, 0)
            ev.use = bool(it and it.checkState() == Qt.CheckState.Checked)
        return saphir_wiki.events_to_import(self._events, self.date.text().strip(), self.T_K, self.P_hPa, self.chamber_name)

    def _preview(self):
        if not self._events:
            self.preview.setPlainText(""); return
        imp = self._current_import()
        lines = [f"Date {imp.date or '?'};  roof open: " + (", ".join(f"{o}–{c}" for o, c in imp.roof_open) or "never")]
        day = saphir_wiki.archive_day_folder(self.archive_root, imp.date) if imp.date else ""
        lines.append("Archive: " + (f"found {day} → measured T/P/H2O/J/NO/NO2/O3 will be used" if day
                                    else "day folder not found → manual conditions, clear-sky J"))
        lines.append("Precursors (MCM): " + (", ".join(imp.precursors) or "none"))
        for inj in imp.injections:
            what = f"{inj.total_ppb or '?'} ppb → {inj.products[:60]}" if inj.products else "recorded only"
            mode = inj.delivery.get("mode", "")
            lines.append(f"  {ACTIONS[inj.action]['label']:42s} {inj.name:10s} {inj.start} +{inj.duration_min:g} min   {what}"
                         + (f"   [{mode}]" if mode else ""))
        for lvl, msg in imp.checks:
            lines.append(("⚠ " if lvl == "warning" else "✖ ") + msg)
        if imp.notes:
            lines.append("Notes kept out of the model: " + "; ".join(n[:60] for n in imp.notes[:8]) + (" …" if len(imp.notes) > 8 else ""))
        self._mix_names = sorted({inj.delivery.get("mixture", "") for inj in imp.injections if inj.delivery.get("mode") == "mixture"} - {""})
        self.sur_btn.setEnabled(bool(self._mix_names))
        if self._mix_names:
            lines.append(f"Mixture recipe(s) used: {', '.join(self._mix_names)} — 'Set surrogates …' opens the recipe to map "
                         "the components the MCM lacks.")
        self.preview.setPlainText("\n".join(lines))

    def _surrogates(self):
        """Open the recipe of the mixture the import used, let the user map components to
        MCM surrogates, save it, and rebuild the preview with it."""
        names = getattr(self, "_mix_names", [])
        if not names:
            return
        mixes = ic.load_mixtures()
        m = mixes.get(names[0])
        if m is None:
            QMessageBox.information(self, "Mixture", f"Recipe '{names[0]}' not found in {ic.MIXTURE_DIR}"); return
        dlg = MixtureDialog(self, cond=None)
        dlg.preload(m)
        dlg.injected.setText(""); dlg.save_box.setChecked(True)
        if dlg.exec() == QDialog.DialogCode.Accepted:
            self._preview()

    def _accept(self):
        if not self._events:
            QMessageBox.information(self, "Wiki import", "Nothing parsed yet — fetch, load or paste the action log first."); return
        self.result_import = self._current_import()
        self.accept()


# ─────────────────────────────────────────────────────────────────────────────
class ExperimentHost:
    """Everything the pages need from their container: the page objects, include/easy
    roots, spec <-> pages mapping, wiki import, output folder and writing.  Mixed into
    the QWizard (ExperimentWizard) and into the step-by-step main window (juice_gui)."""

    include_dir: str = ""
    easy_root: str = ""

    def _init_host(self, include_dir: str, easy_root: str):
        self.include_dir, self.easy_root = include_dir, easy_root or os.getcwd()
        self.wiki_memory: Optional[dict] = None       # last import dialog state of the open experiment (cleared by new / open / template)
        self._base = ExperimentSpec()
        self.page_start = StartPage(self); self.page_cond = ConditionsPage(self); self.page_chem = ChemistryPage(self)
        self.page_inj = InjectionsPage(self); self.page_chamber = ChamberPage(self); self.page_init = InitialsPage(self)
        self.page_chamber.embed_initials(self.page_init)
        self.page_review = ReviewPage(self)
        self.last_result: dict = {}
        self.import_notes: List[str] = []                 # kept with the experiment (ExperimentSpec.import_notes)

    def pages(self):
        """Pages shown one after the other (the initial values sit inside the chamber page)."""
        return (self.page_start, self.page_cond, self.page_chem, self.page_inj, self.page_chamber, self.page_review)

    # ── spec ↔ pages ─────────────────────────────────────────────────────────
    def load_spec(self, spec: ExperimentSpec, keep_name: bool = False):
        self._base = spec
        if not keep_name:
            self.page_start.name_edit.setText(spec.name)
        self.page_cond.set_from(spec.conditions)
        self.page_chem.set_from(spec)
        self.page_inj.set_from(spec)
        self.page_chamber.set_from(spec.chamber)
        self.page_init.set_from(spec)
        self.page_review.overrides.setPlainText("\n".join(spec.overrides))
        self.import_notes = list(spec.import_notes)

    def chamber_selected(self, profile_id: str):
        """The user chose another chamber on 'Where and when': the chamber page takes that profile's
        default switches and DILUTE (strengths and formulas typed so far are kept)."""
        from experiment_spec import chamber_defaults
        self.page_chamber.set_from(chamber_defaults(profile_id, self.page_chamber.get(self._base.chamber)))

    def collect_spec(self) -> ExperimentSpec:
        b = self._base
        spec = ExperimentSpec(name=self.page_start.name_edit.text().strip() or "experiment",
                              title=b.title, mech_incl_suffix=b.mech_incl_suffix, use_anorg=b.use_anorg,
                              use_j_clear_sky_incl=b.use_j_clear_sky_incl, hmax=b.hmax, outstep=b.outstep,
                              extra_includes=list(b.extra_includes), template=b.template)
        spec.conditions = self.page_cond.get()
        self.page_chem.apply(spec)
        spec.injections = self.page_inj.get()
        spec.chamber = self.page_chamber.get(b.chamber)
        self.page_init.apply(spec)
        spec.overrides = [l.strip() for l in self.page_review.overrides.toPlainText().splitlines() if l.strip()]
        spec.import_notes = list(self.import_notes)
        return spec

    # ── wiki import ───────────────────────────────────────────────────────────
    def import_from_wiki(self, parent=None, paste: bool = False):
        cond = self.page_cond.get()
        dlg = WikiImportDialog(parent or self, date=cond.date, archive_root=cond.archive_root,
                               T_K=cond.temperature_K, P_hPa=cond.pressure_hPa, chamber_name=cond.chamber_name)
        dlg.restore(getattr(self, "wiki_memory", None))
        if paste:
            dlg.paste.setFocus()
        accepted = dlg.exec() == QDialog.DialogCode.Accepted
        self.wiki_memory = dlg.memory()                # kept also when cancelled: the text and ticks are work in progress
        if not accepted or dlg.result_import is None:
            return
        self.apply_wiki_import(dlg.result_import, replace=dlg.replace_box.isChecked())

    def apply_wiki_import(self, imp: "saphir_wiki.WikiImport", replace: bool = True):
        """Put a WikiImport into the pages: date + roof + archive mode (Conditions), precursors
        (Chemistry), injections (Injections).  The name defaults to the date if still the placeholder."""
        cond = self.page_cond.get()
        if imp.date:
            cond.date = imp.date
        if imp.roof_open:
            cond.roof_open = list(imp.roof_open)
        day = saphir_wiki.archive_day_folder(cond.archive_root, cond.date)
        if day:
            cond.source, cond.j_source = "archive", "measured"
        else:
            cond.source, cond.j_source = "manual", "clear_sky"
        self.page_cond.set_from(cond)
        chem = self.page_chem
        if imp.precursors:
            chem.mode_compounds.setChecked(True)
            have = {chem.chosen.item(i).text() for i in range(chem.chosen.count())}
            for p in imp.precursors:
                if p not in have:
                    chem.chosen.addItem(p)
        notes = [m for _lvl, m in imp.checks]
        self.import_notes = notes if replace else self.import_notes + notes
        spec = self.collect_spec()
        spec.injections = list(imp.injections) if replace else spec.injections + list(imp.injections)
        self.page_inj.set_from(spec)
        if imp.date and self.page_start.name_edit.text().strip() in ("", "my_experiment", "experiment"):
            self.page_start.name_edit.setText("saphir_" + imp.date.replace("-", ""))
        self.last_wiki_import = imp

    def out_dir(self) -> str:
        return os.path.join(self.page_start.out_edit.text().strip() or self.easy_root,
                            self.page_start.name_edit.text().strip() or "experiment")

    def write_all(self) -> dict:
        """Write <name>.enz, <name>.easy, chemistry subset and <name>.experiment.json into
        out_dir().  Returns write_experiment's result dict (+ 'spec', 'spec_path').  Raises on error."""
        spec = self.collect_spec()
        out = self.out_dir()
        res = write_experiment(spec, out)
        res["spec_path"] = save_spec(spec, os.path.join(out, spec.name + ".experiment.json"))
        res["spec"] = spec
        self.last_result = res
        for n in res.get("notes", []):
            print(f"[experiment] {n}")
        return res


# ─────────────────────────────────────────────────────────────────────────────
class ExperimentWizard(QWizard, ExperimentHost):
    experiment_written = pyqtSignal(str)

    def __init__(self, parent=None, include_dir: str = "", easy_root: str = "", last_easy_dir: str = ""):
        super().__init__(parent)
        self.setWindowTitle("New experiment — JUICE")
        self.resize(1250, 850)
        self.setWizardStyle(QWizard.WizardStyle.ModernStyle)
        self._init_host(include_dir, easy_root or last_easy_dir)
        for p in self.pages():
            self.addPage(p)
        self.load_spec(self._base, keep_name=True)

    # ── finish ───────────────────────────────────────────────────────────────
    def accept(self):
        try:
            res = self.write_all()
        except Exception as e:
            QMessageBox.critical(self, "Write failed", str(e)); return
        # no modal dialog here so that headless use cannot block
        super().accept()
        self.experiment_written.emit(res["easy"])


if __name__ == "__main__":
    import sys
    from PyQt6.QtWidgets import QApplication
    app = QApplication(sys.argv)
    w = ExperimentWizard(include_dir=sys.argv[1] if len(sys.argv) > 1 else "", easy_root=sys.argv[2] if len(sys.argv) > 2 else "")
    w.show()
    sys.exit(app.exec())
