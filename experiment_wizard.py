"""
experiment_wizard.py — guided scenario builder (Step 3 of the JUICE workflow)
============================================================================

A QWizard that asks experiment questions and fills an `ExperimentSpec`
(experiment_spec.py).  On Finish the backend writes <name>.enz, <name>.easy,
the mechanism subset (if compounds were chosen) and <name>.experiment.json,
then hands the .easy to FR_GUI (Open Scenario path).  The user never edits a
file; the review page shows the generated EASY and keeps an "expert overrides"
box and a button to the block-level editor (easy_builder.EasyBuilderDialog).

Pages:  Start (template, name, folder) → Where and when → Chemistry →
        Injections → Chamber effects → Initial conditions → Review
"""

from __future__ import annotations

import os
import re
from typing import Dict, List, Optional, Tuple

from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtWidgets import (QWizard, QWizardPage, QVBoxLayout, QHBoxLayout, QGridLayout, QLabel,
                             QLineEdit, QPushButton, QCheckBox, QComboBox, QListWidget, QListWidgetItem,
                             QTableWidget, QTableWidgetItem, QPlainTextEdit, QFileDialog, QMessageBox,
                             QGroupBox, QHeaderView, QRadioButton, QButtonGroup, QWidget, QSpinBox,
                             QAbstractItemView)

from experiment_spec import (ExperimentSpec, Conditions, InjectionEvent, ChamberEffects,
                             list_templates, load_template, load_spec, save_spec, validate,
                             mechanism_species, render_easy, easy_form, write_experiment,
                             H2O_IMPURITY_DEFAULTS, STANDARD_COLUMNS)
from easy_builder import scan_mechanisms, EasyBuilderDialog

try:
    from mcm_species_db import MCM_DB
except Exception:                               # database optional
    MCM_DB = {}


def _info(text: str) -> QLabel:
    lab = QLabel(text)
    lab.setWordWrap(True)
    lab.setStyleSheet("color: #888;")
    return lab


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
    def __init__(self, wiz: "ExperimentWizard"):
        super().__init__()
        self.wiz = wiz
        self.setTitle("Start")
        self.setSubTitle("Pick a template that resembles your experiment, give it a name and choose where the files go.")
        lay = QVBoxLayout(self)
        self.tpl_list = QListWidget()
        self.templates = list_templates()
        self.tpl_list.addItem(QListWidgetItem("Empty — start from the SAPHIR defaults"))
        for name, path in self.templates.items():
            try:
                title, desc, _ = load_template(path)
            except Exception as e:
                title, desc = name, f"(unreadable: {e})"
            it = QListWidgetItem(f"{title}\n    {desc}")
            it.setData(Qt.ItemDataRole.UserRole, path)
            self.tpl_list.addItem(it)
        self.tpl_list.setCurrentRow(0)
        self.tpl_list.itemSelectionChanged.connect(self._template_chosen)
        lay.addWidget(QLabel("Template"))
        lay.addWidget(self.tpl_list)
        g = QGridLayout()
        g.addWidget(QLabel("Experiment name"), 0, 0)
        self.name_edit = QLineEdit("my_experiment"); g.addWidget(self.name_edit, 0, 1, 1, 2)
        g.addWidget(QLabel("Project folder"), 1, 0)
        self.out_edit = QLineEdit(wiz.easy_root); g.addWidget(self.out_edit, 1, 1)
        b = QPushButton("…"); b.clicked.connect(self._browse); g.addWidget(b, 1, 2)
        lay.addLayout(g)
        lay.addWidget(_info("A sub-folder with the experiment name is created inside the project folder. "
                            "It will hold the generated .enz, .easy, chemistry and the experiment description "
                            "(.experiment.json) that you can reopen here later."))
        hb = QHBoxLayout()
        ob = QPushButton("Open a saved experiment (.experiment.json) …"); ob.clicked.connect(self._open_saved)
        hb.addWidget(ob); hb.addStretch(); lay.addLayout(hb)

    def _browse(self):
        d = QFileDialog.getExistingDirectory(self, "Project folder", self.out_edit.text() or self.wiz.easy_root)
        if d:
            self.out_edit.setText(d)

    def _template_chosen(self):
        it = self.tpl_list.currentItem()
        path = it.data(Qt.ItemDataRole.UserRole) if it else None
        if path:
            try:
                _t, _d, spec = load_template(path)
            except Exception as e:
                QMessageBox.warning(self, "Template", str(e)); return
        else:
            spec = ExperimentSpec()
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
        self.setSubTitle("Date and time window of the experiment, chamber conditions and the roof schedule. "
                         "These become the .enz environment file.")
        g = QGridLayout(self)
        r = 0
        g.addWidget(QLabel("Date (UTC)"), r, 0); self.date = QLineEdit("2024-07-09"); g.addWidget(self.date, r, 1)
        g.addWidget(QLabel("from"), r, 2); self.start = QLineEdit("04:00"); g.addWidget(self.start, r, 3)
        g.addWidget(QLabel("to"), r, 4); self.end = QLineEdit("21:59"); g.addWidget(self.end, r, 5)
        g.addWidget(QLabel("step [s]"), r, 6); self.step = QSpinBox(); self.step.setRange(10, 3600); self.step.setValue(60); g.addWidget(self.step, r, 7)
        r += 1
        self.src_manual = QRadioButton("Manual conditions (constants below)")
        self.src_archive = QRadioButton("SAPHIR archive (NetCDF) — measured T, flow, humidity, J values, NO/NO2/O3")
        self.src_manual.setChecked(True)
        g.addWidget(self.src_manual, r, 0, 1, 4); g.addWidget(self.src_archive, r, 4, 1, 4); r += 1
        g.addWidget(QLabel("Archive root"), r, 0); self.archive = QLineEdit(r"Z:\IEK8-SAPHIR\data"); g.addWidget(self.archive, r, 1, 1, 7); r += 1
        g.addWidget(QLabel("Temperature [K]"), r, 0); self.temp = QLineEdit("298"); g.addWidget(self.temp, r, 1)
        self.temp_mode = QComboBox(); self.temp_mode.addItems(["follow ENZ column (T = T_in)", "hold constant in the model"])
        g.addWidget(self.temp_mode, r, 2, 1, 2)
        g.addWidget(QLabel("Pressure [hPa]"), r, 4); self.press = QLineEdit("1006.5"); g.addWidget(self.press, r, 5); r += 1
        g.addWidget(QLabel("Water vapour [vol-%]"), r, 0); self.h2o = QLineEdit("1.0"); g.addWidget(self.h2o, r, 1)
        g.addWidget(QLabel("Replenishment flow [m³/h]"), r, 2); self.flow = QLineEdit("6.0"); g.addWidget(self.flow, r, 3)
        g.addWidget(QLabel("Photolysis frequencies"), r, 4)
        self.jsrc = QComboBox(); self.jsrc.addItems(["measured (archive J values)", "clear sky from solar zenith angle"])
        g.addWidget(self.jsrc, r, 5, 1, 3); r += 1
        g.addWidget(QLabel("Roof open / close (HH:MM)"), r, 0)
        self.roof = _table(["open", "close"]); self.roof.setMaximumHeight(110); g.addWidget(self.roof, r, 1, 1, 5)
        bb = QVBoxLayout(); a = QPushButton("+"); a.clicked.connect(lambda: _add_row(self.roof, ["09:30", "18:30"])); d = QPushButton("−"); d.clicked.connect(lambda: _del_row(self.roof)); bb.addWidget(a); bb.addWidget(d); g.addLayout(bb, r, 6); r += 1
        g.addWidget(QLabel("Hold species to measured values"), r, 0)
        self.constrain = QLineEdit(""); self.constrain.setPlaceholderText("e.g. NO2, O3  (archive only; the model then does not predict them)")
        g.addWidget(self.constrain, r, 1, 1, 7); r += 1
        g.addWidget(_info("Without archive access choose 'clear sky' photolysis; measured J values need the archive."), r, 0, 1, 8)

    def set_from(self, c: Conditions):
        self.date.setText(c.date); self.start.setText(c.start); self.end.setText(c.end); self.step.setValue(int(c.step_s))
        (self.src_archive if c.source == "archive" else self.src_manual).setChecked(True)
        self.archive.setText(c.archive_root); self.temp.setText(f"{c.temperature_K:g}")
        self.temp_mode.setCurrentIndex(0 if c.temperature_mode == "enz" else 1)
        self.press.setText(f"{c.pressure_hPa:g}"); self.h2o.setText(f"{c.h2o_percent:g}"); self.flow.setText(f"{c.flow_m3h:g}")
        self.jsrc.setCurrentIndex(0 if c.j_source == "measured" else 1)
        _set_rows(self.roof, [[o, cl] for o, cl in c.roof_open])
        self.constrain.setText(", ".join(c.constrain_species))

    def get(self) -> Conditions:
        c = Conditions()
        c.date, c.start, c.end, c.step_s = self.date.text().strip(), self.start.text().strip(), self.end.text().strip(), int(self.step.value())
        c.source = "archive" if self.src_archive.isChecked() else "manual"
        c.archive_root = self.archive.text().strip()
        c.temperature_K = float(self.temp.text() or 298); c.temperature_mode = "enz" if self.temp_mode.currentIndex() == 0 else "const"
        c.pressure_hPa = float(self.press.text() or 1006.5); c.h2o_percent = float(self.h2o.text() or 0); c.flow_m3h = float(self.flow.text() or 0)
        c.j_source = "measured" if self.jsrc.currentIndex() == 0 else "clear_sky"
        c.roof_open = [(o, cl) for o, cl in _rows(self.roof)]
        c.constrain_species = [x.strip() for x in self.constrain.text().split(",") if x.strip()]
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
        self.mode_compounds = QRadioButton("Compounds (MCM subset is generated)")
        self.mode_existing = QRadioButton("Existing mechanism include set")
        self.mode_compounds.setChecked(True)
        hb = QHBoxLayout(); hb.addWidget(self.mode_compounds); hb.addWidget(self.mode_existing); hb.addStretch(); lay.addLayout(hb)
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
        self._filter("")

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

    def _remove_selected(self):
        for it in self.chosen.selectedItems():
            self.chosen.takeItem(self.chosen.row(it))

    def set_from(self, spec: ExperimentSpec):
        self.chosen.clear()
        for p in spec.precursors:
            self.chosen.addItem(p)
        self.mech_combo.setCurrentText(spec.mechanism_stem)
        (self.mode_existing if spec.mechanism_stem else self.mode_compounds).setChecked(True)

    def apply(self, spec: ExperimentSpec):
        if self.mode_existing.isChecked():
            spec.mechanism_stem = self.mech_combo.currentText().strip(); spec.precursors = []
        else:
            spec.mechanism_stem = ""; spec.precursors = [self.chosen.item(i).text() for i in range(self.chosen.count())]

    def validatePage(self) -> bool:
        if self.mode_existing.isChecked() and not self.mech_combo.currentText().strip():
            QMessageBox.warning(self, "Chemistry", "Choose a mechanism include set."); return False
        if self.mode_compounds.isChecked() and self.chosen.count() == 0:
            QMessageBox.warning(self, "Chemistry", "Add at least one compound."); return False
        return True


# ─────────────────────────────────────────────────────────────────────────────
INJ_HEADERS = ["Name", "Product(s) — empty = recorded only", "Total [ppb]", "Start HH:MM", "Duration [min]",
               "ENZ column", "Switch", "Comment"]


class InjectionsPage(QWizardPage):
    def __init__(self, wiz):
        super().__init__()
        self.wiz = wiz
        self.setTitle("Injections")
        self.setSubTitle("What was injected, how much in total, when and over how long. One row per product channel; "
                         "rows with the same name belong to the same injection (e.g. impurities of a humidification).")
        lay = QVBoxLayout(self)
        self.table = _table(INJ_HEADERS); lay.addWidget(self.table)
        hb = QHBoxLayout()
        for label, fn in [("+ compounds from Chemistry", self._add_compounds), ("+ humidification", self._add_h2o),
                          ("+ O3", lambda: self._add_simple("O3", "O3", "50")), ("+ NO", lambda: self._add_simple("NO", "NO", "5")),
                          ("+ NO2", lambda: self._add_simple("NO2", "NO2", "5")), ("+ CO2 tracer", lambda: self._add_simple("CO2", "CO2", "0")),
                          ("+ seed (recorded only)", lambda: _add_row(self.table, ["seed", "", "", "", "10", "Q_seed", "", "recorded only"])),
                          ("+ row", lambda: _add_row(self.table)), ("− row", lambda: _del_row(self.table))]:
            b = QPushButton(label); b.clicked.connect(fn); hb.addWidget(b)
        hb.addStretch(); lay.addLayout(hb)
        lay.addWidget(_info("Amount = total injected in ppb (the backend spreads it over the duration). "
                            "Product(s) may be an MCM species or mole fractions of several (0.9*A+0.1*B). "
                            "Leave the column and switch empty for the standard names (QNO, QO3, QH2O, Q_<name>)."))
        lay.addWidget(QLabel("Follow-up reactions (advanced; e.g. TwoReaction chain: VOCin1 → VOCin_all, rate 1./60)"))
        self.chain = _table(["Reactants", "Products", "Rate expression"]); self.chain.setMaximumHeight(110); lay.addWidget(self.chain)
        hb2 = QHBoxLayout(); a = QPushButton("+"); a.clicked.connect(lambda: _add_row(self.chain)); d = QPushButton("−"); d.clicked.connect(lambda: _del_row(self.chain)); hb2.addWidget(a); hb2.addWidget(d); hb2.addStretch(); lay.addLayout(hb2)

    def _add_compounds(self):
        page = self.wiz.page_chem
        names = [page.chosen.item(i).text() for i in range(page.chosen.count())]
        if not names and page.mech_combo.currentText().strip():
            QMessageBox.information(self, "Compounds", "Add rows by hand: the compounds are defined by the existing mechanism."); return
        for n in names:
            _add_row(self.table, [n, n, "10", "08:00", "2", "", "", ""])

    def _add_h2o(self):
        _add_row(self.table, ["H2O", H2O_IMPURITY_DEFAULTS[0][0], H2O_IMPURITY_DEFAULTS[0][1], "05:30", "45", "QH2O", "qH2O",
                              "humidification: HONO impurity (SAPHIR default)"])
        for prod, ppb in H2O_IMPURITY_DEFAULTS[1:]:
            _add_row(self.table, ["H2O", prod, ppb, "", "", "QH2O", "qH2O", "impurity carried with the water"])

    def _add_simple(self, name, prod, ppb):
        _add_row(self.table, [name, prod, ppb, "09:30", "20" if name == "O3" else "5", "", "", ""])

    def set_from(self, spec: ExperimentSpec):
        rows = []
        for inj in spec.injections:
            rows.append([inj.name, inj.products, inj.total_ppb, inj.start, f"{inj.duration_min:g}", inj.column, inj.switch, inj.comment])
            for prod, ppb in inj.extra_channels:
                rows.append([inj.name, prod, ppb, "", "", inj.column, inj.switch, ""])
        _set_rows(self.table, rows)
        _set_rows(self.chain, [[a, b, c] for inj in spec.injections for a, b, c in inj.chain])

    def get(self) -> List[InjectionEvent]:
        events: Dict[str, InjectionEvent] = {}
        order: List[str] = []
        for name, prod, ppb, start, dur, col, sw, comment in _rows(self.table):
            if not name:
                continue
            ev = events.get(name)
            if ev is None:
                ev = InjectionEvent(name, prod, ppb or "0", start, float(dur or 1), col, sw, comment=comment)
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
     "The chamber is kept at slight overpressure; all species are diluted with the flow (DILUTE from FLOW, 276 m³)."),
    ("dark_hono", "Dark HONO formation from water on the walls",
     "Slow heterogeneous HONO production without light."),
    ("light_source", "NO2 source under illumination (SNO2)",
     "Light-driven NO2 release from the walls, used in some low-NOx experiments."),
]


class ChamberPage(QWizardPage):
    def __init__(self, wiz):
        super().__init__()
        self.wiz = wiz
        self.setTitle("What the chamber adds on its own")
        self.setSubTitle("SAPHIR chamber effects. All are on by default with the values used in the IDL scenarios; "
                         "open 'parameters' to change them.")
        lay = QVBoxLayout(self)
        self.checks: Dict[str, QCheckBox] = {}
        for key, title, desc in CHAMBER_ITEMS:
            cb = QCheckBox(title); cb.setChecked(key != "light_source"); self.checks[key] = cb
            lay.addWidget(cb); lay.addWidget(_info("      " + desc))
        self.params = QGroupBox("Parameters (advanced)"); self.params.setCheckable(True); self.params.setChecked(False)
        g = QGridLayout(self.params)
        self.hono_const = QLineEdit("2.6*0.7"); self.hono_lhs = QLineEdit("SHONO+HV"); self.hono_rhs = QLineEdit("HONO+XNOy")
        g.addWidget(QLabel("HONO source: SHONO = CONST("), 0, 0); g.addWidget(self.hono_const, 0, 1); g.addWidget(QLabel(") reaction"), 0, 2); g.addWidget(self.hono_lhs, 0, 3); g.addWidget(QLabel("-->"), 0, 4); g.addWidget(self.hono_rhs, 0, 5)
        self.hcho_const = QLineEdit("1"); self.hcho_lhs = QLineEdit("qSHCHO+HV"); self.hcho_rhs = QLineEdit("HCHO"); self.hcho_cname = QLineEdit("SHCHO")
        g.addWidget(QLabel("HCHO source:"), 1, 0); g.addWidget(self.hcho_cname, 1, 1); g.addWidget(self.hcho_const, 1, 2); g.addWidget(self.hcho_lhs, 1, 3); g.addWidget(QLabel("-->"), 1, 4); g.addWidget(self.hcho_rhs, 1, 5)
        self.y_ppb = QLineEdit("250"); self.y_name = QLineEdit("qY"); self.y_prod = QLineEdit("CH3O2")
        g.addWidget(QLabel("Background reactivity Y [ppb CO-eq]"), 2, 0); g.addWidget(self.y_ppb, 2, 1); g.addWidget(QLabel("name / product"), 2, 2); g.addWidget(self.y_name, 2, 3); g.addWidget(self.y_prod, 2, 5)
        self.dilute = QLineEdit("CONST((FLOW+0.010*60)/(276.*3600.))"); self.tracer = QLineEdit("XNOy, NA")
        g.addWidget(QLabel("DILUTE ="), 3, 0); g.addWidget(self.dilute, 3, 1, 1, 3); g.addWidget(QLabel("tracer dilution"), 3, 4); g.addWidget(self.tracer, 3, 5)
        self.dark_rate = QLineEdit("1e-012/10"); self.dark_rhs = QLineEdit("HONO+XNOy")
        g.addWidget(QLabel("Dark HONO: H2O --> "), 4, 0); g.addWidget(self.dark_rhs, 4, 1); g.addWidget(QLabel("k = CONST("), 4, 2); g.addWidget(self.dark_rate, 4, 3)
        self.light_rate = QLineEdit("jno2/5e-3/3600*HV_in/1.6*0.55"); self.light_rhs = QLineEdit("NO2+XNOy"); self.light_init = QLineEdit("M*1e-9*3.5*0")
        g.addWidget(QLabel("SNO2 --> "), 5, 0); g.addWidget(self.light_rhs, 5, 1); g.addWidget(QLabel("k = CONST("), 5, 2); g.addWidget(self.light_rate, 5, 3); g.addWidget(QLabel("INITIAL[SNO2]"), 5, 4); g.addWidget(self.light_init, 5, 5)
        lay.addWidget(self.params)

    def set_from(self, ch: ChamberEffects):
        for key, cb in self.checks.items():
            cb.setChecked(bool(getattr(ch, key)))
        self.hono_const.setText(ch.hono_const); self.hono_lhs.setText(ch.hono_lhs); self.hono_rhs.setText(ch.hono_rhs)
        self.hcho_const.setText(ch.hcho_const); self.hcho_lhs.setText(ch.hcho_lhs); self.hcho_rhs.setText(ch.hcho_rhs); self.hcho_cname.setText(ch.hcho_const_name)
        self.y_ppb.setText(ch.y_expr or ch.y_ppb); self.y_name.setText(ch.y_name); self.y_prod.setText(ch.y_product)
        self.dilute.setText(ch.dilute_expr); self.tracer.setText(", ".join(ch.tracer_dilution))
        self.dark_rate.setText(ch.dark_hono_rate); self.dark_rhs.setText(ch.dark_hono_rhs)
        self.light_rate.setText(ch.light_source_rate); self.light_rhs.setText(ch.light_source_rhs); self.light_init.setText(ch.light_source_initial)

    def get(self, base: ChamberEffects) -> ChamberEffects:
        ch = ChamberEffects(**{k: getattr(base, k) for k in ChamberEffects.__dataclass_fields__})
        for key, cb in self.checks.items():
            setattr(ch, key, cb.isChecked())
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
        return ch


# ─────────────────────────────────────────────────────────────────────────────
class InitialsPage(QWizardPage):
    def __init__(self, wiz):
        super().__init__()
        self.wiz = wiz
        self.setTitle("Initial conditions")
        self.setSubTitle("Concentrations at the start of the time window. Everything not listed starts at the floor value.")
        lay = QVBoxLayout(self)
        hb = QHBoxLayout(); hb.addWidget(QLabel("Floor for all species, INITIAL[*] = CONST(")); self.all_edit = QLineEdit("1."); hb.addWidget(self.all_edit); hb.addWidget(QLabel(") molecules/cm³")); hb.addStretch(); lay.addLayout(hb)
        self.table = _table(["Species", "Expression inside CONST( ), e.g. M*1e-9*14.0 for 14 ppb"]); lay.addWidget(self.table)
        hb2 = QHBoxLayout(); a = QPushButton("+ row"); a.clicked.connect(lambda: _add_row(self.table, ["", "M*1e-9*"])); d = QPushButton("− row"); d.clicked.connect(lambda: _del_row(self.table)); hb2.addWidget(a); hb2.addWidget(d); hb2.addStretch(); lay.addLayout(hb2)
        lay.addWidget(_info("Typical: CO 14 ppb (M*1e-9*14.0); with archive data the measured NO2/O3 at the start can be used."))

    def set_from(self, spec: ExperimentSpec):
        self.all_edit.setText(spec.initial_all)
        _set_rows(self.table, [[sp, ex] for sp, ex in spec.initials])

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
        self.summary = QPlainTextEdit(); self.summary.setReadOnly(True); self.summary.setMaximumHeight(170); lay.addWidget(self.summary)
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
                it = QListWidgetItem(("✖ " if level == "error" else "⚠ ") + m)
                self.checks.addItem(it)
                self._has_error |= level == "error"
            if not msgs:
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
class ExperimentWizard(QWizard):
    experiment_written = pyqtSignal(str)

    def __init__(self, parent=None, include_dir: str = "", easy_root: str = "", last_easy_dir: str = ""):
        super().__init__(parent)
        self.setWindowTitle("3. Build scenario — JUICE experiment wizard")
        self.resize(1250, 850)
        self.setWizardStyle(QWizard.WizardStyle.ModernStyle)
        self.include_dir, self.easy_root = include_dir, easy_root or last_easy_dir or os.getcwd()
        self._base = ExperimentSpec()
        self.page_start = StartPage(self); self.page_cond = ConditionsPage(self); self.page_chem = ChemistryPage(self)
        self.page_inj = InjectionsPage(self); self.page_chamber = ChamberPage(self); self.page_init = InitialsPage(self)
        self.page_review = ReviewPage(self)
        for p in (self.page_start, self.page_cond, self.page_chem, self.page_inj, self.page_chamber, self.page_init, self.page_review):
            self.addPage(p)
        self.load_spec(self._base, keep_name=True)

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
        return spec

    def out_dir(self) -> str:
        return os.path.join(self.page_start.out_edit.text().strip() or self.easy_root,
                            self.page_start.name_edit.text().strip() or "experiment")

    # ── finish ───────────────────────────────────────────────────────────────
    def accept(self):
        try:
            spec = self.collect_spec()
            out = self.out_dir()
            res = write_experiment(spec, out)
            save_spec(spec, os.path.join(out, spec.name + ".experiment.json"))
        except Exception as e:
            QMessageBox.critical(self, "Write failed", str(e)); return
        # Notes (e.g. archive columns replaced by constants) are kept for the caller
        # and printed; no modal dialog here so that headless use cannot block.
        self.last_result = res
        for n in res.get("notes", []):
            print(f"[wizard] {n}")
        super().accept()
        self.experiment_written.emit(res["easy"])


if __name__ == "__main__":
    import sys
    from PyQt6.QtWidgets import QApplication
    app = QApplication(sys.argv)
    w = ExperimentWizard(include_dir=sys.argv[1] if len(sys.argv) > 1 else "", easy_root=sys.argv[2] if len(sys.argv) > 2 else "")
    w.show()
    sys.exit(app.exec())
