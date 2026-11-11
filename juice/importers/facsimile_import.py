"""
facsimile_import.py — an experiment set up for FACSIMILE (.easy + ENZ) → ExperimentSpec
========================================================================================

    python -m juice.importers.facsimile_import path\\to\\<name>.easy      # opens the experiment window on it

An experiment modelled with FACSIMILE earlier (a campaign day, an archived scenario) becomes an
ordinary JUICE experiment: the master file is read back into the experiment description
(`import_facsimile`), the ENZ of that run becomes the measured-inputs file, and the experiment
opens in the step-by-step window, where every part can be looked at, changed and run.

What is read back:
    .easy   includes → mechanism include set; WALLDEP / SOURCE / BGREAC / DIL blocks → chamber
            effects; MEASURED block → J-value source, temperature mode, pressure; FILES / FACS;
            INITIAL; the injection blocks (switch = input(COL), channels, delayed release); every
            `k[...] = CONST(...)` on a mechanism reaction → reaction adjustment; everything else
            verbatim as expert overrides (the review shows the regenerated EASY).
    ENZ     date, time window and step; roof times from hv; injection start and duration from the
            pulse columns (value 1/n over n rows); the measured columns (t, H2O, dilute, sz, J values,
            NO/NO2/O3) stay in the file, which the experiment reads as its measured inputs.

The import checks itself: the EASY regenerated from the description is compared statement by
statement with the original, and the ENZ columns regenerated from the description with the
original ones.  Differences come back as import notes (shown under 'Needs attention').

The include files of the mechanism stay where they are: the experiment folder and its
chemistry\\ subfolder are added to the experiment's include folders (Chemistry step).  The EASY
version that made the FACSIMILE run (read from the .fac / .out next to it) is kept with the
experiment, because 5.96b and 6.00 apply a redefined FUNCTION differently.
"""
from __future__ import annotations
from juice._root import ROOT as _JUICE_ROOT

import os
import re
import sys
import tempfile
from datetime import datetime, timedelta
from typing import Dict, List, Optional, Tuple

import numpy as np

HERE = _JUICE_ROOT
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from juice.core.experiment_spec import (ExperimentSpec, Conditions, ChamberEffects, InjectionEvent, ENZ_EPOCH,     # noqa: E402
                             STANDARD_COLUMNS, enz_columns, render_easy)
from juice.core.easy_builder import load_form, EasyForm                                                             # noqa: E402
from juice.core.easy_scene import statements_equivalent                                                             # noqa: E402
from juice.core import measured_inputs as mi                                                                             # noqa: E402

_RX_ADJUST = re.compile(r"^k\[(.+?)\]\s*=\s*CONST\s*\((.*)\)\s*(?:;\s*(.*))?$", re.I)
_RX_CONSTRAIN_IN = re.compile(r"^(\w+)_in\s*=\s*input\s*\(\s*(\w+)\s*\)\s*$", re.I)
_RX_CONSTRAIN = re.compile(r"^(\w+)\s*=\s*CONST\s*\(\s*(\w+)_in\s*\*\s*1e-9\s*\*\s*M\s*\)\s*$", re.I)


def enz_name_in(path: str) -> str:
    """The ENZ file an EASY / elog names in FILES[ENZ] = …  ('' when none)."""
    try:
        with open(path, encoding="latin-1") as fh:
            for line in fh:
                s = line.strip()
                if s.startswith((";", "*")):
                    continue
                m = re.match(r"FILES\s*\[\s*ENZ\s*\]\s*=\s*(.+?)\s*(?:;.*)?$", s.split(";")[0].strip(), re.I)
                if m:
                    return m.group(1).strip().strip('"\'')
    except OSError:
        pass
    return ""


def find_enz(folder: str, stem: str, named: str) -> str:
    """The ENZ to use: the one FILES[ENZ] names, else <stem>.enz, else the only .enz in the folder."""
    for cand in ([named] if named else []) + [stem + ".enz"]:
        p = os.path.join(folder, cand)
        if os.path.isfile(p):
            return p
    enz = [f for f in os.listdir(folder) if f.lower().endswith(".enz")] if os.path.isdir(folder) else []
    return os.path.join(folder, enz[0]) if len(enz) == 1 else ""


def safe_name(stem: str) -> str:
    """An experiment name the GUI accepts (letters, digits, '_', '-', '.')."""
    return re.sub(r"[^A-Za-z0-9_.\-]+", "_", stem).strip("_") or "experiment"


def name_from_column(col: str) -> str:
    """The injection name behind an ENZ switch column: QNO → NO, Q_ethen → ethen, qVCP2 → VCP2,
    QCl2 → Cl2 (the standard columns map back exactly; others drop the leading q / q_)."""
    for name, std in STANDARD_COLUMNS.items():
        if std.upper() == col.upper():
            return {"CL2": "Cl2", "SEED": "seed"}.get(name, name)
    if col[:2].lower() == "q_":
        return col[2:] or col
    if col[:1].lower() == "q" and len(col) > 1:
        return col[1:]
    return col


def _hhmm(seconds_since_2000: float) -> str:
    return (ENZ_EPOCH + timedelta(seconds=float(seconds_since_2000))).strftime("%H:%M")


def pulse_of(values: np.ndarray, times: np.ndarray, step_s: float) -> Tuple[str, float, bool]:
    """(start HH:MM, duration in minutes, exact) of an injection switch column.  Exact when the
    non-zero rows form one block holding 1/n each (the IDL convention: integral 1); otherwise the
    first and last non-zero rows bound the pulse and the caller notes the approximation.
    ('', 0, True) when the column is 0 throughout."""
    nz = np.flatnonzero(np.asarray(values, float) != 0)
    if not len(nz):
        return "", 0.0, True
    i0, i1 = int(nz[0]), int(nz[-1])
    n = i1 - i0 + 1
    block = values[i0:i1 + 1]
    exact = len(nz) == n and np.allclose(block, 1.0 / n, rtol=1e-3, atol=1e-12)   # the file holds 1/n rounded to a few digits
    return _hhmm(times[i0]), n * step_s / 60.0, exact


def _blocks(flag: np.ndarray) -> List[Tuple[int, int]]:
    """[(first index, one past last)] of the runs where flag is true."""
    out, start = [], None
    for i, on in enumerate(list(flag) + [False]):
        if on and start is None:
            start = i
        elif not on and start is not None:
            out.append((start, i)); start = None
    return out


def _const_inner(expr: str) -> str:
    m = re.match(r"CONST\s*\((.*)\)\s*$", (expr or "").strip(), re.I | re.S)
    return m.group(1).strip() if m else (expr or "").strip()


def _chamber_from_form(f: EasyForm, notes: List[str]) -> ChamberEffects:
    ch = ChamberEffects(
        walldep=f.walldep, walldep_value=f.walldep_value or "1.",
        hono=f.hono_source, hono_lhs=f.hono_lhs, hono_rhs=f.hono_rhs, hono_const=f.hono_const, hono_rate=f.hono_rate,
        hcho=f.hcho_source, hcho_lhs=f.hcho_lhs, hcho_rhs=f.hcho_rhs, hcho_const_name=f.hcho_const_name,
        hcho_const=f.hcho_const, hcho_rate=f.hcho_rate,
        bgreac=f.bgreac, y_name=f.y_name, y_product=f.y_product, y_ppb=f.y_ppb, y_expr=f.y_const_expr, y_rate=f.y_rate,
        dilution=f.dilution, dilute_expr=f.dilute_expr, dil_value=f.dil_value or "CONST(1.)", dil_flow_const=f.dil_flow_const,
        tracer_dilution=list(f.tracer_dilution),
        dark_hono=f.dark_hono, dark_hono_rhs=f.dark_hono_rhs, dark_hono_rate=f.dark_hono_rate,
        light_source=f.light_source, light_source_rhs=f.light_source_rhs,
        light_source_rate=f.light_source_rate, light_source_initial=f.light_source_initial)
    if f.dilution and not f.dil_value:
        notes.append("DILUTE without a DIL switch in the EASY: DIL = CONST(1.) assumed")
    return ch


def import_facsimile(easy_path: str, enz_path: str, name: str = "", easy_version: Optional[str] = None,
                     reference_cs: str = "", check: bool = True) -> Tuple[ExperimentSpec, List[str]]:
    """The experiment behind a FACSIMILE master file and its ENZ, as (spec, notes).  `notes` say
    what could not be read back exactly; with `check` the regenerated EASY and ENZ are compared with
    the originals and every difference is noted.  The notes are also stored in spec.import_notes."""
    easy_path = os.path.abspath(easy_path); enz_path = os.path.abspath(enz_path)
    folder, stem = os.path.dirname(easy_path), os.path.splitext(os.path.basename(easy_path))[0]
    notes: List[str] = []
    f = load_form(easy_path)
    times, raw_cols = mi.read_enz(enz_path)
    if len(times) < 2:
        raise ValueError(f"{os.path.basename(enz_path)}: fewer than two rows")
    col_of = {k.upper(): k for k in raw_cols}                  # FACSIMILE reads INPUT(jno2) from a column 'jNO2'

    def enz_col(name: str) -> Optional[np.ndarray]:
        return raw_cols.get(col_of.get(name.upper(), ""))
    step_s = float(np.median(np.diff(times)))
    if not np.allclose(np.diff(times), step_s, atol=1e-6):
        notes.append(f"the ENZ rows are not evenly spaced: step {step_s:g} s assumed (the experiment's grid is regular)")

    # ── where and when ──
    c = Conditions()
    t0 = ENZ_EPOCH + timedelta(seconds=float(times[0])); t1 = ENZ_EPOCH + timedelta(seconds=float(times[-1]))
    c.date, c.start, c.end, c.step_s = t0.strftime("%Y-%m-%d"), t0.strftime("%H:%M"), t1.strftime("%H:%M"), int(round(step_s))
    if t1.date() != t0.date():
        notes.append(f"the ENZ runs past midnight (to {t1.strftime('%Y-%m-%d %H:%M')}): the time window ends at {c.end} of {c.date}")
    c.source, c.input_file = "file", enz_path
    settings = {r.name.upper(): r for r in f.setting_rows}
    t_row, p_row = settings.get("T"), settings.get("P")
    t_val = _const_inner(t_row.value) if (t_row is not None and t_row.mode == "const") else ""
    c.temperature_mode = "const" if t_val else "enz"
    if t_val and re.fullmatch(r"[0-9.eE+\-]+", t_val):
        c.temperature_K = float(t_val)
    elif t_val:
        notes.append(f"T = CONST({t_val}) is not a number: the GUI holds T at {c.temperature_K:g} K (change it on 'Where and when')")
    elif enz_col("t") is not None:
        c.temperature_K = float(np.round(np.nanmean(enz_col("t")), 1))      # the fallback constant only; T follows the ENZ
    p_val = _const_inner(p_row.value) if (p_row is not None and p_row.mode == "const") else ""
    if p_val and re.fullmatch(r"[0-9.eE+\-]+", p_val):
        c.pressure_hPa = float(p_val)
    elif p_row is not None:
        notes.append(f"pressure P = {p_row.value or ('INPUT(' + p_row.column + ')')}: the GUI writes P as a number; {c.pressure_hPa:g} hPa used")
    if enz_col("H2O") is not None:
        c.h2o_percent = float(np.round(np.nanmean(enz_col("H2O")), 3))
    if enz_col("dilute") is not None:
        c.flow_m3h = float(np.round(np.nanmean(enz_col("dilute")), 2))
    hv = enz_col("hv")
    if hv is not None:
        c.roof_open = [(_hhmm(times[a]), _hhmm(times[b]) if b < len(times) else _hhmm(times[-1] + step_s))
                       for a, b in _blocks(np.asarray(hv) > 0.5)]
    c.j_source = "measured" if any(r.mode in ("enz", "input") for r in f.j_rows) else "clear_sky"
    for r in f.setting_rows + f.j_rows:
        if r.mode in ("enz", "input") and r.column and enz_col(r.column) is None:
            notes.append(f"{r.name} = INPUT({r.column}): no column '{r.column}' in the ENZ")
    for other in [r for r in f.setting_rows if r.name.upper() not in ("T", "P", "SZA", "SZ", "HV", "FLOW", "H2O")]:
        notes.append(f"setting {other.name} ({other.mode}) is not one the GUI writes: dropped")

    # ── chemistry and files ──
    spec = ExperimentSpec(name=safe_name(name or stem), title=f.title, conditions=c)
    spec.mechanism_stem, spec.mech_incl_suffix = f.mech_stem, f.mech_incl_suffix
    spec.use_anorg, spec.use_j_clear_sky_incl = f.use_anorg, f.use_j_clear_sky
    spec.extra_includes = list(f.includes_pre) + list(f.includes_other)
    if f.includes_pre:
        notes.append("include(s) before the inorganic ones (" + ", ".join(f.includes_pre) + ") are written after the mechanism includes")
    if not f.mech_stem:
        notes.append("no mechanism include set (<stem>.fac.mech …) in the EASY: choose the chemistry on the Chemistry step")
    spec.hmax, spec.outstep = f.hmax or "30", f.outstep or "100000"
    spec.initial_all = next((_const_inner(ex) for sp, ex, _c in f.initials if sp == "*"), "1.")
    spec.initials = [(sp, _const_inner(ex)) for sp, ex, _c in f.initials if sp != "*"]
    spec.chamber = _chamber_from_form(f, notes)
    if f.jmeasured_include != "saphir.jmeasured_02.incl":
        notes.append("the EASY " + ("has no saphir.jmeasured_02.incl" if not f.jmeasured_include else f"includes {f.jmeasured_include}")
                     + "; the GUI writes saphir.jmeasured_02.incl")
    if f.dil_stem:
        notes.append(f"dilution files of another stem ({f.dil_stem}.dil.*): the GUI writes those of the mechanism ({f.mech_stem}.dil.*)")

    # ── injections: one event per switch, channels from the EASY, timing from the ENZ ──
    by_switch: Dict[str, Tuple[str, list]] = {}
    for r in f.injections:
        by_switch.setdefault(r.switch, (r.column, []))[1].append(r)
    events: List[InjectionEvent] = []
    used_names: set = set()
    for sw, (col, rows) in by_switch.items():
        chans = [r for r in rows if r.products.strip() or r.rate.strip()]
        nm = name_from_column(col)
        k = 2
        while nm in used_names:
            nm = f"{name_from_column(col)}_{k}"; k += 1
        used_names.add(nm)
        ev = InjectionEvent(nm, chans[0].products if chans else "", (chans[0].rate or "0") if chans else "0",
                            column=col, switch=sw, comment=(chans[0].comment if chans else (rows[0].comment if rows else "")),
                            extra_channels=[(r.products, r.rate or "0") for r in chans[1:]])
        ev.action = ev.action_key()
        if enz_col(col) is None:
            notes.append(f"injection {nm}: switch column '{col}' is not in the ENZ — no pulse, start left empty")
        else:
            start, dur, exact = pulse_of(enz_col(col), times, step_s)
            if not start:
                notes.append(f"injection {nm}: column '{col}' is 0 throughout the ENZ (never injected)")
            else:
                ev.start, ev.duration_min = start, dur
                if not exact:
                    notes.append(f"injection {nm}: column '{col}' is not one pulse of 1/n — start {start} and {dur:g} min "
                                 "bound its non-zero rows; the regenerated pulse differs from the file")
        events.append(ev)
    taken = {e.column.upper() for e in events}
    for col, vals in raw_cols.items():                      # ENZ columns with a pulse but no EASY block: recorded only
        if col[:1].upper() != "Q" or col.upper() in taken:
            continue
        start, dur, exact = pulse_of(vals, times, step_s)
        if not start:
            continue
        nm = name_from_column(col); k = 2
        while nm in used_names:
            nm = f"{name_from_column(col)}_{k}"; k += 1
        used_names.add(nm)
        ev = InjectionEvent(nm, "", "0", start, dur, column=col, comment="recorded in the ENZ, not modelled in the EASY")
        ev.action = ev.action_key(); events.append(ev)
        if not exact:
            notes.append(f"recorded event {nm}: column '{col}' is not one pulse of 1/n — timing approximated")
    if f.chain_reactions:
        target = next((e for e in reversed(events) if e.products), None)
        if target is not None:
            target.chain = [tuple(x) for x in f.chain_reactions]
        else:
            notes.append("delayed-release reactions without a modelled injection: kept as expert overrides")
            for a, b, r in f.chain_reactions:
                spec.overrides += [f"{a} --> {b}", f"k[{a} --> {b}] = CONST({r})"]
    spec.injections = events

    # ── everything else: rate adjustments to the Chemistry step, constrained species, the rest verbatim ──
    rest: List[str] = []
    pending_in: Dict[str, str] = {}
    for code in f.extra_statements:
        code = code.strip()
        if not code:
            continue
        m = _RX_CONSTRAIN_IN.match(code)
        if m and m.group(1).upper() == m.group(2).upper():
            pending_in[m.group(1).upper()] = code; continue
        m = _RX_CONSTRAIN.match(code)
        if m and m.group(1).upper() == m.group(2).upper() and m.group(1).upper() in pending_in:
            pending_in.pop(m.group(1).upper()); c.constrain_species.append(m.group(1)); continue
        m = _RX_ADJUST.match(code)
        if m:
            spec.reaction_adjustments.append((m.group(1).strip(), m.group(2).strip(), (m.group(3) or "").strip())); continue
        rest.append(code)
    rest = list(pending_in.values()) + rest
    spec.overrides += rest
    for k_, v in f.extra_facs:
        spec.overrides.append(f"FACS[{k_}] = {v}")

    # ── provenance ──
    spec.include_dirs = [d for d in (folder, os.path.join(folder, "chemistry")) if os.path.isdir(d)]
    if easy_version is None:
        from juice.core.elog_generator import easy_version as _ver
        easy_version = _ver(folder, stem=stem)
    spec.easy_version = easy_version or ""
    spec.reference_cs = reference_cs or (os.path.join(folder, stem + ".cs") if os.path.isfile(os.path.join(folder, stem + ".cs")) else "")
    spec.import_notes = []

    if check:
        notes += fidelity(spec, easy_path, enz_path)
    spec.import_notes = list(notes)
    return spec, notes


def fidelity(spec: ExperimentSpec, easy_path: str, enz_path: str) -> List[str]:
    """How the regenerated inputs differ from the originals: EASY statement by statement, the ENZ
    column by column (grid, switch columns, hv and the measured ones).  Empty when nothing differs."""
    out: List[str] = []
    tmp = tempfile.NamedTemporaryFile("w", suffix=".easy", delete=False, encoding="latin-1", errors="replace", newline="\n")
    try:
        tmp.write(render_easy(spec)); tmp.close()
        ok, diffs = statements_equivalent(easy_path, tmp.name, ignore_duplicate_reactions=True)
    finally:
        try:
            os.unlink(tmp.name)
        except OSError:
            pass
    for d in diffs:
        if "FILES[ENZ]" in d.upper():                  # JUICE names the ENZ after the experiment: always different, never a loss
            continue
        out.append("EASY: " + d)
    try:
        times, raw = mi.read_enz(enz_path)
        grid, cols, _n = enz_columns(spec)
        if len(grid) != len(times) or not np.allclose(grid, times, atol=1e-3):
            out.append(f"ENZ: time grid differs ({len(times)} rows in the file, {len(grid)} regenerated)")
        else:
            for name, vals in raw.items():
                if name not in cols:
                    if np.any(vals != 0) and name[:1].upper() == "Q":
                        out.append(f"ENZ: column {name} is not regenerated")
                    continue
                if not np.allclose(vals, cols[name], rtol=1e-3, atol=1e-9):      # pulse values are 1/n rounded in the file
                    out.append(f"ENZ: column {name} differs (max |Δ| {np.max(np.abs(vals - cols[name])):.3g})")
    except Exception as e:
        out.append(f"ENZ: could not be checked ({e})")
    return out


# ─────────────────────────────────────────────────────────────────────────────
# Qt dialog
# ─────────────────────────────────────────────────────────────────────────────

try:
    from PyQt6.QtCore import QSettings
    from PyQt6.QtWidgets import (QDialog, QDialogButtonBox, QFileDialog, QGridLayout, QHBoxLayout, QLabel, QLineEdit,
                                 QListWidget, QMessageBox, QPlainTextEdit, QPushButton, QVBoxLayout, QWidget)
    _QT = True
except ImportError:
    _QT = False


if _QT:
    def _info(text: str) -> QLabel:
        lab = QLabel(text); lab.setWordWrap(True); lab.setStyleSheet("color: #555555; font-size: 11px;")
        return lab

    class ImportDialog(QDialog):
        """Choose the .easy and ENZ of a FACSIMILE experiment, see what the import makes of them, and
        hand the ExperimentSpec to the experiment window.  The ENZ is copied into the experiment
        folder (<project>/<name>/<name>.measured.enz); nothing in the original folder is changed."""

        def __init__(self, parent=None, path: str = "", project_root: str = "", include_dir: str = ""):
            super().__init__(parent)
            self.setWindowTitle("Import a FACSIMILE experiment")
            self.resize(900, 640)
            self.settings = QSettings("FZJ", "JUICE")
            if not include_dir or not project_root:
                from juice.core.scenario_prep import load_paths
                paths = load_paths("")
                include_dir, project_root = include_dir or paths.get("include_dir", ""), project_root or paths.get("output_root", HERE)
            self.project_root, self.include_dir = project_root, include_dir
            self.result_spec: Optional[ExperimentSpec] = None
            self.result_notes: List[str] = []
            self.folder = self.stem = ""
            v = QVBoxLayout(self)
            v.addWidget(_info("An experiment modelled with FACSIMILE earlier: its master file (.easy) and the ENZ of that run become an "
                              "ordinary JUICE experiment — conditions, chemistry, injections and chamber effects are filled in and can "
                              "be changed before the run. The ENZ is copied into the experiment folder; the original folder is not changed."))
            g = QGridLayout(); r = 0
            g.addWidget(QLabel("EASY file"), r, 0); self.file = QLineEdit(); self.file.setPlaceholderText(".easy of the FACSIMILE run (an .elog is accepted when its .easy lies next to it)")
            self.file.editingFinished.connect(lambda: self._file_changed(self.file.text())); g.addWidget(self.file, r, 1)
            b = QPushButton("Browse …"); b.clicked.connect(self._pick_file); g.addWidget(b, r, 2); r += 1
            self.found = _info(""); g.addWidget(self.found, r, 1, 1, 2); r += 1
            g.addWidget(QLabel("ENZ file"), r, 0); self.enz = QLineEdit(); g.addWidget(self.enz, r, 1)
            b = QPushButton("Browse …"); b.clicked.connect(self._pick_enz); g.addWidget(b, r, 2); r += 1
            g.addWidget(QLabel("Include folders"), r, 0)
            self.dirs = QListWidget(); self.dirs.setMaximumHeight(70)
            self.dirs.setToolTip("searched (with their subfolders) for the #include files of the EASY, before JUICE's own include folder; "
                                 "the experiment folder and its chemistry\\ subfolder are added by themselves")
            for d in self.settings.value("facsimile/include_dirs", []) or []:
                if os.path.isdir(d):
                    self.dirs.addItem(d)
            g.addWidget(self.dirs, r, 1)
            bb = QVBoxLayout(); b = QPushButton("Add …"); b.clicked.connect(self._add_dir); bb.addWidget(b)
            b = QPushButton("Remove"); b.clicked.connect(self._remove_dir); bb.addWidget(b); bb.addStretch(); g.addLayout(bb, r, 2); r += 1
            g.addWidget(QLabel("Experiment name"), r, 0); self.name = QLineEdit(); g.addWidget(self.name, r, 1); r += 1
            g.addWidget(QLabel("Project folder"), r, 0); self.project = QLineEdit(project_root); g.addWidget(self.project, r, 1)
            b = QPushButton("Browse …"); b.clicked.connect(self._pick_project); g.addWidget(b, r, 2); r += 1
            v.addLayout(g)
            hb = QHBoxLayout()
            self.read_btn = QPushButton("Read the files"); self.read_btn.setObjectName("primary"); self.read_btn.clicked.connect(self.read)
            hb.addWidget(self.read_btn); hb.addStretch(); v.addLayout(hb)
            v.addWidget(QLabel("What the import found — ⚠ marks what could not be read back exactly (kept as notes with the experiment)"))
            self.report = QPlainTextEdit(); self.report.setReadOnly(True); v.addWidget(self.report, stretch=1)
            self.buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
            self.buttons.button(QDialogButtonBox.StandardButton.Ok).setText("Open in the experiment window")
            self.buttons.button(QDialogButtonBox.StandardButton.Ok).setEnabled(False)
            self.buttons.accepted.connect(self._accept); self.buttons.rejected.connect(self.reject); v.addWidget(self.buttons)
            if path:
                self.file.setText(path); self._file_changed(path)

        # ── choosing the files ──
        def _start_dir(self) -> str:
            last = self.settings.value("facsimile/last_dir", "")
            return last if last and os.path.isdir(last) else (self.project_root if os.path.isdir(self.project_root) else HERE)

        def _pick_file(self):
            p, _ = QFileDialog.getOpenFileName(self, "FACSIMILE experiment", self._start_dir(), "FACSIMILE experiment (*.easy *.elog);;All files (*)")
            if p:
                self.file.setText(p); self._file_changed(p)

        def _file_changed(self, p: str):
            p = p.strip().strip('"')
            if not p or not os.path.isfile(p) or p == getattr(self, "_loaded", ""):
                return                                   # editingFinished fires on every focus change: only a new path counts
            self._loaded = p
            from juice.core.elog_generator import easy_version, function_scope_for, FUNCTION_SCOPE_USE
            self.settings.setValue("facsimile/last_dir", os.path.dirname(p))
            self.folder = os.path.dirname(p); self.stem = os.path.splitext(os.path.basename(p))[0]
            self.easy_path = os.path.join(self.folder, self.stem + ".easy")
            self.ver = easy_version(self.folder, stem=self.stem)
            self.enz.setText(find_enz(self.folder, self.stem, enz_name_in(p)))
            self.name.setText(safe_name(self.stem))
            parts = []
            if not os.path.isfile(self.easy_path):
                parts.append(f"✖ no {self.stem}.easy next to it — the import reads the master file (an .elog holds every include inline)")
            parts.append(f"EASY {self.ver} ({'5.96b rule: a redefined FUNCTION applies from its line on' if function_scope_for(self.ver) == FUNCTION_SCOPE_USE else '6.00 rule: the last FUNCTION definition applies everywhere'})"
                         if self.ver else "no .fac/.out next to it: EASY version unknown, JUICE's own rules (6.00)")
            cs = os.path.join(self.folder, self.stem + ".cs")
            parts.append("FACSIMILE results (.cs) found: 'Compare with FACSIMILE' on the Run step" if os.path.isfile(cs) else "no FACSIMILE results (.cs) next to it")
            self.found.setText(" · ".join(parts))
            self.result_spec = None; self.buttons.button(QDialogButtonBox.StandardButton.Ok).setEnabled(False); self.report.clear()

        def _pick_enz(self):
            p, _ = QFileDialog.getOpenFileName(self, "ENZ file", self.folder or self._start_dir(), "ENZ (*.enz);;All files (*)")
            if p:
                self.enz.setText(p)

        def _pick_project(self):
            p = QFileDialog.getExistingDirectory(self, "Project folder", self.project.text() or self.project_root)
            if p:
                self.project.setText(p)

        def _add_dir(self):
            p = QFileDialog.getExistingDirectory(self, "Folder with include files", self.folder or self._start_dir())
            if p and not any(self.dirs.item(i).text() == p for i in range(self.dirs.count())):
                self.dirs.addItem(p); self._save_dirs()

        def _remove_dir(self):
            for it in self.dirs.selectedItems():
                self.dirs.takeItem(self.dirs.row(it))
            self._save_dirs()

        def _save_dirs(self):
            self.settings.setValue("facsimile/include_dirs", [self.dirs.item(i).text() for i in range(self.dirs.count())])

        # ── import ──
        def read(self):
            easy = getattr(self, "easy_path", "")
            if not easy or not os.path.isfile(easy):
                QMessageBox.warning(self, "Import", "Choose the .easy of a FACSIMILE experiment first."); return
            enz = self.enz.text().strip()
            if not os.path.isfile(enz):
                QMessageBox.warning(self, "Import", "The ENZ file was not found: choose it."); return
            try:
                spec, notes = import_facsimile(easy, enz, self.name.text().strip(), easy_version=self.ver)
            except Exception as e:
                self.report.setPlainText(f"✖ The files could not be read: {e}"); return
            spec.include_dirs += [self.dirs.item(i).text() for i in range(self.dirs.count()) if self.dirs.item(i).text() not in spec.include_dirs]
            self.result_spec, self.result_notes = spec, notes
            self._missing = self._missing_includes(spec)
            self.report.setPlainText(self._summary(spec, notes))
            self.buttons.button(QDialogButtonBox.StandardButton.Ok).setEnabled(True)

        def _missing_includes(self, spec: ExperimentSpec) -> List[str]:
            """Include files named in the EASY that none of the folders holds (the run would stop on them).
            Looked up as elog_generator does: by file name, ignoring case, in every subfolder."""
            names = []
            try:
                with open(self.easy_path, encoding="latin-1") as fh:
                    for line in fh:
                        m = re.match(r"^\s*#include\s+(.+?)\s*(?:;.*)?$", line)
                        if m:
                            names.append(m.group(1).strip())
            except OSError:
                return []
            have = set()
            for d in list(spec.include_dirs) + [self.include_dir]:
                if d and os.path.isdir(d):
                    for _root, _dirs, files in os.walk(d):
                        have.update(f.lower() for f in files)
            return [n for n in names if os.path.basename(n).lower() not in have]

        def _summary(self, spec: ExperimentSpec, notes: List[str]) -> str:
            c = spec.conditions
            lines = [f"{spec.name}: {c.date} {c.start}–{c.end} UTC, step {c.step_s} s; measured conditions from the ENZ "
                     f"({len([1 for r in mi.read_enz(c.input_file)[1]])} columns)",
                     f"Roof open: {', '.join(f'{o}–{cl}' for o, cl in c.roof_open) or 'never (dark)'};  photolysis: {c.j_source}; "
                     f"T {'from the ENZ' if c.temperature_mode == 'enz' else f'{c.temperature_K:g} K'}, P {c.pressure_hPa:g} hPa",
                     f"Chemistry: {spec.mechanism_stem or '(none)'}" + (f" + {', '.join(spec.extra_includes)}" if spec.extra_includes else "")
                     + f";  {len(spec.reaction_adjustments)} rate adjustment(s);  include folders: {len(spec.include_dirs)}",
                     "Chamber: " + ", ".join(k for k in ("walldep", "hono", "hcho", "bgreac", "dilution", "dark_hono", "light_source") if getattr(spec.chamber, k))
                     + f";  HONO strength {spec.chamber.hono_const}, Y {spec.chamber.y_expr or spec.chamber.y_ppb + ' ppb'}",
                     f"Initial values: {len(spec.initials)};  expert overrides: {len(spec.overrides)} line(s)" + (f";  EASY {spec.easy_version}" if spec.easy_version else ""),
                     f"Injections ({len(spec.injections)}):"]
            for inj in spec.injections:
                what = f"{inj.total_ppb} ppb → {inj.products[:50]}{'…' if len(inj.products) > 50 else ''}" if inj.products else "recorded only"
                lines.append(f"   {inj.name:14s} {inj.start or '--:--'} +{inj.duration_min:g} min   {what}" + (f"  (+{len(inj.extra_channels)} channel(s))" if inj.extra_channels else ""))
            if self._missing:
                lines.append(f"✖ include file(s) not found in any folder: {', '.join(self._missing)} — add the folder that holds them under 'Include folders'")
            lines.append("")
            lines.append("✔ The regenerated EASY and ENZ match the originals statement by statement and column by column." if not notes
                         else f"{len(notes)} note(s):")
            lines += ["⚠ " + n for n in notes]
            return "\n".join(lines)

        def _accept(self):
            spec = self.result_spec
            if spec is None:
                return
            name = safe_name(self.name.text().strip() or spec.name); spec.name = name
            project = self.project.text().strip() or self.project_root
            out = os.path.join(project, name)
            try:
                import shutil
                os.makedirs(out, exist_ok=True)
                dst = os.path.join(out, name + ".measured.enz")
                if os.path.abspath(spec.conditions.input_file) != os.path.abspath(dst):
                    shutil.copy2(spec.conditions.input_file, dst)
                spec.conditions.input_file = dst
            except Exception as e:
                QMessageBox.critical(self, "Import", f"The ENZ could not be copied into {out}:\n{e}"); return
            self.project_root = project
            self.accept()


def main(argv=None) -> int:
    """`python -m juice.importers.facsimile_import <file.easy>`: the experiment window opens its import dialog on the file."""
    argv = list(sys.argv if argv is None else argv)
    from juice.gui import juice_gui
    return juice_gui.main(argv)


if __name__ == "__main__":
    sys.exit(main())
