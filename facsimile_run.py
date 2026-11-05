"""
facsimile_run.py — run an experiment that was set up for FACSIMILE (an .easy or its .elog) in JUICE.

    python facsimile_run.py [<file.easy | file.elog>] [-t toml]

For experiments modelled earlier with FACSIMILE (a campaign day, an archived scenario): JUICE runs
the same inputs and, when the FACSIMILE results (.cs) lie next to them, compares the two runs.
Nothing is written into the original folder — the inputs are copied to a work folder under
output_root/facsimile/<name>.

Two sources:
  .elog  the file EASY wrote when FACSIMILE ran: every include is inline, exactly as used then.
         Preferred when it exists next to the .easy.
  .easy  the master file; its #include files are looked up in the experiment folder, its
         chemistry/ subfolder, the include folders added here (remembered) and JUICE's own
         include folder.

The EASY version that made the run (read from the .fac / .out next to it) decides how a FUNCTION
macro defined more than once applies — EASY 5.96b and 6.00 differ (elog_generator.FUNCTION_SCOPE_*).
"""
from __future__ import annotations

import contextlib
import io
import json
import os
import re
import shutil
import sys
import time
from typing import Optional

from PyQt6.QtCore import Qt, QProcess, QProcessEnvironment, QSettings
from PyQt6.QtWidgets import (QApplication, QCheckBox, QComboBox, QFileDialog, QGridLayout, QGroupBox, QHBoxLayout, QLabel,
                             QLineEdit, QListWidget, QMessageBox, QProgressBar, QPushButton, QRadioButton, QTextEdit,
                             QVBoxLayout, QWidget)

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from scenario_prep import load_paths, pick_toml, solver_profiles, default_profile, solver_from_toml, prepare_scenario  # noqa: E402
from elog_generator import easy_version, function_scope_for, FUNCTION_SCOPE_USE  # noqa: E402

COMPARE_FLOOR = 2.5e4      # molec/cm³ ≈ 1e-3 ppt: as tests/test_solver.py — below this, CPU rounding dominates
COMPARE_TOL = 0.02


def _info(text: str) -> QLabel:
    lab = QLabel(text); lab.setWordWrap(True); lab.setStyleSheet("color: #555555; font-size: 11px;")
    return lab


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


class FacsimileRunWindow(QWidget):
    def __init__(self, toml: str = "", path: str = ""):
        super().__init__()
        import juice_theme as theme
        self.toml = toml
        self.paths = load_paths(toml)
        self.settings = QSettings("FZJ", "JUICE")
        self._proc: Optional[QProcess] = None
        self._t0 = 0.0
        self.json_path = ""
        self.stem = self.folder = self.easy_path = self.elog_path = ""
        self.setWindowTitle("JUICE — Run a FACSIMILE experiment"); self.setWindowIcon(theme.app_icon()); self.setMinimumWidth(900)
        v = QVBoxLayout(self); v.setContentsMargins(16, 12, 16, 12); v.setSpacing(8)
        t = QLabel("Run a FACSIMILE experiment"); t.setObjectName("title"); v.addWidget(t)
        v.addWidget(_info("An experiment modelled with FACSIMILE earlier (its .easy, or the .elog EASY wrote for that run): JUICE runs "
                          "the same inputs in a work folder and, when the FACSIMILE results lie next to them, compares the two runs. "
                          "The original folder is not changed."))

        g = QGroupBox("Experiment"); gl = QGridLayout(g); r = 0
        gl.addWidget(QLabel("File"), r, 0); self.file = QLineEdit(); self.file.setPlaceholderText(".easy or .elog of the FACSIMILE run")
        self.file.editingFinished.connect(lambda: self._file_changed(self.file.text())); gl.addWidget(self.file, r, 1)
        b = QPushButton("Browse …"); b.clicked.connect(self._pick_file); gl.addWidget(b, r, 2); r += 1
        self.found = _info(""); gl.addWidget(self.found, r, 1, 1, 2); r += 1
        gl.addWidget(QLabel("Source"), r, 0)
        hb = QHBoxLayout()
        self.src_elog = QRadioButton("the .elog of that run (every include inline, exactly as used — recommended)")
        self.src_easy = QRadioButton("the .easy, with the include files looked up in the folders below")
        self.src_elog.setChecked(True); hb.addWidget(self.src_elog); hb.addWidget(self.src_easy); hb.addStretch()
        gl.addLayout(hb, r, 1, 1, 2); r += 1
        gl.addWidget(QLabel("ENZ file"), r, 0); self.enz = QLineEdit(); gl.addWidget(self.enz, r, 1)
        b = QPushButton("Browse …"); b.clicked.connect(self._pick_enz); gl.addWidget(b, r, 2); r += 1
        gl.addWidget(QLabel("Include folders"), r, 0)
        self.dirs = QListWidget(); self.dirs.setMaximumHeight(70); self.dirs.setToolTip("searched (with their subfolders) for #include files of the .easy; remembered")
        for d in self.settings.value("facsimile/include_dirs", []) or []:
            if os.path.isdir(d):
                self.dirs.addItem(d)
        gl.addWidget(self.dirs, r, 1)
        bb = QVBoxLayout(); b = QPushButton("Add …"); b.clicked.connect(self._add_dir); bb.addWidget(b)
        b = QPushButton("Remove"); b.clicked.connect(self._remove_dir); bb.addWidget(b); bb.addStretch(); gl.addLayout(bb, r, 2); r += 1
        gl.addWidget(QLabel("Work folder"), r, 0); self.work = QLineEdit(); gl.addWidget(self.work, r, 1)
        b = QPushButton("Browse …"); b.clicked.connect(self._pick_work); gl.addWidget(b, r, 2); r += 1
        v.addWidget(g)

        g = QGroupBox("Run"); gl = QGridLayout(g); r = 0
        gl.addWidget(QLabel("Solver"), r, 0); self.solver = QComboBox()
        self._profiles = solver_profiles(toml)
        for name, p in self._profiles.items():
            self.solver.addItem(p["label"], name)
        dflt = default_profile(toml)
        if dflt in self._profiles:
            self.solver.setCurrentIndex(self.solver.findData(dflt))
        self.solver.setToolTip("\n\n".join(f"{p['label']}:\n{p['summary']}" for p in self._profiles.values()))
        self.solver.currentIndexChanged.connect(self._solver_changed); gl.addWidget(self.solver, r, 1, 1, 2); r += 1
        self.kvpd = QCheckBox("Also write the reaction rates and production/loss budgets (large files; for the process analysis "
                              "and reactivity views in the viewer)"); gl.addWidget(self.kvpd, r, 0, 1, 3); r += 1
        v.addWidget(g)

        hb = QHBoxLayout()
        self.prep_btn = QPushButton("1. Prepare"); self.prep_btn.setObjectName("primary"); self.prep_btn.clicked.connect(self.prepare); hb.addWidget(self.prep_btn)
        self.run_btn = QPushButton("2. Run"); self.run_btn.setObjectName("primary"); self.run_btn.setEnabled(False); self.run_btn.clicked.connect(self.run); hb.addWidget(self.run_btn)
        self.stop_btn = QPushButton("Stop"); self.stop_btn.setEnabled(False); self.stop_btn.clicked.connect(self.stop); hb.addWidget(self.stop_btn)
        self.view_btn = QPushButton("3. Open results"); self.view_btn.setObjectName("primary"); self.view_btn.setEnabled(False); self.view_btn.clicked.connect(self.open_results); hb.addWidget(self.view_btn)
        self.cmp_btn = QPushButton("Compare with FACSIMILE"); self.cmp_btn.setEnabled(False); self.cmp_btn.clicked.connect(self.compare)
        self.cmp_btn.setToolTip(f"max relative deviation over the run, ±{COMPARE_TOL*100:g} %; concentrations measured against at least "
                                f"{COMPARE_FLOOR:g} molec/cm³ (≈ 1e-3 ppt); settings (J values, rate constants) without a floor")
        hb.addWidget(self.cmp_btn); hb.addStretch(); v.addLayout(hb)
        self.bar = QProgressBar(); self.bar.setRange(0, 100); self.bar.setTextVisible(True); self.bar.setFormat("not started"); v.addWidget(self.bar)
        self.log = QTextEdit(); self.log.setReadOnly(True); self.log.setMinimumHeight(220)
        self.log.setStyleSheet("font-family: Consolas, monospace; font-size: 11px; background: #fbfbfb; color: #222222;")
        v.addWidget(self.log, stretch=1)
        if path:
            self.file.setText(path); self._file_changed(path)

    # ── choosing the experiment ──
    def _start_dir(self) -> str:
        last = self.settings.value("facsimile/last_dir", "")
        return last if last and os.path.isdir(last) else next((d for d in (self.paths.get("easy_root"), HERE) if d and os.path.isdir(d)), HERE)

    def _pick_file(self):
        p, _ = QFileDialog.getOpenFileName(self, "FACSIMILE experiment", self._start_dir(), "FACSIMILE experiment (*.easy *.elog);;All files (*)")
        if p:
            self.file.setText(p); self._file_changed(p)

    def _file_changed(self, p: str):
        p = p.strip().strip('"')
        if not p or not os.path.isfile(p):
            return
        self.settings.setValue("facsimile/last_dir", os.path.dirname(p))
        self.folder = os.path.dirname(p); self.stem = os.path.splitext(os.path.basename(p))[0]
        self.easy_path = os.path.join(self.folder, self.stem + ".easy"); self.elog_path = os.path.join(self.folder, self.stem + ".elog")
        has_easy, has_elog = os.path.isfile(self.easy_path), os.path.isfile(self.elog_path)
        self.src_elog.setEnabled(has_elog); self.src_easy.setEnabled(has_easy)
        (self.src_elog if has_elog else self.src_easy).setChecked(True)
        self.ver = easy_version(self.folder, stem=self.stem)
        self.enz.setText(find_enz(self.folder, self.stem, enz_name_in(p)))
        self.fac_cs = os.path.join(self.folder, self.stem + ".cs")
        parts = [f"EASY {self.ver} ({'5.96b rule: a redefined FUNCTION applies from its line on' if function_scope_for(self.ver) == FUNCTION_SCOPE_USE else '6.00 rule: the last FUNCTION definition applies everywhere'})"
                 if self.ver else "no .fac/.out next to it: EASY version unknown, JUICE's own rules (6.00)"]
        parts.append("FACSIMILE results (.cs) found: can be compared after the run" if os.path.isfile(self.fac_cs) else "no FACSIMILE results (.cs) next to it")
        if not has_elog:
            parts.append("no .elog: the include files are looked up")
        self.found.setText(" · ".join(parts))
        self.work.setText(os.path.join(self.paths.get("output_root", HERE), "facsimile", self.stem))
        self.json_path = ""; self.run_btn.setEnabled(False)
        self.view_btn.setEnabled(os.path.isfile(self.cs_path())); self.cmp_btn.setEnabled(os.path.isfile(self.cs_path()) and os.path.isfile(self.fac_cs))

    def _pick_enz(self):
        p, _ = QFileDialog.getOpenFileName(self, "ENZ file", self.folder or self._start_dir(), "ENZ (*.enz);;All files (*)")
        if p:
            self.enz.setText(p)

    def _pick_work(self):
        p = QFileDialog.getExistingDirectory(self, "Work folder", self.work.text() or self.paths.get("output_root", HERE))
        if p:
            self.work.setText(p)

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

    def _solver_changed(self, *_a):
        self.json_path = ""; self.run_btn.setEnabled(False)

    # ── prepare / run / results ──
    def cs_path(self) -> str:
        return os.path.join(self.work.text(), self.stem + ".cs") if self.stem and self.work.text() else ""

    def prepare(self) -> str:
        src = self.elog_path if self.src_elog.isChecked() else self.easy_path
        if not self.stem or not os.path.isfile(src):
            self.log.append("Choose the .easy or .elog of a FACSIMILE experiment first."); return ""
        enz = self.enz.text().strip()
        if not os.path.isfile(enz):
            self.log.append("The ENZ file was not found: choose it."); return ""
        work = self.work.text().strip()
        self.bar.setValue(0); self.bar.setFormat("preparing …"); QApplication.processEvents()
        try:
            os.makedirs(work, exist_ok=True)
            easy_copy = os.path.join(work, self.stem + ".easy")
            shutil.copy2(src, easy_copy)
            # the ENZ keeps the name the file gives in FILES[ENZ], so that line still resolves
            enz_copy = os.path.join(work, enz_name_in(src) or os.path.basename(enz))
            shutil.copy2(enz, enz_copy)
            dirs = [self.folder, os.path.join(self.folder, "chemistry")] + [self.dirs.item(i).text() for i in range(self.dirs.count())]
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                self.json_path = prepare_scenario({"easy": easy_copy, "enz": enz_copy}, self.paths["include_dir"], self.kvpd.isChecked(),
                                                  self.toml, extra_dirs=dirs, solver=solver_from_toml(self.toml, self.solver.currentData()),
                                                  easy_version=self.ver)
            report = buf.getvalue()
            with open(os.path.join(work, self.stem + "_prepare.log"), "w", encoding="utf-8") as fh:
                fh.write(report)
        except Exception as e:
            self.log.append(f"Prepare failed: {e}"); self.bar.setValue(0); self.bar.setFormat("prepare failed"); return ""
        missing = []
        try:
            with open(os.path.join(work, self.stem + "_juice.elog"), encoding="utf-8", errors="replace") as fh:
                missing = sorted({m.group(1).strip() for m in re.finditer(r"^;#include (.+?) ! NOT FOUND", fh.read(), re.M)})
        except OSError:
            pass
        self.log.append(f"Source: {os.path.basename(src)}  ·  ENZ: {os.path.basename(enz)}  ·  EASY {self.ver or '(unknown)'}  ·  "
                        f"solver {self.solver.currentText()}")
        for line in report.splitlines():
            if "[scenario_prep]" in line or "FUNCTION scope" in line:
                self.log.append(line)
        if missing:
            self.log.append(f"✖ {len(missing)} include file(s) not found: " + ", ".join(missing) +
                            "\n  Add the folder that holds them under 'Include folders', or use the .elog of the run.")
        self.log.append(f"Inputs ready in {work}  (details: {self.stem}_prepare.log)")
        self.bar.setValue(100); self.bar.setFormat("inputs ready" if not missing else f"inputs ready — {len(missing)} include(s) missing")
        self.run_btn.setEnabled(True); self.view_btn.setEnabled(os.path.isfile(self.cs_path()))
        return self.json_path

    def run(self):
        if not self.json_path:
            return
        args = ["-u", os.path.join(HERE, "scenario_runner.py"), self.json_path] + (["--toml", self.toml] if self.toml else [])
        self._proc = QProcess(self)
        env = QProcessEnvironment.systemEnvironment(); env.insert("PYTHONIOENCODING", "utf-8"); env.insert("PYTHONUTF8", "1")
        self._proc.setProcessEnvironment(env); self._proc.setProcessChannelMode(QProcess.ProcessChannelMode.MergedChannels)
        self._proc.readyReadStandardOutput.connect(self._on_output); self._proc.finished.connect(self._on_finished)
        self._t0 = time.monotonic(); self.bar.setValue(0); self.bar.setFormat("running   0%")
        self.run_btn.setEnabled(False); self.prep_btn.setEnabled(False); self.stop_btn.setEnabled(True); self.view_btn.setEnabled(False); self.cmp_btn.setEnabled(False)
        self.log.append("$ python scenario_runner.py " + os.path.basename(self.json_path))
        self._proc.start(sys.executable, args)

    def stop(self):
        if self._proc and self._proc.state() != QProcess.ProcessState.NotRunning:
            self._proc.kill()

    def _on_output(self):
        raw = bytes(self._proc.readAllStandardOutput()).decode("utf-8", errors="replace")
        self.log.moveCursor(self.log.textCursor().MoveOperation.End); self.log.insertPlainText(raw); self.log.moveCursor(self.log.textCursor().MoveOperation.End)
        for m in re.finditer(r"Fraction Done:\s*([\d.]+)", raw):
            pct = min(int(float(m.group(1)) * 100), 99); el = int(time.monotonic() - self._t0); mm, ss = divmod(el, 60)
            self.bar.setValue(pct); self.bar.setFormat(f"running   {mm}m {ss:02d}s   {pct}%")

    def _on_finished(self, code, _status):
        el = int(time.monotonic() - self._t0); mm, ss = divmod(el, 60)
        ok = code == 0 and os.path.isfile(self.cs_path())
        try:
            with open(self.json_path[:-5] + ".run.json", "w", encoding="utf-8") as fh:
                json.dump({"exit": code, "elapsed": time.monotonic() - self._t0, "finished": time.strftime("%Y-%m-%d %H:%M:%S"),
                           "cs_bytes": os.path.getsize(self.cs_path()) if ok else 0}, fh, indent=1)
        except Exception:
            pass
        self.bar.setValue(100 if ok else 0); self.bar.setFormat(f"done   {mm}m {ss:02d}s" if ok else f"failed (exit {code})   {mm}m {ss:02d}s")
        self.run_btn.setEnabled(True); self.prep_btn.setEnabled(True); self.stop_btn.setEnabled(False); self.view_btn.setEnabled(ok)
        self.cmp_btn.setEnabled(ok and os.path.isfile(self.fac_cs))
        if ok and os.path.isfile(self.fac_cs):
            self.log.append("FACSIMILE results found next to the inputs: 'Compare with FACSIMILE' sets the two runs against each other.")

    def open_results(self):
        cs = self.cs_path()
        if os.path.isfile(cs):
            p = QProcess(self); p.start(sys.executable, [os.path.join(HERE, "results_viewer.py"), cs]); self._viewer = p

    def compare(self):
        import compare_juice_facsimile as cjf
        work = self.work.text().strip(); log_path = os.path.join(work, f"{self.stem}_difference_JUICE_Facsimile.log")
        self.bar.setFormat("comparing …"); QApplication.processEvents()
        buf = io.StringIO()
        try:
            with contextlib.redirect_stdout(buf):
                s = cjf.compare(self.stem, COMPARE_TOL, no_stop=True, floor=COMPARE_FLOOR, juice_dir_override=work, fac_dir_override=self.folder)
        except Exception as e:
            self.log.append(f"Comparison failed: {e}"); self.bar.setFormat("comparison failed"); return
        with open(log_path, "w", encoding="utf-8") as fh:
            fh.write(buf.getvalue())
        self.log.append(f"JUICE vs FACSIMILE, ±{COMPARE_TOL*100:g} % (concentrations measured against ≥ {COMPARE_FLOOR:g} molec/cm³):")
        for label, g in (s.get("groups") or {}).items():
            worst = ", ".join(f"{n} {r*100:.2g} %" for r, n in g["worst"][:3])
            self.log.append(f"  {label}: {g['n'] - g['n_exceed']}/{g['n']} within  ·  worst {worst}")
        self.log.append(("✔ the two runs agree" if s.get("ok") else "✖ differences above the tolerance") + f"  (full report: {os.path.basename(log_path)})")
        self.bar.setFormat("agree with FACSIMILE" if s.get("ok") else "differences — see the log")


def main(argv=None) -> int:
    argv = list(sys.argv if argv is None else argv)
    import juice_theme as theme
    toml = pick_toml(argv)
    skip = {i + 1 for i, a in enumerate(argv) if a in ("-t", "--toml")}
    files = [a for i, a in enumerate(argv[1:], 1) if i not in skip and not a.startswith("-") and os.path.isfile(a)]
    theme.taskbar_identity()
    app = QApplication(argv); theme.apply_theme(app)
    w = FacsimileRunWindow(toml, files[0] if files else ""); w.show(); theme.force_taskbar_icon(w)
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
