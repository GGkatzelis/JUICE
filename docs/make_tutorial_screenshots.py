"""
make_tutorial_screenshots.py - the screenshots of docs/tutorial.md, rendered offscreen (no desktop needed).

    python docs/make_tutorial_screenshots.py [--work <folder>] [--no-run]

Builds the tutorial's experiments exactly as the text describes them (isoprene SAPHIR day from the
template, the UV-lamp template, the site_day field example), runs the isoprene day once for the results
viewer, and writes docs/img/tutorial/*.png. Rerun it after a GUI change so the tutorial stays current.
The experiments are written to --work (default: a folder in the system temp directory), never into the
repository.
"""
import os, sys, glob, argparse, subprocess, tempfile, shutil
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("QT_QPA_FONTDIR", "C:/Windows/Fonts")
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
OUT = os.path.join(HERE, "img", "tutorial")
sys.path.insert(0, ROOT); os.chdir(ROOT)

ap = argparse.ArgumentParser()
ap.add_argument("--work", default=os.path.join(tempfile.gettempdir(), "juice_tutorial_shots"))
ap.add_argument("--no-run", action="store_true", help="reuse the isoprene run of an earlier call")
args = ap.parse_args()
os.makedirs(OUT, exist_ok=True); os.makedirs(args.work, exist_ok=True)

from PyQt6.QtWidgets import QApplication
from PyQt6.QtCore import Qt
app = QApplication(sys.argv)
try:
    from juice_theme import apply_theme; apply_theme(app)
except Exception:
    pass
W, H = 1400, 900
TOML = os.path.join(ROOT, "FR_GUI.toml")
made = []


def settle(n=10):
    for _ in range(n):
        app.processEvents()


NEUTRAL_FOLDER = r"C:\JUICE\experiments"


def mask_paths(w):
    """No personal paths in public screenshots. Every line edit showing a path under the user's home or the
    temp work folder shows it under a neutral project folder instead, keeping the part below the work folder
    (a table still reads ...\\site_day\\site_day.csv). Returns the originals; grab() restores them right after
    the capture, so the window's state never changes. Signals are blocked while masking."""
    from PyQt6.QtWidgets import QLineEdit, QLabel
    home = os.path.normcase(os.path.expanduser("~"))
    work = os.path.normcase(os.path.abspath(args.work))
    saved = []
    for e in w.findChildren(QLineEdit) + w.findChildren(QLabel):   # labels too: the Run step shows the folder as text
        t = e.text()
        if not t or "\n" in t or "<" in t or len(t) > 400:
            continue
        n = os.path.normcase(os.path.abspath(t)) if (":" in t[:3] or t.startswith(("/", "\\"))) else ""
        if not n or not (n.startswith(home) or n.startswith(work) or os.sep + "users" + os.sep in n):
            continue
        tail = os.path.relpath(n, work) if n.startswith(work) else "."
        e.blockSignals(True); e.setText(NEUTRAL_FOLDER if tail == "." else os.path.join(NEUTRAL_FOLDER, tail))
        e.blockSignals(False); saved.append((e, t))
    return saved


def set_field(w, placeholder_start, text):
    """Fill the line edit whose placeholder starts with the given text (e.g. the experiment name)."""
    from PyQt6.QtWidgets import QLineEdit
    for e in w.findChildren(QLineEdit):
        if e.placeholderText().startswith(placeholder_start):
            e.blockSignals(True); e.setText(text); e.blockSignals(False)


def show_template(w, title):
    """Show the template the tutorial picks in the template chooser, without re-applying it."""
    from PyQt6.QtWidgets import QComboBox
    for cb in w.findChildren(QComboBox):
        i = cb.findText(title, Qt.MatchFlag.MatchStartsWith)
        if i >= 0:
            cb.blockSignals(True); cb.setCurrentIndex(i); cb.blockSignals(False)


def grab(w, name, size=(W, H)):
    """The offscreen platform clamps shown windows to 800x600; a window not shown on screen keeps its size."""
    w.setAttribute(Qt.WidgetAttribute.WA_DontShowOnScreen, True)
    w.show(); w.resize(*size); settle(); saved = mask_paths(w); w.resize(*size); settle()
    path = os.path.join(OUT, name)
    ok = w.grab().save(path); made.append(name)
    for e, t in saved:
        e.blockSignals(True); e.setText(t); e.blockSignals(False); print(f"{name:42s} {'ok' if ok else 'FAILED'} {w.width()}x{w.height()}")


import experiment_spec as es, scenario_prep as sp


def tutorial_isoprene():
    """Section 1 of the tutorial: template 'VOC photo-oxidation at low NOx', typed conditions, isoprene."""
    _, _, s = es.load_template(es.list_templates()["voc_photooxidation_low_nox"])
    s.name = "tutorial_isoprene"
    c = s.conditions
    c.source, c.date, c.start, c.end, c.step_s = "manual", "2024-07-09", "07:00", "16:00", 60
    c.temperature_K, c.temperature_mode, c.pressure_hPa, c.h2o_percent, c.flow_m3h = 298.0, "const", 1013.0, 1.0, 0.5
    c.roof_open = [("09:00", "16:00")]
    c.j_source = "clear_sky"
    s.precursors = ["C5H8"]
    s.injections = [es.InjectionEvent(name="isoprene", start="08:30", products="C5H8", total_ppb="10"),
                    es.InjectionEvent(name="O3", start="09:00", products="O3", total_ppb="40")]
    return s


# ── 0. the two start windows ──
import juice_start
grab(juice_start.FirstChoice(TOML), "00_first_choice.png", (900, 520))
grab(juice_start.StartWindow(TOML), "01_start_menu.png", (1000, 640))

# ── 1. the experiment window, every step, on the tutorial experiment ──
import juice_gui
spec = tutorial_isoprene()
w = juice_gui.JuiceWindow(TOML); w.load_spec(spec, keep_name=True); settle()
set_field(w, "e.g. isoprene", "tutorial_isoprene"); show_template(w, "VOC photo-oxidation at low NOx")
names = ["10_step1_where_and_when.png", "11_step2_chemistry.png", "12_step3_injections.png",
         "13_step4_chamber.png", "14_step5_run.png"]
for row, name in enumerate(names):
    w.steps.setCurrentRow(row); settle(); grab(w, name)
w.close()

# ── run the experiment once for the results viewer ──
exp_dir = os.path.join(args.work, "tutorial_isoprene")
cs = os.path.join(exp_dir, "tutorial_isoprene.cs")
if not (args.no_run and os.path.isfile(cs)):
    shutil.rmtree(exp_dir, ignore_errors=True)
    res = es.write_experiment(spec, exp_dir)
    es.save_spec(spec, os.path.join(exp_dir, "tutorial_isoprene.experiment.json"))   # the viewer's overview reads it
    js = sp.prepare_scenario(res, sp.load_paths(TOML)["include_dir"], export_kv_pd=True, toml_path=TOML,
                             solver=sp.solver_from_toml(TOML, "fast"))
    print("running", js)
    r = subprocess.run([sys.executable, os.path.join(ROOT, "scenario_runner.py"), js], cwd=ROOT,
                       capture_output=True, text=True, env={**os.environ, "PYTHONIOENCODING": "utf-8"})
    if r.returncode:
        print(r.stdout[-2000:], r.stderr[-2000:]); sys.exit("the tutorial run failed")
    cs = next(iter(glob.glob(os.path.join(exp_dir, "*.cs"))), cs)

# ── 2. the results viewer, one picture per tab ──
import results_viewer as rv
v = rv.ResultsViewer(rv.load_column_file(cs), csv_path=cs, toml_path=TOML)
from PyQt6.QtWidgets import QCheckBox


def no_observations():
    """Synthetic day: never overlay the real archive data of that date (Overview and Timeseries have a box each)."""
    for box in v.findChildren(QCheckBox):
        if box.text() in ("SAPHIR observations", "observations") and box.isChecked():
            box.setChecked(False)
    settle(20)


no_observations()
grab(v, "20_viewer_overview.png")
tabs = v.right_tabs
for i in range(tabs.count()):
    label = tabs.tabText(i)
    if label in ("Timeseries", "Process Analysis", "Reactivity", "Budgets", "Chemical Space"):
        tabs.setCurrentIndex(i); settle(20); no_observations()
        if label == "Process Analysis":
            v._update_budget("OH"); settle(20)        # the tutorial text: click OH
        grab(v, f"2{i + 1}_viewer_{label.replace(' ', '_').lower()}.png")
v.close()

# ── 3. the UV-lamp template, first step ──
_, _, lamp = es.load_template(es.list_templates()["uv_chamber_lamps"])
lamp.name = "tutorial_uv_lamps"
lamp.precursors = ["C3H6"]                       # the tutorial: "an alkene and NO, say"
lamp.injections = [es.InjectionEvent(name="NO", start="08:15", products="NO", total_ppb="20"),
                   es.InjectionEvent(name="propene", start="08:30", products="C3H6", total_ppb="50")]
w = juice_gui.JuiceWindow(TOML); w.load_spec(lamp, keep_name=True); settle()
set_field(w, "e.g. isoprene", "tutorial_uv_lamps"); show_template(w, "UV-lamp chamber — photooxidation under blacklights")
w.steps.setCurrentRow(0); grab(w, "30_uv_lamps_where_and_when.png")
w.steps.setCurrentRow(3); grab(w, "31_uv_lamps_chamber.png")
w.close()

# ── 4. the field window on a copy of the site_day example, every step ──
import field_window as fw, field_spec as fs
case = fs.copy_case(os.path.join(ROOT, "field_examples", "site_day", "site_day.field.json"), args.work)
f = fw.FieldWindow(case, TOML); settle(); show_template(f, "site_day")
for row, name in enumerate(["40_field_step1_case.png", "41_field_step2_chemistry.png",
                            "42_field_step3_emissions_constraints.png", "43_field_step4_mixing.png",
                            "44_field_step5_run.png"]):
    f.steps.setCurrentRow(row); settle(); grab(f, name)
f.close()

# ── 4b. the site_day run in the viewer: HCHO modelled against its O_ column, and the radicals ──
fr_dir = os.path.join(args.work, "site_day_run")
fcs = os.path.join(fr_dir, "site_day.cs")
if not (args.no_run and os.path.isfile(fcs)):
    shutil.rmtree(fr_dir, ignore_errors=True)
    r = subprocess.run([sys.executable, os.path.join(ROOT, "field_scenario.py"), case, "--out", fr_dir, "--run"],
                       cwd=ROOT, capture_output=True, text=True, env={**os.environ, "PYTHONIOENCODING": "utf-8"})
    print("field run rc", r.returncode, "".join(l + "\n" for l in r.stdout.splitlines() if "HCHO" in l or "skill" in l))
    if r.returncode:
        print(r.stdout[-2000:], r.stderr[-2000:]); sys.exit("the field run failed")
v = rv.ResultsViewer(rv.load_column_file(fcs), csv_path=fcs, toml_path=TOML)
tabs = v.right_tabs
for i in range(tabs.count()):
    if tabs.tabText(i) == "Timeseries":
        tabs.setCurrentIndex(i); settle(20)
ts = v.ts_tab
ts.clear_all()
ts.add_species(["HCHO"], unit="ppb", panel=1)
ts.add_species(["OH", "HO2"], unit="native", panel=2, log=True)
ts.layout_combo.setCurrentText("2"); settle(30)
grab(v, "45_field_result_hcho.png")
v.close()

# ── 4c. the two other kinds of run: a plume and a Lagrangian air parcel ──
for key, shots in (("plume_ageing", [(0, "46_plume_case.png"), (3, "47_plume_mixing.png")]),
                   ("la_parcel", [(0, "48_parcel_case.png"), (2, "49_parcel_emissions.png")])):
    c2 = fs.copy_case(os.path.join(ROOT, "field_examples", key, f"{key}.field.json"), args.work)
    f = fw.FieldWindow(c2, TOML); settle(); show_template(f, key)
    for row, name in shots:
        f.steps.setCurrentRow(row); settle(); grab(f, name)
    f.close()

# ── 5. the study of Sect. 4: initial NO x isoprene amount on the tutorial experiment, run once ──
import study as S, study_window as sw
st_dir = os.path.join(args.work, "tutorial_isopleth")
st_path = os.path.join(st_dir, "tutorial_isopleth.study.json")
if not (args.no_run and os.path.isfile(st_path)):
    shutil.rmtree(st_dir, ignore_errors=True); os.makedirs(st_dir)
    st = S.new_study("tutorial_isopleth", st_dir, TOML)
    st.mode, st.focus, st.workers = "sweep", "ozone", 4
    st.experiments = [os.path.join(exp_dir, "tutorial_isoprene.experiment.json")]
    st.axes = [S.Axis("initial", "NO", S.parse_values("log 0.1 20 6")),
               S.Axis("injection_scale", "isoprene", S.parse_values("log 0.5 4 6"))]
    S.save_study(st, st_path)
    s0 = sw.StudyWindow(TOML, study_path=st_path); settle(20)
    s0._save(st_path); s0.close()        # saved as the window plans it, so its run names are the window's
    print("running the study (36 runs)")
    r = subprocess.run([sys.executable, os.path.join(ROOT, "study_runner.py"), st_path, "--workers", "4", "--yes"],
                       cwd=ROOT, capture_output=True, text=True, env={**os.environ, "PYTHONIOENCODING": "utf-8"})
    print("study runner rc", r.returncode, r.stdout[-400:])
s = sw.StudyWindow(TOML, study_path=st_path); settle(20)
grab(s, "50_study_setup.png", (1900, 1000))
for i in range(s.tabs.count()):
    label = s.tabs.tabText(i)
    if label in ("Results", "Plot"):
        s.tabs.setCurrentIndex(i); settle(30)
        for m in ("_refresh_results", "_plot"):
            try:
                getattr(s, m)(); settle(30)
            except Exception as e:
                print(m, "failed:", e)
        grab(s, f"5{1 if label == 'Results' else 2}_study_{label.lower()}.png", (1900, 1000))
s.close()

print(f"\n{len(made)} screenshots in {OUT}")
