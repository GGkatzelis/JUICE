"""A field study: one case × several input tables × emission factors, through the study machinery."""
import json
import os

import numpy as np
import pytest

from juice.core import field_spec as fs
from juice.core import field_scenario as fsc
from juice.core import study as S
from juice.core.scenario_prep import load_paths

INCLUDE_DIR = load_paths("")["include_dir"]


def _table(path, e_no):
    lines = ["time,lat,lon,T,p,RH,h,jNO2,E_NO [mol/s],E_C5H8 [mol/s],B_O3 [ppb]"]
    for i in range(3):
        hh = 20 + i
        lines.append(f"2021-08-15 {hh:02d}:00,34.1,-118.2,298,1000,40,800,0.008,{e_no},0.05,40")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return str(path)


@pytest.fixture
def field_study(tmp_path):
    t1, t2 = _table(tmp_path / "traj_01.csv", 2.0), _table(tmp_path / "traj_02.csv", 4.0)
    spec = fs.FieldSpec(name="parcel", run_type="lagrangian", table=t1, area_km2=64.0, default_emission_unit="mol/s", dilution="pbl",
                        chemistry={"base_mechanism": "MCM_v3.3.1", "precursors": ["C5H8"], "inorganic_only": False, "addons": [],
                                   "mechanism_stem": "", "include_dirs": [], "reaction_adjustments": []},
                        initials={"O3": 30.0, "NO2": 2.0})
    case = fs.save_spec(spec, str(tmp_path / "parcel"))
    st = S.new_study("nox_sweep", str(tmp_path / "nox_sweep"), "")
    st.field_case, st.tables = case, [t1, t2]
    st.axes = [S.Axis("emission_scale", "NOX", [0.5, 1.0])]
    st.export_kv_pd = False
    return st, (t1, t2)


def test_plans_and_validation(field_study):
    st, (t1, t2) = field_study
    plans = S.plan_runs(st)
    assert [p.name for p in plans] == ["traj_01__NOX_x-0.5", "traj_01__NOX_x-1", "traj_02__NOX_x-0.5", "traj_02__NOX_x-1"]
    assert plans[0].base_path == t1 and plans[2].base_path == t2
    assert S.bases(st) == [t1, t2]
    msgs = S.validate_study(st)
    assert not [m for lvl, m in msgs if lvl == "error"], msgs
    assert S.nox_axis(st.axes) == 0
    assert "emissions of NOX" in S.axis_description(st, st.axes[0])[0]


def test_validation_catches_the_wrong_things(field_study, tmp_path):
    st, _t = field_study
    st.axes.append(S.Axis("injection_scale", "all", [1, 2]))
    errors = [m for lvl, m in S.validate_study(st) if lvl == "error"]
    assert any("not for a field case" in e for e in errors)
    st.axes = [S.Axis("emission_scale", "NOX", [0.5])]
    st.tables.append(str(tmp_path / "missing.csv"))
    errors = [m for lvl, m in S.validate_study(st) if lvl == "error"]
    assert any("not found" in e for e in errors)
    st.tables = []; st.field_case = str(tmp_path / "nope.field.json")
    errors = [m for lvl, m in S.validate_study(st) if lvl == "error"]
    assert any("cannot read" in e for e in errors)


def test_materialise_writes_each_run_with_its_factor_and_shared_chemistry(field_study):
    st, _t = field_study
    plans = S.plan_runs(st)
    js = S.materialise(st, plans[0], INCLUDE_DIR)
    assert os.path.isfile(js) and js == plans[0].json_path
    easy = open(plans[0].json_path[:-5] + ".easy", encoding="latin-1").read()
    assert "emission scale factors of this run: NOX × 0.5" in easy
    run_spec = fs.load_spec(os.path.join(plans[0].folder, plans[0].name + fs.SPEC_SUFFIX))
    assert run_spec.emission_scale == {"NOX": 0.5} and run_spec.name == plans[0].name
    assert os.path.isfile(os.path.join(st.folder, "chemistry", "chemistry.json"))
    assert not os.path.isdir(os.path.join(plans[0].folder, "chemistry"))       # the subset is shared, not extracted per run
    # the factor reaches the ENZ: NO emission halved, isoprene untouched
    from juice.core.measured_inputs import read_enz
    _times, cols = read_enz(plans[0].json_path[:-5] + ".enz")
    S.materialise(st, plans[1], INCLUDE_DIR)
    _times, cols1 = read_enz(plans[1].json_path[:-5] + ".enz")
    assert cols["QE_NO"][0] == pytest.approx(0.5 * cols1["QE_NO"][0]) and cols["QE_C5H8"][0] == pytest.approx(cols1["QE_C5H8"][0])
    info = json.load(open(plans[0].info_path, encoding="utf-8"))
    assert info["table"].endswith("traj_01.csv") and "emissions NOX × 0.5" in info["notes"]
    assert plans[0].status == "ready"
    rows = S.collect_results(st, plans)
    assert len(rows) == 4 and rows[0]["experiment"] == "traj_01" and rows[0]["NOX_x"] == 0.5 and rows[0]["status"] == "ready"


def test_emission_factor_resolution():
    spec = fs.FieldSpec(emission_scale={"NOX": 0.5, "VOC": 2.0, "CO": 3.0})
    assert fs.emission_factor(spec, "NO") == 0.5 and fs.emission_factor(spec, "NO2") == 0.5
    assert fs.emission_factor(spec, "C5H8") == 2.0 and fs.emission_factor(spec, "CO") == 3.0
    assert fs.emission_factor(fs.FieldSpec(emission_scale={"ALL": 0.1}), "SO2") == 0.1
    assert fs.emission_factor(fs.FieldSpec(), "NO") == 1.0
    assert fs.emission_class("HONO") == "NOX" and fs.emission_class("SO2") == "OTHER" and fs.emission_class("TOLUENE") == "VOC"


def test_study_window_opens_a_field_study(field_study):
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PyQt6.QtWidgets import QApplication
    from juice.gui import study_window as sw
    app = QApplication.instance() or QApplication([])
    st, (t1, t2) = field_study
    os.makedirs(st.folder, exist_ok=True)
    path = S.save_study(st, os.path.join(st.folder, "nox_sweep.study.json"))
    w = sw.StudyWindow("", path)
    assert "Field study" in w.windowTitle()
    assert [w.exps.item(i).text() for i in range(w.exps.count())] == [t1, t2]
    assert "2 table(s) × 2 grid point(s) = 4 runs" in w.count.text()
    assert "Field study: 2 input table(s)" in w.ready.text()
    back = w._to_study()
    assert back.field_case == st.field_case and back.tables == [t1, t2] and [a.kind for a in back.axes] == ["emission_scale"]
    assert w.runs.rowCount() == 4
    w.close()
