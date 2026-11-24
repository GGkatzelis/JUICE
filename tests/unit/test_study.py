"""study.py, the parts that need no run: axis values and names, the grid of runs, the changes a study applies to an
experiment, the results read from a run (on hand-made series), run-time estimates, the ozone ridge."""
import json
import math
import os
import time
from datetime import datetime, timezone

import numpy as np
import pytest

from juice.core import study as sd
from juice.core.experiment_spec import ExperimentSpec, InjectionEvent
from juice.core.study import Axis, Modification, RunPlan, Study


# ── axis values, names, grid ─────────────────────────────────────────────────

def test_values_as_a_list_or_a_range():
    assert sd.parse_values("0.5, 1; 2  4") == [0.5, 1.0, 2.0, 4.0]
    assert sd.parse_values("lin 0 20 5") == [0.0, 5.0, 10.0, 15.0, 20.0]
    assert sd.parse_values("LOG 0.1 10 3") == [0.1, 1.0, 10.0]
    assert sd.parse_values(" log 1 100 5 ")[1] == 3.16228              # six significant digits


def test_log_range_needs_positive_limits_and_text_is_an_error():
    with pytest.raises(ValueError, match="positive"):
        sd.parse_values("log 0 10 3")
    with pytest.raises(ValueError):
        sd.parse_values("one, two")


@pytest.mark.parametrize("axis, name", [
    (Axis("initial", "NO"), "NO_t0_ppb"),
    (Axis("injection_set", "ethanol"), "ethanol_ppb"),
    (Axis("injection_scale", "voc"), "voc_x"),
    (Axis("rate_factor", "C5H8 + OH --> CISOPA"), "C5H8 + OH --> CISOPA_x"),
    (Axis("emission_scale", "NOX"), "NOX_x"),
    (Axis("injection_new", "no|11:00|5"), "NO_added_1100_ppb"),
    (Axis("chamber", "hono_const"), "hono_const"),
    (Axis("initial", "NO", label="my axis"), "my axis"),
    (Axis("injection_scale", "NO", also=[["injection_scale", "NO2"]]), "NO+NO2_x"),
    (Axis("initial", "NO", also=[["initial", "NO2"]]), "NO+NO2"),
])
def test_axis_names(axis, name):
    assert axis.name() == name


def test_grid_is_every_combination_first_axis_slowest():
    axes = [Axis("initial", "NO", [1.0, 2.0]), Axis("injection_scale", "voc", [0.5, 1.0, 2.0])]
    points = sd.grid_points(axes)
    assert len(points) == 6
    assert points[0] == {"NO_t0_ppb": 1.0, "voc_x": 0.5} and points[1] == {"NO_t0_ppb": 1.0, "voc_x": 1.0}
    assert points[-1] == {"NO_t0_ppb": 2.0, "voc_x": 2.0}
    assert sd.grid_points([]) == [{}]                                  # no axis: one run, the experiment as it is


def test_run_and_base_names():
    assert sd.run_name("day", {}) == "day"
    assert sd.run_name("day", {"NO_t0_ppb": 0.5, "my axis": 1e-05}) == "day__NO_t0_ppb-0.5__my_axis-1e-05"
    assert sd.base_name(os.path.join("x", "uv_lamps.experiment.json")) == "uv_lamps"
    assert sd.base_name(os.path.join("x", "traj_01.csv")) == "traj_01"


def test_runs_are_planned_per_experiment_and_point_or_per_table_in_a_field_study(tmp_path):
    st = Study(name="s", folder=str(tmp_path), experiments=["a.experiment.json", "b.experiment.json"],
               axes=[Axis("initial", "NO", [1.0, 2.0])])
    plans = sd.plan_runs(st)
    assert [p.name for p in plans] == ["a__NO_t0_ppb-1", "a__NO_t0_ppb-2", "b__NO_t0_ppb-1", "b__NO_t0_ppb-2"]
    assert plans[0].folder == os.path.join(str(tmp_path), "runs", "a__NO_t0_ppb-1")
    assert plans[0].cs_path.endswith("a__NO_t0_ppb-1.cs") and plans[0].info_path.endswith(".run.json")
    st.field_case, st.tables = "case.field.json", ["t1.csv", "t2.csv", "t3.csv"]
    assert sd.bases(st) == ["t1.csv", "t2.csv", "t3.csv"] and len(sd.plan_runs(st)) == 6


def test_guards_on_the_number_of_runs():
    assert sd.check_guards(50)[0] == "ok" and sd.check_guards(51)[0] == "warn"
    assert sd.check_guards(200)[0] == "warn" and sd.check_guards(201)[0] == "confirm"


def test_study_file_round_trip(tmp_path):
    st = Study(name="iso", experiments=["a.experiment.json"], modifications=[Modification("chamber", "hono_const", "1")],
               axes=[Axis("initial", "NO", [1.0, 2.0], also=[["initial", "NO2"]])], notes="n")
    path = sd.save_study(st, str(tmp_path / "iso.study.json"))
    back = sd.load_study(path)
    assert back.folder == str(tmp_path)                                # a study without folder lives next to its file
    back.folder = ""
    assert back == st


# ── changes applied to an experiment ─────────────────────────────────────────

def experiment():
    spec = ExperimentSpec()
    spec.precursors = ["C2H5OH"]
    spec.injections = [InjectionEvent("NO", "NO", "2", "09:00", 1.0, "QNO", "qNO", action="no"),
                       InjectionEvent("NO2", "NO2", "8", "09:00", 1.0, "QNO2", "qNO2", action="no2"),
                       InjectionEvent("ethanol", "C2H5OH", "200", "10:30", 2.0, action="voc",
                                      extra_channels=[("CH3CHO", "1.5")]),
                       InjectionEvent("O3", "O3", "40+5", "12:00", 20.0, "QO3", "qO3", action="o3")]
    return spec


def test_value_takes_the_type_of_the_field_it_replaces():
    assert sd._coerce(True, "no") is False and sd._coerce(False, "ON") is True
    assert sd._coerce(3, "4.0") == 4 and isinstance(sd._coerce(3, "4.0"), int)
    assert sd._coerce(298.0, "310") == 310.0 and sd._coerce("1.", "2*x") == "2*x"


def test_scaling_a_number_or_an_expression():
    assert sd._scale_expr("200", 0.5) == "100"
    assert sd._scale_expr("2.6*0.7", "0.5") == "(2.6*0.7)*0.5"         # an expression is kept and multiplied
    assert sd._scale_expr("", 3) == "0"


@pytest.mark.parametrize("target, expected", [
    ("NO|11:00|5", ("NO", "11:00", 5.0)),
    ("o3|11:00", ("O3", "11:00", 20.0)),                               # the SAPHIR default for ozone
    ("C5H8|9:05|", ("C5H8", "9:05", 1.0)),
    ("C5H8|soon|x", ("C5H8", "", 1.0)),
    ("", ("", "", 1.0)),
])
def test_new_injection_target(target, expected):
    assert sd.parse_new_injection(target) == expected


def test_injections_are_found_by_name_group_or_product():
    spec = experiment()
    names = lambda target: [e.name for e in sd._events_matching(spec, target)]
    assert names("ethanol") == ["ethanol"] and names("ETHANOL") == ["ethanol"]
    assert names("all") == ["NO", "NO2", "ethanol", "O3"]
    assert names("nox") == ["NO", "NO2"]
    assert names("o3") == ["O3"]
    assert names("ethanol, NO, ethanol") == ["ethanol", "NO"]           # a list, each once
    assert names("C2H5OH") == ["ethanol"]                              # by what it injects
    assert names("limonene") == []


def test_voc_group_takes_single_compounds_and_injected_solutions():
    """Decision of 24 Nov 2026.  Before, 'voc' left the solutions out: the action of that name was matched first
    and the line that adds the mixtures was never reached, while the study window's VOC group had both."""
    spec = experiment()
    spec.injections += [InjectionEvent("VCP", "0.6*LIMONENE+0.4*C2H5OH", "12", "11:00", 2.0, action="mixture"),
                        InjectionEvent("cooking", "0.5*NC12H26+0.5*ACR", "4", "11:30", 2.0, action="mixture")]
    names = lambda target: [e.name for e in sd._events_matching(spec, target)]
    assert names("voc") == ["ethanol", "VCP", "cooking"] and names("VOC") == names("voc")
    assert names("mixture") == ["VCP", "cooking"]                       # the injected solutions alone
    assert names("compound") == ["ethanol"]                             # the single compounds alone
    sd.apply_modification(spec, Modification("injection_scale", "mixture", "2"))
    vcp = spec.injections[4]
    assert vcp.total_ppb == "24" and vcp.products == "0.6*LIMONENE+0.4*C2H5OH"      # the whole solution, same composition
    assert spec.injections[2].total_ppb == "200"                        # ethanol untouched


def test_rate_replaces_an_earlier_adjustment_of_the_same_reaction():
    spec = experiment()
    spec.reaction_adjustments = [("C2H5OH + OH --> CH3CHO", "1e-12", "old")]
    notes = sd.apply_modification(spec, Modification("rate", "c2h5oh+oh-->ch3cho", "0"))
    assert notes == [] and spec.reaction_adjustments == [("c2h5oh+oh-->ch3cho", "0", "study variant")]


def test_rate_factor_multiplies_the_adjusted_rate_or_says_that_the_reaction_is_unknown(monkeypatch):
    spec = experiment()
    spec.reaction_adjustments = [("A + B --> C", "KMT01", "")]
    sd.apply_modification(spec, Modification("rate_factor", "A+B-->C", "0.8"))
    assert spec.reaction_adjustments == [("A+B-->C", "(KMT01)*0.8", "study: × 0.8")]
    monkeypatch.setattr(sd, "mcm_rate_of", lambda reaction, fac="": None)
    notes = sd.apply_modification(spec, Modification("rate_factor", "X --> Y", "2"))
    assert "not found" in notes[0] and len(spec.reaction_adjustments) == 1
    monkeypatch.setattr(sd, "mcm_rate_of", lambda reaction, fac="": "1.2D-11")
    sd.apply_modification(spec, Modification("rate_factor", "X --> Y", "2"))
    assert spec.reaction_adjustments[-1][:2] == ("X --> Y", "(1.2D-11)*2")


def test_addon_and_base_mechanism_ask_for_a_new_subset():
    spec = experiment()
    spec.mechanism_stem = "ethanol"
    sd.apply_modification(spec, Modification("addon", "NOAA_BB"))
    sd.apply_modification(spec, Modification("addon", "NOAA_BB"))      # twice: listed once
    assert spec.addons == ["NOAA_BB"] and spec.mechanism_stem == ""
    spec.mechanism_stem = "ethanol"
    sd.apply_modification(spec, Modification("base_mechanism", "MCM_v3.2"))
    assert spec.base_mechanism == "MCM_v3.2" and spec.mechanism_stem == ""


def test_chamber_and_condition_values():
    spec = experiment()
    assert sd.apply_modification(spec, Modification("chamber_scale", "hono_const", "0.5")) == []
    assert spec.chamber.hono_const == "(2.6*0.7)*0.5"
    sd.apply_modification(spec, Modification("chamber_scale", "walldep_value", "3"))
    assert spec.chamber.walldep_value == "3"                           # '1.' is a number: multiplied out
    sd.apply_modification(spec, Modification("chamber", "hcho_const", "0"))
    sd.apply_modification(spec, Modification("condition", "temperature_K", "310"))
    assert spec.chamber.hcho_const == "0" and spec.conditions.temperature_K == 310.0
    for kind in ("chamber", "chamber_scale", "condition"):
        assert "unknown" in sd.apply_modification(spec, Modification(kind, "no_such_field", "1"))[0]


def test_injected_amounts_are_scaled_with_their_side_channels_or_set():
    spec = experiment()
    sd.apply_modification(spec, Modification("injection_scale", "ethanol", "2"))
    ethanol = spec.injections[2]
    assert ethanol.total_ppb == "400" and ethanol.extra_channels == [("CH3CHO", "3")]
    sd.apply_modification(spec, Modification("injection_scale", "O3", "2"))
    assert spec.injections[3].total_ppb == "(40+5)*2"
    sd.apply_modification(spec, Modification("injection_set", "nox", "1.5"))
    assert [e.total_ppb for e in spec.injections[:2]] == ["1.5", "1.5"]
    assert "no injection matches" in sd.apply_modification(spec, Modification("injection_scale", "limonene", "2"))[0]


def test_initial_value_replaces_an_earlier_one_of_that_species():
    spec = experiment()
    spec.initials = [("no", "M*1e-9*3"), ("O3", "M*1e-9*10")]
    sd.apply_modification(spec, Modification("initial", "NO", "0.5"))
    assert spec.initials == [("O3", "M*1e-9*10"), ("NO", "M*1e-9*0.5")]


def test_unknown_kind_is_a_note_not_a_crash():
    assert sd.apply_modification(experiment(), Modification("magic", "x", "1")) == ["unknown modification 'magic'"]


def test_new_inorganic_injection_gets_the_saphir_column_or_its_own_when_the_name_is_taken():
    spec = ExperimentSpec()
    assert sd.add_injection(spec, "CO|11:00|", "500") == []
    co = spec.injections[0]
    assert (co.name, co.products, co.total_ppb, co.start, co.duration_min, co.column, co.switch, co.action) == \
        ("CO", "CO", "500", "11:00", 1.0, "QCO", "qCO", "co")
    spec = experiment()
    sd.add_injection(spec, "NO|11:00|5", "3")
    added = spec.injections[-1]
    assert (added.name, added.column, added.switch, added.duration_min) == ("NO_add", "Q_NO_add", "q_NO_add", 5.0)


def test_new_injection_is_set_not_added_again_at_the_next_grid_point():
    spec = experiment()
    sd.add_injection(spec, "NO|11:00|5", "3")
    sd.add_injection(spec, "NO|11:00|5", "6")
    assert len(spec.injections) == 5 and spec.injections[-1].total_ppb == "6"
    sd.add_injection(spec, "NO|13:00|5", "1")                           # another time: another injection
    assert len(spec.injections) == 6


def test_new_injection_of_a_compound_outside_the_chemistry_is_said_and_a_bad_target_refused():
    spec = experiment()
    assert "not among the experiment's compounds" in sd.add_injection(spec, "C5H8|11:00|2", "5")[0]
    assert spec.injections[-1].column == "Q_C5H8" and spec.injections[-1].action == "voc"
    assert sd.add_injection(spec, "C2H5OH|11:00|2", "5") == []          # ethanol is among them
    before = len(spec.injections)
    assert "give species|HH:MM|minutes" in sd.add_injection(spec, "C5H8", "5")[0] and len(spec.injections) == before


def test_a_grid_point_is_applied_to_every_member_of_its_axis():
    spec = experiment()
    axis = Axis("injection_set", "NO", [1.0], also=[["injection_set", "NO2"]])
    assert sd.apply_point(spec, [axis], {axis.name(): 4.0}) == []
    assert [e.total_ppb for e in spec.injections[:2]] == ["4", "4"]


# ── checks before anything is written ────────────────────────────────────────

def test_empty_study_and_missing_files_are_errors(tmp_path):
    messages = sd.validate_study(Study(folder=""))
    assert ("error", "No base experiment chosen.") in messages and ("error", "No study folder.") in messages
    st = Study(folder=str(tmp_path), experiments=[str(tmp_path / "gone.experiment.json")],
               axes=[Axis("initial", "NO", [])])
    texts = [text for level, text in sd.validate_study(st) if level == "error"]
    assert any("not found" in t for t in texts) and any("has no values" in t for t in texts)


def test_many_runs_are_a_warning(tmp_path):
    st = Study(folder=str(tmp_path), experiments=[str(tmp_path / "gone.experiment.json")],
               axes=[Axis("initial", "NO", list(range(1, 9))), Axis("initial", "NO2", list(range(1, 9)))])
    assert ("warning", "64 runs — above 50; check the estimate before starting.") in sd.validate_study(st)


# ── finished runs ────────────────────────────────────────────────────────────

def finished(plan, exit_code=0, **info):
    os.makedirs(plan.folder, exist_ok=True)
    open(plan.json_path, "w").write("{}")
    open(plan.cs_path, "w").write("cs")
    now = time.time()
    os.utime(plan.json_path, (now - 10, now - 10))
    sd._write_info(plan, {"exit": exit_code, **info})


def test_a_run_is_done_only_with_a_result_newer_than_its_input_and_exit_zero(tmp_path):
    plan = RunPlan("a", "a", {}, "a", str(tmp_path / "a"))
    assert not sd.is_done(plan)
    finished(plan)
    assert sd.is_done(plan)
    sd._write_info(plan, {"exit": 1})
    assert not sd.is_done(plan)                                         # the solver failed
    sd._write_info(plan, {})
    assert not sd.is_done(plan)                                         # no exit code recorded
    sd._write_info(plan, {"exit": 0})
    later = time.time() + 10
    os.utime(plan.json_path, (later, later))
    assert not sd.is_done(plan)                                         # inputs written again after the run


def test_progress_is_the_last_fraction_in_the_log(tmp_path):
    log = tmp_path / "run.log"
    log.write_text("start\nFraction Done: 0.25\nmore\nFraction Done: 0.5\n", encoding="utf-8")
    assert sd._tail_fraction(str(log)) == 0.5
    log.write_text("nothing yet", encoding="utf-8")
    assert sd._tail_fraction(str(log)) == 0.0 and sd._tail_fraction(str(tmp_path / "none.log")) == 0.0


def test_results_table_and_file(tmp_path):
    st = Study(name="s", folder=str(tmp_path), experiments=["a.experiment.json"], axes=[Axis("initial", "NO", [1.0, 2.0])],
               metrics=["o3_max", "nox_end"])
    plans = sd.plan_runs(st)
    finished(plans[0], elapsed=12.5, metrics={k: 1.0 for k in sd.METRICS} | {"o3_max": 61.5, "nox_end": 0.2})
    rows = sd.collect_results(st, plans)
    assert rows[0] == {"run": "a__NO_t0_ppb-1", "experiment": "a", "NO_t0_ppb": 1.0, "status": "done", "elapsed_s": 12.5,
                       "cs": plans[0].cs_path, "o3_max": 61.5, "nox_end": 0.2}
    assert rows[1]["status"] == "planned" and rows[1]["o3_max"] is None and rows[1]["cs"] == ""
    path = sd.write_results(st, rows)
    lines = open(path, encoding="utf-8").read().splitlines()
    assert lines[0] == "run,experiment,NO_t0_ppb,status,elapsed_s,cs,o3_max,nox_end" and len(lines) == 3


# ── run-time estimates ───────────────────────────────────────────────────────

def test_estimate_of_the_remaining_time():
    assert sd.parallel_factor(1) == 1.0 and sd.parallel_factor(11) == pytest.approx(2.8)
    assert sd.estimate_seconds(10, 60.0, 1) == 600.0
    assert sd.estimate_seconds(10, 60.0, 4, alone=False) == 180.0       # three rounds of four
    assert sd.estimate_seconds(10, 60.0, 4) == pytest.approx(180.0 * 1.54)
    assert sd.estimate_seconds(2, 60.0, 8, alone=False) == 60.0         # never more workers than runs
    assert sd.estimate_seconds(0, 60.0, 4) == 0.0 and sd.estimate_seconds(5, float("nan"), 4) == 0.0


def test_durations_are_written_for_people():
    assert [sd.fmt_duration(s) for s in (0, float("nan"), 42.4, 125, 3600, 5400)] == \
        ["—", "—", "42 s", "2 min 05 s", "1 h 00 min", "1 h 30 min"]


def test_prior_run_time_is_the_median_of_the_experiments_that_have_run(tmp_path):
    paths = []
    for name, elapsed in (("a", 10.0), ("b", 30.0), ("c", 0.0), ("d", None)):
        paths.append(str(tmp_path / f"{name}.experiment.json"))
        if elapsed is not None:
            (tmp_path / f"{name}.run.json").write_text(json.dumps({"elapsed": elapsed}), encoding="utf-8")
    assert sd.prior_seconds(Study(experiments=paths)) == 20.0
    assert sd.prior_seconds(Study(experiments=[paths[3]])) == 0.0


# ── results from the series of a run ─────────────────────────────────────────

M = 2.5e19                       # molecules cm-3: 1 ppb = 2.5e10
PPB = 2.5e10
T0 = (datetime(2024, 7, 9, 10, 0, tzinfo=timezone.utc) - datetime(2000, 1, 1, tzinfo=timezone.utc)).total_seconds()
TIME = T0 + 60.0 * np.arange(4)  # 10:00 … 10:03
LIT = np.array([False, True, True, False])


def series(**columns):
    d = {"M": np.full(4, M), "T": np.full(4, 300.0)}
    d.update({k: np.asarray(v, float) for k, v in columns.items()})
    return d


def derived(key, d, lit=LIT):
    return sd._derived(key, TIME, d, lit, d.get("M"), "no.cs")


def test_oh_exposure_is_the_time_integral_of_oh():
    assert derived("oh_exposure", series(OH=[1e6, 1e6, 2e6, 2e6])) == pytest.approx((1e6 + 2e6 + 2e6) * 60.0)
    assert math.isnan(derived("oh_exposure", series()))


def test_ho2_to_ro2_ratio_uses_the_lit_part_only():
    d = series(HO2=[9e9, 2e8, 4e8, 9e9], RO2=[1.0, 1e8, 2e8, 1.0])
    assert derived("ho2_ro2", d) == pytest.approx(2.0)
    assert math.isnan(derived("ho2_ro2", series(HO2=[1, 1, 1, 1], RO2=[0, 0, 0, 0])))


def test_ozone_production_from_ho2_and_no():
    # 1 ppb NO and 1e8 HO2 at 300 K: k = 3.45e-12·e^0.9 = 8.486e-12 → 8.486e-4 ppb/s = 3.055 ppb/h
    d = series(HO2=[1e8] * 4, NO=[PPB] * 4)
    assert derived("po3_mean_lit", d) == pytest.approx(3.055, rel=1e-3)
    assert derived("po3_lit_int", d) == pytest.approx(2 * 60 * 8.486e-4, rel=1e-3)      # two lit minutes
    with_ro2 = derived("po3_mean_lit", series(HO2=[1e8] * 4, NO=[PPB] * 4, RO2=[1e8] * 4))
    assert with_ro2 == pytest.approx(3.055 + 2.7e-12 * math.exp(1.2) * 1e8 * 3600, rel=1e-3)
    assert math.isnan(derived("po3_mean_lit", series(HO2=[1e8] * 4)))                   # no NO in the run


def test_kmt08_is_the_mcm_value_and_stays_below_its_high_pressure_limit():
    k = sd._kmt08(np.array([300.0]), np.array([2.46e19]))[0]
    assert k == pytest.approx(9.75e-12, rel=0.01)                      # worked out by hand from the MCM expression
    low, high = sd._kmt08(np.array([300.0, 300.0]), np.array([2.46e17, 2.46e21]))
    assert low < k < high < 3.0e-11


def test_regime_indicator_scales_as_ho2_squared_over_oh_times_no2():
    base = dict(HO2=[1e8] * 4, OH=[1e6] * 4, NO2=[PPB] * 4, H2O=[0.0] * 4)
    r = derived("regime", series(**base))
    assert r > 0
    assert derived("regime", series(**dict(base, HO2=[2e8] * 4))) == pytest.approx(4 * r)
    assert derived("regime", series(**dict(base, NO2=[2 * PPB] * 4))) == pytest.approx(r / 2)
    assert derived("regime", series(**dict(base, H2O=[0.02 * M] * 4))) > r              # water speeds up HO2 + HO2
    assert math.isnan(derived("regime", series(**dict(base, OH=[0.0] * 4))))            # no HNO3 made


def test_ozone_production_efficiency_counts_nitrogen_atoms(monkeypatch):
    monkeypatch.setattr(sd, "_nitrogen_species", lambda cs: ["HNO3", "N2O5", "ABSENT"])
    d = series(O3=[0, 10, 30, 30], NO2=[5, 5, 5, 5], HNO3=[0, 1, 3, 4], N2O5=[0, 0, 1, 1])
    # while lit: Ox 15 → 35, NOz 1 → 3 + 2·1 = 5
    assert derived("ope", d) == pytest.approx(20.0 / 4.0)
    assert math.isnan(derived("ope", series(O3=[0, 10, 30, 30], NO2=[5] * 4, HNO3=[1] * 4)))   # no NOz made


def test_voc_consumed_is_relative_to_the_largest_amount_present(monkeypatch):
    monkeypatch.setattr(sd, "_precursors_for", lambda cs: ["C5H8", "APINENE", "NOT_IN_RUN"])
    d = series(C5H8=[10, 8, 4, 1], APINENE=[0, 2, 2, 1])
    assert derived("voc_consumed", d) == pytest.approx(80.0)            # 10 at most, 2 at the end
    monkeypatch.setattr(sd, "_precursors_for", lambda cs: [])
    assert math.isnan(derived("voc_consumed", d))


def test_results_need_the_air_density():
    d = series(HO2=[1e8] * 4, NO=[PPB] * 4)
    assert math.isnan(sd._derived("po3_mean_lit", TIME, d, LIT, None, "no.cs"))
    assert math.isnan(derived("no_such_result", d))


def metrics(monkeypatch, keys, d):
    from juice.gui import experiment_plots
    monkeypatch.setattr(experiment_plots, "read_columns", lambda path, names: (TIME, d))
    return sd.compute_metrics("no.cs", keys)


def test_metrics_of_a_run(monkeypatch):
    d = series(HV=[0, 1, 1, 0], O3=np.array([10, 20, 50, 40]) * PPB, NO=np.array([4, 3, 2, 1]) * PPB,
               NO2=np.array([1, 2, 3, 4]) * PPB, OH=[0, 2e6, 4e6, 0], HONO=np.array([0.1, 0.4, 0.2, 0.1]) * PPB)
    m = metrics(monkeypatch, ["o3_max", "o3_end", "o3_time_max", "o3_lit_delta", "ox_lit_delta", "oh_mean_lit", "oh_max",
                              "nox_end", "no2_max", "hono_max", "pan_max", "not_a_metric"], d)
    assert m["o3_max"] == pytest.approx(50.0) and m["o3_end"] == pytest.approx(40.0)
    assert m["o3_time_max"] == pytest.approx(10.0 + 2.0 / 60.0)         # 10:02 UTC
    assert m["o3_lit_delta"] == pytest.approx(30.0) and m["ox_lit_delta"] == pytest.approx(31.0)
    assert m["oh_mean_lit"] == pytest.approx(3e6) and m["oh_max"] == 4e6
    assert m["nox_end"] == pytest.approx(5.0) and m["no2_max"] == pytest.approx(4.0) and m["hono_max"] == pytest.approx(0.4)
    assert math.isnan(m["pan_max"]) and math.isnan(m["not_a_metric"])   # not in the run / not a result


def test_a_run_without_light_column_or_without_light_counts_as_lit_throughout(monkeypatch):
    d = series(O3=np.array([10, 20, 50, 40]) * PPB)
    assert metrics(monkeypatch, ["o3_lit_delta"], d)["o3_lit_delta"] == pytest.approx(30.0)
    d["HV"] = np.zeros(4)
    assert metrics(monkeypatch, ["o3_lit_delta"], d)["o3_lit_delta"] == pytest.approx(30.0)


# ── the NOx axis, the ridge, sensitivities ───────────────────────────────────

@pytest.mark.parametrize("axes, index", [
    ([Axis("injection_scale", "voc"), Axis("initial", "NO2")], 1),
    ([Axis("injection_set", "NO,NO2"), Axis("injection_scale", "voc")], 0),
    ([Axis("injection_new", "no|11:00|5")], 0),
    ([Axis("emission_scale", "NOX")], 0),
    ([Axis("rate_factor", "NO + O3 --> NO2"), Axis("initial", "O3")], -1),
    ([Axis("injection_scale", "NO, ethanol")], 0),                      # named first
    ([Axis("injection_scale", "ethanol, NO")], -1),
    ([Axis("chamber", "hono_const", label="NOx_source")], 0),           # by its label
    ([], -1),
])
def test_which_axis_varies_nox(axes, index):
    assert sd.nox_axis(axes) == index


def rows_of(level, pairs):
    return [{"nox": x, "voc": level, "o3": z} for x, z in pairs]


def test_ridge_is_the_vertex_of_the_parabola_through_the_three_points_around_the_maximum():
    # through (1, 1), (2, 4), (3, 3): z = -2x² + 9x - 6, vertex at x = 2.25
    ridge = sd.ozone_ridge(rows_of(1.0, [(3, 3.0), (1, 1.0), (2, 4.0)]), "nox", "voc", "o3", log=False)
    assert ridge == [{"level": 1.0, "ridge": pytest.approx(2.25), "status": "inside", "n": 3}]


def test_ridge_on_a_logarithmic_nox_axis():
    ridge = sd.ozone_ridge(rows_of(1.0, [(1, 1.0), (10, 4.0), (100, 3.0)]), "nox", "voc", "o3", log=True)
    assert ridge[0]["ridge"] == pytest.approx(10 ** 1.25)


def test_levels_where_ozone_only_rises_or_only_falls_and_levels_without_enough_points():
    rows = (rows_of(1.0, [(1, 1.0), (2, 2.0), (3, 3.0)]) + rows_of(2.0, [(1, 3.0), (2, 2.0), (3, 1.0)])
            + rows_of(3.0, [(1, 1.0)]) + rows_of(4.0, [(1, 1.0), (2, float("nan")), (3, None)]))
    ridge = sd.ozone_ridge(rows, "nox", "voc", "o3", log=False)
    assert [(r["level"], r["status"], r["ridge"]) for r in ridge] == [(1.0, "above", None), (2.0, "below", None)]


def test_regime_map_gives_the_slope_along_the_nox_axis():
    x, y = np.meshgrid(np.linspace(0, 4, 5), np.linspace(0, 2, 5))
    x, y = x.ravel(), y.ravel()
    X, Y, Z, slope = sd.regime_map(x, y, 3 * x - y, nox_on_x=True, n=9)
    assert Z.shape == (9, 9) and np.nanmax(np.abs(slope - 3.0)) < 1e-9
    assert np.nanmax(np.abs(sd.regime_map(x, y, 3 * x - y, nox_on_x=False, n=9)[3] + 1.0)) < 1e-9


def test_normalised_sensitivity_is_the_exponent_of_a_power_law():
    xs = [1.0, 2.0, 4.0, 8.0]
    assert sd.normalised_sensitivity(xs, [x ** 2 for x in xs]) == pytest.approx([2.0, 2.0, 2.0])
    out = sd.normalised_sensitivity([1.0, 2.0, 2.0, 4.0], [1.0, -1.0, 3.0, None])
    assert all(math.isnan(v) for v in out)                              # negative value, same x, missing value


def test_axis_descriptions_for_the_plots(tmp_path):
    st = Study(folder=str(tmp_path))
    describe = lambda axis: sd.axis_description(st, axis)
    assert describe(Axis("initial", "NO")) == ("NO at the start [ppb]", 0.0, "NO at the start [ppb]")
    assert describe(Axis("injection_new", "NO|11:00|5")) == ("NO added at 11:00 over 5 min [ppb]", 0.0, "NO added [ppb]")
    assert describe(Axis("rate_factor", "A --> B"))[0] == "k(A --> B) × factor (1 = MCM)"
    assert describe(Axis("injection_scale", "voc")) == ("injections voc × factor (1 = as in the experiment)", 0.0, "× factor")
    assert describe(Axis("initial", "NO", also=[["initial", "NO2"]]))[0].endswith("— together with initial NO2")
