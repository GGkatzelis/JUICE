"""A study that has really run (the `small_study` fixture: a UV-lamp ethanol experiment over initial NO2 x injected
ethanol, 2 x 2): what each run folder holds, that the base experiment is untouched, what the results say about
the chemistry, and that a finished study is not run again."""
import json
import os

import numpy as np
import pytest

from juice.core import study as S
from juice.core.experiment_spec import load_spec


@pytest.fixture(scope="module")
def study(small_study):
    return S.load_study(small_study)


def plan_at(st, no2, ethanol):
    return next(p for p in S.plan_runs(st) if p.point == {"NO2_t0_ppb": no2, "ethanol_x": ethanol})


def test_every_grid_point_is_a_finished_run_in_its_own_folder(study):
    plans = S.plan_runs(study)
    assert sorted(p.name for p in plans) == ["uv_ethanol__NO2_t0_ppb-50__ethanol_x-0.5", "uv_ethanol__NO2_t0_ppb-50__ethanol_x-2",
                                              "uv_ethanol__NO2_t0_ppb-5__ethanol_x-0.5", "uv_ethanol__NO2_t0_ppb-5__ethanol_x-2"]
    for p in plans:
        assert S.is_done(p)
        for ext in (".easy", ".enz", ".json", ".cs", ".experiment.json", ".run.json"):
            assert os.path.isfile(os.path.join(p.folder, p.name + ext)), ext


@pytest.mark.parametrize("no2, ethanol", [(5.0, 0.5), (50.0, 2.0)])
def test_the_run_easy_carries_the_grid_point_and_where_it_came_from(study, no2, ethanol):
    p = plan_at(study, no2, ethanol)
    easy = open(p.json_path[:-5] + ".easy", encoding="latin-1").read()
    assert f"INITIAL[NO2] = CONST(M*1e-9*{no2:g})" in easy
    assert f"k[Q_ethanol-->C2H5OH]=CONST(M*1e-9/60*{200 * ethanol:g})" in easy     # 200 ppb in the base experiment
    assert f"; STUDY: iso  base experiment: uv_ethanol  run: {p.name}" in easy
    assert f"; STUDY: point: NO2_t0_ppb = {no2:g}, ethanol_x = {ethanol:g}" in easy
    spec = load_spec(os.path.join(p.folder, p.name + ".experiment.json"))
    assert spec.name == p.name and float(spec.injections[0].total_ppb) == 200 * ethanol


def test_the_base_experiment_is_not_touched(study):
    base = load_spec(study.experiments[0])
    assert base.injections[0].total_ppb == "200" and ("NO2", "M*1e-9*20") in [tuple(x) for x in base.initials]


def test_the_results_show_the_two_regimes(study):
    rows = {(r["NO2_t0_ppb"], r["ethanol_x"]): r for r in S.collect_results(study)}
    assert all(r["status"] == "done" and r["elapsed_s"] > 0 for r in rows.values())
    low, high = rows[(5.0, 0.5)], rows[(50.0, 0.5)]
    # little NO2: ozone does not care about four times more ethanol (NOx-limited), the indicator is high
    assert rows[(5.0, 2.0)]["o3_max"] == pytest.approx(low["o3_max"], rel=0.02)
    assert low["regime"] > 3 and rows[(5.0, 2.0)]["regime"] > 3
    # ten times the NO2: about three times the ozone, the indicator low (VOC-limited side)
    assert 2 < high["o3_max"] / low["o3_max"] < 4
    assert high["regime"] < 0.5 and rows[(50.0, 2.0)]["regime"] < 0.5
    assert rows[(50.0, 2.0)]["o3_max"] > high["o3_max"]                 # and there more ethanol does help


def test_each_run_records_its_metrics_and_run_time(study):
    info = json.load(open(plan_at(study, 50.0, 2.0).info_path, encoding="utf-8"))
    assert info["exit"] == 0 and info["elapsed"] > 0 and info["point"] == {"NO2_t0_ppb": 50.0, "ethanol_x": 2.0}
    assert set(S.METRICS) <= set(info["metrics"])
    assert info["metrics"]["no2_max"] == pytest.approx(50.0, rel=1e-6)            # the initial value is the maximum


def test_a_finished_study_is_not_run_again(study):
    runner = S.StudyRunner(study)
    before = {p.name: os.path.getmtime(p.cs_path) for p in runner.plans}
    events = []
    runner.run_all(on_event=lambda kind, plan, _r: events.append((kind, plan.name if plan else None)), workers=2)
    assert sorted(k for k, _n in events) == ["all_done"] + ["skipped"] * 4
    assert {p.name: os.path.getmtime(p.cs_path) for p in runner.plans} == before
    assert runner.remaining() == 0 and runner.estimate() == (0.0, "—") and runner.per_run_seconds() > 0
    assert S.materialise(study, runner.plans[0], runner.include_dir) == runner.plans[0].json_path   # inputs kept as they are


def test_extending_an_axis_keeps_the_finished_runs(study):
    longer = S.Study(**{**study.__dict__, "axes": [S.Axis("initial", "NO2", [5.0, 50.0, 100.0]), study.axes[1]]})
    plans = S.plan_runs(longer)
    assert len(plans) == 6 and sum(S.is_done(p) for p in plans) == 4


def test_time_series_of_all_runs_for_the_isopleth_video(study):
    data = S.time_series(study, S.plan_runs(study), step_min=5.0)
    assert sorted(data["runs"]) == sorted(p.name for p in S.plan_runs(study))
    assert data["O3"].shape == (4, len(data["hours"])) and data["hours"][0] == pytest.approx(8.0)
    assert data["hours"][-1] == pytest.approx(12.0) and len(data["hours"]) == 49
    lit = data["hv"] > 0.5
    assert not lit[data["hours"] < 9.0].any() and lit[(data["hours"] > 9.05) & (data["hours"] < 11.95)].all()
    assert np.all(np.isfinite(data["O3"])) and np.all(data["NOx"] >= 0)
    assert sorted(data["axes"]) == ["NO2_t0_ppb", "ethanol_x"]
    assert os.path.isfile(os.path.join(study.folder, "iso.timeseries.npz"))
    again = S.time_series(study, S.plan_runs(study), step_min=5.0)                # from the cache file
    assert np.array_equal(again["O3"], data["O3"]) and again["runs"] == data["runs"]
