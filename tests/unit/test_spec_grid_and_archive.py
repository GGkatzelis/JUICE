"""experiment_spec.py: the time grid, the roof and pulse columns, and the reading of the SAPHIR archive.  The archive
is a network drive that no test machine has, so a small one is built here: <root>/YYYY/MM/DD/<date>.<file>.nc."""
import os

import numpy as np
import pytest

from juice.core import experiment_spec as es
from juice.core.experiment_spec import Conditions

DATE = "2024-07-09"
DAY = es._day_seconds(DATE)


def conditions(root="", **kw):
    c = Conditions()
    c.date, c.start, c.end, c.step_s, c.archive_root = DATE, "10:00", "10:05", 60, str(root)
    for k, v in kw.items():
        setattr(c, k, v)
    return c


@pytest.fixture(autouse=True)
def fresh_archive_cache():
    es.clear_archive_cache()
    yield
    es.clear_archive_cache()


# ── time grid, roof, pulses ──────────────────────────────────────────────────

def test_time_helpers():
    assert es._day_seconds("2000-01-02") == 86400.0
    assert es._hhmm_s("01:30") == 5400 and es._hhmm_s(" 9:05 ") == 32700


def test_grid_runs_from_start_to_end_inclusive():
    grid = es.time_grid(conditions())
    assert len(grid) == 6 and grid[0] == DAY + 36000 and grid[-1] == DAY + 36300
    assert len(es.time_grid(conditions(step_s=120))) == 3              # 10:00, 10:02, 10:04: the end is not overshot


def test_roof_column_is_one_from_opening_up_to_but_not_at_closing():
    c = conditions(roof_open=[("10:01", "10:03"), ("10:05", "10:30")])
    assert list(es._roof(es.time_grid(c), c)) == [0, 1, 1, 0, 0, 1]


def test_pulse_has_integral_one_over_its_duration():
    c = conditions()
    grid = es.time_grid(c)
    assert list(es.pulse_column(grid, c, "10:01", 2.0)) == [0, 0.5, 0.5, 0, 0, 0]
    assert list(es.pulse_column(grid, c, "10:02", 1.0)) == [0, 0, 1.0, 0, 0, 0]
    assert list(es.pulse_column(grid, c, "10:02", 0.1)) == [0, 0, 1.0, 0, 0, 0]       # shorter than a step: one row
    assert list(es.pulse_column(grid, c, "", 2.0)) == [0] * 6                         # no start time: no pulse


def test_pulse_that_runs_past_the_end_is_cut_and_one_after_the_end_is_absent():
    c = conditions()
    grid = es.time_grid(c)
    assert list(es.pulse_column(grid, c, "10:04", 4.0)) == [0, 0, 0, 0, 0.25, 0.25]   # half of it is injected
    assert list(es.pulse_column(grid, c, "11:00", 1.0)) == [0] * 6


def test_pulse_timed_before_the_run_enters_in_the_first_step_and_the_review_says_so():
    """Decision of 24 Nov 2026: the ENZ column stays as it was (the validated files do not change); the review of
    the experiment now warns about an injection before the start and about one after the end."""
    c = conditions()
    grid = es.time_grid(c)
    assert list(es.pulse_column(grid, c, "09:00", 2.0)) == [0.5, 0.5, 0, 0, 0, 0]
    spec = es.ExperimentSpec()
    spec.conditions = c
    spec.inorganic_only = True
    spec.injections = [es.InjectionEvent("NO", "NO", "2", "09:00", 1.0, "QNO", "qNO", action="no"),
                       es.InjectionEvent("O3", "O3", "40", "12:00", 20.0, "QO3", "qO3", action="o3"),
                       es.InjectionEvent("NO2", "NO2", "5", "10:02", 1.0, "QNO2", "qNO2", action="no2")]
    warnings = [text for level, text in es.validate(spec) if level == "warning"]
    assert any("'NO' at 09:00 is before the start of the run (10:00)" in t and "first step" in t for t in warnings)
    assert any("'O3' at 12:00 is after the end of the run (10:05)" in t and "not injected" in t for t in warnings)
    assert not any("'NO2'" in t and "the run" in t for t in warnings)


# ── the archive ──────────────────────────────────────────────────────────────

def write_nc(root, suffix, variables, minutes=(0, 2, 4, 6), start="10:00"):
    """One archive file with a time axis in minutes after `start` and the given variables."""
    nc4 = pytest.importorskip("netCDF4")
    y, m, d = DATE.split("-")
    folder = os.path.join(str(root), y, m, d)
    os.makedirs(folder, exist_ok=True)
    path = os.path.join(folder, f"{DATE}.{suffix}")
    with nc4.Dataset(path, "w") as ds:
        ds.createDimension("time", len(minutes))
        t = ds.createVariable("time", "f8", ("time",))
        t.units = f"seconds since {DATE} {start}:00"
        t[:] = [60.0 * x for x in minutes]
        for name, values in variables.items():
            ds.createVariable(name, "f8", ("time",))[:] = values
    return path


def test_where_a_column_is_looked_for(tmp_path):
    sources = es.archive_sources(conditions(tmp_path), "jNO2")
    day = os.path.join(str(tmp_path), "2024", "07", "09")
    assert sources[0] == (os.path.join(day, f"{DATE}.SAPHIR.SR.jno2.nc"), ["jNO2", "jno2", "JNO2", "j_no2"])
    assert sources[1] == (os.path.join(day, f"{DATE}.saphir_photolysis.jno2.nc"), ["jNO2"])     # the naming of 2017
    assert es.archive_sources(conditions(tmp_path), "not_a_column") == []


def test_column_is_read_and_interpolated_onto_the_grid(tmp_path):
    write_nc(tmp_path, "SAPHIR.SR.jno2.nc", {"jNO2": [0.0, 0.002, 0.004, 0.006]})
    c = conditions(tmp_path)
    col = es._archive_column(c, "jNO2", es.time_grid(c), 0.0)
    assert col == pytest.approx([0.0, 0.001, 0.002, 0.003, 0.004, 0.005])


def test_grid_times_outside_the_data_get_the_fill_value(tmp_path):
    write_nc(tmp_path, "SAPHIR.collected.all_param.nc", {"TK_M_SAPHIR": [290.0, 292.0]}, minutes=(2, 4))
    c = conditions(tmp_path)
    assert list(es._archive_column(c, "t", es.time_grid(c), 284.0)) == [284.0, 284.0, 290.0, 291.0, 292.0, 284.0]


def test_variable_under_another_spelling_and_file_under_another_case(tmp_path):
    write_nc(tmp_path, "SAPHIR.SR.jNO2.nc", {"jno2": [0.001] * 4})      # the file name as some days have it
    c = conditions(tmp_path)
    assert list(es._archive_column(c, "jNO2", es.time_grid(c), 0.0)) == [0.001] * 6


def test_old_archive_names_are_used_when_the_current_file_is_missing(tmp_path):
    write_nc(tmp_path, "NOxContainer.o3_uv.nc", {"O3": [40.0] * 4})
    c = conditions(tmp_path)
    assert list(es._archive_column(c, "O3", es.time_grid(c), 0.0)) == [40.0] * 6


def test_current_file_without_the_variable_falls_through_to_the_old_one(tmp_path):
    write_nc(tmp_path, "SAPHIR.TR43.NO_NO2_NOx.nc", {"NO2": [5.0] * 4})                # no NO in it
    write_nc(tmp_path, "NOxContainer.no.nc", {"NO": [1.5] * 4})
    c = conditions(tmp_path)
    assert list(es._archive_column(c, "NO", es.time_grid(c), 0.0)) == [1.5] * 6


def test_column_that_the_day_does_not_have_is_none(tmp_path):
    c = conditions(tmp_path)
    assert es._archive_column(c, "jNO2", es.time_grid(c), 0.0) is None                 # no day folder at all
    write_nc(tmp_path, "SAPHIR.SR.jo1d.nc", {"jO1D": [1e-5] * 4})
    es.clear_archive_cache()
    assert es._archive_column(c, "jNO2", es.time_grid(c), 0.0) is None                 # folder, but not this file
    assert es._archive_column(c, "H2O", es.time_grid(c), 0.0) is None                  # no water file, no collected file


def test_water_from_relative_humidity_when_there_is_no_water_file(tmp_path):
    # 50 % RH at 20 °C and 1000 hPa: saturation 23.4 hPa (Magnus) → 0.5 · 23.4 / 1000 = 1.17 vol-%
    write_nc(tmp_path, "SAPHIR.collected.all_param.nc",
             {"RH_SAPHIR_VAISALA2": [50.0, 50.0, -2.0, 50.0], "TK_M_SAPHIR": [293.15] * 4, "p_SAPHIR": [1000.0] * 4})
    c = conditions(tmp_path)
    h2o = es._archive_column(c, "H2O", es.time_grid(c), 0.0)
    assert h2o[0] == pytest.approx(1.17, rel=0.005) and h2o[2] == pytest.approx(1.17, rel=0.005)
    assert h2o[4] == 0.0                                                # a negative reading of a dry chamber is 0
    assert 0.0 < h2o[3] < h2o[2]                                        # on the way down to it


def test_water_from_humidity_uses_the_stated_pressure_without_a_pressure_series(tmp_path):
    write_nc(tmp_path, "SAPHIR.collected.all_param.nc", {"MT401B_RH_SAPHIR": [50.0] * 4, "TK_M_SAPHIR": [293.15] * 4})
    c = conditions(tmp_path, pressure_hPa=500.0)
    assert es._archive_column(c, "H2O", es.time_grid(c), 0.0)[0] == pytest.approx(2.34, rel=0.005)


def test_water_from_humidity_needs_humidity_and_temperature(tmp_path):
    write_nc(tmp_path, "SAPHIR.collected.all_param.nc", {"TK_M_SAPHIR": [293.15] * 4})
    c = conditions(tmp_path)
    assert es._archive_column(c, "H2O", es.time_grid(c), 0.0) is None


def test_the_water_file_wins_over_the_humidity(tmp_path):
    write_nc(tmp_path, "SAPHIR.PICARRO_CFKADS2270.CO_CO2_CH4_H2O.nc", {"H2O": [0.8] * 4})
    write_nc(tmp_path, "SAPHIR.collected.all_param.nc", {"RH_SAPHIR_VAISALA2": [50.0] * 4, "TK_M_SAPHIR": [293.15] * 4})
    c = conditions(tmp_path)
    assert list(es._archive_column(c, "H2O", es.time_grid(c), 0.0)) == [0.8] * 6


def test_a_file_is_listed_and_read_once_per_session(tmp_path):
    path = write_nc(tmp_path, "SAPHIR.SR.jno2.nc", {"jNO2": [0.002] * 4})
    c = conditions(tmp_path)
    grid = es.time_grid(c)
    assert es._archive_column(c, "jNO2", grid, 0.0) is not None
    os.remove(path)
    assert es._archive_column(c, "jNO2", grid, 0.0) is not None         # still from the session's memory
    es.clear_archive_cache()
    assert es._archive_column(c, "jNO2", grid, 0.0) is None             # forgotten: the file is looked for again
