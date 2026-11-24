"""enz_builder.py, the functions without a window: the solar zenith angle, NetCDF time and variables,
interpolation onto the ENZ grid, the step, roof and pulse columns, and the ENZ writer.  The experiment
description, the field mode and the measured inputs all build their ENZ files from these."""
from datetime import datetime, timezone
from types import SimpleNamespace as NS

import numpy as np
import pytest

from juice.core import enz_builder as eb

DAY = "2024-07-09"


def seconds(text):
    """'2024-07-09 10:00' → seconds since the ENZ epoch (1 Jan 2000 UTC)."""
    dt = datetime.strptime(text, "%Y-%m-%d %H:%M").replace(tzinfo=timezone.utc)
    return (dt - eb.ENZ_EPOCH).total_seconds()


def minute_grid(start, n):
    return seconds(start) + 60.0 * np.arange(n)


def test_the_epoch_helper_of_this_file():
    assert seconds("2000-01-01 00:00") == 0.0 and seconds("2000-01-02 00:01") == 86460.0


# ── solar zenith angle ───────────────────────────────────────────────────────

def whole_day(date):
    return seconds(f"{date} 00:00") + 60.0 * np.arange(1440)


# The routine takes one declination per calendar day (that of 00 UTC).  At the solstices this does not matter; on
# 20 March 2023 the sun crossed the equator at 21:24 UTC, the routine uses the declination of midnight before,
# -0.4°, where the true value at midday was -0.15°: the zenith angle is up to about 0.4° off near the equinoxes.
@pytest.mark.parametrize("date, declination", [("2023-06-21", 23.44), ("2023-12-21", -23.44), ("2023-03-20", -0.4)])
def test_noon_and_midnight_zenith_follow_latitude_and_declination(date, declination):
    sza = eb.compute_solar_zenith(whole_day(date))
    assert sza.min() == pytest.approx(eb.SAPHIR_LAT - declination, abs=0.2)             # sun highest
    assert sza.max() == pytest.approx(180.0 - (eb.SAPHIR_LAT + declination), abs=0.2)   # sun lowest


def test_solar_noon_in_juelich_is_before_12_utc():
    # 6.36° E is 25.5 min ahead of Greenwich; the equation of time on 21 June is about −2 min → about 11:36 UTC
    sza = eb.compute_solar_zenith(whole_day("2023-06-21"))
    assert abs(int(sza.argmin()) - (11 * 60 + 36)) <= 3


def test_zenith_angle_moves_with_the_site():
    noon = [seconds("2023-03-20 12:00")]
    equator = eb.compute_solar_zenith(noon, lat=0.0, lon=0.0)[0]
    assert equator < 3.0                                              # overhead at the equator at the equinox
    assert eb.compute_solar_zenith(noon, lat=60.0, lon=0.0)[0] == pytest.approx(equator + 60.0, abs=2.0)
    # six hours of longitude to the east it is evening at the same UTC time
    assert eb.compute_solar_zenith(noon, lat=0.0, lon=90.0)[0] == pytest.approx(90.0, abs=3.0)


def test_zenith_angle_is_symmetric_about_solar_noon():
    sza = eb.compute_solar_zenith(whole_day("2023-06-21"))
    k = int(sza.argmin())
    assert sza[k - 120] == pytest.approx(sza[k + 120], abs=0.3)


# ── NetCDF time and variables ────────────────────────────────────────────────

def test_time_units_are_converted_to_enz_seconds():
    convert = eb._convert_nc_time_to_enz_seconds
    assert list(convert([0, 60], "seconds since 2000-01-01 00:00:00 UTC")) == [0.0, 60.0]
    assert list(convert([0, 1], "days since 2000-01-02")) == [86400.0, 172800.0]
    assert list(convert([1.5], "hours since 2000-1-1 01:00")) == [9000.0]
    assert list(convert([2], "Minutes since 2000/01/01 00:00:30")) == [150.0]
    assert convert([0], "seconds since 2024-07-09 00:00:00")[0] == seconds("2024-07-09 00:00")


def test_time_without_a_reference_is_taken_as_enz_seconds():
    assert list(eb._convert_nc_time_to_enz_seconds([5, 6], "s")) == [5.0, 6.0]


@pytest.fixture
def nc_file(tmp_path):
    nc4 = pytest.importorskip("netCDF4")
    path = str(tmp_path / "day.nc")
    with nc4.Dataset(path, "w") as ds:
        ds.createDimension("time", 4)
        ds.createDimension("short", 3)
        t = ds.createVariable("time", "f8", ("time",))
        t.units = "seconds since 2024-07-09 00:00:00"
        t[:] = [0, 60, 120, 180]
        no2 = ds.createVariable("no2", "f8", ("time",), fill_value=-999.0)     # the lower-case spelling
        no2[:] = np.ma.masked_values([1.0, -999.0, 3.0, 4.0], -999.0)
        ds.createVariable("O3", "f8", ("short",))[:] = [10.0, 20.0, 30.0]
    return path


def test_variable_is_read_with_enz_times_and_gaps_as_nan(nc_file):
    times, values = eb.read_nc_variable(nc_file, "NO2")               # found under its alternative name 'no2'
    assert times[0] == seconds("2024-07-09 00:00") and times[3] - times[0] == 180.0
    assert values[0] == 1.0 and np.isnan(values[1]) and list(values[2:]) == [3.0, 4.0]


def test_time_and_values_are_cut_to_the_shorter_of_the_two(nc_file):
    times, values = eb.read_nc_variable(nc_file, "O3")
    assert len(times) == len(values) == 3


def test_missing_variable_lists_what_the_file_has(nc_file):
    with pytest.raises(KeyError, match="'jNO2' not found.*no2"):
        eb.read_nc_variable(nc_file, "jNO2")


def test_missing_file_and_file_without_time(tmp_path):
    nc4 = pytest.importorskip("netCDF4")
    with pytest.raises(FileNotFoundError):
        eb.read_nc_variable(str(tmp_path / "nothing.nc"), "NO2")
    path = str(tmp_path / "no_time.nc")
    with nc4.Dataset(path, "w") as ds:
        ds.createDimension("n", 2)
        ds.createVariable("NO2", "f8", ("n",))[:] = [1.0, 2.0]
    with pytest.raises(KeyError, match="No time variable"):
        eb.read_nc_variable(path, "NO2")


# ── onto the ENZ grid ────────────────────────────────────────────────────────

def test_interpolation_is_linear_and_fills_outside_the_data():
    out = eb.interpolate_to_grid(np.array([10.0, 20.0]), np.array([1.0, 3.0]), np.array([0.0, 10.0, 15.0, 20.0, 30.0]),
                                 fill_value=-1.0)
    assert list(out) == [-1.0, 1.0, 2.0, 3.0, -1.0]


def test_interpolation_bridges_gaps_in_the_data():
    out = eb.interpolate_to_grid(np.array([0.0, 10.0, 20.0]), np.array([1.0, np.nan, 3.0]), np.array([10.0]))
    assert list(out) == [2.0]                                         # straight line across the missing value


def test_interpolation_of_data_without_one_valid_value_is_the_fill_value():
    out = eb.interpolate_to_grid(np.array([0.0, 10.0]), np.array([np.nan, np.nan]), np.arange(3.0), fill_value=284.0)
    assert list(out) == [284.0, 284.0, 284.0]


def test_step_function_holds_each_value_until_the_next():
    grid = minute_grid(f"{DAY} 11:47", 9)                             # 11:47 … 11:55
    out = eb.step_to_grid("11:54,0; 11:49,0.2", grid, DAY)            # unsorted on purpose
    assert list(out) == [0, 0, 0.2, 0.2, 0.2, 0.2, 0.2, 0, 0]


def test_step_function_skips_what_it_cannot_read():
    grid = minute_grid(f"{DAY} 10:00", 3)
    assert list(eb.step_to_grid("10:01,5; nonsense; 10:xx,7; ;10:02", grid, DAY)) == [0, 5, 5]
    assert list(eb.step_to_grid("", grid, DAY)) == [0, 0, 0]


def test_roof_is_open_from_the_opening_minute_up_to_but_not_at_the_closing_minute():
    grid = minute_grid(f"{DAY} 08:58", 6)                             # 08:58 … 09:03
    out = eb.hv_intervals_to_grid([("09:00", "09:02")], grid, DAY)
    assert list(out) == [0, 0, 1, 1, 0, 0]


def test_several_roof_intervals_and_an_unreadable_one():
    grid = minute_grid(f"{DAY} 09:00", 6)
    out = eb.hv_intervals_to_grid([("09:00", "09:01"), ("open", "close"), ("09:03", "09:05")], grid, DAY)
    assert list(out) == [1, 0, 0, 1, 1, 0]


def test_roof_interval_of_another_day_leaves_the_grid_closed():
    grid = minute_grid("2024-07-10 09:00", 3)
    assert list(eb.hv_intervals_to_grid([("09:00", "09:02")], grid, DAY)) == [0, 0, 0]


def test_pulses_add_up_where_they_overlap():
    grid = minute_grid(f"{DAY} 10:00", 6)
    out = eb.pulse_to_grid([("10:01", "10:04", 2.0), ("10:03", "10:05", 0.5)], grid, DAY)
    assert list(out) == [0, 2.0, 2.0, 2.5, 0.5, 0]


def test_pulse_that_ends_before_it_starts_or_cannot_be_read_is_left_out():
    grid = minute_grid(f"{DAY} 10:00", 4)
    out = eb.pulse_to_grid([("10:02", "10:01", 9.0), ("10:02", "10:02", 9.0), (None, "10:03", 9.0),
                            ("10:01", "10:02", "much"), ("10:03", "10:04", "1.5")], grid, DAY)
    assert list(out) == [0, 0, 0, 1.5]


# ── the ENZ file ─────────────────────────────────────────────────────────────

def test_enz_file_header_columns_and_number_formats(tmp_path):
    path = str(tmp_path / "kammer.enz")
    eb.write_enz(path, np.array([772000000.0, 772000060.0]),
                 {"hv": np.array([0.0, 1.0]), "jNO2": np.array([0.0052, float("nan")]), "t": np.array([-3.5, float("inf")])})
    lines = open(path).read().splitlines()
    assert lines[:6] == ["COLUMN 1=TIME / SEC SINCE 1.1.2000 00:00:00", "COLUMN 2=hv", "COLUMN 3=jNO2", "COLUMN 4=t",
                         "NUMBER OF COLUMNS=4", "&&&&&&&&&&&&&&&&&&&&&"]
    assert lines[6] == " 772000000.00000000, 0, 5.20000000E-03, -3.50000000"       # small numbers in E notation
    assert lines[7] == " 772000060.00000000, 1.00000000, 0.00000000, 0.00000000"   # NaN and infinity become zero
    assert len(lines) == 8


def test_written_enz_reads_back(tmp_path):
    path = str(tmp_path / "k.enz")
    grid = minute_grid(f"{DAY} 10:00", 5)
    columns = {"sz": eb.compute_solar_zenith(grid), "jO1D": np.array([0, 1.2e-5, 2.4e-5, 1e-7, 0.0])}
    eb.write_enz(path, grid, columns)
    rows = np.array([[float(x) for x in line.split(",")] for line in open(path).read().splitlines()[5:]])
    assert np.array_equal(rows[:, 0], grid)
    assert rows[:, 1] == pytest.approx(columns["sz"], abs=1e-8)
    assert rows[:, 2] == pytest.approx(columns["jO1D"], rel=1e-8)


# ── the rest ─────────────────────────────────────────────────────────────────

def test_chamber_volume_is_saphirs_by_default_and_for_an_unknown_name():
    assert eb._chamber_volume() == 276.0
    assert eb._chamber_volume("SAPHIR") == 276.0


def test_wiki_action_log_lines_are_found_in_the_page(monkeypatch):
    requests = pytest.importorskip("requests")
    page = ("<html><body><h1>2024-07-09</h1>\n<pre>\n"
            "09:00 R  <b>roof</b> open\n"
            "10:15 X 20 &amp; 5 µl   limonene&nbsp;injected\n"
            "not an action line\n"
            "10:30 XY two letters are not a code\n"
            "14:30 r roof closed\n</pre></body></html>")
    seen = {}

    def fake_get(url, **kwargs):
        seen["url"] = url
        return NS(text=page, raise_for_status=lambda: None)

    monkeypatch.setattr(requests, "get", fake_get)
    events = eb.fetch_wiki_events("2024-07-09")
    assert seen["url"].endswith("/2024-07-09")
    assert events == [{"time": "09:00", "code": "R", "description": "roof open"},
                      {"time": "10:15", "code": "X", "description": "20 & 5 µl limonene injected"},
                      {"time": "14:30", "code": "r", "description": "roof closed"}]
