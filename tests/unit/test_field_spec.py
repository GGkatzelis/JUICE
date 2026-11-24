"""Field cases: the spec file, the input table, unit conversion and the checks (field_spec.py)."""
import os

import numpy as np
import pytest

from juice.core import field_spec as fs


def _table(tmp_path, text, name="t.csv"):
    p = tmp_path / name
    p.write_text(text, encoding="utf-8")
    return str(p)


BASE = ("time,lat,lon,T,p,RH,jNO2,h,{extra_h}\n"
        "2021-08-15 12:00,34.1,-118.1,298,1000,50,0.008,1000,{v}\n"
        "2021-08-15 12:15,34.2,-118.0,299,1000,50,0.008,2000,{v}\n")


def _spec(path, **kw):
    kw.setdefault("chemistry", {"mechanism_stem": "x"})
    return fs.FieldSpec(table=path, **kw)


@pytest.mark.parametrize("run_type", list(fs.RUN_TYPES))
def test_each_template_is_a_valid_table(tmp_path, run_type):
    p = fs.write_template(str(tmp_path / "t.csv"), run_type)
    t, r = fs.check(_spec(p, run_type=run_type, area_km2=64.0, dilution=fs.DEFAULT_DILUTION[run_type]))
    assert r.ok, r.text()
    assert t.n == 3 and "NO2" in t.photolysis


def test_columns_land_in_their_roles(tmp_path):
    p = _table(tmp_path, "time,T,p,H2O,sz,jNO2,J_PANa,JO31D,E_no [ppb/s],b_O3,C_NO [ppt],mystery\n"
                         "2021-08-15 12:00,25,1000,1,30,0.008,1e-5,2e-5,0.1,40,500,7\n"
                         "2021-08-15 13:00,26,1000,1,30,0.008,1e-5,2e-5,0.1,40,1500,7\n")
    t = fs.read_table(p)
    assert t.conditions["t"][0] == pytest.approx(298.15)
    assert set(t.photolysis) == {"NO2", "PANa", "O31D"}
    assert t.emission_units == {"NO": "ppb/s"}
    assert t.background["O3"][0] == 40
    assert t.constraints["NO"].tolist() == [0.5, 1.5]
    assert t.ignored == ["mystery"]


def test_rows_are_sorted_with_their_values(tmp_path):
    p = _table(tmp_path, "time,T,p,RH,C_NO\n2021-08-15 13:00,300,1000,50,2\n2021-08-15 12:00,290,1000,50,1\n")
    t = fs.read_table(p)
    assert t.conditions["t"].tolist() == [290, 300] and t.constraints["NO"].tolist() == [1, 2]
    assert any("sorted" in n for n in t.notes)


def test_german_excel_csv_with_comment_header(tmp_path):
    p = _table(tmp_path, "# a comment, with commas, and more, commas\n# again, again, again\n# x, y\n# x, y\n# x, y\n"
                         "time;T;p;RH;C_NO\n2021-08-15 12:00;298,5;1000;50,5;1,25\n2021-08-15 13:00;299,5;1000;50,5;1,5\n")
    t = fs.read_table(p)
    assert t.conditions["t"].tolist() == [298.5, 299.5] and t.constraints["NO"].tolist() == [1.25, 1.5]


def test_semicolon_table_with_decimal_points(tmp_path):
    p = _table(tmp_path, "time;T;p;RH;C_NO\n2021-08-15 12:00;298.5;1000;50;1.25\n2021-08-15 13:00;299.5;1000;50;1.5\n")
    assert fs.read_table(p).constraints["NO"].tolist() == [1.25, 1.5]


def test_excel_table(tmp_path):
    pd = pytest.importorskip("pandas")
    pytest.importorskip("openpyxl")
    p = str(tmp_path / "t.xlsx")
    pd.DataFrame({"time": pd.to_datetime(["2021-08-15 12:00", "2021-08-15 13:00"]), "T": [298, 299], "p": [1000, 1000],
                  "RH": [50, 50], "C_NO [ppb]": [1.0, 2.0]}).to_excel(p, index=False)
    t = fs.read_table(p)
    assert t.n == 2 and t.constraints["NO"].tolist() == [1.0, 2.0]


def test_emission_units_agree(tmp_path):
    """The same physical source written in each unit gives the same molec cm-3 s-1."""
    T, p_hpa, h_m, area_km2 = 298.0, 1000.0, 1000.0, 64.0
    M = p_hpa / (T * fs.DENSITY_CONST)
    volume_rate = 1.0e6                                     # molec cm-3 s-1
    h_cm, area_cm2 = h_m * 100.0, area_km2 * 1e10
    flux = volume_rate * h_cm                               # molec cm-2 s-1
    per_second = flux * area_cm2                            # molec s-1
    given = {"molec/cm3/s": volume_rate, "ppb/s": volume_rate / M * 1e9, "molec/cm2/s": flux, "mol/km2/s": flux / fs.NA * 1e10,
             "mol/km2/h": flux / fs.NA * 1e10 * 3600.0, "molec/s": per_second, "mol/s": per_second / fs.NA}
    assert set(given) == set(fs.EMISSION_UNITS)
    cols = ",".join(f"E_S{i} [{u}]" for i, u in enumerate(given))
    vals = ",".join(repr(v) for v in given.values())
    p = _table(tmp_path, f"time,T,p,RH,h,lat,lon,jNO2,{cols}\n"
                         f"2021-08-15 12:00,{T},{p_hpa},50,{h_m},34,-118,0.008,{vals}\n"
                         f"2021-08-15 12:15,{T},{p_hpa},50,{h_m},34,-118,0.008,{vals}\n")
    spec = _spec(p, run_type="lagrangian", dilution="pbl", area_km2=area_km2)
    t, r = fs.check(spec)
    assert r.ok, r.text()
    for sp, rate in fs.emission_rates(t, spec).items():
        assert rate == pytest.approx(volume_rate, rel=1e-12), sp


def test_flux_follows_the_boundary_layer(tmp_path):
    p = _table(tmp_path, BASE.format(extra_h="E_NO [molec/cm2/s]", v="1e12"))
    spec = _spec(p, run_type="lagrangian", dilution="pbl")
    t, r = fs.check(spec)
    assert r.ok, r.text()
    assert fs.emission_rates(t, spec)["NO"] == pytest.approx([1e7, 5e6])


def test_mixing_ratio_units(tmp_path):
    p = _table(tmp_path, "time,T,p,RH,C_A [ppm],C_B [ppt],C_C [molec/cm3]\n"
                         "2021-08-15 12:00,298,1000,50,2,500,2.4e10\n2021-08-15 13:00,298,1000,50,2,500,2.4e10\n")
    t = fs.read_table(p)
    M = 1000 / (298 * fs.DENSITY_CONST)
    assert t.constraints["A"][0] == 2000 and t.constraints["B"][0] == 0.5
    assert t.constraints["C"][0] == pytest.approx(2.4e10 / M * 1e9)


@pytest.mark.parametrize("fragment, kwargs, text", [
    ("E_NO [mol/s]", {}, "footprint area"),
    ("E_NO [bananas]", {}, "unknown emission unit"),
])
def test_emission_unit_errors(tmp_path, fragment, kwargs, text):
    p = _table(tmp_path, BASE.format(extra_h=fragment, v="1"))
    _t, r = fs.check(_spec(p, run_type="lagrangian", dilution="pbl", **kwargs))
    assert not r.ok and any(text in e for e in r.errors), r.text()


def test_flux_without_height_is_an_error(tmp_path):
    p = _table(tmp_path, "time,T,p,RH,jNO2,lat,lon,E_NO\n2021-08-15 12:00,298,1000,50,0.008,34,-118,1\n"
                         "2021-08-15 13:00,298,1000,50,0.008,34,-118,1\n")
    _t, r = fs.check(_spec(p, run_type="plume", dilution="constant", k_dil_per_h=0.1))
    assert any("boundary-layer height" in e for e in r.errors), r.text()


def test_each_run_type_asks_for_what_it_needs(tmp_path):
    p = _table(tmp_path, "time,lat,lon,T,p,RH,jNO2\n2021-08-15 12:00,34,-118,298,1000,50,0.008\n"
                         "2021-08-15 13:00,34,-118,298,1000,50,0.008\n")
    _t, r = fs.check(_spec(p, run_type="constrained"))
    assert any("measured species" in e for e in r.errors)
    _t, r = fs.check(_spec(p, run_type="plume", dilution="table"))
    assert any("kdil" in e for e in r.errors)
    _t, r = fs.check(_spec(p, run_type="lagrangian", dilution="pbl"))
    assert any("emission" in e for e in r.errors) and any("column h" in e for e in r.errors)


def test_tracer_must_be_in_the_table(tmp_path):
    p = _table(tmp_path, "time,lat,lon,T,p,RH,jNO2,C_CO\n2021-08-15 12:00,34,-118,298,1000,50,0.008,100\n"
                         "2021-08-15 13:00,34,-118,298,1000,50,0.008,90\n")
    _t, r = fs.check(_spec(p, run_type="plume", dilution="tracer", tracer="co"))
    assert r.ok, r.text()
    _t, r = fs.check(_spec(p, run_type="plume", dilution="tracer", tracer="SO2"))
    assert any("C_SO2" in e for e in r.errors)


def test_bad_values_are_reported_not_raised(tmp_path):
    p = _table(tmp_path, "time,lat,lon,T,p,RH,jNO2,h,E_NO [molec/cm2/s],B_O3,C_NO\n"
                         "2021-08-15 12:00,34,-118,298,1000,50,0.008,-5,-1,,n/a\n"
                         "2021-08-15 12:00,34,-118,298,1000,50,0.008,500,1,40,1\n")
    t, r = fs.check(_spec(p, run_type="lagrangian", dilution="pbl"))
    text = r.text()
    for piece in ("repeat an earlier time", "height h must be above zero", "E_NO has negative", "B_O3"):
        assert piece in text, text
    assert t.unreadable == {"C_NO": 1} or "C_NO" in text


@pytest.mark.parametrize("content", ["", "just one line\n", "a,b\n"])
def test_empty_files_are_reported(tmp_path, content):
    p = _table(tmp_path, content)
    t, r = fs.check(_spec(p))
    assert t is None and not r.ok


def test_missing_file_and_unreadable_time(tmp_path):
    _t, r = fs.check(_spec(str(tmp_path / "nope.csv")))
    assert any("not found" in e for e in r.errors)
    p = _table(tmp_path, "time,T\nyesterday,298\nlater,299\n")
    t, r = fs.check(_spec(p))
    assert t is None and any("time column" in e for e in r.errors)


def test_hours_of_day_need_a_date(tmp_path):
    p = _table(tmp_path, "time,T,p,RH,C_NO\n12.0,298,1000,50,1\n13.0,298,1000,50,1\n")
    t, r = fs.check(_spec(p))
    assert t is None and any("date" in e for e in r.errors)
    t, _r = fs.check(_spec(p, date="2021-08-15"))
    assert t is not None and t.times_s[1] - t.times_s[0] == 3600


def test_repeated_species_is_an_error(tmp_path):
    p = _table(tmp_path, "time,T,p,RH,C_NO,c_no [ppt]\n2021-08-15 12:00,298,1000,50,1,1\n2021-08-15 13:00,298,1000,50,1,1\n")
    _t, r = fs.check(_spec(p))
    assert any("given twice" in e for e in r.errors)


def test_spec_round_trip_keeps_table_next_to_it(tmp_path):
    tp = fs.write_template(str(tmp_path / "t.csv"), "lagrangian")
    spec = fs.FieldSpec(name="la", run_type="lagrangian", table=tp, lat=34.1, lon=-118.1, chemistry={"mechanism_stem": "RACM2B-VCP"}, area_km2=64.0,
                        dilution="pbl", initials={"O3": 20.0})
    path = fs.save_spec(spec, str(tmp_path / "la"))
    assert path.endswith(fs.SPEC_SUFFIX)
    assert '"table": "t.csv"' in open(path, encoding="utf-8").read()
    back = fs.load_spec(path)
    assert os.path.samefile(back.table, tp)
    back.table = spec.table
    assert back == spec


def test_table_outside_the_folder_stays_absolute(tmp_path):
    (tmp_path / "a").mkdir(); (tmp_path / "b").mkdir()
    tp = fs.write_template(str(tmp_path / "b" / "t.csv"), "plume")
    path = fs.save_spec(fs.FieldSpec(table=tp), str(tmp_path / "a" / "x.field.json"))
    assert os.path.samefile(fs.load_spec(path).table, tp)


def test_spec_file_with_unknown_key_and_bad_values(tmp_path):
    p = tmp_path / "x.field.json"
    p.write_text('{"name": "x", "from_the_future": 1}', encoding="utf-8")
    assert fs.load_spec(str(p)).name == "x"
    p.write_text('{"area_km2": "big"}', encoding="utf-8")
    with pytest.raises(ValueError):
        fs.load_spec(str(p))
    p.write_text('{"version": 99}', encoding="utf-8")
    with pytest.raises(ValueError):
        fs.load_spec(str(p))
    p.write_text('[1, 2]', encoding="utf-8")
    with pytest.raises(ValueError):
        fs.load_spec(str(p))


def test_initials_text():
    assert fs.parse_initials("O3=20, no = 0.5; PAN=3e-1") == {"O3": 20.0, "NO": 0.5, "PAN": 0.3}
    assert fs.parse_initials("  ") == {}
    assert fs.format_initials({"O3": 20.0, "NO": 0.5}) == "O3=20, NO=0.5"
    for bad in ("O3", "O3=x", "O3=-1", "O3=nan"):
        with pytest.raises(ValueError):
            fs.parse_initials(bad)


def test_step_above_a_minute_warns(tmp_path):
    p = fs.write_template(str(tmp_path / "t.csv"), "constrained")
    _t, r = fs.check(_spec(p, step_s=900))
    assert r.ok and any("RO2" in w for w in r.warnings)


def _two_rows(tmp_path, h="1000,2000"):
    a, b = h.split(",")
    return _table(tmp_path, f"time,lat,lon,T,p,RH,jNO2,h,C_CO,B_CO,kdil\n2021-08-15 12:00,34,-118,298,1000,50,0.008,{a},100,20,1e-4\n"
                            f"2021-08-15 12:15,34,-118,298,1000,50,0.008,{b},60,20,2e-4\n")


def test_dilution_from_each_source(tmp_path):
    t, _r = fs.check(_spec(_two_rows(tmp_path)))
    n = t.n
    assert fs.dilution_rate(t, _spec("", dilution="none")).tolist() == [0.0] * n
    assert fs.dilution_rate(t, _spec("", dilution="constant", k_dil_per_h=0.36)) == pytest.approx([1e-4, 1e-4])
    assert fs.dilution_rate(t, _spec("", dilution="table")) == pytest.approx([1e-4, 2e-4])
    growth = np.log(2.0) / 900.0
    assert fs.dilution_rate(t, _spec("", dilution="pbl")) == pytest.approx([growth, growth])
    decay = np.log(80.0 / 40.0) / 900.0                # above the background of 20 ppb
    assert fs.dilution_rate(t, _spec("", dilution="tracer", tracer="CO")) == pytest.approx([decay, decay])


def test_collapsing_layer_uses_the_collapse_rate(tmp_path):
    t, _r = fs.check(_spec(_two_rows(tmp_path, "2000,1000")))
    assert fs.dilution_rate(t, _spec("", dilution="pbl", k_collapse_per_h=0.36)) == pytest.approx([1e-4, 1e-4])


def test_deposition_rates(tmp_path):
    t, _r = fs.check(_spec(_two_rows(tmp_path)))
    d = fs.deposition_rates(t, _spec("", deposition_velocity_cm_s={"O3": 0.5}, deposition_lifetime_h={"HNO3": 2.0}))
    assert d["O3"] == pytest.approx([0.5 / 1e5, 0.5 / 2e5])
    assert d["HNO3"] == pytest.approx([1 / 7200.0] * 2)


def test_deposition_is_checked(tmp_path):
    p = _table(tmp_path, "time,T,p,RH,sz,C_NO\n2021-08-15 12:00,298,1000,50,30,1\n2021-08-15 13:00,298,1000,50,30,1\n")
    _t, r = fs.check(_spec(p, deposition_velocity_cm_s={"O3": 0.5}))
    assert any(st == fs.STEP_MIXING and "height column h" in m for _l, st, m in r.issues), r.text()
    _t, r = fs.check(_spec(p, deposition_velocity_cm_s={"O3": 0.5}, deposition_lifetime_h={"O3": 1.0}))
    assert any("not both" in e for e in r.errors)
    _t, r = fs.check(_spec(p, deposition_lifetime_h={"O3": 0.0}))
    assert any("above zero" in e for e in r.errors)


def test_every_check_names_its_step(tmp_path):
    p = _table(tmp_path, "time,lat,lon,T,p,RH,jNO2,E_NO [bananas],B_O3,kdil\n2021-08-15 12:00,34,-118,298,1000,50,0.008,1,-1,-1\n"
                         "2021-08-15 13:00,34,-118,298,1000,50,0.008,1,40,1\n")
    _t, r = fs.check(fs.FieldSpec(table=p, name="bad name", run_type="lagrangian", dilution="pbl"))
    step = {m: st for _l, st, m in r.issues}
    assert len(r.issues) == len(r.errors) + len(r.warnings)
    assert step["no chemistry chosen yet"] == fs.STEP_CHEMISTRY
    assert step["B_O3 has negative values"] == fs.STEP_MIXING
    assert next(st for m, st in step.items() if "bananas" in m) == fs.STEP_SOURCES
    assert next(st for m, st in step.items() if "case name" in m) == fs.STEP_CASE
    assert r.worst(fs.STEP_CHEMISTRY) == "error" and r.worst(fs.STEP_RUN) == ""


def test_no_chemistry_is_an_error(tmp_path):
    p = fs.write_template(str(tmp_path / "t.csv"), "constrained")
    _t, r = fs.check(fs.FieldSpec(table=p))
    assert "no chemistry chosen yet" in r.errors
    for chem in ({"precursors": ["C5H8"]}, {"mechanism_stem": "RACM2"}, {"inorganic_only": True}):
        _t, r = fs.check(fs.FieldSpec(table=p, chemistry=chem))
        assert r.ok, r.text()


def test_observation_columns_are_read_and_not_held(tmp_path):
    p = _table(tmp_path, "time,T,p,RH,C_NO,O_HCHO [ppt],O_O3\n2021-08-15 12:00,298,1000,50,1,2500,30\n2021-08-15 13:00,298,1000,50,1,,40\n")
    t = fs.read_table(p)
    assert t.observations["HCHO"][0] == pytest.approx(2.5) and np.isnan(t.observations["HCHO"][1])
    assert t.observations["O3"].tolist() == [30.0, 40.0] and "O3" not in t.constraints
    _t, r = fs.check(_spec(p, lat=34.0, lon=-118.0))
    assert r.ok and any("2 observed for comparison" in i for i in r.info), r.text()
    p = _table(tmp_path, "time,T,p,RH,C_NO,O_NO\n2021-08-15 12:00,298,1000,50,1,1\n2021-08-15 13:00,298,1000,50,1,1\n", "u.csv")
    _t, r = fs.check(_spec(p, lat=34.0, lon=-118.0))
    assert any("both held" in w for w in r.warnings)


def test_emission_scale_is_applied_and_reported(tmp_path):
    p = _table(tmp_path, BASE.format(extra_h="E_NO [molec/cm2/s],E_C5H8 [molec/cm2/s]", v="1e12,1e12"))
    spec = _spec(p, run_type="lagrangian", dilution="pbl", emission_scale={"NOX": 0.5})
    t, r = fs.check(spec)
    assert r.ok and any("NO × 0.5" in i for i in r.info), r.text()
    rates = fs.emission_rates(t, spec)
    assert rates["NO"] == pytest.approx(0.5 * rates["C5H8"])
    _t, r = fs.check(_spec(p, run_type="lagrangian", dilution="pbl", emission_scale={"VOC": -1}))
    assert any("emission scale" in e for e in r.errors)
