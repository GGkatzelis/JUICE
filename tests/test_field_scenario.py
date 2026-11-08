"""Field case → ENZ + EASY → scenario (field_scenario.py), and the solver's held species.

The quick tests check the generated files and the parsed scenario.  The slow one runs the solver
on a dark night and compares with analytic solutions: emission + dilution towards a background,
a species held to a ramp, first-order deposition."""
import json
import os

import numpy as np
import pytest

import field_scenario as fsc
import field_spec as fs
from scenario_prep import load_paths, solver_from_toml

INCLUDE_DIR = load_paths("")["include_dir"]
INORGANIC = {"inorganic_only": True, "base_mechanism": "MCM_v3.3.1"}


def night_table(path, rows=13, h="500"):
    """15-min rows, 2021-08-15 06:00 UTC onwards (night in Los Angeles): T, p, RH, h, jNO2 = 0, an emitted and
    diluted CO with a background, NO held to a ramp."""
    lines = ["time,lat,lon,T,p,RH,h,jNO2,E_CO [molec/cm3/s],B_CO [ppb],C_NO [ppb]"]
    for i in range(rows):
        hh, mm = divmod(6 * 60 + 15 * i, 60)
        lines.append(f"2021-08-15 {hh:02d}:{mm:02d},34.14,-118.13,298,1000,50,{h},0,1e7,100,{1 + 4 * i / (rows - 1):.4f}")
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    return str(path)


def night_spec(table, **kw):
    kw.setdefault("name", "night"); kw.setdefault("run_type", "plume"); kw.setdefault("chemistry", INORGANIC)
    kw.setdefault("dilution", "constant"); kw.setdefault("k_dil_per_h", 0.36)
    kw.setdefault("initials", {"CO": 200.0, "HNO3": 10.0}); kw.setdefault("deposition_lifetime_h", {"HNO3": 2.0})
    return fs.FieldSpec(table=table, **kw)


def _checked(spec):
    table, rep = fs.check(spec)
    assert rep.ok, rep.text()
    return table


def test_grid_and_columns(tmp_path):
    spec = night_spec(night_table(tmp_path / "n.csv"))
    table = _checked(spec)
    grid, cols, notes = fsc.enz_columns(spec, table)
    assert len(grid) == 181 and grid[1] - grid[0] == 60 and grid[0] == table.times_s[0]
    assert set(cols) >= {"sz", "hv", "t", "P", "H2O", "KDIL", "HPBL", "QE_CO", "B_CO", "C_NO", "jNO2"}
    assert np.all(cols["sz"] > 90) and np.all(cols["hv"] == 1)          # night
    assert cols["KDIL"] == pytest.approx(np.full(181, 1e-4))
    assert cols["QE_CO"] == pytest.approx(np.full(181, 1e7)) and cols["B_CO"] == pytest.approx(np.full(181, 100.0))
    assert cols["C_NO"][0] == 1.0 and cols["C_NO"][-1] == pytest.approx(5.0) and cols["C_NO"][15] == pytest.approx(1 + 4 / 12, abs=1e-4)
    assert not notes


def test_emissions_are_held_between_rows_and_follow_the_layer(tmp_path):
    p = tmp_path / "e.csv"
    p.write_text("time,T,p,RH,lat,lon,h,jNO2,E_NO [molec/cm2/s]\n2021-08-15 12:00,298,1000,50,34,-118,1000,0.008,1e12\n"
                 "2021-08-15 12:15,298,1000,50,34,-118,2000,0.008,1e12\n2021-08-15 12:30,298,1000,50,34,-118,2000,0.008,2e12\n",
                 encoding="utf-8")
    spec = fs.FieldSpec(name="e", run_type="lagrangian", table=str(p), chemistry=INORGANIC, dilution="pbl")
    table = _checked(spec)
    grid, cols, _notes = fsc.enz_columns(spec, table)
    e = cols["QE_NO"]
    assert e[0] == pytest.approx(1e7) and e[14] == pytest.approx(1e7) and e[15] == pytest.approx(5e6) and e[-1] == pytest.approx(1e7)
    assert cols["HPBL"][7] == pytest.approx(1000 + 7 * 1000 / 15)       # conditions interpolate linearly
    growth = np.log(2.0) / 900.0
    assert cols["KDIL"][0] == pytest.approx(growth) and cols["KDIL"][20] == pytest.approx(0.3 / 3600.0)


def test_spin_up_repeats_the_table(tmp_path):
    p = tmp_path / "c.csv"
    lines = ["time,T,p,RH,lat,lon,jNO2,C_NO"] + [f"2021-08-15 {h:02d}:00,298,1000,50,34,-118,0.001,{h}" for h in range(24)]
    p.write_text("\n".join(lines) + "\n", encoding="utf-8")
    spec = fs.FieldSpec(name="c", run_type="constrained", table=str(p), chemistry=INORGANIC, spin_up_days=2, step_s=900)
    table = _checked(spec)
    grid, cols, notes = fsc.enz_columns(spec, table)
    assert fsc.spin_up_period(table) == 86400.0
    assert grid[0] == table.times_s[0] - 2 * 86400 and grid[-1] == table.times_s[-1]
    day = 96
    assert cols["C_NO"][0] == 0 and cols["C_NO"][day] == 0 and cols["C_NO"][2 * day] == 0
    assert cols["C_NO"][4] == pytest.approx(1.0) and cols["C_NO"][day + 4] == pytest.approx(1.0)
    assert cols["C_NO"][day - 1] == pytest.approx(23.0)                 # the gap before the next copy holds the last row
    assert any("spin-up" in n for n in notes)


def test_gaps_are_bridged_and_reported(tmp_path):
    p = tmp_path / "g.csv"
    p.write_text("time,T,p,RH,lat,lon,jNO2,C_NO\n2021-08-15 12:00,298,1000,50,34,-118,0.001,1\n2021-08-15 12:15,,1000,50,34,-118,0.001,\n"
                 "2021-08-15 12:30,300,1000,50,34,-118,0.001,3\n", encoding="utf-8")
    spec = fs.FieldSpec(name="g", run_type="constrained", table=str(p), chemistry=INORGANIC)
    table = _checked(spec)
    _grid, cols, notes = fsc.enz_columns(spec, table)
    assert cols["t"][15] == pytest.approx(299.0) and cols["C_NO"][15] == pytest.approx(2.0)
    assert any(n.startswith("T: 1 empty cell") and "30 min" in n for n in notes) and any(n.startswith("C_NO: 1 empty") for n in notes)


def test_written_files_and_parsed_scenario(tmp_path):
    spec = night_spec(night_table(tmp_path / "n.csv"), deposition_velocity_cm_s={"O3": 0.4})
    table = _checked(spec)
    path, res = fsc.prepare(spec, table, str(tmp_path / "out"), INCLUDE_DIR, "", solver_from_toml("", "fast"), export_kv_pd=False)
    easy = open(res["easy"], encoding="latin-1").read()
    for line in ("KDIL = INPUT(KDIL)", "HPBL = INPUT(HPBL)", "DILUTE = CONST(KDIL+0.000000e+00)", "QE_CO = INPUT(QE_CO)", "QE_CO --> CO",
                 "k[QE_CO --> CO] = CONST(1.)", "B_CO --> CO", "k[B_CO --> CO] = CONST(KDIL*1e-9*M)", "C_NO = INPUT(C_NO)",
                 "NO = CONST(C_NO*1e-9*M)", "k[HNO3 -->] = CONST(1./7200)", "k[O3 -->] = CONST(0.4/(HPBL*100.))", "P = INPUT(P)",
                 "INITIAL[CO] = CONST(M*1e-9*200)", "FILES[ENZ]   = night.enz", "jno2 = INPUT(jNO2)",
                 "k[NO2 + HV --> O3 + NO] = CONST(jno2)"):
        assert line in easy, line
    assert "PHTFKT_CS(1.165D-02" not in easy                      # measured jNO2 scales the clear-sky frequencies
    for absent in ("SHONO", "WALLDEP", "qY", "XNOy"):
        assert absent not in easy, absent
    assert os.path.isfile(res["spec_path"]) and os.path.isfile(os.path.join(tmp_path, "out", "night.enz"))
    cfg = json.load(open(path, encoding="utf-8"))
    assert cfg["held_species"] == {"NO": "C_NO*1e-9*M"}
    assert cfg["env_map"]["PRESS"] == {"type": "SAPHIR File", "value": "P"} and cfg["env_map"]["TEMP"]["type"] == "SAPHIR File"
    sources = {(r["reactants"], r["products"]): r["rate"] for r in cfg["injection_reactions"]}
    assert sources[("QE_CO", "CO")] == "(1.) * QE_CO" and sources[("B_CO", "CO")] == "(KDIL*1e-9*M) * B_CO"
    assert cfg["idl_parameters"]["DILUTE"] == "KDIL+0.000000e+00"
    custom = {r["reactants"]: r["rate"] for r in cfg["custom_reactions"] if r["products"] == ""}
    assert custom["HNO3"] == "1./7200" and "HPBL*100." in custom["O3"]
    assert not [n for n in res["notes"] if "empty" in n]


def test_held_species_are_not_misread_in_a_chamber_scenario(tmp_path):
    """A SAPHIR scenario defines settings from INPUTs (SHONO, H2O, T) but holds no species."""
    import shutil
    from scenario_prep import prepare_scenario
    from conftest import IDL_INPUT
    src = os.path.join(IDL_INPUT, "scene_05")
    for ext in (".easy", ".enz"):
        shutil.copy(os.path.join(src, "scene_05" + ext), tmp_path / ("scene_05" + ext))
    path = prepare_scenario({"easy": str(tmp_path / "scene_05.easy"), "enz": str(tmp_path / "scene_05.enz")}, INCLUDE_DIR, False, "")
    cfg = json.load(open(path, encoding="utf-8"))
    assert cfg["held_species"] == {}
    assert cfg["env_map"]["PRESS"] == {"type": "Fixed Constant", "value": "1006.5"}


@pytest.mark.slow
@pytest.mark.parametrize("profile", ["fast", "precise"])
def test_night_matches_the_analytic_solutions(tmp_path, profile):
    """Dark inorganic chemistry: CO = B + E/k + (C0 - B - E/k) e^-kt, HNO3 decays with 1/tau + k, NO follows the
    held ramp with FACSIMILE's row convention (the value of a row is the one read at its start)."""
    import run_benchmark as rb
    spec = night_spec(night_table(tmp_path / "n.csv"))
    table = _checked(spec)
    path, _res = fsc.prepare(spec, table, str(tmp_path / profile), INCLUDE_DIR, "", solver_from_toml("", profile), export_kv_pd=False)
    assert fsc.run(path, "", str(tmp_path / profile / "run.log")) == 0, open(tmp_path / profile / "run.log", encoding="utf-8", errors="replace").read()[-3000:]
    df = rb.parse_cs_file(str(tmp_path / profile / "night.cs"))
    assert len(df) == 181
    t = df["TIME"].to_numpy() - df["TIME"].to_numpy()[0]
    M = 1000.0 / (298.0 * fs.DENSITY_CONST)
    E, k, B, C0 = 1e7, 0.36 / 3600.0, 100e-9 * M, 200e-9 * M
    co = B + E / k + (C0 - B - E / k) * np.exp(-k * t)
    assert np.max(np.abs(df["CO"].to_numpy() / co - 1)) < 1e-3
    hno3 = 10e-9 * M * np.exp(-(1 / 7200.0 + k) * t)
    assert np.max(np.abs(df["HNO3"].to_numpy() / hno3 - 1)) < 1e-2
    held = np.interp(df["TIME"].to_numpy() - 60.0, table.times_s, table.constraints["NO"]) * 1e-9 * M
    assert np.max(np.abs(df["NO"].to_numpy()[1:] / held[1:] - 1)) < 1e-5
    assert df["KDIL"].iloc[5] == pytest.approx(1e-4) and df["P"].iloc[5] == pytest.approx(1000.0)


def test_without_jno2_the_clear_sky_frequencies_are_used(tmp_path):
    p = tmp_path / "d.csv"
    p.write_text("time,T,p,RH,lat,lon,C_NO\n2021-08-15 20:00,298,1000,50,34,-118,1\n2021-08-15 21:00,298,1000,50,34,-118,1\n",
                 encoding="utf-8")
    spec = fs.FieldSpec(name="d", run_type="constrained", table=str(p), chemistry=INORGANIC)
    table = _checked(spec)
    res = fsc.write_case(spec, table, str(tmp_path / "out"))
    easy = open(res["easy"], encoding="latin-1").read()
    assert "jno2 = PHTFKT_CS(1.165D-02,0.244,-0.267)" in easy and "INPUT(jNO2)" not in easy
    assert "jNO2" not in res["cols"] and 15 < res["cols"]["sz"][30] < 45          # 20:30 UTC: early afternoon in Los Angeles
