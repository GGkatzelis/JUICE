"""Particles of the gas-particle-wall partitioning: the condensation sink, the size-distribution file, the EASY lines
of a GECKO-A run, and the wall terms scaled to the chamber."""
import os

import numpy as np
import pytest

import juice
from juice.core import particles as pt
from juice.core.experiment_spec import enz_columns, partitioning_lines, render_easy, validate


def test_condensation_sink_values():
    assert pt.condensation_sink(160, 1800) == pytest.approx(4.33e-3, rel=0.01)      # the SAPHIR seed of the include
    assert pt.condensation_sink(100, 1e4) == pytest.approx(1.10e-2, rel=0.01)
    assert pt.condensation_sink(100, 0) == 0.0
    # continuum limit (large particles, Kn -> 0): beta -> 1
    assert pt.fuchs_sutugin(np.array([1e-3]))[0] == pytest.approx(1.0, abs=1e-3)


def test_size_distribution_file(tmp_path):
    f = tmp_path / "smps.csv"
    f.write_text("# made up\ntime;Dp_50 nm;Dp_100 nm;Dp_200 nm;flag\n09:00;2000;10000;2000;ok\n12:00;20000;100000;20000;ok\n")
    sd = pt.read_size_distribution(str(f), "2024-07-09")
    assert list(sd.d_nm) == [50, 100, 200]
    assert any("ignored" in n for n in sd.notes)                                     # the flag column
    w = pt._dlogdp(sd.d_nm)
    assert sd.n_cm3[0] == pytest.approx(np.array([2000, 10000, 2000]) * w)         # dN/dlogDp -> per bin
    cs = sd.cs()
    assert cs[1] == pytest.approx(10 * cs[0])
    grid = sd.times_s[0] + np.array([-600.0, 0.0, 5400.0, 10800.0, 12000.0])
    on = pt.cs_on_grid(str(f), "2024-07-09", grid)
    assert on[0] == pytest.approx(cs[0]) and on[-1] == pytest.approx(cs[1])        # held beyond the table
    assert on[2] == pytest.approx((cs[0] + cs[1]) / 2)


def _gecko(**kw):
    exp = juice.new_experiment("prt", kw.pop("chamber", "SAPHIR"))
    exp.base_mechanism = "GECKOA_ethanol"
    c = exp.conditions
    c.date, c.start, c.end = "2024-07-09", "09:00", "10:00"
    for k, v in kw.items():
        setattr(exp, k, v)
    return exp


def test_default_is_no_particles():
    exp = _gecko()
    assert exp.particles == "none"
    assert "KCOND = CONST(0.0000E+00)" in "\n".join(partitioning_lines(exp))
    assert "KCOND = CONST(0.0000E+00)" in render_easy(exp)


def test_one_size_and_file_lines(tmp_path):
    exp = _gecko(particles="mono", particle_d_nm=100.0, particle_n_cm3=1e4)
    cs = pt.condensation_sink(100, 1e4, exp.conditions.temperature_K)
    assert f"KCOND = CONST({cs:.4E})" in "\n".join(partitioning_lines(exp))
    f = tmp_path / "smps.csv"
    f.write_text("time,100\n09:00,30000\n10:00,30000\n")
    exp = _gecko(particles="file", particle_file=str(f))
    text = "\n".join(partitioning_lines(exp))
    assert "CS_in = input(CS)" in text and "KCOND = CONST(CS_in)" in text
    _grid, cols, _notes = enz_columns(exp)
    assert np.allclose(cols["CS"], pt.condensation_sink(100, 30000 * pt._dlogdp(np.array([100.0]))[0], exp.conditions.temperature_K))


def test_missing_file_is_an_error():
    exp = _gecko(particles="file", particle_file="no_such_file.csv")
    assert any(lvl == "error" and "size-distribution file" in m for lvl, m in validate(exp))


def test_walls_scale_with_the_chamber():
    assert not any(l.startswith(("KFWALL", "CW_SCALE")) for l in partitioning_lines(_gecko()))      # SAPHIR: as written
    lines = partitioning_lines(_gecko(chamber="CHARME"))
    s = (22.0 / 9.2) / (320.0 / 276.0)
    assert f"KFWALL = CONST({1e-3 * s:.4E})" in lines and f"CW_SCALE = CONST({s:.4g})" in lines
    exp = _gecko(chamber="AIDA")                                                                      # wall area unknown
    assert not any(l.startswith("KFWALL") for l in partitioning_lines(exp))
    assert any(lvl == "info" and "wall area" in m for lvl, m in validate(exp))


def test_mcm_runs_unchanged():
    exp = juice.new_experiment("prt", "SAPHIR")
    exp.precursors = ["C5H8"]
    exp.particles = "mono"
    assert "KCOND" not in render_easy(exp)


def test_particles_survive_a_script(tmp_path):
    exp = _gecko(particles="mono", particle_d_nm=150.0, particle_n_cm3=5e3)
    path = juice.write_script(exp, str(tmp_path / "prt_run.py"))
    back = juice.build_from_script(path)
    assert (back.particles, back.particle_d_nm, back.particle_n_cm3) == ("mono", 150.0, 5e3)
