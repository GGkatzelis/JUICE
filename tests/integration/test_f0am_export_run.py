"""f0am_export on a real JUICE run: the exact-answer experiment (tests/exact_small.py) recorded, written as F0AM files
and integrated again from those files.

The experiment has what the export must fold into F0AM's rate table: a light switch (HV), an injection input (QJ,
a zero-order reaction in F0AM), explicit dilution (DIL x FLOW) and a temperature step (k(T) of N --> O).  The files
F0AM would read must reproduce the closed-form answers, not only JUICE's own numbers."""
import os
import shutil
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
import exact_small  # noqa: E402

from juice.importers import f0am_export as fx  # noqa: E402

STEM = "exact_small"
pytestmark = pytest.mark.slow            # one precise run in-process (~30 s)


@pytest.fixture(scope="module")
def export(tmp_path_factory):
    d = tmp_path_factory.mktemp("f0am_export")
    js = exact_small.prepare(d / "work", "precise")
    rec = fx.record(js, STEM)
    return rec, fx.write(rec, str(d / "out"))


def test_the_rate_table_carries_every_setting(export):
    rec, _d = export
    assert rec.species == list("ABCDEFGHIJLNO") and rec.ro2 == []
    assert rec.rxn_names == ["A=B", "C=D", "D=C", "E+F=G", "H=I", "QJ=J", "L+DIL=", "N=O"]
    assert rec.rxn_reactants[5] == [] and rec.rxn_stoich[5] == {"J": 1.0}            # injection: zero order
    assert rec.rxn_reactants[6] == ["L"]                                              # DIL folded into k
    t = rec.t0 - rec.t0[0]
    k = dict(zip(rec.rxn_names, rec.k.T))
    lit = (t >= 1800) & (t < 5400)
    assert np.all(k["H=I"][lit] == 5e-4) and np.all(k["H=I"][~lit] == 0.0)          # HV
    inj = k["QJ=J"] > 0
    assert list(t[inj]) == [600.0, 660.0]                                            # QJ rows of the ENZ
    assert k["QJ=J"][inj] == pytest.approx(exact_small.M * 1e-9 / 60 * 8.0 * 0.5, rel=1e-6)
    assert np.all(k["L+DIL="][t < 3600] == 0.0) and k["L+DIL="][-1] == pytest.approx(27.6 / (276.0 * 3600.0))
    assert k["N=O"][0] == pytest.approx(1e-2 * np.exp(-4.0)) and k["N=O"][-1] == pytest.approx(1e-2 * np.exp(-2.0))
    assert np.all(rec.k[:, :4] == rec.k[0, :4])                                      # constants stay constant
    assert rec.y0["A"] == pytest.approx(1e9) and rec.y0["J"] == 0.0


def test_the_written_files_give_the_exact_answers(export):
    rec, d = export
    for f in (f"{STEM}_JUICE_mech.m", f"{STEM}_JUICE_k.mat", f"{STEM}_JUICE_inputs.mat", f"run_{STEM}.m",
              f"{STEM}_JUICE_J.m", f"{STEM}_JUICE.csv", f"{STEM}_export.txt"):
        assert os.path.isfile(os.path.join(d, f)), f
    m = fx.read_mech(os.path.join(d, f"{STEM}_JUICE_mech.m"))
    # the ENZ pressure has six decimals, so M changes by ~6e-10 per step: one first-order loss per species carries it
    assert rec.k_loss is not None and np.abs(rec.k_loss).max() < 1e-11
    assert m["species"] == rec.species and len(m["rx"]) == 8 + 13
    assert [r[1] for r in m["rx"][8:]] == [[s] for s in rec.species]
    sp, Y = fx.emulate(d, STEM)
    worst = 0.0
    for i, t in enumerate(rec.t0 + rec.dt - rec.t0[0]):
        for s, v in exact_small.exact(t).items():
            if v > 1e6:
                worst = max(worst, abs(Y[i, sp.index(s)] / v - 1))
    assert worst < 1e-4, worst


def test_check_and_compare(export, tmp_path):
    _rec, d0 = export
    d = str(tmp_path / STEM)
    shutil.copytree(d0, d)
    text = fx.check(d, STEM)
    assert "within 0.1 %: 13/13" in text
    assert "within 0.1 %: 13/13" in fx.check(d, STEM, ro2_frozen=True)
    assert fx.compare_outputs(d, STEM).endswith("(run run_exact_small.m in MATLAB)")
    shutil.copy(os.path.join(d, f"{STEM}_JUICE.csv"), os.path.join(d, f"{STEM}_F0AM.csv"))   # a perfect F0AM
    out = fx.compare_outputs(d, STEM)
    assert "F0AM vs JUICE" in out and "within 0.1 %: 13/13" in out and "FACSIMILE" not in out
    # the command line: --compare reads, --scripts-only rewrites the MATLAB scripts of an existing export
    os.remove(os.path.join(d, f"run_{STEM}.m"))
    assert fx.main([STEM, "--out", str(tmp_path), "--scripts-only"]) == 0
    assert os.path.isfile(os.path.join(d, f"run_{STEM}.m"))
    assert fx.main(["no_such_scene", "--out", str(tmp_path), "--scripts-only"]) == 1
    assert fx.main([STEM, "--out", str(tmp_path), "--compare"]) == 0


def test_deviation_table_uses_the_benchmark_floor():
    a = np.array([[1.0e5, 10.0], [2.0e5, 20.0]])
    b = np.array([[1.0e5, 0.0], [1.0e5, 0.0]])
    dev = fx.deviation_table(["X", "Y"], a, ["Y", "X", "Z"], b[:, ::-1])
    assert dev[0][1] == "X" and dev[0][0] == pytest.approx(1.0)                      # 1e5 off on a 1e5 peak
    assert dev[1] == (pytest.approx(20.0 / fx.COMPARE_FLOOR), "Y", 0.0)             # zero reference: the floor


def test_an_unknown_scene_fails_with_exit_code_1(tmp_path, capsys):
    assert fx.main(["no_such_scene", "--out", str(tmp_path), "--profile", "fast"]) == 1
    assert "no_such_scene: FAILED" in capsys.readouterr().out
