"""f0am_export: the written F0AM files describe exactly the recorded mechanism (no JUICE run needed)."""
import os

import numpy as np

from juice.importers import f0am_export as fx


def _toy(n=5, dt=60.0):
    # A -> B (constant), B + RO2 -> C (RO2 reaction, time-varying), source -> A (zero order)
    k = np.zeros((n, 3))
    k[:, 0] = 1e-3
    k[:, 1] = np.linspace(1e-12, 2e-12, n)
    k[:, 2] = [0, 1e5, 1e5, 0, 0]
    return fx.Recording(
        stem="toy", species=["A", "B", "C", "RA"], ro2=["RA"],
        rxn_names=["A=B", "B+RO2=C", "Q=A"], rxn_reactants=[["A"], ["B", "RO2"], []],
        rxn_stoich=[{"A": -1.0, "B": 1.0}, {"B": -1.0, "C": 1.0}, {"A": 1.0}],
        k=k, t0=np.arange(n) * dt + 1e8, dt=dt, P=np.full(n, 1013.0), T=np.full(n, 298.0),
        H2O=np.full(n, 3e17), M=np.full(n, 2.46e19), y0={"A": 1e9, "B": 0.0, "C": 0.0, "RA": 1e8},
        k_loss=None, juice_t=np.arange(n + 1) * dt + 1e8,
        juice_y={s: np.zeros(n + 1) for s in ("A", "B", "C", "RA")})


def test_written_mechanism_reads_back(tmp_path):
    rec = _toy()
    d = fx.write(rec, str(tmp_path))
    m = fx.read_mech(os.path.join(d, "toy_JUICE_mech.m"))
    assert m["species"] == rec.species and m["ro2"] == ["RA"]
    assert len(m["rx"]) == 3
    assert m["rx"][0][0] == repr(1e-3) and m["rx"][0][1] == ["A"]
    assert m["rx"][1][0].startswith("JUICEk_(:,") and m["rx"][1][1] == ["B", "RO2"]
    assert m["rx"][2][1] == [] and m["rx"][2][2] == {"A": 1.0}
    for f in ("run_toy.m", "toy_JUICE_J.m", "toy_JUICE_k.mat", "toy_JUICE_inputs.mat", "toy_JUICE.csv"):
        assert os.path.isfile(os.path.join(d, f))
    run = open(os.path.join(d, "run_toy.m")).read()
    assert "'toy_JUICE_J(Met)'" in run          # F0AM needs a J function in ChemFiles{2}
    assert "function J = toy_JUICE_J(Met)" in open(os.path.join(d, "toy_JUICE_J.m")).read()


def test_emulation_follows_the_rate_constants(tmp_path):
    rec = _toy()
    d = fx.write(rec, str(tmp_path))
    sp, Y = fx.emulate(d, "toy")
    a = Y[:, sp.index("A")]
    # step 1: pure decay; step 2: + 1e5 molec/cm3/s source; the source switches off after step 3
    exp1 = 1e9 * np.exp(-1e-3 * 60)
    assert abs(a[0] / exp1 - 1) < 1e-5
    exp2 = exp1 * np.exp(-0.06) + 1e5 / 1e-3 * (1 - np.exp(-0.06))
    assert abs(a[1] / exp2 - 1) < 1e-5
