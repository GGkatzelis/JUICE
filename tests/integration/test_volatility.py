"""Volatility methods and the gas-particle-wall partitioning of an MCM experiment."""
import math
import os

import pytest

import juice
from juice.chemistry import volatility as vol
from juice.chemistry.species_db import load_species_db
from juice.core import experiment_spec as es

M = vol.Molecule


def test_formula_methods_by_hand():
    pinic = M("PINIC", 9, 14, 4, 0, 186.2)                       # C9H14O4
    x = 9 * 4 / 13
    assert math.log10(vol.cstar_donahue(pinic)) == pytest.approx((25 - 9) * 0.475 - 4 * 2.3 + 2 * x * 0.3)
    # Li et al. (2016) Table 1: bCO = -0.7790 for CHO and -0.07790 for CHON (a factor 10 apart; checked in the paper)
    assert math.log10(vol.cstar_li(pinic)) == pytest.approx((22.66 - 9) * 0.4481 - 4 * 1.656 + 2 * x * 0.7790)
    assert math.log10(vol.cstar_li(M("C10H16", 10, 16, 0, 0, 136.2))) == pytest.approx((23.80 - 10) * 0.4861)   # CH class
    assert math.log10(vol.cstar_li(M("amine", 4, 11, 0, 1, 73.1))) == pytest.approx((24.59 - 4) * 0.4066 - 0.9619)  # CHN
    nit = M("N", 10, 17, 5, 1, 231.2)                            # an organic nitrate: the CHON coefficients
    x = 10 * 5 / 15
    assert math.log10(vol.cstar_li(nit)) == pytest.approx((24.13 - 10) * 0.3667 - 5 * 0.7732 + 2 * x * 0.0779 - 1.114)
    assert [vol.volatility_class(c) for c in (1e-5, 1e-2, 10.0, 1e5)] == ["ELVOC", "LVOC", "SVOC", "VOC"]


def test_simpol_from_the_structure():
    def log10p(smiles):
        from rdkit.Chem import Descriptors
        c = vol.cstar_simpol(M("x", 0, 0, 0, 0, 0.0, smiles))
        return math.log10(c * vol.R_ATM * vol.T_REF / 1e6 / Descriptors.MolWt(vol._mol(smiles)))
    assert log10p("CC(=O)OON(=O)=O") == pytest.approx(-1.33, abs=0.02)            # PAN: volatile
    assert log10p("c1ccc(O)c(O)c1OC") == pytest.approx(-7.04, abs=0.02)           # 3-methoxycatechol
    assert log10p("Oc1ccccc1O") < log10p("COc1ccccc1O")                           # catechol below guaiacol
    assert vol.simpol_groups("CCO") == {0: 1, 1: 2, 7: 1}                         # ethanol: zeroeth group, 2 C, 1 OH
    with pytest.raises(ValueError):
        vol.cstar("simpol", M("bad", 1, 1, 1, 0, 30.0, "not a smiles"))


def test_which_species_partition():
    db = load_species_db()
    names = ["PINONIC", "PINIC", "PAN", "APINENE", "HCHO", "CO", "CO2", "CH4", "APINBO2", "CH3CO3", "OH", "NOTASPECIES"]
    r = vol.partitioning_species(names, db, "simpol")
    got = [s[0] for s in r.species]
    assert got == ["PINONIC", "PINIC", "PAN", "APINENE", "HCHO"]                  # no radicals, no CO / CO2 / CH4, no inorganics
    cst = {s[0]: float(s[5]) for s in r.species}
    assert cst["PINIC"] < cst["PINONIC"] < cst["APINENE"] < cst["HCHO"]
    limited = vol.partitioning_species(names, db, "simpol", cstar_max=1e4)
    assert [s[0] for s in limited.species] == ["PINONIC", "PINIC"] and limited.above_limit == 3
    assert "2 species partition" in limited.summary() and "3 above the C* limit" in limited.summary()
    # the formula methods make PAN far stickier than the structure methods do (F. Rohrer's nitrate finding)
    pan = {m: float({s[0]: s[5] for s in vol.partitioning_species(["PAN"], db, m).species}["PAN"]) for m in ("simpol", "li2016", "donahue2011", "li_evap")}
    assert pan["simpol"] > 1e7 and pan["donahue2011"] < 1e2 and pan["li2016"] < pan["li_evap"]
    with pytest.raises(KeyError):
        vol.partitioning_species(names, db, "no_such_method")


def test_every_method_is_offered_or_says_why():
    assert vol.DEFAULT_METHOD in vol.METHODS and {"simpol", "nannoolal", "evaporation", "li_evap", "li2016", "donahue2011"} <= set(vol.METHODS)
    for m in vol.METHODS.values():
        ok, why = m.available()
        assert ok or ("UManSysProp" in why and m.backend == "umansysprop")
        assert m.reference and m.label
    if not vol.umansysprop_backend()[0]:
        with pytest.raises(RuntimeError, match="UManSysProp"):
            vol.partitioning_species(["PINIC"], load_species_db(), "nannoolal")


def _exp(name="part", method="li2016", particles="mono", on=True):
    e = juice.new_experiment(name, chamber="SAPHIR")
    c = e.conditions
    c.source, c.date, c.start, c.end = "manual", "2024-07-09", "08:00", "09:30"
    c.temperature_mode, c.j_source, c.roof_open = "const", "clear_sky", [("08:10", "09:30")]
    e.precursors = ["APINENE"]
    e.injections = [juice.Injection(name="O3", products="O3", total_ppb="60", start="08:05"),
                    juice.Injection(name="apinene", products="APINENE", total_ppb="10", start="08:15")]
    e.partitioning, e.volatility_method, e.particles = on, method, particles
    return e


def test_off_by_default_and_written_when_asked(tmp_path):
    off = _exp("off", on=False)
    assert es.ExperimentSpec().partitioning is False and es.partitioning_kind(off) == ""
    res = es.write_experiment(off, str(tmp_path / "off"))
    assert "SOA.include" not in open(res["easy"], encoding="latin-1").read() and "partitioning" not in res["mechanism"]
    # switched on for an experiment whose subset exists already: only the include is added
    off.partitioning = True
    res2 = es.write_experiment(off, str(tmp_path / "off"))
    assert res2["mechanism"]["partitioning"]["species"] > 50
    on = _exp("on")
    res = es.write_experiment(on, str(tmp_path / "on"))
    easy = open(res["easy"], encoding="latin-1").read()
    assert "#include on.SOA.include" in easy and "KCOND = CONST(" in easy and "; PARTICLES: 100 nm" in easy
    inc = open(res["mechanism"]["partitioning"]["include"], encoding="latin-1").read()
    assert "Li et al. (2016)" in inc.splitlines()[1] and "k[SOA_PINIC-->PINIC]=CONST(KCOND*" in inc and "SOA_SEED" not in inc
    assert "k[PINIC-->WALL_PINIC]=CONST(KFWALL)" in inc and "APINBO2 -->" not in inc
    assert not [m for lvl, m in juice.check(on) if lvl == "error"]
    bad = _exp("bad", method="nope")
    assert any("unknown volatility method" in m for lvl, m in juice.check(bad) if lvl == "error")
    lumped = _exp("lumped"); lumped.base_mechanism = "RACM2"
    if "RACM2" in __import__("juice.chemistry.mechanism_registry", fromlist=["x"]).load_bases():
        assert es.partitioning_kind(lumped) == "" and any("available only for" in m for _l, m in juice.check(lumped))


def test_alpha_pinene_forms_particles_and_loses_vapour_to_the_wall(tmp_path):
    import numpy as np
    run = juice.run(_exp(), str(tmp_path / "run"), solver="fast", kv_pd=False, script=False)
    df = run.results().df
    soa = [c for c in df.columns if c.startswith("SOA_") and c.split("_")[1] not in ("SVOC", "LVOC", "ELVOC")]
    wall = [c for c in df.columns if c.startswith("WALL_") and c.split("_")[1] not in ("PL", "SVOC", "LVOC", "ELVOC")]
    assert len(soa) > 50 and len(soa) == len(wall)
    assert float(df["COA"].values[0]) == pytest.approx(0.5, abs=0.01) and float(df["COA"].values[-1]) > 1.0   # µg m-3
    assert df[wall].values[-1].sum() > 0 and df[soa].values[-1].sum() > 0
    assert np.isfinite(df.values.astype(float)).all()
    # without particles nothing condenses, the walls still take vapour
    run0 = juice.run(_exp("walls", particles="none"), str(tmp_path / "walls"), solver="fast", kv_pd=False, script=False)
    d0 = run0.results().df
    # (every species starts at the INITIAL[*] floor of 1 molecule cm-3, so 'nothing' is about one molecule per species)
    assert d0[[c for c in d0.columns if c in soa]].values[-1].sum() < 2.0 * len(soa)
    assert float(d0["COA"].values[-1]) == pytest.approx(0.5, abs=1e-3)
    assert d0[[c for c in d0.columns if c in wall]].values[-1].sum() > 0


def test_chemistry_page_carries_the_choice():
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PyQt6.QtWidgets import QApplication
    from juice.gui import experiment_wizard as ew
    from juice.core.scenario_prep import load_paths
    app = QApplication.instance() or QApplication([])                # noqa: F841

    class W:
        include_dir, easy_root = load_paths("")["include_dir"], ""
    s = es.ExperimentSpec(name="x"); s.precursors = ["APINENE"]
    s.partitioning, s.volatility_method, s.partition_cstar_max, s.particles = True, "donahue2011", 1e5, "mono"
    page = ew.ChemistryPage(W()); page.set_from(s)
    p = page.partition
    assert p._kind == "computed" and p.enable.isChecked() and p.method.currentData() == "donahue2011"
    assert [p.method.itemData(i) for i in range(p.method.count())] == list(vol.METHODS)
    for i, m in enumerate(vol.METHODS.values()):                      # a method without its package is listed, not selectable
        assert p.method.model().item(i).isEnabled() == m.available()[0]
    out = es.ExperimentSpec(name="y"); page.apply(out)
    assert (out.partitioning, out.volatility_method, out.partition_cstar_max, out.particles) == (True, "donahue2011", 1e5, "mono")
    assert out.gecko_partitioning is True                             # the GECKO-A tick is kept apart
    p.spec_provider = lambda: (lambda sp: (page.apply(sp), sp)[1])(es.ExperimentSpec(name="z"))
    p.preview()
    assert p.species.rowCount() > 50 and "species partition" in p.method_info.text()
    p.enable.setChecked(False); page.apply(out)
    assert out.partitioning is False and out.gecko_partitioning is True
    page.close()


@pytest.mark.skipif(not vol.umansysprop_backend()[0], reason="UManSysProp / Open Babel not installed")
def test_umansysprop_methods_agree_with_simpol_for_aliphatics():
    pan = M("PAN", 2, 3, 5, 1, 121.05, "CC(=O)OON(=O)=O")
    pinic = M("PINIC", 9, 14, 4, 0, 186.2, "CC1(C)C(CC1C(O)=O)CC(O)=O")
    for m in (pan, pinic):
        s = math.log10(vol.cstar_simpol(m))
        assert abs(math.log10(vol.cstar_nannoolal(m)) - s) < 1.0 and abs(math.log10(vol.cstar_evaporation(m)) - s) < 1.0
    # EVAPORATION is not for aromatic rings (catechol about 4 orders above SIMPOL; F. Rohrer saw the same)
    cat = M("catechol", 6, 6, 2, 0, 110.1, "Oc1ccccc1O")
    assert math.log10(vol.cstar_evaporation(cat)) - math.log10(vol.cstar_simpol(cat)) > 3
