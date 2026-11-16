"""Photolysis reactions that the light switch does not close are reported by the experiment check."""
import os
import re

import pytest

import juice
from juice.chemistry import mechanism_registry as mreg
from juice.chemistry import photolysis_check as pc
from juice.core import experiment_spec as es
from juice.importers import f0am_import as fi


def _exp(base="", addons=(), j="measured", chamber="SAPHIR"):
    exp = juice.new_experiment("pc_test", chamber=chamber)
    c = exp.conditions
    c.source, c.date, c.start, c.end = "manual", "2024-07-09", "08:00", "10:00"
    c.temperature_mode, c.j_source, c.roof_open = "const", j, [("09:00", "10:00")]
    if base:
        exp.base_mechanism, exp.precursors = base, []
    else:
        exp.precursors = ["TOLUENE"]
    exp.addons = list(addons)
    exp.injections = [juice.Injection(name="NO", products="NO", total_ppb="5", start="08:05")]
    return exp


def test_find_on_plain_text():
    easy = "\n".join([
        "#include saphir.jmeasured_02.incl", "jno2 = INPUT(jno2)", "HV_in = INPUT(hv)", "HV = CONST(HV_in)",
        "JA = PHTFKT(1e-5,1,-0.3)", "JB = CONST(jno2*0.1*HV)", "SRC = CONST(HV_in*2.0)",
        "k[A --> B] = CONST(JA)",                    # PHTFKT of the measured-J include: no switch
        "k[A + hv --> C] = CONST(JA)",               # hv as a reactant
        "k[D --> E] = CONST(JB*0.5)",                # the switch inside the J
        "k[F --> G] = CONST(jno2*0.3)",              # the measured column itself
        "k[H --> I] = CONST(jno2*0.3*HV_in)",
        "k[SRC --> HONO] = CONST(4.7e13*jno2)",      # a source whose reactant constant carries the switch
        "k[L + OH --> M] = CONST(1.2D-11)"])         # no photolysis
    hits = pc.find(easy, [])
    assert [h[0] for h in hits] == ["A --> B", "F --> G"]
    assert [h[0] for h in pc.find(easy.replace("#include saphir.jmeasured_02.incl", "#include j_clear_sky.incl"), [])] == ["F --> G"]
    assert [h[0] for h in pc.find(easy, [], inputs_follow_switch=True)] == ["A --> B"]      # lamps: jNO2 is the lamps' own


def test_shipped_chemistry_is_clean():
    assert not pc.ungated(es.spec_scene_05()) and not pc.ungated(es.spec_scene_09())
    for base in ("RACM2", "CRACMM2", "CRACMM1_aq"):
        if base not in mreg.load_bases():
            continue
        for j, chamber in (("measured", "SAPHIR"), ("clear_sky", "SAPHIR"), ("clear_sky", "generic_lamps")):
            assert not pc.ungated(_exp(base, j=j, chamber=chamber)), (base, j, chamber)
    if "NOAA_BB" in mreg.load_addons():
        for j in ("measured", "clear_sky"):
            assert not pc.ungated(_exp(addons=["NOAA_BB"], j=j)), j


def test_every_installed_addon_and_lumped_base_is_clean():
    """Plan C for chemistry that arrives later: whatever is in chemistry/mechanisms is checked, not a fixed list."""
    assert pc.check_installed() == []


def test_a_new_addon_without_the_switch_is_caught(tmp_path, monkeypatch):
    import json
    import shutil
    addons = tmp_path / "addons"
    shutil.copytree(mreg.ADDON_DIR, addons)
    new = addons / "NEW_TEST"
    new.mkdir()
    (new / "new_test.incl").write_text("JNEW = CONST(J41*0.5)\nk[NEWX --> NEWY] = CONST(JNEW)\n"
                                       "k[NEWZ + hv --> NEWY] = CONST(J41)\n", encoding="utf-8")
    (new / "manifest.json").write_text(json.dumps({"kind": "addon", "id": "NEW_TEST", "name": "test", "version": "0",
                                                   "files": ["new_test.incl"], "adds_species": {"NEWX": "C2H4O"}}),
                                       encoding="utf-8")
    monkeypatch.setattr(mreg, "ADDON_DIR", str(addons))
    hits = [h for h in pc.check_installed(lights=(("measured", "SAPHIR"),)) if "NEW_TEST" in h]
    assert len(hits) == 1 and "NEWX --> NEWY" in hits[0] and "new_test.incl" in hits[0]


def test_hand_written_reaction_is_reported():
    exp = _exp()
    assert not [m for _l, m in juice.check(exp) if "light switch" in m]
    exp.overrides = ["XTEST --> YTEST", "k[XTEST --> YTEST]=CONST(J41)", "k[ZTEST --> YTEST]=CONST(J41*HV_in)"]
    msgs = [m for lvl, m in juice.check(exp) if lvl == "warning" and "light switch" in m]
    assert len(msgs) == 1 and "XTEST --> YTEST" in msgs[0] and "ZTEST" not in msgs[0] and "'+ hv'" in msgs[0]


def test_mechanism_j_names_do_not_take_a_measured_input_name():
    """JUICE keeps one variable per name whatever its case: RACM2's JNO2 was the measured input jno2, so its NO2
    photolysis read the building-roof measurement past the light switch."""
    for base in ("RACM2", "CRACMM2", "CRACMM1_aq"):
        if base not in mreg.load_bases():
            continue
        names = re.findall(r"^(\w+) = CONST\(", open(mreg.base(base).path("phot")).read(), re.M)
        assert names and not [n for n in names if n.upper() in fi.RESERVED_J], base
    if "RACM2" in mreg.load_bases():
        b = mreg.base("RACM2")
        assert "k[NO2 --> O3P + NO] = CONST(JNO2_MECH)" in open(b.path("chemistry")).read()
        assert "JNO2_MECH = CONST(J4)" in open(b.path("phot")).read() and b.status == "reviewed"


@pytest.mark.skipif("RACM2" not in mreg.load_bases(), reason="RACM2 not installed")
def test_racm2_photolysis_is_off_with_the_roof_closed(tmp_path):
    import numpy as np
    csv = tmp_path / "meas.csv"
    csv.write_text("time,jNO2,jO1D\n07:30,5e-3,1e-5\n10:30,5e-3,1e-5\n")       # the roof instrument never reads zero
    exp = _exp("RACM2")
    exp.conditions.source, exp.conditions.input_file = "file", str(csv)
    run = juice.run(exp, str(tmp_path / "run"), solver="fast", kv_pd=False, script=False)
    df = run.results().df
    frac = np.linspace(0.0, 2.0, len(df))                                         # hours since 08:00; roof opens at 09:00
    for col in ("JNO2_MECH", "JO1D_MECH"):
        v = df[col].values.astype(float)
        assert v[frac < 0.95].max() == 0.0 and v[frac > 1.05].max() > 0.0, col
