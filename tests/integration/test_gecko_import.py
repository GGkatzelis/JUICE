"""The GECKO-A importer (F. Rohrer's gecko_mechanism / gecko_photolysis / gecko_to_easy) on the ethanol generator
output in tests/fixtures/gecko_ethanol: what the reader finds, what the converter writes, and that the written
chemistry goes through JUICE's own EASY parser."""
import os
import re

import pytest

from juice.importers import gecko_mechanism as gm
from juice.importers import gecko_photolysis as gp
from juice.importers import gecko_to_easy as g
from juice.core.evaluator import fortran_d_to_e

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ETHANOL = os.path.join(HERE, "fixtures", "gecko_ethanol")


@pytest.fixture(scope="module")
def mech():
    return gm.read_mechanism(ETHANOL)


def _write_all(mech, folder, name="Ethanol_test"):
    g.MCM_NAMES.setdefault("O02000", "C2H5OH")          # what build() does with the template's injected VOC
    conv = g.convert(mech)
    n_an = g.write_anorg_include(conv, mech, os.path.join(folder, g.ANORG_INCLUDE))
    g.write_alias_include(os.path.join(folder, g.ALIAS_INCLUDE))
    g.write_names_include(mech, os.path.join(folder, f"{name}.geckoa.names.incl"))
    n_org = g.write_organic_include(conv, mech, os.path.join(folder, f"{name}.geckoa.organic.incl"))
    g.write_dil_include(conv.gas_species, os.path.join(folder, f"{name}.dil.incl"))
    soa = g.write_soa_include(mech, os.path.join(folder, f"{name}.SOA.include"), g.CSTAR_MAX, "krechmer")
    n_tsv = g.write_species_tsv(mech, soa, os.path.join(folder, f"{name}.species.tsv"))
    return conv, n_an, n_org, soa, n_tsv


def test_reader_finds_the_generator_output(mech):
    assert mech.primaries == ["O02000"]                               # CH3CH2(OH)
    assert len(mech.species) == 66 and len(mech.reactions) == 289     # gas <-> particle/wall transfer lines are not reactions
    assert mech.species["O02000"].notation == "CH3CH2(OH)" and mech.species["O02000"].formula == "C2H6O"
    assert mech.species["2O2000"].is_peroxy and mech.species["2O2000"].radical and not mech.species["OO2000"].radical
    assert gm.strip_phase("GOO2000") == "OO2000" and gm.strip_phase("AOO2000") == "OO2000"
    assert {r.kind for r in mech.reactions} >= {"arrh", "tbody", "falloff", "extra", "hv", "pero"}
    assert mech.pero[8] == ["2O2000", "2D2000", "2A2000", "2G2000", "2P2000"]
    assert mech.tref == 298.0 and mech.pvap["OO2000"] == (pytest.approx(9.596e-05), pytest.approx(66.1))
    hv = [r for r in mech.reactions if r.kind == "hv"]
    assert hv and all(r.aux for r in hv)                              # every HV reaction carries its label line


def test_converter_writes_the_includes(mech, tmp_path):
    conv, n_an, n_org, soa, n_tsv = _write_all(mech, str(tmp_path))
    assert (n_an, n_org, len(soa)) == (42, 183, 19) and n_tsv >= 66      # ethanol scheme: 289 generator reactions, duplicates merged
    org = (tmp_path / "Ethanol_test.geckoa.organic.incl").read_text(encoding="latin-1")
    assert "GTROE = FUNCTION[" in org and "GPHOT = FUNCTION[" in org and "C2H5OH" in org
    keys = re.findall(r"^k\[(.+?)\]\s*=\s*(.+)$", org, re.M)
    assert len(keys) == n_org
    # rate expressions contain no Fortran D exponents and no name the D->E conversion would touch
    assert all(fortran_d_to_e(e) == e for _, e in keys)
    tsv = (tmp_path / "Ethanol_test.species.tsv").read_text(encoding="utf-8").splitlines()
    assert tsv[0].startswith("Name\tSmiles")
    row = next(l.split("\t") for l in tsv[1:] if l.startswith("C2H5OH\t"))
    assert row[1] == "CCO" and row[4] == "C2H6O"                      # SMILES from rdkit, formula from the dictionary


def test_photolysis_labels_are_fitted():
    labels = gp.known_labels()
    assert len(labels) == 137 and gp.is_measured(2) and "jo1d" in gp.easy_rate(2)
    other = next(l for l in labels if not gp.is_measured(l))
    assert "GPHOT(" in gp.easy_rate(other)


def test_generated_chemistry_parses_in_juice(mech, tmp_path):
    """A minimal EASY that includes the generated files goes through elog_generator with every reaction kept."""
    from juice.core.elog_generator import generate_both
    _conv, n_an, n_org, _soa, _n = _write_all(mech, str(tmp_path))
    easy = tmp_path / "Ethanol_test.easy"
    easy.write_text("\n".join([
        ";TITLE=GECKO-A ethanol import test", "CONST=FUNCTION[%1]", "DENSITY=FUNCTION[(%1)/((%2)*1.379E-19)]",
        "M =DENSITY(P,T)", "O2=CONST(0.20946*M)", "N2=CONST(0.78084*M)", "T = CONST(298.)", "P = CONST(1000.)",
        "SZ = CONST(40.)", "jno2 = CONST(8e-3)", "jo1d = CONST(2e-5)", "jhono = CONST(1.5e-3)", "jh2o2 = CONST(6e-6)",
        "jhcho_m = CONST(4e-5)", "jhcho_r = CONST(3e-5)", "H2O = CONST(3e17)", "HV = CONST(1.)",
        f"#include {g.ALIAS_INCLUDE}", f"#include {g.ANORG_INCLUDE}", "#include Ethanol_test.geckoa.names.incl",
        "#include Ethanol_test.geckoa.organic.incl", ""]), encoding="latin-1")
    _raw, canon = generate_both(str(easy), search_dirs=[str(tmp_path)], output_dir=str(tmp_path), function_scope="last")
    text = open(canon, encoding="latin-1").read()
    n_rxn = sum(1 for l in text.splitlines() if "-->" in l and not l.strip().startswith(("k[", ";")))
    assert n_rxn == n_an + n_org, (n_rxn, n_an, n_org)
    rxn_lines = [l for l in text.splitlines() if "-->" in l and not l.strip().startswith(";")]
    assert not any("GOO2000" in l for l in rxn_lines) and any("ETHGLY" in l for l in rxn_lines)   # codes -> chosen names
    # the ALIAS lines take effect when the runner's parser reads the elog: GECKO-A's HO, HNO2, O3P, SULF are the
    # MCM's OH, HONO, O, SA - one species each, never two
    from juice.core.easy_parser import EasyParser
    m = EasyParser(base_path=str(tmp_path)).parse(canon)
    names = {getattr(s, "name", s) for s in m.species}
    names = {n.upper() if isinstance(n, str) else str(n).upper() for n in names}
    assert {"OH", "HONO", "O", "SA", "C2H5OH", "ETHGLY"} <= names, sorted(names)[:20]
    assert not ({"HO", "HNO2", "O3P", "SULF", "GOO2000"} & names)
    assert len(m.reactions) == n_an + n_org


# ── the command line: a whole scenario from a SAPHIR template (25 Nov 2026) ─────────────────────────────────────────
# scene_05 is an ethanol experiment whose injection is a delayed release (qVOC --> 0.9*VOCin1 + 0.1*VOCin2, the two
# intermediates --> VOCin_all --> C2H5OH).  Until 25 Nov 2026 build() took the intermediates for MCM-only species,
# disabled the whole chain and wrote a scenario into which no ethanol was ever injected.
SCENE_05 = os.path.join(os.path.dirname(HERE), "validation", "IDL Input", "scene_05", "scene_05.easy")


def test_template_voc_follows_a_delayed_release():
    tl = open(SCENE_05, encoding="latin-1").read().splitlines()
    assert g.template_voc(tl) == "C2H5OH"
    assert g.template_voc(["qVOC --> C2H5OH", "k[qVOC --> C2H5OH] = CONST(1.)"]) == "C2H5OH"
    assert g.template_voc(["qVOC --> 0.5*A + 0.5*B"]) is None                     # two ends: no single VOC
    assert g.template_voc(["NO --> NO2"]) is None
    assert g._chain_intermediates(tl, {"QVOC", "C2H5OH"}) == {"VOCIN1", "VOCIN2", "VOCIN_ALL"}
    assert g._chain_intermediates(tl, {"QVOC"}) == set()                         # the chain ends on an unknown species


def _active(text):
    return [l.strip() for l in text.splitlines() if l.strip() and not l.lstrip().startswith((";", "*"))]


@pytest.mark.parametrize("mode", ["mcmcore", "pure"])
def test_build_writes_a_scenario_juice_can_expand(mode, tmp_path, monkeypatch):
    from juice.core.elog_generator import generate_both
    monkeypatch.setattr(g, "MCM_NAMES", dict(g.MCM_NAMES))                       # build() names the primary in it
    args = [ETHANOL, "--template", SCENE_05, "--name", "EtG", "--mode", mode, "--out", str(tmp_path / "easy"),
            "--juice-dir", str(tmp_path / "juice")] + (["--maxgen", "2"] if mode == "pure" else [])
    g.main(args)
    out = tmp_path / "easy" / "EtG"
    assert sorted(os.listdir(out)) == sorted(os.listdir(tmp_path / "juice" / "EtG"))
    assert (out / "scene_05.enz").is_file()
    text = (out / "EtG.easy").read_text(encoding="latin-1")
    act = _active(text)
    assert "qVOC-->0.9*VOCin1+0.1*VOCin2" in act and "VOCin_all-->C2H5OH" in act   # the injection is kept
    assert not any(".fac." in l for l in act if l.startswith("#include"))            # the MCM converter files are gone
    assert ("#include anorg.geckoa.incl" in " ".join(act)) == (mode == "pure")
    assert any(l.startswith("#include anorg.mcm331.dark.rate") for l in act) == (mode == "mcmcore")
    if mode == "pure":
        assert "maxgen=2 (static generation cut" in text
        assert not any(l.startswith("k[NO2 + HV") for l in act)                    # measured-J overrides of the MCM set
    raw, _canon = generate_both(str(out / "EtG.easy"), search_dirs=[str(out), os.path.join(os.path.dirname(SCENE_05), "..", "include files")],
                                output_dir=str(tmp_path / "x"))
    elog = open(raw, encoding="latin-1").read()
    assert "VOCin_all-->C2H5OH" in elog.replace(" ", "") and "ETHGLY" in elog


def test_build_refuses_an_unknown_mode_and_a_missing_enz(tmp_path):
    with pytest.raises(ValueError):
        g.build(ETHANOL, SCENE_05, "x", "both", str(tmp_path))
    t = tmp_path / "t.easy"
    t.write_text("FILES[ENZ] = nowhere.enz\nqVOC --> C2H5OH\n", encoding="latin-1")
    with pytest.raises(FileNotFoundError):
        g.build(ETHANOL, str(t), "x", "mcmcore", str(tmp_path), log=lambda *_: None)
