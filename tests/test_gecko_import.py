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

HERE = os.path.dirname(os.path.abspath(__file__))
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
