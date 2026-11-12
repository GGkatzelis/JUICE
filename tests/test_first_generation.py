"""First-generation chemistry from every kind of base mechanism, and the partitioning values the panel reads."""
import pytest

from juice.chemistry.first_generation import first_generation


@pytest.mark.parametrize("species, base, addons, name, partners", [
    ("C5H8", "MCM_v3.3.1", (), "C5H8", {"OH", "O3", "NO3"}),
    ("FURAN", "MCM_v3.3.1", ("NOAA_BB",), "FURAN", {"OH"}),           # an add-on compound: its EASY include
    ("APINENE", "RACM2", (), "API", {"OH", "O3", "NO3"}),              # MCM name -> RACM2 name via ALIAS
    ("HCHO", "RACM2", (), "HCHO", {"OH", "NO3", "hv"}),                # photolysis written without an hv reactant
    ("C5H8", "CRACMM2", (), "ISO", {"OH", "O3", "NO3"}),
    ("C2H5OH", "GECKOA_ethanol", (), "C2H5OH", {"OH", "NO3"}),         # GECKO-A's HO is the MCM's OH
])
def test_first_generation_in_every_base(species, base, addons, name, partners):
    got, rxns = first_generation(species, base, addons)
    assert got == name
    assert {r["partner"] for r in rxns} == partners
    assert all(r["easy_rate"] for r in rxns), "every reaction has its rate"


def test_partitioning_steps_are_not_chemistry():
    _n, rxns = first_generation("C2H5OH", "GECKOA_ethanol")
    assert not any("SOA_" in r["easy"] or "WALL_" in r["easy"] for r in rxns)


def test_partitioning_values_come_from_the_include():
    from juice.chemistry import mechanism_registry as mreg
    from juice.gui.partitioning_panel import read_soa_include
    v = read_soa_include(mreg.base("GECKOA_ethanol").path("soa"))
    assert (v["KCOND"], v["KFWALL"], v["KPWALL"], v["CW_SCALE"], v["COA0"]) == (5.5e-3, 1e-3, 5e-5, 1.0, 0.5)
    assert len(v["species"]) == 19
    cstar = dict((n, (c, w)) for n, c, w in v["species"])
    assert cstar["ETHGLY"] == (2.437e5, 9101.0)
