"""The check that every injected species is in the chosen mechanism, for each kind of base mechanism."""
import pytest

import juice
from juice.core.experiment_spec import InjectionEvent


def _flagged(base, injected, precursors=()):
    exp = juice.new_experiment("probe", "SAPHIR")
    exp.base_mechanism = base
    exp.precursors = list(precursors)
    exp.injections = [InjectionEvent(sp, sp, "5", "10:00", 2) for sp in injected]
    return {m.split("'")[3]: m for _lvl, m in juice.check(exp) if "is not a species" in m}


@pytest.mark.parametrize("base, inside, outside", [
    ("GECKOA_ethanol", ["C2H5OH", "HCHO", "NO", "O3"], ["TOLUENE", "C5H8", "APINENE"]),
    ("GECKOA_toluene", ["TOLUENE", "C2H5OH"], ["C5H8", "LIMONENE"]),
    ("RACM2", ["TOLUENE", "C5H8", "APINENE"], ["NOTASPECIES"]),    # MCM names of its lumped species
    ("CRACMM2", ["TOLUENE", "HCHO"], ["CO2", "H2"]),               # not in CRACMM2, though common elsewhere
])
def test_complete_base_checks_its_own_species(base, inside, outside):
    from juice.chemistry import mechanism_registry as mreg
    if base not in mreg.load_bases():                  # GECKO-A imports other than ethanol are git-ignored (size)
        pytest.skip(f"{base} is not imported on this machine")
    flagged = _flagged(base, inside + outside)
    assert set(flagged) == set(outside)
    assert all("used as a whole" in m for m in flagged.values())


def test_precursors_do_not_excuse_a_complete_base():
    # precursors are ignored for a GECKO-A base, so naming one there must not hide a missing species
    assert set(_flagged("GECKOA_ethanol", ["TOLUENE"], precursors=["TOLUENE"])) == {"TOLUENE"}


def test_mcm_compounds_mode_unchanged():
    flagged = _flagged("MCM_v3.3.1", ["C5H8", "TOLUENE", "CO2", "NOTASPECIES"], precursors=["C5H8"])
    assert set(flagged) == {"TOLUENE", "NOTASPECIES"}
    assert not any("used as a whole" in m for m in flagged.values())


def test_chamber_tracer_is_known():
    exp = juice.new_experiment("probe", "SAPHIR")
    exp.base_mechanism = "CRACMM2"
    exp.chamber.tracer_dilution = ["CO2"]
    exp.injections = [InjectionEvent("CO2", "CO2", "5", "10:00", 2)]
    assert not [m for _lvl, m in juice.check(exp) if "is not a species" in m]


@pytest.fixture
def imported_gecko(tmp_path, monkeypatch):
    """A GECKO-A base imported here, as a user imports one: the ethanol generator output of tests/fixtures under an id
    of its own, in a temporary registry.  Stands in for the large imports (GECKOA_toluene, ...), which are git-ignored
    and so missing on every machine but the one that made them."""
    import os
    import shutil
    from juice._root import ROOT
    from juice.chemistry import mechanism_registry as mreg
    from juice.importers import gecko_import as gi
    base_dir = tmp_path / "base"
    base_dir.mkdir()
    shutil.copytree(os.path.join(ROOT, "chemistry", "mechanisms", "base", "MCM_v3.3.1"), base_dir / "MCM_v3.3.1")
    monkeypatch.setattr(mreg, "BASE_DIR", str(base_dir))
    folder = os.path.join(ROOT, "tests", "fixtures", "gecko_ethanol")
    return gi.import_generator_folder(folder, "GECKOA_imported", "C2H5OH", log=lambda *a: None).id


def test_an_imported_gecko_base_checks_its_own_species(imported_gecko):
    flagged = _flagged(imported_gecko, ["C2H5OH", "HCHO", "NO", "O3", "TOLUENE", "C5H8"])
    assert set(flagged) == {"TOLUENE", "C5H8"}
    assert all("used as a whole" in m for m in flagged.values())
