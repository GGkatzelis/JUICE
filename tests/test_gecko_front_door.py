"""GECKO-A through the front door: a generator folder becomes a base mechanism of the registry, an experiment (chamber
or field) chooses it like the MCM, the writer emits the includes in the right order for both inorganic cores, and the
Chemistry step shows it.  Uses the ethanol generator output in tests/fixtures/gecko_ethanol and a temporary registry."""
import os
import shutil

import pytest

import experiment_spec as es
import gecko_import as gi
import mechanism_registry as mreg

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
ETHANOL = os.path.join(HERE, "fixtures", "gecko_ethanol")


@pytest.fixture
def registry(tmp_path, monkeypatch):
    """A registry holding the MCM download (unchanged, for the 'mcm' core) and the imported ethanol mechanism."""
    base_dir = tmp_path / "base"
    base_dir.mkdir()
    shutil.copytree(os.path.join(ROOT, "mechanisms", "base", "MCM_v3.3.1"), base_dir / "MCM_v3.3.1")
    monkeypatch.setattr(mreg, "BASE_DIR", str(base_dir))
    b = gi.import_generator_folder(ETHANOL, "GECKOA_ethanol", "C2H5OH", log=lambda *a: None)
    return b


def _spec(base_id: str, core: str, name: str = "eth") -> es.ExperimentSpec:
    s = es.ExperimentSpec(name=name)
    es.apply_chamber(s, "SAPHIR")
    c = s.conditions
    c.source, c.date, c.start, c.end, c.step_s = "manual", "2024-07-09", "08:00", "09:00", 60
    c.temperature_K, c.temperature_mode, c.pressure_hPa, c.h2o_percent, c.flow_m3h = 298.0, "const", 1013.0, 1.0, 0.5
    c.roof_open = [("08:30", "09:00")]
    c.j_source = "clear_sky"
    s.base_mechanism, s.gecko_core, s.gecko_partitioning = base_id, core, True
    s.precursors = ["C2H5OH"]
    s.injections = [es.InjectionEvent(name="ethanol", start="08:40", products="C2H5OH", total_ppb="20")]
    return s


def test_import_registers_a_gecko_base(registry):
    b = registry
    assert b.is_gecko and b.kind == "gecko" and mreg.is_gecko(b.id) and not mreg.is_gecko("MCM_v3.3.1")
    assert b.gecko["precursor"] == "C2H5OH" and b.gecko["primary_code"] == "O02000" and b.gecko["n_species"] == 66
    assert (b.gecko["n_reactions_inorganic"], b.gecko["n_reactions_organic"], b.gecko["n_partitioning"]) == (42, 183, 19)
    for role in ("physics", "anorg", "aliases", "names", "organic", "dil_mech", "dil_rate", "soa", "species_tsv"):
        assert os.path.isfile(b.path(role)), role
    assert not b.verify()                                      # checksums match what was written
    assert b.source["generator_files_sha256"]["dictionary.out"] and b.source["options"]["wall"] == "krechmer"
    names = mreg.species_names(b.id)
    assert names["C2H5OH"] == "C2H6O" and "ETHGLY" in names and "HOCH2CH2O2" in names
    # the dilution pair is the IDL layout: reactions in .dil.mech, k[] lines in .dil.rate
    mech = open(b.path("dil_mech"), encoding="latin-1").read(); rate = open(b.path("dil_rate"), encoding="latin-1").read()
    code = lambda text: [l.split(";", 1)[0].strip() for l in text.splitlines() if l.split(";", 1)[0].strip()]
    assert any(l.startswith("C2H5OH + DIL") for l in code(mech)) and not any(l.startswith("k[") for l in code(mech))
    assert all(l.startswith("k[") for l in code(rate)) and any("C2H5OH + DIL" in l and "CONST(DILUTE)" in l for l in code(rate))
    assert any("GECKO-A generator output" in l for l in mreg.provenance_lines(b.id, []))
    assert "MCM_v3.3.1" in mreg.load_bases() and b.id in mreg.load_bases()


@pytest.mark.parametrize("core", ["gecko", "mcm"])
def test_experiment_writer_emits_the_gecko_includes(registry, tmp_path, core):
    b = registry
    s = _spec(b.id, core)
    assert not [m for lvl, m in es.validate(s) if lvl == "error"], es.validate(s)
    res = es.write_experiment(s, str(tmp_path / core))
    chem = sorted(os.listdir(tmp_path / core / "chemistry"))
    assert b.files["organic"] in chem and b.files["dil_mech"] in chem and b.files["aliases"] in chem
    easy = open(res["easy"], encoding="latin-1").read()
    inc = [l.split(";")[0].split(None, 1)[1].strip() for l in easy.splitlines() if l.strip().startswith("#include")]
    assert inc[:2] == ["geckoa.physics.incl", "geckoa.aliases.incl"]   # M/O2/N2 first, names resolved before any reaction
    assert ("anorg.geckoa.incl" in inc) == (core == "gecko")
    assert ("anorg.mcm331.dark.mech" in inc) == (core == "mcm")
    assert ("saphir.jmeasured_02.incl" in inc) == (core == "mcm")
    assert "j_clear_sky.incl" in inc                            # the anchors need its J<n>
    assert inc.index(b.files["names"]) < inc.index(b.files["organic"]) < inc.index(b.files["soa"])
    assert f"{b.id}.dil.mech" in inc and f"{b.id}.dil.rate" in inc
    assert "jo1d = PHTFKT(6.073D-05,1.743,-0.474)" in easy and "jhono = PHTFKT(2.644D-03" in easy and easy.count("photolysis anchor") == 5
    assert "CONST(J1)" not in easy                                  # J<n> exist only with the MCM phot include
    assert ("k[NO2 + HV --> O3 + NO] = CONST(jno2)" in easy) == (core == "mcm")   # the MCM rate overrides stay with the MCM core only
    assert "CHEMISTRY: base GECKOA_ethanol" in easy
    # partitioning off: no SOA include
    s2 = _spec(b.id, core, "eth2"); s2.gecko_partitioning = False
    res2 = es.write_experiment(s2, str(tmp_path / (core + "2")))
    assert b.files["soa"] not in open(res2["easy"], encoding="latin-1").read()


@pytest.mark.slow
def test_gecko_experiment_runs_in_juice(registry, tmp_path):
    """The written experiment goes through EASY -> JSON and the solver; GECKO-A names are merged into the MCM ones."""
    import subprocess, sys
    import run_benchmark as rb
    import scenario_prep as sp
    b = registry
    s = _spec(b.id, "gecko", "run_gecko")
    res = es.write_experiment(s, str(tmp_path / "run"))
    js = sp.prepare_scenario(res, sp.load_paths("")["include_dir"], export_kv_pd=False, toml_path="", solver=sp.solver_from_toml("", "twostep"))
    r = subprocess.run([sys.executable, os.path.join(ROOT, "scenario_runner.py"), js], capture_output=True, text=True, cwd=ROOT,
                       env={**os.environ, "PYTHONIOENCODING": "utf-8"})
    assert r.returncode == 0, r.stdout[-2000:]
    df = rb.parse_cs_file(str(tmp_path / "run" / "run_gecko.cs"))
    names = {c.split("<")[0] for c in df.columns}
    assert {"C2H5OH", "CH3CHO", "OH", "HO2", "O3", "HONO", "O", "SA", "ETHGLY"} <= names
    assert not ({"HO", "HNO2", "O3P", "SULF", "GOO2000"} & names)
    eth = df[next(c for c in df.columns if c.startswith("C2H5OH"))]
    assert eth.max() > 1e11 and eth.iloc[-1] < eth.max()        # injected, then consumed


@pytest.mark.slow
def test_gecko_mcm_core_follows_the_mcm_subset_control(registry, tmp_path):
    """The same synthetic SAPHIR day with the MCM v3.3.1 subset for ethanol and with the GECKO-A ethanol mechanism on
    the MCM core: inorganic chemistry, chamber sources and photolysis are then identical, only the organic chemistry
    differs, so O3, OH and HONO must agree closely.  This is the check that caught the missing physics lines (every
    termolecular rate 0, OH 1e9) and the CONST(J<n>) anchors (an undefined J1 evaluated as 1 s-1)."""
    import subprocess, sys
    import run_benchmark as rb
    import scenario_prep as sp
    b = registry
    out = {}
    for tag in ("control", "gecko"):
        s = _spec(b.id if tag == "gecko" else "MCM_v3.3.1", "mcm", tag)
        if tag == "control":
            s.gecko_partitioning = True
        res = es.write_experiment(s, str(tmp_path / tag))
        js = sp.prepare_scenario(res, sp.load_paths("")["include_dir"], export_kv_pd=False, toml_path="", solver=sp.solver_from_toml("", "twostep"))
        r = subprocess.run([sys.executable, os.path.join(ROOT, "scenario_runner.py"), js], capture_output=True, text=True, cwd=ROOT,
                           env={**os.environ, "PYTHONIOENCODING": "utf-8"})
        assert r.returncode == 0, r.stdout[-2000:]
        df = rb.parse_cs_file(str(tmp_path / tag / f"{tag}.cs"))
        out[tag] = {n: df[next(c for c in df.columns if c.split("<")[0] == n)].to_numpy() for n in ("O3", "OH", "HONO", "NO2", "C2H5OH")}
    for n, tol in (("O3", 0.05), ("OH", 0.10), ("HONO", 0.10), ("NO2", 0.05), ("C2H5OH", 0.05)):
        a, g = out["control"][n][-1], out["gecko"][n][-1]
        assert abs(g / a - 1) < tol, f"{n}: MCM subset {a:.3e}, GECKO-A on the MCM core {g:.3e}"


def test_field_case_can_use_a_gecko_base(registry, tmp_path):
    import numpy as np
    import field_spec as fs
    import field_scenario as fsc
    from scenario_prep import load_paths
    b = registry
    lines = ["time,lat,lon,T,p,RH,h,jNO2,E_C2H5OH [mol/s],B_O3 [ppb]"]
    for i in range(3):
        lines.append(f"2021-08-15 {12 + i:02d}:00,50.9,6.4,298,1000,40,800,0.008,0.05,40")
    table = tmp_path / "t.csv"; table.write_text("\n".join(lines) + "\n", encoding="utf-8")
    chem = {"base_mechanism": b.id, "addons": [], "precursors": ["C2H5OH"], "mechanism_stem": "", "inorganic_only": False,
            "include_dirs": [], "reaction_adjustments": [], "gecko_core": "mcm", "gecko_partitioning": False}
    spec = fs.FieldSpec(name="parcel", run_type="lagrangian", table=str(table), area_km2=64.0, default_emission_unit="mol/s",
                        dilution="pbl", chemistry=chem)
    assert fs.chemistry_chosen(spec) and "GECKOA_ethanol" in fs.chemistry_label(spec) and "MCM inorganic core" in fs.chemistry_label(spec)
    t, report = fs.check(spec)
    assert report.ok, report.text()
    res = fsc.write_case(spec, t, str(tmp_path / "out"))
    easy = open(res["easy"], encoding="latin-1").read()
    assert "#include geckoa.physics.incl" in easy and "#include geckoa.aliases.incl" in easy and b.files["organic"] in easy and b.files["soa"] not in easy
    assert "QE_C2H5OH --> C2H5OH" in easy and "jo1d = PHTFKT(6.073D-05" in easy
    assert os.path.isfile(os.path.join(tmp_path / "out", "chemistry", b.files["organic"]))


def test_chemistry_page_shows_a_gecko_base(registry):
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PyQt6.QtWidgets import QApplication
    import experiment_wizard as ew
    from scenario_prep import load_paths
    app = QApplication.instance() or QApplication([])
    b = registry

    class W:
        include_dir, easy_root = load_paths("")["include_dir"], ""
    page = ew.ChemistryPage(W())
    labels = [page.base_combo.itemText(i) for i in range(page.base_combo.count())]
    assert any(l.startswith("GECKO-A:") for l in labels) and any("Master Chemical" in l for l in labels)
    page.base_combo.setCurrentIndex(page.base_combo.findData(b.id))
    assert page.is_gecko() and page.gecko_box.isVisibleTo(page) and not page.mcm_box.isVisibleTo(page)
    assert "C2H5OH" in page.gecko_info.text() and "183 organic" in page.gecko_info.text()
    page.gecko_core.setCurrentIndex(page.gecko_core.findData("mcm")); page.gecko_part.setChecked(False)
    s = es.ExperimentSpec(name="x"); page.apply(s)
    assert (s.base_mechanism, s.precursors, s.gecko_core, s.gecko_partitioning, s.mechanism_stem) == (b.id, ["C2H5OH"], "mcm", False, "")
    assert page.validatePage()
    # back again: set_from restores the two options, the MCM page is hidden
    page2 = ew.ChemistryPage(W()); page2.set_from(s)
    assert page2.is_gecko() and page2.gecko_core.currentData() == "mcm" and not page2.gecko_part.isChecked()
    # the reaction search works on the GECKO-A include
    page2.rxn_search.setText("C2H5OH"); page2._find_reactions()
    assert page2.found.count() >= 1 and "C2H5OH" in page2.found.item(0).text()
    page2.base_combo.setCurrentIndex(page2.base_combo.findData("MCM_v3.3.1"))
    assert not page2.is_gecko() and page2.mcm_box.isVisibleTo(page2)
    page.close(); page2.close()
