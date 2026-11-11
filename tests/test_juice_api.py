"""juice.py, JUICE from code: the experiment scripts rebuild exactly what the windows hold.

The promise of <name>_run.py is that build() returns the same experiment, field by field, and that
the files it writes are byte-identical to the window's.  Checked on the golden SAPHIR scenes, the
lamp-chamber check, every template, every field example, and through the two windows themselves."""
import filecmp
import os
from dataclasses import asdict

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import pytest

import experiment_spec as es
import field_spec as fs
import juice
from check_lamp_chamber import lamp_spec


def _chamber_specs():
    out = {"scene_05": es.spec_scene_05(), "scene_09": es.spec_scene_09(), "lamps": lamp_spec()}
    for n in juice.templates():
        out["template_" + n] = juice.new_experiment("t_" + n, template=n)
    return out


def _field_specs():
    return {n: juice.load(p) for n, (p, _d) in fs.list_examples().items()}


def _roundtrip(spec, tmp_path, key):
    path = juice.write_script(spec, str(tmp_path / f"{key}_run.py"), solver="fast")
    return juice.build_from_script(path)


@pytest.mark.parametrize("key", sorted(_chamber_specs()))
def test_chamber_script_rebuilds_the_experiment(tmp_path, key):
    spec = _chamber_specs()[key]
    assert asdict(_roundtrip(spec, tmp_path, key)) == asdict(spec)


@pytest.mark.parametrize("key", sorted(_field_specs()))
def test_field_script_rebuilds_the_case(tmp_path, key):
    spec = _field_specs()[key]
    assert asdict(_roundtrip(spec, tmp_path, key)) == asdict(spec)


@pytest.mark.parametrize("key", ["lamps", "site_day"])
def test_script_writes_byte_identical_files(tmp_path, key):
    spec = {**_chamber_specs(), **_field_specs()}[key]
    back = _roundtrip(spec, tmp_path, key)
    a, b = juice.write(spec, str(tmp_path / "a")), juice.write(back, str(tmp_path / "b"))
    for k in ("enz", "easy"):
        assert filecmp.cmp(a[k], b[k], shallow=False), k


def test_scripts_explain_the_lamps_and_compare_a_field_case(tmp_path):
    lamp = juice.script_text(lamp_spec())
    assert "chambers/generic_lamps.chamber.json" in lamp and "jNO2 at full power" in lamp
    assert "chamber.json" not in juice.script_text(es.spec_scene_05())         # SAPHIR: the sun, nothing to explain
    field = juice.script_text(_field_specs()["site_day"])
    assert "juice.skill(case, cs_path)" in field and "# STEP 4 — MIXING AND DEPOSITION" in field
    assert field.index("case.initials") > field.index("# STEP 4")               # initial values on step 4, as in the window


def test_examples_build_what_they_describe():
    import importlib.util
    def build(name):
        s = importlib.util.spec_from_file_location("_ex_" + name, os.path.join(juice.HERE, "examples", name + ".py"))
        m = importlib.util.module_from_spec(s); s.loader.exec_module(m)
        return m.build()
    site = build("field_site")
    ref = juice.load(fs.list_examples()["site_day"][0])
    assert asdict(site) == asdict(ref)                                          # the shipped case, written out
    for name in ("saphir_isoprene", "uv_lamp_chamber"):
        assert not [m for lvl, m in juice.check(build(name)) if lvl == "error"], name


def test_new_experiment_is_what_picking_the_chamber_does():
    spec = juice.new_experiment("x", chamber="generic_lamps")
    ref = es.ExperimentSpec(name="x")
    es.apply_chamber(ref, "generic_lamps")
    assert asdict(spec) == asdict(ref)
    assert asdict(juice.new_experiment("y")) == asdict(es.ExperimentSpec(name="y"))   # SAPHIR = the dataclass defaults


def test_catalogues_list_what_can_be_chosen():
    assert "SAPHIR" in juice.chambers() and "generic_lamps" in juice.chambers()
    assert "MCM_v3.3.1" in juice.mechanisms()["bases"]
    assert set(juice.solvers()) == {"fast", "precise", "twostep"} == set(juice.OPTIONS["solver"])


def test_write_refuses_an_experiment_with_errors(tmp_path):
    spec = juice.new_experiment("broken")
    spec.conditions.date = ""
    assert any(lvl == "error" for lvl, _m in juice.check(spec))
    with pytest.raises(ValueError, match="errors"):
        juice.write(spec, str(tmp_path))


@pytest.fixture(scope="module")
def app():
    from PyQt6.QtWidgets import QApplication
    yield QApplication.instance() or QApplication([])


def test_chamber_window_prepare_writes_the_script(app, tmp_path):
    import juice_gui
    w = juice_gui.JuiceWindow("")
    w.load_spec(lamp_spec(), keep_name=True)
    w.page_start.out_edit.setText(str(tmp_path))
    w.page_run.kvpd.setChecked(False)
    assert w.page_run.prepare()
    held = w.collect_spec()
    script = tmp_path / held.name / (held.name + "_run.py")
    assert script.is_file()
    back = juice.build_from_script(str(script))
    assert asdict(back) == asdict(held)
    again = juice.write(back, str(tmp_path / "again"))
    for k, ext in (("enz", ".enz"), ("easy", ".easy")):
        assert filecmp.cmp(again[k], str(tmp_path / held.name / (held.name + ext)), shallow=False), k
    text = script.read_text(encoding="utf-8")
    assert f"solver = {w.page_run.solver.currentData()!r}" in text and "kv_pd = False" in text
    assert "# STEP 1 — EXPERIMENT" in text and "# STEP 5 — RUN" in text
    w.close()


def test_field_window_prepare_writes_the_script(app, tmp_path):
    import field_window as fw
    w = fw.FieldWindow()
    w.folder.setText(str(tmp_path))
    w.examples.setCurrentIndex(w.examples.findData("site_day"))
    w.dirty = False
    w._load_example()
    w.page_run.kvpd.setChecked(False)
    assert w.page_run.prepare()
    spec = w.case[0]
    script = os.path.join(os.path.dirname(w.page_run.json_path), spec.name + "_run.py")
    assert os.path.isfile(script)
    assert asdict(juice.build_from_script(script)) == asdict(spec)
    w.dirty = False
    w.close()


def test_catalogue_lists_every_planned_chemistry():
    import mechanism_registry as mreg
    cat = mreg.catalog()
    by_id = {e["id"]: e for e in cat["bases"]}
    assert by_id["MCM_v3.3.1"]["status"] == "available"
    assert by_id["CRACMM2"]["status"] in ("planned", "available")
    assert {"RACM2", "CB6r2", "CB7", "SAPRC07B", "GEOSCHEM_v12.07", "CB05", "CRACMM1_aq", "MCM_v3.2"} <= set(by_id)
    assert "NOAA_BB" in {e["id"] for e in cat["addons"]}
    for e in cat["bases"] + cat["addons"]:
        assert e.get("citation") and e.get("source"), e["id"]


def test_chemistry_step_shows_planned_entries_greyed_out(app):
    import juice_gui
    import mechanism_registry as mreg
    w = juice_gui.JuiceWindow("")
    combo = w.page_chem.base_combo if hasattr(w, "page_chem") else None
    if combo is None:
        from experiment_wizard import ChemistryPage
        combo = w.findChildren(ChemistryPage)[0].base_combo
    labels = [combo.itemText(i) for i in range(combo.count())]
    for e in mreg.planned("bases"):
        i = next(k for k, t in enumerate(labels) if t.startswith(e["name"]))
        assert not combo.model().item(i).isEnabled()
    assert combo.currentData() == "MCM_v3.3.1"
    w.close()


def test_species_search_lists_addon_species_and_ticks_the_addon(app):
    import juice_gui
    from experiment_wizard import ChemistryPage
    w = juice_gui.JuiceWindow("")
    page = w.findChildren(ChemistryPage)[0]
    page.search.setText("FURAN")
    items = [page.candidates.item(i).text() for i in range(page.candidates.count())]
    assert any(t.startswith("FURAN ") and "[add-on NOAA_BB]" in t for t in items)
    row = next(i for i, t in enumerate(items) if t.startswith("FURAN "))
    page.candidates.item(row).setSelected(True)
    page._add_selected()
    assert "NOAA_BB" in page.addon_ids()
    assert "FURAN" in [page.chosen.item(i).text() for i in range(page.chosen.count())]
    w.close()


def test_species_lookup_is_per_base_mechanism():
    import species_db as sdb
    assert sdb.load_species_db().get("G2T70GE") is None                 # GECKO-A codes are not global names
    import mechanism_registry as mreg
    bases = mreg.load_bases()
    if {"GECKOA_toluene", "GECKOA_apinene"} <= set(bases):
        assert sdb.load_species_db(base_id="GECKOA_toluene").get("G2T70GE").formula == "C7H8NO10"
        assert sdb.load_species_db(base_id="GECKOA_apinene").get("G2T70GE").formula == "C7H7N2O11"
    if "CRACMM2" in bases:
        assert sdb.load_species_db(base_id="CRACMM2").get("ISO").formula == "C5H8"
    assert sdb.load_species_db().get("FURAN").formula == "C4H4O"            # add-on species from its manifest
