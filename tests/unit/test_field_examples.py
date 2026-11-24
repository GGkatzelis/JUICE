"""The shipped synthetic field examples: complete, reproducible, and loaded as a copy."""
import os

import pytest

from juice.core import field_scenario as fsc
from juice.core import field_spec as fs
from juice.core.scenario_prep import load_paths

EXAMPLES = fs.list_examples()


def test_three_examples_one_per_kind():
    assert set(EXAMPLES) == {"site_day", "plume_ageing", "la_parcel"}
    kinds = {fs.load_spec(p).run_type for p, _d in EXAMPLES.values()}
    assert kinds == set(fs.RUN_TYPES)


@pytest.mark.parametrize("key", sorted(EXAMPLES))
def test_example_is_complete_and_prepares(key, tmp_path):
    spec = fs.load_spec(EXAMPLES[key][0])
    table, rep = fs.check(spec)
    assert rep.ok, rep.text()
    assert not rep.warnings or all("jNO2" not in w for w in rep.warnings), rep.text()
    path, res = fsc.prepare(spec, table, str(tmp_path / key), load_paths("")["include_dir"], "", export_kv_pd=False)
    assert os.path.isfile(path) and os.path.isfile(res["enz"]) and os.path.isfile(res["easy"])
    assert not [n for n in res["notes"] if "empty" in n or "not used" in n], res["notes"]


def test_examples_are_what_the_script_writes(tmp_path):
    """make_examples.py reproduces the shipped files byte for byte."""
    import importlib.util
    spec = importlib.util.spec_from_file_location("make_examples", os.path.join(fs.EXAMPLES_DIR, "make_examples.py"))
    mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    for key, make in mod.EXAMPLES.items():
        made = make(str(tmp_path / key))
        for name in (key + ".csv", key + fs.SPEC_SUFFIX):
            shipped = open(os.path.join(fs.EXAMPLES_DIR, key, name), encoding="utf-8", newline=None).read()
            written = open(os.path.join(tmp_path, key, name), encoding="utf-8", newline=None).read()
            assert shipped == written, name                 # universal newlines: a CRLF checkout is the same content


def test_copy_case_leaves_the_example_alone(tmp_path):
    src = EXAMPLES["plume_ageing"][0]
    before = open(src, encoding="utf-8").read()
    p = fs.copy_case(src, str(tmp_path))
    assert p == str(tmp_path / "plume_ageing" / "plume_ageing.field.json")
    spec = fs.load_spec(p)
    assert spec.table == str(tmp_path / "plume_ageing" / "plume_ageing.csv") and os.path.isfile(spec.table)
    assert open(src, encoding="utf-8").read() == before
    assert fs.is_example(src) and not fs.is_example(p)
