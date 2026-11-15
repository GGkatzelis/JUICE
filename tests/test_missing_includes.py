"""A missing include file stops the preparation (it used to run a fragment of the mechanism and report success), and a
study rebuilds a chemistry folder that a failed rebuild left empty behind its marker."""
import json
import os
import shutil

import pytest

import juice
from juice.core import scenario_prep as sp
from juice.core import study


def _exp(name="inc"):
    e = juice.new_experiment(name, chamber="SAPHIR")
    c = e.conditions
    c.source, c.date, c.start, c.end = "manual", "2024-07-09", "08:00", "08:30"
    c.temperature_mode, c.j_source, c.roof_open = "const", "clear_sky", [("08:10", "08:30")]
    e.precursors = ["C2H5OH"]
    e.injections = [juice.Injection(name="ethanol", products="C2H5OH", total_ppb="10", start="08:10")]
    return e


def test_missing_include_stops_the_preparation(tmp_path):
    e = _exp()
    files = juice.write(e, str(tmp_path))
    easy = files["easy"]
    text = open(easy, encoding="latin-1").read()
    open(easy, "w", encoding="latin-1").write(text.replace("#include saphir.walldep.incl", "#include no_such_file.incl"))
    with pytest.raises(FileNotFoundError, match="no_such_file.incl"):
        sp.prepare_scenario(files, sp.load_paths("")["include_dir"], export_kv_pd=False)
    open(easy, "w", encoding="latin-1").write(text)                     # the original prepares
    assert os.path.isfile(sp.prepare_scenario(files, sp.load_paths("")["include_dir"], export_kv_pd=False))


def test_study_rebuilds_an_empty_chemistry_folder(tmp_path):
    e = _exp("base")
    st = study.Study(name="s", folder=str(tmp_path / "study"))
    out = os.path.join(st.folder, "chemistry", e.name)
    os.makedirs(os.path.join(out, "chemistry"))
    json.dump({"mechanism_stem": "gone", "mech_incl_suffix": ".incl", "use_anorg": False, "extra_includes": [],
               "addons": []}, open(os.path.join(out, "chemistry.json"), "w"))          # a marker with nothing behind it
    chem, fields = study._chemistry_for(st, str(tmp_path / "base.experiment.json"), e)
    assert fields["mechanism_stem"] != "gone" and os.path.isfile(os.path.join(chem, fields["mechanism_stem"] + ".mech"))
    chem2, fields2 = study._chemistry_for(st, str(tmp_path / "base.experiment.json"), e)  # complete now: reused
    assert fields2 == fields
