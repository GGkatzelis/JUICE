"""An experiment set up for FACSIMILE (.easy + ENZ) read back into an ExperimentSpec: the EASY the
GUI regenerates from it is statement-equivalent to the original and the ENZ columns are identical."""
import os

import numpy as np
import pytest

from juice.importers import facsimile_import as fi
from juice.core import measured_inputs as mi
from juice.core.experiment_spec import enz_columns, render_easy, spec_from_dict, spec_to_dict
from juice.core.easy_scene import statements_equivalent

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
IDL = os.path.join(HERE, "validation", "IDL Input")


def _scene(name):
    return os.path.join(IDL, name, name + ".easy"), os.path.join(IDL, name, name + ".enz")


@pytest.mark.parametrize("scene", ["scene_05", "scene_09"])
def test_golden_scenes_read_back_without_loss(scene, tmp_path):
    easy, enz = _scene(scene)
    spec, notes = fi.import_facsimile(easy, enz)
    # nothing the import could not express — only notes about the inputs themselves (a 0 pulse column)
    assert not [n for n in notes if n.startswith(("EASY:", "ENZ:"))], notes
    out = tmp_path / "re.easy"; out.write_text(render_easy(spec), encoding="latin-1")
    ok, diffs = statements_equivalent(easy, str(out), ignore_duplicate_reactions=True)
    assert ok, diffs
    times, raw = mi.read_enz(enz)
    grid, cols, _n = enz_columns(spec)
    assert np.allclose(grid, times)
    for name, vals in raw.items():
        if name in cols:
            assert np.allclose(vals, cols[name], atol=1e-9), name
    assert spec.conditions.source == "file" and spec.conditions.input_file == os.path.abspath(enz)
    assert spec.easy_version == "6.00"                          # from the .fac next to the scene
    assert spec.reference_cs.endswith(scene + ".cs")
    assert os.path.join(IDL, scene) in spec.include_dirs


def test_scene_09_details():
    spec, _notes = fi.import_facsimile(*_scene("scene_09"))
    c = spec.conditions
    assert (c.date, c.start, c.end, c.step_s) == ("2024-08-05", "04:00", "21:59", 60)
    assert c.roof_open == [("09:33", "18:37")] and c.j_source == "measured" and c.temperature_mode == "enz"
    assert spec.mechanism_stem == "chanel.all_mixes.fac" and spec.chamber.light_source
    by = {i.name: i for i in spec.injections}
    assert by["H2O"].start == "05:26" and by["H2O"].duration_min == 45 and len(by["H2O"].extra_channels) == 3
    assert by["O3"].total_ppb == "52" and by["O3"].start == "09:36"
    assert by["seed"].products == "" and by["seed"].start == "06:28"           # recorded in the ENZ only
    assert by["ethanol"].action == "voc" and by["H2O"].action == "humidification"
    assert ("CO", "M*1e-9*14.0") in spec.initials


def test_scene_05_delayed_release_and_constant_temperature():
    spec, notes = fi.import_facsimile(*_scene("scene_05"))
    assert spec.conditions.temperature_mode == "const" and spec.conditions.temperature_K == 300
    vcp2 = next(i for i in spec.injections if i.name == "VCP2")
    assert vcp2.switch == "qVOC" and vcp2.column == "qVCP2" and len(vcp2.chain) == 3
    assert any("QO3" in n and "0 throughout" in n for n in notes)              # O3 = 113*0, never pulsed


def test_pulse_of():
    t = np.arange(0, 600, 60.0)
    v = np.zeros(10); v[3:5] = 0.5
    assert fi.pulse_of(v, t, 60) == ("00:03", 2.0, True)
    v2 = np.zeros(10); v2[3] = 1.0; v2[6] = 1.0                                 # two spikes: not one pulse
    assert fi.pulse_of(v2, t, 60) == ("00:03", 4.0, False)
    assert fi.pulse_of(np.zeros(10), t, 60) == ("", 0.0, True)


def test_name_from_column():
    assert [fi.name_from_column(c) for c in ("QNO", "Q_ethen", "qVCP2", "QH2O", "QCl2", "Q_seed")] == \
           ["NO", "ethen", "VCP2", "H2O", "Cl2", "seed"]


def test_enz_as_measured_inputs_file():
    _easy, enz = _scene("scene_09")
    inp = mi.read_inputs_file(enz)
    assert {"t", "H2O", "dilute", "hv", "sz"} <= set(inp.cols)
    assert "injection switch" in inp.summary()


def test_spec_round_trips_through_json():
    spec, _ = fi.import_facsimile(*_scene("scene_09"))
    back = spec_from_dict(spec_to_dict(spec))
    assert back.include_dirs == spec.include_dirs and back.easy_version == spec.easy_version
    assert back.reference_cs == spec.reference_cs and [i.name for i in back.injections] == [i.name for i in spec.injections]
