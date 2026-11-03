"""Input side: the MCM download, the EASY writer and the experiment description → ENZ + EASY.

These are the self-tests that already existed as scripts (python <module> --selftest); here they
run under pytest, one test per reference scene where the script loops over scenes.
"""
import contextlib
import io

import pytest

from conftest import SCENES


def test_mcm_download_unmodified():
    """Every base mechanism file matches the sha256 recorded when it was downloaded."""
    import mechanism_registry as mreg
    bases = mreg.load_bases()
    assert bases, "no base mechanism installed under mechanisms/base"
    problems = [p for b in bases.values() for p in b.verify()]
    assert not problems, "\n".join(problems)


def test_addons_declare_a_compatible_base():
    import mechanism_registry as mreg
    bases = mreg.load_bases()
    for a in mreg.load_addons().values():
        assert any(a.compatible(b) for b in bases), f"add-on {a.id} fits no installed base"
        assert a.status in ("draft", "reviewed"), f"add-on {a.id}: unknown status {a.status!r}"


@pytest.mark.parametrize("scene", SCENES)
def test_easy_roundtrip(scene, tmp_path):
    """parse → render is byte-identical; extract → write is statement-equivalent; same scenario JSON."""
    import check_easy_roundtrip as rt
    msg = rt.check(scene, str(tmp_path))
    assert msg.startswith("OK"), msg


@pytest.mark.parametrize("scene", SCENES)
def test_easy_builder_form(scene):
    """Loading a scene into the EASY builder form and writing it back is statement-equivalent."""
    import easy_builder as eb
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        rc = eb._selftest([scene])
    assert rc == 0, out.getvalue()


def test_experiment_spec_goldens():
    """scene_05 and scene_09 rebuilt from experiment descriptions give the IDL-written EASY and ENZ."""
    import experiment_spec as es
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        rc = es._selftest()
    assert rc == 0, out.getvalue()


def test_chamber_profiles_load():
    """Every chambers/*.chamber.json parses into a profile with a positive volume."""
    import glob
    import os
    import chambers
    files = glob.glob(os.path.join("chambers", "*.chamber.json"))
    loaded = chambers.load_dir()
    assert files and len(loaded) >= len(files), f"{len(files)} files, {len(loaded)} profiles loaded"
    for cid, p in loaded.items():
        assert p.volume_m3 > 0, f"{cid}: volume {p.volume_m3}"


def test_chambers_selftest():
    import chambers
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        rc = chambers._selftest()
    assert rc == 0, out.getvalue()
