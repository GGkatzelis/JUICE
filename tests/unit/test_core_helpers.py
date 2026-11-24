"""Small helpers of the core: injection rates of the EASY builder, the ENZ column list, the fall-off function,
the species search for surrogates, saving a chamber profile, and the parse cache."""
import copy
import json
import os

import numpy as np
import pytest

from juice.chemistry import mcm_species_db as db
from juice.core import chambers, easy_builder as ebd, mcm_kinetics, parse_cache
from juice.core.mechanism import Mechanism


# ── EASY builder ─────────────────────────────────────────────────────────────

@pytest.mark.parametrize("entry, written", [
    ("8.1", "CONST(M*1e-9/60*8.1)"),                       # a number is ppb per minute
    (" 575*180/175 ", "CONST(M*1e-9/60*575*180/175)"),     # also as arithmetic
    ("M*1e-9/60*qNO", "CONST(M*1e-9/60*qNO)"),             # anything with a name is the expression itself
    ("CONST(1./60)", "CONST(1./60)"),                      # already written out
    ("const (2)", "const (2)"),
])
def test_injection_rate_as_written_into_the_easy_file(entry, written):
    assert ebd._rate_expr(entry) == written


def test_chain_reaction_rate_is_always_the_expression():
    assert ebd._rate_expr("1./60", ppb=False) == "CONST(1./60)"


@pytest.mark.parametrize("entry", ["8.1", "575*180/175", "M*1e-9/60*qNO", "jno2/5e-3/3600*HV_in"])
def test_injection_rate_reads_back_as_it_was_entered(entry):
    assert ebd._ppb_from_rate(ebd._rate_expr(entry)) == entry


def test_reading_a_rate_back():
    assert ebd._ppb_from_rate("CONST( M*1e-9 / 60 * 113*0 )") == "113*0"
    assert ebd._ppb_from_rate("CONST(M*1e-9/60*8.1)", ppb=False) == "M*1e-9/60*8.1"
    assert ebd._ppb_from_rate("") == "" and ebd._ppb_from_rate(None) == ""
    assert ebd._ppb_from_rate("PHTFKT(1,2,3)") == "PHTFKT(1,2,3)"       # not a CONST: left as it is
    assert ebd._const_inner(" CONST( 2.8e6 ) ") == "2.8e6" and ebd._const_inner("INPUT(t)") == "INPUT(t)"


def test_enz_columns_without_time_and_without_what_follows_the_header(tmp_path):
    path = tmp_path / "k.enz"
    path.write_text("COLUMN 1=TIME / SEC SINCE 1.1.2000 00:00:00\nCOLUMN 2=sz\nCOLUMN 3=jNO2 / 1/s\nNUMBER OF COLUMNS=3\n"
                    "&&&&&&&\nCOLUMN 4=not a column, this is data\n", encoding="latin-1")
    assert ebd.read_enz_columns(str(path)) == ["sz", "jNO2"]
    assert ebd.read_enz_columns(str(tmp_path / "missing.enz")) == []


def test_mechanisms_are_found_by_their_mech_file(tmp_path):
    (tmp_path / "sub").mkdir()
    for name in ("ethanol.fac.mech", "ethanol.fac.rate", "sub/limonene.fac.mech", "notes.mech"):
        (tmp_path / name).write_text("", encoding="utf-8")
    found = ebd.scan_mechanisms(str(tmp_path))
    assert found == {"ethanol.fac": str(tmp_path), "limonene.fac": str(tmp_path / "sub")}
    assert ebd.scan_mechanisms("") == {} and ebd.scan_mechanisms(str(tmp_path / "none")) == {}


# ── fall-off function ────────────────────────────────────────────────────────

def test_comprc_is_the_iupac_fall_off_expression():
    # OH + NO2 at 300 K and M = 2.46e19 (MCM KMT08): k0·M = 7.872e-11, k∞ = 3e-11, Fc = 0.41 → 9.75e-12 by hand
    k = mcm_kinetics.COMPRC(3.2e-30 * 2.46e19, 3.0e-11, 0.41)
    assert k == pytest.approx(9.75e-12, rel=0.005)


def test_comprc_limits():
    assert mcm_kinetics.COMPRC(1e-20, 3e-11, 0.41) == pytest.approx(1e-20, rel=1e-3)     # low pressure: k0·M
    assert mcm_kinetics.COMPRC(1e-2, 3e-11, 0.41) == pytest.approx(3e-11, rel=1e-3)      # high pressure: k∞
    assert mcm_kinetics.COMPRC(0.0, 3e-11, 0.41) == 0.0
    both = mcm_kinetics.COMPRC(np.array([1e-20, 1e-2]), 3e-11, 0.41)
    assert both.shape == (2,)                                                             # works on series


def test_generic_ro2_no_rate():
    assert mcm_kinetics.KRO2(350.0) == pytest.approx(2.9e-12 * np.e)


# ── species search for surrogates ────────────────────────────────────────────

needs_db = pytest.mark.skipif(not db.MCM_DB, reason="the MCM species table is not there")


@needs_db
def test_formula_finds_the_species_of_that_formula_first():
    hits = db.search_species("C10H16")
    names = [sp.name for sp, _why in hits]
    assert {"APINENE", "BPINENE", "LIMONENE"} <= set(names[:6])
    assert all(why == "same formula" for sp, why in hits if sp.name in ("APINENE", "BPINENE", "LIMONENE"))


@needs_db
def test_compound_outside_the_mcm_gets_the_nearest_chain_of_the_same_kind():
    hits = db.search_species("C14H30")                                  # tetradecane: the MCM stops at dodecane
    assert hits[0][0].name == "NC12H26" and hits[0][1] == "nearest carbon number with the same O"


@needs_db
def test_part_of_a_name_and_radicals_last():
    hits = db.search_species("pinene")
    names = [sp.name for sp, _why in hits]
    assert names[:2] == ["APINENE", "BPINENE"] and hits[0][1] == "name match"
    radicals = [i for i, (sp, why) in enumerate(hits) if "radical" in why]
    stable = [i for i, (sp, why) in enumerate(hits) if "radical" not in why]
    assert not radicals or max(stable) < min(radicals)


@needs_db
def test_empty_query_and_the_limit():
    assert db.search_species("") == [] and db.search_species(None) == []
    assert len(db.search_species("C5H8O2", limit=3)) <= 3


# ── chamber profile ──────────────────────────────────────────────────────────

def test_profile_is_saved_under_a_safe_file_name_and_the_registry_is_left_alone(tmp_path):
    before = set(chambers.CHAMBERS)
    profile = copy.deepcopy(chambers.chamber("SAPHIR"))
    profile.id, profile.name, profile.volume_m3 = "", "My test chamber #1", 12.5
    path = chambers.save_profile(profile, str(tmp_path))
    assert os.path.dirname(path) == str(tmp_path) and path.endswith(chambers.SUFFIX)
    assert profile.id == chambers.safe_id("My test chamber #1") and os.path.basename(path) == profile.id + chambers.SUFFIX
    saved = json.load(open(path, encoding="utf-8"))
    assert saved["volume_m3"] == 12.5 and saved["name"] == "My test chamber #1"
    assert set(chambers.CHAMBERS) == before                             # another folder: not among the installed ones
    assert chambers.chamber("no such chamber").id == chambers.DEFAULT_CHAMBER


# ── parse cache ──────────────────────────────────────────────────────────────

class FakeParser:
    """Counts how often the text is read and parsed; the 'mechanism' is the number of lines.  The counts are kept
    outside the object's own state, because a cache hit replaces that state with the one that was stored."""

    def __init__(self):
        self.mech = None

    reads = property(lambda self: FakeParser.counts[id(self)]["reads"])
    parses = property(lambda self: FakeParser.counts[id(self)]["parses"])
    counts = {}

    def __new__(cls):
        obj = super().__new__(cls)
        FakeParser.counts[id(obj)] = {"reads": 0, "parses": 0}
        return obj

    def preprocess(self, path):
        FakeParser.counts[id(self)]["reads"] += 1
        return open(path, encoding="utf-8").read()

    def parse(self, path):
        FakeParser.counts[id(self)]["parses"] += 1
        print("parser says: one warning")
        text = self.preprocess(path)
        self.mech = Mechanism(species=[f"S{text.count(chr(10))}"])
        return self.mech


@pytest.fixture
def cache_home(tmp_path, monkeypatch):
    monkeypatch.setattr(os.path, "expanduser", lambda p: str(tmp_path / "home") if p == "~" else p)
    return tmp_path / "home" / ".juice" / "parse_cache"


def elog(tmp_path, megabytes, name="big.elog"):
    path = tmp_path / name
    path.write_text("A --> B\n" * int(megabytes * 1e6 / 8), encoding="utf-8")
    return str(path)


SOURCE = os.path.join(os.path.dirname(parse_cache.__file__), "scenario_runner.py")


def test_small_file_is_parsed_without_a_cache(tmp_path, cache_home):
    parser, messages = FakeParser(), []
    mech = parse_cache.parse_cached(parser, elog(tmp_path, 0.01), log=lambda *a, **k: messages.append(a[0]), source_file=SOURCE)
    assert parser.parses == 1 and mech is parser.mech and not cache_home.exists()


def test_large_file_is_parsed_once_and_then_taken_from_the_cache(tmp_path, cache_home):
    path = elog(tmp_path, 5.1)
    first, messages = FakeParser(), []
    log = lambda *a, **k: messages.append(a[0])
    mech = parse_cache.parse_cached(first, path, log=log, source_file=SOURCE)
    assert first.parses == 1 and first.reads == 1                       # the text is read once, not twice
    assert len(list(cache_home.glob("*.pkl"))) == 1 and any("parse cache written" in m for m in messages)
    assert any("parser says: one warning" in m for m in messages)
    again, messages2 = FakeParser(), []
    mech2 = parse_cache.parse_cached(again, path, log=lambda *a, **k: messages2.append(a[0]), source_file=SOURCE)
    assert again.parses == 0 and mech2.species == mech.species and again.mech is mech2
    assert any("taken from the parse cache" in m for m in messages2)
    assert any("parser says: one warning" in m for m in messages2)      # the parser's messages come with the entry


def test_changed_text_is_parsed_again_and_the_cache_can_be_switched_off(tmp_path, cache_home):
    path = elog(tmp_path, 5.1)
    parse_cache.parse_cached(FakeParser(), path, log=lambda *a, **k: None, source_file=SOURCE)
    with open(path, "a", encoding="utf-8") as fh:
        fh.write("C --> D\n")
    changed = FakeParser()
    parse_cache.parse_cached(changed, path, log=lambda *a, **k: None, source_file=SOURCE)
    assert changed.parses == 1 and len(list(cache_home.glob("*.pkl"))) == 2
    off = FakeParser()
    parse_cache.parse_cached(off, path, enabled=False, log=lambda *a, **k: None, source_file=SOURCE)
    assert off.parses == 1 and len(list(cache_home.glob("*.pkl"))) == 2


def test_broken_cache_entry_falls_back_to_parsing(tmp_path, cache_home):
    path = elog(tmp_path, 5.1)
    parse_cache.parse_cached(FakeParser(), path, log=lambda *a, **k: None, source_file=SOURCE)
    (entry,) = cache_home.glob("*.pkl")
    entry.write_bytes(b"not a pickle")
    parser, messages = FakeParser(), []
    mech = parse_cache.parse_cached(parser, path, log=lambda *a, **k: messages.append(a[0]), source_file=SOURCE)
    assert parser.parses == 1 and mech is parser.mech
    assert any("parse cache not used" in m for m in messages)
