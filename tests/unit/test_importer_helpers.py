"""Readers of other programs' files, function by function: the FACSIMILE result files of the comparison tool, the
rate-constant and photolysis files of F0AM, the ENZ a FACSIMILE experiment names, and the generation count of a
GECKO-A mechanism."""
import pytest

from juice.importers import f0am_export, f0am_import, facsimile_import
from juice.importers import gecko_mechanism as gm
from juice.tools import compare_juice_facsimile as cmp


# ── comparison tool: reaction names ──────────────────────────────────────────

@pytest.mark.parametrize("written, canonical", [
    ("O3+NO-->NO2", "NO+O3-->NO2"),
    ("no + o3 --> no2", "NO+O3-->NO2"),
    ("N2O5-->NA+NA", "N2O5-->2*NA"),
    ("N2O5-->2*NA", "N2O5-->2*NA"),
    ("A-->0.5*B+0.5*B+0.3*C", "A-->B+0.3*C"),
    ("X+DIL-->", "DIL+X-->"),
    ("X + DIL --> (NONE)", "DIL+X-->"),
    ("JSCALE", "JSCALE"),
])
def test_reaction_names_are_made_comparable(written, canonical):
    assert cmp._canon_rxn_str(written) == canonical


def test_rate_constant_column_labels():
    assert cmp._canon_kv_label("k_R00006:O3+NO-->NO2") == "NO+O3-->NO2"
    assert cmp._canon_kv_label("v_O3+NO-->NO2") == "NO+O3-->NO2"


@pytest.mark.parametrize("name, setting", [
    ("T", True), ("m", True), ("KMT08", True), ("KRO2NO", True), ("J4", True), ("J11", True), ("jno2_cs", True),
    ("DILUTE", True), ("OH", False), ("JUICE", False), ("NO2", False), ("J", False), ("C5H8", False),
])
def test_settings_are_told_from_species(name, setting):
    assert cmp._is_setting(name) is setting


# ── comparison tool: FACSIMILE files ─────────────────────────────────────────

CS = """ * Facsimile-Output
 COLUMN 1=TIME / SEC SINCE 1.1.2000 00:00:00
 COLUMN 2=OH<0>
 COLUMN 3=T<0>
 NUMBER OF COLUMNS=3
 &&&&&&&&
 773812860.0      1.5E+06  2.9276631E+02
 773812920.0      1.540849-169  293.0
 773812920.0      9.9E+06  0.0
"""
KV = """ COLUMN 1=TIME
 COLUMN 2=k<1,0> !UO3+NO-->NO2!N
 COLUMN 3=v<1,0> !UO3+NO-->NO2!N
 COLUMN 4=k<2,0> !UN2O5-->NA+NA!N
 NUMBER OF COLUMNS=4
 &&&&&&&&
 773812860.0  1.8E-14  2.5E+06  3.0E-04
 773812920.0  1.9E-14  2.6E+06  3.1E-04
"""
PD = """ COLUMN 1=TIME
 COLUMN 2=p<0> oh
 COLUMN 3=d<0> oh
 NUMBER OF COLUMNS=3
 &&&&&&&&
 773812860.0  1.0E+06  1.1E+06
 773812920.0  2.0E+06  2.1E+06
"""


def written(tmp_path, text, name):
    path = tmp_path / name
    path.write_text(text, encoding="latin-1")
    return str(path)


def test_facsimile_concentration_file(tmp_path):
    df = cmp.load_facsimile_cs(written(tmp_path, CS, "s.cs"))
    assert list(df.columns) == ["OH", "T"] and list(df.index) == [773812860.0, 773812920.0]
    assert df.loc[773812920.0, "OH"] == 1.540849e-169                   # the E that FACSIMILE leaves out
    assert df.loc[773812920.0, "T"] == 293.0                            # of the repeated last time, the first row


def test_facsimile_rate_constant_file(tmp_path):
    k, v = cmp.load_facsimile_kv(written(tmp_path, KV, "s.kv"))
    assert list(k.columns) == ["NO+O3-->NO2", "N2O5-->2*NA"] and list(v.columns) == ["NO+O3-->NO2"]
    assert k.loc[773812920.0, "NO+O3-->NO2"] == 1.9e-14 and v.loc[773812860.0, "NO+O3-->NO2"] == 2.5e6


def test_facsimile_production_and_destruction_file(tmp_path):
    df = cmp.load_facsimile_pd(written(tmp_path, PD, "s.pd"))
    assert list(df.columns) == ["P_OH", "D_OH"] and df.loc[773812920.0, "D_OH"] == 2.1e6


def test_facsimile_files_without_data_are_refused(tmp_path):
    header = " COLUMN 1=TIME\n COLUMN 2=k<1,0> !UA-->B!N\n &&&&\n"
    with pytest.raises(ValueError):
        cmp.load_facsimile_kv(written(tmp_path, header, "e.kv"))
    with pytest.raises(ValueError):
        cmp.load_facsimile_pd(written(tmp_path, header, "e.pd"))
    with pytest.raises(ValueError):
        cmp.load_facsimile_cs(written(tmp_path, "nothing here\n", "e.cs"))


# ── F0AM ─────────────────────────────────────────────────────────────────────

@pytest.mark.parametrize("matlab, easy", [
    ("2.7e-12.*exp(360./T);", "2.7e-12*EXP(360/T)"),
    ("k0 .* (T./300).^-2.6", "k0*(T/300)@(-2.6)"),
    ("1e-11*(300./T)^1.5", "1e-11*(300/T)@1.5"),
    ("sqrt(a) + log10(b)", "SQRT(a)+LOG10(b)"),
])
def test_matlab_rate_in_easy_notation(matlab, easy):
    assert f0am_import.matlab_to_easy(matlab) == easy


def test_helper_variables_are_substituted_in_brackets():
    local = {"k0": "3e-30*M", "kinf": "(k0)*2"}
    assert f0am_import._substitute("k0/kinf+k0x", local) == "(3e-30*M)/(((3e-30*M))*2)+k0x"


@pytest.mark.parametrize("rate, zero", [("0", True), ("0.0", True), ("(0.0)", True), ("a*1.0*0.0", True), ("0.0*b", True),
                                         ("1e-12", False), ("10.0*a", False), ("0.05", False)])
def test_a_rate_that_is_zero(rate, zero):
    assert f0am_import._is_zero(rate) is zero


K_FILE = """% rate constants
i = 0;
i=i+1;
Knames{i} = 'KRO2NO';
krx(:,i) = 2.7e-12.*exp(360./T);   % MCM

i = i + 1;
Knames{i} = 'KTROE';
k0 = 3.2e-30.*M.*(T./300).^-4.5; kinf = 3.0e-11;
krx(:,i) = k0.*kinf./(k0+kinf);

i=i+1;
% Knames{i} = 'COMMENTED_OUT';
krx(:,i) = 1;
"""


def test_named_rate_constants_of_an_f0am_file(tmp_path):
    k = f0am_import.read_k(written(tmp_path, K_FILE, "M_K.m"))
    assert list(k) == ["KRO2NO", "KTROE"]
    assert k["KRO2NO"] == "2.7e-12*EXP(360/T)"
    assert k["KTROE"] == "(3.2e-30*M*(T/300)@(-4.5))*(3.0e-11)/((3.2e-30*M*(T/300)@(-4.5))+(3.0e-11))"
    assert f0am_import.read_k(str(tmp_path / "missing_K.m")) == {}


def test_photolysis_names_of_an_f0am_file(tmp_path):
    text = "J.NO2 = Jmcm.J4;\n  J.O3_O1D = 0.9 .* Jmcm.J1; % comment\n%J.OLD = Jmcm.J9;\nJ.HONO2 = Jmcm.Jn12 + Jmcm.J8;\n"
    assert f0am_import.read_j(written(tmp_path, text, "M_J.m")) == {"NO2": "J4", "O3_O1D": "0.9*J1", "HONO2": "JN12+J8"}


def test_reaction_sides_for_the_f0am_export():
    assert f0am_export._parse_side(["oh", " 2*HO2 ", "", "0.5*ch3o2"]) == [(1.0, "OH"), (2.0, "HO2"), (0.5, "CH3O2")]


# ── a FACSIMILE experiment's ENZ ─────────────────────────────────────────────

def test_enz_named_in_the_easy_file(tmp_path):
    easy = written(tmp_path, "; FILES[ENZ] = old.enz\n*FILES[ENZ] = older.enz\nT = INPUT(t)\n files [ enz ]  =  \"day 1.enz\"  ; today\n", "a.easy")
    assert facsimile_import.enz_name_in(easy) == "day 1.enz"
    assert facsimile_import.enz_name_in(written(tmp_path, "T = INPUT(t)\n", "b.easy")) == ""
    assert facsimile_import.enz_name_in(str(tmp_path / "missing.easy")) == ""


def test_which_enz_is_used(tmp_path):
    for name in ("named.enz", "scene.enz"):
        (tmp_path / name).write_text("", encoding="utf-8")
    folder = str(tmp_path)
    assert facsimile_import.find_enz(folder, "scene", "named.enz").endswith("named.enz")       # the one the EASY names
    assert facsimile_import.find_enz(folder, "scene", "gone.enz").endswith("scene.enz")        # else <stem>.enz
    assert facsimile_import.find_enz(folder, "other", "") == ""                                # two candidates: none
    (tmp_path / "named.enz").unlink()
    assert facsimile_import.find_enz(folder, "other", "").endswith("scene.enz")                # the only one there
    assert facsimile_import.find_enz(str(tmp_path / "nowhere"), "x", "") == ""


# ── GECKO-A generations ──────────────────────────────────────────────────────

def species(code, radical=False, carbons=1):
    return gm.GSpecies(code, code, 30.0, radical, carbons, 2, 0, 1)


def reaction(reactants, products, line, kind="arrh"):
    return gm.GReaction([(1.0, c) for c in reactants], [(1.0, c) for c in products], 1e-12, 0.0, 0.0, kind, line_no=line)


def small_mechanism():
    """VOC -(OH)-> RO2 -> P1 -(OH)-> RO2B -> P2 -(OH)-> P3 ;  P1 also straight to P3 ;  ORPHAN -> P2 is unreachable."""
    sp = {c: species(c) for c in ("VOC", "P1", "P2", "P3", "ORPHAN")}
    sp.update({c: species(c, radical=True) for c in ("RO2", "RO2B")})
    sp["HO"] = gm.GSpecies("HO", "HO", 17.0, True, 0, 1, 0, 1)
    sp["CO"] = gm.GSpecies("CO", "CO", 28.0, False, 1, 0, 0, 1)
    rx = [reaction(["VOC", "HO"], ["RO2"], 1), reaction(["RO2"], ["P1", "CO"], 2), reaction(["P1", "HO"], ["RO2B"], 3),
          reaction(["RO2B"], ["P2"], 4), reaction(["P2", "HO"], ["P3"], 5), reaction(["P1"], ["P3"], 6),
          reaction(["ORPHAN"], ["P2"], 7), reaction(["HO", "HO"], ["HO"], 8), reaction(["P1"], ["P1"], 9, kind="transfer")]
    return gm.GeckoMechanism("", sp, rx, {"P1": (1e-5, 80.0), "ORPHAN": (1e-9, 100.0)}, 298.0, {1: ["RO2", "RO2B", "ORPHAN"]}, ["VOC"])


def test_generation_is_counted_at_each_stable_species():
    gen = gm.generations(small_mechanism())
    assert gen == {"VOC": 0, "RO2": 1, "P1": 1, "RO2B": 2, "P2": 2, "P3": 2}       # P3 by the shortest way; CO not counted
    assert "ORPHAN" not in gen


def test_cut_stops_the_stable_species_beyond_the_generation_and_drops_what_cannot_be_reached():
    cut, gen = gm.generation_cut(small_mechanism(), 1)
    assert [r.line_no for r in cut.reactions] == [1, 2, 3, 4, 6, 8, 9]              # P2 (generation 2) reacts no further
    assert set(cut.species) == {"VOC", "RO2", "P1", "RO2B", "P2", "P3", "HO", "CO"} and "ORPHAN" not in cut.species
    assert cut.pvap == {"P1": (1e-5, 80.0)} and cut.pero == {1: ["RO2", "RO2B"]} and gen["P2"] == 2
    cut0, _ = gm.generation_cut(small_mechanism(), 0)
    assert [r.line_no for r in cut0.reactions] == [1, 2, 8]                         # P1 stays as a product
    everything, _ = gm.generation_cut(small_mechanism(), 9)
    assert [r.line_no for r in everything.reactions] == [1, 2, 3, 4, 5, 6, 8, 9]    # only the orphan's reaction is gone
