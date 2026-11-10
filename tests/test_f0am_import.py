"""f0am_import.py: F0AM mechanisms as JUICE base mechanisms (CRACMM2 first).

The converter is checked on a small synthetic F0AM file (every construct the library uses), on the
installed CRACMM2 base (counts, the J-map fix, no undefined name in any rate), on the EASY files an
experiment with it gets (sun and lamps), and with one run in which an MCM-named injection (C5H8)
reaches the mechanism's own species (ISO) through the alias file.  The conversion itself was checked
against F0AM running its own cracmm2_*.m (dark and lit with the same J: median 0.1 % / 0.05 %)."""
import os
import re

import pytest

import experiment_spec as es
import f0am_import as fi
import juice
import mechanism_registry as mreg

SYNTH = r"""
SpeciesToAdd = {'A'; 'B'; 'C'; 'HO'; 'D'};
AddSpecies

%   1, <R001>
i=i+1;
Rnames{   1} = 'A + HO = 0.5*B + 0.5*C ';
k(:,i) = (  2.0000E-12.*exp(  3.0000E+02./T) );
Gstr{i,   1}='A';Gstr{i,   2}='HO';
fA(i)=fA(i)-1.0;fHO(i)=fHO(i)-1.0;
fB(i)=fB(i)+  0.5000;fC(i)=fC(i)+  0.5000;

%   2, <R002>
i=i+1;
Rnames{   2} = 'B + HO = D ';
xko =   1.8000E-30.*M.*exp(  0.0000E+00./T).*(T./300).^ -3.0000E+00;
xkinf =   2.8000E-11.*exp(  0.0000E+00./T).*(T./300).^  0.0000E+00;
xn =   1.0000E+00;
F =   6.0000E-01;
G=1.0./(1.0+(log10(xko./xkinf)./xn).^2);
k(:,i) = (xko./( 1.0+xko./xkinf).*F.^G );
Gstr{i,   1}='B';Gstr{i,   2}='HO';
fB(i)=fB(i)-1.0;fHO(i)=fHO(i)-1.0;
fD(i)=fD(i)+  1.0000;

%   3, <R003>
i=i+1;
Rnames{   3} = 'D = B + HO ';
k(:,i) = (  1.0E+16.*exp( -1.0000E+04./T).*k(:,2) );
Gstr{i,   1}='D';
fD(i)=fD(i)-1.0;
fB(i)=fB(i)+  1.0000;fHO(i)=fHO(i)+  1.0000;

%   4, <R004>
i=i+1;
Rnames{   4} = 'D = B + HO ';
k(:,i) = (JD_TEST );
Gstr{i,   1}='D';
fD(i)=fD(i)-1.0;
fB(i)=fB(i)+  1.0000;fHO(i)=fHO(i)+  1.0000;

%   5, <R005>
i=i+1;
Rnames{   5} = 'C + HO + H2O = HO + A ';
k(:,i) = (  1.0E-30 ).*H2O;
Gstr{i,   1}='C';Gstr{i,   2}='HO';
fC(i)=fC(i)-1.0;fHO(i)=fHO(i)-1.0;
fHO(i)=fHO(i)+  1.0000;fA(i)=fA(i)+  1.0000;

%   6, <R006>
i=i+1;
Rnames{   6} = 'A =';
SEAWATER = 0.0;
a = 1.0E-5;
a(a>2.0) = 2.0;
k(:,i) = a.*SEAWATER;
Gstr{i,   1}='A';
fA(i)=fA(i)-1.0;
"""


def test_matlab_expressions_become_easy():
    assert fi.matlab_to_easy("(  6.1E-34.*(T./300).^( -2.4) ).*O2.*M;") == "(6.1E-34*(T/300)@(-2.4))*O2*M"
    assert fi.matlab_to_easy("1.8E-30.*M.*(T./300).^ -3.0E+00") == "1.8E-30*M*(T/300)@(-3.0E+00)"
    assert fi.matlab_to_easy("2.3E-12.*exp( -1.0E+3./T)") == "2.3E-12*EXP(-1.0E+3/T)"
    assert "LOG10(" in fi.matlab_to_easy("log10(x./y)")


def test_reading_handles_every_construct(tmp_path):
    p = tmp_path / "synth_AllRxns.m"
    p.write_text(SYNTH)
    species, rx, dropped = fi.read_allrxns(str(p), fixed={"SEAWATER": 0.0})
    by = {fi.easy_reaction(r): r for r in rx}
    assert species == ["A", "B", "C", "HO", "D"]
    assert "A + HO --> 0.5*B + 0.5*C" in by                              # fractional yields
    troe = by["B + HO --> D"].rate
    assert "xko" not in troe and "LOG10(" in troe and "@(-3.0000E+00)" in troe   # helpers substituted
    back = [r.rate for r in rx if fi.easy_reaction(r) == "D --> B + HO"]
    assert len(back) == 2 and any(troe in b and b.startswith("(1.0E+16") for b in back)   # k(:,2) = reaction 2's rate
    assert by["C + HO --> HO + A"].rate.endswith("*H2O")                  # catalytic HO stays a product
    assert len(dropped) == 1 and "rate zero" in dropped[0]                # the sea-water reaction
    merged = {}
    for r in rx:                                                           # read() merges duplicates
        merged.setdefault(fi.easy_reaction(r), []).append(r)
    assert len(merged["D --> B + HO"]) == 2


def test_wrap_never_breaks_a_number():
    line = "k[X --> Y] = CONST(" + "*".join(["1.4000E+04*EXP(-1.0E+03/T)"] * 12) + ")"
    out = fi.wrap(line, 80)
    assert all(l.endswith(" $") for l in out.split("\n")[:-1])
    assert out.replace(" $\n     ", "") == line


@pytest.fixture(scope="module")
def cracmm2():
    if "CRACMM2" not in mreg.load_bases():
        pytest.skip("CRACMM2 not installed")
    return mreg.base("CRACMM2")


def test_installed_cracmm2(cracmm2):
    d = cracmm2.lumped
    assert cracmm2.is_lumped and cracmm2.is_complete and not cracmm2.is_gecko
    assert d["n_species"] == 193 and d["n_reactions"] == 526
    assert 21 in d["j_f0am"] and 22 in d["j_f0am"]
    assert not cracmm2.verify()
    phot = open(cracmm2.path("phot")).read()
    assert "JHNO4_RACM2 = CONST(JN21+JN22)" in phot                       # F0AM's J4 for HNO4 corrected
    defined = set(re.findall(r"^(\w+)\s*=", phot, re.M)) | set(re.findall(r"^(\w+)\s*=", open(cracmm2.path("chemistry")).read(), re.M))
    text = open(cracmm2.path("chemistry")).read().replace(" $\n", "")
    used = set()
    for line in text.splitlines():
        if line.startswith("k["):
            used |= fi._identifiers(line.split("= CONST(", 1)[1])
    assert used - defined - {"T", "M", "O2", "N2", "H2", "H2O"} == set()
    assert "ISO" in mreg.species_names("CRACMM2") and "C5H8" in mreg.species_names("CRACMM2")


def _spec(chamber="SAPHIR"):
    exp = juice.new_experiment("cracmm2_test", chamber=chamber)
    c = exp.conditions
    c.source, c.date, c.start, c.end = "manual", "2024-07-09", "08:00", "10:00"
    c.temperature_mode, c.j_source, c.roof_open = "const", "clear_sky", [("08:00", "10:00")]
    exp.base_mechanism, exp.precursors = "CRACMM2", []
    exp.injections = [juice.Injection(name="NO", products="NO", total_ppb="5", start="08:05"),
                      juice.Injection(name="isoprene", products="C5H8", total_ppb="10", start="08:10")]
    return exp


def test_experiment_easy_under_the_sun(cracmm2):
    exp = _spec()
    assert not [m for lvl, m in juice.check(exp) if lvl != "info"]
    text = es.render_easy(exp)
    for fn in ("cracmm2.physics.incl", "cracmm2.aliases.incl", "cracmm2.chem.incl", "cracmm2.phot.incl", "CRACMM2.dil.mech"):
        assert f"#include {fn}" in text
    assert "anorg.mcm331" not in text and "saphir.jmeasured" not in text
    assert re.search(r"^J4 = PHTFKT\(1\.1650e-02,0\.244,-0\.267\)", text, re.M)
    assert re.search(r"^JN21 = PHTFKT\(", text, re.M)


def test_experiment_easy_under_lamps(cracmm2):
    exp = _spec("generic_lamps")
    text = es.render_easy(exp)
    lines = text.splitlines()
    assert not any(l.startswith("k[O3 + HV") or l.startswith("k[NO2 + HV") for l in lines)   # no MCM inorganic k[]
    j4 = next(i for i, l in enumerate(lines) if l.startswith("J4 = CONST("))
    jn = next(i for i, l in enumerate(lines) if l.startswith("JN21 = CONST(J3*"))
    assert j4 < jn


def test_cracmm2_runs_and_mcm_names_reach_its_species(cracmm2, tmp_path):
    run = juice.run(_spec(), str(tmp_path / "run"), solver="fast", kv_pd=False, script=False)
    r = run.results()
    assert "ISO" in r.species and "OH" in r.species and "HO" not in r.species
    iso = r.ppb("ISO")
    assert iso.max() > 5.0 and iso[-1] < 0.5 * iso.max()                  # injected as C5H8, oxidised as ISO
    assert r.ppb("O3")[-1] > 1.0
