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


# ── the NOAA biomass-burning add-on (F0AM Chem/MCMv331/NOAA_BiomassBurning) ─────────────────────────────────────
# Built as the difference of F0AM's NOAA MCM to F0AM's own MCM v3.3.1 (78 reactions, as in the maintained .m);
# of NOAA's five rate changes only the two BENZAL photolyses differ from JUICE's MCM download (the CISOPCO2 + RO2
# rates NOAA corrects are already 2.0e-12 in the MCM v3.3.1 download).

@pytest.fixture(scope="module")
def noaa_bb():
    if "NOAA_BB" not in mreg.load_addons():
        pytest.skip("NOAA_BB not installed")
    return mreg.load_addons()["NOAA_BB"]


def test_noaa_bb_addon_contents(noaa_bb):
    text = open(noaa_bb.paths()[0], encoding="latin-1").read()
    mech = [l for l in text.splitlines() if "-->" in l and not l.startswith(("k[", ";"))]
    rates = [l for l in text.splitlines() if l.startswith("k[")]
    assert len(mech) == 78 and len(rates) == 80                            # 78 new + 2 re-rated
    assert noaa_bb.overrides == ["BENZAL + hv  -->  C6H5CO3 + HO2", "BENZAL + hv  -->  C6H5O2 + CO + HO2"]
    assert "k[BENZAL+hv-->C6H5CO3+HO2]=CONST(JN2*0.5)" in text
    furf = [l for l in mech if l.startswith("FURFURAL") and "OH" not in l]
    assert any("FURFC3H4" in l for l in furf) and not any(re.search(r"\bC3H4\b", l) for l in furf)   # the .m product, not allene
    assert noaa_bb.j_f0am == [2, 52] and re.search(r"^JN52 = PHTFKT\(", text, re.M)
    assert "HYDFURANO2" in noaa_bb.ro2 and len(noaa_bb.ro2) == 12
    assert {"FURAN", "FURFURAL", "GUAIACOL", "SYRINGOL"} <= set(noaa_bb.precursors)
    assert noaa_bb.adds_species["GUAIACOL"] == "C7H8O2"


def _bb_spec(chamber="SAPHIR"):
    exp = juice.new_experiment("bb_test", chamber=chamber)
    c = exp.conditions
    c.source, c.date, c.start, c.end = "manual", "2024-07-09", "08:00", "10:00"
    c.temperature_mode, c.j_source, c.roof_open = "const", "clear_sky", [("08:10", "10:00")]
    exp.precursors, exp.addons = ["TOLUENE"], ["NOAA_BB"]
    exp.injections = [juice.Injection(name="NO", products="NO", total_ppb="5", start="08:05"),
                      juice.Injection(name="furan", products="FURAN", total_ppb="10", start="08:15"),
                      juice.Injection(name="guaiacol", products="GUAIACOL", total_ppb="5", start="08:15")]
    return exp


def test_noaa_bb_under_lamps_redefines_its_channels(noaa_bb):
    text = es.render_easy(_bb_spec("generic_lamps"))
    lines = text.splitlines()
    assert any(l.startswith("JN52 = CONST(J18*") for l in lines)


def test_noaa_bb_runs_with_the_mcm(noaa_bb, tmp_path):
    exp = _bb_spec()
    assert not [m for lvl, m in juice.check(exp) if lvl == "error"]
    run = juice.run(exp, str(tmp_path / "bb"), solver="fast", kv_pd=False, script=False)
    r = run.results()
    assert r.ppb("FURAN").max() > 5 and r.ppb("FURAN")[-1] < 0.5 * r.ppb("FURAN").max()
    assert r.ppb("NGUAIACOL")[-1] > 0.01
    chem = os.path.join(run.out_dir, "chemistry")
    ro2 = next(f for f in os.listdir(chem) if f.endswith(".ro2.incl"))
    assert "HYDFURANO2" in open(os.path.join(chem, ro2)).read()
    elog = open(os.path.join(run.out_dir, "bb_test_juice.canonical.elog"), encoding="utf-8", errors="replace").read()
    assert "k[BENZAL+hv-->C6H5CO3+HO2]=CONST(JN2*0.5)" in elog


def test_installed_cracmm1():
    if "CRACMM1_aq" not in mreg.load_bases():
        pytest.skip("CRACMM1_aq not installed")
    b = mreg.base("CRACMM1_aq")
    assert b.is_lumped and b.lumped["n_species"] == 178 and b.lumped["n_reactions"] == 501
    assert "CH4=CONST(1.85E-6*M)" in open(b.path("physics")).read()          # held at 1850 ppb, as in CMAQ
    assert "JHNO4_RACM2 = CONST(JN21+JN22)" in open(b.path("phot")).read()
    assert any("OPEN_OCEAN" in d for d in b.lumped["dropped"])                 # F0AM's marine O3 sink: not in a chamber
    text = open(b.path("chemistry")).read().replace(" $\n", "")
    defined = set(re.findall(r"^(\w+)\s*=", open(b.path("phot")).read() + "\n" + text, re.M))
    used = set()
    for line in text.splitlines():
        if line.startswith("k["):
            used |= fi._identifiers(line.split("= CONST(", 1)[1])
    assert used - defined - {"T", "M", "O2", "N2", "H2", "H2O", "CH4"} == set()


def test_installed_racm2():
    if "RACM2" not in mreg.load_bases():
        pytest.skip("RACM2 not installed")
    b = mreg.base("RACM2")
    assert b.is_lumped and b.lumped["n_species"] == 116 and b.lumped["n_reactions"] == 360
    chem = open(b.path("chemistry")).read()
    assert "RCO3 + MO2 --> 0.9*MO2 + 0.9*HO2 + 0.4*CO2 + HCHO + 0.1*ORA2" in chem      # F0AM's R268 consumed ACO3
    assert "INCONSISTENT" not in chem
    assert "H2=CONST" not in open(b.path("physics")).read()                          # H2 is a RACM2 species
    assert "JHO2NO2 = CONST(JN21+JN22)" in open(b.path("phot")).read()
    k_lines = re.findall(r"^(K_\w+) = CONST\(", chem, re.M)
    assert len(k_lines) == 16                                                      # the K file's fall-offs, written out


def test_comment_lines_of_a_block_are_ignored(tmp_path):
    src = SYNTH.replace("fB(i)=fB(i)-1.0;fHO(i)=fHO(i)-1.0;\nfD(i)=fD(i)+  1.0000;",
                        "fB(i)=fB(i)-1.0;fHO(i)=fHO(i)-1.0;\n% Gstr{i,3}='C'; fC(i)=fC(i)-1.0;\nfD(i)=fD(i)+  1.0000;")
    assert src != SYNTH
    p = tmp_path / "c_AllRxns.m"
    p.write_text(src)
    _sp, rx, _dr = fi.read_allrxns(str(p), fixed={"SEAWATER": 0.0})
    assert "B + HO --> D" in {fi.easy_reaction(r) for r in rx}


def test_arithmetic_stoichiometry():
    assert fi._arith("-1+0.9") == pytest.approx(-0.1) and fi._arith("+0.5*2") == 1.0
    with pytest.raises(ValueError):
        fi._arith("+x")
