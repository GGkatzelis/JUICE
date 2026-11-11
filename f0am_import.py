"""
f0am_import.py — a mechanism of the F0AM chemistry library as a JUICE base mechanism
====================================================================================

F0AM (Wolfe et al., 2016; github.com/AirChem/F0AM, GPL-3.0 like JUICE) ships its mechanisms as three
MATLAB files with the same layout for every mechanism:

    <mech>_AllRxns.m   one block per reaction:  Rnames{i} = 'A + B = 0.5*C + D';  k(:,i) = <rate>;
                       Gstr{i,1} = 'A'; ...  fA(i)=fA(i)-1; fC(i)=fC(i)+0.5; ...   (local helper variables
                       such as xko/xkinf/F/G of a fall-off are assigned inside the block before k)
    <mech>_K.m         named rate constants the reactions refer to (heterogeneous: 1e-60 in F0AM)
    <mech>_J.m         named photolysis frequencies as MCM J numbers (Jmcm.J4) or F0AM's own channels
                       without an MCM number (Jmcm.Jn14)

This module turns them into EASY include files, so the mechanism runs in JUICE and in FACSIMILE alike:

    <id>.physics.incl  M, O2, N2, H2 (F0AM's constant species; H2O comes from the experiment)
    <id>.aliases.incl  ALIAS lines: the standard SAPHIR names (OH = HO, O = O3P, ...) and the MCM names
                       of the mechanism's explicit species (C5H8 -> ISO, ...), so templates and chamber
                       sources written for the MCM work unchanged
    <id>.chem.incl     every reaction and its rate, written out
    <id>.phot.incl     the mechanism's J names in terms of J<n> (MCM numbers) and JN<k> (F0AM channels);
                       the experiment defines J<n>/JN<k> for its light (sun, measured, lamps)
    <id>.dil.mech/.rate  dilution of every species
    <id>.species.tsv   the species

Reactants are F0AM's Gstr (the species solved for); the rate already carries the constant species
(O2, M, H2O, H2) as factors, exactly as in F0AM.  Products are the net stoichiometry of the f-lines plus
the reactants, so catalytic species and fractional yields come out right.  Duplicate reactions (the same
reactants and products) are merged into one with the summed rate.  Expressions are translated from MATLAB
(.* ./ .^ exp log10) to EASY (* / @ EXP LOG10) with the block's helper variables substituted.

    python f0am_import.py CRACMM2 <folder with cracmm2_AllRxns.m, cracmm2_K.m, cracmm2_J.m>
"""
from __future__ import annotations

import datetime as _dt
import hashlib
import json
import os
import re
import sys
from collections import OrderedDict
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

HERE = os.path.dirname(os.path.abspath(__file__))
BASE_DIR = os.path.join(HERE, "mechanisms", "base")

# What a mechanism needs besides its files: the MCM names of its explicit species (canonical name =
# the mechanism's own), corrections to the F0AM J map, and the citation.
PROFILES: Dict[str, dict] = {
    "CRACMM2": {
        "name": "CRACMM 2",
        "files": ("cracmm2_AllRxns.m", "cracmm2_K.m", "cracmm2_J.m"),
        "citation": "Skipper et al. (2024); F0AM version of the CMAQ mech.def by N. Skipper (EPA)",
        "f0am_path": "Chem/CRACMM2",
        "family": "lumped",
        # MCM name -> mechanism species (only where the species is one compound, not a lumped class)
        "mcm_names": {"C5H8": "ISO", "TOLUENE": "TOL", "BENZENE": "BEN", "C2H5OH": "EOH", "CH3OH": "MOH",
                      "C2H6": "ETH", "C2H4": "ETE", "APINENE": "API", "LIMONENE": "LIM", "C2H2": "ACE",
                      "CH3CHO": "ACD", "CH3COCH3": "ACT", "HCOOH": "ORA1", "CH3CO2H": "ORA2", "GLYOX": "GLY",
                      "MGLYOX": "MGLY", "STYRENE": "STY", "C4H6": "BDE13", "ACR": "ACRO", "CH3O2": "MO2",
                      "CH3CO3": "ACO3", "CH3CO3H": "PAA", "C2H5O2": "ETHP"},
        # F0AM's J map has JHNO4_RACM2 = Jmcm.J4 (the NO2 photolysis frequency, ~1e3 too fast); HNO4 photolysis is
        # the sum of F0AM's own HO2NO2 channels Jn21 (HO2 + NO2) and Jn22 (HO + NO3).
        "j_fixes": {"JHNO4_RACM2": ("JN21+JN22", "F0AM maps JHNO4_RACM2 to J4 (jNO2); HO2NO2 photolysis is Jn21 + Jn22")},
        # CMAQ-only switches inside rate blocks (no ocean in a chamber or box): value used
        "fixed": {"SEAWATER": 0.0, "OPEN_OCEAN": 0.0, "SURF_ZONE": 0.0, "SEA_ICE": 0.0, "ILLUMINATED": 1.0},
    },
}
PROFILES["CRACMM1_aq"] = {
    **PROFILES["CRACMM2"],
    "name": "CRACMM 1 (with aqueous chemistry)",
    "files": ("CRACMM1_aq_AllRxns.m", "CRACMM1_aq_K.m", "CRACMM1_aq_J.m"),
    "citation": "Place et al. (2023); F0AM version of the CMAQ mech.def by B. Place (EPA)",
    "f0am_path": "Chem/CRACMM1",
    # CMAQ (and F0AM's CRACMM1) hold CH4 at 1850 ppb and H2 at 560 ppb; CH4 is not a species of the mechanism
    "constants": {"CH4": "1.85E-6*M"},
    "j_fixes": {"JHNO4_RACM2": ("JN21+JN22", "F0AM maps JHNO4_RACM2 to J3 (H2O2); HO2NO2 photolysis is Jn21 + Jn22, as for CRACMM2")},
}
PROFILES["RACM2"] = {
    "name": "RACM2",
    "files": ("RACM2_AllRxns.m", "RACM2_K.m", "RACM2_J.m"),
    "citation": "Goliff et al. (2013), Atmos. Environ.; F0AM version (S. R., M. M., G. M. Wolfe, 2015-2016, with the corrections listed in its header)",
    "f0am_path": "Chem/RACM2",
    "family": "lumped",
    "mcm_names": PROFILES["CRACMM2"]["mcm_names"],          # only the names of species RACM2 has are written
    "j_fixes": {"JHO2NO2": ("JN21+JN22", "F0AM maps JHO2NO2 to J3 (H2O2); HO2NO2 photolysis is Jn21 + Jn22, as for CRACMM")},
    # R268 (RCO3 + MO2) in F0AM's file consumes ACO3 instead of RCO3 (copied from R267): its reactants are RCO3 + MO2
    "block_fixes": {"RCO3+MO2=0.9HO2+0.9MO2+0.4CO2+HCHO+0.1ORA2": [("fACO3(i)=fACO3(i)-1; fMO2", "fRCO3(i)=fRCO3(i)-1; fMO2")]},
    "fixed": {},
}

# the standard SAPHIR aliases (as geckoa.aliases.incl): canonical = what the SAPHIR includes and the outputs use
STANDARD_ALIASES = [("HV", "hv"), ("OH", "HO"), ("HCHO", "H2CO CH2O"), ("O", "O3P"), ("CH3OOH", "CH3O2H MOOH"),
                    ("HONO", "HNO2"), ("HNO3", "HONO2"), ("HNO4", "HO2NO2"), ("SA", "SULF H2SO4")]
CONSTANT_SPECIES = ("M", "O2", "N2", "H2", "H2O")
# names a lumped mechanism shares with the MCM for the SAME compound (inorganics and a few explicit organics):
# their formula is taken from the MCM species table; every other species of a lumped mechanism has none
SAME_AS_MCM = {"O3", "O", "O1D", "OH", "HO2", "H2O2", "NO", "NO2", "NO3", "N2O5", "HONO", "HNO3", "HNO4", "HO2NO2",
               "CO", "CO2", "SO2", "H2", "CH4", "HCHO", "PAN", "MVK", "MACR", "MEK", "MPAN"}
EASY_VARIABLES = {"T", "M", "O2", "N2", "H2", "H2O", "P"}


@dataclass
class Reaction:
    label: str                        # F0AM's <R001> label
    reactants: List[str]
    products: "OrderedDict[str, float]"
    rate: str                         # EASY expression
    notes: List[str] = field(default_factory=list)


@dataclass
class Mechanism:
    species: List[str]
    reactions: List[Reaction]
    k_names: Dict[str, str]           # K file: name -> EASY value
    j_map: Dict[str, str]             # J name -> EASY expression in J<n> / JN<k>
    dropped: List[str] = field(default_factory=list)
    notes: List[str] = field(default_factory=list)


# ── MATLAB -> EASY ───────────────────────────────────────────────────────────────────────────────
_NUM = r"[0-9]*\.?[0-9]+(?:[EeDd][-+]?[0-9]+)?"


def matlab_to_easy(expr: str) -> str:
    """One MATLAB rate expression in EASY notation (no helper variables left)."""
    e = expr.strip().rstrip(";").strip()
    e = re.sub(r"\.\^\s*(-\s*" + _NUM + r")", r"@(\1)", e)            # x.^ -2.6  ->  x@(-2.6)
    e = re.sub(r"\^\s*(-\s*" + _NUM + r")", r"@(\1)", e)
    e = e.replace(".*", "*").replace("./", "/").replace(".^", "@").replace("^", "@")
    e = re.sub(r"\bexp\(", "EXP(", e)
    e = re.sub(r"\blog10\(", "LOG10(", e)
    e = re.sub(r"\bsqrt\(", "SQRT(", e)
    e = re.sub(r"\s+", "", e)
    return e


def _identifiers(e: str) -> set:
    """Names in an EASY expression (variables, species, functions other than EXP/LOG10/SQRT)."""
    return set(re.findall(r"(?<![0-9.A-Za-z_])([A-Za-z_][A-Za-z0-9_]*)", e)) - {"EXP", "LOG10", "SQRT"}


def _substitute(expr: str, local: Dict[str, str]) -> str:
    """Replace helper variables by their (already substituted) definitions, in parentheses."""
    def rep(m):
        n = m.group(1)
        return f"({local[n]})" if n in local else n
    prev = None
    while prev != expr:
        prev = expr
        expr = re.sub(r"(?<![0-9.A-Za-z_])([A-Za-z_][A-Za-z0-9_]*)(?![A-Za-z0-9_])", rep, expr)
    return expr


def _num(text: str) -> float:
    return float(text.replace("D", "E").replace("d", "e"))


# ── reading the three files ─────────────────────────────────────────────────────────────────────
_ARITH = re.compile(r"^[0-9.EeDd+\-*/ ()]+$")


def _arith(text: str) -> float:
    """A stoichiometric increment written as arithmetic ('1+0.9', '-1', '0.5*2'): numbers and + - * / only."""
    t = text.strip().replace("D", "E").replace("d", "e")
    if not _ARITH.match(t):
        raise ValueError(f"not a number: {text!r}")
    return float(eval(t, {"__builtins__": {}}, {}))


def _strip_comments(blk: str) -> str:
    return "\n".join(l.split("%", 1)[0] for l in blk.split("\n"))


def read_allrxns(path: str, fixed: Optional[Dict[str, float]] = None,
                 block_fixes: Optional[Dict[str, List[Tuple[str, str]]]] = None) -> Tuple[List[str], List[Reaction], List[str]]:
    """(species, reactions, dropped) of an F0AM <mech>_AllRxns.m.  Comment text is ignored (F0AM files keep
    superseded stoichiometry as % lines).  block_fixes: {Rnames text: [(old, new), ...]} replacements in that block,
    for errors of the source file (each one is listed in the reaction's notes)."""
    fixed = dict(fixed or {})
    block_fixes = dict(block_fixes or {})
    text = open(path, encoding="latin-1").read()
    m = re.search(r"SpeciesToAdd\s*=\s*\{(.*?)\};", text, re.S)
    species = re.findall(r"'([^']+)'", m.group(1)) if m else []
    blocks = re.split(r"\n\s*i\s*=\s*i\s*\+\s*1\s*;", text)[1:]
    reactions: List[Reaction] = []
    dropped: List[str] = []
    rate_by_index: Dict[int, str] = {}                 # F0AM reaction number -> EASY rate (k(:,n) refers to it)
    for index, raw_blk in enumerate(blocks, start=1):
        label_m = re.search(r"Rnames\{\s*([^}]*?)\s*\}\s*=\s*'([^']*)'", raw_blk)
        name = label_m.group(2).strip() if label_m else "?"
        if label_m and label_m.group(1).isdigit():
            index = int(label_m.group(1))
        fixes_here = block_fixes.get(name, [])
        for old, new in fixes_here:
            if old not in raw_blk:
                raise ValueError(f"block fix for '{name}' does not apply: {old!r} not found")
            raw_blk = raw_blk.replace(old, new)
        tag = re.match(r"\s*%\s*(R\w+)", raw_blk)          # RACM2 style: i=i+1; %R266
        if tag:
            name = f"{tag.group(1)}: {name}"
        blk = _strip_comments(raw_blk)
        local: Dict[str, str] = {}
        rate = None
        notes: List[str] = []
        unsupported = []
        for raw in blk.split("\n"):
            line = raw.split("%", 1)[0].strip()
            if not line or line.startswith(("Rnames", "Gstr")) or re.match(r"f[A-Za-z0-9_]+\(i\)", line):
                continue
            mk = re.match(r"k\(:,i\)\s*=\s*(.+?);?\s*$", line)
            if mk:
                rate = mk.group(1); continue
            ma = re.match(r"([A-Za-z_][A-Za-z0-9_]*)\s*=\s*(.+?);\s*$", line)
            if ma:
                var, val = ma.group(1), ma.group(2)
                if var in fixed:
                    continue
                if re.search(r"[<>]|==|~=", val):                         # a logical mask (SZA > 0)
                    unsupported.append(line); continue
                local[var] = _substitute(matlab_to_easy(val), local)
                continue
            if re.match(r"[A-Za-z_][A-Za-z0-9_]*\(.*\)\s*=", line):           # masked assignment a(a>b) = b
                unsupported.append(line); continue
        if rate is None:
            continue
        # constants of the profile (CMAQ ocean switches) first, then helpers
        consts = {k: repr(float(v)) for k, v in fixed.items()}
        rate = re.sub(r"k\(:,\s*(\d+)\)", lambda m_: "(" + rate_by_index[int(m_.group(1))] + ")", rate)
        e = _substitute(matlab_to_easy(rate), {**consts, **local})
        rate_by_index[index] = e
        reac = re.findall(r"Gstr\{i,\s*\d+\}\s*=\s*'([^']+)'", blk)
        net: "OrderedDict[str, float]" = OrderedDict()
        for sp_, val in re.findall(r"f([A-Za-z0-9_]+)\(i\)\s*=\s*f\1\(i\)\s*([-+][0-9.EeDd+\-*/ ]*?)\s*(?:;|$)", blk, re.M):
            net[sp_] = net.get(sp_, 0.0) + _arith(val)
        notes += [f"source corrected: {old.strip()} -> {new.strip()}" for old, new in fixes_here]
        prods: "OrderedDict[str, float]" = OrderedDict()
        for sp_, v in net.items():
            c = v + reac.count(sp_)
            if abs(c) > 1e-12:
                prods[sp_] = round(c, 6)
        for sp_ in reac:
            if net.get(sp_, 0.0) + reac.count(sp_) < -1e-12:
                notes.append(f"{sp_}: consumed more than the reactants say")
        bad = [s for s, c in prods.items() if c < 0 and not s.startswith("DELTA_")]   # CMAQ's DELTA_C/N: bookkeeping, negative by design
        if bad:                                        # the f-lines consume a species that is not a reactant
            notes.append("INCONSISTENT: consumes " + ", ".join(bad) + " which the reactants (Gstr) do not list")
        if _is_zero(e):
            dropped.append(f"{name}: rate zero in a chamber or box ({rate.strip()})"); continue
        if unsupported:
            dropped.append(f"{name}: not translatable ({'; '.join(unsupported)})"); continue
        reactions.append(Reaction(label=name, reactants=reac, products=prods, rate=e, notes=notes))
    return species, reactions, dropped


def _is_zero(e: str) -> bool:
    """A rate that is zero: a literal 0, or a product with a zero factor such as a*1.0*0.0 (the CMAQ sea-water switch)."""
    return bool(re.fullmatch(r"\(?0(\.0*)?\)?", e) or re.search(r"(^|\*)\(?0\.0\)?(\*|$)", e))


def read_k(path: str) -> Dict[str, str]:
    """Named rate constants of <mech>_K.m: one block per constant (Knames{i} = 'X'; helper assignments;
    krx(:,i) = expression), translated like a reaction block."""
    out: Dict[str, str] = {}
    if not os.path.isfile(path):
        return out
    text = open(path, encoding="latin-1").read()
    for blk in re.split(r"\n\s*i\s*=\s*i\s*\+\s*1\s*;", text)[1:]:
        blk = _strip_comments(blk)
        mname = re.search(r"Knames\{i\}\s*=\s*'([^']+)'", blk)
        if not mname:
            continue
        local: Dict[str, str] = {}
        value = None
        for line in blk.split("\n"):
            for stmt in (s.strip() for s in line.split(";")):
                if not stmt or stmt.startswith("Knames"):
                    continue
                mk = re.match(r"krx\(:,i\)\s*=\s*(.+)$", stmt)
                if mk:
                    value = _substitute(matlab_to_easy(mk.group(1)), local); continue
                ma = re.match(r"([A-Za-z_][A-Za-z0-9_]*)\s*=\s*(.+)$", stmt)
                if ma:
                    local[ma.group(1)] = _substitute(matlab_to_easy(ma.group(2)), local)
        if value is not None:
            out[mname.group(1)] = value
    return out


def read_j(path: str) -> Dict[str, str]:
    """J name -> EASY expression in J<n> / JN<k> of <mech>_J.m."""
    out: Dict[str, str] = {}
    for name, val in re.findall(r"^\s*J\.([A-Za-z0-9_]+)\s*=\s*([^;%]+);", open(path, encoding="latin-1").read(), re.M):
        e = re.sub(r"Jmcm\.Jn(\d+)", r"JN\1", val)
        e = re.sub(r"Jmcm\.J(\d+)", r"J\1", e)
        out[name] = re.sub(r"\s+", "", e)
    return out


def read(folder: str, profile: dict) -> Mechanism:
    rx, k, j = (os.path.join(folder, f) for f in profile["files"])
    species, reactions, dropped = read_allrxns(rx, profile.get("fixed"), profile.get("block_fixes"))
    jm = read_j(j)
    notes = []
    for name, (expr, why) in (profile.get("j_fixes") or {}).items():
        if name in jm:
            notes.append(f"{name} = {expr} instead of F0AM's {jm[name]}: {why}")
            jm[name] = expr
    # merge duplicates (same reactants and products): summed rate
    merged: "OrderedDict[str, Reaction]" = OrderedDict()
    for r in reactions:
        key = easy_reaction(r)
        if key in merged:
            merged[key].rate = f"{merged[key].rate}+{r.rate}"
            merged[key].notes.append(f"merged with {r.label}")
        else:
            merged[key] = r
    return Mechanism(species=species, reactions=list(merged.values()), k_names=read_k(k), j_map=jm,
                     dropped=dropped, notes=notes)


# ── writing ──────────────────────────────────────────────────────────────────────────────────────
def _coef(c: float) -> str:
    return "" if abs(c - 1.0) < 1e-12 else f"{c:g}*"


def wrap(line: str, width: int = 110) -> str:
    """A long statement over several lines with FACSIMILE's $ continuation, broken only after '*' or ')+'
    (never inside a number such as 1.4E+04)."""
    if len(line) <= width:
        return line
    out, cur = [], ""
    i = 0
    while i < len(line):
        cur += line[i]
        brk = line[i] == "*" or (line[i] == "+" and i > 0 and line[i - 1] == ")")
        if brk and len(cur) >= width - 10 and i < len(line) - 1:
            out.append(cur + " $"); cur = "     "
        i += 1
    out.append(cur)
    return chr(10).join(out)


def easy_reaction(r: Reaction) -> str:
    lhs = " + ".join(r.reactants)
    rhs = " + ".join(f"{_coef(c)}{s}" for s, c in r.products.items())
    return f"{lhs} --> {rhs}".rstrip()


def mcm_j_numbers(mech: Mechanism) -> Tuple[List[int], List[int]]:
    """(MCM J numbers, F0AM Jn numbers) the J map refers to."""
    used = " ".join(mech.j_map[n] for n in mech.j_map)
    return (sorted({int(x) for x in re.findall(r"(?<![A-Za-z0-9_])J(\d+)", used)}),
            sorted({int(x) for x in re.findall(r"\bJN(\d+)", used)}))


def write_base(base_id: str, folder: str, out_root: str = BASE_DIR, profile: Optional[dict] = None) -> str:
    """Convert the F0AM files in `folder` into mechanisms/base/<base_id>/ (overwritten).  Returns the folder."""
    prof = profile or PROFILES[base_id]
    mech = read(folder, prof)
    out = os.path.join(out_root, base_id)
    os.makedirs(out, exist_ok=True)
    stem = base_id.lower()
    files = {"physics": f"{stem}.physics.incl", "aliases": f"{stem}.aliases.incl", "chemistry": f"{stem}.chem.incl",
             "phot": f"{stem}.phot.incl", "dil_mech": f"{base_id}.dil.mech", "dil_rate": f"{base_id}.dil.rate",
             "species_tsv": f"{stem}.species.tsv"}
    src = {f: hashlib.sha256(open(os.path.join(folder, f), "rb").read()).hexdigest() for f in prof["files"]}
    head = [f"; {{file}} -- {prof['name']} from the F0AM chemistry library ({prof['f0am_path']}), converted by f0am_import.py "
            f"{_dt.date.today().isoformat()}",
            f"; {prof['citation']}.  F0AM: Wolfe et al. (2016), github.com/AirChem/F0AM, GPL-3.0.",
            "; source files sha256[:16]: " + ", ".join(f"{k} {v[:16]}" for k, v in src.items())]

    def header(name):
        return [h.replace("{file}", name) for h in head]

    used_species = []
    for r in mech.reactions:
        for s in list(r.reactants) + list(r.products):
            if s not in used_species:
                used_species.append(s)
    species = [s for s in mech.species if s in used_species] + [s for s in used_species if s not in mech.species]

    with open(os.path.join(out, files["physics"]), "w", encoding="ascii", errors="replace", newline="\n") as fh:
        fh.write("\n".join(header(files["physics"]) + [
            "; physics: air density from the pressure and temperature settings, F0AM's constant species",
            "CONST=FUNCTION[%1]", "DENSITY=FUNCTION[(%1)/((%2)*1.379E-19)]", "M =DENSITY(P,T)",
            "O2=CONST(0.2095*M)", "N2=CONST(0.7808*M)"]
            # H2 is a constant (560 ppb) where the mechanism has no H2 species (CRACMM); RACM2 solves for H2
            + ([] if "H2" in species else ["H2=CONST(5.6E-7*M)"])
            + [f"{n}=CONST({v})   ; held constant, as in F0AM and CMAQ" for n, v in prof.get("constants", {}).items()
               if n not in species] + [""]))

    with open(os.path.join(out, files["aliases"]), "w", encoding="ascii", newline="\n") as fh:
        L = header(files["aliases"]) + [
            "; ALIAS[canonical] = alternatives: the parsers replace every alternative by the canonical name in reactions",
            ";   and k[] keys. First the standard SAPHIR names, then the MCM names of this mechanism's explicit species",
            ";   (the canonical name is the mechanism's own, which the output carries)."]
        L += [f"ALIAS[{c}] = {alt}" for c, alt in STANDARD_ALIASES]
        taken = {c for c, _ in STANDARD_ALIASES} | {a for _, alts in STANDARD_ALIASES for a in alts.split()}
        for mcm, own in prof.get("mcm_names", {}).items():
            if own in species and mcm not in taken and own not in taken:
                L.append(f"ALIAS[{own}] = {mcm}")
        fh.write("\n".join(L) + "\n")

    with open(os.path.join(out, files["chemistry"]), "w", encoding="ascii", errors="replace", newline="\n") as fh:
        L = header(files["chemistry"]) + [
            f"; {len(species)} species, {len(mech.reactions)} reactions; rates as in F0AM (constant species O2, M, H2O, H2",
            ";   are factors of the rate); duplicates merged into one reaction with the summed rate."]
        if mech.dropped:
            L.append("; not converted: " + " | ".join(mech.dropped))
        L.append("")
        if mech.k_names:
            L.append("; named rate constants (F0AM " + prof["files"][1] + "; heterogeneous uptake is off in F0AM: 1e-60)")
            for n, v in mech.k_names.items():
                L.append(wrap(f"{n} = CONST({v})"))
            L.append("")
        for r in mech.reactions:
            rx = easy_reaction(r)
            L.append(f"; {r.label}" + (f"  ({'; '.join(r.notes)})" if r.notes else ""))
            L.append(rx)
            L.append(wrap(f"k[{rx}] = CONST({r.rate})"))
        fh.write("\n".join(L) + "\n")

    jm, jn = mcm_j_numbers(mech)
    with open(os.path.join(out, files["phot"]), "w", encoding="ascii", newline="\n") as fh:
        L = header(files["phot"]) + [
            "; the mechanism's photolysis frequencies as F0AM maps them: J<n> = MCM v3.3.1 J number n, JN<k> = F0AM's channel",
            ";   Jn<k> (no MCM number). The experiment defines J<n> and JN<k> for its light (sun, measured J, lamps)."]
        L += [f"; {n}" for n in mech.notes]
        for name, expr in mech.j_map.items():
            L.append(f"{name} = CONST({expr})")
        fh.write("\n".join(L) + "\n")

    with open(os.path.join(out, files["dil_mech"]), "w", encoding="ascii", newline="\n") as fh:
        fh.write("\n".join([f"; {files['dil_mech']} -- dilution of every species (f0am_import.py)"] + [f"{s} + DIL --> " for s in species]) + "\n")
    with open(os.path.join(out, files["dil_rate"]), "w", encoding="ascii", newline="\n") as fh:
        fh.write("\n".join([f"; {files['dil_rate']} -- dilution of every species (f0am_import.py)"] + [f"k[{s} + DIL --> ] = CONST(DILUTE)" for s in species]) + "\n")
    with open(os.path.join(out, files["species_tsv"]), "w", encoding="utf-8", newline="\n") as fh:
        rev = {own: mcm for mcm, own in prof.get("mcm_names", {}).items()}
        alias_of = {a: c for c, alts in STANDARD_ALIASES for a in alts.split()}
        try:                                           # formulas where the species is ONE compound: its MCM name,
            import species_db as _sdb                  #   or an inorganic / unambiguous name shared with the MCM
            mcm_db = _sdb.load_species_db()
        except Exception:
            mcm_db = None

        def formula_of(s):
            if mcm_db is None:
                return ""
            name = rev.get(s) or alias_of.get(s) or (s if s in SAME_AS_MCM else "")
            return mcm_db.formula(name) if name else ""
        fh.write("Name\tFormula\tMCM name\tPeroxyRadical\n"
                 + "".join(f"{s}\t{formula_of(s)}\t{rev.get(s, '')}\tfalse\n" for s in species))

    import shutil                                      # the F0AM originals, for provenance (GPL-3.0, as JUICE)
    os.makedirs(os.path.join(out, "f0am_source"), exist_ok=True)
    for fn in prof["files"]:
        shutil.copy2(os.path.join(folder, fn), os.path.join(out, "f0am_source", fn))
    checks = {fn: hashlib.sha256(open(os.path.join(out, fn), "rb").read()).hexdigest() for fn in files.values()}
    manifest = {
        "kind": "lumped", "id": base_id, "name": prof["name"],
        "version": f"F0AM {prof['f0am_path']}, converted {_dt.date.today().isoformat()}",
        "files": files, "checksums": checks,
        "source": {"f0am_path": prof["f0am_path"], "url": "https://github.com/AirChem/F0AM/tree/master/" + prof["f0am_path"],
                   "copy": "f0am_source/", "files_sha256": src,
                   "converted": _dt.date.today().isoformat(), "converter": "f0am_import.py",
                   "license": "GPL-3.0 (F0AM)"},
        "citation": prof["citation"], "status": "draft",
        "notes": "; ".join(mech.notes + mech.dropped),
        "lumped": {"family": prof.get("family", "lumped"), "n_species": len(species), "n_reactions": len(mech.reactions),
                   "j_mcm": jm, "j_f0am": jn, "mcm_names": prof.get("mcm_names", {}), "dropped": mech.dropped},
    }
    with open(os.path.join(out, "manifest.json"), "w", encoding="utf-8", newline="\n") as fh:
        json.dump(manifest, fh, indent=1)
    return out



# == Add-ons: what an F0AM variant of the MCM changes, as an EASY include on top of JUICE's MCM ==============
# F0AM ships some MCM extensions as a complete modified MCM (FACSIMILE .fac and MATLAB .m).  The add-on is the
# difference to F0AM's own MCM v3.3.1: new reactions, and k[...] = CONST(...) overrides for re-rated base
# reactions (written on JUICE's reaction keys, after the base: last definition wins).  The reaction text comes
# from the .fac (MCM conventions: RO2 in the rate, J<n> labels) and is converted with mcm_converter, exactly as
# the base; corrections F0AM made later in the .m are applied by the profile and checked against the .m.
ADDON_PROFILES: Dict[str, dict] = {
    "NOAA_BB": {
        "name": "NOAA biomass burning: furans, furfurals, phenolics",
        "f0am_path": "Chem/MCMv331/NOAA_BiomassBurning",
        "citation": "Decker et al. (2019), ES&T; Coggon et al. (MCMUpdates_BiomassBurning_DeckerCoggon.pptx); F0AM version G. M. Wolfe",
        "f0am_base": ("MCMv331_AllRxns.fac", "MCMv331_AllRxns.m"),
        "f0am_mod": ("MCMv331_AllRxns_NOAABB.fac", "MCMv331_AllRxns_NOAABB.m"),
        "compatible_with": ["MCM_v3.3.1"],
        # the .fac predates two corrections of the .m (GMW 2019): BENZAL photolysis on F0AM's Jn2 (IUPAC cross
        # sections) instead of the placeholder J<59>; furfural photolysis on Jn52 with the product FURFC3H4
        # (C3H4 is allene in the MCM)
        "j_rename": {"J<59>": "JN2", "J<58>": "JN52"},
        "product_rename": {("FURFURAL", "C3H4"): "FURFC3H4"},
        "formulas": {"FURAN": "C4H4O", "M2F": "C5H6O", "M3F": "C5H6O", "DIMEFURAN": "C6H8O", "FURFURAL": "C5H4O2",
                     "MEFURFURAL": "C6H6O2", "GUAIACOL": "C7H8O2", "MGUAIACOL": "C8H10O2", "SYRINGOL": "C8H10O3",
                     "ACETONIT": "C2H3N"},
        "precursors": ["FURAN", "M2F", "M3F", "DIMEFURAN", "FURFURAL", "MEFURFURAL", "GUAIACOL", "MGUAIACOL", "SYRINGOL",
                       "DBUTANOL"],
    },
}


def _fac_reactions(path: str) -> List[Tuple[str, str]]:
    """(rate, reaction) of every `% rate : A + B = C ;` of a FACSIMILE file (continuations joined)."""
    from pathlib import Path
    import mcm_converter as mc
    return mc._extract_reactions(mc._combine_continuations(mc._read_lines(Path(path))))


def _fac_key(reaction: str) -> Tuple[tuple, tuple]:
    lhs, rhs = reaction.split("=", 1)

    def norm(side):
        return tuple(sorted(re.sub(r"\s+", "", t) for t in side.split("+") if t.strip()))
    return norm(lhs), norm(rhs)


def _multiset(rx: List[Tuple[str, str]]) -> Dict[tuple, List[str]]:
    out: Dict[tuple, List[str]] = {}
    for rate, r in rx:
        if "=" in r:
            out.setdefault(_fac_key(r), []).append(re.sub(r"\s+", "", rate))
    return out


def _species_of(reactions: List[str], sep: str) -> List[str]:
    out: List[str] = []
    for r in reactions:
        for t in re.split(r"\+|" + re.escape(sep), r):
            t = re.sub(r"^\s*[0-9.]+\s*\*?\s*", "", t.strip())
            if t and t.lower() != "hv" and t not in out:
                out.append(t)
    return out


def _m_ro2(path: str) -> List[str]:
    text = open(path, encoding="latin-1").read()
    m = re.search(r"RO2ToAdd\s*=\s*\{(.*?)\};", text, re.S)
    return re.findall(r"'([^']+)'", m.group(1)) if m else []


def build_addon(addon_id: str, folder: str, juice_base_fac: str, out_root: str = "") -> str:
    """Write mechanisms/addons/<addon_id>/ from the F0AM files in `folder` (see ADDON_PROFILES)."""
    import csv
    import shutil
    import mcm_converter as mc
    prof = ADDON_PROFILES[addon_id]
    out = os.path.join(out_root or os.path.join(HERE, "mechanisms", "addons"), addon_id)
    os.makedirs(out, exist_ok=True)
    base_fac = os.path.join(folder, prof["f0am_base"][0])
    mod_fac, mod_m = (os.path.join(folder, f) for f in prof["f0am_mod"])
    mod_rx = _fac_reactions(mod_fac)
    juice_all = _fac_reactions(juice_base_fac)
    base, mod, juice = _multiset(_fac_reactions(base_fac)), _multiset(mod_rx), _multiset(juice_all)
    juice_rx: Dict[tuple, list] = {}
    for rate, r in juice_all:
        if "=" in r:
            juice_rx.setdefault(_fac_key(r), []).append((rate, r))

    def fix(rate: str, reaction: str) -> Tuple[str, str]:
        for a, b in prof.get("j_rename", {}).items():
            rate = rate.replace(a, b)
        lhs, rhs = reaction.split("=", 1)
        for (reac, old), new in prof.get("product_rename", {}).items():
            if reac in [t.strip() for t in lhs.split("+")]:
                rhs = re.sub(r"(?<![A-Za-z0-9_])" + re.escape(old) + r"(?![A-Za-z0-9_])", new, rhs)
        return rate, f"{lhs}={rhs}"

    added, changed, notes = [], [], []
    for key, rates in mod.items():
        if key not in base:
            added += [fix(rt, r) for rt, r in mod_rx if "=" in r and _fac_key(r) == key]
        elif sorted(rates) != sorted(base[key]):
            if key not in juice_rx or len(juice_rx[key]) != 1 or len(rates) != 1:
                raise ValueError(f"{addon_id}: cannot place the re-rated reaction {key} on JUICE's MCM")
            if sorted(juice[key]) != sorted(base[key]):
                notes.append(f"{' + '.join(key[0])}: JUICE's MCM rate {juice[key]} differs from F0AM's {base[key]}")
            j_rate, j_rxn = juice_rx[key][0]
            new_rate, _r = fix(rates[0], j_rxn)
            if re.sub(r"\s+", "", new_rate) == re.sub(r"\s+", "", j_rate):
                notes[-1:] = [f"{' + '.join(key[0])} = {' + '.join(key[1])}: the extension's rate {new_rate} is already "
                              f"JUICE's MCM v3.3.1 rate (F0AM's MCM copy has {base[key][0]}); no override"]
                continue
            mech, rate_def = mc._convert_reaction(j_rate, j_rxn)[0]
            rate_def = rate_def.split("=CONST(", 1)[0] + "=CONST(" + mc._normalise_rate_expr(new_rate) + ")"
            changed.append((mech.strip(), rate_def, re.sub(r"\s+", "", j_rate), re.sub(r"\s+", "", new_rate)))
    conv = [mc._convert_reaction(rt, r)[0] for rt, r in added]

    juice_species = set(_species_of([r for _rt, r in juice_all if "=" in r], "="))
    used = _species_of([m for m, _ in conv], "-->")
    new_species = [s for s in used if s not in juice_species]
    needs = [s for s in used if s in juice_species]
    ro2 = [s for s in _m_ro2(mod_m) if s in new_species]
    jn = sorted({int(x) for _m, rd in conv for x in re.findall(r"\bJN(\d+)\b", rd)}
                | {int(x) for c in changed for x in re.findall(r"\bJN(\d+)\b", c[1])})
    jtab = {}
    with open(os.path.join(HERE, "photolysis", "f0am_jn_clear_sky.csv"), encoding="utf-8") as fh:
        for r in csv.DictReader(l for l in fh if not l.startswith("#")):
            jtab[r["J"]] = r
    fn = f"{addon_id.lower()}.incl"
    src = {os.path.basename(p): hashlib.sha256(open(p, "rb").read()).hexdigest() for p in (base_fac, mod_fac, mod_m)}
    L = [f"; {fn} -- {prof['name']}: the F0AM extension of the MCM v3.3.1 ({prof['f0am_path']}) as an add-on,",
         f";   built by f0am_import.py {_dt.date.today().isoformat()} from the difference to F0AM's own MCM v3.3.1.",
         f"; {prof['citation']}.  F0AM: github.com/AirChem/F0AM, GPL-3.0.",
         "; source files sha256[:16]: " + ", ".join(f"{k} {v[:16]}" for k, v in src.items()),
         f"; {len(conv)} new reactions, {len(changed)} re-rated MCM reactions, {len(new_species)} new species ({len(ro2)} in RO2).",
         "; corrections of the later F0AM .m applied: " + "; ".join(f"{a} -> {b}" for a, b in prof.get("j_rename", {}).items())
         + "; " + "; ".join(f"{k[1]} -> {v} (from {k[0]})" for k, v in prof.get("product_rename", {}).items())]
    L += [f"; {n}" for n in notes]
    L += ["", "; F0AM photolysis channels without an MCM number (sun: clear-sky fit, photolysis/f0am_jn_clear_sky.csv;",
          ";   a lamp chamber redefines them after this file)"]
    for k in jn:
        r = jtab[f"Jn{k}"]
        L.append(f"JN{k} = PHTFKT({float(r['l']):.4e},{float(r['m']):g},{-float(r['n']):g})   ; {r['reaction'][:60]}")
    L += ["", "; ---- reactions ----"] + [m.strip() for m, _ in conv]
    L += ["", "; ---- rates ----"] + [rd for _m, rd in conv]
    L += ["", "; ---- MCM reactions re-rated by this extension (last definition wins) ----"]
    for mech, rd, old, new in changed:
        L += [f"; was {old}", rd]
    with open(os.path.join(out, fn), "w", encoding="latin-1", newline="\n") as fh:
        fh.write("\n".join(L) + "\n")
    manifest = {
        "kind": "addon", "id": addon_id, "name": prof["name"],
        "version": f"F0AM {prof['f0am_path']}, built {_dt.date.today().isoformat()}", "status": "draft",
        "compatible_with": prof["compatible_with"], "files": [fn],
        "description": (f"{len(conv)} reactions of {', '.join(prof['precursors'])} and their first products; "
                        f"{len(changed)} MCM rates changed (BENZAL photolysis on IUPAC cross sections, CISOPCO2 + RO2). "
                        f"Source: F0AM {prof['f0am_path']}."),
        "adds_species": {s: prof["formulas"].get(s, "") for s in new_species},
        "needs_base_species": needs, "precursors": [p for p in prof["precursors"] if p in new_species],
        "overrides": [c[0] for c in changed], "references": [prof["citation"]],
        "notes": "; ".join(notes), "ro2": ro2, "j_f0am": jn,
        "source": {"f0am_path": prof["f0am_path"], "url": "https://github.com/AirChem/F0AM/tree/master/" + prof["f0am_path"],
                   "files_sha256": src},
    }
    with open(os.path.join(out, "manifest.json"), "w", encoding="utf-8", newline="\n") as fh:
        json.dump(manifest, fh, indent=1)
    return out


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if len(argv) == 2 and argv[0] in ADDON_PROFILES:      # an add-on: python f0am_import.py NOAA_BB <folder with the F0AM files>
        import mechanism_registry as mreg
        out = build_addon(argv[0], argv[1], mreg.base(mreg.DEFAULT_BASE).fac)
        m = json.load(open(os.path.join(out, "manifest.json")))
        print(f"{m['id']}: {len(m['adds_species'])} new species, {len(m['overrides'])} re-rated MCM reactions -> {out}")
        if m["notes"]:
            print("  notes:", m["notes"])
        return 0
    if len(argv) != 2 or argv[0] not in PROFILES:
        print(__doc__); print("bases:", ", ".join(PROFILES), "| add-ons:", ", ".join(ADDON_PROFILES)); return 2
    out = write_base(argv[0], argv[1])
    m = json.load(open(os.path.join(out, "manifest.json")))
    print(f"{m['id']}: {m['lumped']['n_species']} species, {m['lumped']['n_reactions']} reactions -> {out}")
    for d in m["lumped"]["dropped"]:
        print("  not converted:", d)
    if m["notes"]:
        print("  notes:", m["notes"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
