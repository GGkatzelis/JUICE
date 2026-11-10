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

# the standard SAPHIR aliases (as geckoa.aliases.incl): canonical = what the SAPHIR includes and the outputs use
STANDARD_ALIASES = [("HV", "hv"), ("OH", "HO"), ("HCHO", "H2CO CH2O"), ("O", "O3P"), ("CH3OOH", "CH3O2H MOOH"),
                    ("HONO", "HNO2"), ("HNO3", "HONO2"), ("HNO4", "HO2NO2"), ("SA", "SULF H2SO4")]
CONSTANT_SPECIES = ("M", "O2", "N2", "H2", "H2O")
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
def read_allrxns(path: str, fixed: Optional[Dict[str, float]] = None) -> Tuple[List[str], List[Reaction], List[str]]:
    """(species, reactions, dropped) of an F0AM <mech>_AllRxns.m."""
    fixed = dict(fixed or {})
    text = open(path, encoding="latin-1").read()
    m = re.search(r"SpeciesToAdd\s*=\s*\{(.*?)\};", text, re.S)
    species = re.findall(r"'([^']+)'", m.group(1)) if m else []
    blocks = re.split(r"\n\s*i\s*=\s*i\s*\+\s*1\s*;", text)[1:]
    reactions: List[Reaction] = []
    dropped: List[str] = []
    rate_by_index: Dict[int, str] = {}                 # F0AM reaction number -> EASY rate (k(:,n) refers to it)
    for index, blk in enumerate(blocks, start=1):
        label_m = re.search(r"Rnames\{\s*(\d+)\s*\}\s*=\s*'([^']*)'", blk)
        name = label_m.group(2).strip() if label_m else "?"
        if label_m:
            index = int(label_m.group(1))
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
        for sp_, sign, val in re.findall(r"f([A-Za-z0-9_]+)\(i\)\s*=\s*f\1\(i\)\s*([-+])\s*([0-9.Ee+-]+)", blk):
            net[sp_] = net.get(sp_, 0.0) + (1 if sign == "+" else -1) * float(val)
        prods: "OrderedDict[str, float]" = OrderedDict()
        for sp_, v in net.items():
            c = v + reac.count(sp_)
            if abs(c) > 1e-12:
                prods[sp_] = round(c, 6)
        for sp_ in reac:
            if net.get(sp_, 0.0) + reac.count(sp_) < -1e-12:
                notes.append(f"{sp_}: consumed more than the reactants say")
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
    """Named rate constants of <mech>_K.m (constant values only)."""
    out: Dict[str, str] = {}
    if not os.path.isfile(path):
        return out
    text = open(path, encoding="latin-1").read()
    for name, val in re.findall(r"Knames\{i\}\s*=\s*'([^']+)';\s*\n\s*krx\(:,i\)\s*=\s*([^;]+);", text):
        out[name] = matlab_to_easy(val)
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
    species, reactions, dropped = read_allrxns(rx, profile.get("fixed"))
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
            "O2=CONST(0.2095*M)", "N2=CONST(0.7808*M)", "H2=CONST(5.6E-7*M)", ""]))

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
                L.append(f"{n} = CONST({v})")
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
        fh.write("Name\tFormula\tMCM name\n" + "".join(f"{s}\t\t{rev.get(s, '')}\n" for s in species))

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


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if len(argv) != 2 or argv[0] not in PROFILES:
        print(__doc__); print("known:", ", ".join(PROFILES)); return 2
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
