"""gecko_import.py -- a GECKO-A generator output folder becomes a base mechanism of the registry.

The front door for GECKO-A chemistry: `import_generator_folder()` runs F. Rohrer's converter
(gecko_mechanism / gecko_to_easy) once on an OUT/<mech>/ folder and files the result under
mechanisms/base/<ID>/ with a manifest of kind "gecko", so the Chemistry step offers it next to the
MCM download and every experiment that chooses it records what it used (generator files and their
SHA-256, generator depth, import options).  Nothing of GECKO-A's own code is involved: the import
reads the generator's output files.

Files written (the names the experiment writer copies into <run>/chemistry/):
    geckoa.physics.incl             CONST, DENSITY, M = DENSITY(P,T), O2, N2 - the lines an MCM include set carries in its
                                    .saphir.incl; every termolecular rate (M, O2) is zero without them.  Included first.
    anorg.geckoa.incl               GECKO-A's inorganic core          (used when the core choice is "gecko")
    geckoa.aliases.incl             GECKO-A names -> EASY/MCM names   (always loaded first)
    <ID>.geckoa.names.incl          GECKO code -> chosen name synonyms
    <ID>.geckoa.organic.incl        the organic chemistry
    <ID>.dil.mech / <ID>.dil.rate   X + DIL --> for every gas species, in the IDL layout the EASY writer expects
    <ID>.SOA.include                gas <-> particle <-> wall partitioning (Nannoolal C*), optional per experiment
    <ID>.species.tsv                species table with SMILES

    python gecko_import.py <generator folder> --id GECKOA_isoprene --voc C5H8 [--maxgen N] [--cstar-max X] [--wall krechmer|gecko]
"""
from __future__ import annotations

import argparse
import os
import re
from datetime import date
from typing import Dict, List, Optional

import mechanism_registry as mreg

GECKO_PREFIX = "GECKOA_"
CITATION = ("GECKO-A: Aumont, Szopa & Madronich, ACP 5, 2497 (2005); Valorso et al., ACP 11, 6895 (2011); "
            "GECKO-A v1.0 (2026), https://gitlab.in2p3.fr/ (GPL-3.0). JUICE reads the generator output; "
            "it ships none of its code.")
# primary codes GECKO-A names explicitly (the generic C<n> codes such as O02000 carry no name): the MCM name to suggest
PRIMARY_NAMES = {"ISOPRE": "C5H8", "APINEN": "APINENE", "BPINEN": "BPINENE", "LIMONE": "LIMONENE", "TOLUEN": "TOLUENE",
                 "BENZEN": "BENZENE", "OXYLEN": "OXYL", "MXYLEN": "MXYL", "PXYLEN": "PXYL"}
NOTATION_NAMES = {"CH3CH2(OH)": "C2H5OH", "CH3(OH)": "CH3OH", "CH3CHO": "CH3CHO", "CH3COCH3": "CH3COCH3",
                  "CH2=CH2": "C2H4", "CH3CH=CH2": "C3H6", "CH3CH2CH2CH2CH2CH2CH2CH3": "NC8H18", "CH3CH2CH2CH3": "NC4H10"}


def generator_summary(folder: str) -> dict:
    """What the generator produced, from its small bookkeeping files (no mechanism read)."""
    out = {"folder": folder, "species": None, "reactions": None, "maxgen": None, "primaries": []}
    p = os.path.join(folder, "size.dum")
    if os.path.isfile(p):
        for line in open(p, encoding="latin-1", errors="replace"):
            m = re.match(r"\s*(\d+)\s*!\s*(.*)", line)
            if m and "total number of species" in m.group(2):
                out["species"] = int(m.group(1))
            elif m and "total number of reactions" in m.group(2):
                out["reactions"] = int(m.group(1))
    p = os.path.join(folder, "scheme.log")
    if os.path.isfile(p):
        for line in open(p, encoding="latin-1", errors="replace"):
            m = re.match(r"\s*maxgen\s*=\s*(\d+)", line)
            if m:
                out["maxgen"] = int(m.group(1))
    p = os.path.join(folder, "listprimary.dat")
    if os.path.isfile(p):
        for line in open(p, encoding="latin-1", errors="replace"):
            toks = line.split()
            if toks:
                out["primaries"].append((toks[0], toks[1] if len(toks) > 1 else ""))
    return out


def is_generator_folder(folder: str) -> bool:
    return all(os.path.isfile(os.path.join(folder, f)) for f in ("dictionary.out", "reactionswithcom.dum"))


def suggest_voc_name(folder: str) -> str:
    """MCM name for the generator's primary species, when it is known; '' otherwise (the user types it)."""
    import gecko_to_easy as g
    for code, notation in generator_summary(folder)["primaries"]:
        if code in g.MCM_NAMES:
            return g.MCM_NAMES[code]
        if code in PRIMARY_NAMES:
            return PRIMARY_NAMES[code]
        if notation in NOTATION_NAMES:
            return NOTATION_NAMES[notation]
    return ""


def suggest_id(folder: str, voc: str = "") -> str:
    tag = (voc or os.path.basename(os.path.normpath(folder)) or "mechanism")
    return GECKO_PREFIX + re.sub(r"[^A-Za-z0-9_]+", "_", tag).strip("_")


PHYSICS_INCLUDE = "geckoa.physics.incl"
PHYSICS_LINES = [
    "; geckoa.physics.incl -- physics, not chemistry: what an MCM include set defines in its <stem>.saphir.incl and a",
    ";   GECKO-A mechanism has to bring itself. Air density from the pressure and temperature settings, O2 and N2 for",
    ";   the termolecular rates (GECKO-A TBODY/OXYGEN, MCM KMT fall-off). Written by gecko_import.py.",
    "CONST=FUNCTION[%1]",
    "DENSITY=FUNCTION[(%1)/((%2)*1.379E-19)]",
    "M =DENSITY(P,T)",
    "O2=CONST(0.20946*M)",
    "N2=CONST(0.78084*M)",
]


def write_physics_include(path: str) -> str:
    with open(path, "w", encoding="latin-1", newline="\n") as fh:
        fh.write("\n".join(PHYSICS_LINES) + "\n")
    return path


def _split_dil(dil_incl: str, mech_path: str, rate_path: str) -> None:
    """The converter's one dilution include -> the .dil.mech / .dil.rate pair the EASY writer references."""
    head, rxn, rate = [], [], []
    for line in open(dil_incl, encoding="latin-1"):
        s = line.rstrip("\r\n")
        code = s.split(";", 1)[0].strip()
        if not code:
            if not rxn and not rate:
                head.append(s)
            continue
        (rate if code.startswith("k[") else rxn).append(s)
    with open(mech_path, "w", encoding="latin-1", newline="\n") as fh:
        fh.write("\n".join(head + rxn) + "\n")
    with open(rate_path, "w", encoding="latin-1", newline="\n") as fh:
        fh.write("\n".join(head + rate) + "\n")


def import_generator_folder(folder: str, base_id: str = "", voc: str = "", name: str = "", maxgen: Optional[int] = None,
                            cstar_max: Optional[float] = None, wall: str = "krechmer", base_dir: str = "",
                            log=print) -> mreg.BaseMechanism:
    """Convert a GECKO-A generator folder and register it as base mechanism `base_id`.

    voc      MCM name given to the generator's primary species (what the Injections step will call it)
    maxgen   static generation cut (None = the complete generator output)
    cstar_max  partition only species with C* <= this (ug/m3); None = every species GECKO-A partitions
    wall     'krechmer' (project scheme) or 'gecko' (GECKO-A's own wall parameters)"""
    import gecko_to_easy as g
    from gecko_mechanism import read_mechanism, generation_cut
    folder = os.path.abspath(folder)
    if not is_generator_folder(folder):
        raise FileNotFoundError(f"{folder} is not a GECKO-A generator output folder (dictionary.out, reactionswithcom.dum)")
    base_id = base_id or suggest_id(folder, voc)
    if not re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", base_id):
        raise ValueError(f"mechanism id '{base_id}' may contain letters, digits and _ only")
    out = os.path.join(base_dir or mreg.BASE_DIR, base_id)
    os.makedirs(out, exist_ok=True)
    mech = read_mechanism(folder)
    n_full = (len(mech.species), len(mech.reactions))
    primary = mech.primaries[0] if len(mech.primaries) == 1 else ""
    voc = (voc or "").strip().upper()
    if voc and primary and primary not in g.MCM_NAMES:
        g.MCM_NAMES[primary] = voc                    # the generator's primary gets the experiment's name for it
    if maxgen is not None:
        mech, _gens = generation_cut(mech, maxgen)
        log(f"[gecko_import] generation cut maxgen={maxgen}: {n_full[0]} -> {len(mech.species)} species, "
            f"{n_full[1]} -> {len(mech.reactions)} reactions")
    conv = g.convert(mech)
    files = {"physics": PHYSICS_INCLUDE, "anorg": g.ANORG_INCLUDE, "aliases": g.ALIAS_INCLUDE, "names": f"{base_id}.geckoa.names.incl",
             "organic": f"{base_id}.geckoa.organic.incl", "dil_mech": f"{base_id}.dil.mech", "dil_rate": f"{base_id}.dil.rate",
             "soa": f"{base_id}.SOA.include", "species_tsv": f"{base_id}.species.tsv"}
    write_physics_include(os.path.join(out, files["physics"]))
    n_an = g.write_anorg_include(conv, mech, os.path.join(out, files["anorg"]))
    g.write_alias_include(os.path.join(out, files["aliases"]))
    n_names = g.write_names_include(mech, os.path.join(out, files["names"]))
    n_org = g.write_organic_include(conv, mech, os.path.join(out, files["organic"]))
    dil_incl = os.path.join(out, f"{base_id}.dil.incl")
    g.write_dil_include(conv.gas_species, dil_incl)
    _split_dil(dil_incl, os.path.join(out, files["dil_mech"]), os.path.join(out, files["dil_rate"]))
    os.remove(dil_incl)
    soa = g.write_soa_include(mech, os.path.join(out, files["soa"]), cstar_max, wall)
    n_tsv = g.write_species_tsv(mech, soa, os.path.join(out, files["species_tsv"]))
    gen_depth = g._generator_maxgen(folder)
    primary_name = g.easy_name(primary) if primary else ""
    source_files = {f: mreg.sha256_of(os.path.join(folder, f)) for f in ("dictionary.out", "reactionswithcom.dum", "pvap.nan.dat", "scheme.log")
                    if os.path.isfile(os.path.join(folder, f))}
    manifest = {
        "kind": "gecko",
        "id": base_id,
        "name": name or f"GECKO-A {voc or primary_name or os.path.basename(folder)}"
                        + (f" (generation cut {maxgen})" if maxgen is not None else " (complete generator output)"),
        "version": f"generator maxgen {gen_depth if gen_depth is not None else '?'}, imported {date.today().isoformat()}",
        "files": files,
        "checksums": {fn: mreg.sha256_of(os.path.join(out, fn)) for fn in files.values()},
        "source": {"folder": folder, "generator_files_sha256": source_files, "generator_maxgen": gen_depth,
                   "imported": date.today().isoformat(),
                   "options": {"maxgen_cut": maxgen, "cstar_max": cstar_max, "wall": wall}},
        "gecko": {"precursor": voc or primary_name, "primary_code": primary,
                  "primary_notation": mech.species[primary].notation if primary in mech.species else "",
                  "n_species": len(mech.species), "n_species_generator": n_full[0], "n_reactions_generator": n_full[1],
                  "n_reactions_organic": n_org, "n_reactions_inorganic": n_an, "n_partitioning": len(soa), "n_names": n_names,
                  "n_species_tsv": n_tsv, "dropped": dict(conv.dropped), "rescaled": len(conv.rescaled)},
        "citation": CITATION,
        "status": "imported",
        "notes": "Written by gecko_import.py from the generator output; regenerate rather than edit. "
                 "Experiments choose the inorganic core (GECKO-A's own or MCM v3.3.1 with aliases) and whether the "
                 "partitioning include is loaded.",
    }
    mreg.write_manifest(out, manifest)
    log(f"[gecko_import] {base_id}: {len(mech.species)} species, {n_an} inorganic + {n_org} organic reactions, "
        f"{len(soa)} partitioning species -> {out}")
    bases = mreg.load_bases(base_dir or "")
    return bases[base_id]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("folder", help="GECKO-A generator output folder")
    ap.add_argument("--id", default="", help="mechanism id (default GECKOA_<voc>)")
    ap.add_argument("--voc", default="", help="MCM name of the primary species (default: suggested from the generator files)")
    ap.add_argument("--name", default="", help="display name")
    ap.add_argument("--maxgen", type=int, default=None, help="static generation cut")
    ap.add_argument("--cstar-max", type=float, default=None, help="partition only species with C* <= this (ug/m3)")
    ap.add_argument("--wall", choices=("krechmer", "gecko"), default="krechmer")
    ap.add_argument("--summary", action="store_true", help="only print what the generator produced")
    a = ap.parse_args(argv)
    if a.summary:
        s = generator_summary(a.folder)
        print(f"{s['species']} species, {s['reactions']} reactions, maxgen {s['maxgen']}, primaries {s['primaries']}, "
              f"suggested name {suggest_voc_name(a.folder) or '(none)'}")
        return
    b = import_generator_folder(a.folder, a.id, a.voc or suggest_voc_name(a.folder), a.name, a.maxgen, a.cstar_max, a.wall)
    print(f"registered {b.id}: {b.name} ({b.version})")


if __name__ == "__main__":
    main()
