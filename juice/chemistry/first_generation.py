"""First-generation chemistry of a compound in any base mechanism: the reactions that consume it.

The MCM base is read from its .fac (mcm_species_db.reactions_of). Every other mechanism JUICE ships is EASY
include text, a reaction line followed (directly, or in a block further down) by its rate line:

    API + OH --> APIP
    k[API + OH --> APIP] = CONST(1.21E-11*EXP(440/T))

That covers GECKO-A, the F0AM-library mechanisms (CRACMM, RACM2) and the add-ons (NOAA_BB). Names are matched
through the mechanism's ALIAS lines, so the MCM name works everywhere: APINENE finds RACM2's API, OH finds
GECKO-A's HO. The files are scanned line by line, never parsed whole, so a GECKO-A mechanism of a million
reactions costs one pass over the file.

Each reaction is a dict like mcm_species_db.reactions_of gives: easy, easy_rate, reactants, products,
photolysis, partner ('OH', 'O3', 'NO3', 'hv', another species, or '' for unimolecular).
"""
from __future__ import annotations

import os
import re
from typing import Dict, Iterable, List, Tuple

from juice.chemistry import mechanism_registry as mreg

_ALIAS = re.compile(r"^\s*ALIAS\[\s*([^\]]+?)\s*\]\s*=\s*(.+)$", re.I)
_RATE = re.compile(r"^\s*k\[(.+?)\]\s*=\s*(.*)$", re.I)
_STOICH = re.compile(r"^[0-9.eEdD+\-]*\*")
_CONST = re.compile(r"^CONST\((.*)\)$", re.I | re.S)
ORDER = ("OH", "O3", "NO3", "hv", "")
_ALIAS_CACHE: Dict[Tuple[str, ...], Dict[str, str]] = {}


def _logical_lines(path: str) -> Iterable[str]:
    """Lines of an EASY file without comments; a line ending in '$' continues on the next."""
    buf = ""
    with open(path, encoding="latin-1") as fh:
        for raw in fh:
            code = raw.split(";", 1)[0].rstrip()
            if code.endswith("$"):
                buf += code[:-1] + " "
                continue
            line, buf = buf + code, ""
            if line.strip():
                yield line


def _name(term: str) -> str:
    return _STOICH.sub("", term.strip()).upper()


def _key(reaction: str) -> str:
    return re.sub(r"\s+", "", reaction).upper()


def aliases(paths: List[str]) -> Dict[str, str]:
    """synonym -> name used in the mechanism, from the ALIAS[name] = syn1 syn2 lines of these files."""
    k = tuple(paths)
    if k not in _ALIAS_CACHE:
        out: Dict[str, str] = {}
        for p in paths:
            if not os.path.isfile(p):
                continue
            with open(p, encoding="latin-1") as fh:
                for raw in fh:
                    m = _ALIAS.match(raw.split(";", 1)[0])
                    if m:
                        for syn in m.group(2).split():
                            out.setdefault(syn.upper(), m.group(1).strip().upper())
        _ALIAS_CACHE[k] = out
    return _ALIAS_CACHE[k]


def easy_reactions_of(species: str, paths: List[str], alias: Dict[str, str], limit: int = 0) -> List[dict]:
    """Reactions in these EASY files with `species` (any synonym) as a reactant."""
    canon = lambda n: alias.get(n, n)
    want = canon(species.strip().upper())
    names_of_want = {want} | {s for s, c in alias.items() if c == want}
    found: List[dict] = []
    rates: Dict[str, str] = {}
    for p in paths:
        if not os.path.isfile(p):
            continue
        for line in _logical_lines(p):
            up = line.upper()
            if not any(n in up for n in names_of_want):        # cheap filter before any parsing
                continue
            m = _RATE.match(line)
            if m:
                rates[_key(m.group(1))] = re.sub(r"\s+", " ", m.group(2).strip())
                continue
            if "-->" not in line or "=" in line:
                continue
            lhs, rhs = line.split("-->", 1)
            reactants = [_name(t) for t in lhs.split("+") if t.strip()]
            if want not in [canon(r) for r in reactants]:
                continue
            found.append({"easy": re.sub(r"\s+", " ", line.strip()), "reactants": reactants,
                          "products": [t.strip() for t in rhs.split("+") if t.strip()]})
            if limit and len(found) >= limit:
                break
    out = []
    for r in found:
        rate = rates.get(_key(r["easy"]), "")
        m = _CONST.match(rate)
        r["easy_rate"] = m.group(1).strip() if m else rate
        others = [canon(n) for n in r["reactants"] if canon(n) != want]
        r["photolysis"] = any(o in ("HV",) for o in others) or (not others and bool(re.match(r"\s*J", r["easy_rate"], re.I)))
        others = [o for o in others if o != "HV"]
        r["partner"] = "hv" if r["photolysis"] else (others[0] if others else "")
        out.append(r)
    out.sort(key=lambda r: (ORDER.index(r["partner"]) if r["partner"] in ORDER else 9, r["easy"]))
    return out


def _files(base: "mreg.BaseMechanism", addon_ids: List[str]) -> Tuple[List[str], List[str]]:
    """(EASY files holding reactions, files holding ALIAS lines) of a base and its add-ons."""
    rx = []
    if not base.fac:                                   # a complete base: every include of its folder
        rx = [os.path.join(base.folder, f) for f in sorted(os.listdir(base.folder))
              if f.lower().endswith((".incl", ".include"))
              and not any(s in f.lower() for s in ("physics", ".soa."))]   # physics and gas-particle-wall steps
    for a in mreg.addons_for(base.id, list(addon_ids)) if addon_ids else []:
        rx += a.paths()
    return rx, rx


_CACHE: Dict[tuple, Tuple[str, List[dict]]] = {}


def first_generation(species: str, base_id: str = "", addon_ids: Iterable[str] = (), limit: int = 0) -> Tuple[str, List[dict]]:
    """(name in the mechanism, reactions consuming it) for `species` in base `base_id` with the given add-ons."""
    key = (species.strip().upper(), base_id, tuple(addon_ids), limit)
    if key not in _CACHE:
        _CACHE[key] = _first_generation(*key)
    return _CACHE[key]


def _first_generation(species: str, base_id: str, addon_ids: Tuple[str, ...], limit: int) -> Tuple[str, List[dict]]:
    b = mreg.base(base_id)
    rx_files, alias_files = _files(b, list(addon_ids))
    alias = aliases(alias_files)
    name = alias.get(species.strip().upper(), species.strip().upper())
    out: List[dict] = []
    if b.fac:                                          # the MCM: its .fac, then the add-ons' EASY reactions
        from juice.chemistry.mcm_species_db import reactions_of
        out = reactions_of(name, b.fac)
    if rx_files:
        out += easy_reactions_of(name, rx_files, alias, limit=limit)
        out.sort(key=lambda r: (ORDER.index(r["partner"]) if r["partner"] in ORDER else 9, r["easy"]))
    return name, out
