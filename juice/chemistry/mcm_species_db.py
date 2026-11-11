"""mcm_species_db.py — Loader for the MCM v3.3.1 species database.

Parses mcm_export_species.tsv (downloaded directly from http://mcm.york.ac.uk/)
and exposes a single dict:

    MCM_DB: dict[str, MCMSpecies]

Keys are uppercase MCM species names.  Access via:

    from mcm_species_db import MCM_DB
    entry = MCM_DB.get("TOLUENE")   # → MCMSpecies or None
    n_c   = entry.n_C               # 7
    mass  = entry.mass              # 92.14
"""
from __future__ import annotations
from juice._root import ROOT as _JUICE_ROOT

import os
import re
from dataclasses import dataclass, field
from typing import Optional

# ---------------------------------------------------------------------------
# Path to the TSV (same directory as this file)
# ---------------------------------------------------------------------------
_DB_DIR = os.path.join(_JUICE_ROOT, "chemistry", "mechanisms", "base", "MCM_v3.3.1")
_TSV_FILE = os.path.join(_DB_DIR, "mcm_export_species.tsv")
if not os.path.isfile(_TSV_FILE):                               # pre-registry layout
    _DB_DIR = os.path.join(_JUICE_ROOT, "MCM_database")
    _TSV_FILE = os.path.join(_DB_DIR, "mcm_export_species.tsv")


# ---------------------------------------------------------------------------
# Data container
# ---------------------------------------------------------------------------
@dataclass
class MCMSpecies:
    name:       str
    formula:    str          # e.g. "C7H8"
    smiles:     str
    mass:       float        # g/mol
    is_peroxy:  bool         # True = RO2 (PeroxyRadical column)
    synonyms:   str
    # Atom counts derived from formula
    n_C: int = 0
    n_H: int = 0
    n_N: int = 0
    n_O: int = 0
    n_Cl: int = 0
    n_Br: int = 0
    n_S: int = 0
    n_F: int = 0


# ---------------------------------------------------------------------------
# Formula parser
# ---------------------------------------------------------------------------
_ATOM_RE = re.compile(
    r'(Cl|Br|Si|[CHNOFSP])(\d*)'   # two-char atoms first, then one-char
)
_ATOM_MAP = {
    'C':  'n_C', 'H': 'n_H', 'N': 'n_N', 'O': 'n_O',
    'Cl': 'n_Cl', 'Br': 'n_Br', 'S': 'n_S', 'F': 'n_F',
}


def _parse_formula(formula: str) -> dict[str, int]:
    """Return {n_C, n_H, …} atom counts from a molecular formula string."""
    counts: dict[str, int] = {}
    if not formula or formula in ('-', 'X', ''):
        return counts
    for m in _ATOM_RE.finditer(formula):
        atom, num = m.group(1), m.group(2)
        key = _ATOM_MAP.get(atom)
        if key:
            counts[key] = counts.get(key, 0) + (int(num) if num else 1)
    return counts


# ---------------------------------------------------------------------------
# TSV parser
# ---------------------------------------------------------------------------
def _find_col(parts_lower: list[str], *candidates: str) -> int:
    """Return the first column index matching any of the candidate substrings."""
    for cand in candidates:
        for i, p in enumerate(parts_lower):
            if cand in p:
                return i
    return -1


def _load_tsv(path: str) -> dict[str, MCMSpecies]:
    db: dict[str, MCMSpecies] = {}
    if not os.path.isfile(path):
        return db

    with open(path, encoding='utf-8', errors='replace') as fh:
        header_found = False
        col_name = col_smiles = col_formula = col_mass = col_peroxy = col_synonyms = -1

        for raw in fh:
            line = raw.rstrip('\n')

            # Skip blank or pure-comment lines (but still scan for header)
            if not line or line.startswith('*') or line.startswith(';'):
                continue

            parts = line.split('\t')
            parts_lower = [p.strip().lower() for p in parts]

            if not header_found:
                # Detect header row: first cell must be exactly "name"
                if parts_lower[0] == 'name':
                    col_name     = _find_col(parts_lower, 'name')
                    col_smiles   = _find_col(parts_lower, 'smiles')
                    col_formula  = _find_col(parts_lower, 'formula')
                    col_mass     = _find_col(parts_lower, 'mass')
                    col_peroxy   = _find_col(parts_lower, 'peroxyradical')  # matches "excited peroxyradical"
                    col_synonyms = _find_col(parts_lower, 'synonyms')
                    header_found = col_name >= 0 and col_formula >= 0 and col_mass >= 0
                continue

            # Data row
            parts = line.split('\t')
            if len(parts) <= max(col_name, col_formula, col_mass, col_peroxy):
                continue

            name = parts[col_name].strip()
            if not name or name.startswith('*') or name.startswith(';'):
                continue

            formula  = parts[col_formula].strip()  if col_formula  < len(parts) else ''
            smiles   = parts[col_smiles].strip()   if col_smiles   < len(parts) else ''
            synonyms = parts[col_synonyms].strip() if col_synonyms < len(parts) else ''
            try:
                mass = float(parts[col_mass].strip())
            except (ValueError, IndexError):
                mass = 0.0
            is_peroxy = parts[col_peroxy].strip().lower() == 'true' if col_peroxy < len(parts) else False

            atoms = _parse_formula(formula)
            entry = MCMSpecies(
                name=name.upper(),
                formula=formula,
                smiles=smiles,
                mass=mass,
                is_peroxy=is_peroxy,
                synonyms=synonyms,
                **atoms,
            )
            db[entry.name] = entry

    return db


# ---------------------------------------------------------------------------
# Module-level singleton — loaded once on first import
# ---------------------------------------------------------------------------
MCM_DB: dict[str, MCMSpecies] = _load_tsv(_TSV_FILE)


# ---------------------------------------------------------------------------
# Formula → species name reverse index (built lazily on first use)
# ---------------------------------------------------------------------------
_FORMULA_INDEX: dict | None = None


def _normalise_formula(formula: str) -> str:
    """Canonicalise a chemical formula for lookup.

    Parses into atom counts then re-emits in Hill order:
    C first, H second, then other elements alphabetically. So 'C5H8',
    'H8C5', 'C5 H8', 'c5h8' all canonicalise to 'C5H8'.
    """
    # Normalise case: capitalise atom symbols (h8c5 → H8C5)
    # Two-letter atoms (Cl, Br, Si) need special handling
    s = formula.strip()
    s = re.sub(r'\bcl\b', 'Cl', s, flags=re.I)
    s = re.sub(r'\bbr\b', 'Br', s, flags=re.I)
    # Then upper-case single letters
    s = re.sub(r'(?<![a-zA-Z])([cChHnNoOsSfFpP])', lambda m: m.group(1).upper(), s)
    counts = _parse_formula(s)
    if not counts:
        return ''
    out = []
    if 'n_C' in counts:
        out.append('C' + (str(counts['n_C']) if counts['n_C'] > 1 else ''))
    if 'n_H' in counts:
        out.append('H' + (str(counts['n_H']) if counts['n_H'] > 1 else ''))
    rest = sorted(k for k in counts if k not in ('n_C', 'n_H'))
    for k in rest:
        sym = k.replace('n_', '')
        out.append(sym + (str(counts[k]) if counts[k] > 1 else ''))
    return ''.join(out)


def lookup_by_formula(formula: str) -> list[MCMSpecies]:
    """Return all MCM species matching the given chemical formula.

    The formula is canonicalised in Hill order before lookup, so 'C5H8',
    'c5h8', 'H8C5', 'C5 H8' all return the same species list. Results
    are sorted by name.
    """
    global _FORMULA_INDEX
    if _FORMULA_INDEX is None:
        _FORMULA_INDEX = {}
        for sp in MCM_DB.values():
            key = _normalise_formula(sp.formula)
            if key:
                _FORMULA_INDEX.setdefault(key, []).append(sp.name)
        for k in _FORMULA_INDEX:
            _FORMULA_INDEX[k].sort()

    canon = _normalise_formula(formula)
    if not canon:
        return []
    names = _FORMULA_INDEX.get(canon, [])
    return [MCM_DB[n] for n in names]


# ---------------------------------------------------------------------------
# Known zero-carbon inorganic/boundary species not in the MCM TSV export.
# These are SAPHIR measurement columns, model constraints, or atmospheric
# background species — definitively carbon-free, never flagged as "unknown".
# ---------------------------------------------------------------------------
KNOWN_ZERO_CARBON: frozenset[str] = frozenset({
    # NOx / HONO / HNOx
    'NO', 'NO2', 'NO3', 'N2O5', 'N2O', 'HONO', 'HNO3', 'HNO4',
    # Oxygen species
    'O', 'O1D', 'O3', 'OH', 'HO2', 'H2O2',
    # Water / hydrogen
    'H2O', 'H2', 'H',
    # Sulfur (inorganic)
    'SO2', 'SO3', 'H2SO4',
    # Halogens
    'CL', 'CL2', 'HCL', 'CLNO', 'CLNO2', 'CLNO3', 'BR', 'BR2', 'HBR',
    # Atmospheric diatomics / noble gases
    'N2', 'O2', 'AR', 'HE', 'N', 'NH3',
    # SAPHIR model internals
    'M', 'RO2', 'RO2T',
    # Common aliases / measurement columns
    'NA', 'SA',
})


# ---------------------------------------------------------------------------
# Reactions of a species (first-generation products) from the full MCM .fac
# ---------------------------------------------------------------------------
_FAC_CACHE: dict = {}
_FAC_LINE = re.compile(r"^\s*%\s*(?P<rate>.+?)\s*:\s*(?P<lhs>.*?)\s*=\s*(?P<rhs>.*?)\s*;")
_PARTNERS = ("OH", "O3", "NO3", "HV", "NO", "HO2", "NO2", "CL")


def _fac_species(side: str) -> list[str]:
    """'HCHO + 2 HO2' → ['HCHO', '2*HO2']  (EASY coefficient style)."""
    out = []
    for tok in side.split("+"):
        tok = tok.strip()
        if not tok:
            continue
        m = re.match(r"^([0-9.]+)\s+(\S+)$", tok)
        out.append(f"{m.group(1)}*{m.group(2)}" if m else tok)
    return out


def fac_to_easy_rate(rate: str) -> str:
    """MCM .fac rate text → the form the EASY include files use (TEMP → T, J<n> → Jn)."""
    r = re.sub(r"J<(\d+)>", r"J\1", rate)
    return r.replace("TEMP", "T")


def load_fac_reactions(fac_path: str) -> list[dict]:
    """All reactions of a .fac file: {'rate', 'reactants', 'products', 'easy'} (cached)."""
    if fac_path in _FAC_CACHE:
        return _FAC_CACHE[fac_path]
    rx = []
    with open(fac_path, encoding="latin-1") as fh:
        for line in fh:
            m = _FAC_LINE.match(line)
            if not m:
                continue
            lhs, rhs = _fac_species(m.group("lhs")), _fac_species(m.group("rhs"))
            rate = m.group("rate").strip()
            photo = rate.startswith("J<")
            easy_lhs = lhs + (["HV"] if photo else [])
            rx.append({"rate": rate, "reactants": lhs, "products": rhs, "photolysis": photo,
                       "easy": " + ".join(easy_lhs) + " --> " + " + ".join(rhs),
                       "easy_rate": fac_to_easy_rate(rate)})
    _FAC_CACHE[fac_path] = rx
    return rx


def reactions_of(species: str, fac_path: str) -> list[dict]:
    """Reactions consuming `species`, each with a 'partner' (OH, O3, NO3, hv, …, or '' for unimolecular)."""
    sp = species.strip().upper()
    out = []
    for r in load_fac_reactions(fac_path):
        names = [t.split("*")[-1].upper() for t in r["reactants"]]
        if sp not in names:
            continue
        others = [n for n in names if n != sp]
        partner = "hv" if r["photolysis"] else (others[0] if others else "")
        out.append(dict(r, partner=partner))
    order = {p: i for i, p in enumerate(("OH", "O3", "NO3", "hv", ""))}
    out.sort(key=lambda r: (order.get(r["partner"], 9), r["easy"]))
    return out


# ---------------------------------------------------------------------------
# Species search for surrogates: exact formula, then near formula, then name
# ---------------------------------------------------------------------------
def _counts(formula: str) -> dict[str, int]:
    try:
        return _parse_formula(formula)
    except Exception:
        return {}


def search_species(query: str, limit: int = 40) -> list[tuple["MCMSpecies", str]]:
    """Candidates for `query` (a formula such as C6H14O3, or part of an MCM name).
    Returns [(species, why)] best first: exact formula → same C and O count →
    same C count → name match.  Only non-radical, closed-shell species are listed
    first; peroxy radicals come last."""
    q = (query or "").strip()
    if not q:
        return []
    out: list[tuple[int, str, "MCMSpecies"]] = []
    looks_formula = bool(re.fullmatch(r"[A-Z][A-Za-z0-9]*", q)) and any(ch.isdigit() for ch in q)
    qn = _normalise_formula(q) if looks_formula else ""
    qc = _counts(q) if looks_formula else {}
    qu = q.upper()
    for name, sp in MCM_DB.items():
        rank, why = None, ""
        if qn and _normalise_formula(sp.formula or "") == qn:
            rank, why = 0, f"same formula {sp.formula}"
        elif qc and sp.n_C == qc.get("n_C", -1) and sp.n_O == qc.get("n_O", -2) and sp.n_N == qc.get("n_N", 0):
            rank, why = 1, f"same C and O count ({sp.formula})"
        elif qc and sp.n_C == qc.get("n_C", -1) and abs(sp.n_O - qc.get("n_O", 0)) <= 1 and sp.n_N == qc.get("n_N", 0):
            rank, why = 2, f"same C count, O ±1 ({sp.formula})"
        if rank in (1, 2) and qc:               # a C14 alkane should not become a C14 ketone with 3 double bonds
            d_dbe = abs((sp.n_C - sp.n_H / 2.0 + sp.n_N / 2.0) - (qc.get("n_C", 0) - qc.get("n_H", 0) / 2.0 + qc.get("n_N", 0) / 2.0))
            if d_dbe > 2:
                rank = None
        elif qu in name:
            rank, why = 3, f"name contains {qu}"
        if rank is None:
            continue
        if any(tag in name for tag in ("OOH", "NO3", "PAN", "OOA", "OOB")):
            rank += 5                       # oxidation products: possible, but an emitted compound is the better surrogate
        if sp.is_peroxy or name.endswith("O2") or name.endswith("O"):
            rank += 10                      # radicals are never sensible surrogates
        out.append((rank, name, sp))
    if qc and not any(r <= 2 for r, _n, _s in out):   # no close match: same O and N, same saturation, nearest chain
        def dbe(c, h, n):                   # double-bond equivalents: an alkane stays an alkane
            return c - h / 2.0 + n / 2.0 + 1
        q_dbe = dbe(qc.get("n_C", 0), qc.get("n_H", 0), qc.get("n_N", 0))
        cands = [((abs(dbe(sp.n_C, sp.n_H, sp.n_N) - q_dbe), abs(sp.n_C - qc.get("n_C", 0))), name, sp)
                 for name, sp in MCM_DB.items()
                 if sp.n_O == qc.get("n_O", 0) and sp.n_N == qc.get("n_N", 0) and not sp.is_peroxy
                 and not name.endswith("O2") and not name.endswith("O")]
        if cands:
            best = min(c[0] for c in cands)
            out += [(4, name, sp) for d, name, sp in cands if d == best]
    out.sort(key=lambda t: (t[0], t[1]))
    return [(sp, _why(rank, sp, qu, qn)) for rank, name, sp in out[:limit]]


def _why(rank: int, sp: "MCMSpecies", qu: str, qn: str) -> str:
    base = {0: "same formula", 1: "same C and O count", 2: "same C count, O ±1", 3: "name match",
            4: "nearest carbon number with the same O"}.get(rank % 5 if rank < 10 else (rank - 10) % 5, "")
    tags = []
    if 5 <= rank < 10 or rank >= 15:
        tags.append("oxidation product")
    if rank >= 10:
        tags.append("radical")
    return base + (f" ({', '.join(tags)})" if tags else "")
