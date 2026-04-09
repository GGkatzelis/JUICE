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

import os
import re
from dataclasses import dataclass, field
from typing import Optional

# ---------------------------------------------------------------------------
# Path to the TSV (same directory as this file)
# ---------------------------------------------------------------------------
_DB_DIR = os.path.join(os.path.dirname(__file__), "MCM species formula")
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
