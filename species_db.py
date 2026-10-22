"""species_db.py -- ONE loader for ALL species TSVs; the single source of truth for
formula / atom counts / molar mass / RO2-flag / SMILES across the whole project.

RULE (user, 2026-10-05): "whenever you need TSV data -- really everywhere -- take the
object from this subroutine." So every make_*/plot_* script that needs 'is X an RO2?',
'what is X's formula / nC / nO / nN / nH / MW?' must get it from here, not re-parse TSVs.

Usage:
    from species_db import load_species_db
    db = load_species_db()                 # cached singleton (reads every TSV once)
    db.is_ro2('APINAO2')                    # True / False (PeroxyRadical column)
    db.formula('C729CHO')                   # 'C7H12O2'  (TSV truth, NOT the MCM-code name)
    db.nC('SOA_PINAL'), db.mw('WALL_PINAL') # phase prefixes SOA_/WALL_/WALL_PL_ are stripped
    db.nO('P_C10H15O8')                     # systematic P_CxHyOz parsed if not tabulated
    sp = db.get('TOLUENE')                  # Species dataclass (or None)
    db.nc, db.no, db.nn, db.nh             # raw {NAME:count} dicts (tabulated species only)

Accessors (nC/nO/nN/nH/mw/formula/is_ro2) apply, in order:
  1. exact tabulated value (TSV, first file that defines the name wins -- setdefault);
  2. after stripping a SOA_/WALL_PL_/WALL_ phase prefix, the tabulated value of the base;
  3. a systematic formula name (optionally P_-prefixed), e.g. P_C10H15O8 / C11H14O10 -> parsed;
  4. else 0 (counts) / 0.0 (mw) / '' (formula) / False (is_ro2).
This reproduces the old per-script `_fromname`/`ro2form` logic, in one place.
"""
from __future__ import annotations
import os, re, sys
from dataclasses import dataclass, field

# Model internals / diagnostic columns that are legitimately NOT species in a TSV (so they must not
# trigger an "unknown molecule" halt). Add to this list rather than to a TSV when a genuine non-species
# name surfaces. Phase-accumulator group sums (SOA_SVOC_C ...) are matched by _GRPSUM below.
NON_MOLECULE = frozenset({
    'M', 'P', 'T', 'TEMP', 'SZ', 'SZA', 'SZANGLE', 'HV', 'RH', 'FLOW', 'TIME', 'DUMMY', 'J',
    'JNO2', 'JO1D', 'JHONO', 'JH2O2', 'JHCHO_R', 'JHCHO_M',
    'DIL', 'DILUTE', 'RO2', 'RO2T', 'COA', 'CW', 'WALLDEP',
    'QSEED', 'Q_SEED', 'CW_SCALE', 'KCOND', 'KFWALL', 'KPWALL', 'KROOR_PRAM', 'KOLIG', 'SCALE_PAN',
    # inorganic/alias columns absent from the carbon species TSVs, and SAPHIR env diagnostics:
    'N2', 'NA', 'HNO4', 'DV', 'DVO3', 'COSX', 'XNOY', 'LNO2', 'MNO2', 'NNO2', 'QL', 'QY',
})
_GRPSUM = re.compile(r'^(SOA|WALL|WALL_PL)_(SVOC|LVOC|ELVOC)_[CON]$')
# Systematic model-internal column classes that are never species (rate constants, photolysis J-values,
# PRAM/yield knobs, Q injection sources, *_IN forcing inputs, bare group element sums). Centrally
# maintained so a non-species column never triggers the unknown-molecule halt.
_NONMOL_RE = re.compile(
    r'^(J\d+(_CS)?'
    r'|K[A-Z0-9_]*'                                   # ALL rate coefficients (KMT.., KRO2.., KAUTOOX.., KQ_n, KR_n, KAT_NC12_.., KBPPN, KUNI_.., KROOR_.., K298.., K14..); MCM/SIM-HOM has no species starting with K
    r'|Y(HOM|ACC|O3|OH|NO3|O2)[A-Z0-9_]*|Y[A-Z0-9_]*_(PRAM|LIMO)'   # HOM/accretion yield knobs
    r'|Q[A-Z0-9_]*|INJECTION_VOLUME_[A-Z0-9]+'        # injection sources / volumes
    r'|S(H[A-Z0-9]*|NO2|O3|HCHO|HONO)'                # SAPHIR S-sources (SHONO, SHCHO, SNO2, SO3...)
    r'|X_[A-Z0-9]+|RES\d*'                            # inert injection tracers (X_limonene), reservoir diagnostics
    r'|(SVOC|LVOC|ELVOC)_[CON]'                       # bare group element sums
    r'|[A-Z][A-Z0-9]*_IN)$')                          # forcing INPUT columns (H2O_IN, T_IN, HV_IN)

# ---------------------------------------------------------------------------
# Canonical TSV list. First file to define a NAME wins (setdefault). Order =
# most-authoritative first. Missing files are skipped silently. Add paths via
# load_species_db(extra_paths=[...]) rather than editing call sites.
# ---------------------------------------------------------------------------
_MCMd = "O:/IDL_drive_O/fr_lib/easy/Mechanismen/MCM"
_SH   = "E:/Python Code/JUICE archive/EASY/include files/SIM_HOM"
_MCM331 = "E:/Python Code/JUICE archive/EASY/include files/MCM331/mcm_3-3-1_windows"
# Each source is (path, mode). 'clean' = well-formed tab+header TSV (Name..Formula..PeroxyRadical,
# constant column count) -> authoritative for formula AND the RO2 flag. 'positional' = messy/IDL
# MCM-export file with empty padding cells and inconsistent tabs -> read formula from the 5th
# WHITESPACE-collapsed token only (as the historical budget-plot loader did); it carries no usable
# RO2 column, so RO2 is never taken from it (that is what caused NO2/CO2/H2O2 to be mis-flagged).
_HERE = os.path.dirname(os.path.abspath(__file__))
_REG = os.path.join(_HERE, "mechanisms")


def _registry_sources():
    """Species tables shipped with the repository (mechanism registry): every machine has these.
    base/<ID>/*.tsv are MCM exports (positional), addons/<ID>/*.tsv are clean tables."""
    out = []
    for kind, mode in (("base", 'positional'), ("addons", 'clean')):
        d = os.path.join(_REG, kind)
        if not os.path.isdir(d):
            continue
        for mid in sorted(os.listdir(d)):
            md = os.path.join(d, mid)
            if not os.path.isdir(md):
                continue
            for fn in sorted(os.listdir(md), key=lambda f: f != "mcm_export_species.tsv"):   # the clean export first
                if fn.lower().endswith(".tsv"):
                    out.append((os.path.join(md, fn), 'clean' if fn == "mcm_export_species.tsv" else mode))
    return out


DEFAULT_SOURCES = [
    # clean, authoritative (formula + RO2) -- first, so HOM/isoprene formulas are the radical formulas:
    (_MCMd + "/CHANEL.all_mixes_gas_species.tsv",        'clean'),
    (_MCMd + "/SIMHOM2026_gas_species.FULL.tsv",         'clean'),
    (_MCMd + "/SIMHOM2026_HOMprototype_gas_species.tsv", 'clean'),
    # positional MCM base/extra (formula only; fills gaps the clean files don't cover):
    (_MCMd + "/mcm_export_species.tsv",                  'positional'),
    (_MCMd + "/additional_species.tsv",                  'positional'),
    (_MCM331 + "/mcm_3-3-1_species_complete.tsv",        'positional'),
    (_SH + "/guaiacol_creosol_gas.tsv",                  'positional'),
    # condenser formulas used by make_soa_condensation to build the SOA include group-sums -- fill
    # the closed-shell condensers (incl. LIMO/BPIN HOM) the clean files don't carry (formula @ col 1):
    (_SH + "/SIMHOM_volatility.param",                   'volparam'),
    # last-resort fill (messy padding -> only trusted for gaps; clean files already won the radicals):
    (_SH + "/SIMHOM2026_gas_species.tsv",               'positional'),
] + _registry_sources()   # mechanisms/base + mechanisms/addons of this checkout (the O:/E: drives are Franz's)

_AW = {'C': 12.011, 'H': 1.008, 'O': 15.999, 'N': 14.007, 'Cl': 35.453,
       'Br': 79.904, 'S': 32.06, 'F': 18.998}
_ATOM_RE = re.compile(r'(Cl|Br|Si|[CHNOFSP])(\d*)')
_PHASE_PREFIXES = ('SOA_', 'WALL_PL_', 'WALL_')
_SYS_RE = re.compile(r'^C\d+H\d+(?:O\d*)?(?:N\d*)?$')   # systematic CxHyOzNw name (after optional P_)


def parse_formula(formula: str) -> dict:
    """{'C':n,'H':n,...} from a molecular-formula string ('' / '-' -> {})."""
    d = {}
    if not formula or formula in ('-', 'X', ''):
        return d
    for atom, num in _ATOM_RE.findall(formula):
        if atom in _AW:
            d[atom] = d.get(atom, 0) + (int(num) if num else 1)
    return d


@dataclass
class Species:
    name: str
    formula: str = ''
    smiles: str = ''
    mass: float = 0.0          # TSV Mass column (0 if absent)
    is_ro2: bool = False       # PeroxyRadical == true
    source: str = ''           # basename of the TSV it came from
    nC: int = 0; nH: int = 0; nO: int = 0; nN: int = 0
    # --- aliases for drop-in compatibility with the old mcm_species_db.MCMSpecies ---
    @property
    def n_C(self): return self.nC
    @property
    def n_H(self): return self.nH
    @property
    def n_O(self): return self.nO
    @property
    def n_N(self): return self.nN
    @property
    def is_peroxy(self): return self.is_ro2


def _merge(db, name, formula, mass, is_ro2, smiles, src):
    existing = db.get(name)
    if existing is None:                                            # first definition: formula is authoritative
        d = parse_formula(formula)
        db[name] = Species(name=name, formula=(formula if formula not in ('', '-') else ''),
                           smiles=(smiles if smiles not in ('', '-') else ''), mass=mass, is_ro2=is_ro2, source=src,
                           nC=d.get('C', 0), nH=d.get('H', 0), nO=d.get('O', 0), nN=d.get('N', 0))
    else:                                                          # merge: OR the RO2 flag, backfill a missing formula
        if is_ro2: existing.is_ro2 = True
        if not existing.formula and formula not in ('', '-'):
            d = parse_formula(formula)
            existing.formula = formula; existing.source = src
            existing.nC, existing.nH, existing.nO, existing.nN = d.get('C', 0), d.get('H', 0), d.get('O', 0), d.get('N', 0)
        if not existing.mass and mass: existing.mass = mass
        if not existing.smiles and smiles not in ('', '-'): existing.smiles = smiles


def _load_one(path: str, mode: str, db: dict, quiet: bool):
    if not os.path.isfile(path):
        if not quiet: print("  species_db: missing", path)
        return 0
    src = os.path.basename(path); n0 = len(db)
    try:
        lines = open(path, encoding='latin-1', errors='replace').read().splitlines()
    except Exception as e:
        if not quiet: print("  species_db: cannot read", path, e)
        return 0
    if mode == 'clean':
        # tab-delimited with a 'Name .. Formula .. PeroxyRadical' header and a constant column count.
        iN = iF = iM = iP = iS = -1; have = False
        for ln in lines:
            if not ln.strip() or ln.lstrip().startswith((';', '*')):
                continue
            cells = ln.split('\t')
            low = [c.strip().lower() for c in cells]
            if not have:
                if 'name' in low and 'formula' in low:
                    iN = low.index('name'); iF = low.index('formula')
                    iM = low.index('mass') if 'mass' in low else -1
                    iP = next((i for i, c in enumerate(low) if 'peroxyradical' in c), -1)
                    iS = low.index('smiles') if 'smiles' in low else -1
                    have = True
                continue
            if len(cells) <= iF: continue
            name = cells[iN].strip().upper()
            if not name or name.startswith(('*', ';')): continue
            formula = cells[iF].strip()
            is_ro2  = (0 <= iP < len(cells) and cells[iP].strip().lower() == 'true')
            smiles  = cells[iS].strip() if (0 <= iS < len(cells)) else ''
            mass = 0.0
            if 0 <= iM < len(cells):
                try: mass = float(cells[iM].strip())
                except ValueError: mass = 0.0
            _merge(db, name, formula, mass, is_ro2, smiles, src)
    elif mode == 'volparam':
        # SIMHOM_volatility.param: tab-delimited, header 'name formula MW smiles nC nO ...'; formula @ col 1; no RO2.
        have = False
        for ln in lines:
            if not ln.strip() or ln.lstrip().startswith(';'): continue
            cells = ln.split('\t')
            low = [c.strip().lower() for c in cells]
            if not have:
                if 'name' in low and 'formula' in low: have = True
                continue
            if len(cells) < 2: continue
            name = cells[0].strip().upper()
            if not name: continue
            _merge(db, name, cells[1].strip(), 0.0, False, '', src)
    else:
        # positional: formula = 5th WHITESPACE-collapsed token (historical budget-plot behaviour); no RO2.
        for ln in lines:
            s = ln.strip()
            if not s or s.startswith((';', '*')): continue
            t = [x for x in re.split(r'[ \t]+', s) if x != '']
            if not t or t[0].lower() == 'name' or len(t) < 5: continue
            name = t[0].strip().upper()
            if name.startswith(('*', ';')): continue
            mass = 0.0
            if len(t) > 5:
                try: mass = float(t[5])
                except ValueError: mass = 0.0
            _merge(db, name, t[4], mass, False, (t[1] if len(t) > 1 else ''), src)
    return len(db) - n0


class UnknownSpecies(KeyError):
    """Raised by SpeciesDB.require() for a molecule name absent from every species TSV."""


class SpeciesDB:
    """Merged species table + accessors with phase-prefix strip and systematic-name fallback.
    Every accessor NOTES a name that resolves to nothing (genuine unknown molecule) in self.unknown;
    call check()/report_unknown() to halt+report so a missing species is fixed centrally (one TSV)."""
    def __init__(self, table: dict):
        self._t = table
        self.nc = {k: v.nC for k, v in table.items()}
        self.no = {k: v.nO for k, v in table.items()}
        self.nn = {k: v.nN for k, v in table.items()}
        self.nh = {k: v.nH for k, v in table.items()}
        self.unknown = set()                 # molecule names queried but found in no TSV (for report_unknown)

    def __len__(self): return len(self._t)
    def __contains__(self, name): return name.upper() in self._t
    def __getitem__(self, name): return self._t[name.upper()]

    @staticmethod
    def _strip(u):
        for p in _PHASE_PREFIXES:
            if u.startswith(p): return u[len(p):]
        return u

    def resolvable(self, name):
        """True unless `name` is an UNKNOWN molecule: known species (after phase strip), systematic
        CxHyOzNw name, model-internal (NON_MOLECULE) or a group-sum column all count as resolvable."""
        u = name.upper(); b = self._strip(u)
        if u in self._t or b in self._t: return True
        y = b[2:] if b.startswith('P_') else b
        if _SYS_RE.match(y): return True
        if not re.fullmatch(r'[A-Za-z0-9_]+', b): return True          # labels/headers (spaces, '/') are not molecules
        if u in NON_MOLECULE or b in NON_MOLECULE: return True
        if _GRPSUM.match(u) or _GRPSUM.match(b): return True
        if _NONMOL_RE.match(u) or _NONMOL_RE.match(b): return True     # rate constants / J-values / Q-sources / *_IN / group sums
        return False

    def _note(self, name):
        if not self.resolvable(name): self.unknown.add(name)

    def get(self, name):
        """Species for name (phase prefix stripped); None if untabulated."""
        u = name.upper()
        return self._t.get(u) or self._t.get(self._strip(u))

    def _count(self, name, el):
        """el in 'CHON'. Tabulated -> prefix-stripped tabulated -> systematic-name parse -> 0 (noted)."""
        self._note(name)
        u = name.upper()
        mp = {'C': self.nc, 'O': self.no, 'N': self.nn, 'H': self.nh}[el]
        if mp.get(u, 0) > 0: return mp[u]
        b = self._strip(u)
        if b != u and mp.get(b, 0) > 0: return mp[b]
        y = b[2:] if b.startswith('P_') else b
        if _SYS_RE.match(y):
            m = re.search(el + r'(\d+)', y)
            if m: return int(m.group(1))
            return 1 if el == 'C' else 0
        return 0

    def nC(self, name): return self._count(name, 'C')
    def nH(self, name): return self._count(name, 'H')
    def nO(self, name): return self._count(name, 'O')
    def nN(self, name): return self._count(name, 'N')

    def mw(self, name):
        """Molar mass from atom counts (nC/nH/nO/nN already note unknowns). 0.0 if unknown."""
        return (self.nC(name) * _AW['C'] + self.nH(name) * _AW['H']
                + self.nO(name) * _AW['O'] + self.nN(name) * _AW['N'])

    def formula(self, name):
        """TSV formula (prefix stripped); systematic names keep their own string; '' if unknown (noted)."""
        self._note(name)
        sp = self.get(name)
        if sp is not None and sp.formula not in ('', '-'): return sp.formula
        b = self._strip(name.upper())
        y = b[2:] if b.startswith('P_') else b
        return y if _SYS_RE.match(y) else ''

    def smiles(self, name):
        """SMILES from the TSV (prefix stripped); '' if none/unknown (noted)."""
        self._note(name)
        sp = self.get(name)
        return sp.smiles if (sp is not None and sp.smiles not in ('', '-')) else ''

    def mass(self, name):
        """TSV Mass column if present, else computed MW (noted)."""
        self._note(name)
        sp = self.get(name)
        if sp is not None and sp.mass > 0: return sp.mass
        return self.mw(name)

    def is_ro2(self, name):
        """PeroxyRadical flag (prefix stripped). The single truth for 'is X an RO2?' (noted)."""
        self._note(name)
        sp = self.get(name)
        return bool(sp.is_ro2) if sp is not None else False

    def ro2_names(self):
        return [k for k, v in self._t.items() if v.is_ro2]

    def names_from_source(self, src):
        """Names whose FIRST (authoritative) definition came from the TSV basename `src`, in load
        order. For a source listed first in DEFAULT_SOURCES (e.g. CHANEL.all_mixes_gas_species.tsv)
        this is exactly that file's species set, in its original row order -- a drop-in for 'iterate
        the species in that TSV' without re-reading the file."""
        s = os.path.basename(src)
        return [k for k, v in self._t.items() if v.source == s]

    def require(self, name):
        """Return the Species, or raise UnknownSpecies if `name` is an unknown molecule."""
        sp = self.get(name)
        if sp is None and not self.resolvable(name):
            raise UnknownSpecies(name)
        return sp

    def check(self, names, context='', halt=True):
        """Note every name in `names`, then report+halt if any are unknown molecules. Call this with
        a consumer's full species-name list up front to surface missing species centrally."""
        for n in names: self._note(n)
        return self.report_unknown(context=context, halt=halt)

    def report_unknown(self, context='', halt=True):
        """If any unknown molecule names were noted, print them (+ how to fix) and, if halt, SystemExit."""
        if not self.unknown: return False
        names = sorted(self.unknown)
        head = f"[species_db] {len(names)} UNKNOWN molecule name(s) not in any species TSV"
        if context: head += f"  (context: {context})"
        body = "  " + ", ".join(names[:60]) + ("" if len(names) <= 60 else f"  ... (+{len(names)-60} more)")
        hint = ("  -> add each to a species TSV (Name/Formula/Mass/.../PeroxyRadical), or, if it is a model\n"
                "     internal/diagnostic column, to species_db.NON_MOLECULE; then rerun (single source of truth).")
        print(head + "\n" + body + "\n" + hint, file=sys.stderr)
        if halt: raise SystemExit(2)
        return True


_CACHE = {}
def load_species_db(extra_paths=None, quiet=True, reload=False) -> SpeciesDB:
    """Read every source TSV once (cached) and return the shared SpeciesDB. extra_paths appends
    scenario-specific sources: each is a path (parsed 'clean') or a (path, mode) tuple.
    reload=True rebuilds the cache."""
    extras = [(p if isinstance(p, tuple) else (p, 'clean')) for p in (extra_paths or [])]
    key = tuple(extras)
    if reload: _CACHE.pop(key, None)
    if key in _CACHE: return _CACHE[key]
    table = {}
    for path, mode in DEFAULT_SOURCES + extras:
        _load_one(path, mode, table, quiet)
    if not quiet:
        nro2 = sum(1 for v in table.values() if v.is_ro2)
        print(f"species_db: {len(table)} species ({nro2} RO2) from {len(DEFAULT_SOURCES)+len(extras)} sources")
    db = SpeciesDB(table)
    _CACHE[key] = db
    return db


if __name__ == "__main__":
    db = load_species_db(quiet=False)
    for t in ["TOLUENE", "APINAO2", "C729CHO", "PINAL", "P_C10H15O8", "SOA_PINAL", "WALL_C2H5OH", "C2H5OH"]:
        print(f"  {t:14} formula={db.formula(t):12} nC={db.nC(t):2} nO={db.nO(t):2} nN={db.nN(t)} MW={db.mw(t):7.2f} RO2={db.is_ro2(t)}")
