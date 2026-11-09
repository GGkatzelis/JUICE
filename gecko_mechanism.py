"""Read a GECKO-A generator output folder (OUT/<mech>/): species dictionary, Nannoolal vapour
pressures, RO2 counter classes and the reaction list with its kinetics keywords.

Pure data layer (no Qt, no EASY) -- the EASY conversion is gecko_to_easy.py.

Kinetics conventions (verified in boxmodel4gecko LIBSRC/ratetool.f90, 2026-11-07):
  plain line  A n Ea : k = A * T^n * exp(-Ea/T)            (T in K -- NOT T/300)
  +TBODY             : k *= M
  +OXYGEN            : k *= 0.2*M
  +MEPERO            : k *= [CH3O2]              (counter, not consumed)
  +PEROj             : k *= sum(RO2 in class j)  (counter, not consumed; pero<j>.dat lists the class)
  +FALLOFF / k0A k0n k0Ea Fc flag . . /
                     : kinf = A*(T/300)^n*exp(-Ea/T), k0M = k0A*(T/300)^k0n*exp(-k0Ea/T)*M
                       flag 0 -> k = k0M/(1+k0M/kinf) * Fc^(1/(1+log10(k0M/kinf)^2))   (standard Troe)
                       flag 1 -> MCM broadening N = 0.75-1.27*log10(Fc)
  +EXTRA / label c2 c3 c4 c5 c6 c7 /
                     : 100 k*0.2*M*M | 500 k*H2O | 501 k*H2O + c2*exp(-c4/T)*H2O*M
                       550 k + k3*M/(1+k3*M/k2) with k2 = c2*T^c3*exp(-c4/T), k3 = c5*T^c6*exp(-c7/T)
  +HV / label hvfact / : photolysis J(label)*hvfact (J table = jfile.phot, see gecko_photolysis.py)
  AIN/AOU/WIN/WOU    : gas<->particle/wall transfer pseudo-reactions (reactions.dum only) -> kind 'transfer'
Names in the reaction/pvap/pero files carry a phase prefix (G gas, A particle, W wall); the bare
dictionary code is the name without that first character.
"""
from __future__ import annotations
import os
import re
from dataclasses import dataclass, field

KEYWORDS = {'HV', 'TBODY', 'OXYGEN', 'MEPERO', 'FALLOFF', 'EXTRA', 'ISOM', 'AIN', 'AOU', 'WIN', 'WOU'}
_PERO = re.compile(r'^PERO(\d)$')
_AUX = re.compile(r'^\s*(HV|FALLOFF|EXTRA|ISOM|WIN|WOU)\s*/\s*(.*?)\s*/\s*$')
# ISOM (ratetool.f90 'isomerisation reaction'): k = A*T^n*exp(-Ea/T) * (c1*T^4 + c2*T^3 + c3*T^2 + c4*T + c5),
# the five c_i on the `ISOM / c1 c2 c3 c4 c5 /` line (Vereecken-type T-polynomial correction).
_TERM = re.compile(r'^\s*((?:[0-9]+\.?[0-9]*|\.[0-9]+)(?:[eE][-+]?[0-9]+)?)?\s*([A-Za-z][A-Za-z0-9_]*)\s*$')
_ELEMENT_ORDER = ('C', 'H', 'N', 'O')


class GeckoReadError(Exception):
    pass


@dataclass
class GSpecies:
    code: str          # bare dictionary code: 'O02000', 'HO', 'CH2O'
    notation: str      # GECKO structural notation, or the fixed name for inorganics/C1
    mw: float          # g/mol (0.0 for the fixed inorganic names)
    radical: bool
    nC: int
    nH: int
    nN: int
    nO: int

    @property
    def formula(self) -> str:
        counts = {'C': self.nC, 'H': self.nH, 'N': self.nN, 'O': self.nO}
        out = ''
        for el in _ELEMENT_ORDER:
            n = counts[el]
            if n == 1:
                out += el
            elif n > 1:
                out += f"{el}{n}"
        return out

    @property
    def is_peroxy(self) -> bool:
        return '(OO.)' in self.notation


@dataclass
class GReaction:
    reactants: list          # [(coef, code)]
    products: list           # [(coef, code)]
    A: float
    n: float
    Ea: float
    kind: str = 'arrh'       # arrh|tbody|oxygen|mepero|pero|falloff|extra|hv|transfer
    aux: list = field(default_factory=list)
    pero_class: int = 0
    comment: str = ''
    line_no: int = 0

    @property
    def species(self) -> set:
        return {c for _, c in self.reactants} | {c for _, c in self.products}


@dataclass
class GeckoMechanism:
    folder: str
    species: dict            # code -> GSpecies
    reactions: list          # [GReaction]
    pvap: dict               # code -> (Psat_atm, dHvap_kJ_mol)
    tref: float
    pero: dict               # class -> [codes]
    primaries: list          # [codes]

    def rxn_species(self) -> set:
        s = set()
        for r in self.reactions:
            s |= r.species
        return s


def strip_phase(name: str) -> str:
    if name and name[0] in 'GAW':
        return name[1:]
    raise GeckoReadError(f"species name without phase prefix: '{name}'")


def read_dictionary(path: str) -> dict:
    species = {}
    with open(path, encoding='latin-1') as fh:
        for ln in fh:
            s = ln.strip()
            if not s or s.startswith('!') or s.upper().startswith('END'):
                continue
            toks = s.split()
            if len(toks) < 12:
                continue                      # the record-count header line
            try:
                mw = float(toks[-10])
                ints = [int(x) for x in toks[-9:]]
            except ValueError:
                continue
            code, notation = toks[0], toks[1]
            species[code] = GSpecies(code, notation, mw, ints[0] == 1, ints[1], ints[2], ints[3], ints[4])
    if not species:
        raise GeckoReadError(f"no species read from {path}")
    return species


def read_pvap(path: str) -> tuple:
    """-> (tref_K, {code: (Psat_atm, dHvap_kJ_mol)}); empty dict if the file is missing."""
    if not os.path.isfile(path):
        return 298.0, {}
    tref, pv = 298.0, {}
    with open(path, encoding='latin-1') as fh:
        for ln in fh:
            s = ln.strip()
            if not s or s.startswith('!') or s.upper().startswith('END'):
                continue
            toks = s.split()
            if toks[0].upper() == 'TREF':
                tref = float(toks[1]); continue
            if len(toks) < 3 or toks[1] == '!':
                continue                      # record-count line
            try:
                pv[strip_phase(toks[0])] = (float(toks[1]), float(toks[2]))
            except ValueError:
                continue
    return tref, pv


def read_pero(folder: str) -> dict:
    classes = {}
    for j in range(1, 10):
        p = os.path.join(folder, f"pero{j}.dat")
        if not os.path.isfile(p):
            continue
        names = []
        with open(p, encoding='latin-1') as fh:
            first = True
            for ln in fh:
                s = ln.strip()
                if first:
                    first = False; continue   # record count
                if not s or s.upper().startswith('END'):
                    continue
                names.append(strip_phase(s.split()[0]))
        classes[j] = names
    return classes


def read_primaries(folder: str) -> list:
    p = os.path.join(folder, "listprimary.dat")
    if not os.path.isfile(p):
        return []
    out = []
    with open(p, encoding='latin-1') as fh:
        for ln in fh:
            toks = ln.split()
            if toks:
                out.append(toks[0])
    return out


def _parse_side(side: str, is_lhs: bool, rxn: GReaction, line_no: int):
    terms = [t.strip() for t in side.split('+') if t.strip()]
    for t in terms:
        if is_lhs:
            up = t.upper()
            m = _PERO.match(up)
            if m or up in KEYWORDS:
                kw = 'pero' if m else up.lower()
                if up in ('AIN', 'AOU', 'WIN', 'WOU'):
                    kw = 'transfer'
                if rxn.kind != 'arrh':
                    raise GeckoReadError(f"line {line_no}: two kinetics keywords ({rxn.kind}, {kw})")
                rxn.kind = kw
                if m:
                    rxn.pero_class = int(m.group(1))
                continue
            rxn.reactants.append((1.0, strip_phase(t)))
        else:
            if t.upper() == 'NOTHING':
                continue
            m = _TERM.match(t)
            if not m:
                raise GeckoReadError(f"line {line_no}: cannot parse product term '{t}'")
            coef = float(m.group(1)) if m.group(1) else 1.0
            rxn.products.append((coef, strip_phase(m.group(2))))


def read_reactions(path: str) -> list:
    rxns = []
    cur = None
    with open(path, encoding='latin-1') as fh:
        raw_lines = fh.readlines()
    # a reaction whose product list ends with '+' continues on the next line
    lines, carry, carry_no = [], '', 0
    for no, raw in enumerate(raw_lines, 1):
        body = raw.rstrip('\r\n')
        if carry:
            body = carry + ' ' + body.strip()
            no = carry_no
            carry = ''
        if body.partition(';')[0].rstrip().endswith('+'):
            carry, carry_no = body.partition(';')[0].rstrip(), no
            continue
        lines.append((no, body))
    for i, line in lines:
        s = line.strip()
        if not s or s.startswith('!'):
            continue
        if s.upper() == 'END':
            break
        if s.upper() in ('REACTIONS', 'SPECIES', 'END SPECIES'):
            continue
        if '=>' not in line:
            m = _AUX.match(line)
            if not m:
                raise GeckoReadError(f"line {i}: unexpected line '{s}'")
            if cur is None:
                raise GeckoReadError(f"line {i}: auxiliary line before any reaction")
            kw, inner = m.group(1), m.group(2)
            vals = [float(x) for x in inner.split()]
            if kw == 'HV':
                cur.aux = [int(round(vals[0])), vals[1] if len(vals) > 1 else 1.0]
            elif kw in ('FALLOFF', 'EXTRA', 'ISOM'):
                cur.aux = vals
            continue
        body, _, comment = line.partition(';')
        toks = body.split()
        try:
            A, n, Ea = float(toks[-3]), float(toks[-2]), float(toks[-1])
        except (ValueError, IndexError):
            raise GeckoReadError(f"line {i}: no A n Ea at end of '{s}'")
        eq = ' '.join(toks[:-3])
        lhs, rhs = eq.split('=>', 1)
        cur = GReaction([], [], A, n, Ea, comment=' '.join(comment.split()), line_no=i)
        _parse_side(lhs, True, cur, i)
        _parse_side(rhs, False, cur, i)
        rxns.append(cur)
    return rxns


def read_mechanism(folder: str, reactions_file: str | None = None) -> GeckoMechanism:
    """reactions_file defaults to reactionswithcom.dum (kinetics provenance tags as comments, no
    transfer pseudo-reactions); reactions.dum works too (transfer reactions come back as kind 'transfer')."""
    folder = os.path.abspath(folder)
    species = read_dictionary(os.path.join(folder, "dictionary.out"))
    if reactions_file is None:
        reactions_file = os.path.join(folder, "reactionswithcom.dum")
        if not os.path.isfile(reactions_file):
            reactions_file = os.path.join(folder, "reactions.dum")
    rxns = read_reactions(reactions_file)
    tref, pv = read_pvap(os.path.join(folder, "pvap.nan.dat"))
    pero = read_pero(folder)
    mech = GeckoMechanism(folder, species, rxns, pv, tref, pero, read_primaries(folder))
    unknown = mech.rxn_species() - set(species)
    if unknown:
        raise GeckoReadError(f"reaction species missing from dictionary.out: {sorted(unknown)}")
    return mech


def is_radical(sp: GSpecies) -> bool:
    return sp.radical or '.' in sp.notation


def is_organic(sp: GSpecies) -> bool:
    return sp.nC > 0 and sp.code not in ('CO', 'CO2')


def generations(mech: GeckoMechanism) -> dict:
    """Static oxidation-generation estimate (code -> generation) by breadth-first search from the primaries:
    a closed-shell organic reactant starts the next generation, a radical reactant stays in its own;
    gen(product) = minimum over all producing reactions. Primaries = 0, their first stable products = 1."""
    from collections import deque
    gen = {c: 0 for c in mech.primaries}
    by_reactant = {}
    for r in mech.reactions:
        if r.kind == 'transfer':
            continue
        for _, c in r.reactants:
            sp = mech.species.get(c)
            if sp is not None and is_organic(sp):
                by_reactant.setdefault(c, []).append(r)
    q = deque(mech.primaries)
    while q:
        c = q.popleft()
        g = gen[c] + (0 if is_radical(mech.species[c]) else 1)
        for r in by_reactant.get(c, []):
            for _, p in r.products:
                sp = mech.species.get(p)
                if sp is None or not is_organic(sp):
                    continue
                if p not in gen or g < gen[p]:
                    gen[p] = g
                    q.append(p)
    return gen


def generation_cut(mech: GeckoMechanism, maxgen: int):
    """Reduced mechanism in which closed-shell organics of generation > maxgen do not react any further
    (they stay as terminal products, like GECKO-A's own maxgen); radicals always react. -> (mech, gen)."""
    gen = generations(mech)

    def may_react(c: str) -> bool:
        sp = mech.species[c]
        if not is_organic(sp):
            return True
        if c not in gen:
            return False                       # unreachable from the primaries
        return is_radical(sp) or gen[c] <= maxgen

    # forward closure from the primaries: a reaction is kept once all its organic reactants are reachable
    # and allowed to react; its products become reachable. Inorganic-only reactions are always kept.
    from collections import deque
    by_reactant, rxns, kept_ids = {}, [], set()
    for r in mech.reactions:
        orgs = [c for _, c in r.reactants if c in mech.species and is_organic(mech.species[c])]
        if not orgs:
            rxns.append(r)
            kept_ids.add(id(r))
        for c in orgs:
            by_reactant.setdefault(c, []).append(r)
    reach = set(mech.primaries)
    q = deque(mech.primaries)
    while q:
        c = q.popleft()
        if not may_react(c):
            continue
        for r in by_reactant.get(c, []):
            if id(r) in kept_ids:
                continue
            orgs = [x for _, x in r.reactants if x in mech.species and is_organic(mech.species[x])]
            if any(x not in reach or not may_react(x) for x in orgs):
                continue
            rxns.append(r)
            kept_ids.add(id(r))
            for _, p in r.products:
                sp = mech.species.get(p)
                if sp is not None and is_organic(sp) and p not in reach:
                    reach.add(p)
                    q.append(p)
    rxns.sort(key=lambda r: r.line_no)
    used = reach | {c for c, s in mech.species.items() if not is_organic(s)}
    species = {c: s for c, s in mech.species.items() if c in used}
    pvap = {c: v for c, v in mech.pvap.items() if c in species}
    pero = {k: [c for c in v if c in species] for k, v in mech.pero.items()}
    return GeckoMechanism(mech.folder, species, rxns, pvap, mech.tref, pero, list(mech.primaries)), gen


if __name__ == '__main__':
    import sys
    from collections import Counter
    m = read_mechanism(sys.argv[1] if len(sys.argv) > 1 else r"Z:\TEMP\g.gkatzelis\Franz\GECKO-A\ethanol")
    print(f"{m.folder}: {len(m.species)} species, {len(m.reactions)} reactions, "
          f"{len(m.pvap)} with Nannoolal Psat @ {m.tref} K, primaries {m.primaries}")
    print("kinds:", dict(Counter(r.kind for r in m.reactions)))
    print("RO2 classes:", {k: v for k, v in m.pero.items() if v})
    print("HV labels:", sorted({r.aux[0] for r in m.reactions if r.kind == 'hv'}))
    print("EXTRA labels:", sorted({int(r.aux[0]) for r in m.reactions if r.kind == 'extra'}))
    fo_flags = {tuple(r.aux[4:]) for r in m.reactions if r.kind == 'falloff'}
    print("FALLOFF flag/extra columns:", fo_flags)
