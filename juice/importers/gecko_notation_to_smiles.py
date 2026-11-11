"""Translate GECKO-A's structural notation (dictionary.out formula column) into SMILES.

GECKO-A grammar (from the official GitLab wiki, home/input-files):
  - linear chain carbon: C / CH / CH2 / CH3 (explicit H count, 0 if bare 'C')
  - double-bonded chain carbon: 'Cd' (both atoms of a C=C are 'Cd'; the '=' is written
    immediately before the double-bond partner)
  - aromatic ring carbon: lowercase 'c' / 'cH'
  - ring closure: a single digit right after an atom (before its H-count) opens/closes a ring
  - branch: '(...)' holding either a sub-chain or a terminal functional group:
    (OH) (OOH) (NO2) (ONO2) (OONO2) (O.) (OO.)
  - 'CHO' = terminal aldehyde carbon, 'CO' = carbonyl carbon (both carry an implicit =O)
  - '-O-' = inline ether oxygen in the main chain

Only the subset actually needed so far is implemented robustly (linear chains, CHO/CO,
the six functional-group branches, radicals). Rings/Cd/aromatic/ether raise NotImplementedError
so a bad parse is never silently wrong -- extend when the next molecule needs it.

Self-check built in: for every heavy atom, GECKO-A's explicit H-count in the text must equal
RDKit's standard-valence H-count once the skeleton (bonds, radical electrons) is built. A
mismatch raises instead of returning a wrong structure.
"""
from __future__ import annotations
import re
from rdkit import Chem
from rdkit.Chem import RWMol, Atom, BondType


class GeckoParseError(Exception):
    pass


# (RDKit SMARTS-free group builder) name -> list of atoms to append with (element, bonds_from_prev, extra_info)
# each group is built as its own little chain attached to the owner atom by a single bond.
_GROUP_RE = re.compile(r'^(OONO2|ONO2|OOH|OH|OO\.|O\.|NO2)')

_MAIN_ATOM_RE = re.compile(r'^(CHO|CO|C)(d)?(H([23])?)?')


def _add_functional_group(mol: RWMol, owner_idx: int, name: str) -> None:
    """Attach one of the 6 terminal functional groups to `owner_idx` by a single bond."""
    if name == 'OH':
        o = mol.AddAtom(Atom(8)); mol.GetAtomWithIdx(o).SetNoImplicit(False)
        mol.AddBond(owner_idx, o, BondType.SINGLE)
    elif name == 'OOH':
        o1 = mol.AddAtom(Atom(8)); o2 = mol.AddAtom(Atom(8))
        mol.AddBond(owner_idx, o1, BondType.SINGLE)
        mol.AddBond(o1, o2, BondType.SINGLE)
    elif name == 'O.':
        o = mol.AddAtom(Atom(8))
        mol.GetAtomWithIdx(o).SetNumRadicalElectrons(1)
        mol.AddBond(owner_idx, o, BondType.SINGLE)
    elif name == 'OO.':
        o1 = mol.AddAtom(Atom(8)); o2 = mol.AddAtom(Atom(8))
        mol.GetAtomWithIdx(o2).SetNumRadicalElectrons(1)
        mol.AddBond(owner_idx, o1, BondType.SINGLE)
        mol.AddBond(o1, o2, BondType.SINGLE)
    elif name == 'NO2':
        n = mol.AddAtom(Atom(7)); mol.GetAtomWithIdx(n).SetFormalCharge(1)
        o1 = mol.AddAtom(Atom(8)); o2 = mol.AddAtom(Atom(8))
        mol.GetAtomWithIdx(o2).SetFormalCharge(-1)
        mol.AddBond(owner_idx, n, BondType.SINGLE)
        mol.AddBond(n, o1, BondType.DOUBLE)
        mol.AddBond(n, o2, BondType.SINGLE)
    elif name == 'ONO2':
        o0 = mol.AddAtom(Atom(8))
        mol.AddBond(owner_idx, o0, BondType.SINGLE)
        _add_functional_group(mol, o0, 'NO2')  # reuse: bonds the N(O)(=O) onto the ester O
    elif name == 'OONO2':
        o0 = mol.AddAtom(Atom(8)); o1 = mol.AddAtom(Atom(8))
        mol.AddBond(owner_idx, o0, BondType.SINGLE)
        mol.AddBond(o0, o1, BondType.SINGLE)
        _add_functional_group(mol, o1, 'NO2')
    else:
        raise GeckoParseError(f"unknown functional group '{name}'")


def _parse_chain(mol: RWMol, s: str, pos: int, prev_idx: int | None,
                  checks: list) -> tuple[int, int | None]:
    """Parse a chain starting at s[pos:]; return (new_pos, last_main_chain_atom_idx).

    `checks` collects (idx, explicit_h, token) for every main-chain atom; validation is
    deferred to after the WHOLE molecule (all chain + branch bonds) has been built, since
    an atom's bond to its chain successor is only added once the next token is parsed."""
    last_idx = prev_idx
    n = len(s)
    while pos < n:
        c = s[pos]
        if c == ')':
            break
        if c in '-' and s[pos:pos + 3] == '-O-':
            raise NotImplementedError("inline ether '-O-' not yet supported")
        if c in 'cC':
            if c == 'c':
                raise NotImplementedError("aromatic 'c' not yet supported")
            m = _MAIN_ATOM_RE.match(s[pos:])
            if not m:
                raise GeckoParseError(f"cannot parse atom at {pos}: ...{s[pos:pos+10]}...")
            token, is_d, hpart, hdigit = m.groups()
            pos += len(m.group(0))
            if is_d:
                raise NotImplementedError("'Cd' (double-bond chain carbon) not yet supported")
            # ring-closure digit (only relevant once rings are supported)
            if pos < n and s[pos].isdigit():
                raise NotImplementedError("ring closure not yet supported")
            explicit_h = 0 if hpart is None else (1 if hdigit is None else int(hdigit))

            idx = mol.AddAtom(Atom(6))
            if token in ('CHO', 'CO'):
                o = mol.AddAtom(Atom(8))
                mol.AddBond(idx, o, BondType.DOUBLE)
                if token == 'CHO':
                    explicit_h = 1  # the aldehyde H (not written as 'H' in the text)
            if last_idx is not None:
                mol.AddBond(last_idx, idx, BondType.SINGLE)

            # consume branches immediately following this atom
            while pos < n and s[pos] == '(':
                close = _find_matching_paren(s, pos)
                inner = s[pos + 1:close]
                mg = _GROUP_RE.match(inner)
                if mg and mg.group(0) == inner:
                    _add_functional_group(mol, idx, inner)
                else:
                    # sub-chain branch, e.g. '(CH3)' or '(CH2CH3)'
                    _, _ = _parse_chain(mol, inner, 0, idx, checks)
                pos = close + 1

            checks.append((idx, explicit_h, token))
            last_idx = idx
            continue
        raise GeckoParseError(f"unexpected character '{c}' at {pos} in '{s}'")
    return pos, last_idx


def _find_matching_paren(s: str, open_pos: int) -> int:
    depth = 0
    for i in range(open_pos, len(s)):
        if s[i] == '(':
            depth += 1
        elif s[i] == ')':
            depth -= 1
            if depth == 0:
                return i
    raise GeckoParseError(f"unbalanced parenthesis in '{s}'")


def _check_valence(mol: RWMol, idx: int, explicit_h: int, token: str) -> None:
    """GECKO-A text gives an explicit H count; verify it matches what the skeleton implies."""
    atom = mol.GetAtomWithIdx(idx)
    bond_order_sum = sum(b.GetBondTypeAsDouble() for b in atom.GetBonds())
    radical = atom.GetNumRadicalElectrons()
    implied_h = 4 - int(bond_order_sum) - radical
    if implied_h != explicit_h:
        raise GeckoParseError(
            f"H-count mismatch for '{token}' (idx {idx}): text says {explicit_h}H, "
            f"skeleton implies {implied_h}H (bonds={bond_order_sum}, radical={radical})")
    atom.SetNoImplicit(True)
    atom.SetNumExplicitHs(explicit_h)


def gecko_to_mol(formula: str) -> Chem.Mol:
    formula = formula.strip()
    mol = RWMol()
    checks: list = []
    pos, last = _parse_chain(mol, formula, 0, None, checks)
    if pos != len(formula):
        raise GeckoParseError(f"trailing unparsed text: '{formula[pos:]}'")
    for idx, explicit_h, token in checks:
        _check_valence(mol, idx, explicit_h, token)
    m = mol.GetMol()
    Chem.SanitizeMol(m)
    return m


def gecko_to_smiles(formula: str) -> str:
    return Chem.MolToSmiles(gecko_to_mol(formula))


def gecko_formula_counts(formula: str) -> dict:
    """nC/nH/nN/nO/... from the parsed structure (for cross-check against dictionary.out)."""
    m = gecko_to_mol(formula)
    m = Chem.AddHs(m)
    counts: dict[str, int] = {}
    for atom in m.GetAtoms():
        counts[atom.GetSymbol()] = counts.get(atom.GetSymbol(), 0) + 1
    return counts


if __name__ == '__main__':
    tests = [
        ("CH3CH2(OH)", "ethanol"),
        ("CH3CHO", "acetaldehyde"),
        ("CH3CO(OH)", "acetic acid"),
        ("CH3CO(OOH)", "peracetic acid"),
        ("CH3CO(OONO2)", "PAN"),
        ("CH3(ONO2)", "methyl nitrate"),
        ("CH2(OH)CH2(OO.)", "hydroxyethyl peroxy"),
        ("CHO(OH)", "formic acid"),
        ("CHOCHO", "glyoxal"),
        ("CH3(OO.)", "methylperoxy"),
    ]
    for formula, name in tests:
        try:
            smi = gecko_to_smiles(formula)
            print(f"{name:22s} {formula:20s} -> {smi}")
        except Exception as e:
            print(f"{name:22s} {formula:20s} -> FAILED: {e}")
