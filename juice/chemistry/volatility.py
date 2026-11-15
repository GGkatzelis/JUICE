"""Volatility (C*) of the species of a mechanism, by a method the user chooses.

C* is the saturation mass concentration in ug m-3 at 298.15 K.  The gas <-> particle <-> wall scheme
(importers/make_soa_condensation.emit_soa) needs one C* per condensable species; where it comes from is a
modelling choice with orders of magnitude between the answers, so every method is one entry of METHODS and the
Chemistry step offers them in a drop-down.  A new method is a function (Molecule) -> C* and one METHODS line.

    structure based (need the SMILES)
      simpol        SIMPOL.1, Pankow and Asher (2008): group contributions                       RDKit
      nannoolal     Nannoolal et al. (2008) vapour pressure with the Nannoolal et al. (2004)
                    boiling point: what GECKO-A uses for its own species                         UManSysProp
      evaporation   EVAPORATION, Compernolle et al. (2011): aliphatic oxygenates, not aromatics  UManSysProp
    formula based (volatility basis set parameterisations; carbon, oxygen and nitrogen numbers only)
      li2016        Li et al. (2016), molecular corridors, CHO and CHON coefficients
      donahue2011   Donahue et al. (2011), 2-D volatility basis set
    mixed
      li_evap       Li et al. (2016) for CHO; organic nitrogen by EVAPORATION where the structure allows
                    (F. Rohrer's routing for organonitrates and PANs, which the formula methods make far too sticky)

SIMPOL coefficients and SMARTS, the Li/Donahue forms and the li_evap routing are taken from F. Rohrer's
make_simpol.py, compare_volatility_methods.py and b2_volatility.py (MCM_SAPHIR_Tool, branch easy-start-gui), which
transcribe SIMPOL from UManSysProp (Topping et al., 2016).  His comparison on 387 limonene products: SIMPOL,
EVAPORATION and Nannoolal agree within 0.1 orders of magnitude in the median.  The Li and Donahue coefficients
were checked against the papers on 2026-11-15 and differ from his scripts in two places: Donahue's bCO is -0.3
(his scripts: +0.3), Li's CHO bCO is -0.7790 (his scripts: -0.07790, the CHON value).  With the papers' values,
on the 571 closed-shell products of alpha-pinene + limonene in the MCM, the medians against SIMPOL are: Li +0.7
orders (CHO) and -2.4 (with N); Donahue -2.2 (CHO) and -7.5 (with N: it counts the nitrate oxygens and has no
nitrogen term, so organic nitrates and PANs become condensable).

UManSysProp and Open Babel are optional: `python -m pip install openbabel-wheel "setuptools<81"` (UManSysProp reads
its tables with pkg_resources, which setuptools 81 removed) and a copy of github.com/loftytopping/UManSysProp_public
(it has no installer): unpacked as GitHub/_external/UManSysProp_public next to JUICE, or anywhere with the
environment variable UMANSYS pointing to it.
"""
from __future__ import annotations

import math
import os
import re
import sys
from collections import OrderedDict
from dataclasses import dataclass
from typing import Callable, Dict, Iterable, List, Optional, Tuple

R_ATM = 8.206e-5          # m3 atm mol-1 K-1
T_REF = 298.15
NOT_PARTITIONED = {"CO", "CO2", "CH4"}
GROUP_LIMITS = (3e-4, 0.3, 300.0)        # ELVOC | LVOC | SVOC | (IVOC and VOC)


@dataclass
class Molecule:
    name: str
    nC: int
    nH: int
    nO: int
    nN: int
    mass: float
    smiles: str = ""


def volatility_class(cstar: float) -> str:
    return "ELVOC" if cstar < GROUP_LIMITS[0] else ("LVOC" if cstar < GROUP_LIMITS[1] else ("SVOC" if cstar < GROUP_LIMITS[2] else "VOC"))


# ── formula based ──────────────────────────────────────────────────────────────────────────────────────────────

def _corridor(nC, nO, n0, bC, bO, bCO):
    return 10 ** ((n0 - nC) * bC - nO * bO - 2 * (nC * nO / max(1, nC + nO)) * bCO)


def cstar_donahue(m: Molecule) -> float:
    """Donahue et al. (2011), ACP 11, 3303, Eq. (3): log10 C0 = (nC0 - nC) bC - nO bO - 2 nC nO/(nC + nO) bCO with
    nC0 = 25, bC = 0.475, bO = 2.3 and bCO = -0.3 (their Fig. 3), so the cross term raises C0."""
    return _corridor(m.nC, m.nO, 25.0, 0.475, 2.3, -0.3) if m.nC > 0 else 1e9


# Li et al. (2016), ACP 16, 3327, Eq. (2) and Table 1 (checked against the paper 2026-11-15):
#   log10 C0 = (n0C - nC) bC - nO bO - 2 nC nO/(nC + nO) bCO - nN bN - nS bS      class: (n0C, bC, bO, bCO, bN)
LI_2016 = {"CH": (23.80, 0.4861, 0.0, 0.0, 0.0), "CHO": (22.66, 0.4481, 1.656, -0.7790, 0.0),
           "CHN": (24.59, 0.4066, 0.0, 0.0, 0.9619), "CHON": (24.13, 0.3667, 0.7732, -0.07790, 1.114)}


def cstar_li(m: Molecule) -> float:
    """Li et al. (2016) with the coefficients of the compound's class (CH, CHO, CHN or CHON)."""
    if m.nC <= 0:
        return 1e9
    cls = "CH" + ("O" if m.nO > 0 else "") + ("N" if m.nN > 0 else "")
    n0, bC, bO, bCO, bN = LI_2016[cls]
    return _corridor(m.nC, m.nO, n0, bC, bO, bCO) * 10 ** (-m.nN * bN)


def _li_cho(nC, nO):
    return cstar_li(Molecule("", nC, 0, nO, 0, 0.0)) if nC > 0 else 1e9


# ── SIMPOL.1 (RDKit) ───────────────────────────────────────────────────────────────────────────────────────────

_SIMPOL_COEF = {
    0: (-4.26938E+02, 2.89223E-01, 4.42057E-03, 2.92846E-01), 1: (-4.11248E+02, 8.96919E-01, -2.48607E-03, 1.40312E-01),
    2: (-1.46442E+02, 1.54528E+00, 1.71021E-03, -2.78291E-01), 3: (3.50262E+01, -9.20839E-01, 2.24399E-03, -9.36300E-02),
    4: (-8.72770E+01, 1.78059E+00, -3.07187E-03, -1.04341E-01), 5: (5.73335E+00, 1.69764E-02, -6.28957E-04, 7.55434E-03),
    6: (-2.61268E+02, -7.63282E-01, -1.68213E-03, 2.89038E-01), 7: (-7.25373E+02, 8.26326E-01, 2.50957E-03, -2.32304E-01),
    8: (-7.29501E+02, 9.86017E-01, -2.92664E-03, 1.78077E-01), 9: (-1.37456E+01, 5.23486E-01, 5.50298E-04, -2.76950E-01),
    10: (-7.98796E+02, -1.09436E+00, 5.24132E-03, -2.28040E-01), 11: (-3.93345E+02, -9.51778E-01, -2.19071E-03, 3.05843E-01),
    12: (-1.44334E+02, -1.85617E+00, -2.37491E-05, 2.88290E-01), 13: (4.05265E+01, -2.43780E+00, 3.60133E-03, 9.86422E-02),
    14: (-7.07406E+01, -1.06674E+00, 3.73104E-03, -1.44003E-01), 15: (-7.83648E+02, -1.03439E+00, -1.07148E-03, 3.15535E-01),
    16: (-5.63872E+02, -7.18416E-01, 2.63016E-03, -4.99470E-02), 17: (-4.53961E+02, -3.26105E-01, -1.39780E-04, -3.93916E-02),
    18: (3.71375E+01, -2.66753E+00, 1.01483E-03, 2.14233E-01), 19: (-5.03710E+02, 1.04092E+00, -4.12746E-03, 1.82790E-01),
    20: (-3.59763E+01, -4.08458E-01, 1.67264E-03, -9.98919E-02), 21: (-6.09432E+02, 1.50436E+00, -9.09024E-04, -1.35495E-01),
    22: (-1.02367E+02, -7.16253E-01, -2.90670E-04, -5.88556E-01), 23: (-1.93802E+03, 6.48262E-01, 1.73245E-03, 3.47940E-02),
    24: (-5.26919E+00, 3.06435E-01, 3.25397E-03, -6.81506E-01), 25: (-2.84042E+02, -6.25424E-01, -8.22474E-04, -8.80549E-02),
    26: (1.50093E+02, 2.39875E-02, -3.37969E-03, 1.52789E-02), 27: (-2.03387E+01, -5.48718E+00, 8.39075E-03, 1.07884E-01),
    28: (-8.38064E+02, -1.09600E+00, -4.24385E-04, 2.81812E-01), 29: (-5.27934E+01, -4.63689E-01, -5.11647E-03, 3.84965E-01),
    30: (-1.61520E+03, 9.01669E-01, 1.44536E-03, 2.66889E-01),
}
_SIMPOL_SMARTS = {
    1: '[#6]', 2: '[#6][!$([NX3][CX3]=O),$([CX3](=O)[NX3])]',
    301: '[#6,#7,#8;r3;a]', 302: '[#6,#7,#8;r4;a]', 303: '[#6,#7,#8;r5;a]', 304: '[#6,#7,#8;r6;a]', 305: '[#6,#7,#8;r7;a]', 306: '[#6,#7,#8;r8;a]',
    401: '[#6,#7,#8;r3;A]', 402: '[#6,#7,#8;r4;A]', 403: '[#6,#7,#8;r5;A]', 404: '[#6,#7,#8;r6;A]', 405: '[#6,#7,#8;r7;A]', 406: '[#6,#7,#8;r8;A]',
    5: '[CX3]=[CX3]', 6: '[CX3;r]=[CX3][CX3;r]=[OX1]', 7: '[OX2;H1][CX4]', 8: '[CX3;H1](=O)[#6]', 9: '[#6][CX3](=O)[#6]',
    10: '[#6][CX3](=[OX1])[OX2;H1]',
    11: '[CX4;!$([NX3](=O)=O),!$([NX3+](=O)[O-])][CX3](=[OX1])[OX2][CX4;!$([CX3]=[#8])]',
    12: '[C;!$([CX3]=[#8])][OX2;!R][C;!$([CX3]=[#8])]', 13: '[CX4;!$([CX3]=[#8])][OX2;R][CX4;!$([CX3]=[#8])]',
    14: '[OX2;!$([CX3]=[O])]([#6;a;R])[#6]',
    15: '[$([OX2;!$(O-O)][NX3](=[OX1])(=[OX1])),$([OX2;!$(O-O)][NX3+]([OX1-])(=[OX1]))]',
    1601: '[#6][$([NX3](=O)=O)]', 1602: '[#6][$([NX3+](=O)[O-])]', 17: '[OX2H][cX3]:[c]',
    18: '[NX3;H2][CX4]', 19: '[CX4][NX3;H1][CX4]', 20: '[CX4][NX3]([CX4])[CX4]',
    21: '[c][NX3;!$(N~O);!$(N~N);!$(NC~O);!$(NC~N)]',
    22: '[CX3]([#6])(=O)[NX3;H2]', 23: '[CX3]([#6])(=O)[NX3;H1][#6]', 24: '[CX3]([#6])(=O)[NX3]([#6])[#6]',
    2501: '[$([CX3]([#6])(=O)[OX2][OX2][NX3](=O)=O)]', 2502: '[$([CX3]([#6])(=O)[OX2][OX2][NX3+](=O)[O-])]',
    26: '[#6][OX2][OX2][#6]', 27: '[#6;!$(C=O)][OX2][OX2;H1]', 28: '[#6][CX3](=O)[OX2][OX2;H1]',
    2901: '[NX3](=O)(=O)cc[OX2;H1]', 2902: '[NX3](=O)(=O)ccc[OX2;H1]', 2903: '[NX3](=O)(=O)cccc[OX2;H1]',
    2904: '[NX3+](=O)([O-])cc[OX2;H1]', 2905: '[NX3+](=O)([O-])ccc[OX2;H1]', 2906: '[NX3+](=O)([O-])cccc[OX2;H1]',
    3001: '[$([#6][OX2][CX3](=O)[CX4;H2][NX3](=O)=O)]', 3002: '[$([#6][OX2][CX3](=O)[CX4;H2][NX3+](=O)[O-])]',
}
_SIMPOL_PATT: Dict[int, object] = {}


def _mol(smiles: str):
    from rdkit import Chem
    from rdkit import RDLogger
    RDLogger.DisableLog("rdApp.*")
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("SMILES not readable")
    return mol


def simpol_groups(smiles: str) -> Dict[int, float]:
    """SIMPOL.1 group counts of a SMILES (group number -> count), as UManSysProp's groups.simpol combines them."""
    from rdkit import Chem
    if not _SIMPOL_PATT:
        for k, s in _SIMPOL_SMARTS.items():
            _SIMPOL_PATT[k] = Chem.MolFromSmarts(s)
    mol = _mol(smiles)
    m = {k: len(mol.GetSubstructMatches(p, uniquify=True)) for k, p in _SIMPOL_PATT.items()}
    r: Dict[int, float] = {0: 1}
    for k in (1, 2, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 18, 19, 20, 21, 22, 23, 24, 26, 27, 28):
        r[k] = m[k]
    r[3] = math.ceil(m[301] / 3 + m[302] / 4 + m[303] / 5 + m[304] / 6 + m[305] / 7 + m[306] / 8)
    r[4] = math.ceil(m[401] / 3 + m[402] / 4 + m[403] / 5 + m[404] / 6 + m[405] / 7 + m[406] / 8)
    r[16] = m[1601] + m[1602]
    r[17] = m[17] / 2
    if r[22] + r[23] + r[24] == 0:
        r[2] = 0
    r[25] = m[2501] + m[2502]
    r[29] = m[2901] + m[2902] + m[2903] / 2 + m[2904] + m[2905] + m[2906] / 2
    if r[29] > 0:
        r[17] = 0
    r[30] = m[3001] + m[3002]
    if r[30] > 0:
        r[11] = 0
    return {k: v for k, v in r.items() if v}


def _mass(m: Molecule) -> float:
    if m.mass:
        return m.mass
    from rdkit.Chem import Descriptors
    return Descriptors.MolWt(_mol(m.smiles))


def cstar_simpol(m: Molecule, T: float = T_REF) -> float:
    """SIMPOL.1, Pankow and Asher (2008), ACP 8, 2773: log10(p/atm) = sum_k nu_k (B1/T + B2 + B3 T + B4 ln T)."""
    log10p = sum(n * (c[0] / T + c[1] + c[2] * T + c[3] * math.log(T))
                 for k, n in simpol_groups(m.smiles).items() for c in [_SIMPOL_COEF[k]])
    return 1e6 * _mass(m) * 10 ** log10p / (R_ATM * T)


def is_aromatic(smiles: str) -> bool:
    try:
        return any(a.GetIsAromatic() for a in _mol(smiles).GetAtoms())
    except Exception:
        return False


# ── UManSysProp (Nannoolal, EVAPORATION) ───────────────────────────────────────────────────────────────────────

_BACKEND: Dict[str, object] = {}


def umansysprop_backend() -> Tuple[bool, str]:
    """(available, why not).  Loaded once; the checkout named by UMANSYS is tried when the package is not installed."""
    if "ok" in _BACKEND:
        return bool(_BACKEND["ok"]), str(_BACKEND.get("why", ""))
    try:
        from juice._root import ROOT
        for path in (os.environ.get("UMANSYS", ""),                       # a checkout named by the user
                     os.path.join(os.path.dirname(ROOT), "_external", "UManSysProp_public")):   # next to JUICE
            if path and os.path.isfile(os.path.join(path, "umansysprop", "vapour_pressures.py")):
                if path not in sys.path:
                    sys.path.append(path)
                break
        from openbabel import pybel
        from umansysprop import boiling_points, vapour_pressures
        _BACKEND.update(ok=True, pybel=pybel, vp=vapour_pressures, bp=boiling_points)
    except Exception as e:
        hint = ('python -m pip install "setuptools<81" (UManSysProp needs pkg_resources)' if "pkg_resources" in str(e) else
                "python -m pip install openbabel-wheel" if "openbabel" in str(e) else
                "unpack github.com/loftytopping/UManSysProp_public as GitHub/_external/UManSysProp_public, or set UMANSYS")
        _BACKEND.update(ok=False, why=f"needs Open Babel and UManSysProp ({type(e).__name__}: {str(e)[:60]}); {hint}")
    return bool(_BACKEND["ok"]), str(_BACKEND.get("why", ""))


def _need_backend():
    ok, why = umansysprop_backend()
    if not ok:
        raise RuntimeError(why)


def cstar_nannoolal(m: Molecule, T: float = T_REF) -> float:
    """Nannoolal et al. (2008) vapour pressure with the Nannoolal et al. (2004) boiling point (UManSysProp)."""
    _need_backend()
    mol = _BACKEND["pybel"].readstring("smi", m.smiles)
    lp = _BACKEND["vp"].nannoolal(mol, T, _BACKEND["bp"].nannoolal(mol))
    return 1e6 * (m.mass or mol.molwt) * 10 ** lp / (R_ATM * T)


def cstar_evaporation(m: Molecule, T: float = T_REF) -> float:
    """EVAPORATION, Compernolle et al. (2011), ACP 11, 9431 (UManSysProp). Not for aromatic rings."""
    _need_backend()
    mol = _BACKEND["pybel"].readstring("smi", m.smiles)
    return 1e6 * (m.mass or mol.molwt) * 10 ** _BACKEND["vp"].evaporation(mol, T) / (R_ATM * T)


LI_N_OFFSET = 2.977     # log10 offset of the nitrate fallback, calibrated against EVAPORATION (F. Rohrer, 2026-10-12)


def _li_n_fallback(nC, nO, nN):
    """F. Rohrer's fit for an organic nitrate without a usable structure: the CHO corridor on the oxygens that are not
    in the nitrate group, with bCO = -0.0779 (not Li's -0.779 for CHO) and LI_N_OFFSET. The two numbers were fitted
    together against EVAPORATION, so they stay as fitted; this is an empirical fit, not Li et al. (2016)."""
    return _corridor(nC, max(1, nO - 2 * nN), 22.66, 0.4481, 1.656, -0.07790) * 10 ** LI_N_OFFSET


def cstar_li_evap(m: Molecule, T: float = T_REF) -> float:
    """Li et al. (2016) for CHO; with nitrogen: Li (CHON) for an aromatic ring, EVAPORATION for a structure, else Li
    with the nitrate oxygens de-weighted and the EVAPORATION-calibrated offset (F. Rohrer's b2_volatility routing)."""
    if m.nN <= 0:
        return _li_cho(m.nC, m.nO)
    if m.smiles and is_aromatic(m.smiles):
        return cstar_li(m)
    if m.smiles and umansysprop_backend()[0]:
        try:
            c = cstar_evaporation(m, T)
            if c > 0 and c == c:
                return c
        except Exception:
            pass
    return _li_n_fallback(m.nC, m.nO, m.nN)


# ── the methods offered ────────────────────────────────────────────────────────────────────────────────────────

@dataclass
class Method:
    id: str
    label: str
    needs: str                      # "formula" | "smiles"
    reference: str
    fn: Callable[[Molecule], float]
    backend: str = ""               # "" | "umansysprop"

    def available(self) -> Tuple[bool, str]:
        if self.backend == "umansysprop":
            return umansysprop_backend()
        if self.needs == "smiles":
            try:
                import rdkit  # noqa: F401
            except Exception as e:
                return False, f"needs RDKit ({e})"
        return True, ""


METHODS: "OrderedDict[str, Method]" = OrderedDict((m.id, m) for m in (
    Method("simpol", "SIMPOL.1, group contributions from the structure", "smiles",
           "Pankow and Asher (2008), Atmos. Chem. Phys. 8, 2773", cstar_simpol),
    Method("nannoolal", "Nannoolal, from the structure (as GECKO-A)", "smiles",
           "Nannoolal et al. (2008), Fluid Phase Equilib. 269, 117; boiling point Nannoolal et al. (2004)", cstar_nannoolal, "umansysprop"),
    Method("evaporation", "EVAPORATION, from the structure (not for aromatics)", "smiles",
           "Compernolle et al. (2011), Atmos. Chem. Phys. 11, 9431", cstar_evaporation, "umansysprop"),
    Method("li_evap", "Li et al. (2016), with EVAPORATION for organic nitrogen", "formula",
           "Li et al. (2016), Atmos. Chem. Phys. 16, 3327; Compernolle et al. (2011); routing F. Rohrer", cstar_li_evap),
    Method("li2016", "Li et al. (2016), volatility basis set from the formula", "formula",
           "Li et al. (2016), Atmos. Chem. Phys. 16, 3327", cstar_li),
    Method("donahue2011", "Donahue et al. (2011), 2-D volatility basis set from the formula", "formula",
           "Donahue et al. (2011), Atmos. Chem. Phys. 11, 3303", cstar_donahue),
))
DEFAULT_METHOD = "simpol"


def method(method_id: str) -> Method:
    if method_id not in METHODS:
        raise KeyError(f"unknown volatility method '{method_id}' (choose {', '.join(METHODS)})")
    return METHODS[method_id]


def cstar(method_id: str, m: Molecule) -> float:
    c = float(method(method_id).fn(m))
    if not (c > 0) or c != c or math.isinf(c):
        raise ValueError("no finite C*")
    return c


# ── the species of a mechanism ─────────────────────────────────────────────────────────────────────────────────

_TOKEN = re.compile(r"^(?:[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?\s*\*?\s*)?([A-Za-z_][A-Za-z0-9_]*)$")


def mechanism_species(texts: Iterable[str]) -> List[str]:
    """Species (upper case, first appearance first) in the reactions of EASY text: bare lines and k[...] keys."""
    seen, order = set(), []
    for text in texts:
        for raw in text.splitlines():
            line = raw.split(";")[0].strip()
            if "-->" not in line or "-/->" in line:
                continue
            if line.lower().startswith("k["):
                mm = re.match(r"k\[(.*?)\]", line, re.I)
                if not mm:
                    continue
                eq = mm.group(1)
            else:
                eq = re.sub(r"\{[^}]*\}", "", line)
            for side in eq.split("-->"):
                for term in side.split("+"):
                    mm = _TOKEN.match(term.strip())
                    if mm and mm.group(1).upper() not in seen:
                        seen.add(mm.group(1).upper()); order.append(mm.group(1).upper())
    return order


def _is_radical(smiles: str) -> bool:
    try:
        from rdkit.Chem import Descriptors
        return Descriptors.NumRadicalElectrons(_mol(smiles)) > 0
    except Exception:
        return "[O]" in smiles or "[C]" in smiles or "[CH" in smiles


@dataclass
class Result:
    method: str
    species: List[tuple]            # (name, MW, nC, nO, nN, 'C*', class) as emit_soa takes them
    skipped: Dict[str, str]         # name -> why it is not partitioned
    above_limit: int = 0            # species left out because C* is above the limit

    def summary(self) -> str:
        n = {g: sum(1 for s in self.species if volatility_class(float(s[5])) == g) for g in ("ELVOC", "LVOC", "SVOC", "VOC")}
        return (f"{len(self.species)} species partition ({n['ELVOC']} ELVOC, {n['LVOC']} LVOC, {n['SVOC']} SVOC, "
                f"{n['VOC']} more volatile)"
                + (f"; {self.above_limit} above the C* limit" if self.above_limit else "")
                + (f"; {len(self.skipped)} without a C*" if self.skipped else ""))


def partitioning_species(names: Iterable[str], db, method_id: str = DEFAULT_METHOD, cstar_max: float = 0.0) -> Result:
    """C* by `method_id` for the closed-shell organic species among `names` (`db`: a species_db.SpeciesDB).
    Radicals, CO, CO2 and CH4 are left out; a species the method cannot do is listed in `skipped`.
    `cstar_max` > 0 leaves out everything more volatile (fewer reactions; 6 are written per species)."""
    meth = method(method_id)
    ok, why = meth.available()
    if not ok:
        raise RuntimeError(f"volatility method '{method_id}': {why}")
    out, skipped, above = [], {}, 0
    for n in names:
        if n not in db or n in NOT_PARTITIONED:
            continue
        sp = db[n]
        if sp.nC < 1 or sp.is_ro2:
            continue
        smiles = sp.smiles if sp.smiles not in ("-", "") else ""
        if smiles and _is_radical(smiles):
            continue
        if not smiles and (n.endswith("O2") or n.endswith("O3") or (n.endswith("O") and not n.endswith("OOH"))):
            continue                                     # a radical by its MCM name, where no structure says so
        if meth.needs == "smiles" and not smiles:
            skipped[n] = "no structure (SMILES)"
            continue
        mol = Molecule(n, sp.nC, sp.nH, sp.nO, sp.nN, float(sp.mass or 0.0), smiles)
        if not mol.mass:
            mol.mass = 12.011 * sp.nC + 1.008 * sp.nH + 15.999 * sp.nO + 14.007 * sp.nN
        try:
            c = cstar(method_id, mol)
        except Exception as e:
            skipped[n] = f"{type(e).__name__}: {str(e)[:50]}"
            continue
        if cstar_max and c > cstar_max:
            above += 1
            continue
        g = volatility_class(c)
        out.append((n, mol.mass, sp.nC, sp.nO, sp.nN, f"{c:.3E}", "SVOC" if g == "VOC" else g))
    return Result(method_id, out, skipped, above)
