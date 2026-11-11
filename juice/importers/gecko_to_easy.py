"""GECKO-A mechanism (OUT/<mech>/) -> EASY/JUICE scenario with gas/particle/wall partitioning.

Writes into <out>/<name>/ (the chemistry as include files, the way EASY scenarios are built):
  anorg.geckoa.incl            GECKO-A's inorganic core (every reaction between core species only), GECKO-A's
                               own species names (HO, HNO2, HNO4, O3P, SULF) -- the counterpart of the
                               anorg.mcm331.* includes; a scenario loads ONE of the two sets
  geckoa.aliases.incl          ALIAS[canonical] = alt ... : maps the GECKO-A names onto the EASY/MCM names
                               (OH, HONO, HNO4, O, SA, HCHO); extend here for further species
  <name>.geckoa.names.incl     GECKO-A dictionary code -> name used in the organic include, as ALIAS lines
                               (every code is accepted as a synonym; override/extend here)
  <name>.geckoa.organic.incl   the organic GECKO-A chemistry (GECKO-A kinetics 1:1; photolysis scaled to the
                               SAPHIR-measured J values via gecko_photolysis) -- identical in both modes
  <name>.dil.incl              X + DIL --> for every gas species (same form as the IDL <stem>.fac.dil.mech/.rate)
  <name>.SOA.include           our 6-reaction Krechmer partitioning (make_soa_condensation.emit_soa) for the
                               species GECKO-A itself partitions (pvap.nan.dat), C* from GECKO-A's Nannoolal Psat,
                               KCOND/KFWALL/KPWALL/CW_SCALE = project constants (fixed SAPHIR aerosol)
  <name>.species.tsv           species table (IDL-reader / species_db 'clean' layout) with SMILES
  <name>.easy                  the scenario; the ENZ named in the template's FILES[ENZ] is copied next to it

Modes (they differ ONLY in the inorganic include set and the chamber lines; swap by (un)commenting):
  mcmcore  SAPHIR template (Ethanol.easy) with its anorg.mcm331 core, j_clear_sky, measured J, wall
           deposition, chamber sources and injections; the template's <stem>.fac.* includes are replaced
           by the GECKO-A organic include (+ names/dil/SOA); anorg.geckoa.incl is written but commented.
  pure     self-contained: anorg.geckoa.incl + the same organic include, no old chemistry include; the
           chamber set-up (ENZ, injections, HONO/HCHO sources, dilution, wall deposition) is written out.

Species names: inorganic core species keep GECKO-A's names inside the GECKO-A includes (the alias include
makes them the EASY names); organics get MCM names where the structure is unambiguous (C2H5OH, CH3CHO,
PAN, HOCH2CH2O2, ...; table MCM_NAMES), otherwise G + dictionary code (e.g. G1A2000). The GECKO code and
structure are kept in the species TSV, the names include and as reaction comments.

Usage:
  python -m juice.importers.gecko_to_easy "Z:\\TEMP\\g.gkatzelis\\Franz\\GECKO-A\\ethanol" --template "O:\\2023\\23-07-14\\EASY\\scenario\\Ethanol\\Ethanol.easy" --name Ethanol_GECKOA_MCMcore --mode mcmcore --out "E:\\Python Code\\JUICE archive\\EASY" [--juice-dir "E:\\Python Code\\JUICE archive\\JUICE"]
"""
from __future__ import annotations
import argparse
import hashlib
import os
import re
import shutil
from collections import OrderedDict
from dataclasses import dataclass, field
from datetime import datetime

from juice.importers.gecko_mechanism import GeckoMechanism, GReaction, read_mechanism, generation_cut
from juice.importers.gecko_photolysis import easy_rate as photolysis_rate, GPHOT_LINE, GPHOT_COMMENT, HELPER_LINES
from juice.importers.make_soa_condensation import emit_soa

R_ATM = 8.20573e-5          # atm m3 / (mol K)
_AW = {'C': 12.011, 'H': 1.008, 'N': 14.007, 'O': 15.999, 'S': 32.06}

# GECKO-A codes whose reactions are the inorganic core (anorg.mcm331 covers them in mcmcore mode)
INORG_CORE = {'O3P', 'O1D', 'O3', 'O2', 'NO', 'NO2', 'NO3', 'N2O5', 'HO', 'HO2', 'H2O2', 'HNO2', 'HNO3',
              'HNO4', 'H2', 'H2O', 'CO', 'CO2', 'SO2', 'SULF'}
ENV_CODES = {'H2O', 'O2', 'N2'}      # EASY settings, never species/products
SKIP_CODES = {'XCLOST'}              # GECKO-A carbon-loss bookkeeping (unused in the ethanol scheme)

MCM_NAMES = {
    'HO': 'OH', 'HNO2': 'HONO', 'HNO4': 'HO2NO2', 'O3P': 'O', 'SULF': 'SA', 'CH2O': 'HCHO',
    'N01001': 'CH3NO3',
    'O02000': 'C2H5OH', 'D02000': 'CH3CHO', 'A02000': 'CH3CO2H', 'G02000': 'CH3CO3H', 'P02000': 'PAN',
    '3K2000': 'CH3CO3', 'DD2000': 'GLYOX',
    '2O2000': 'HOCH2CH2O2', '1O2000': 'HOCH2CH2O', 'OO2000': 'ETHGLY', 'HO2000': 'HYETHO2H',
    'DO2000': 'HOCH2CHO', 'AO2000': 'HOCH2CO2H', '3K2001': 'HOCH2CO3', 'GO2000': 'HOCH2CO3H', 'PO2000': 'PHAN',
    '2D2000': 'HCOCH2O2', '1D2000': 'HCOCH2O', 'HD2000': 'HCOCH2OOH', '3D2000': 'HCOCO3', 'GD2000': 'HCOCO3H',
    'AD2000': 'HCOCO2H',
}
IDENTITY = {'H2', 'H2O', 'H2O2', 'HNO3', 'HO2', 'N2O5', 'NO', 'NO2', 'NO3', 'O1D', 'O2', 'O3', 'SO2',
            'CH3O', 'CH3OH', 'CH3O2', 'CH3OOH', 'CH4', 'HCOOH', 'CO', 'CO2', 'XCLOST'}
# fixed (non-structural) dictionary names: formula, SMILES
FIXED = {
    'H2': ('H2', '[H][H]'), 'H2O': ('H2O', 'O'), 'H2O2': ('H2O2', 'OO'), 'HNO2': ('HNO2', 'O=NO'),
    'HNO3': ('HNO3', 'O=[N+]([O-])O'), 'HNO4': ('HNO4', 'O=[N+]([O-])OO'), 'HO': ('HO', '[OH]'),
    'HO2': ('HO2', '[O]O'), 'N2O5': ('N2O5', 'O=[N+]([O-])O[N+](=O)[O-]'), 'NO': ('NO', '[N]=O'),
    'NO2': ('NO2', '[O]N=O'), 'NO3': ('NO3', '[O][N+](=O)[O-]'), 'O3P': ('O', '[O]'), 'O1D': ('O', '[O]'),
    'O2': ('O2', 'O=O'), 'O3': ('O3', '[O-][O+]=O'), 'SO2': ('O2S', 'O=S=O'), 'SULF': ('H2O4S', 'OS(=O)(=O)O'),
    'XCLOST': ('X', '-'), 'CH2O': ('CH2O', 'C=O'), 'CH4': ('CH4', 'C'), 'CO': ('CO', '[C-]#[O+]'),
    'CO2': ('CO2', 'O=C=O'),
}

BASE_SETTINGS = [
    "; ---- base settings: physics, not chemistry (the lines the IDL MCM converter wrote into",
    ";      <stem>.fac.saphir.incl; identical to the original Ethanol.elog). Species aliases: geckoa.aliases.incl ----",
    "CONST=FUNCTION[%1]",
    "DENSITY=FUNCTION[(%1)/((%2)*1.379E-19)]",
    "M =DENSITY(P,T)",
    "O2=CONST(0.20946*M)",
    "N2=CONST(0.78084*M)",
]
# ALIAS[canonical] = alt1 alt2 ...  (EASY convention, as in general_settings.incl / <stem>.fac.saphir.incl);
# the parsers substitute alternatives in reactions and k[] keys, the canonical name is what the output carries.
ALIAS_LINES = [
    "ALIAS[HV]     = hv",
    "ALIAS[OH]     = HO",
    "ALIAS[HCHO]   = H2CO CH2O",
    "ALIAS[O]      = O3P",
    "ALIAS[CH3OOH] = CH3O2H MOOH",
    "ALIAS[HONO]   = HNO2",
    "ALIAS[HNO3]   = HONO2",
    "ALIAS[HNO4]   = HO2NO2",
    "ALIAS[SA]     = SULF H2SO4",
]
ALIAS_INCLUDE = "geckoa.aliases.incl"
ANORG_INCLUDE = "anorg.geckoa.incl"
MCM_ANORG_INCLUDES = ["anorg.mcm331.dark.rate", "anorg.mcm331.dark.mech", "anorg.mcm331.phot.rate",
                      "anorg.mcm331.phot.mech", "j_clear_sky.incl"]
# FUNCTION definitions and parameter lines carry NO trailing comment: the EASY parsers' FUNCTION regex needs
# ']' at the line end, and a line comment made them fall through as (uncompilable) parameters -> rate 0.
GTROE_COMMENT = "; GTROE(k0*M, kinf, Fc) = standard Troe form of boxmodel4gecko forate (flag 0)"
GTROE_LINE = "GTROE = FUNCTION[((%1)/(1.+(%1)/(%2)))*((%3)@(1./(1.+(LOG10((%1)/(%2)))@2.)))]"


def code_name(code: str) -> str:
    """G<code> as an EASY identifier. GECKO-A codes are case-sensitive (G2T700k and G2T700K are different
    species) while EASY/JUICE upper-case every name, so a lower-case letter x becomes '_X' (injective:
    upper-case letters never carry an underscore)."""
    return 'G' + ''.join('_' + ch.upper() if ch.islower() else ch for ch in code)


def easy_name(code: str) -> str:
    """Canonical EASY name (what the output files carry): MCM name, identity, or G<code> (code_name)."""
    if code in MCM_NAMES:
        return MCM_NAMES[code]
    if code in IDENTITY:
        return code
    return code_name(code)


def gecko_name(code: str) -> str:
    """Name written inside the GECKO-A includes: inorganic core species keep GECKO-A's own name
    (HO, HNO2, HNO4, O3P, SULF -> mapped by geckoa.aliases.incl), everything else as easy_name."""
    if code in INORG_CORE:
        return code
    return easy_name(code)


def is_inorganic(r: GReaction) -> bool:
    return r.species <= INORG_CORE


# ----------------------------------------------------------------------------- number formatting
def _real(x: float) -> str:
    """FACSIMILE-safe real literal (always carries a decimal point or exponent)."""
    if x == 0:
        return '0.'
    if float(x).is_integer() and abs(x) < 1e7:
        return f"{int(x)}."
    s = f"{x:.6g}"
    if '.' not in s and 'e' not in s and 'E' not in s:
        s += '.'
    return s.upper()


def _coef(c: float) -> str:
    if abs(c - 1.0) < 1e-12:
        return ''
    return f"{c:g}*"


def _arrh(A: float, n: float, Ea: float, t300: bool = False) -> str:
    """k = A * T^n * exp(-Ea/T)  (t300: (T/300)^n, the falloff convention)."""
    parts = [f"{A:.3E}"]
    if n != 0.0:
        parts.append(f"(T/300.)@({_real(n)})" if t300 else f"T@({_real(n)})")
    if Ea != 0.0:
        parts.append(f"EXP({_real(-Ea)}/T)")
    return '*'.join(parts)


# ----------------------------------------------------------------------------- conversion
@dataclass
class EasyRxn:
    reactants: list                 # display names, pseudo reactants (HV, MEPERO, PEROj) last
    products: list                  # [(coef, name)]
    terms: list = field(default_factory=list)   # CONST-summable rate terms
    bare: str = ''                  # bare FUNCTION-call rate (photolysis GPHOT)
    comments: list = field(default_factory=list)

    @property
    def key(self):
        return (tuple(sorted(self.reactants)), tuple(sorted((round(c, 6), n) for c, n in self.products)))

    def equation(self) -> str:
        lhs = ' + '.join(self.reactants)
        rhs = ' + '.join(f"{_coef(c)}{n}" for c, n in self.products)
        return f"{lhs} --> {rhs}" if rhs else f"{lhs} --> "

    def rate(self) -> str:
        if self.bare:
            if self.terms:
                raise ValueError(f"cannot merge a bare rate with CONST terms: {self.equation()}")
            return self.bare
        return f"CONST({' + '.join(self.terms)})"


@dataclass
class Group:
    """One include's worth of reactions (inorganic core or organic chemistry)."""
    reactions: list = field(default_factory=list)   # merged EasyRxn in first-appearance order
    params: list = field(default_factory=list)      # parameter lines (Troe, EXTRA helpers, GPHOT)
    merged: int = 0


@dataclass
class Converted:
    inorganic: Group                # -> anorg.geckoa.incl
    organic: Group                  # -> <name>.geckoa.organic.incl
    counters: list                  # MEPERO/PEROj definition lines (organic include)
    dropped: dict                   # reason -> count
    gas_species: list               # canonical names of all gas species in the scenario (incl. core)
    rescaled: dict = field(default_factory=dict)   # code -> lifetime factor f (< 1)


# ---- lifetime rescaling of ultrafast intermediates --------------------------------------------------------
# GECKO-A writes barrier-less alkoxy decompositions as k = A*T^1.7 (~2e13 s-1 at 298 K). Together with rates of
# 1e-5 s-1 the Jacobian spans >1e18 and the BDF Newton/LU fails at ENZ discontinuities. MCM's convention is
# KDEC = 1e6 s-1. We scale EVERY loss channel of such a species by the same factor f = KMAX/kmax (unimolecular,
# O2/M, NO/HO2/..., counters, HV): branching ratios are exactly preserved at all conditions, only the lifetime of
# the transient (< 1 us) becomes ~1/KMAX. Reference conditions for the pseudo-first-order estimate below.
KMAX = 1.0e6                     # s-1, MCM KDEC
REF = {'T': 298.0, 'M': 2.46e19, 'H2O': 3.7e17, 'RO2': 1e9, 'NO': 2.5e10, 'NO2': 2.5e10, 'HO2': 2.5e8, 'OH': 2.5e6,
       'O3': 2.5e12, 'NO3': 2.5e7, 'CO': 2.5e12, 'H2O2': 2.5e10, 'HONO': 2.5e9, 'HNO3': 2.5e10, 'SO2': 2.5e9,
       'O3P': 1e5, 'O1D': 1.0, 'H2': 1.3e13, 'NO3_': 0.0}


def k_ref(r: GReaction) -> float:
    """Pseudo-first-order loss rate (s-1) of the organic reactant of r at REF conditions (0 for photolysis)."""
    import math
    T, M = REF['T'], REF['M']
    k = r.A * T ** r.n * math.exp(-r.Ea / T)
    kind = r.kind
    if kind == 'hv':
        return 0.0
    if kind == 'tbody':
        k *= M
    elif kind == 'oxygen':
        k *= 0.2 * M
    elif kind in ('mepero', 'pero'):
        k *= REF['RO2']
    elif kind == 'isom':
        c = r.aux
        k *= (c[0] * T ** 4 + c[1] * T ** 3 + c[2] * T ** 2 + c[3] * T + c[4])
    elif kind == 'falloff':
        c = r.aux
        k0 = c[0] * (T / 300) ** c[1] * math.exp(-c[2] / T) * M
        ki = r.A * (T / 300) ** r.n * math.exp(-r.Ea / T)
        x = k0 / ki
        k = k0 / (1 + x) * c[3] ** (1 / (1 + math.log10(x) ** 2))
    elif kind == 'extra':
        lab = int(round(r.aux[0]))
        if lab == 100:
            k *= 0.2 * M * M
        elif lab in (500, 501):
            k *= REF['H2O']
        elif lab == 502:
            kd = 4.7856e-4 * math.exp(1851 / T - 5.10485e-3 * T) * 8.314 * T * 1e6 / (1.01325e5 * 6.02e23)
            k *= kd * REF['H2O'] ** 2
    for _, c in r.reactants:
        if c in INORG_CORE:
            k *= REF.get({'HO': 'OH', 'HNO2': 'HONO', 'HNO4': 'HNO3', 'SULF': 'SO2'}.get(c, c), 1e9)
    return k


def lifetime_scales(mech: GeckoMechanism, kmax: float = KMAX) -> dict:
    """code -> f (< 1) for every organic species whose fastest loss channel exceeds kmax at REF conditions."""
    fastest = {}
    for r in mech.reactions:
        if r.kind == 'transfer':
            continue
        orgs = [c for _, c in r.reactants if c in mech.species and c not in INORG_CORE]
        if len(orgs) != 1:
            continue
        k = k_ref(r)
        if k > fastest.get(orgs[0], 0.0):
            fastest[orgs[0]] = k
    return {c: kmax / k for c, k in fastest.items() if k > kmax}


def _convert_one(r: GReaction, mech: GeckoMechanism, state: dict, prefix: str, scale: dict | None = None) -> tuple:
    """-> (EasyRxn | None, [param lines], drop_reason | '').  state: running counters for parameter names,
    prefix: 'A' (anorg include) or 'G' (organic include) so parameter names never clash between the two;
    scale: lifetime_scales() factors applied to every channel of the rescaled species."""
    if r.kind == 'transfer':
        return None, [], 'gas<->particle/wall pseudo-reaction (replaced by the SOA include)'
    if any(c in SKIP_CODES for _, c in r.reactants):
        return None, [], 'XCLOST bookkeeping species as reactant'
    reac = [gecko_name(c) for _, c in r.reactants]
    # XCLOST (GECKO-A's lost-carbon counter) is dropped from the products, the reaction itself is kept
    prod = [(c, gecko_name(code)) for c, code in r.products if code not in ENV_CODES and code not in SKIP_CODES]
    if any(code in SKIP_CODES for _, code in r.products):
        state['xclost'] = state.get('xclost', 0) + 1
    e = EasyRxn(reac, prod)
    params = []
    base = _arrh(r.A, r.n, r.Ea)
    src = r.comment.replace(';', '/').strip()
    e.comments.append(f"GECKO-A l.{r.line_no}" + (f" {src}" if src else ''))
    k = r.kind
    if k == 'arrh':
        e.terms = [base]
    elif k == 'tbody':
        e.terms = [f"{base}*M"]
    elif k == 'oxygen':
        e.terms = [f"{base}*0.2*M"]
    elif k == 'mepero':
        e.reactants.append('MEPERO'); e.terms = [base]
    elif k == 'pero':
        if not mech.pero.get(r.pero_class):
            return None, [], f'RO2 counter class PERO{r.pero_class} empty in this mechanism (rate 0)'
        e.reactants.append(f'PERO{r.pero_class}'); e.terms = [base]
    elif k == 'isom':
        c = r.aux
        if len(c) < 5:
            raise ValueError(f"ISOM line {r.line_no}: expected 5 polynomial coefficients, got {c}")
        poly = ' + '.join(f"({ci:.4E})*T@{p}" if p > 1 else (f"({ci:.4E})*T" if p == 1 else f"({ci:.4E})")
                          for ci, p in zip(c[:5], (4, 3, 2, 1, 0)) if ci != 0.0) or '0.'
        e.terms = [f"{base}*({poly})"]
        e.comments.append("ISOM T-polynomial")
    elif k == 'falloff':
        c = r.aux
        if len(c) < 5:
            raise ValueError(f"FALLOFF line {r.line_no}: expected 7 coefficients, got {c}")
        if c[4] != 0.0:
            raise NotImplementedError(f"FALLOFF flag {c[4]} (MCM broadening) at line {r.line_no} not supported")
        state['fo'] += 1
        pname = f"KFO_{prefix}{state['fo']}"
        k0 = _arrh(c[0], c[1], c[2], t300=True) + "*M"
        ki = _arrh(r.A, r.n, r.Ea, t300=True)
        params.append(f"; {pname}: {e.equation()}")
        params.append(f"{pname} = GTROE({k0}, {ki}, {_real(c[3])})")
        e.terms = [pname]
    elif k == 'extra':
        c = r.aux
        label = int(round(c[0]))
        if label == 100:
            e.terms = [f"{base}*0.2*M*M"]
        elif label == 500:
            e.terms = [f"{base}*H2O"]
        elif label == 501:
            e.terms = [f"{base}*H2O", f"{c[1]:.3E}*EXP({_real(-c[3])}/T)*H2O*M"]
        elif label == 502:
            # water dimer: k * Kd(T) * [H2O]^2, Kd = 4.7856e-4*exp(1851/T - 5.10485e-3*T) atm-1 -> cm3 (Scribano 2006)
            e.terms = [f"{base}*(4.7856E-4*EXP(1851./T-5.10485E-3*T)*8.314*T*1.E6/(1.01325E5*6.02E23))*H2O*H2O"]
        elif label == 550:
            state['ex'] += 1
            p2, p3 = f"KEX_{prefix}{state['ex']}_K2", f"KEX_{prefix}{state['ex']}_K3"
            params.append(f"; {p2}, {p3}: EXTRA 550 k2/k3 of {e.equation()}")
            params.append(f"{p2} = CONST({_arrh(c[1], c[2], c[3])})")
            params.append(f"{p3} = CONST({_arrh(c[4], c[5], c[6])})")
            e.terms = [base, f"{p3}*M/(1.+{p3}*M/{p2})"]
        else:
            raise NotImplementedError(f"EXTRA label {label} at line {r.line_no} not supported")
        e.comments.append(f"EXTRA {label}")
    elif k == 'hv':
        label, hvf = r.aux
        e.reactants.append('HV')
        rate = photolysis_rate(label, 1.0)
        if rate.startswith('CONST('):
            inner = rate[6:-1]                       # measured J straight from the ENZ
        else:
            # like MCM's J<n> = PHTFKT(...): the GPHOT call becomes a named parameter (evaluated in the
            # settings pass, where COSX/jno2 live); a bare call inside k[...] is not evaluated there.
            pname = f"KHV_{prefix}{label}"
            if pname not in state['hv']:
                state['hv'][pname] = rate
                params.append(f"; {pname}: GECKO-A HV label {label}")
                params.append(f"{pname} = {rate}")
            inner = pname
        if hvf != 1.0:
            inner = f"{inner}*{hvf:g}"
        e.terms = [inner]
        e.comments.append(f"HV label {label}")
    else:
        raise NotImplementedError(f"kinetics kind '{k}' at line {r.line_no}")
    if scale:
        orgs = [c for _, c in r.reactants if c in mech.species and c not in INORG_CORE]
        f = min((scale.get(c, 1.0) for c in orgs), default=1.0)
        if f < 1.0:
            e.terms = [f"({t})*{f:.3E}" for t in e.terms]
            e.comments.append(f"lifetime rescaled x{f:.1E} (all channels of {'/'.join(gecko_name(c) for c in orgs if c in scale)}; k_max {KMAX:.0E} s-1)")
    return e, params, ''


SUM_CHUNK = 600     # species per CONST(...) sum line (keeps the compiled expressions short)


def _sum_setting(pname: str, names: list, comment: str) -> list:
    """`pname = CONST( 0.0 + a + b + ... )`, split into pname_1, pname_2, ... partial sums when long."""
    if len(names) <= SUM_CHUNK:
        return [comment, f"{pname} = CONST( 0.0 + {' + '.join(names)} )"]
    L, subs = [comment], []
    for k in range(0, len(names), SUM_CHUNK):
        sub = f"{pname}_{k // SUM_CHUNK + 1}"
        L.append(f"{sub} = CONST( 0.0 + {' + '.join(names[k:k + SUM_CHUNK])} )")
        subs.append(sub)
    L.append(f"{pname} = CONST( {' + '.join(subs)} )")
    return L


def convert(mech: GeckoMechanism, kmax: float | None = KMAX) -> Converted:
    """Split the mechanism into the inorganic core and the organic chemistry; merge duplicate channels
    (same reactants + products) into one summed rate within each group (the runner dedups reactions);
    rescale the lifetimes of ultrafast intermediates (kmax, None = off)."""
    groups = {'A': Group(), 'G': Group()}
    merged = {'A': OrderedDict(), 'G': OrderedDict()}
    states = {'A': {'fo': 0, 'ex': 0, 'hv': {}}, 'G': {'fo': 0, 'ex': 0, 'hv': {}}}
    dropped = {}
    scale = lifetime_scales(mech, kmax) if kmax else {}
    for r in mech.reactions:
        key = 'A' if is_inorganic(r) else 'G'
        e, p, why = _convert_one(r, mech, states[key], key, scale)
        if e is None:
            dropped[why] = dropped.get(why, 0) + 1
            continue
        groups[key].params += p
        m = merged[key]
        if e.key in m:
            tgt = m[e.key]
            if e.bare or tgt.bare:
                raise ValueError(f"duplicate photolysis channel cannot be merged: {e.equation()}")
            tgt.terms += e.terms
            tgt.comments += e.comments
            groups[key].merged += 1
        else:
            m[e.key] = e
    for key in groups:
        groups[key].reactions = list(merged[key].values())
    n_xclost = states['G'].get('xclost', 0) + states['A'].get('xclost', 0)
    if n_xclost:
        dropped[f'XCLOST product removed (lost-carbon counter; reaction kept) in'] = n_xclost
    org = groups['G'].reactions
    counters = []
    if any('MEPERO' in e.reactants for e in org):
        counters.append("; MEPERO = GECKO-A CH3O2 counter (X + MEPERO: X consumed, CH3O2 enters through the rate)")
        counters.append("MEPERO = CONST( 0.0 + CH3O2 )")
    for j in sorted(mech.pero):
        if mech.pero[j] and any(f'PERO{j}' in e.reactants for e in org):
            counters += _sum_setting(f"PERO{j}", [easy_name(c) for c in mech.pero[j]],
                                     f"; PERO{j} = GECKO-A RO2 counter class {j} ({len(mech.pero[j])} species)")
    gas = [easy_name(c) for c in mech.species if c not in ENV_CODES | SKIP_CODES]
    return Converted(groups['A'], groups['G'], counters, dropped, gas, scale)


# ----------------------------------------------------------------------------- writers
def _sha(path: str) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()[:16]


def _write(path: str, lines: list):
    with open(path, 'w', encoding='latin-1', newline='\n') as fh:
        fh.write('\n'.join(lines) + '\n')


def _provenance(mech: GeckoMechanism) -> list:
    files = ['dictionary.out', 'reactionswithcom.dum', 'pvap.nan.dat', 'scheme.log']
    out = [f"; GECKO-A output: {mech.folder}"]
    for f in files:
        p = os.path.join(mech.folder, f)
        if os.path.isfile(p):
            out.append(f";   {f:<22s} sha256[:16] {_sha(p)}")
    return out


KINETICS_NOTES = [
    "; GECKO-A kinetics 1:1 (boxmodel4gecko ratetool.f90 conventions: k = A*T^n*exp(-Ea/T); FALLOFF = standard Troe;",
    ";   TBODY = *M; OXYGEN = *0.2*M; EXTRA 100/500/501/550 written out; MEPERO/PEROj = RO2 counters as settings).",
    "; Photolysis: SAPHIR-measured J (jno2, jo1d, jh2o2, jhono, jhcho_m, jhcho_r) where measured; every other",
    ";   GECKO-A label as GPHOT(l_ratio,dm,-dn,anchor,hvfact) = GECKO-A's own J(SZA) ratio to jno2/jo1d (gecko_photolysis.py).",
]
FUNCTION_BLOCK = ["", "; ---- functions (the same copy stands in every GECKO-A include; identical, last definition wins) ----"] \
    + HELPER_LINES + [GTROE_COMMENT, GTROE_LINE, GPHOT_COMMENT, GPHOT_LINE]


def _rxn_lines(g: Group, what: str) -> list:
    L = []
    if g.params:
        L += ["", f"; ---- rate parameters ({what}: Troe falloff, EXTRA helpers, GPHOT photolysis) ----"] + g.params
    L += ["", f"; ---- reactions ({what}) ----"]
    for e in g.reactions:
        eq = e.equation()
        L.append(f"; {' | '.join(e.comments)}")
        L.append(eq)
        L.append(f"k[{eq}] = {e.rate()}")
    return L


def write_anorg_include(conv: Converted, mech: GeckoMechanism, path: str) -> int:
    g = conv.inorganic
    L = [f"; {os.path.basename(path)} -- GECKO-A inorganic core, GENERATED by gecko_to_easy.py {datetime.now():%Y-%m-%d %H:%M}",
         "; The counterpart of the MCM set anorg.mcm331.dark/phot.rate/.mech (+ j_clear_sky.incl): a scenario loads ONE",
         ";   of the two. Species carry GECKO-A's own names (HO, HNO2, HNO4, O3P, SULF, O1D ...); geckoa.aliases.incl maps",
         ";   them onto the EASY names (OH, HONO, HNO4, O, SA). Needs the measured J inputs (jno2, jo1d, jh2o2, jhono)."]
    L += KINETICS_NOTES + _provenance(mech)
    L += [f"; reactions: {len(g.reactions)} ({g.merged} multi-channel duplicates merged into summed rates)"]
    L += FUNCTION_BLOCK + _rxn_lines(g, 'inorganic core')
    _write(path, L)
    return len(g.reactions)


def write_organic_include(conv: Converted, mech: GeckoMechanism, path: str) -> int:
    g = conv.organic
    L = [f"; {os.path.basename(path)} -- GECKO-A organic chemistry, GENERATED by gecko_to_easy.py {datetime.now():%Y-%m-%d %H:%M}",
         "; Every GECKO-A reaction with at least one non-core species (the inorganic core is anorg.geckoa.incl or anorg.mcm331.*).",
         "; Names: MCM names where unambiguous, else G<GECKO code> (<name>.geckoa.names.incl lists the mapping as ALIAS lines);",
         ";   inorganic partners by their GECKO-A names (HO, HNO2, ...; geckoa.aliases.incl)."]
    L += KINETICS_NOTES + _provenance(mech)
    L += [f"; reactions: {len(g.reactions)} ({g.merged} multi-channel duplicates merged into summed rates)"]
    for why, n in conv.dropped.items():
        L.append(f";   dropped {n:3d}: {why}")
    if conv.rescaled:
        fmin = min(conv.rescaled.values())
        L += [f"; lifetime rescaling: {len(conv.rescaled)} transient species whose fastest loss channel exceeds {KMAX:.0E} s-1 at",
              f";   298 K (barrier-less decompositions up to {KMAX / fmin:.1E} s-1) have ALL their loss channels multiplied by one",
              f";   factor f = {KMAX:.0E}/k_max (f >= {fmin:.1E}): branching ratios unchanged, lifetime ~1/{KMAX:.0E} s (MCM KDEC convention)."]
    L += FUNCTION_BLOCK + ["", "; ---- RO2 counters ----"] + conv.counters + _rxn_lines(g, 'organic')
    _write(path, L)
    return len(g.reactions)


def write_alias_include(path: str):
    L = [f"; {os.path.basename(path)} -- species name aliases GECKO-A <-> EASY/MCM (gecko_to_easy.py)",
         "; ALIAS[canonical] = alt1 alt2 ...  The EASY/JUICE parsers replace every alternative by the canonical name in",
         ";   reactions and k[] keys (NOT inside CONST(...) expressions: settings must use the canonical name). The canonical",
         ";   name is what the output files carry. Extend this file for further species (one line per canonical name).",
         "; First block = the standard SAPHIR aliases (general_settings.incl / <stem>.fac.saphir.incl), then GECKO-A specific."]
    L += ALIAS_LINES
    _write(path, L)


def write_names_include(mech: GeckoMechanism, path: str) -> int:
    L = [f"; {os.path.basename(path)} -- GECKO-A dictionary code -> name used in the organic include (gecko_to_easy.py)",
         "; ALIAS[name] = G<code>: the GECKO-A code (G-prefixed, as a non-MCM species would be called) is accepted as a synonym",
         ";   of the MCM name chosen by the converter. Override or extend here; the comment carries the GECKO-A notation.",
         "; Inorganic core names are in geckoa.aliases.incl; species without an MCM name are called G<code> already."]
    n = 0
    for code, sp in mech.species.items():
        if code in SKIP_CODES or code in INORG_CORE or code not in MCM_NAMES or code == 'CH2O':
            continue
        L.append(f"ALIAS[{MCM_NAMES[code]:<12s}] = {code_name(code):<10s}; {sp.notation}")
        n += 1
    _write(path, L)
    return n


def write_dil_include(names: list, path: str):
    L = [f"; {os.path.basename(path)} -- GENERATED by gecko_to_easy.py: dilution of every gas species (as <stem>.fac.dil.mech/.rate)"]
    for n in names:
        L.append(f"{n} + DIL --> ")
    L.append("")
    for n in names:
        L.append(f"k[{n} + DIL --> ] = CONST(DILUTE)")
    _write(path, L)


CSTAR_MAX = None     # ug/m3; None = partition every species GECKO-A partitions (user decision for ethanol, 2026-11-07)


def soa_species(mech: GeckoMechanism, cstar_max: float | None = None) -> list:
    """(name, MW, nC, nO, nN, Cstar_str, group) for the species GECKO-A partitions (pvap.nan.dat), optionally
    only those with C* <= cstar_max (evaporation k = KCOND*C*/COA reaches 1e7..1e10 s-1 for C* > 1e9 -> stiffness
    failures of the warm BDF in large mechanisms; a species with C* > 1e6 holds < 1e-4 of its mass in the particle
    and < 6 % on the wall at equilibrium)."""
    out = []
    for code, (psat, _dh) in mech.pvap.items():
        sp = mech.species[code]
        cstar = psat * sp.mw * 1e6 / (R_ATM * mech.tref)
        if cstar_max is not None and cstar > cstar_max:
            continue
        grp = 'ELVOC' if cstar < 3e-4 else ('LVOC' if cstar < 0.3 else 'SVOC')
        out.append((easy_name(code), sp.mw, sp.nC, sp.nO, sp.nN, f"{cstar:.3E}", grp))
    return out


def write_soa_include(mech: GeckoMechanism, path: str, cstar_max: float | None = None, wall: str = 'krechmer') -> list:
    """wall='krechmer' = project scheme (Cw = 16/3*C*^0.6, KFWALL 1e-3); wall='gecko' = GECKO-A's own gas-wall
    parameters (k_in 1.111e-3 s-1, Cw = 1.2e-4 mol/m3 * MW for every species, Matsunaga & Ziemann bag) -- sensitivity."""
    sp = soa_species(mech, cstar_max)
    if wall == 'gecko':
        from juice.importers.make_soa_condensation import GECKO_KFWALL
        emit_soa(sp, path, seed=(), cw_mode='gecko', kfwall=GECKO_KFWALL)
    else:
        emit_soa(sp, path, seed=())
    with open(path, encoding='utf-8') as fh:
        txt = fh.read()
    head = (f"; GENERATED by gecko_to_easy.py {datetime.now():%Y-%m-%d %H:%M} via make_soa_condensation.emit_soa\n"
            f"; Species = the {len(sp)} species GECKO-A itself partitions (pvap.nan.dat, {len(mech.pvap)} total"
            + (f", C* <= {cstar_max:g} ug/m3 kept" if cstar_max is not None else '') + f"); C* = Nannoolal Psat @ {mech.tref:g} K\n"
            f";   x MW / (R T) (GECKO-A's own vapour-pressure SAR). Rate constants = project constants (fixed SAPHIR aerosol\n"
            f";   d~160 nm, N~1800 cm-3 -> KCOND; Krechmer Cw). Group labels are the standard VBS bins (all ethanol species\n"
            f";   are above the SVOC range -> no SOA mass expected; the scheme is wired for completeness).\n")
    with open(path, 'w', encoding='latin-1', newline='\n') as fh:
        fh.write(head + txt)
    return sp


def _smiles(sp) -> str:
    if sp.code in FIXED:
        return FIXED[sp.code][1]
    try:
        from juice.importers.gecko_notation_to_smiles import gecko_to_smiles
        return gecko_to_smiles(sp.notation)
    except Exception:
        return '-'


def _mass(formula: str) -> float:
    m = 0.0
    for el, n in re.findall(r'([A-Z][a-z]?)(\d*)', formula):
        m += _AW.get(el, 0.0) * (int(n) if n else 1)
    return m


def write_species_tsv(mech: GeckoMechanism, soa: list, path: str):
    rows = ["Name\tSmiles\tInchi\tInchiKey\tFormula\tMass\tExcited\tPeroxyRadical\tSynonyms"]
    for code, sp in mech.species.items():
        if code in SKIP_CODES:
            continue
        formula = FIXED[code][0] if code in FIXED else sp.formula
        mass = sp.mw if sp.mw > 0 else _mass(formula)
        rows.append('\t'.join([easy_name(code), _smiles(sp), '-', '-', formula or 'X', f"{mass:.4f}",
                               'true' if code == 'O1D' else 'false', 'true' if sp.is_peroxy else 'false',
                               f"GECKO:{code}|{sp.notation}"]))
    for name, *_ in soa:
        for pre in ('SOA_', 'WALL_', 'WALL_PL_'):
            rows.append('\t'.join([pre + name, '-', '-', '-', 'C', '12.011', 'false', 'false', f"{pre}phase of {name}"]))
    with open(path, 'w', encoding='ascii', newline='\r\n') as fh:
        fh.write('\n'.join(rows) + '\n')
    return len(rows) - 1


# ----------------------------------------------------------------------------- scenarios
def _enz_name(template_lines: list) -> str:
    for ln in template_lines:
        m = re.match(r'^\s*FILES\s*\[\s*ENZ\s*\]\s*=\s*(\S+)', ln, re.I)
        if m:
            return m.group(1)
    raise ValueError("template has no FILES[ENZ] line")


def _gecko_block(name: str, mode: str) -> list:
    """The chemistry include lines. mode decides which inorganic core is active; the other stays as a
    commented line so a scenario is switched by (un)commenting."""
    L = ["",
         "*********************************************************************",
         "; ---- chemistry (gecko_to_easy.py). The inorganic core is ONE of two include sets: MCM v3.3.1",
         ";      (anorg.mcm331.dark/phot.rate/.mech + j_clear_sky.incl) or GECKO-A (anorg.geckoa.incl) -- swap by",
         ";      (un)commenting. geckoa.aliases.incl maps the GECKO-A names (HO, HNO2, HNO4, O3P, SULF) onto the EASY names.",
         f"#include {ALIAS_INCLUDE}"]
    if mode == 'pure':
        L += [f"#include {ANORG_INCLUDE}            ; GECKO-A inorganic core (active)"]
        L += [f";#include {inc}" for inc in MCM_ANORG_INCLUDES]
        L += [";#include saphir.jmeasured_02.incl   ; (MCM set only: PHTFKT + measured J definitions)"]
    else:
        L += [f";#include {ANORG_INCLUDE}           ; GECKO-A inorganic core instead of the anorg.mcm331.* set above"]
    L += [f"#include {name}.geckoa.names.incl    ; GECKO-A code -> MCM name synonyms (ALIAS)",
          f"#include {name}.geckoa.organic.incl  ; GECKO-A organic chemistry (same file in both modes)",
          f"#include {name}.dil.incl",
          f"#include {name}.SOA.include   ; gas<->particle<->wall partitioning (Nannoolal C*)",
          ""]
    return L


MCM_CHEM_INCLUDE = re.compile(r'^\s*#include\s+\S+\.fac\.', re.I)      # IDL MCM-converter files (.fac.mech/.rate/.ro2/.saphir/.dil)
PURE_DROP_INCLUDE = re.compile(r'^\s*#include\s+(anorg\.mcm331\.|j_clear_sky\.incl|saphir\.jmeasured)', re.I)
CHAMBER_WORDS = {'HV', 'DIL', 'M', 'NA', 'O2', 'N2', 'H2O', 'NOTHING', 'HV_IN'}
INORG_EASY = {'OH', 'HO', 'HONO', 'HNO2', 'HO2NO2', 'HNO4', 'O', 'O3P', 'O1D', 'O3', 'NO', 'NO2', 'NO3', 'N2O5', 'HO2',
              'H2O2', 'HNO3', 'CO', 'CO2', 'SO2', 'SA', 'SULF', 'H2', 'HCHO', 'CH2O', 'H2SO4'}
_TERM_RE = re.compile(r'^(?:[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?\s*\*?\s*)?([A-Za-z_]\w*)$')


def _eq_sides(eq: str):
    eq = re.sub(r'\{[^}]*\}', '', eq)
    sides = []
    for side in eq.split('-->'):
        names = set()
        for term in side.split('+'):
            m = _TERM_RE.match(term.strip())
            if m:
                names.add(m.group(1).upper())
        sides.append(frozenset(names))
    return sides


def _eq_species(eq: str) -> set:
    return set().union(*_eq_sides(eq))


def _eq_key(eq: str):
    sides = _eq_sides(eq)
    return (sides[0], sides[1] if len(sides) > 1 else frozenset())


def write_scenario_easy(template: str, name: str, out_dir: str, header: list, mode: str, known: set) -> tuple:
    """The scenario .easy from the SAPHIR template: the IDL MCM-converter includes (*.fac.*) are replaced by the
    GECKO-A include block; in `pure` mode the MCM inorganic set (anorg.mcm331.*, j_clear_sky, jmeasured) goes too.
    Template reaction/k[] lines that name MCM-only species are commented out (the chamber lines stay), and in
    `pure` mode so are k[] overrides of reactions that live in the dropped MCM includes (measured-J overrides:
    the GECKO-A core uses the measured J directly). -> (path, [notes])."""
    with open(template, encoding='latin-1') as fh:
        tl = fh.read().splitlines()
    assigned, bare = set(), set()
    for ln in tl:
        s = ln.split(';')[0].strip()
        m = re.match(r'^(?:INITIAL\s*\[\s*)?([A-Za-z_]\w*)\s*\]?\s*=', s)
        if m:
            assigned.add(m.group(1).upper())
        if s and '-->' in s and '-/->' not in s and not s.lower().startswith('k['):
            bare.add(_eq_key(s))
    allowed = set(known) | INORG_EASY | CHAMBER_WORDS | assigned
    label = {'mcmcore': '(GECKO-A organics + MCM core)', 'pure': '(complete GECKO-A mechanism)'}[mode]
    out, notes, base_done = [], [], False
    for i, ln in enumerate(tl):
        s = ln.split(';')[0].strip()
        if i == 0 and ln.upper().startswith(';TITLE'):
            out.append(re.sub(r'Modell:[^!]*', f'Modell:{name} {label}', ln))
            out += header
            continue
        if MCM_CHEM_INCLUDE.match(ln):
            if not base_done:
                out += BASE_SETTINGS
                base_done = True
            out.append(f"; [gecko_to_easy] replaced by the GECKO-A includes below: {ln.strip()}")
            continue
        if mode == 'pure' and PURE_DROP_INCLUDE.match(ln):
            out.append(f"; [gecko_to_easy pure] MCM inorganic set not loaded: {ln.strip()}")
            continue
        if s and '-->' in s and '-/->' not in s:
            is_k = s.lower().startswith('k[')
            if is_k:
                m = re.match(r'k\[(.*?)\]', s, re.I)
                eq = m.group(1) if m else s
            else:
                eq = s
            unknown = _eq_species(eq) - allowed
            if unknown:
                out.append(f"; [gecko_to_easy] MCM-only species {' '.join(sorted(unknown))} -- line disabled: {ln.strip()}")
                notes.append(f"disabled (MCM-only species {', '.join(sorted(unknown))}): {s[:70]}")
                continue
            if mode == 'pure' and is_k and _eq_key(eq) not in bare:
                out.append(f"; [gecko_to_easy pure] rate override of an MCM-include reaction (GECKO-A core has its own): {ln.strip()}")
                notes.append(f"disabled (pure: MCM-include rate override): {s[:70]}")
                continue
        out.append(ln)
    if not base_done:
        out[1:1] = BASE_SETTINGS
    out += _gecko_block(name, mode)
    path = os.path.join(out_dir, f"{name}.easy")
    _write(path, out)
    return path, notes


# ----------------------------------------------------------------------------- driver
def _generator_maxgen(gecko_dir: str) -> int | None:
    """Generation depth the GECKO-A generator was run with (`maxgen = N` in scheme.log); None if not found."""
    p = os.path.join(gecko_dir, 'scheme.log')
    if not os.path.isfile(p):
        return None
    with open(p, encoding='latin-1', errors='replace') as fh:
        for line in fh:
            m = re.match(r'\s*maxgen\s*=\s*(\d+)', line)
            if m:
                return int(m.group(1))
    return None


def template_voc(template_lines: list) -> str | None:
    for ln in template_lines:
        m = re.match(r'^\s*qVOC\s*-->\s*([A-Za-z_]\w*)', ln, re.I)
        if m:
            return m.group(1).upper()
    return None


def build(gecko_dir: str, template: str, name: str, mode: str, out_root: str,
          enz_path: str | None = None, juice_root: str | None = None, log=print,
          maxgen: int | None = None, voc: str | None = None, cstar_max: float | None = CSTAR_MAX,
          wall: str = 'krechmer') -> dict:
    if mode not in ('mcmcore', 'pure'):
        raise ValueError(mode)
    with open(template, encoding='latin-1') as fh:
        tl = fh.read().splitlines()
    mech = read_mechanism(gecko_dir)
    n_full = (len(mech.species), len(mech.reactions))
    # the injected VOC of the template is GECKO-A's primary species: give it the template's name
    voc = (voc or template_voc(tl) or '').upper() or None
    if voc and len(mech.primaries) == 1 and mech.primaries[0] not in MCM_NAMES:
        MCM_NAMES[mech.primaries[0]] = voc
    gens = None
    if maxgen is not None:
        mech, gens = generation_cut(mech, maxgen)
        log(f"[gecko_to_easy] generation cut maxgen={maxgen}: {n_full[0]} -> {len(mech.species)} species, "
            f"{n_full[1]} -> {len(mech.reactions)} reactions (closed-shell organics of generation > {maxgen} are terminal)")
    conv = convert(mech)
    out_dir = os.path.join(out_root, name)
    os.makedirs(out_dir, exist_ok=True)
    for stale in (f"{name}.chem.incl",):            # file names of earlier converter versions
        if os.path.isfile(os.path.join(out_dir, stale)):
            os.remove(os.path.join(out_dir, stale))
    enz_name = _enz_name(tl)
    enz_src = enz_path or os.path.join(os.path.dirname(os.path.abspath(template)), enz_name)
    if not os.path.isfile(enz_src):
        raise FileNotFoundError(enz_src)
    shutil.copy2(enz_src, os.path.join(out_dir, enz_name))

    gen_depth = _generator_maxgen(gecko_dir)          # generator depth from scheme.log ("full" alone says nothing about it)
    gen_txt = f"GECKO-A generator maxgen={gen_depth}" if gen_depth is not None else "GECKO-A generator depth unknown"
    header = [f"; GENERATED by gecko_to_easy.py {datetime.now():%Y-%m-%d %H:%M}  mode={mode}"
              + (f"  maxgen={maxgen} (static generation cut of the {gen_txt} output: {n_full[0]} -> {len(mech.species)} species)"
                 if maxgen is not None else f"  full = complete {gen_txt} output ({n_full[0]} species, no static cut)"),
              f"; template (chamber set-up): {template}",
              f"; ENZ: {enz_name} (copied from {os.path.dirname(enz_src)})"
              + (f"; injected VOC {voc} = GECKO-A primary {mech.primaries}" if voc else '')] + _provenance(mech)
    n_an = write_anorg_include(conv, mech, os.path.join(out_dir, ANORG_INCLUDE))
    write_alias_include(os.path.join(out_dir, ALIAS_INCLUDE))
    n_names = write_names_include(mech, os.path.join(out_dir, f"{name}.geckoa.names.incl"))
    n_org = write_organic_include(conv, mech, os.path.join(out_dir, f"{name}.geckoa.organic.incl"))
    write_dil_include(conv.gas_species, os.path.join(out_dir, f"{name}.dil.incl"))
    soa = write_soa_include(mech, os.path.join(out_dir, f"{name}.SOA.include"), cstar_max, wall)
    n_tsv = write_species_tsv(mech, soa, os.path.join(out_dir, f"{name}.species.tsv"))
    known = set(conv.gas_species) | {gecko_name(c) for c in mech.species} | {easy_name(c) for c in mech.species}
    easy, notes = write_scenario_easy(template, name, out_dir, header, mode, known)
    for n in notes:
        log(f"    template: {n}")
    nparam = lambda g: sum(1 for p in g.params if not p.startswith(';'))
    if conv.rescaled:
        log(f"[gecko_to_easy] lifetime rescaling: {len(conv.rescaled)} transient species (fastest channel > {KMAX:.0E} s-1, "
            f"min factor {min(conv.rescaled.values()):.1E}), all their channels scaled, branching preserved")
    log(f"[gecko_to_easy] {name} ({mode}): {ANORG_INCLUDE} {n_an} inorganic reactions ({conv.inorganic.merged} merged, "
        f"{nparam(conv.inorganic)} parameters){' [written, not active]' if mode == 'mcmcore' else ''}; "
        f"organic include {n_org} reactions ({conv.organic.merged} merged, {nparam(conv.organic)} parameters), "
        f"counters {[c.split('=')[0].strip() for c in conv.counters if not c.startswith(';')]}; "
        f"{n_names} name aliases, {len(soa)} partitioning species, {n_tsv} TSV rows")
    for why, n in conv.dropped.items():
        log(f"    dropped {n:3d}: {why}")
    log(f"    -> {easy}")
    juice_dir = None
    if juice_root:
        juice_dir = os.path.join(juice_root, name)
        os.makedirs(juice_dir, exist_ok=True)
        for f in os.listdir(out_dir):
            shutil.copy2(os.path.join(out_dir, f), os.path.join(juice_dir, f))
        log(f"    copied inputs -> {juice_dir}")
    return {'easy': easy, 'out_dir': out_dir, 'juice_dir': juice_dir, 'converted': conv, 'mechanism': mech,
            'generations': gens, 'notes': notes}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('gecko_dir', help='GECKO-A output folder (dictionary.out, reactionswithcom.dum, pvap.nan.dat, pero*.dat)')
    ap.add_argument('--template', required=True, help='SAPHIR scenario .easy providing the chamber set-up')
    ap.add_argument('--name', required=True, help='scenario name (= folder and file stem)')
    ap.add_argument('--mode', choices=('mcmcore', 'pure'), default='mcmcore')
    ap.add_argument('--out', required=True, help='EASY root; the scenario goes to <out>/<name>/')
    ap.add_argument('--enz', help='ENZ file (default: FILES[ENZ] of the template, next to it)')
    ap.add_argument('--juice-dir', help='JUICE root; all inputs are copied to <juice-dir>/<name>/ as well')
    ap.add_argument('--maxgen', type=int, help='static generation cut: closed-shell organics beyond this oxidation '
                                               'generation do not react further (size reduction for JUICE)')
    ap.add_argument('--voc', help="name of the injected VOC (default: the template's qVOC --> X line); GECKO-A's primary gets it")
    ap.add_argument('--cstar-max', type=float, default=CSTAR_MAX,
                    help='partition only species with C* <= this (ug/m3); default: all GECKO-A partitions')
    ap.add_argument('--wall', choices=('krechmer', 'gecko'), default='krechmer',
                    help="gas-wall parameters: project scheme (default) or GECKO-A's own (k_in 1.111e-3, Cw = 1.2e-4 mol/m3 * MW)")
    a = ap.parse_args(argv)
    build(a.gecko_dir, a.template, a.name, a.mode, a.out, a.enz, a.juice_dir, maxgen=a.maxgen, voc=a.voc, cstar_max=a.cstar_max,
          wall=a.wall)


if __name__ == '__main__':
    main()
