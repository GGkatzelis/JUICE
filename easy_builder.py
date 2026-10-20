"""
easy_builder.py — Step 3 of the JUICE workflow: build the scene.easy file
=========================================================================

The dialog collects everything the IDL workflow used to put into a SAPHIR
EASY master file by hand (see the scene_03..09 survey, 2026-10-20) and writes
the file through easy_scene.write_scene, so that FACSIMILE and JUICE receive
identical input:

    mechanism include set   → INCLUDES block
    chamber options         → WALLDEP / SOURCE (HONO, HCHO) / BGREAC (Y) / DIL
    measured inputs         → MEASURED block (J values, T, P, SZ, HV, FLOW, H2O)
    FACSIMILE parameters    → FILES[ENZ], FACS[HMAX], FACS[OUTSTEP]
    initial concentrations  → INITIAL[...]
    chamber sources         → tracer dilution, dark HONO, SNO2 light source
    injections              → switch = input(ENZ column) + channel reactions
    extra statements        → verbatim

The form ↔ SceneSpec mapping lives in pure functions (`spec_from_form`,
`form_from_spec`) so it can be tested without Qt.  Loading an existing
scene.easy into the dialog and writing it back is statement-equivalent for the
scene_01..09 test set (check: python easy_builder.py --selftest).
"""

from __future__ import annotations

import os
import re
import sys
import shutil
from dataclasses import dataclass, field, asdict
from typing import Dict, List, Optional, Tuple

from easy_scene import (SceneSpec, Setting, ReactionDef, Injection, Stmt,
                        parse_document, extract_scene, write_scene,
                        statements_equivalent, normalise)

# ─────────────────────────────────────────────────────────────────────────────
# Defaults (from the hand-made scene files; newest convention = scene_06..09)
# ─────────────────────────────────────────────────────────────────────────────

ANORG_INCLUDES = [
    ("anorg.mcm331.dark.rate", "dark anorganic reactions"),
    ("anorg.mcm331.dark.mech", "rate constants of dark reactions"),
    ("anorg.mcm331.phot.rate", "anorganic photolysis reactions"),
    ("anorg.mcm331.phot.mech", "anorganic photolysis frequencies"),
]
MECH_SUFFIXES = [
    (".saphir.incl", "MCM3-spezific constants and functions"),
    (".ro2.incl", "RO2 species"),
    (".mech", "reaction mechanism for given VOC"),
    (".rate", "rate constants for given VOC reactions"),
]
WVC2RH_DEF = "FUNCTION[%1*100./(exp(21.36469-5339.66/T)/P*M)]"
HONO_RATE  = "CONST( 4.7e13 * (1. + (RH/11.6)@2) * exp(-3950./T) * jno2 )"
HCHO_RATE  = "CONST(2.8e6)"
Y_RATE     = "CONST(2.4e-13)"
DILUTE_DEF = "CONST((FLOW+0.010*60)/(276.*3600.))"

# measured photolysis inputs: (variable, default ENZ column, reaction whose k is overridden)
J_INPUTS = [
    ("jno2",    "jno2",    "NO2 + HV --> O3 + NO"),
    ("jo1d",    "jo1d",    "O3 + HV --> O1D "),
    ("jhcho_m", "jhcho_m", "HCHO + HV --> H2 + CO "),
    ("jhcho_r", "jhcho_r", "HCHO + HV --> 2*HO2 + CO "),
    ("jh2o2",   "jh2o2",   "H2O2 + HV --> 2*HO"),
    ("jhono",   "jhono",   "HONO + HV --> OH + NO"),
]
# physical settings: name → (input alias, default ENZ column, comment)
SETTINGS = [
    ("T",    "T_in",   "t",      "temperature in K"),
    ("P",    "",       "",       "pressure in hPa"),
    ("SZA",  "",       "sz",     "solar zenith angle"),
    ("SZ",   "",       "sz",     "solar zenith angle"),
    ("HV",   "HV_in",  "hv",     "switch for light"),
    ("FLOW", "",       "dilute", "flow into the chamber in m^3/h"),
    ("H2O",  "H2O_in", "H2O",    "mixing ratio of water"),
]
SETTING_CONST_FROM_INPUT = {"T": "CONST(T_in)", "HV": "CONST(HV_in)", "H2O": "CONST(1e-2*M*H2O_in)"}


# ─────────────────────────────────────────────────────────────────────────────
# Form model (what the dialog edits)
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class MeasuredRow:
    name: str                 # T, P, SZA, SZ, HV, FLOW, H2O  or jno2, ...
    mode: str = "enz"         # enz (alias + CONST(alias)) | input (NAME = INPUT(col)) | const | off
    column: str = ""          # ENZ column (mode enz)
    value: str = ""           # expression inside CONST() (mode const)
    comment: str = ""


@dataclass
class InjectionRow:
    switch: str               # qNO, Q_ethanol, ...
    column: str               # ENZ column name
    products: str             # "NO" | "NO+XNOy" | "0.9*VOCin1+0.1*VOCin2"
    rate: str                 # ppb/min number  OR  full expression (contains letters)
    comment: str = ""


@dataclass
class EasyForm:
    title: str = ""
    header_comments: List[str] = field(default_factory=list)
    # includes
    mech_stem: str = ""                 # e.g. ethanol.fac
    use_anorg: bool = True
    use_j_clear_sky: bool = True
    includes_pre: List[str] = field(default_factory=list)
    includes_other: List[str] = field(default_factory=list)
    # chamber blocks
    walldep: bool = True
    walldep_value: str = "1."
    hono_source: bool = True
    hono_lhs: str = "SHONO+HV"
    hono_rhs: str = "HONO+XNOy"
    hono_const: str = "2.6*0.7"
    hono_rate: str = HONO_RATE
    hcho_source: bool = True
    hcho_lhs: str = "qSHCHO+HV"
    hcho_rhs: str = "HCHO"
    hcho_const: str = "1"
    hcho_const_name: str = "SHCHO"      # constant may be named differently from the reactant (scene_09)
    hcho_rate: str = HCHO_RATE
    bgreac: bool = True
    y_name: str = "qY"
    y_product: str = "CH3O2"
    y_ppb: str = "250"
    y_const_expr: str = ""              # full expression overrides y_ppb when set
    y_rate: str = Y_RATE
    dilution: bool = True
    dilute_expr: str = DILUTE_DEF
    dil_stem: str = ""                  # stem of the .dil.mech/.dil.rate pair (default: mechanism stem)
    dil_flow_const: str = ""            # FLOW = CONST(6.) inside the DIL block (optional)
    dil_value: str = "CONST(1.)"
    # measured
    jmeasured_include: str = "saphir.jmeasured_02.incl"
    j_rows: List[MeasuredRow] = field(default_factory=list)
    j_overrides: Dict[str, str] = field(default_factory=dict)   # jname → reaction text
    setting_rows: List[MeasuredRow] = field(default_factory=list)
    # files / facs
    enz_file: str = "scene.enz"
    hmax: str = "30"
    outstep: str = "100000"
    extra_facs: List[Tuple[str, str]] = field(default_factory=list)
    # initials: (species, expression, comment); '*' allowed
    initials: List[Tuple[str, str, str]] = field(default_factory=list)
    # chamber sources
    tracer_dilution: List[str] = field(default_factory=lambda: ["XNOy", "NA"])
    dark_hono: bool = True
    dark_hono_rhs: str = "HONO+XNOy"
    dark_hono_rate: str = "1e-012/10"
    light_source: bool = False
    light_source_rhs: str = "NO2+XNOy"
    light_source_rate: str = "jno2/5e-3/3600*HV_in/1.6*0.55"
    light_source_initial: str = "M*1e-9*3.5*0"
    # injections
    injections: List[InjectionRow] = field(default_factory=list)
    chain_reactions: List[Tuple[str, str, str]] = field(default_factory=list)   # (lhs, rhs, rate expr)
    # anything else
    extra_statements: List[str] = field(default_factory=list)


def default_form() -> EasyForm:
    f = EasyForm()
    f.j_rows = [MeasuredRow(j, "enz", col) for j, col, _ in J_INPUTS]
    f.j_overrides = {j: rxn for j, _, rxn in J_INPUTS}
    f.setting_rows = [
        MeasuredRow("T", "enz", "t", comment="temperature in K"),
        MeasuredRow("P", "const", "", "1006.5", "pressure in hPa"),
        MeasuredRow("SZA", "enz", "sz", comment="solar zenith angle"),
        MeasuredRow("SZ", "enz", "sz", comment="solar zenith angle"),
        MeasuredRow("HV", "enz", "hv", comment="switch for light"),
        MeasuredRow("FLOW", "enz", "dilute", comment="flow into the chamber in m^3/h"),
        MeasuredRow("H2O", "enz", "H2O", comment="mixing ratio of water"),
    ]
    f.initials = [("*", "CONST(1.)", "")]
    return f


# ─────────────────────────────────────────────────────────────────────────────
# helpers
# ─────────────────────────────────────────────────────────────────────────────

def read_enz_columns(path: str) -> List[str]:
    cols = []
    try:
        with open(path, encoding="latin-1") as fh:
            for line in fh:
                if line.startswith("COLUMN"):
                    cols.append(line.split("=", 1)[1].split("/")[0].strip())
                elif "&&&" in line:
                    break
    except OSError:
        pass
    return [c for c in cols if c.upper() != "TIME"]


def scan_mechanisms(include_dir: str) -> Dict[str, str]:
    """stem (e.g. 'ethanol.fac') → directory containing <stem>.mech."""
    found = {}
    if include_dir and os.path.isdir(include_dir):
        for root, _dirs, files in os.walk(include_dir):
            for f in files:
                if f.lower().endswith(".fac.mech"):
                    found[f[:-len(".mech")]] = root
    return dict(sorted(found.items()))


_NUMERIC = r"[0-9.eE+\-*/() ]+"


def _rate_expr(rate: str, ppb: bool = True) -> str:
    """Injection channel rate: a purely numeric entry means ppb/min and becomes
    CONST(M*1e-9/60*<n>); anything else is the expression inside CONST( ).
    With ppb=False (chain reactions) the entry is always the expression."""
    r = rate.strip()
    if re.match(r"CONST\s*\(", r, re.I):
        return r
    if ppb and re.fullmatch(_NUMERIC, r):
        return f"CONST(M*1e-9/60*{r})"
    return f"CONST({r})"


def _ppb_from_rate(rate_raw: Optional[str], ppb: bool = True) -> str:
    """Inverse of _rate_expr: CONST(M*1e-9/60*<numeric>) → <numeric>; otherwise the
    full expression inside CONST( ) (so that variables such as qNO or QL survive)."""
    if not rate_raw:
        return ""
    if ppb:
        m = re.match(r"CONST\s*\(\s*M\s*\*\s*1e-9\s*/\s*60\s*\*\s*(" + _NUMERIC + r")\s*\)\s*$", rate_raw, re.I)
        if m:
            return m.group(1).strip()
    m = re.match(r"CONST\s*\((.*)\)\s*$", rate_raw.strip(), re.I | re.S)
    return m.group(1).strip() if m else rate_raw


# ─────────────────────────────────────────────────────────────────────────────
# form → SceneSpec
# ─────────────────────────────────────────────────────────────────────────────

def spec_from_form(f: EasyForm) -> SceneSpec:
    s = SceneSpec()
    s.title = f.title
    s.header_comments = list(f.header_comments)
    s.includes_pre = [Stmt(f"#include {n}") for n in f.includes_pre]
    if f.use_anorg:
        s.anorg_includes = [Stmt(f"#include {n}", c) for n, c in ANORG_INCLUDES]
    if f.mech_stem:
        s.mech_stem = f.mech_stem
        s.mech_includes = [Stmt(f"#include {f.mech_stem}{suf}", c) for suf, c in MECH_SUFFIXES]
    if f.use_j_clear_sky:
        s.j_clear_sky = Stmt("#include j_clear_sky.incl")
    s.includes_other = [Stmt(f"#include {n}") for n in f.includes_other]

    if f.walldep:
        s.walldep_include = Stmt("#include saphir.walldep.incl")
        s.walldep = Setting("WALLDEP", f"CONST({f.walldep_value})", "Switch fuer die Wandverluste")

    if f.hono_source or f.hcho_source:
        s.wvc2rh = Setting("WVC2RH", WVC2RH_DEF)
        s.rh = Setting("RH", "WVC2RH(H2O)")
    if f.hono_source:
        s.hono_source = ReactionDef(f.hono_lhs, f.hono_rhs, rate_raw=f.hono_rate,
                                    text_raw=f"{f.hono_lhs} --> {f.hono_rhs}")
        s.hono_const = Setting("SHONO", f"CONST({f.hono_const})")
    if f.hcho_source:
        s.hcho_source = ReactionDef(f.hcho_lhs, f.hcho_rhs, rate_raw=f.hcho_rate,
                                    text_raw=f"{f.hcho_lhs} --> {f.hcho_rhs}")
        s.hcho_const = Setting(f.hcho_const_name or f.hcho_lhs.split("+")[0].strip(), f"CONST({f.hcho_const})")

    if f.bgreac:
        s.y_reaction = ReactionDef(f"{f.y_name} + OH", f.y_product, rate_raw=f.y_rate,
                                   text_raw=f"{f.y_name} + OH --> {f.y_product}")
        expr = f.y_const_expr or f"{f.y_ppb}e-9*M"
        s.y_const = Setting(f.y_name, f"CONST({expr})", "Y in CO Aequivalenten")

    _dstem = f.dil_stem or f.mech_stem
    if f.dilution and _dstem:
        s.dil_includes = [Stmt(f"#include {_dstem}.dil.mech", "dilution mechanism"),
                          Stmt(f"#include {_dstem}.dil.rate", "dilution rate")]
    if f.dilution:
        s.dilute = Setting("DILUTE", f.dilute_expr)
        if f.dil_flow_const:
            s.dil_flow = Setting("FLOW", f"CONST({f.dil_flow_const})", "Experimetal flow in m^3/h")
        s.dil = Setting("DIL", f.dil_value, "Switch for dilution")

    if f.jmeasured_include:
        s.jmeasured_include = Stmt(f"#include {f.jmeasured_include}",
                                   "definitions for further use of jno2 and jo1d as scaling factors")
    for r in f.j_rows:
        if r.mode == "off":
            continue
        if r.mode in ("enz", "input"):
            st = Setting(r.name, f"INPUT({r.column})", r.comment)
        else:
            st = Setting(r.name, r.value if re.match(r"\w+\s*\(", r.value) else f"CONST({r.value})", r.comment)
        s.j_inputs.append(st); s.measured_items.append(st)
        rxn = f.j_overrides.get(r.name, "")
        if rxn:
            rd = ReactionDef(rxn.split("-->")[0].strip(), rxn.split("-->", 1)[1].strip(),
                             rate_raw=f"CONST({r.name})")
            rd._rate_only_text = rxn                # type: ignore[attr-defined]
            s.j_overrides.append(rd); s.measured_items.append(rd)
    for r in f.setting_rows:
        if r.mode == "off":
            continue
        alias = dict((n, a) for n, a, _c, _cm in SETTINGS).get(r.name, "")
        if r.mode == "input":
            st = Setting(r.name, f"INPUT({r.column})", r.comment)
        elif r.mode == "enz":
            if alias:
                st_in = Setting(alias, f"INPUT({r.column})", r.comment)
                s.settings.append(st_in); s.measured_items.append(st_in)
                st = Setting(r.name, SETTING_CONST_FROM_INPUT.get(r.name, f"CONST({alias})"))
            else:
                st = Setting(r.name, f"INPUT({r.column})", r.comment)
        else:
            if alias and r.column:                   # keep the raw ENZ alias for the output
                st_in = Setting(alias, f"INPUT({r.column})", r.comment)
                s.settings.append(st_in); s.measured_items.append(st_in)
            st = Setting(r.name, r.value if re.match(r"\w+\s*\(", r.value) else f"CONST({r.value})",
                         r.comment if not alias else "")
        s.settings.append(st); s.measured_items.append(st)

    s.files = [("ENZ", f.enz_file)]
    s.facs = [("HMAX", f.hmax), ("OUTSTEP", f.outstep)] + list(f.extra_facs)
    s.initials = [(sp, ex, c) for sp, ex, c in f.initials]

    for sp in f.tracer_dilution:
        s.tracer_dilution.append(ReactionDef(f"{sp}+DIL", "", rate_raw="CONST(DILUTE)", text_raw=f"{sp}+DIL-->"))
    if f.dark_hono:
        s.dark_hono = ReactionDef("H2O", f.dark_hono_rhs, rate_raw=f"CONST({f.dark_hono_rate})",
                                  text_raw=f"H2O-->{f.dark_hono_rhs}")
    if f.light_source:
        s.light_sources.append(ReactionDef("SNO2", f.light_source_rhs, rate_raw=f"CONST({f.light_source_rate})",
                                           text_raw=f"SNO2 --> {f.light_source_rhs}"))
        s.initials.append(("SNO2", f"CONST({f.light_source_initial})", ""))

    inj_by_switch: Dict[str, Injection] = {}
    for r in f.injections:
        if not r.switch or not r.column:
            continue
        inj = inj_by_switch.get(r.switch)
        if inj is None:
            inj = Injection(r.switch, r.column, comment=r.comment)
            inj_by_switch[r.switch] = inj; s.injections.append(inj); s.tail.append(inj)
        if r.products.strip() or r.rate.strip():
            inj.channels.append(ReactionDef(r.switch, r.products, rate_raw=_rate_expr(r.rate),
                                            text_raw=f"{r.switch}-->{r.products}"))
    for lhs, rhs, rate in f.chain_reactions:
        rd = ReactionDef(lhs, rhs, rate_raw=_rate_expr(rate, ppb=False), text_raw=f"{lhs}-->{rhs}")
        s.extra_reactions.append(rd); s.tail.append(rd)
    for code in f.extra_statements:
        if code.strip():
            st = Stmt(code.strip()); s.extra.append(st); s.tail.append(st)
    return s


# ─────────────────────────────────────────────────────────────────────────────
# SceneSpec → form (loading an existing scene.easy into the dialog)
# ─────────────────────────────────────────────────────────────────────────────

def _const_inner(rhs: str) -> str:
    m = re.match(r"CONST\s*\((.*)\)\s*$", rhs.strip(), re.I | re.S)
    return m.group(1).strip() if m else rhs.strip()


def form_from_spec(s: SceneSpec) -> EasyForm:
    f = default_form()
    f.title = s.title
    f.header_comments = list(s.header_comments)
    f.includes_pre = [st.code.split(None, 1)[1] for st in s.includes_pre]
    f.use_anorg = bool(s.anorg_includes)
    f.mech_stem = s.mech_stem
    f.use_j_clear_sky = s.j_clear_sky is not None
    f.includes_other = [st.code.split(None, 1)[1] for st in s.includes_other]
    extra: List[str] = []
    # includes we do not regenerate from the stem (e.g. a custom anorg set) → keep verbatim
    std_anorg = {n for n, _ in ANORG_INCLUDES}
    for st in s.anorg_includes:
        if st.code.split(None, 1)[1] not in std_anorg:
            extra.append(st.code)
    if s.mech_stem:
        std_mech = {f"{s.mech_stem}{suf}".lower() for suf, _ in MECH_SUFFIXES}
        for st in s.mech_includes:
            if st.code.split(None, 1)[1].lower() not in std_mech:
                extra.append(st.code)

    f.walldep = s.walldep is not None or s.walldep_include is not None
    if s.walldep:
        f.walldep_value = _const_inner(s.walldep.rhs)
    f.hono_source = s.hono_source is not None
    if s.hono_source:
        f.hono_lhs, f.hono_rhs = s.hono_source.lhs, s.hono_source.rhs
        f.hono_rate = s.hono_source.rate_raw or HONO_RATE
        f.hono_const = _const_inner(s.hono_const.rhs) if s.hono_const else "1."
    f.hcho_source = s.hcho_source is not None
    if s.hcho_source:
        f.hcho_lhs, f.hcho_rhs = s.hcho_source.lhs, s.hcho_source.rhs
        f.hcho_rate = s.hcho_source.rate_raw or HCHO_RATE
        f.hcho_const = _const_inner(s.hcho_const.rhs) if s.hcho_const else "1"
        f.hcho_const_name = s.hcho_const.name if s.hcho_const else f.hcho_lhs.split("+")[0].strip()
    if s.wvc2rh and s.wvc2rh.rhs.replace(" ", "") != WVC2RH_DEF.replace(" ", ""):
        extra.append(s.wvc2rh.render())
    if (s.wvc2rh or s.rh) and not (f.hono_source or f.hcho_source):
        if s.wvc2rh: extra.append(s.wvc2rh.render())
        if s.rh: extra.append(s.rh.render())
    f.bgreac = s.y_reaction is not None
    if s.y_reaction:
        f.y_name = s.y_reaction.lhs.split("+")[0].strip()
        f.y_product = s.y_reaction.rhs
        f.y_rate = s.y_reaction.rate_raw or Y_RATE
    if s.y_const:
        inner = _const_inner(s.y_const.rhs)
        m = re.fullmatch(r"\s*([0-9.]+)e-9\s*\*\s*M\s*", inner, re.I)
        if m:
            f.y_ppb, f.y_const_expr = m.group(1), ""
        else:
            f.y_const_expr = inner
    f.dilution = s.dilute is not None or bool(s.dil_includes) or s.dil is not None
    if s.dilute:
        f.dilute_expr = s.dilute.rhs
    _dstems = {re.sub(r"\.dil\.(mech|rate)$", "", st.code.split(None, 1)[1], flags=re.I) for st in s.dil_includes}
    if len(_dstems) == 1:
        _d = _dstems.pop()
        f.dil_stem = "" if _d.lower() == (s.mech_stem or "").lower() else _d
    elif _dstems:
        extra.extend(st.code for st in s.dil_includes)
    f.dil_flow_const = _const_inner(s.dil_flow.rhs) if s.dil_flow else ""
    f.dil_value = s.dil.rhs if s.dil else ""
    if s.dil is None:
        f.dil_value = ""
    f.jmeasured_include = s.jmeasured_include.code.split(None, 1)[1] if s.jmeasured_include else ""

    # measured: rebuild rows from the encounter-ordered items
    f.j_rows, f.setting_rows, f.j_overrides = [], [], {}
    alias_of = {a: n for n, a, _c, _cm in SETTINGS if a}
    pending_alias: Dict[str, Setting] = {}
    names_seen = set()
    for item in s.measured_items:
        if isinstance(item, Setting):
            up = item.name.upper()
            m_in = re.match(r"INPUT\s*\(\s*(\w+)\s*\)\s*$", item.rhs, re.I)
            if up.startswith("J") and up not in {"J"}:
                if m_in:
                    f.j_rows.append(MeasuredRow(item.name, "enz", m_in.group(1), comment=item.comment))
                else:
                    f.j_rows.append(MeasuredRow(item.name, "const", "", item.rhs, item.comment))
                names_seen.add(up)
            elif up in alias_of:                      # T_in / HV_in / H2O_in
                pending_alias[alias_of[up]] = item
            elif up in {n for n, _a, _c, _cm in SETTINGS}:
                name = next(n for n, _a, _c, _cm in SETTINGS if n == up)
                al = pending_alias.pop(name, None)
                if m_in:
                    f.setting_rows.append(MeasuredRow(name, "input", m_in.group(1), comment=item.comment))
                    if al is not None:
                        extra.append(al.render())
                elif al is not None and item.rhs.replace(" ", "") == SETTING_CONST_FROM_INPUT.get(name, "").replace(" ", ""):
                    col = re.match(r"INPUT\s*\(\s*(\w+)\s*\)", al.rhs, re.I).group(1)
                    f.setting_rows.append(MeasuredRow(name, "enz", col, comment=al.comment))
                else:
                    col = re.match(r"INPUT\s*\(\s*(\w+)\s*\)", al.rhs, re.I).group(1) if al is not None else ""
                    f.setting_rows.append(MeasuredRow(name, "const", col, item.rhs, item.comment or (al.comment if al else "")))
                names_seen.add(up)
            else:
                extra.append(item.render())
        else:                                          # rate-only J override
            m = re.match(r"CONST\s*\(\s*(\w+)\s*\)", item.rate_raw or "", re.I)
            if m:
                f.j_overrides[m.group(1)] = item._rate_only_text
            else:
                extra.append(f"k[{item._rate_only_text}] = {item.rate_raw}")
    for name, al in pending_alias.items():             # alias without a main definition
        extra.append(al.render())
    # j overrides keyed case-insensitively to the j row names
    jmap = {k.upper(): v for k, v in f.j_overrides.items()}
    f.j_overrides = {r.name: jmap.get(r.name.upper(), "") for r in f.j_rows}

    f.enz_file = next((v for k, v in s.files if k == "ENZ"), "")
    f.extra_facs = [(k, v) for k, v in s.facs if k not in ("HMAX", "OUTSTEP")]
    f.hmax = next((v for k, v in s.facs if k == "HMAX"), "")
    f.outstep = next((v for k, v in s.facs if k == "OUTSTEP"), "")
    for k, v in s.files:
        if k != "ENZ":
            extra.append(f"FILES[{k}] = {v}")

    f.initials = [(sp, ex, c) for sp, ex, c in s.initials if sp.upper() != "SNO2" or not s.light_sources]
    f.tracer_dilution = [rd.lhs.split("+")[0].strip() for rd in s.tracer_dilution]
    f.dark_hono = s.dark_hono is not None
    if s.dark_hono:
        f.dark_hono_rhs = s.dark_hono.rhs
        f.dark_hono_rate = _const_inner(s.dark_hono.rate_raw or "")
    f.light_source = bool(s.light_sources)
    if s.light_sources:
        ls = s.light_sources[0]
        f.light_source_rhs, f.light_source_rate = ls.rhs, _const_inner(ls.rate_raw or "")
        init = next((ex for sp, ex, _c in s.initials if sp.upper() == "SNO2"), "")
        f.light_source_initial = _const_inner(init) if init else "0"
        for ls2 in s.light_sources[1:]:
            extra.extend(ls2.render())

    f.injections, f.chain_reactions = [], []
    for item in s.tail:
        if isinstance(item, Injection):
            if not item.channels:
                f.injections.append(InjectionRow(item.switch, item.enz_column, "", "", item.comment))
            for rd in item.channels:
                f.injections.append(InjectionRow(item.switch, item.enz_column, rd.rhs,
                                                 _ppb_from_rate(rd.rate_raw), item.comment or rd.comment))
            for rd in item.chain:
                f.chain_reactions.append((rd.lhs, rd.rhs, _ppb_from_rate(rd.rate_raw, ppb=False)))
        elif isinstance(item, tuple):
            extra.append(item[1].code)
        elif isinstance(item, Stmt):
            extra.append(item.code)
        elif hasattr(item, "_rate_only_text"):
            extra.append(f"k[{item._rate_only_text}] = {item.rate_raw}")
        else:
            if item.arrow == "-->" and item.rate_raw:
                f.chain_reactions.append((item.lhs, item.rhs, _ppb_from_rate(item.rate_raw, ppb=False)))
            else:
                extra.extend(item.render())
    f.extra_statements = extra
    return f


def load_form(easy_path: str) -> EasyForm:
    return form_from_spec(extract_scene(parse_document(easy_path)))


def render_form(f: EasyForm) -> str:
    return write_scene(spec_from_form(f))


# ─────────────────────────────────────────────────────────────────────────────
# Qt dialog
# ─────────────────────────────────────────────────────────────────────────────

try:
    from PyQt6.QtCore import pyqtSignal, Qt
    from PyQt6.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QGridLayout, QLabel, QLineEdit,
                                 QPushButton, QCheckBox, QComboBox, QTabWidget, QWidget, QTableWidget,
                                 QTableWidgetItem, QPlainTextEdit, QFileDialog, QMessageBox, QGroupBox,
                                 QHeaderView, QSplitter)
    _QT = True
except ImportError:                                   # pure-function use without Qt
    _QT = False


if _QT:
    class EasyBuilderDialog(QDialog):
        """Step 3: build scene.easy from mechanism + ENZ + chamber options + injections."""
        easy_written = pyqtSignal(str)

        def __init__(self, parent=None, include_dir: str = "", easy_root: str = "",
                     enz_path: str = "", last_easy_dir: str = ""):
            super().__init__(parent)
            self.setWindowTitle("3. Build EASY scenario file")
            self.resize(1400, 900)
            self.include_dir, self.easy_root = include_dir, easy_root
            self.last_easy_dir = last_easy_dir or easy_root
            self.enz_path = enz_path
            self.enz_columns: List[str] = read_enz_columns(enz_path) if enz_path else []
            self.mechs = scan_mechanisms(include_dir)
            self._build_ui()
            self.set_form(default_form())
            if enz_path:
                self.enz_edit.setText(enz_path)
                self._enz_changed()

        # ── UI ───────────────────────────────────────────────────────────────
        def _build_ui(self):
            root = QVBoxLayout(self)
            top = QGridLayout()
            top.addWidget(QLabel("Title"), 0, 0)
            self.title_edit = QLineEdit(); top.addWidget(self.title_edit, 0, 1, 1, 5)
            top.addWidget(QLabel("Mechanism"), 1, 0)
            self.mech_combo = QComboBox(); self.mech_combo.setEditable(True)
            self.mech_combo.addItems(list(self.mechs)); top.addWidget(self.mech_combo, 1, 1)
            self.anorg_chk = QCheckBox("anorg.mcm331 includes"); self.anorg_chk.setChecked(True); top.addWidget(self.anorg_chk, 1, 2)
            self.jcs_chk = QCheckBox("j_clear_sky.incl"); self.jcs_chk.setChecked(True); top.addWidget(self.jcs_chk, 1, 3)
            top.addWidget(QLabel("ENZ file"), 2, 0)
            self.enz_edit = QLineEdit(); top.addWidget(self.enz_edit, 2, 1, 1, 3)
            b = QPushButton("…"); b.clicked.connect(self._browse_enz); top.addWidget(b, 2, 4)
            b2 = QPushButton("Load existing .easy …"); b2.clicked.connect(self._load_easy); top.addWidget(b2, 2, 5)
            root.addLayout(top)

            tabs = QTabWidget(); root.addWidget(tabs, stretch=1)
            tabs.addTab(self._tab_chamber(), "Chamber")
            tabs.addTab(self._tab_measured(), "Measured inputs")
            tabs.addTab(self._tab_initials(), "Initials && FACSIMILE")
            tabs.addTab(self._tab_injections(), "Injections")
            tabs.addTab(self._tab_extra(), "Extra && Preview")
            self.tabs = tabs

            bottom = QHBoxLayout()
            self.out_edit = QLineEdit(); self.out_edit.setPlaceholderText("output scene.easy path")
            bottom.addWidget(QLabel("Write to")); bottom.addWidget(self.out_edit, stretch=1)
            b = QPushButton("…"); b.clicked.connect(self._browse_out); bottom.addWidget(b)
            pv = QPushButton("Preview"); pv.clicked.connect(self._preview); bottom.addWidget(pv)
            wr = QPushButton("Write .easy"); wr.clicked.connect(lambda: self._write(open_after=False)); bottom.addWidget(wr)
            wo = QPushButton("Write && open in Step 2"); wo.setStyleSheet("font-weight: bold;")
            wo.clicked.connect(lambda: self._write(open_after=True)); bottom.addWidget(wo)
            root.addLayout(bottom)

        def _tab_chamber(self):
            w = QWidget(); g = QGridLayout(w); r = 0
            def row(label, *widgets):
                nonlocal r
                g.addWidget(QLabel(label), r, 0)
                for i, wd in enumerate(widgets, 1):
                    g.addWidget(wd, r, i)
                r += 1
            self.walldep_chk = QCheckBox("saphir.walldep.incl + WALLDEP = CONST("); self.walldep_val = QLineEdit("1.")
            row("Wall deposition", self.walldep_chk, self.walldep_val, QLabel(")"))
            self.hono_chk = QCheckBox("HONO source"); self.hono_lhs = QLineEdit(); self.hono_rhs = QLineEdit()
            self.hono_const = QLineEdit(); self.hono_rate = QLineEdit()
            row("HONO (Rohrer et al.)", self.hono_chk, self.hono_lhs, QLabel("-->"), self.hono_rhs, QLabel("SHONO = CONST("), self.hono_const)
            row("   k[HONO]", self.hono_rate)
            self.hcho_chk = QCheckBox("HCHO source"); self.hcho_lhs = QLineEdit(); self.hcho_rhs = QLineEdit()
            self.hcho_const = QLineEdit(); self.hcho_rate = QLineEdit(); self.hcho_cname = QLineEdit("SHCHO")
            row("HCHO", self.hcho_chk, self.hcho_lhs, QLabel("-->"), self.hcho_rhs, self.hcho_cname, QLabel("= CONST("), self.hcho_const)
            row("   k[HCHO]", self.hcho_rate)
            self.y_chk = QCheckBox("Background reactivity"); self.y_name = QLineEdit("qY"); self.y_prod = QLineEdit("CH3O2")
            self.y_ppb = QLineEdit("250"); self.y_rate = QLineEdit(Y_RATE)
            row("Y + OH -->", self.y_chk, self.y_name, self.y_prod, QLabel("Y [ppb CO-eq]"), self.y_ppb, self.y_rate)
            self.dil_chk = QCheckBox("dilution includes + DILUTE"); self.dilute_expr = QLineEdit(DILUTE_DEF)
            self.dil_flow = QLineEdit(); self.dil_flow.setPlaceholderText("FLOW = CONST(…) in DIL block (optional)")
            self.dil_val = QLineEdit("CONST(1.)"); self.dil_stem = QLineEdit(); self.dil_stem.setPlaceholderText(".dil stem (default: mechanism)")
            row("Dilution", self.dil_chk, self.dilute_expr, self.dil_flow, QLabel("DIL ="), self.dil_val, self.dil_stem)
            self.tracer_edit = QLineEdit("XNOy, NA")
            row("Tracer dilution (X+DIL-->)", self.tracer_edit)
            self.dark_chk = QCheckBox("H2O --> "); self.dark_rhs = QLineEdit("HONO+XNOy"); self.dark_rate = QLineEdit("1e-012/10")
            row("Dark HONO", self.dark_chk, self.dark_rhs, QLabel("k = CONST("), self.dark_rate, QLabel(")"))
            self.light_chk = QCheckBox("SNO2 --> "); self.light_rhs = QLineEdit("NO2+XNOy")
            self.light_rate = QLineEdit(); self.light_init = QLineEdit()
            row("NO2 light source", self.light_chk, self.light_rhs, QLabel("k = CONST("), self.light_rate, QLabel("INITIAL[SNO2] = CONST("), self.light_init)
            g.setRowStretch(r, 1)
            return w

        def _tab_measured(self):
            w = QWidget(); v = QVBoxLayout(w)
            h = QHBoxLayout(); h.addWidget(QLabel("J-measured include")); self.jmeas_edit = QLineEdit("saphir.jmeasured_02.incl"); h.addWidget(self.jmeas_edit); v.addLayout(h)
            self.j_table = QTableWidget(0, 5); self.j_table.setHorizontalHeaderLabels(["Variable", "Mode", "ENZ column", "Constant / expression", "k[] override (reaction)"])
            self.j_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
            v.addWidget(QLabel("Photolysis frequencies")); v.addWidget(self.j_table)
            self.s_table = QTableWidget(0, 5); self.s_table.setHorizontalHeaderLabels(["Setting", "Mode (enz | input | const | off)", "ENZ column", "Constant / expression", "Comment"])
            self.s_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
            v.addWidget(QLabel("Physical settings (T, P, SZ, HV, FLOW, H2O)")); v.addWidget(self.s_table)
            return w

        def _tab_initials(self):
            w = QWidget(); g = QGridLayout(w)
            g.addWidget(QLabel("FILES[ENZ]"), 0, 0); self.files_enz = QLineEdit(); g.addWidget(self.files_enz, 0, 1)
            g.addWidget(QLabel("FACS[HMAX]"), 0, 2); self.hmax_edit = QLineEdit("30"); g.addWidget(self.hmax_edit, 0, 3)
            g.addWidget(QLabel("FACS[OUTSTEP]"), 0, 4); self.outstep_edit = QLineEdit("100000"); g.addWidget(self.outstep_edit, 0, 5)
            self.init_table = QTableWidget(0, 3); self.init_table.setHorizontalHeaderLabels(["Species ('*' = all)", "Expression (e.g. M*1e-9*14.0)", "Comment"])
            self.init_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
            g.addWidget(self.init_table, 1, 0, 1, 6)
            hb = QHBoxLayout(); a = QPushButton("+ row"); a.clicked.connect(lambda: self._add_row(self.init_table)); d = QPushButton("− row"); d.clicked.connect(lambda: self._del_row(self.init_table)); hb.addWidget(a); hb.addWidget(d); hb.addStretch()
            g.addLayout(hb, 2, 0, 1, 6)
            return w

        def _tab_injections(self):
            w = QWidget(); v = QVBoxLayout(w)
            v.addWidget(QLabel("One row per channel: switch = input(ENZ column); switch --> products with k = CONST(M*1e-9/60*ppb_per_min). "
                               "Rate accepts a number (ppb/min) or a full expression."))
            self.inj_table = QTableWidget(0, 5); self.inj_table.setHorizontalHeaderLabels(["Switch", "ENZ column", "Products", "Rate (ppb/min or expr)", "Comment"])
            self.inj_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
            v.addWidget(self.inj_table)
            hb = QHBoxLayout()
            a = QPushButton("+ row"); a.clicked.connect(lambda: self._add_row(self.inj_table)); hb.addWidget(a)
            d = QPushButton("− row"); d.clicked.connect(lambda: self._del_row(self.inj_table)); hb.addWidget(d)
            f = QPushButton("Add rows from ENZ Q-columns"); f.clicked.connect(self._fill_from_enz); hb.addWidget(f)
            hb.addStretch(); v.addLayout(hb)
            v.addWidget(QLabel("Follow-up reactions (e.g. TwoReaction chain VOCin1 --> VOCin_all, k = CONST(1./60))"))
            self.chain_table = QTableWidget(0, 3); self.chain_table.setHorizontalHeaderLabels(["Reactants", "Products", "Rate expression"])
            self.chain_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
            v.addWidget(self.chain_table)
            hb2 = QHBoxLayout(); a2 = QPushButton("+ row"); a2.clicked.connect(lambda: self._add_row(self.chain_table)); d2 = QPushButton("− row"); d2.clicked.connect(lambda: self._del_row(self.chain_table)); hb2.addWidget(a2); hb2.addWidget(d2); hb2.addStretch(); v.addLayout(hb2)
            return w

        def _tab_extra(self):
            w = QWidget(); v = QVBoxLayout(w)
            sp = QSplitter(Qt.Orientation.Horizontal)
            left = QWidget(); lv = QVBoxLayout(left); lv.addWidget(QLabel("Extra statements (verbatim, one per line)"))
            self.extra_edit = QPlainTextEdit(); lv.addWidget(self.extra_edit)
            lv.addWidget(QLabel("Extra includes (one per line, e.g. VOC_Quellen_CHANEL2024.incl)"))
            self.inc_other_edit = QPlainTextEdit(); self.inc_other_edit.setMaximumHeight(90); lv.addWidget(self.inc_other_edit)
            sp.addWidget(left)
            right = QWidget(); rv = QVBoxLayout(right); rv.addWidget(QLabel("Preview of scene.easy"))
            self.preview = QPlainTextEdit(); self.preview.setReadOnly(True); rv.addWidget(self.preview)
            sp.addWidget(right); sp.setSizes([500, 900])
            v.addWidget(sp)
            return w

        # ── table helpers ────────────────────────────────────────────────────
        @staticmethod
        def _set_rows(table: QTableWidget, rows: List[List[str]]):
            table.setRowCount(0)
            for vals in rows:
                r = table.rowCount(); table.insertRow(r)
                for c, v in enumerate(vals):
                    table.setItem(r, c, QTableWidgetItem(str(v)))

        @staticmethod
        def _get_rows(table: QTableWidget) -> List[List[str]]:
            out = []
            for r in range(table.rowCount()):
                vals = []
                for c in range(table.columnCount()):
                    it = table.item(r, c); vals.append(it.text().strip() if it else "")
                if any(vals):
                    out.append(vals)
            return out

        @staticmethod
        def _add_row(table: QTableWidget):
            table.insertRow(table.rowCount())

        @staticmethod
        def _del_row(table: QTableWidget):
            r = table.currentRow()
            if r >= 0:
                table.removeRow(r)

        # ── form ↔ widgets ───────────────────────────────────────────────────
        def set_form(self, f: EasyForm):
            self.title_edit.setText(f.title)
            self.mech_combo.setCurrentText(f.mech_stem)
            self.anorg_chk.setChecked(f.use_anorg); self.jcs_chk.setChecked(f.use_j_clear_sky)
            self.walldep_chk.setChecked(f.walldep); self.walldep_val.setText(f.walldep_value)
            self.hono_chk.setChecked(f.hono_source); self.hono_lhs.setText(f.hono_lhs); self.hono_rhs.setText(f.hono_rhs)
            self.hono_const.setText(f.hono_const); self.hono_rate.setText(f.hono_rate)
            self.hcho_chk.setChecked(f.hcho_source); self.hcho_lhs.setText(f.hcho_lhs); self.hcho_rhs.setText(f.hcho_rhs)
            self.hcho_const.setText(f.hcho_const); self.hcho_rate.setText(f.hcho_rate); self.hcho_cname.setText(f.hcho_const_name)
            self.y_chk.setChecked(f.bgreac); self.y_name.setText(f.y_name); self.y_prod.setText(f.y_product)
            self.y_ppb.setText(f.y_const_expr or f.y_ppb); self.y_rate.setText(f.y_rate)
            self.dil_chk.setChecked(f.dilution); self.dilute_expr.setText(f.dilute_expr); self.dil_flow.setText(f.dil_flow_const); self.dil_val.setText(f.dil_value); self.dil_stem.setText(f.dil_stem)
            self.tracer_edit.setText(", ".join(f.tracer_dilution))
            self.dark_chk.setChecked(f.dark_hono); self.dark_rhs.setText(f.dark_hono_rhs); self.dark_rate.setText(f.dark_hono_rate)
            self.light_chk.setChecked(f.light_source); self.light_rhs.setText(f.light_source_rhs); self.light_rate.setText(f.light_source_rate); self.light_init.setText(f.light_source_initial)
            self.jmeas_edit.setText(f.jmeasured_include)
            self._set_rows(self.j_table, [[r.name, r.mode, r.column, r.value, f.j_overrides.get(r.name, "")] for r in f.j_rows])
            self._set_rows(self.s_table, [[r.name, r.mode, r.column, r.value, r.comment] for r in f.setting_rows])
            self.files_enz.setText(f.enz_file); self.hmax_edit.setText(f.hmax); self.outstep_edit.setText(f.outstep)
            self._set_rows(self.init_table, [[sp, ex, c] for sp, ex, c in f.initials])
            self._set_rows(self.inj_table, [[r.switch, r.column, r.products, r.rate, r.comment] for r in f.injections])
            self._set_rows(self.chain_table, [[a, b, c] for a, b, c in f.chain_reactions])
            self.extra_edit.setPlainText("\n".join(f.extra_statements + [f"FACS[{k}] = {v}" for k, v in f.extra_facs]))
            self.inc_other_edit.setPlainText("\n".join(f.includes_pre + f.includes_other))
            self._form_header = list(f.header_comments)

        def get_form(self) -> EasyForm:
            f = default_form()
            f.title = self.title_edit.text().strip()
            f.header_comments = list(getattr(self, "_form_header", []))
            f.mech_stem = self.mech_combo.currentText().strip()
            f.use_anorg = self.anorg_chk.isChecked(); f.use_j_clear_sky = self.jcs_chk.isChecked()
            f.walldep = self.walldep_chk.isChecked(); f.walldep_value = self.walldep_val.text().strip()
            f.hono_source = self.hono_chk.isChecked(); f.hono_lhs = self.hono_lhs.text().strip(); f.hono_rhs = self.hono_rhs.text().strip()
            f.hono_const = self.hono_const.text().strip(); f.hono_rate = self.hono_rate.text().strip()
            f.hcho_source = self.hcho_chk.isChecked(); f.hcho_lhs = self.hcho_lhs.text().strip(); f.hcho_rhs = self.hcho_rhs.text().strip()
            f.hcho_const = self.hcho_const.text().strip(); f.hcho_rate = self.hcho_rate.text().strip(); f.hcho_const_name = self.hcho_cname.text().strip()
            f.bgreac = self.y_chk.isChecked(); f.y_name = self.y_name.text().strip(); f.y_product = self.y_prod.text().strip()
            yv = self.y_ppb.text().strip()
            if re.fullmatch(r"[0-9.]+", yv):
                f.y_ppb, f.y_const_expr = yv, ""
            else:
                f.y_const_expr = yv
            f.y_rate = self.y_rate.text().strip()
            f.dilution = self.dil_chk.isChecked(); f.dilute_expr = self.dilute_expr.text().strip()
            f.dil_flow_const = self.dil_flow.text().strip(); f.dil_value = self.dil_val.text().strip(); f.dil_stem = self.dil_stem.text().strip()
            f.tracer_dilution = [t.strip() for t in self.tracer_edit.text().split(",") if t.strip()]
            f.dark_hono = self.dark_chk.isChecked(); f.dark_hono_rhs = self.dark_rhs.text().strip(); f.dark_hono_rate = self.dark_rate.text().strip()
            f.light_source = self.light_chk.isChecked(); f.light_source_rhs = self.light_rhs.text().strip()
            f.light_source_rate = self.light_rate.text().strip(); f.light_source_initial = self.light_init.text().strip()
            f.jmeasured_include = self.jmeas_edit.text().strip()
            f.j_rows, f.j_overrides = [], {}
            for name, mode, col, val, ov in self._get_rows(self.j_table):
                f.j_rows.append(MeasuredRow(name, mode or "enz", col, val)); f.j_overrides[name] = ov
            f.setting_rows = [MeasuredRow(n, m or "enz", c, v, cm) for n, m, c, v, cm in self._get_rows(self.s_table)]
            f.enz_file = self.files_enz.text().strip(); f.hmax = self.hmax_edit.text().strip(); f.outstep = self.outstep_edit.text().strip()
            f.initials = [(sp, ex, c) for sp, ex, c in self._get_rows(self.init_table)]
            f.injections = [InjectionRow(sw, col, pr, rt, cm) for sw, col, pr, rt, cm in self._get_rows(self.inj_table)]
            f.chain_reactions = [(a, b, c) for a, b, c in self._get_rows(self.chain_table)]
            f.extra_facs, f.extra_statements = [], []
            for line in self.extra_edit.toPlainText().splitlines():
                m = re.match(r"FACS\s*\[\s*(\w+)\s*\]\s*=\s*(.+)$", line.strip(), re.I)
                if m:
                    f.extra_facs.append((m.group(1).upper(), m.group(2).strip()))
                elif line.strip():
                    f.extra_statements.append(line.strip())
            incs = [l.strip() for l in self.inc_other_edit.toPlainText().splitlines() if l.strip()]
            f.includes_pre = [i for i in incs if i.lower().startswith("default_")]
            f.includes_other = [i for i in incs if not i.lower().startswith("default_")]
            return f

        # ── actions ──────────────────────────────────────────────────────────
        def _browse_enz(self):
            p, _ = QFileDialog.getOpenFileName(self, "Select ENZ file", self.last_easy_dir, "ENZ files (*.enz);;All files (*)")
            if p:
                self.enz_edit.setText(p); self._enz_changed()

        def _enz_changed(self):
            self.enz_path = self.enz_edit.text().strip()
            self.enz_columns = read_enz_columns(self.enz_path) if self.enz_path else []
            if self.enz_path:
                self.files_enz.setText(os.path.basename(self.enz_path))
                if not self.out_edit.text():
                    stem = os.path.splitext(os.path.basename(self.enz_path))[0]
                    self.out_edit.setText(os.path.join(os.path.dirname(self.enz_path), stem + ".easy"))

        def _fill_from_enz(self):
            """Add one injection row per Q*/q* ENZ column that is not yet used."""
            used = {r[1].upper() for r in self._get_rows(self.inj_table)}
            for col in self.enz_columns:
                if col.upper().startswith("Q") and col.upper() not in used:
                    sw = col if col.lower().startswith("q_") else "q" + col[1:]
                    r = self.inj_table.rowCount(); self.inj_table.insertRow(r)
                    for c, v in enumerate([sw, col, "", "", "from ENZ"]):
                        self.inj_table.setItem(r, c, QTableWidgetItem(v))

        def _load_easy(self):
            p, _ = QFileDialog.getOpenFileName(self, "Load existing EASY file", self.last_easy_dir, "EASY files (*.easy);;All files (*)")
            if not p:
                return
            try:
                self.set_form(load_form(p))
            except Exception as e:
                QMessageBox.critical(self, "Load failed", str(e)); return
            self.last_easy_dir = os.path.dirname(p)
            enz = os.path.join(os.path.dirname(p), self.files_enz.text())
            if os.path.isfile(enz):
                self.enz_edit.setText(enz); self._enz_changed()
            self.out_edit.setText(os.path.splitext(p)[0] + ".easy")
            self._preview()

        def _browse_out(self):
            p, _ = QFileDialog.getSaveFileName(self, "Write EASY file", self.out_edit.text() or self.last_easy_dir, "EASY files (*.easy)")
            if p:
                self.out_edit.setText(p)

        def _preview(self):
            try:
                self.preview.setPlainText(render_form(self.get_form()))
                self.tabs.setCurrentIndex(self.tabs.count() - 1)
            except Exception as e:
                QMessageBox.critical(self, "Preview failed", str(e))

        def _write(self, open_after: bool):
            out = self.out_edit.text().strip()
            if not out:
                QMessageBox.warning(self, "No output path", "Choose where to write the .easy file."); return
            try:
                text = render_form(self.get_form())
                os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
                with open(out, "w", encoding="latin-1", errors="replace", newline="\n") as fh:
                    fh.write(text)
                enz = self.enz_edit.text().strip()
                if enz and os.path.isfile(enz):
                    dest = os.path.join(os.path.dirname(out), os.path.basename(enz))
                    if os.path.abspath(dest) != os.path.abspath(enz):
                        shutil.copy2(enz, dest)
            except Exception as e:
                QMessageBox.critical(self, "Write failed", str(e)); return
            self.preview.setPlainText(text)
            if open_after:
                self.easy_written.emit(out)
                self.accept()
            else:
                QMessageBox.information(self, "Written", f"EASY file written:\n{out}")


# ─────────────────────────────────────────────────────────────────────────────
# self-test: load every scene → form → write → statement-equivalent?
# ─────────────────────────────────────────────────────────────────────────────

def _selftest(scenes: List[str]) -> int:
    import tempfile
    here = os.path.dirname(os.path.abspath(__file__))
    idl = os.path.join(here, "MCM_tested_scenarios", "IDL Input")
    scenes = scenes or sorted(d for d in os.listdir(idl) if d.startswith("scene_"))
    tmp = tempfile.mkdtemp(prefix="easy_builder_")
    bad = 0
    for sc in scenes:
        p = os.path.join(idl, sc, sc + ".easy")
        if not os.path.isfile(p):
            continue
        form = load_form(p)
        text = render_form(form)
        out = os.path.join(tmp, sc + ".easy")
        with open(out, "w", encoding="latin-1", errors="replace", newline="\n") as fh:
            fh.write(text)
        ok, diffs = statements_equivalent(p, out)
        bad += not ok
        print(f"{'OK  ' if ok else 'FAIL'} {sc}: injections={len(form.injections)} chain={len(form.chain_reactions)} "
              f"extra={len(form.extra_statements)}")
        for d in diffs[:6]:
            print("       " + d[:150])
    print(f"{len(scenes) - bad}/{len(scenes)} pass  (files in {tmp})")
    return 1 if bad else 0


if __name__ == "__main__":
    if "--selftest" in sys.argv:
        sys.exit(_selftest([a for a in sys.argv[1:] if not a.startswith("--")]))
    if _QT:
        from PyQt6.QtWidgets import QApplication
        app = QApplication(sys.argv)
        dlg = EasyBuilderDialog(include_dir=sys.argv[1] if len(sys.argv) > 1 else "")
        dlg.show()
        sys.exit(app.exec())
