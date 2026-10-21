"""
injection_calc.py — how much ended up in the chamber, and what the MCM calls it.

Three things the experiment wizard needs that are chamber- and lab-specific:

1. Chamber profiles   (`CHAMBERS`, `ChamberProfile`)  — volume and the standard
   conditions used to turn an injected amount into a mixing ratio.  SAPHIR is the
   only profile today; the IDL/FACSIMILE scripts use 276 m³ (DILUTE formula), so the
   same number is used here for consistency.  Other chambers: add a profile.

2. Delivery calculations (`resolve_delivery`) — the operator notes one of
       ppb          "7.5 ppb"                       expected mixing ratio, used as is
       liquid       "20 µl isoprene by syringe"     µl × density / molar mass → mol
       gas_syringe  "2×3.2 ml from the 1 % bottle"  ml of bottle gas × fraction → mol
       mfc          "50 s at 5 slpm, 100 ppm bottle" sccm × s × fraction → mol
   Each returns (ppb, explanation) so the wizard can SHOW how the number was derived
   and the user can overrule it.  Moles of air = P·V/(R·T) with the chamber's T, P.

3. Compound resolution and mixtures (`resolve_compound`, `Mixture`, `mixture_event`) —
   solutions of several compounds (VCP, cooking, …) are injected as ONE pulse; every
   component is looked up in the MCM (common-name table → MCM name; else formula;
   else the name as typed).  Components without an MCM species need a surrogate
   (user decision, stored in the mixture file) or stay unmodelled with an alert.
   The result is one InjectionEvent whose product string carries the mole fractions
   (`0.31*LIMONENE+0.42*C2H5OH+…`, the IDL convention for VCP2) and whose total ppb
   is the modelled sum.  Mixture recipes live in `mixtures/*.mixture.json`.

Nothing here needs Qt.
"""
from __future__ import annotations

import json
import os
import re
from dataclasses import dataclass, field, asdict
from typing import Dict, List, Optional, Tuple

R_GAS = 8.314462618            # J/(mol K)

# ─────────────────────────────────────────────────────────────────────────────
# Chamber profiles
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class ChamberProfile:
    name: str
    volume_m3: float
    note: str = ""
    mfc_std_T_K: float = 273.15          # "standard" of the mass-flow controllers (sccm / slpm)
    mfc_std_P_hPa: float = 1013.25
    syringe_T_K: float = 293.15          # gas syringes are filled at lab temperature


CHAMBERS: Dict[str, ChamberProfile] = {
    "SAPHIR": ChamberProfile("SAPHIR", 276.0,
                             "Jülich SAPHIR, 270 m³ nominal; 276 m³ is the value in the IDL DILUTE formula "
                             "and is kept for consistency with the FACSIMILE scenarios."),
}
DEFAULT_CHAMBER = "SAPHIR"


def chamber(name: str = "") -> ChamberProfile:
    return CHAMBERS.get(name or DEFAULT_CHAMBER, CHAMBERS[DEFAULT_CHAMBER])


def moles_air(volume_m3: float, T_K: float, P_hPa: float) -> float:
    return P_hPa * 100.0 * volume_m3 / (R_GAS * T_K)


# ─────────────────────────────────────────────────────────────────────────────
# Delivery calculations — every function returns (ppb, explanation)
# ─────────────────────────────────────────────────────────────────────────────

def ppb_liquid(vol_ul: float, mw: float, density: float, T_K: float, P_hPa: float, ch: ChamberProfile) -> Tuple[float, str]:
    n = vol_ul * 1e-3 * density / mw
    n_air = moles_air(ch.volume_m3, T_K, P_hPa)
    ppb = n / n_air * 1e9
    return ppb, (f"liquid syringe: {vol_ul:g} µl × {density:g} g/ml ÷ {mw:g} g/mol = {n * 1e6:.3g} µmol; "
                 f"air in {ch.volume_m3:g} m³ at {T_K:g} K, {P_hPa:g} hPa = {n_air:.4g} mol → {ppb:.3g} ppb")


def ppb_gas_syringe(vol_ml: float, fraction: float, T_K: float, P_hPa: float, ch: ChamberProfile) -> Tuple[float, str]:
    """`fraction` = mole fraction in the bottle (1 for the pure gas, 1e-6 per ppm)."""
    n = P_hPa * 100.0 * vol_ml * 1e-6 / (R_GAS * ch.syringe_T_K) * fraction
    n_air = moles_air(ch.volume_m3, T_K, P_hPa)
    ppb = n / n_air * 1e9
    return ppb, (f"gas syringe: {vol_ml:g} ml of bottle gas (fraction {fraction:g}) at {ch.syringe_T_K:g} K, {P_hPa:g} hPa "
                 f"= {n * 1e6:.3g} µmol; air {n_air:.4g} mol → {ppb:.3g} ppb")


def ppb_mfc(flow_sccm: float, duration_s: float, fraction: float, T_K: float, P_hPa: float, ch: ChamberProfile) -> Tuple[float, str]:
    v_std_cm3 = flow_sccm * duration_s / 60.0
    n = ch.mfc_std_P_hPa * 100.0 * v_std_cm3 * 1e-6 / (R_GAS * ch.mfc_std_T_K) * fraction
    n_air = moles_air(ch.volume_m3, T_K, P_hPa)
    ppb = n / n_air * 1e9
    return ppb, (f"MFC: {flow_sccm:g} sccm × {duration_s:g} s = {v_std_cm3:g} std-cm³ "
                 f"({ch.mfc_std_T_K:g} K, {ch.mfc_std_P_hPa:g} hPa) × fraction {fraction:g} = {n * 1e6:.3g} µmol; "
                 f"air {n_air:.4g} mol → {ppb:.3g} ppb")


def bottle_fraction(value: float, unit: str) -> float:
    """Bottle concentration → mole fraction.  unit: ppm | ppb | % | fraction | pure."""
    u = (unit or "fraction").lower()
    if u == "ppm":
        return value * 1e-6
    if u == "ppb":
        return value * 1e-9
    if u in ("%", "percent"):
        return value / 100.0
    if u == "pure":
        return 1.0
    return value


def resolve_delivery(d: Dict[str, str], T_K: float, P_hPa: float, ch: ChamberProfile) -> Tuple[Optional[float], str]:
    """Turn an InjectionEvent.delivery dict into (ppb, explanation).  Keys by mode:
         ppb:          value
         liquid:       vol_ul, mw, density
         gas_syringe:  vol_ml, bottle_value, bottle_unit
         mfc:          flow_sccm (or flow_slpm), duration_s, bottle_value, bottle_unit
       Returns (None, reason) when a field is missing."""
    mode = (d.get("mode") or "ppb").lower()
    try:
        if mode == "ppb":
            v = float(d["value"]); return v, f"{v:g} ppb as stated"
        if mode == "liquid":
            return ppb_liquid(float(d["vol_ul"]), float(d["mw"]), float(d["density"]), T_K, P_hPa, ch)
        if mode == "gas_syringe":
            return ppb_gas_syringe(float(d["vol_ml"]), bottle_fraction(float(d.get("bottle_value", 1)), d.get("bottle_unit", "pure")),
                                   T_K, P_hPa, ch)
        if mode == "mfc":
            flow = float(d["flow_sccm"]) if d.get("flow_sccm") else float(d["flow_slpm"]) * 1000.0
            return ppb_mfc(flow, float(d["duration_s"]), bottle_fraction(float(d.get("bottle_value", 1)), d.get("bottle_unit", "pure")),
                           T_K, P_hPa, ch)
    except (KeyError, ValueError, TypeError) as e:
        return None, f"{mode}: missing or invalid field ({e})"
    return None, f"unknown delivery mode '{mode}'"


# ─────────────────────────────────────────────────────────────────────────────
# Compounds: common name → MCM name, molar mass, liquid density
# ─────────────────────────────────────────────────────────────────────────────

# lower-case keyword → (MCM name or "" if the MCM has no such species, molar mass g/mol, density g/ml or None)
COMPOUNDS: Dict[str, Tuple[str, float, Optional[float]]] = {
    # terpenes / biogenics
    "alpha-pinene": ("APINENE", 136.23, 0.858), "a-pinene": ("APINENE", 136.23, 0.858), "α-pinene": ("APINENE", 136.23, 0.858),
    "apinene": ("APINENE", 136.23, 0.858),
    "beta-pinene": ("BPINENE", 136.23, 0.872), "b-pinene": ("BPINENE", 136.23, 0.872), "β-pinene": ("BPINENE", 136.23, 0.872),
    "bpinene": ("BPINENE", 136.23, 0.872),
    "d-limonene": ("LIMONENE", 136.23, 0.841), "limonene": ("LIMONENE", 136.23, 0.841),
    "beta-caryophyllene": ("BCARY", 204.36, 0.905), "caryophyllene": ("BCARY", 204.36, 0.905),
    "isoprene": ("C5H8", 68.12, 0.681), "pinonaldehyde": ("PINAL", 168.23, 0.98),
    # aromatics
    "1,3,5-trimethylbenzene": ("TM135B", 120.19, 0.865), "1,2,4-trimethylbenzene": ("TM124B", 120.19, 0.876),
    "1,2,3-trimethylbenzene": ("TM123B", 120.19, 0.894),
    "ethylbenzene": ("EBENZ", 106.17, 0.867), "o-xylene": ("OXYL", 106.17, 0.880), "m-xylene": ("MXYL", 106.17, 0.860),
    "p-xylene": ("PXYL", 106.17, 0.861), "toluene": ("TOLUENE", 92.14, 0.867), "benzene": ("BENZENE", 78.11, 0.879),
    "styrene": ("STYRENE", 104.15, 0.909), "phenol": ("PHENOL", 94.11, 1.07), "cresol": ("CRESOL", 108.14, 1.03),
    # oxygenates
    "methanol": ("CH3OH", 32.04, 0.792), "ethanol": ("C2H5OH", 46.07, 0.789),
    "isopropanol": ("IPROPOL", 60.10, 0.786), "i-propyl alcohol": ("IPROPOL", 60.10, 0.786), "2-propanol": ("IPROPOL", 60.10, 0.786),
    "isopropyl alcohol": ("IPROPOL", 60.10, 0.786), "n-propanol": ("NPROPOL", 60.10, 0.803), "1-propanol": ("NPROPOL", 60.10, 0.803),
    "n-butanol": ("NBUTOL", 74.12, 0.810), "1-butanol": ("NBUTOL", 74.12, 0.810), "butanol": ("NBUTOL", 74.12, 0.810),
    "tert-butanol": ("TBUTOL", 74.12, 0.781),
    "propylene glycol": ("PROPGLY", 76.09, 1.036), "1,2-propanediol": ("PROPGLY", 76.09, 1.036),
    "ethylene glycol": ("ETHGLY", 62.07, 1.113), "2-butoxyethanol": ("BUOX2ETOH", 118.17, 0.901),
    "butoxyethanol": ("BUOX2ETOH", 118.17, 0.901), "butyl glycol": ("BUOX2ETOH", 118.17, 0.901),
    "2-ethoxyethanol": ("", 90.12, 0.930), "ethoxyethanol": ("", 90.12, 0.930),
    "acetaldehyde": ("CH3CHO", 44.05, 0.784), "formaldehyde": ("HCHO", 30.03, None),
    "acetone": ("CH3COCH3", 58.08, 0.784), "mek": ("MEK", 72.11, 0.805), "methyl ethyl ketone": ("MEK", 72.11, 0.805),
    "2-butanone": ("MEK", 72.11, 0.805), "hydroxyacetone": ("ACETOL", 74.08, 1.08),
    "glyoxal": ("GLYOX", 58.04, None), "methylglyoxal": ("MGLYOX", 72.06, None),
    "ethyl acetate": ("ETHACET", 88.11, 0.902), "methyl acetate": ("METHACET", 74.08, 0.932),
    "acetic acid": ("CH3CO2H", 60.05, 1.049), "formic acid": ("HCOOH", 46.03, 1.22),
    "dimethyl ether": ("CH3OCH3", 46.07, None), "diethyl ether": ("DIETETHER", 74.12, 0.713),
    # alkanes
    "dodecane": ("NC12H26", 170.33, 0.750), "n-dodecane": ("NC12H26", 170.33, 0.750),
    "undecane": ("NC11H24", 156.31, 0.740), "decane": ("NC10H22", 142.28, 0.730), "n-decane": ("NC10H22", 142.28, 0.730),
    "nonane": ("NC9H20", 128.26, 0.718), "octane": ("NC8H18", 114.23, 0.703), "n-octane": ("NC8H18", 114.23, 0.703),
    "heptane": ("NC7H16", 100.20, 0.684), "n-heptane": ("NC7H16", 100.20, 0.684),
    "hexane": ("NC6H14", 86.18, 0.659), "n-hexane": ("NC6H14", 86.18, 0.659), "cyclohexane": ("CHEX", 84.16, 0.779),
    "pentane": ("NC5H12", 72.15, 0.626), "n-pentane": ("NC5H12", 72.15, 0.626), "isopentane": ("IC5H12", 72.15, 0.620),
    "n-butane": ("NC4H10", 58.12, None), "butane": ("NC4H10", 58.12, None), "isobutane": ("IC4H10", 58.12, None),
    "i-butane": ("IC4H10", 58.12, None), "propane": ("C3H8", 44.10, None), "ethane": ("C2H6", 30.07, None),
    "methane": ("CH4", 16.04, None),
    # alkenes
    "1-butene": ("BUT1ENE", 56.11, None), "trans-2-butene": ("TBUT2ENE", 56.11, None), "cis-2-butene": ("CBUT2ENE", 56.11, None),
    "isobutene": ("MEPROPENE", 56.11, None), "isobutylene": ("MEPROPENE", 56.11, None),
    "1,3-butadiene": ("C4H6", 54.09, None), "butadiene": ("C4H6", 54.09, None),
    "propene": ("C3H6", 42.08, None), "propylene": ("C3H6", 42.08, None),
    "ethene": ("C2H4", 28.05, None), "ethylene": ("C2H4", 28.05, None),
    # inorganic / tracers
    "carbon monoxide": ("CO", 28.01, None), "carbon dioxide": ("CO2", 44.01, None),
    "nitric oxide": ("NO", 30.01, None), "nitrogen dioxide": ("NO2", 46.01, None), "ozone": ("O3", 48.00, None),
    # things people inject that the MCM does not have (name known → MW/density known, MCM "" → surrogate needed)
    "d5 siloxane": ("", 370.77, 0.959), "decamethylcyclopentasiloxane": ("", 370.77, 0.959), "d5": ("", 370.77, 0.959),
    "triethanolamine": ("", 149.19, 1.124), "ethanolamine": ("", 61.08, 1.012), "monoethanolamine": ("", 61.08, 1.012),
    "carbitol": ("", 134.17, 0.999), "diethylene glycol monoethyl ether": ("", 134.17, 0.999),
    "tetrahydrofuran": ("", 72.11, 0.889), "thf": ("", 72.11, 0.889),
    "texanol": ("", 216.32, 0.95), "benzyl alcohol": ("", 108.14, 1.045),
    "tetradecane": ("", 198.39, 0.763), "n-tetradecane": ("", 198.39, 0.763), "hexadecane": ("", 226.45, 0.773),
    "n-hexadecane": ("", 226.45, 0.773), "octadecane": ("", 254.50, 0.777), "n-octadecane": ("", 254.50, 0.777),
}
_COMPOUND_KEYS = sorted(COMPOUNDS, key=len, reverse=True)
NON_MCM_OK = {"CO2"}                     # tracers that are fine as EASY species without an MCM entry

try:
    from mcm_species_db import MCM_DB, lookup_by_formula
except Exception:                                   # pragma: no cover
    MCM_DB, lookup_by_formula = {}, None


@dataclass
class ResolvedCompound:
    query: str
    mcm: str = ""                 # MCM species name ("" = none)
    mw: Optional[float] = None
    density: Optional[float] = None
    status: str = "unresolved"    # mcm | no_mcm | ambiguous | unresolved
    note: str = ""
    candidates: List[str] = field(default_factory=list)


def resolve_compound(text: str) -> ResolvedCompound:
    """Common name, MCM name or formula → ResolvedCompound.  Never raises."""
    q = (text or "").strip()
    low = q.lower()
    rc = ResolvedCompound(q)
    if not q:
        rc.note = "empty"; return rc
    # 1) name table (longest keyword first; whole-word match only — 'octadecane' must not hit 'decane')
    for key in _COMPOUND_KEYS:
        if re.search(rf"(?<![a-z0-9]){re.escape(key)}(?![a-z0-9])", low):
            mcm, mw, dens = COMPOUNDS[key]
            rc.mw, rc.density = mw, dens
            if mcm and (not MCM_DB or mcm in MCM_DB or mcm in NON_MCM_OK):
                rc.mcm, rc.status, rc.note = mcm, "mcm", f"'{key}' → MCM {mcm}"
            else:
                rc.status, rc.note = "no_mcm", f"'{key}' is not an MCM species — choose a surrogate"
            return rc
    # 2) typed MCM name
    up = q.upper()
    if up in MCM_DB:
        e = MCM_DB[up]
        rc.mcm, rc.mw, rc.status, rc.note = up, float(e.mass) if e.mass else None, "mcm", f"MCM species {up}"
        return rc
    # 3) formula
    if lookup_by_formula and re.fullmatch(r"[A-Z][A-Za-z0-9]*\d[A-Za-z0-9]*", q):
        hits = lookup_by_formula(q)
        if len(hits) == 1:
            rc.mcm, rc.mw, rc.status, rc.note = hits[0].name, float(hits[0].mass), "mcm", f"formula {q} → {hits[0].name}"
            return rc
        if hits:
            rc.candidates = [h.name for h in hits[:12]]
            rc.status, rc.note = "ambiguous", f"formula {q} matches {len(hits)} MCM species: " + ", ".join(rc.candidates)
            return rc
    rc.note = f"'{q}' not found (name table, MCM name, formula) — give the MCM name or a surrogate"
    return rc


# ─────────────────────────────────────────────────────────────────────────────
# Mixtures
# ─────────────────────────────────────────────────────────────────────────────

MIXTURE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "mixtures")
AMOUNT_UNITS = ("ul", "ml", "vol%", "mg", "g", "wt%", "mol%")


@dataclass
class MixtureComponent:
    name: str                       # as on the bottle / in the wiki
    amount: float                   # recipe amount
    unit: str = "ul"                # ul | ml | vol% | mg | g | wt% | mol%
    mcm: str = ""                   # MCM species (resolved or typed)
    surrogate: str = ""             # MCM species to use when `mcm` is empty ("" = not modelled)
    mw: Optional[float] = None
    density: Optional[float] = None
    note: str = ""

    def species(self) -> str:
        return self.mcm or self.surrogate


@dataclass
class Mixture:
    name: str
    description: str = ""
    aliases: List[str] = field(default_factory=list)     # words in a wiki line that mean this mixture ("VCP")
    components: List[MixtureComponent] = field(default_factory=list)
    source: str = ""

    def path(self) -> str:
        return os.path.join(MIXTURE_DIR, f"{safe_name(self.name)}.mixture.json")


def safe_name(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9_]+", "_", s.strip()).strip("_") or "mixture"


def save_mixture(m: Mixture, path: str = "") -> str:
    path = path or m.path()
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(asdict(m), fh, indent=1, ensure_ascii=False)
    return path


def load_mixture(path: str) -> Mixture:
    with open(path, encoding="utf-8") as fh:
        d = json.load(fh)
    comps = [MixtureComponent(**{k: v for k, v in c.items() if k in MixtureComponent.__dataclass_fields__})
             for c in d.get("components", [])]
    m = Mixture(**{k: v for k, v in d.items() if k in Mixture.__dataclass_fields__ and k != "components"})
    m.components = comps
    return m


def load_mixtures(directory: str = "") -> Dict[str, Mixture]:
    out: Dict[str, Mixture] = {}
    directory = directory or MIXTURE_DIR
    if os.path.isdir(directory):
        for f in sorted(os.listdir(directory)):
            if f.endswith(".mixture.json"):
                try:
                    m = load_mixture(os.path.join(directory, f)); out[m.name] = m
                except Exception:
                    pass
    return out


def find_mixture(text: str, mixtures: Dict[str, Mixture]) -> Optional[Mixture]:
    """Mixture whose name or alias occurs as a word in `text` (e.g. the wiki line)."""
    low = (text or "").lower()
    for m in mixtures.values():
        for word in [m.name] + list(m.aliases):
            if word and re.search(rf"(?<![a-z0-9]){re.escape(word.lower())}(?![a-z0-9])", low):
                return m
    return None


_RECIPE_LINE = re.compile(r"^\s*[\*\-•]?\s*(.+?)[\s,;:\t]+(\d+(?:[.,]\d+)?)\s*(µl|μl|ul|ml|vol%|mg|g|wt%|mol%|%)?\s*$", re.IGNORECASE)


def parse_mixture_text(text: str, default_unit: str = "ul") -> List[MixtureComponent]:
    """Recipe lines such as '* Limonene 5.61 μl', 'ethanol, 59.19, ul', 'toluene<TAB>6.08'."""
    comps: List[MixtureComponent] = []
    for line in text.splitlines():
        m = _RECIPE_LINE.match(line)
        if not m:
            continue
        name, val, unit = m.group(1).strip(), float(m.group(2).replace(",", ".")), (m.group(3) or default_unit).lower()
        unit = {"µl": "ul", "μl": "ul", "%": "vol%"}.get(unit, unit)
        rc = resolve_compound(name)
        comps.append(MixtureComponent(name, val, unit, rc.mcm, "", rc.mw, rc.density, rc.note))
    return comps


def resolve_mixture(m: Mixture) -> List[Tuple[str, str]]:
    """Fill mcm/mw/density of every component from the tables (typed values are kept).
    Returns alerts [(level, message)]."""
    alerts: List[Tuple[str, str]] = []
    for c in m.components:
        rc = resolve_compound(c.name)
        if not c.mcm and rc.mcm:
            c.mcm = rc.mcm
        c.mw = c.mw or rc.mw
        c.density = c.density or rc.density
        if c.mcm and MCM_DB and c.mcm not in MCM_DB:
            alerts.append(("error", f"{c.name}: '{c.mcm}' is not an MCM species"))
        if not c.mcm:
            if c.surrogate:
                if MCM_DB and c.surrogate not in MCM_DB:
                    alerts.append(("error", f"{c.name}: surrogate '{c.surrogate}' is not an MCM species"))
                else:
                    alerts.append(("warning", f"{c.name}: no MCM species — modelled as surrogate {c.surrogate}"))
            else:
                alerts.append(("warning", f"{c.name}: no MCM species and no surrogate — injected but NOT modelled ({rc.note})"))
        if c.unit in ("ul", "ml", "vol%") and not c.density:
            alerts.append(("error", f"{c.name}: density needed to convert {c.unit} to moles"))
        if not c.mw:
            alerts.append(("error", f"{c.name}: molar mass needed"))
    return alerts


def mixture_moles(m: Mixture, injected_ul: Optional[float]) -> Tuple[Dict[str, float], float, str]:
    """Moles of every component actually injected.  Recipe amounts are converted to
    (mass, volume) per component; if `injected_ul` is given the recipe is scaled so the
    total liquid volume equals it (mol% recipes need no scaling by volume — then
    injected_ul is treated as the total volume of the resulting liquid).
    Returns ({component name: mol}, recipe volume µl, explanation)."""
    mass: Dict[str, float] = {}
    vol: Dict[str, float] = {}
    for c in m.components:
        if not c.mw:
            continue
        u = c.unit
        if u == "ul":
            v = c.amount
        elif u == "ml":
            v = c.amount * 1000.0
        elif u == "vol%":
            v = c.amount
        elif u in ("mg", "g", "wt%"):
            mg = c.amount * (1000.0 if u == "g" else 1.0)
            v = mg * 1e-3 / c.density if c.density else 0.0
            mass[c.name] = mg * 1e-3
            vol[c.name] = v
            continue
        elif u == "mol%":
            mg = c.amount * c.mw               # relative moles → relative mass
            v = mg * 1e-3 / c.density if c.density else 0.0
            mass[c.name] = mg * 1e-3; vol[c.name] = v
            continue
        else:
            continue
        vol[c.name] = v
        mass[c.name] = v * 1e-3 * (c.density or 0.0)
    v_recipe = sum(vol.values())
    scale = (injected_ul / v_recipe) if (injected_ul and v_recipe > 0) else 1.0
    moles = {}
    for c in m.components:
        if c.name in mass and c.mw:
            moles[c.name] = mass[c.name] * scale / c.mw
    expl = (f"recipe volume {v_recipe:.3g} µl" + (f", injected {injected_ul:g} µl → scale {scale:.3g}" if injected_ul else ", used as is"))
    return moles, v_recipe, expl


def mixture_event(m: Mixture, injected_ul: Optional[float], start: str, duration_min: float,
                  T_K: float, P_hPa: float, ch: ChamberProfile, name: str = ""):
    """One InjectionEvent for the whole mixture: products = mole-fraction expression of the
    modelled components, total_ppb = modelled sum; plus precursor list and alerts."""
    from experiment_spec import InjectionEvent
    alerts = resolve_mixture(m)
    moles, v_recipe, expl = mixture_moles(m, injected_ul)
    n_air = moles_air(ch.volume_m3, T_K, P_hPa)
    ppb = {k: v / n_air * 1e9 for k, v in moles.items()}
    modelled: Dict[str, float] = {}
    unmodelled: List[str] = []
    for c in m.components:
        sp = c.species()
        if c.name not in ppb:
            unmodelled.append(c.name); continue
        if sp:
            modelled[sp] = modelled.get(sp, 0.0) + ppb[c.name]
        else:
            unmodelled.append(c.name)
    total = sum(modelled.values())
    total_all = sum(ppb.values())
    if total > 0:
        terms = [f"{v / total:.4g}*{sp}" for sp, v in sorted(modelled.items(), key=lambda kv: -kv[1])]
        products = "+".join(terms)
    else:
        products = ""
        alerts.append(("error", f"{m.name}: no component can be modelled"))
    ev_name = name or safe_name(m.name)
    lines = [f"mixture {m.name}: {expl}; air {n_air:.4g} mol at {T_K:g} K, {P_hPa:g} hPa",
             f"total injected {total_all:.3g} ppb, modelled {total:.3g} ppb in {len(modelled)} MCM species"]
    for c in m.components:
        if c.name in ppb:
            lines.append(f"  {c.name}: {ppb[c.name]:.3g} ppb → {c.species() or 'NOT MODELLED'}")
    if unmodelled:
        lines.append("  not modelled: " + ", ".join(unmodelled))
    ev = InjectionEvent(ev_name, products, f"{total:.4g}", start, duration_min, f"Q_{ev_name}", "",
                        comment=f"mixture {m.name}, {total_all:.3g} ppb injected, {total:.3g} ppb modelled (see .experiment.json)",
                        action="mixture",
                        delivery={"mode": "mixture", "mixture": m.name, "injected_ul": f"{injected_ul:g}" if injected_ul else "",
                                  "calc": "\n".join(lines)})
    return ev, list(modelled.keys()), alerts


# ─────────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    ch = chamber("SAPHIR")
    print(ppb_liquid(20, 68.12, 0.681, 298, 1013, ch))
    print(ppb_mfc(5000, 50, bottle_fraction(100, "ppm"), 298, 1013, ch))
    print(ppb_gas_syringe(3.2 * 2 + 0.32, 1.0, 298, 1013, ch))
    for q in ("limonene", "a-pinene", "D5 siloxane", "C10H16", "NC12H26", "carbitol", "foo"):
        r = resolve_compound(q); print(f"{q:12s} → {r.status:10s} {r.mcm:10s} {r.note}")
