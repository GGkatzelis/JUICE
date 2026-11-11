"""
chambers.py — the chamber as data: one profile per simulation chamber (no Qt).

JUICE grew up on SAPHIR, and SAPHIR's numbers used to be constants in the code: 276 m³ and a
0.6 m³/h leak in the DILUTE formula, the Jülich sun for the solar zenith angle, the Z: archive,
the Teflon-film wall sources.  They now live in a profile, `setups/chambers/<id>.chamber.json`, so another
chamber is a file, not a code change, and users can add their own.

    import chambers
    chambers.chamber("SAPHIR")          # ChamberProfile (the default when the id is unknown or empty)
    chambers.CHAMBERS                   # id → profile, every *.chamber.json in setups/chambers/ (+ the built-in SAPHIR)
    chambers.save_profile(profile)      # writes setups/chambers/<id>.chamber.json and reloads
    chambers.render_dilute(profile)     # the DILUTE expression for the EASY file

What a profile holds (see the JSON files for examples; every field has a default):

    id, name, institution, location, note, references, verified
    volume_m3, surface_m2, material
    leak_flow_m3h   flow entering the chamber besides the FLOW column (SAPHIR: 0.010*60 = 0.6 m³/h)
    dilute_expr     DILUTE verbatim; empty → CONST((FLOW+leak)/(volume*3600.)).  SAPHIR keeps the
                    IDL string byte for byte so the FACSIMILE scenarios stay identical.
    light           {"type": "sun" | "lamps" | "dark", "lat", "lon" (sun), "lamps": text,
                     "jno2_full_s1": jNO2 at full lamp power}.  Lamp photolysis itself is the next step;
                    until then every chamber runs with the sunlight parametrisation (PHTFKT).
    archive         {"kind": "saphir_nc", "root": path} for a chamber with a data archive JUICE can read; {} otherwise
    effects         default switches of experiment_spec.ChamberEffects (walldep, hono, hcho, bgreac,
                    dilution, dark_hono, light_source) — SAPHIR's wall sources are SAPHIR's; other chambers start with dilution only
    wall_loss_s1    {"O3": 2.1e-6, ...} first-order wall loss per species (stored now, used with the lamp step)
    mfc_std_T_K, mfc_std_P_hPa, syringe_T_K   standards of the injection calculator

Built-in SAPHIR: `_SAPHIR` below is the same profile as setups/chambers/SAPHIR.chamber.json, so JUICE works
(and the golden tests pass) when the folder is missing; `python -m juice.core.chambers --selftest` checks the two agree.
"""
from __future__ import annotations
from juice._root import ROOT as _JUICE_ROOT

import json
import os
import re
from dataclasses import dataclass, field, asdict
from typing import Dict, List, Optional

CHAMBER_DIR = os.path.join(_JUICE_ROOT, "setups", "chambers")
SUFFIX = ".chamber.json"
DEFAULT_CHAMBER = "SAPHIR"
SAPHIR_LAT, SAPHIR_LON = 50.9086, 6.3639          # Jülich Research Centre, degrees N / E
SAPHIR_DILUTE = "CONST((FLOW+0.010*60)/(276.*3600.))"   # the IDL formula, verbatim
EFFECT_KEYS = ("walldep", "hono", "hcho", "bgreac", "dilution", "dark_hono", "light_source")
LIGHT_TYPES = ("sun", "lamps", "dark")


@dataclass
class ChamberProfile:
    name: str
    volume_m3: float
    id: str = ""
    institution: str = ""
    location: str = ""
    note: str = ""
    references: List[str] = field(default_factory=list)
    verified: bool = False                      # numbers checked against the chamber's own publication / operators
    surface_m2: float = 0.0                     # 0 = not known
    material: str = ""
    leak_flow_m3h: float = 0.0
    dilute_expr: str = ""
    light: Dict = field(default_factory=lambda: {"type": "sun", "lat": SAPHIR_LAT, "lon": SAPHIR_LON})
    archive: Dict = field(default_factory=dict)
    effects: Dict[str, bool] = field(default_factory=lambda: {k: k == "dilution" for k in EFFECT_KEYS})
    wall_loss_s1: Dict[str, float] = field(default_factory=dict)
    mfc_std_T_K: float = 273.15                 # "standard" of the mass-flow controllers (sccm / slpm)
    mfc_std_P_hPa: float = 1013.25
    syringe_T_K: float = 293.15                 # gas syringes are filled at lab temperature

    # ── derived ──
    @property
    def light_type(self) -> str:
        t = str(self.light.get("type", "sun")).lower()
        return t if t in LIGHT_TYPES else "sun"

    @property
    def lat_lon(self):
        """(lat, lon) for the solar geometry; Jülich when the profile gives none (lamp chambers:
        the sz column is then only decorative)."""
        try:
            return float(self.light.get("lat", SAPHIR_LAT)), float(self.light.get("lon", SAPHIR_LON))
        except (TypeError, ValueError):
            return SAPHIR_LAT, SAPHIR_LON

    @property
    def archive_root(self) -> str:
        return str(self.archive.get("root", "")) if self.archive else ""

    @property
    def has_archive(self) -> bool:
        return bool(self.archive_root)

    def effect(self, key: str) -> bool:
        return bool(self.effects.get(key, key == "dilution"))

    def summary(self) -> str:
        light = {"sun": "sunlight", "lamps": "lamps" + (f" ({self.light.get('lamps')})" if self.light.get("lamps") else ""),
                 "dark": "no light source"}[self.light_type]
        s = f"{self.volume_m3:g} m³" + (f", {self.material}" if self.material else "") + f"; {light}"
        if self.light_type == "lamps" and self.light.get("jno2_full_s1"):
            s += f", jNO2 ≈ {float(self.light['jno2_full_s1']):.2g} 1/s"
        s += "; data archive readable" if self.has_archive else "; conditions from a file or typed"
        return s


def render_dilute(p: ChamberProfile) -> str:
    """DILUTE = (FLOW + leak) / (V · 3600): the fraction of the chamber replaced per second."""
    if p.dilute_expr:
        return p.dilute_expr
    return f"CONST((FLOW+{p.leak_flow_m3h:g})/({float(p.volume_m3)}*3600.))"


def safe_id(text: str) -> str:
    s = re.sub(r"[^A-Za-z0-9_.-]+", "_", (text or "").strip()).strip("_")
    return s or "chamber"


# The built-in SAPHIR — identical to setups/chambers/SAPHIR.chamber.json (checked by --selftest).
_SAPHIR = ChamberProfile(
    id="SAPHIR", name="SAPHIR", institution="Forschungszentrum Jülich, IEK-8 / ICE-3", location="Jülich, Germany",
    note="Jülich SAPHIR, 270 m³ nominal; 276 m³ is the value in the IDL DILUTE formula and is kept for consistency "
         "with the FACSIMILE scenarios.",
    references=["Rohrer et al. 2005, Atmos. Chem. Phys. 5, 2189", "Bohn et al. 2005, Atmos. Chem. Phys. 5, 493"],
    verified=True, volume_m3=276.0, surface_m2=320.0, material="double-wall FEP film",
    leak_flow_m3h=0.6, dilute_expr=SAPHIR_DILUTE,
    light={"type": "sun", "lat": SAPHIR_LAT, "lon": SAPHIR_LON, "lamps": "", "jno2_full_s1": 0.0},
    archive={"kind": "saphir_nc", "root": r"Z:\IEK8-SAPHIR\data"},
    effects={"walldep": True, "hono": True, "hcho": True, "bgreac": True, "dilution": True, "dark_hono": True, "light_source": False},
    wall_loss_s1={}, mfc_std_T_K=273.15, mfc_std_P_hPa=1013.25, syringe_T_K=293.15,
)


def profile_from_dict(d: Dict, fallback_id: str = "") -> ChamberProfile:
    d = dict(d)
    known = {k: v for k, v in d.items() if k in ChamberProfile.__dataclass_fields__}
    p = ChamberProfile(**known) if "name" in known and "volume_m3" in known else ChamberProfile(name=known.get("name", fallback_id), volume_m3=float(known.get("volume_m3", 0) or 0), **{k: v for k, v in known.items() if k not in ("name", "volume_m3")})
    p.id = p.id or fallback_id or safe_id(p.name)
    p.effects = {k: bool(p.effects.get(k, k == "dilution")) for k in EFFECT_KEYS}
    p.light = normalise_light(p.light)
    return p


LAMP_LIGHT_KEYS = ("spectrum", "spectrum_builtin", "solar_like", "effective_sza", "j_ratios")


def normalise_light(light: Optional[Dict]) -> Dict:
    """The light dict with its full key set for its type (so profiles compare equal whether they
    came from disk, the built-in or the editor): sun → lat, lon; lamps → the lamp J keys
    (see lamp_photolysis.py); every type → lamps text and jno2_full_s1."""
    L = dict(light or {})
    t = str(L.get("type", "sun")).lower()
    L["type"] = t if t in LIGHT_TYPES else "sun"
    if L["type"] == "sun":
        L.setdefault("lat", SAPHIR_LAT); L.setdefault("lon", SAPHIR_LON)
    L.setdefault("lamps", ""); L.setdefault("jno2_full_s1", 0.0)
    if L["type"] == "lamps":
        L.setdefault("spectrum", ""); L.setdefault("spectrum_builtin", ""); L.setdefault("solar_like", False)
        L.setdefault("effective_sza", 30.0); L.setdefault("j_ratios", {})
    return L


def profile_to_dict(p: ChamberProfile) -> Dict:
    return asdict(p)


def load_dir(directory: str = CHAMBER_DIR) -> Dict[str, ChamberProfile]:
    out: Dict[str, ChamberProfile] = {}
    if os.path.isdir(directory):
        for f in sorted(os.listdir(directory)):
            if not f.endswith(SUFFIX):
                continue
            stem = f[:-len(SUFFIX)]
            try:
                with open(os.path.join(directory, f), encoding="utf-8") as fh:
                    p = profile_from_dict(json.load(fh), stem)
                out[p.id] = p
            except Exception as e:                      # a broken user file must not take the GUI down
                print(f"[chambers] {f} not read: {e}")
    return out


CHAMBERS: Dict[str, ChamberProfile] = {}


def reload() -> Dict[str, ChamberProfile]:
    """Re-read setups/chambers/ into CHAMBERS (same dict object: importers keep their reference).
    SAPHIR comes first and is always present."""
    found = load_dir()
    CHAMBERS.clear()
    CHAMBERS[DEFAULT_CHAMBER] = found.pop(DEFAULT_CHAMBER, _SAPHIR)
    CHAMBERS.update(found)
    return CHAMBERS


def chamber(name: str = "") -> ChamberProfile:
    return CHAMBERS.get(name or DEFAULT_CHAMBER, CHAMBERS[DEFAULT_CHAMBER])


def save_profile(p: ChamberProfile, directory: str = CHAMBER_DIR) -> str:
    p.id = safe_id(p.id or p.name)
    os.makedirs(directory, exist_ok=True)
    path = os.path.join(directory, p.id + SUFFIX)
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(profile_to_dict(p), fh, indent=1, ensure_ascii=False)
        fh.write("\n")
    reload()
    return path


reload()


def _selftest() -> int:
    bad = 0
    disk = load_dir().get(DEFAULT_CHAMBER)
    if disk is None:
        print("FAIL setups/chambers/SAPHIR.chamber.json missing"); bad += 1
    elif profile_to_dict(disk) != profile_to_dict(_SAPHIR):
        a, b = profile_to_dict(disk), profile_to_dict(_SAPHIR)
        print("FAIL SAPHIR on disk differs from the built-in:", [k for k in a if a[k] != b.get(k)]); bad += 1
    else:
        print("OK   SAPHIR.chamber.json == built-in SAPHIR")
    if render_dilute(chamber("SAPHIR")) != SAPHIR_DILUTE:
        print("FAIL SAPHIR DILUTE changed:", render_dilute(chamber("SAPHIR"))); bad += 1
    else:
        print("OK   SAPHIR DILUTE verbatim:", SAPHIR_DILUTE)
    g = ChamberProfile(name="g", volume_m3=10, leak_flow_m3h=0.6)
    print("OK   rendered DILUTE for a 10 m³ chamber:", render_dilute(g))
    print(f"{len(CHAMBERS)} profiles: " + ", ".join(CHAMBERS))
    return 1 if bad else 0


if __name__ == "__main__":
    import sys
    sys.exit(_selftest() if "--selftest" in sys.argv else print(__doc__))
