"""
lamp_photolysis.py — photolysis frequencies for lamp-lit chambers (no Qt)

Sunlit chambers (SAPHIR, EUPHORE) run the mechanism's own clear-sky parametrisation scaled by the
measured jNO2 (PHTFKT in j_clear_sky.incl).  That shape is the solar spectrum at the real zenith
angle; it is wrong for lamps (zero at night, capped at clear sky, and a blacklight has no 300-320 nm
tail).  Every chamber study that modelled lamp experiments with the MCM or SAPRC did instead what
Carter's chamber software does (docs/lamp_photolysis_practice.md):

    J_i = jNO2 · ∫ σ_i(λ) φ_i(λ) F(λ) dλ / ∫ σ_NO2(λ) φ_NO2(λ) F(λ) dλ

with F the relative lamp spectrum and jNO2 from NO2 actinometry.  This module does that:

    lamp_j(profile)            → LampJ: ratio J_i / jNO2 for every MCM J<n>, with provenance
    easy_lines(lampj, ...)     → the EASY statements that replace the PHTFKT definitions
    spectrum_from_file(path)   → (λ, F) from a two-column text/CSV file
    clear_sky_ratios(sza)      → the MCM l,m,n ratios (what a solar-like source should give)

σ·φ at 298 K / 1 atm per MCM J number: chemistry/photolysis/sigma_phi_298K.csv (built from the NCAR tuv-x
data set by chemistry/photolysis/build_sigma_phi.py; sources in sigma_phi_298K.sources.md).

How a profile says what its lamps do — keys of setups/chambers/<id>.chamber.json "light" (type "lamps"):
    jno2_full_s1      jNO2 at full lamp power [1/s] (actinometry) — the absolute anchor
    spectrum          relative lamp spectrum: a file (path relative to setups/chambers/, or absolute) with
                      two columns wavelength_nm, intensity; header line '# units: energy|photons'
                      (default energy: W m-2 nm-1 or counts of a spectroradiometer → multiplied by λ)
    spectrum_builtin  "blacklight_350" | "blacklight_365" — synthetic Gaussian stand-ins, flagged
    solar_like        true for a filtered Xe/Ar arc: MCM clear-sky ratios at `effective_sza` (deg)
    j_ratios          {"J1": 1.0e-3, ...} measured ratios that override whatever the above gives
Priority: j_ratios > spectrum > spectrum_builtin > solar_like > (nothing: clear-sky ratios, loud note).
J numbers without σ·φ data are scaled from the covered J with the closest MCM (m, n) shape.
"""
from __future__ import annotations
from juice._root import ROOT as _JUICE_ROOT

import csv
import math
import os
import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np

HERE = _JUICE_ROOT
PHOTOLYSIS_DIR = os.path.join(HERE, "chemistry", "photolysis")
SIGMA_PHI_CSV = os.path.join(PHOTOLYSIS_DIR, "sigma_phi_298K.csv")
MCM_J_CSV = os.path.join(PHOTOLYSIS_DIR, "mcm_photolysis.csv")

# MCM v3.3.1 clear-sky parameters J = l·cos(χ)^m·exp(-n·sec χ) (mcm_converter.py / MCM web site)
MCM_LMN: Dict[int, Tuple[float, float, float]] = {
    1: (6.073e-5, 1.743, 0.474), 2: (4.775e-4, 0.298, 0.080), 3: (1.041e-5, 0.723, 0.279), 4: (1.165e-2, 0.244, 0.267),
    5: (2.485e-2, 0.168, 0.108), 6: (1.747e-1, 0.155, 0.125), 7: (2.644e-3, 0.261, 0.288), 8: (9.312e-7, 1.230, 0.307),
    11: (4.642e-5, 0.762, 0.353), 12: (6.853e-5, 0.477, 0.323), 13: (7.344e-6, 1.202, 0.417), 14: (2.879e-5, 1.067, 0.358),
    15: (2.792e-5, 0.805, 0.338), 16: (1.675e-5, 0.805, 0.338), 17: (7.914e-5, 0.764, 0.364), 18: (1.482e-6, 0.396, 0.298),
    19: (1.482e-6, 0.396, 0.298), 20: (7.600e-4, 0.396, 0.298), 21: (7.992e-7, 1.578, 0.271), 22: (5.804e-6, 1.092, 0.377),
    23: (2.4246e-6, 0.395, 0.296), 24: (2.424e-6, 0.395, 0.296), 31: (6.845e-5, 0.130, 0.201), 32: (1.032e-5, 0.130, 0.201),
    33: (3.802e-5, 0.644, 0.312), 34: (1.537e-4, 0.170, 0.208), 35: (3.326e-4, 0.148, 0.215), 41: (7.649e-6, 0.682, 0.279),
    51: (1.588e-6, 1.154, 0.318), 52: (1.907e-6, 1.244, 0.335), 53: (2.485e-6, 1.196, 0.328), 54: (4.095e-6, 1.111, 0.316),
    55: (1.135e-5, 0.974, 0.309), 56: (4.365e-5, 1.089, 0.323), 61: (7.537e-4, 0.499, 0.266),
}
J_NO2 = 4
# J numbers without σ·φ data in chemistry/photolysis/sigma_phi_298K.csv and the chemically closest covered one
# (n-/i-butanal → propanal); anything else falls back to the closest MCM (m, n) shape
MISSING_ANALOGUE: Dict[int, int] = {15: 14, 16: 14, 17: 14}

# the inorganic photolysis reactions as anorg.mcm331.phot.rate / the converter write them
# (the EASY k[...] override needs the reaction text; the organic ones use the J<n> variable)
ANORG_REACTIONS: Dict[int, str] = {
    1: "O3+hv-->O1D", 2: "O3+hv-->O", 3: "H2O2+hv-->OH+OH", 4: "NO2+hv-->NO+O",
    5: "NO3+hv-->NO", 6: "NO3+hv-->NO2+O", 7: "HONO+hv-->OH+NO", 8: "HNO3+hv-->OH+NO2",
}
# measured-J variable of the ENZ/EASY (experiment_spec.J_INPUTS) → J number it replaces
MEASURED_J: Dict[str, int] = {"jno2": 4, "jo1d": 1, "jhono": 7, "jh2o2": 3, "jhcho_r": 11, "jhcho_m": 12}

BUILTIN_SPECTRA = {                 # name → (centre nm, FWHM nm, note); Gaussian stand-ins until a measured spectrum exists
    "blacklight_350": (352.0, 40.0, "synthetic Gaussian for a UV-A 350 nm blacklight (Sylvania F40/350BL class); "
                                   "the 300-320 nm tail that sets j(O1D) and j(H2O2) is NOT known — measure the lamps"),
    "blacklight_365": (368.0, 25.0, "synthetic Gaussian for a 365 nm BLB / TL-K 40W/05 lamp; no emission below 340 nm"),
}


@dataclass
class LampJ:
    ratios: Dict[int, float] = field(default_factory=dict)   # J number → J / jNO2
    source: str = ""                                          # one line of provenance
    notes: List[str] = field(default_factory=list)
    method: str = ""                                          # "measured ratios" | "spectrum" | "builtin spectrum" | "solar-like" | "clear-sky fallback"
    scaled_from: Dict[int, int] = field(default_factory=dict) # J without σφ data → J it was scaled from

    def ratio(self, j: int) -> float:
        return float(self.ratios.get(j, 0.0))


# ── σ·φ table ────────────────────────────────────────────────────────────────
_SIGMA_PHI: Optional[Tuple[np.ndarray, Dict[int, np.ndarray]]] = None


def sigma_phi(path: str = SIGMA_PHI_CSV) -> Tuple[np.ndarray, Dict[int, np.ndarray]]:
    """(wavelength nm, {J number: σ·φ in cm²}) from chemistry/photolysis/sigma_phi_298K.csv; J columns that are
    empty in the file are left out."""
    global _SIGMA_PHI
    if _SIGMA_PHI is not None and path == SIGMA_PHI_CSV:
        return _SIGMA_PHI
    if not os.path.isfile(path):
        raise FileNotFoundError(f"{path} — run chemistry/photolysis/build_sigma_phi.py <tuv-x checkout> to build it")
    with open(path, encoding="utf-8") as fh:
        rows = list(csv.reader(fh))
    head = [h.strip() for h in rows[0]]
    data = np.array([[float(v) if v.strip() else np.nan for v in r] for r in rows[1:] if r and r[0].strip()])
    lam = data[:, 0]
    out: Dict[int, np.ndarray] = {}
    for i, h in enumerate(head[1:], start=1):
        m = re.fullmatch(r"J(\d+)", h)
        if m and not np.all(np.isnan(data[:, i])):
            out[int(m.group(1))] = np.nan_to_num(data[:, i], nan=0.0)
    if path == SIGMA_PHI_CSV:
        _SIGMA_PHI = (lam, out)
    return lam, out


def sigma_phi_available() -> bool:
    return os.path.isfile(SIGMA_PHI_CSV)


# ── spectra ──────────────────────────────────────────────────────────────────
def spectrum_from_file(path: str, units: str = "auto") -> Tuple[np.ndarray, np.ndarray, str]:
    """Two columns wavelength [nm], intensity; comment lines start with # or ;.  Returns
    (λ, F as relative photon flux, note).  units: 'photons' = the file is a photon flux already;
    'energy' = W m⁻² nm⁻¹ or spectroradiometer counts → F_photons ∝ F_energy·λ; 'auto' reads a
    header line '# units: photons|energy' and defaults to energy."""
    lam, val, hdr_units = [], [], ""
    with open(path, encoding="utf-8-sig", errors="replace") as fh:
        for line in fh:
            s = line.strip()
            if not s:
                continue
            if s[0] in "#;":
                m = re.search(r"units?\s*[:=]\s*(photons?|energy|w|counts)", s, re.I)
                if m:
                    hdr_units = "photons" if m.group(1).lower().startswith("photon") else "energy"
                continue
            parts = [p for p in re.split(r"[,\s;\t]+", s) if p]
            try:
                x, y = float(parts[0]), float(parts[1])
            except (ValueError, IndexError):
                continue                                  # header row with names
            lam.append(x); val.append(y)
    if len(lam) < 5:
        raise ValueError(f"{path}: fewer than 5 numeric rows (expected wavelength_nm, intensity)")
    lam_a, val_a = np.asarray(lam, float), np.asarray(val, float)
    order = np.argsort(lam_a); lam_a, val_a = lam_a[order], np.clip(val_a[order], 0, None)
    u = units if units != "auto" else (hdr_units or "energy")
    if u == "energy":
        val_a = val_a * lam_a
    note = f"spectrum {os.path.basename(path)}: {lam_a.min():g}-{lam_a.max():g} nm, {len(lam_a)} points, read as {u} units" + \
           ("" if hdr_units or units != "auto" else " (no '# units:' header - energy assumed, i.e. multiplied by wavelength)")
    return lam_a, val_a, note


def builtin_spectrum(name: str, sigma_short: float = 0.0) -> Tuple[np.ndarray, np.ndarray, str]:
    """Gaussian stand-in; `sigma_short` narrows the short-wavelength side (the phosphor tail that
    sets j(O1D) and j(H2O2) — unknown without a measurement)."""
    centre, fwhm, note = BUILTIN_SPECTRA[name]
    lam = np.arange(200.0, 750.5, 1.0)
    sig = fwhm / 2.3548
    s_lo = sigma_short or sig
    F = np.where(lam < centre, np.exp(-0.5 * ((lam - centre) / s_lo) ** 2), np.exp(-0.5 * ((lam - centre) / sig) ** 2))
    return lam, F, f"built-in {name}: {note}"


def tail_sensitivity(name: str) -> str:
    """How much j(O1D), j(H2O2) and j(HCHO_r) of a stand-in move when its short-wavelength side is
    30 % narrower — the honest error bar of a synthetic spectrum."""
    sig = BUILTIN_SPECTRA[name][1] / 2.3548
    lam, F, _ = builtin_spectrum(name); wide, _n = convolve_ratios(lam, F)
    lam, F, _ = builtin_spectrum(name, sigma_short=0.7 * sig); narrow, _n = convolve_ratios(lam, F)
    return ("stand-in sensitivity to the unknown <320 nm tail (symmetric vs 30 % narrower short side): "
            + ", ".join(f"J{j} {wide[j]:.1e} vs {narrow[j]:.1e}" for j in (1, 3, 11) if j in wide)
            + " - measure the lamp spectrum to pin these down")


# ── the convolution ──────────────────────────────────────────────────────────
def convolve_ratios(lam_f: np.ndarray, F: np.ndarray) -> Tuple[Dict[int, float], List[str]]:
    """J_i / J_NO2 for every J with σ·φ data, on the σ·φ wavelength grid (F interpolated, 0 outside)."""
    lam, table = sigma_phi()
    Fi = np.interp(lam, lam_f, F, left=0.0, right=0.0)
    notes: List[str] = []
    if lam_f.min() > 300 or lam_f.max() < 400:
        notes.append(f"spectrum covers only {lam_f.min():g}-{lam_f.max():g} nm; photolysis outside is taken as zero")
    denom = float(np.trapezoid(table[J_NO2] * Fi, lam))
    if denom <= 0:
        raise ValueError("the spectrum has no overlap with the NO2 absorption (300-420 nm): cannot normalise to jNO2")
    out = {j: float(np.trapezoid(sp * Fi, lam)) / denom for j, sp in table.items()}
    return out, notes


def clear_sky_ratios(sza_deg: float = 30.0) -> Dict[int, float]:
    """MCM l,m,n ratios J_i / J4 at one zenith angle (the MCM's own clear-sky spectrum)."""
    c = max(math.cos(math.radians(sza_deg)), 0.01)
    def J(n):
        l, m, k = MCM_LMN[n]
        return l * c ** m * math.exp(-k / c)
    j4 = J(J_NO2)
    return {n: J(n) / j4 for n in MCM_LMN}


def _fill_missing(ratios: Dict[int, float], sza_deg: float = 30.0) -> Dict[int, int]:
    """J numbers absent from `ratios` (no σ·φ data): scale from the covered J with the closest MCM
    (m, n) shape — identical (m, n) means the MCM author used the same spectrum assumption — by the
    ratio of their clear-sky values.  Returns {filled J: reference J}."""
    cs = clear_sky_ratios(sza_deg)
    covered = [j for j in ratios if j in MCM_LMN]
    filled: Dict[int, int] = {}
    for j, (l, m, n) in MCM_LMN.items():
        if j in ratios or not covered:
            continue
        ref = MISSING_ANALOGUE.get(j)
        if ref not in covered:
            ref = min(covered, key=lambda k: (MCM_LMN[k][1] - m) ** 2 + (MCM_LMN[k][2] - n) ** 2)
        ratios[j] = ratios[ref] * cs[j] / cs[ref]
        filled[j] = ref
    return filled


# ── profile → LampJ ──────────────────────────────────────────────────────────
def _resolve_spectrum_path(name: str, chamber_dir: str) -> str:
    if os.path.isabs(name) and os.path.isfile(name):
        return name
    for base in (chamber_dir, os.path.join(chamber_dir, "spectra"), PHOTOLYSIS_DIR, os.getcwd()):
        p = os.path.join(base, name)
        if os.path.isfile(p):
            return p
    raise FileNotFoundError(f"lamp spectrum file not found: {name} (looked in {chamber_dir})")


def lamp_j(profile, chamber_dir: str = "") -> LampJ:
    """Resolve the J_i / jNO2 table of a lamp profile (see the module doc for the keys)."""
    from juice.core import chambers
    chamber_dir = chamber_dir or chambers.CHAMBER_DIR
    light = dict(getattr(profile, "light", {}) or {})
    lj = LampJ()
    overrides: Dict[int, float] = {}
    for k, v in (light.get("j_ratios") or {}).items():
        m = re.fullmatch(r"J?(\d+)", str(k).strip(), re.I)
        if m and v is not None:
            overrides[int(m.group(1))] = float(v)
    try:
        if light.get("spectrum"):
            path = _resolve_spectrum_path(str(light["spectrum"]), chamber_dir)
            lam, F, note = spectrum_from_file(path, str(light.get("spectrum_units", "auto")))
            lj.ratios, notes = convolve_ratios(lam, F)
            lj.method, lj.source = "spectrum", f"lamp spectrum x sigma*phi (tuv-x/JPL/IUPAC at 298 K): {os.path.basename(path)}"
            lj.notes += [note] + notes
        elif light.get("spectrum_builtin"):
            name = str(light["spectrum_builtin"])
            lam, F, note = builtin_spectrum(name)
            lj.ratios, notes = convolve_ratios(lam, F)
            lj.method, lj.source = "builtin spectrum", f"SYNTHETIC {name} spectrum x sigma*phi - a stand-in, not a measurement"
            lj.notes += [note] + notes + [tail_sensitivity(name)]
    except (FileNotFoundError, ValueError) as e:
        lj.notes.append(f"lamp spectrum not usable: {e}")
    if not lj.ratios:
        sza = float(light.get("effective_sza", 30.0) or 30.0)
        lj.ratios = clear_sky_ratios(sza)
        if light.get("solar_like"):
            lj.method, lj.source = "solar-like", f"solar-like source: MCM clear-sky ratios at {sza:g} deg zenith angle (filtered Xe/Ar arc)"
        else:
            lj.method, lj.source = "clear-sky fallback", f"NO lamp spectrum in the profile: MCM clear-sky ratios at {sza:g} deg used - wrong for blacklights (j(O1D), j(H2O2) too high)"
            lj.notes.append("give the profile a lamp spectrum (light.spectrum), a built-in stand-in (light.spectrum_builtin) or measured ratios (light.j_ratios)")
    lj.scaled_from = _fill_missing(lj.ratios)
    if overrides:
        lj.ratios.update(overrides)
        for j in overrides:
            lj.scaled_from.pop(j, None)
        lj.method = "measured ratios" if lj.method in ("clear-sky fallback",) and set(overrides) >= {1, 3, 7, 11, 12} else lj.method
        lj.source += f"; measured ratios for J{', J'.join(str(j) for j in sorted(overrides))}"
    lj.ratios[J_NO2] = 1.0
    return lj


# ── EASY output ──────────────────────────────────────────────────────────────
def easy_lines(lj: LampJ, jno2_var: str = "jno2", measured: Optional[Dict[str, int]] = None,
               profile_id: str = "", compare_sza: float = 30.0) -> List[str]:
    """The statements that make the mechanism use the lamp J values:  J<n> = CONST(jno2*ratio)
    for every J number (the organic reactions use these variables) and k[...] = CONST(...) for the
    inorganic reactions (their rates are bound to the PHTFKT definitions).  J values the ENZ carries
    as measurements (`measured`: variable → J number) are tied to that variable instead."""
    measured = measured or {}
    by_j = {j: var for var, j in measured.items()}
    cs = clear_sky_ratios(compare_sza)
    L = [f"; LAMP PHOTOLYSIS ({profile_id or 'lamp chamber'}): J<n> = {jno2_var} x ratio, {lj.source}",
         f";   ratios j(X)/j(NO2) of this light source (clear sky at {compare_sza:g} deg in brackets): "
         + ", ".join(f"J{j} {lj.ratio(j):.3g} ({cs[j]:.3g})" for j in (1, 3, 7, 11, 12, 13, 21, 41) if j in lj.ratios)]
    for n in lj.notes:
        L.append(f";   note: {n}")
    if lj.scaled_from:
        L.append(";   no sigma*phi data, scaled from the J with the closest MCM shape: "
                 + ", ".join(f"J{j}<-J{r}" for j, r in sorted(lj.scaled_from.items())))
    for j in sorted(lj.ratios):
        if j == J_NO2:
            expr = jno2_var
        elif j in by_j:
            expr = by_j[j]
        else:
            expr = f"{jno2_var}*{lj.ratios[j]:.4e}"
        L.append(f"J{j} = CONST({expr})")
        if j in ANORG_REACTIONS and j not in by_j and j != J_NO2:
            L.append(f"k[{ANORG_REACTIONS[j]}] = CONST({expr})")
    return L


def summary(lj: LampJ) -> str:
    cs = clear_sky_ratios(30.0)
    parts = [f"j(O1D)/j(NO2) = {lj.ratio(1):.2g} (clear sky {cs[1]:.2g})", f"j(HONO)/j(NO2) = {lj.ratio(7):.2g} ({cs[7]:.2g})",
             f"j(H2O2)/j(NO2) = {lj.ratio(3):.2g} ({cs[3]:.2g})"]
    return f"{lj.method}: " + ", ".join(parts)


if __name__ == "__main__":
    import sys
    if "--selftest" in sys.argv:
        from juice.core import chambers
        lam, tab = sigma_phi()
        print(f"sigma*phi table: {len(lam)} wavelengths {lam.min():g}-{lam.max():g} nm, J numbers: {sorted(tab)}")
        for name in BUILTIN_SPECTRA:
            l, F, _ = builtin_spectrum(name); r, _n = convolve_ratios(l, F)
            print(f"{name:16s} " + "  ".join(f"J{j}/J4={r[j]:.2e}" for j in (1, 3, 7, 11, 12, 13, 21, 41) if j in r))
        cs = clear_sky_ratios(30)
        print(f"{'MCM clear sky 30':16s} " + "  ".join(f"J{j}/J4={cs[j]:.2e}" for j in (1, 3, 7, 11, 12, 13, 21, 41)))
        for pid, p in chambers.CHAMBERS.items():
            if p.light_type == "lamps":
                print(f"{pid:14s} {summary(lamp_j(p))}")
    else:
        print(__doc__)
