"""Particles in the chamber: the condensation sink CS (s-1) of the gas-particle partitioning, and how the chamber's
size sets the wall terms.

Three sources of particles (ExperimentSpec.particles):
  "none"  no particles (the default): CS = 0, the condensable species exchange with the walls only
  "mono"  one size and number, e.g. 100 nm and 1e4 cm-3, held constant
  "file"  a measured size distribution (SMPS or similar): CS at every time of the table, as an ENZ column

CS = 2 pi D_g sum_i d_i N_i beta(Kn_i), with the Fuchs-Sutugin transition correction
  beta = (1 + Kn) / (1 + (4/(3 alpha) + 0.377) Kn + 4/(3 alpha) Kn^2),   Kn = 2 lambda / d,   lambda = 3 D_g / c_bar
for a generic low-volatility organic vapour: D_g = 5e-6 m2 s-1, M = 200 g mol-1, accommodation alpha = 1.
For the SAPHIR seed the partitioning include was written for (160 nm, 1800 cm-3) this gives 4.3e-3 s-1, against the
5.5e-3 s-1 written in it.

The size-distribution file. Teams write these differently, so the reader is lenient:
  * CSV, tab, semicolon or Excel; lines starting with # are skipped;
  * the first column is the time, in any form the measured-inputs file takes (YYYY-MM-DD HH:MM, HH:MM, hours of the
    day, seconds since 2000 or unix seconds);
  * every other column is one size bin, its header the bin's mid diameter in nm ("100", "100.5 nm", "Dp_100.5");
    a header in µm (all diameters below 5) is converted; columns whose header has no number are ignored;
  * values are dN/dlogDp in cm-3 (the usual SMPS output), or the number per bin in cm-3 (dn = "dN").
"""
from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np

D_GAS = 5e-6            # m2 s-1, gas-phase diffusion coefficient of a C10 oxygenated organic
M_GAS = 0.200           # kg mol-1
ALPHA = 1.0             # mass accommodation coefficient
R = 8.314462618

# the wall terms of the partitioning include (make_soa_condensation) were set for SAPHIR: these are its area and volume
SAPHIR_SURFACE_M2, SAPHIR_VOLUME_M3 = 320.0, 276.0


def mean_speed(T: float = 298.15, M: float = M_GAS) -> float:
    return math.sqrt(8.0 * R * T / (math.pi * M))


def fuchs_sutugin(d_m: np.ndarray, T: float = 298.15, D: float = D_GAS, M: float = M_GAS, alpha: float = ALPHA) -> np.ndarray:
    lam = 3.0 * D / mean_speed(T, M)
    kn = 2.0 * lam / np.asarray(d_m, float)
    a = 4.0 / (3.0 * alpha)
    return (1.0 + kn) / (1.0 + (a + 0.377) * kn + a * kn ** 2)


def condensation_sink(d_nm, n_cm3, T: float = 298.15, D: float = D_GAS, M: float = M_GAS, alpha: float = ALPHA) -> float:
    """CS (s-1) of particles with diameters d_nm and number concentrations n_cm3 (per bin, cm-3)."""
    d = np.atleast_1d(np.asarray(d_nm, float)) * 1e-9
    n = np.atleast_1d(np.asarray(n_cm3, float)) * 1e6
    ok = np.isfinite(d) & np.isfinite(n) & (d > 0)
    return float(2.0 * math.pi * D * np.sum(d[ok] * np.clip(n[ok], 0, None) * fuchs_sutugin(d[ok], T, D, M, alpha)))


def _dlogdp(d_nm: np.ndarray) -> np.ndarray:
    """Width in log10 of each bin, edges at the geometric midpoints (the outer bins mirror their neighbour)."""
    lg = np.log10(d_nm)
    if len(lg) == 1:
        return np.ones(1)
    edges = np.concatenate([[lg[0] - (lg[1] - lg[0]) / 2], (lg[1:] + lg[:-1]) / 2, [lg[-1] + (lg[-1] - lg[-2]) / 2]])
    return np.diff(edges)


@dataclass
class SizeDistribution:
    path: str
    times_s: np.ndarray                   # seconds since 2000-01-01 (ENZ time)
    d_nm: np.ndarray                      # bin mid diameters
    n_cm3: np.ndarray                     # number per bin, cm-3, shape (time, bin)
    notes: List[str] = field(default_factory=list)

    def cs(self, T: float = 298.15) -> np.ndarray:
        return np.array([condensation_sink(self.d_nm, row, T) for row in self.n_cm3])

    def summary(self) -> str:
        tot = np.nansum(self.n_cm3, axis=1)
        return (f"{len(self.times_s)} times, {len(self.d_nm)} bins {self.d_nm.min():.3g}–{self.d_nm.max():.3g} nm, "
                f"N {np.nanmin(tot):.3g}–{np.nanmax(tot):.3g} cm⁻³")


def _diameter(header) -> Optional[float]:
    m = re.search(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", str(header))
    return float(m.group(0)) if m else None


def read_size_distribution(path: str, date: str = "", dn: str = "dNdlogDp") -> SizeDistribution:
    from juice.core import measured_inputs as mi
    df = mi._load_table(path)
    tcol = df.columns[0]
    times = mi._parse_time(df[tcol], date)
    bins: List[Tuple[float, object]] = [(d, c) for c in df.columns[1:] if (d := _diameter(c)) is not None and d > 0]
    if not bins:
        raise ValueError("no size-bin columns: the headers after the time column must be diameters in nm")
    bins.sort()
    notes = []
    d = np.array([b[0] for b in bins])
    if d.max() < 5:
        d = d * 1000.0; notes.append("diameters read as µm → nm")
    vals = np.array([df[c].to_numpy(dtype=float) for _d, c in bins]).T
    if dn.lower().startswith("dndlog"):
        vals = vals * _dlogdp(d)[None, :]
    order = np.argsort(times)
    ok = np.isfinite(times[order])
    skipped = len(df.columns) - 1 - len(bins)
    if skipped:
        notes.append(f"{skipped} column(s) without a diameter in the header ignored")
    return SizeDistribution(path, times[order][ok], d, vals[order][ok], notes)


def cs_on_grid(path: str, date: str, grid_s: np.ndarray, T: float = 298.15, dn: str = "dNdlogDp") -> np.ndarray:
    """CS (s-1) on the ENZ grid; before the first and after the last time of the table the nearest value holds."""
    sd = read_size_distribution(path, date, dn)
    cs = sd.cs(T)
    good = np.isfinite(cs)
    return np.interp(grid_s, sd.times_s[good], cs[good])


def soa_defaults() -> Dict[str, float]:
    from juice.importers import make_soa_condensation as soa
    return {"KCOND": soa.KCOND, "KFWALL": soa.KFWALL, "KPWALL": soa.KPWALL, "CW_SCALE": soa.CW_SCALE,
            "COA0": 0.5, "CW_PRE": soa.CW0 * soa.SV_RATIO, "CW_N": soa.CW_N}


def read_soa_include(path: str) -> Dict:
    """Constants and per-species (name, C*, C_w) of a partitioning include as make_soa_condensation.emit_soa writes
    it; the module defaults for what the file does not say (or with no file)."""
    import os
    out: Dict = dict(soa_defaults(), species=[], path=path)
    if not path or not os.path.isfile(path):
        return out
    text = open(path, encoding="latin-1").read()
    for k in ("KCOND", "KFWALL", "KPWALL", "CW_SCALE"):
        m = re.search(rf"^\s*{k}\s*=\s*CONST\(\s*([0-9.eEdD+\-]+)\s*\)", text, re.M)
        if m:
            out[k] = float(m.group(1).replace("D", "E").replace("d", "e"))
    m = re.search(r"^\s*COA\s*=\s*CONST\(\s*([0-9.eE+\-]+)", text, re.M)
    if m:
        out["COA0"] = float(m.group(1))
    cw = {n: float(v) for n, v in re.findall(r"k\[WALL_(\w+)-->\1\]=CONST\(KFWALL\*[0-9.eE+\-]+/\(([0-9.eE+\-]+)\*CW_SCALE\)\)", text)}
    for name, cstar in re.findall(r"k\[SOA_(\w+)-->\1\]=CONST\(KCOND\*([0-9.eE+\-]+)/COA\)", text):
        out["species"].append((name, float(cstar), cw.get(name, float("nan"))))
    return out


def wall_scale(surface_m2: float, volume_m3: float) -> Optional[float]:
    """(A/V of the chamber) / (A/V of SAPHIR), or None when the chamber's wall area is not known."""
    if not surface_m2 or not volume_m3:
        return None
    return (surface_m2 / volume_m3) / (SAPHIR_SURFACE_M2 / SAPHIR_VOLUME_M3)
