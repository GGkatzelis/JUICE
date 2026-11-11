#!/usr/bin/env python
"""Build ``chemistry/photolysis/sigma_phi_298K.csv`` from an NCAR tuv-x checkout.

The table holds sigma(lambda, 298 K) * phi(lambda, 298 K, 1 atm) in cm^2 on a
1-nm grid (200-750 nm) for the MCM v3.3.1 photolysis reactions J1-J61.  JUICE
convolves it with a lamp spectrum F(lambda) to obtain J_i / J_NO2 ratios for
lamp-lit chambers (the "bottom-up" class-A method of Carter / SAPRC / F0AM).

Usage
-----
    python build_sigma_phi.py <path-to-tuv-x-checkout> [--out DIR] [--validate]

The module is importable without side effects; ``build_table()`` returns the
table as a dict of numpy arrays, ``write_csv()`` writes it, ``validate()`` runs
the spot checks and the clear-sky bracket test described in
``sigma_phi_298K.sources.md``.

Data source: NCAR tuv-x (https://github.com/NCAR/tuv-x), Apache-2.0.
Every parametrised cross section / quantum yield below is a transcription of
the corresponding ``src/cross_sections/*.F90`` or ``src/quantum_yields/*.F90``
routine, evaluated at T = 298 K, p = 1013.25 hPa, M = 2.46e19 cm^-3.  Where
tuv-x offers only a lumped channel (MVK, MACR, MEK, biacetyl, NOA, HPALD) the
MCM channel split / constant quantum yield is applied instead and documented
in the sources file.
"""
from __future__ import annotations

import csv
import math
import os
import sys
from typing import Dict, List, Optional, Tuple

import numpy as np

try:  # netCDF4 is only needed to *build*; keep import failure informative
    import netCDF4  # type: ignore
except ImportError as exc:  # pragma: no cover
    netCDF4 = None
    _NETCDF_IMPORT_ERROR = exc

# ---------------------------------------------------------------------------
# Evaluation conditions and grid
# ---------------------------------------------------------------------------
T_K = 298.0                      # temperature
P_HPA = 1013.25                  # pressure
M_AIR = 2.46e19                  # air number density, cm^-3 (1 atm, 298 K)
GRID = np.round(np.arange(200.0, 750.0 + 0.5, 1.0), 1)   # nm, inclusive

# Order of the columns in the CSV (numeric J order).
J_ORDER = ["J1", "J2", "J3", "J4", "J5", "J6", "J7", "J8",
           "J11", "J12", "J13", "J14", "J15", "J16", "J17", "J18", "J19", "J20",
           "J21", "J22", "J23", "J24",
           "J31", "J32", "J33", "J34", "J35",
           "J41",
           "J51", "J52", "J53", "J54", "J55", "J56",
           "J61"]

# MEK (J22) quantum yield: "tuvx" = IUPAC pressure-dependent yield as implemented in
# tuv-x (0.400 at 1 atm; reproduces the MCM l,m,n in the bracket test), or "mcm_page"
# = the constant 0.16 quoted on the MCM photolysis page.  See sources file.
PHI_MEK_SOURCE = "tuvx"
PHI_MEK_MCM_PAGE = 0.16

# J numbers for which tuv-x has no data (left empty on purpose).
EMPTY_J = {
    "J15": "n-C3H7CHO -> n-C3H7 + HCO: no n-butanal data in tuv-x",
    "J16": "n-C3H7CHO -> C2H4 + CH3CHO: no n-butanal data in tuv-x",
    "J17": "i-C3H7CHO -> products: no i-butanal data in tuv-x",
    "J61": "unidentified MCM J (7.537e-4, 0.499, 0.266); not on the MCM photolysis page",
}

# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------
H_PLANCK = 6.62607015e-34
C_LIGHT = 2.99792458e8


def _nc(path: str):
    if netCDF4 is None:  # pragma: no cover
        raise RuntimeError("netCDF4 is required to read tuv-x data files") from _NETCDF_IMPORT_ERROR
    return netCDF4.Dataset(path)


def read_nc(root: str, rel: str) -> Tuple[np.ndarray, np.ndarray, Optional[np.ndarray]]:
    """Return (wavelength[nm], params[nparams, nbins], temperatures or None)."""
    ds = _nc(os.path.join(root, rel))
    try:
        wl = np.array(ds.variables["wavelength"][:], dtype=float)
        var = "cross_section_parameters" if "cross_section_parameters" in ds.variables \
            else "quantum_yield_parameters"
        p = np.array(ds.variables[var][:], dtype=float)
        if p.ndim == 1:
            p = p[None, :]
        temps = None
        if "temperature" in ds.variables and len(ds.variables["temperature"][:]) > 0:
            temps = np.array(ds.variables["temperature"][:], dtype=float)
    finally:
        ds.close()
    return wl, p, temps


def tint(temps: np.ndarray, params: np.ndarray, T: float) -> np.ndarray:
    """tuv-x 'tint' rule: clamp T to the tabulated range, then interpolate
    linearly between the bracketing temperature rows (src/cross_sections/tint.F90,
    src/quantum_yields/tint.F90)."""
    order = np.argsort(temps)
    temps = temps[order]
    params = params[order]
    Tadj = min(max(T, temps[0]), temps[-1])
    if len(temps) == 1:
        return params[0]
    k = int(np.searchsorted(temps, Tadj, side="left"))
    k = max(1, min(k, len(temps) - 1))
    tstar = (Tadj - temps[k - 1]) / (temps[k] - temps[k - 1])
    return params[k - 1] + tstar * (params[k] - params[k - 1])


def on_grid(wl: np.ndarray, y: np.ndarray, lower: str = "zero", upper: str = "zero",
            lower_value: float = 0.0, upper_value: float = 0.0) -> np.ndarray:
    """Linear interpolation onto GRID; outside the data: 'zero', 'boundary'
    (hold the end value) or 'constant' (use *_value) - the tuv-x
    'lower/upper extrapolation' options."""
    order = np.argsort(wl)
    wl, y = wl[order], y[order]
    left = {"zero": 0.0, "boundary": y[0], "constant": lower_value}[lower]
    right = {"zero": 0.0, "boundary": y[-1], "constant": upper_value}[upper]
    return np.interp(GRID, wl, y, left=left, right=right)


def bin_average(wl: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Mean of a finely resolved source over each 1-nm bin [l-0.5, l+0.5)
    (used only for the 0.01-nm O3 data and the solar spectrum)."""
    order = np.argsort(wl)
    wl, y = wl[order], y[order]
    edges = np.concatenate([GRID - 0.5, [GRID[-1] + 0.5]])
    idx = np.digitize(wl, edges) - 1
    out = np.zeros_like(GRID)
    good = (idx >= 0) & (idx < len(GRID))
    counts = np.bincount(idx[good], minlength=len(GRID))
    sums = np.bincount(idx[good], weights=y[good], minlength=len(GRID))
    nz = counts > 0
    out[nz] = sums[nz] / counts[nz]
    # bins without any fine sample (none for the sources used) -> interpolate
    if (~nz).any():
        out[~nz] = np.interp(GRID[~nz], wl, y, left=0.0, right=0.0)
    return out


def coverage(wl: np.ndarray) -> str:
    return f"{wl.min():.1f}-{wl.max():.1f} nm"


# ---------------------------------------------------------------------------
# Parametrised cross sections (transcribed from tuv-x src/cross_sections/)
# ---------------------------------------------------------------------------
def xs_o3(root: str, meta: dict) -> np.ndarray:
    """O3 at 295 K: Malicet et al. 1995 / Brion et al. 1998, 195-830 nm,
    0.01-nm steps (data/cross_sections/O3_1.nc).  Bin-averaged to 1 nm."""
    wl, p, temps = read_nc(root, "data/cross_sections/O3_1.nc")
    meta["O3"] = {"file": "data/cross_sections/O3_1.nc (DATAE1NEW/O3/1995Malicet_1998Brion_O3_295.txt)",
                  "T": f"{temps[0]:.0f} K data", "coverage": coverage(wl), "method": "1-nm bin average"}
    return bin_average(wl, p[0])


def xs_h2o2(root: str, meta: dict) -> np.ndarray:
    """H2O2: JPL-94 / Kahan table (H2O2_1.nc) outside 260-350 nm; inside,
    the JPL temperature-dependent polynomial of src/cross_sections/h2o2-oh_oh.F90
    (Nicovich & Wine 1988 as recommended by JPL)."""
    wl, p, _ = read_nc(root, "data/cross_sections/H2O2_1.nc")
    sig = on_grid(wl, p[0])
    A = [6.4761E+04, -9.2170972E+02, 4.535649, -4.4589016E-03, -4.035101E-05,
         1.6878206E-07, -2.652014E-10, 1.5534675E-13]
    B = [6.8123E+03, -5.1351E+01, 1.1522E-01, -3.0493E-05, -1.0924E-07]
    lam = GRID
    sumA = np.polyval(A[::-1], lam)
    sumB = np.polyval(B[::-1], lam)
    t = min(max(T_K, 200.0), 400.0)
    chi = 1.0 / (1.0 + math.exp(-1265.0 / t))
    param = (chi * sumA + (1.0 - chi) * sumB) * 1.0e-21
    band = (lam >= 260.0) & (lam < 350.0)
    sig = np.where(band, param, sig)
    meta["H2O2"] = {"file": "data/cross_sections/H2O2_1.nc (DATAJ1NEW/ABS/H2O2_jpl94_Kahan.abs) + "
                            "src/cross_sections/h2o2-oh_oh.F90 polynomial for 260-350 nm",
                    "coverage": coverage(wl), "T": f"{T_K:.0f} K (chi = {chi:.4f})"}
    return sig


def xs_no2(root: str, meta: dict) -> np.ndarray:
    """NO2: NO2_1.nc has rows at 220 and 294 K ('NO2 tint' type, no2_tint.F90);
    298 K clamps to the 294 K row (Vandaele et al. 1998 / JPL)."""
    wl, p, temps = read_nc(root, "data/cross_sections/NO2_1.nc")
    sig = tint(temps, p, T_K)
    meta["NO2"] = {"file": "data/cross_sections/NO2_1.nc (rows 220/294 K; 'NO2 tint')",
                   "T": "294 K row (298 K clamped to table range)", "coverage": coverage(wl)}
    return on_grid(wl, sig)


def xs_hno3(root: str, meta: dict) -> np.ndarray:
    """HNO3: sigma(298) * exp(B (T-298)) (hno3-oh_no2.F90); Burkholder et al. 1993."""
    wl, p, _ = read_nc(root, "data/cross_sections/HNO3_1.nc")
    sig = p[0] * np.exp(p[1] * (T_K - 298.0))
    meta["HNO3"] = {"file": "data/cross_sections/HNO3_1.nc (DATAJ1NEW/ABS/HNO3_burk.abs; col1 sigma298, col2 B)",
                    "coverage": coverage(wl)}
    return on_grid(wl, sig)


def xs_ch2o(root: str, meta: dict) -> np.ndarray:
    """CH2O: sigma(298) + slope*(T-298) (ch2o.F90); JPL-11 (Meller & Moortgat 2000)."""
    wl, p, _ = read_nc(root, "data/cross_sections/CH2O_1.nc")
    sig = p[0] + p[1] * (T_K - 298.0)
    meta["CH2O"] = {"file": "data/cross_sections/CH2O_1.nc (DATAJ1NEW/CH2O/CH2O_jpl11.abs; col1 sigma298, col2 dsigma/dT)",
                    "coverage": coverage(wl)}
    return on_grid(wl, sig)


def xs_acetone(root: str, meta: dict) -> np.ndarray:
    """Acetone: sigma = c1 (1 + T(c2 + T(c3 + T c4))), T clamped to 235-298 K
    (acetone-ch3co_ch3.F90); JPL-11 / Gierczak et al. 1998."""
    wl, p, _ = read_nc(root, "data/cross_sections/CH3COCH3_1.nc")
    T = min(298.0, max(235.0, T_K))
    sig = p[0] * (1.0 + T * (p[1] + T * (p[2] + T * p[3])))
    meta["CH3COCH3"] = {"file": "data/cross_sections/CH3COCH3_1.nc (DATAJ1NEW/CH3COCH3/CH3COCH3_jpl11.abs; 4-coefficient T polynomial)",
                        "coverage": coverage(wl)}
    return on_grid(wl, sig)


def xs_rono2_exp(root: str, rel: str, key: str, meta: dict, origin: str) -> np.ndarray:
    """R-ONO2 with sigma(298) exp(B (T-298)) (rono2.F90 / ch3ono2-ch3o_no2.F90)."""
    wl, p, _ = read_nc(root, rel)
    sig = p[0] * np.exp(p[1] * (T_K - 298.0)) if p.shape[0] >= 2 else p[0]
    meta[key] = {"file": f"{rel} ({origin})", "coverage": coverage(wl)}
    return on_grid(wl, sig)


def xs_simple(root: str, rel: str, key: str, meta: dict, origin: str) -> np.ndarray:
    wl, p, temps = read_nc(root, rel)
    sig = p[0] if temps is None else tint(temps, p, T_K)
    meta[key] = {"file": f"{rel} ({origin})", "coverage": coverage(wl)}
    return on_grid(wl, sig)


def xs_t_butyl_nitrate(meta: dict) -> np.ndarray:
    """t-C4H9ONO2: exp(c + lam (b + a lam)) for 270-330 nm (t_butyl_nitrate.F90;
    fit to Roberts & Fajer 1989)."""
    a, b, c = -0.993E-3, 0.5307, -115.5
    lam = GRID
    sig = np.where((lam >= 270.0) & (lam <= 330.0), np.exp(c + lam * (b + a * lam)), 0.0)
    meta["tC4H9ONO2"] = {"file": "src/cross_sections/t_butyl_nitrate.F90 (analytic fit, Roberts & Fajer 1989)",
                         "coverage": "270.0-330.0 nm (fit range)"}
    return sig


def xs_nitroxy_acetone(meta: dict) -> np.ndarray:
    """NOA (CH3C(O)CH2ONO2): exp(c + lam (b + a lam)) for 284-335 nm
    (nitroxy_acetone.F90; fit to Roberts & Fajer 1989 / Barnes et al. 1993)."""
    a, b, c = -1.365E-3, 0.7834, -156.8
    lam = GRID
    sig = np.where((lam >= 284.0) & (lam <= 335.0), np.exp(c + lam * (b + a * lam)), 0.0)
    meta["NOA"] = {"file": "src/cross_sections/nitroxy_acetone.F90 (analytic fit, Roberts & Fajer 1989)",
                   "coverage": "284.0-335.0 nm (fit range)"}
    return sig


# ---------------------------------------------------------------------------
# Parametrised quantum yields (transcribed from tuv-x src/quantum_yields/)
# ---------------------------------------------------------------------------
def qy_o3_o1d() -> np.ndarray:
    """Matsumi et al. 2002 (JPL) O(1D) yield, src/quantum_yields/o3-o2_o1d.F90."""
    a = (0.8036, 8.9061, 0.1192)
    x = (304.225, 314.957, 310.737)
    om = (5.576, 6.601, 2.187)
    kt = 0.695 * T_K
    q1 = 1.0
    q2 = math.exp(-825.518 / kt)
    qfac1 = q1 / (q1 + q2)
    qfac2 = q2 / (q1 + q2)
    T300 = T_K / 300.0
    lam = GRID
    qy = np.zeros_like(lam)
    qy[lam <= 305.0] = 0.90
    mid = (lam > 305.0) & (lam <= 328.0)
    l = lam[mid]
    qy[mid] = (0.0765
               + a[0] * qfac1 * np.exp(-((x[0] - l) / om[0]) ** 4)
               + a[1] * T300 ** 2 * qfac2 * np.exp(-((x[1] - l) / om[1]) ** 2)
               + a[2] * T300 ** 1.5 * np.exp(-((x[2] - l) / om[2]) ** 2))
    qy[(lam > 328.0) & (lam <= 340.0)] = 0.08
    return qy


def qy_acetone_ch3co() -> np.ndarray:
    """Blitz et al. 2004 pressure/temperature-dependent CH3CO + CH3 yield,
    src/quantum_yields/acetone-ch3co_ch3.F90 (branch 'CH3CO', defaults:
    0.95 below 279 nm, 0 above 327 nm)."""
    Tadj = T_K / 295.0
    M = M_AIR
    lam = GRID
    qy = np.zeros_like(lam)
    for i, w in enumerate(lam):
        if w < 279.0:
            qy[i] = 0.95
        elif w > 327.0:
            qy[i] = 0.0
        else:
            a0 = 0.350 * Tadj ** (-1.28)
            b0 = 0.068 * Tadj ** (-2.65)
            dumexp = b0 * (w - 248.0)
            cA0 = 5.0e34 if dumexp > 80.0 else math.exp(dumexp) * a0 / (1.0 - a0)
            fco = 1.0 / (1.0 + cA0)
            wadj = 1.0e7 / w
            if w < 302.0:
                a1 = 1.600E-19 * Tadj ** (-2.38)
                b1 = 0.55E-3 * Tadj ** (-3.19)
                cA1 = a1 * math.exp(-b1 * (wadj - 33113.0))
                fac = (1.0 - fco) / (1.0 + cA1 * M)
            else:
                a2 = 1.62E-17 * Tadj ** (-10.03)
                b2 = 1.79E-3 * Tadj ** (-1.364)
                cA2 = a2 * math.exp(-b2 * (wadj - 30488.0))
                a3 = 26.29 * Tadj ** (-6.59)
                b3 = 5.72E-7 * Tadj ** (-2.93)
                c3 = 30006.0 * Tadj ** (-0.064)
                cA3 = a3 * math.exp(-b3 * (wadj - c3) ** 2)
                a4 = 1.67E-15 * Tadj ** (-7.25)
                b4 = 2.08E-3 * Tadj ** (-1.16)
                cA4 = a4 * math.exp(-b4 * (wadj - 30488.0))
                fac = (1.0 - fco) * (1.0 + cA3 + cA4 * M) / ((1.0 + cA3 + cA2 * M) * (1.0 + cA4 * M))
            qy[i] = fac
    return qy


def qy_mvk_total() -> np.ndarray:
    """Gierczak et al. 1997 pressure-dependent total MVK yield, src/quantum_yields/mvk.F90."""
    divisor = 5.5 + 9.2e-19 * M_AIR
    return np.minimum(np.exp(-0.055 * (GRID - 308.0)) / divisor, 1.0)


def qy_ch3cocho() -> np.ndarray:
    """Methylglyoxal CH3CO + HCO yield (JPL-11 / Chen et al. 2000 form),
    src/quantum_yields/ch3cocho.F90."""
    lam = GRID
    airfac = M_AIR * 760.0 / 2.456E19
    phi0 = np.clip(1.0 - (lam - 380.0) / 60.0, 0.0, 1.0)
    kq = 1.36e8 * np.exp(-8793.0 / lam)
    qy = np.where((phi0 > 0) & (lam >= 380.0) & (lam <= 440.0), phi0 / (phi0 + kq * airfac), phi0)
    qy[phi0 <= 0] = 0.0
    return qy


def qy_mek_tuvx() -> float:
    """tuv-x MEK yield 1/(0.96 + 2.22e-3 p[torr]) - reported for comparison only."""
    ptorr = 760.0 * M_AIR / 2.69e19
    return min(1.0, 1.0 / (0.96 + 2.22E-3 * ptorr))


def qy_ch2o_channels(root: str, meta: dict) -> Tuple[np.ndarray, np.ndarray]:
    """CH2O yields from data/quantum_yields/CH2O_1.nc (JPL-11): col1 = H + HCO,
    col2 = H2 + CO at 1 atm.  The H2+CO channel is pressure/temperature
    corrected for 330-360 nm as in src/quantum_yields/ch2o.F90.  Lower
    extrapolation 'boundary' as configured in photolysis_rate_constants.json."""
    wl, p, _ = read_nc(root, "data/quantum_yields/CH2O_1.nc")
    phi1 = on_grid(wl, p[0], lower="boundary")
    phi2 = on_grid(wl, p[1], lower="boundary")
    tfac = (300.0 - T_K) / 80.0
    lam = GRID
    band = (lam >= 330.0) & (lam < 360.0) & (phi2 > 0)
    tmp = 1.0 - phi1
    with np.errstate(divide="ignore", invalid="ignore"):
        wrk = (1.0 - (phi1 + phi2)) / (2.45e19 * phi2 * tmp)
        wrk = wrk * (1.0 + 0.05 * (lam - 329.0) * tfac)
        corr = 1.0 / (1.0 / tmp + wrk * M_AIR)
    phi2_p = np.where(band, corr, phi2)
    meta["CH2O_qy"] = {"file": "data/quantum_yields/CH2O_1.nc (DATAJ1/CH2O/CH2O_jpl11.yld; col1 H+HCO, col2 H2+CO at 1 atm)",
                       "coverage": coverage(wl), "extrapolation": "boundary below 250 nm, 0 above 361 nm"}
    return phi1, phi2_p


def qy_ch3cho(root: str, meta: dict) -> np.ndarray:
    """CH3CHO -> CH3 + HCO: data/quantum_yields/CH3CHO_1.nc (IUPAC): col1 =
    phi(CH4 + CO), col2 = phi(CH3 + HCO) at 1 atm; pressure correction as in
    src/quantum_yields/ch3cho-ch3_hco.F90.  No extrapolation configured in
    tuv-x -> 0 outside 260-330 nm."""
    wl, p, _ = read_nc(root, "data/quantum_yields/CH3CHO_1.nc")
    chnl1 = on_grid(wl, p[1])          # CH3 + HCO at 1 atm
    chnl2 = on_grid(wl, 1.0 - p[0])    # 1 - phi(CH4+CO)  (zero-pressure radical yield)
    with np.errstate(divide="ignore", invalid="ignore"):
        wrk = np.where(chnl1 > 0, chnl2 / chnl1 - 1.0, 0.0)
    qy = chnl1 * (1.0 + wrk) / (1.0 + wrk * M_AIR / 2.465e19)
    qy = np.clip(np.nan_to_num(qy), 0.0, 1.0)
    qy[(GRID < wl.min()) | (GRID > wl.max())] = 0.0
    meta["CH3CHO_qy"] = {"file": "data/quantum_yields/CH3CHO_1.nc (DATAJ1/CH3CHO/CH3CHO_iup.yld; col1 CH4+CO, col2 CH3+HCO)",
                         "coverage": coverage(wl), "extrapolation": "0 outside (tuv-x default)"}
    return qy


def qy_c2h5cho(root: str, meta: dict) -> np.ndarray:
    """C2H5CHO -> C2H5 + HCO: data/quantum_yields/C2H5CHO_1.nc (IUPAC, 1 atm),
    Stern-Volmer pressure scaling as in src/quantum_yields/c2h5cho.F90.
    No extrapolation configured -> 0 outside 294-340 nm."""
    wl, p, _ = read_nc(root, "data/quantum_yields/C2H5CHO_1.nc")
    phi = on_grid(wl, p[0])
    fac = M_AIR / 2.45e19
    with np.errstate(divide="ignore", invalid="ignore"):
        qy = np.where(phi > 1e-35, 1.0 / (1.0 + (1.0 / phi - 1.0) * fac), 0.0)
    qy = np.minimum(np.nan_to_num(qy), 1.0)
    meta["C2H5CHO_qy"] = {"file": "data/quantum_yields/C2H5CHO_1.nc (DATAJ1NEW/C2H5CHO/C2H5CHO_iup.yld)",
                          "coverage": coverage(wl), "extrapolation": "0 outside (tuv-x default)"}
    return qy


# ---------------------------------------------------------------------------
# Assemble the table
# ---------------------------------------------------------------------------
def build_table(root: str) -> Tuple[Dict[str, Optional[np.ndarray]], dict]:
    """Return ({J: sigma*phi array or None}, metadata)."""
    meta: dict = {}
    tab: Dict[str, Optional[np.ndarray]] = {j: None for j in J_ORDER}

    # --- O3 ---------------------------------------------------------------
    s_o3 = xs_o3(root, meta)
    phi_o1d = qy_o3_o1d()
    tab["J1"] = s_o3 * phi_o1d
    tab["J2"] = s_o3 * (1.0 - phi_o1d)
    meta["O3_qy"] = {"source": "src/quantum_yields/o3-o2_o1d.F90 (Matsumi et al. 2002); O(3P) = 1 - O(1D)",
                     "note": "phi(O1D) = 0.90 <=305 nm, Matsumi 305-328 nm, 0.08 328-340 nm, 0 above 340 nm"}

    # --- H2O2, NO2, NO3, HONO, HNO3 ---------------------------------------
    tab["J3"] = xs_h2o2(root, meta)                     # phi = 1
    s_no2 = xs_no2(root, meta)
    wl, p, temps = read_nc(root, "data/quantum_yields/NO2_1.nc")
    phi_no2 = on_grid(wl, tint(temps, p, T_K), lower="boundary")
    meta["NO2_qy"] = {"file": "data/quantum_yields/NO2_1.nc (DATAJ1/YLD/NO2_jpl11.yld; rows 298/248 K, 298 K row used)",
                      "coverage": coverage(wl), "extrapolation": "boundary (=1) below 300 nm, 0 above 422 nm (data end at 0.01)"}
    tab["J4"] = s_no2 * phi_no2

    s_no3 = xs_simple(root, "data/cross_sections/NO3_1.nc", "NO3", meta, "DATAJ1NEW/ABS/NO3_jpl11.abs")
    wl, p, temps = read_nc(root, "data/quantum_yields/NO3-NO+O2_1.nc")
    phi_no3_no = on_grid(wl, tint(temps, p, T_K))                     # 0 below 585 nm (data) and above 640 nm
    meta["NO3_qy_NO"] = {"file": "data/quantum_yields/NO3-NO+O2_1.nc (DATAJ1NEW/YLD/NO3_jpl2011.qy1; 298 K row)",
                         "coverage": coverage(wl), "extrapolation": "0 outside"}
    wl, p, temps = read_nc(root, "data/quantum_yields/NO3-NO2+O(3P)_1.nc")
    phi_no3_no2 = on_grid(wl, tint(temps, p, T_K), lower="constant", lower_value=1.0)
    meta["NO3_qy_NO2"] = {"file": "data/quantum_yields/NO3-NO2+O(3P)_1.nc (DATAJ1NEW/YLD/NO3_jpl2011.qy2; 298 K row)",
                          "coverage": coverage(wl), "extrapolation": "1.0 below 585 nm (configured), 0 above 640 nm"}
    tab["J5"] = s_no3 * phi_no3_no
    tab["J6"] = s_no3 * phi_no3_no2

    tab["J7"] = xs_simple(root, "data/cross_sections/HNO2_1.nc", "HONO", meta, "DATAJ1NEW/ABS/HONO_jpl11.abs")
    tab["J8"] = xs_hno3(root, meta)

    # --- carbonyls ---------------------------------------------------------
    s_ch2o = xs_ch2o(root, meta)
    phi1, phi2 = qy_ch2o_channels(root, meta)
    tab["J11"] = s_ch2o * phi1
    tab["J12"] = s_ch2o * phi2

    s_ch3cho = xs_simple(root, "data/cross_sections/CH3CHO_1.nc", "CH3CHO", meta, "DATAJ1NEW/CH3CHO/CH3CHO_jpl11.abs")
    tab["J13"] = s_ch3cho * qy_ch3cho(root, meta)

    s_c2h5cho = xs_simple(root, "data/cross_sections/C2H5CHO_1.nc", "C2H5CHO", meta, "DATAJ1NEW/C2H5CHO/C2H5CHO_iup.abs")
    tab["J14"] = s_c2h5cho * qy_c2h5cho(root, meta)

    # J15, J16, J17: no data in tuv-x -> None

    s_macr = xs_simple(root, "data/cross_sections/Methacrolein_1.nc", "MACR", meta, "DATAJ1NEW/ABS/Methacrolein_jpl11.abs")
    PHI_MACR = 1.95e-3            # MCM: per channel (Pinho et al. 2005 optimisation)
    tab["J18"] = s_macr * PHI_MACR
    tab["J19"] = s_macr * PHI_MACR
    tab["J20"] = s_macr * 1.0     # MCM definition of J20: sigma(MACR), phi = 1 (Peeters et al. 2009)
    meta["MACR_qy"] = {"source": "MCM: phi = 1.95e-3 for each of J18/J19; tuv-x uses a lumped 0.01",
                       "J20": "MCM prescribes sigma(MACR, Gierczak 1997) with phi = 1 for HPALD (proxy, not HPALD data)"}

    s_acet = xs_acetone(root, meta)
    tab["J21"] = s_acet * qy_acetone_ch3co()
    meta["CH3COCH3_qy"] = {"source": "src/quantum_yields/acetone-ch3co_ch3.F90 (Blitz et al. 2004), branch CH3CO",
                           "conditions": f"T = {T_K:.0f} K, M = {M_AIR:.3g} cm^-3; 0.95 below 279 nm, 0 above 327 nm"}

    s_mek = xs_simple(root, "data/cross_sections/CH3COCH2CH3_1.nc", "MEK", meta, "DATAJ1NEW/ABS/Martinez.abs")
    # The MCM page quotes phi = 0.16 (Pinho et al. 2005) but its l,m,n for J22
    # reproduce the clear-sky bracket test only with the IUPAC pressure-dependent
    # yield used by tuv-x (0.400 at 1 atm, src/quantum_yields/ch3coch2ch3-ch3co_ch2ch3.F90);
    # with 0.16 the computed J22/J4 is 0.56 (O3 only) / 0.33 (O3 + Rayleigh) of the MCM.
    phi_mek = qy_mek_tuvx() if PHI_MEK_SOURCE == "tuvx" else PHI_MEK_MCM_PAGE
    tab["J22"] = s_mek * phi_mek
    meta["MEK_qy"] = {"source": f"phi = {phi_mek:.3f} ({PHI_MEK_SOURCE}); tuv-x 1/(0.96 + 2.22e-3 p_torr) = {qy_mek_tuvx():.3f} "
                                f"at 1 atm, MCM page quotes {PHI_MEK_MCM_PAGE} (Pinho et al. 2005)"}

    s_mvk = xs_simple(root, "data/cross_sections/MVK_1.nc", "MVK", meta, "DATAJ1NEW/ABS/MVK_jpl11.abs")
    phi_mvk = qy_mvk_total()
    tab["J23"] = s_mvk * phi_mvk * 0.5
    tab["J24"] = s_mvk * phi_mvk * 0.5
    meta["MVK_qy"] = {"source": "src/quantum_yields/mvk.F90 (Gierczak et al. 1997) total yield at M = 2.46e19, "
                                "split 50/50 between J23 and J24 (MCM: equal l,m,n)"}

    # --- dicarbonyls -------------------------------------------------------
    s_gly = xs_simple(root, "data/cross_sections/CHOCHO_1.nc", "CHOCHO", meta, "DATAJ1NEW/CHOCHO/glyoxal_jpl11.abs")
    for j, rel in (("J31", "data/quantum_yields/CHOCHO-H2_CO_CO_1.nc"),
                   ("J32", "data/quantum_yields/CHOCHO-CH2O_CO_1.nc"),
                   ("J33", "data/quantum_yields/CHOCHO-HCO_HCO_1.nc")):
        wl, p, _ = read_nc(root, rel)
        tab[j] = s_gly * on_grid(wl, p[0], lower="boundary")
        meta[f"{j}_qy"] = {"file": f"{rel} (JPL-11 glyoxal channel yields)", "coverage": coverage(wl),
                           "extrapolation": "boundary below 225 nm, 0 above 420 nm"}

    s_mgly = xs_simple(root, "data/cross_sections/CH3COCHO_1.nc", "CH3COCHO", meta, "DATAJ1NEW/CH3COCHO/CH3COCHO_jpl11.abs")
    tab["J34"] = s_mgly * qy_ch3cocho()
    meta["CH3COCHO_qy"] = {"source": "src/quantum_yields/ch3cocho.F90 (JPL-11 / Chen et al. 2000) at M = 2.46e19"}

    s_biac = xs_simple(root, "data/cross_sections/CH3COCOCH3_1.nc", "BIACET", meta, "DATAJ1NEW/CH3COCOCH3/biacetyl_horowitz.abs")
    tab["J35"] = s_biac * 0.158   # MCM and tuv-x: Plum et al. 1983
    meta["BIACET_qy"] = {"source": "constant 0.158 (Plum et al. 1983) - same in MCM and tuv-x"}

    # --- peroxide and nitrates --------------------------------------------
    tab["J41"] = xs_simple(root, "data/cross_sections/CH3OOH_1.nc", "CH3OOH", meta, "DATAJ1NEW/CH3OOH/CH3OOH_jpl11.abs")
    tab["J51"] = xs_rono2_exp(root, "data/cross_sections/CH3ONO2_1.nc", "CH3ONO2", meta, "DATAJ1NEW/RONO2/CH3ONO2_jpl11.abs")
    tab["J52"] = xs_rono2_exp(root, "data/cross_sections/C2H5ONO2_1.nc", "C2H5ONO2", meta, "DATAJ1NEW/RONO2/C2H5ONO2_iup2006.abs")
    tab["J53"] = xs_simple(root, "data/cross_sections/nC3H7ONO2_1.nc", "nC3H7ONO2", meta, "DATAJ1NEW/RONO2/nC3H7ONO2_iup2006.abs")
    tab["J54"] = xs_rono2_exp(root, "data/cross_sections/CH3CHONO2CH3_1.nc", "iC3H7ONO2", meta, "DATAJ1NEW/RONO2/isopropyl_talukdar.abs")
    tab["J55"] = xs_t_butyl_nitrate(meta)
    tab["J56"] = xs_nitroxy_acetone(meta) * 0.9      # MCM phi = 0.9 (Mueller et al. 2014); tuv-x 1.0
    meta["NOA_qy"] = {"source": "MCM constant 0.9 (Mueller et al. 2014); tuv-x assumes 1.0"}

    # J61 unidentified -> None
    return tab, meta


def write_csv(tab: Dict[str, Optional[np.ndarray]], path: str) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["wavelength_nm"] + J_ORDER)
        for i, lam in enumerate(GRID):
            row = [f"{lam:.1f}"]
            for j in J_ORDER:
                col = tab[j]
                if col is None:
                    row.append("")
                else:
                    v = float(col[i])
                    row.append("0" if v == 0.0 else f"{v:.4e}")
            w.writerow(row)


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------
MCM_LMN = {  # l, m, n (n as on the MCM page, J = l cos^m(X) exp(-n sec X))
    "J1": (6.073e-5, 1.743, 0.474), "J2": (4.775e-4, 0.298, 0.080), "J3": (1.041e-5, 0.723, 0.279),
    "J4": (1.165e-2, 0.244, 0.267), "J5": (2.485e-2, 0.168, 0.108), "J6": (1.747e-1, 0.155, 0.125),
    "J7": (2.644e-3, 0.261, 0.288), "J8": (9.312e-7, 1.230, 0.307),
    "J11": (4.642e-5, 0.762, 0.353), "J12": (6.853e-5, 0.477, 0.323), "J13": (7.344e-6, 1.202, 0.417),
    "J14": (2.879e-5, 1.067, 0.358), "J15": (2.792e-5, 0.805, 0.338), "J16": (1.675e-5, 0.805, 0.338),
    "J17": (7.914e-5, 0.764, 0.364), "J18": (1.482e-6, 0.396, 0.298), "J19": (1.482e-6, 0.396, 0.298),
    "J20": (7.600e-4, 0.396, 0.298), "J21": (7.992e-7, 1.578, 0.271), "J22": (5.804e-6, 1.092, 0.377),
    "J23": (2.4246e-6, 0.395, 0.296), "J24": (2.424e-6, 0.395, 0.296),
    "J31": (6.845e-5, 0.130, 0.201), "J32": (1.032e-5, 0.130, 0.201), "J33": (3.802e-5, 0.644, 0.312),
    "J34": (1.537e-4, 0.170, 0.208), "J35": (3.326e-4, 0.148, 0.215), "J41": (7.649e-6, 0.682, 0.279),
    "J51": (1.588e-6, 1.154, 0.318), "J52": (1.907e-6, 1.244, 0.335), "J53": (2.485e-6, 1.196, 0.328),
    "J54": (4.095e-6, 1.111, 0.316), "J55": (1.135e-5, 0.974, 0.309), "J56": (4.365e-5, 1.089, 0.323),
    "J61": (7.537e-4, 0.499, 0.266),
}


def mcm_j(j: str, sza_deg: float) -> float:
    l, m, n = MCM_LMN[j]
    x = math.radians(sza_deg)
    return l * math.cos(x) ** m * math.exp(-n / math.cos(x))


def solar_photon_flux(root: str) -> np.ndarray:
    """Extraterrestrial spectrum (Chance & Kurucz 2010, data/profiles/solar/
    sao2010.solref.converted, W m-2 nm-1 at 0.01 nm) -> photons cm-2 s-1 nm-1
    averaged over the 1-nm bins of GRID."""
    path = os.path.join(root, "data/profiles/solar/sao2010.solref.converted")
    data = np.loadtxt(path, comments="#")
    wl, irr = data[:, 0], data[:, 1]
    ph = irr * 1.0e-4 * (wl * 1.0e-9) / (H_PLANCK * C_LIGHT)   # W m-2 nm-1 -> photon cm-2 s-1 nm-1
    return bin_average(wl, ph)


def validate(tab: Dict[str, Optional[np.ndarray]], root: str, sza_deg: float = 30.0,
             ozone_du: float = 300.0) -> Tuple[List[str], List[Tuple[str, float, float, float]]]:
    """Spot checks and the clear-sky bracket test.  Returns (lines, table)."""
    lines: List[str] = []
    g = GRID
    s_o3_tot = tab["J1"] + tab["J2"]              # sigma(O3) since phi1 + phi2 = 1

    def at(j, lam):
        return float(tab[j][int(round(lam - g[0]))])

    lines.append(f"sigma*phi(NO2) @254 nm = {at('J4', 254):.3e} cm^2 (expect ~1.05e-20, phi=1)")
    imax = int(np.argmax(tab["J4"]))
    lines.append(f"sigma*phi(NO2) max = {tab['J4'][imax]:.3e} cm^2 at {g[imax]:.0f} nm; "
                 f"@420 = {at('J4', 420):.2e}; @422 = {at('J4', 422):.2e}; @423 = {at('J4', 423):.2e} "
                 f"(0 above 422 nm because the JPL-11 yield data end there)")
    lines.append(f"sigma(O3) @305 nm = {float(s_o3_tot[305 - 200]):.3e}; phi(O1D) @305 = "
                 f"{tab['J1'][305 - 200] / s_o3_tot[305 - 200]:.3f}; sigma*phi(O1D) @305 = {at('J1', 305):.3e} cm^2")
    lines.append(f"sigma(H2O2) @254 nm = {at('J3', 254):.3e} cm^2 (expect ~7e-20)")
    lines.append(f"sigma(HONO) @354 nm = {at('J7', 354):.3e}, @368 nm = {at('J7', 368):.3e} cm^2 (expect ~2e-19)")
    lines.append(f"sigma*phi(HCHO->H+HCO) @305 nm = {at('J11', 305):.3e}, (->H2+CO) @330 nm = {at('J12', 330):.3e} cm^2")
    phi_ac = qy_acetone_ch3co()
    lines.append(f"acetone: phi(CH3CO+CH3) @290/300/310/320 nm = {phi_ac[90]:.3f}/{phi_ac[100]:.3f}/{phi_ac[110]:.3f}/{phi_ac[120]:.3f} "
                 f"(Blitz et al. 2004 at 298 K, 1 atm); sigma*phi @300 nm = {at('J21', 300):.3e} cm^2")
    phi_mg = qy_ch3cocho()
    lines.append(f"methylglyoxal: phi @380/400/420/440 nm = {phi_mg[180]:.3f}/{phi_mg[200]:.3f}/{phi_mg[220]:.3f}/{phi_mg[240]:.3f}")
    phi_mvk = qy_mvk_total()
    lines.append(f"MVK total phi @300/320/340/360 nm = {phi_mvk[100]:.3f}/{phi_mvk[120]:.3f}/{phi_mvk[140]:.3f}/{phi_mvk[160]:.3f}; "
                 f"MEK phi used = {qy_mek_tuvx() if PHI_MEK_SOURCE == 'tuvx' else PHI_MEK_MCM_PAGE:.3f} ({PHI_MEK_SOURCE})")

    # --- bracket test ----------------------------------------------------
    F0 = solar_photon_flux(root)
    sec = 1.0 / math.cos(math.radians(sza_deg))
    n_o3 = ozone_du * 2.6867e16                       # molec cm^-2 (300 DU = 8.06e18)
    tau_o3 = s_o3_tot * n_o3 * sec
    tau_r = 0.00877 * (g / 1000.0) ** (-4.05) * sec
    F_i = F0 * np.exp(-tau_o3)
    F_ii = F0 * np.exp(-tau_o3 - tau_r)

    def ratio(j, F):
        return float(np.sum(tab[j] * F) / np.sum(tab["J4"] * F))

    table = []
    for j in J_ORDER:
        if tab[j] is None or j == "J4":
            continue
        r_mcm = mcm_j(j, sza_deg) / mcm_j("J4", sza_deg)
        table.append((j, r_mcm, ratio(j, F_i), ratio(j, F_ii)))
    return lines, table


def format_table(table) -> str:
    out = ["| J | MCM J_i/J4 | case i (O3 only) | case ii (O3+Rayleigh) | (i)/MCM | (ii)/MCM |",
           "|---|---|---|---|---|---|"]
    for j, r_mcm, r_i, r_ii in table:
        out.append(f"| {j} | {r_mcm:.3e} | {r_i:.3e} | {r_ii:.3e} | {r_i / r_mcm:.2f} | {r_ii / r_mcm:.2f} |")
    return "\n".join(out)


# ---------------------------------------------------------------------------
def main(argv: Optional[List[str]] = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if not argv or argv[0] in ("-h", "--help"):
        print(__doc__)
        return 0
    root = argv[0]
    out_dir = os.path.dirname(os.path.abspath(__file__))
    if "--out" in argv:
        out_dir = argv[argv.index("--out") + 1]
    tab, meta = build_table(root)
    csv_path = os.path.join(out_dir, "sigma_phi_298K.csv")
    write_csv(tab, csv_path)
    print(f"wrote {csv_path}  ({len(GRID)} rows x {len(J_ORDER)} J columns; empty: {', '.join(EMPTY_J)})")
    if "--validate" in argv:
        lines, table = validate(tab, root)
        print("\n".join(lines))
        print()
        print(format_table(table))
        print()
        print("coverage / sources:")
        for k, v in meta.items():
            print(f"  {k}: {v}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
