"""GECKO-A photolysis label -> EASY/JUICE rate expression, scaled to SAPHIR's measured J-values.

GECKO-A's box model (boxmodel4gecko) has no closed-form J(SZA); it only tabulates J at 11
SZA points per label (`INPUT/solarlight.phot`, copy in gecko_data/) and cubic-spline-interpolates
at runtime. For the SAPHIR-measured channels we use the real measurement directly. For every
other GECKO-A label we fit its native table to the functional form MCM uses (`j_clear_sky.incl`:
J = l*cos(SZA)^m*EXP(-n*sec(SZA))), then divide two such fits -- l, m, n of the two channels
combine algebraically into one SZA-dependent ratio, exactly the way `Jx(SZA)/JNO2(SZA)*JNO2(observed)`
was always done, just with GECKO-A's own l,m,n instead of MCM's. The anchor (jno2 vs jo1d) is picked
per label as whichever gives the flatter ratio over SZA 0-60 deg. All 137 labels of the table are
fitted at import (user decision 2026-11-07: SZA-dependent ratios, not fixed factors).

EASY form (a named parameter, evaluated in the settings pass like MCM's J<n> = PHTFKT(...)):
    KHV_G1700 = GPHOT(l_ratio, dm, -dn, jo1d, hvfact)
    k[X + HV --> ...] = CONST(KHV_G1700)
GPHOT needs COSX, AMAX (j_clear_sky.incl, or emit `HELPER_LINES`).
"""
from __future__ import annotations

import os
import re

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PHOT_TABLE = os.path.join(HERE, 'gecko_data', 'solarlight.phot')

# SAPHIR directly measures these -- use the real observation, no scaling.
MEASURED = {
    2: 'jo1d',       # O3 -> O1D
    4: 'jno2',       # NO2 -> NO + O3P
    11: 'jh2o2',     # H2O2 -> 2 OH
    12: 'jhono',     # HONO -> OH + NO
    17: 'jhcho_m',   # HCHO -> H2 + CO   (molecular channel)
    110: 'jhcho_r',  # HCHO -> 2 HO2 + CO (radical channel)
}
ANCHORS = {'jno2': 4, 'jo1d': 2}

# EASY definitions a scenario needs before the GPHOT(...) calls (COSX/AMAX come from
# j_clear_sky.incl in SAPHIR scenarios; the GECKO-A includes emit HELPER_LINES themselves).
HELPER_LINES = [
    "AMAX =FUNCTION[AMAX(%1,%2)]",
    "AMIN =FUNCTION[AMIN(%1,%2)]",
    "COSX=CONST(COS(SZ*0.0174533))",
]
# FUNCTION definitions must stand alone on their line (no trailing comment: the parsers' regex needs ']' last)
GPHOT_COMMENT = "; GPHOT(l_ratio, dm, -dn, anchor, hvfact) = GECKO-A J ratio l_ratio*cos^dm*exp(-dn*sec) x measured anchor x hvfact"
GPHOT_LINE = "GPHOT = FUNCTION[((%1)*((AMAX(COSX,0.01))@(%2))*EXP((%3)/AMAX(COSX,0.01)))*(%4)*(%5)]"


class GeckoPhotolysisError(Exception):
    pass


def load_tables(path: str = PHOT_TABLE) -> dict:
    """label -> (sza[deg], J[s-1]) from solarlight.phot (`PHOT <label> <n>` followed by n rows)."""
    tables = {}
    with open(path, encoding='latin-1') as fh:
        lines = fh.readlines()
    i = 0
    while i < len(lines):
        m = re.match(r'\s*PHOT\s+(\d+)\s+(\d+)', lines[i])
        if m:
            lab, n = int(m.group(1)), int(m.group(2))
            rows = [lines[i + 1 + k].split() for k in range(n)]
            tables[lab] = (np.array([float(r[0]) for r in rows]), np.array([float(r[1]) for r in rows]))
            i += 1 + n
        else:
            i += 1
    return tables


def fit_lmn(sza, j):
    """Least-squares fit of log J = log l + m*log cos(SZA) - n*sec(SZA) over the points with J > 0."""
    mask = (j > 0) & (sza < 90)
    th = np.radians(sza[mask])
    cos = np.cos(th)
    X = np.column_stack([np.ones_like(th), np.log(cos), -1.0 / cos])
    coef, *_ = np.linalg.lstsq(X, np.log(j[mask]), rcond=None)
    return float(np.exp(coef[0])), float(coef[1]), float(coef[2])


_TABLES = None
_FITS = None
_RATIO = None


def _ensure():
    global _TABLES, _FITS, _RATIO
    if _RATIO is not None:
        return
    _TABLES = load_tables()
    _FITS = {lab: fit_lmn(*v) for lab, v in _TABLES.items()}
    _RATIO = {}
    jref = {a: _TABLES[l][1] for a, l in ANCHORS.items()}

    def spread(js, ref):                      # flatness of the raw ratio over the first 7 SZA points (0..60 deg)
        r = js[:7] / ref[:7]
        return (r.max() - r.min()) / r.mean() if r.mean() > 0 else np.inf

    for lab, (sza, js) in _TABLES.items():
        if lab in MEASURED:
            continue
        anchor = min(ANCHORS, key=lambda a: spread(js, jref[a]))
        lr, mr, nr = _FITS[ANCHORS[anchor]]
        lx, mx, nx = _FITS[lab]
        _RATIO[lab] = (anchor, lx / lr, mx - mr, nx - nr)


def ratio_fit(label: int):
    """(anchor, l_ratio, dm, dn) with ratio(SZA) = l_ratio * cos(SZA)^dm * exp(-dn*sec(SZA))."""
    _ensure()
    if label not in _RATIO:
        raise GeckoPhotolysisError(f"HV label {label} not in {PHOT_TABLE}")
    return _RATIO[label]


def known_labels():
    _ensure()
    return sorted(set(MEASURED) | set(_RATIO))


def easy_rate(label: int, hvfact: float = 1.0) -> str:
    """Right-hand side of `KHV_G<label> = ...` for a GECKO-A HV label."""
    if label in MEASURED:
        j = MEASURED[label]
        return f"CONST({j})" if hvfact == 1.0 else f"CONST({j}*{hvfact:g})"
    anchor, l_ratio, dm, dn = ratio_fit(label)
    return f"GPHOT({l_ratio:.4g},{dm:.3f},{-dn:.3f},{anchor},{hvfact:g})"


def is_measured(label: int) -> bool:
    return label in MEASURED


def fit_quality(label: int):
    """(fit/raw at SZA 30, fit/raw at SZA 60) of the ratio expression against the raw tables."""
    _ensure()
    anchor, l_ratio, dm, dn = ratio_fit(label)
    sza, js = _TABLES[label]
    jref = _TABLES[ANCHORS[anchor]][1]
    out = []
    for target in (30.0, 60.0):
        k = int(np.argmin(np.abs(sza - target)))
        c = np.cos(np.radians(sza[k]))
        fitted = l_ratio * c ** dm * np.exp(-dn / c)
        raw = js[k] / jref[k] if jref[k] > 0 else np.nan
        out.append(fitted / raw if raw else np.nan)
    return tuple(out)


if __name__ == '__main__':
    _ensure()
    print(f"# {len(_TABLES)} labels in {PHOT_TABLE}; anchors: " +
          ', '.join(f"{a}=label {l} l={_FITS[l][0]:.5g} m={_FITS[l][1]:.4f} n={_FITS[l][2]:.4f}" for a, l in ANCHORS.items()))
    for lab in known_labels():
        if lab in MEASURED:
            print(f"{lab:>6}   measured   {easy_rate(lab)}")
        else:
            q30, q60 = fit_quality(lab)
            print(f"{lab:>6}  ratio-fit   {easy_rate(lab):50s}  fit/raw 30deg {q30:.3f}  60deg {q60:.3f}")
