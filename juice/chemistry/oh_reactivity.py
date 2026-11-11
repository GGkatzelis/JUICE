"""
oh_reactivity.py — OH reactivity an injection adds, from the MCM's own rate expressions (no Qt).

    k_oh(species, T, p)            → Σ k of all 'species + OH' reactions of the base MCM [cm3 s-1]
    injection_koh(products, ppb, T, p) → OH reactivity of the injected amount [s-1], as if fully mixed
                                       and nothing had reacted yet (the reactivity just after the injection)

The complex coefficients (KMT05, KMT07, KMT15 …) are evaluated from the 'Generic Rate Coefficients'
block of the same .fac file, so the numbers are the ones the model uses.  Species without an OH
reaction in the base MCM (add-on species, tracers such as XNOY) contribute 0 and are reported.
"""
from __future__ import annotations

import math
import re
from functools import lru_cache
from typing import Dict, List, Optional, Tuple

PPB = 1e-9
_DEF = re.compile(r"^\s*([A-Z][A-Z0-9_]*)\s*=\s*(.+?)\s*;\s*$")
_FUNCS = {"EXP": math.exp, "LOG10": math.log10, "LOG": math.log, "SQRT": math.sqrt, "ABS": abs,
          "exp": math.exp, "log10": math.log10, "log": math.log, "sqrt": math.sqrt}


def _py(expr: str) -> str:
    e = re.sub(r"(\d)[Dd]([+-]?\d)", r"\1e\2", expr.strip())
    return e.replace("@", "**")


def number_density(T: float, p_hPa: float) -> float:
    """Air molecules per cm3."""
    return p_hPa * 100.0 / (1.380649e-23 * T) * 1e-6


@lru_cache(maxsize=4)
def _generic_defs(fac_path: str) -> Tuple[Tuple[str, str], ...]:
    """NAME = expression lines of the .fac before the first reaction (the generic/complex coefficients)."""
    out = []
    with open(fac_path, encoding="latin-1") as fh:
        for line in fh:
            if line.lstrip().startswith("%"):
                break
            m = _DEF.match(line)
            if m:
                out.append((m.group(1), _py(m.group(2))))
    return tuple(out)


@lru_cache(maxsize=64)
def _context(fac_path: str, T: float, p_hPa: float) -> Dict[str, float]:
    M = number_density(T, p_hPa)
    ctx: Dict[str, float] = dict(_FUNCS)
    ctx.update({"TEMP": T, "T": T, "M": M, "O2": 0.2095 * M, "N2": 0.7809 * M, "H2O": 0.01 * M, "RO2": 0.0})
    for name, expr in _generic_defs(fac_path):
        try:
            ctx[name] = float(eval(expr, {"__builtins__": {}}, ctx))
        except Exception:
            pass                               # e.g. a coefficient that needs a species concentration
    return ctx


@lru_cache(maxsize=4096)
def k_oh(species: str, T: float = 298.0, p_hPa: float = 1013.25, fac_path: str = "") -> Optional[float]:
    """Total OH rate coefficient of a species in the base MCM [cm3 s-1]; None if it has no OH reaction there."""
    from juice.chemistry.mcm_species_db import reactions_of
    if not fac_path:
        from juice.core.experiment_spec import base_fac
        fac_path = base_fac()
    ctx = _context(fac_path, round(T, 1), round(p_hPa, 1))
    total, found = 0.0, False
    for r in reactions_of(species, fac_path):
        if r.get("partner") != "OH":
            continue
        try:
            total += float(eval(_py(r["rate"]), {"__builtins__": {}}, ctx)); found = True
        except Exception:
            continue
    return total if found else None


def parse_products(products: str) -> List[Tuple[float, str]]:
    """'0.9*A+0.1*B' / 'A' / 'HONO+XNOy' → [(fraction, species)]; a term without a number counts 1."""
    out = []
    for term in re.split(r"\s*\+\s*", (products or "").strip()):
        if not term:
            continue
        m = re.match(r"^\(?\s*([0-9.eE+-]+)\s*\)?\s*\*\s*(.+)$", term)
        frac, sp = (float(m.group(1)), m.group(2)) if m else (1.0, term)
        out.append((frac, sp.strip().strip("()").upper()))
    return out


def injection_koh(products: str, ppb, T: float = 298.0, p_hPa: float = 1013.25, fac_path: str = "") -> Tuple[Optional[float], List[str]]:
    """(OH reactivity of the injected amount [s-1], species without an OH reaction in the MCM)."""
    try:
        amount = float(str(ppb).strip())
    except (TypeError, ValueError):
        return None, []
    terms = parse_products(products)
    if not terms:
        return None, []
    conc = amount * PPB * number_density(T, p_hPa)
    total, missing, any_k = 0.0, [], False
    for frac, sp in terms:
        k = k_oh(sp, T, p_hPa, fac_path)
        if k is None:
            if sp not in ("XNOY", "NA", "SA", "H2O", "CO2"):
                missing.append(sp)
            continue
        total += frac * k * conc; any_k = True
    return (total if any_k else 0.0), missing
