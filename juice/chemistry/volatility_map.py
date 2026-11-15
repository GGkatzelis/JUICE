"""Where the species of an experiment sit in volatility and oxidation state, by each volatility method.

For every closed-shell organic species of the experiment's mechanism (the species the partitioning would use) the
C* of every method in volatility.METHODS, the carbon oxidation state and a structural group.  Drawn as clouds that
hold 80 % of the species (or of their weight, after a run), one per method, and split into groups for one method.
No Qt here: the window is juice/gui/volatility_window.py.

    OSc = 2 O/C - H/C - 5 N/C      Kroll et al. (2011), Nat. Chem. 3, 133; nitrogen counted as in a nitrate (+5),
                                   which holds for the organic nitrates and PANs that carry nearly all MCM nitrogen
    weight = (peak concentration in the run) ** p       p = 0 (each species once), 0.5 (default), 1 (by amount)
"""
from __future__ import annotations

import copy
import io
import contextlib
import math
import os
import shutil
import tempfile
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Tuple

import numpy as np

from juice.chemistry import volatility as vol

GROUPS = ["PANs", "organic nitrates", "aromatics", "hydroperoxides and peracids", "acids", "carbonyls",
          "alcohols and other oxygenates", "hydrocarbons", "other"]
_GROUP_SMARTS = [                                   # first match wins (the order of GROUPS)
    ("PANs", ["[CX3](=O)[OX2][OX2][NX3](=O)=O", "[CX3](=O)[OX2][OX2][N+](=O)[O-]"]),
    ("organic nitrates", ["[#6][OX2][NX3](=O)=O", "[#6][OX2][N+](=O)[O-]"]),
    ("aromatics", ["a"]),
    ("hydroperoxides and peracids", ["[OX2H][OX2]"]),
    ("acids", ["[CX3](=O)[OX2H1]"]),
    ("carbonyls", ["[CX3]=[OX1]"]),
]
METHOD_COLOURS = {"simpol": "#1f77b4", "nannoolal": "#2ca02c", "evaporation": "#17becf", "li_evap": "#9467bd",
                  "li2016": "#ff7f0e", "donahue2011": "#d62728"}
GROUP_COLOURS = {"PANs": "#d62728", "organic nitrates": "#ff7f0e", "aromatics": "#8c564b",
                 "hydroperoxides and peracids": "#e377c2", "acids": "#2ca02c", "carbonyls": "#1f77b4",
                 "alcohols and other oxygenates": "#17becf", "hydrocarbons": "#7f7f7f", "other": "#bcbd22"}
CLASS_BANDS = [(-99, math.log10(3e-4), "ELVOC", "#d9d9d9"), (math.log10(3e-4), math.log10(0.3), "LVOC", "#e8e8e8"),
               (math.log10(0.3), math.log10(300), "SVOC", "#f2f2f2"), (math.log10(300), math.log10(3e6), "IVOC", "#fafafa")]
MIN_CLOUD = 10                                      # fewer species than this: points, no cloud


def osc(nC: int, nH: int, nO: int, nN: int = 0) -> float:
    return (2.0 * nO - nH - 5.0 * nN) / nC if nC else float("nan")


_PATT: Dict[str, list] = {}


def group_of(smiles: str, nC: int, nO: int, nN: int) -> str:
    """The structural group of a species (GROUPS); from the formula alone when there is no structure."""
    if smiles:
        try:
            from rdkit import Chem
            mol = vol._mol(smiles)
            if not _PATT:
                for g, pats in _GROUP_SMARTS:
                    _PATT[g] = [Chem.MolFromSmarts(p) for p in pats]
            for g, _ in _GROUP_SMARTS:
                if any(mol.HasSubstructMatch(p) for p in _PATT[g]):
                    return g
        except Exception:
            pass
    if nO == 0 and nN == 0:
        return "hydrocarbons"
    return "alcohols and other oxygenates" if (nO and not nN and smiles) else "other"


@dataclass
class VolatilityMap:
    names: List[str]
    osc: np.ndarray
    groups: List[str]
    log_cstar: Dict[str, np.ndarray]            # method id -> log10 C* per species (nan where the method gives none)
    unavailable: Dict[str, str] = field(default_factory=dict)    # method id -> why it is not in the plot
    peak: Optional[np.ndarray] = None           # peak concentration of each species in a run (molec cm-3), or None
    run_file: str = ""
    coa_max: Optional[float] = None             # organic aerosol mass of the run at its maximum, ug m-3

    def weights(self, power: float) -> np.ndarray:
        if self.peak is None or power == 0:
            return np.ones(len(self.names))
        w = np.clip(np.nan_to_num(self.peak, nan=0.0), 0.0, None) ** power
        return w if w.sum() > 0 else np.ones(len(self.names))


def build(spec, run_file: str = "", progress: Optional[Callable[[str, float], None]] = None) -> VolatilityMap:
    """The map for an ExperimentSpec (MCM with chosen compounds): its mechanism subset is built in a temporary folder.
    `run_file`: the .cs of a run of this experiment, for the weights (optional)."""
    from juice.core import experiment_spec as es
    from juice.chemistry.species_db import load_species_db
    tmp = tempfile.mkdtemp(prefix="juice_volmap_")
    try:
        s = copy.deepcopy(spec)
        s.mechanism_stem, s.partitioning = "", False
        if progress:
            progress("Building the mechanism of the chosen compounds", 0.0)
        with contextlib.redirect_stdout(io.StringIO()):
            es.ensure_mechanism(s, tmp)
        chem = os.path.join(tmp, "chemistry")
        if not os.path.isdir(chem):
            raise ValueError("no mechanism: choose compounds on the Chemistry step")
        texts = [open(os.path.join(chem, f), encoding="latin-1").read() for f in sorted(os.listdir(chem))
                 if f.lower().endswith((".mech", ".incl")) and ".dil." not in f.lower()]
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    names_all = vol.mechanism_species(texts)
    db = load_species_db(base_id=spec.base_mechanism)
    return from_species(names_all, db, run_file, progress)


def from_species(names_all: List[str], db, run_file: str = "",
                 progress: Optional[Callable[[str, float], None]] = None) -> VolatilityMap:
    per: Dict[str, Dict[str, float]] = {}
    unavailable: Dict[str, str] = {}
    for k, (mid, m) in enumerate(vol.METHODS.items()):
        if progress:
            progress(f"C* by {m.label}" + (" (loading Open Babel and UManSysProp the first time)"
                                            if m.backend == "umansysprop" and "ok" not in vol._BACKEND else ""),
                     0.1 + 0.8 * k / len(vol.METHODS))
        ok, why = m.available()
        if not ok:
            unavailable[mid] = why
            continue
        r = vol.partitioning_species(names_all, db, mid, cstar_max=0.0)
        per[mid] = {s[0]: math.log10(float(s[5])) for s in r.species}
    names = sorted(set().union(*(set(v) for v in per.values()))) if per else []
    osc_v, groups = [], []
    for n in names:
        sp = db[n]
        smi = sp.smiles if sp.smiles not in ("-", "") else ""
        osc_v.append(osc(sp.nC, sp.nH, sp.nO, sp.nN))
        groups.append(group_of(smi, sp.nC, sp.nO, sp.nN))
    log_cstar = {mid: np.array([v.get(n, np.nan) for n in names]) for mid, v in per.items()}
    vm = VolatilityMap(names, np.array(osc_v), groups, log_cstar, unavailable)
    if run_file:
        load_run(vm, run_file)
    return vm


def load_run(vm: VolatilityMap, run_file: str) -> None:
    """Peak concentration of every species in a run (.cs), for the weights; COA if the run had partitioning."""
    from juice.core.column_file_reader import load_column_file
    df = load_column_file(run_file)
    df = df[0] if isinstance(df, tuple) else df
    cols = {c.upper(): c for c in df.columns}
    vm.peak = np.array([float(np.nanmax(df[cols[n]].values)) if n in cols else np.nan for n in vm.names])
    vm.run_file = run_file
    vm.coa_max = float(np.nanmax(df[cols["COA"]].values)) if "COA" in cols else None


# ── drawing (matplotlib, no Qt) ────────────────────────────────────────────────────────────────────────────────



def _hdr_level(z: np.ndarray, frac: float) -> float:
    """The density level whose contour holds `frac` of the density on the grid (highest-density region)."""
    v = np.sort(z.ravel())[::-1]
    c = np.cumsum(v)
    return float(v[min(np.searchsorted(c, frac * c[-1]), len(v) - 1)])


def density(x, y, w, xlim, ylim, grid=200):
    """Weighted Gaussian kernel density on a grid: the points binned onto it and smoothed with Scott's bandwidth
    (per axis, n_eff = (sum w)^2 / sum w^2), which is gaussian_kde without the cross term, at a cost independent of
    the number of species.  Returns (X, Y, Z) or None for too few or degenerate points."""
    from scipy.ndimage import gaussian_filter
    ok = np.isfinite(x) & np.isfinite(y) & (w > 0)
    x, y, w = x[ok], y[ok], w[ok]
    if len(x) < MIN_CLOUD or np.ptp(x) == 0 or np.ptp(y) == 0:
        return None
    p = w / w.sum()
    n_eff = 1.0 / np.sum(p ** 2)
    f = n_eff ** (-1.0 / 6.0)
    sx = f * math.sqrt(np.sum(p * (x - np.sum(p * x)) ** 2))
    sy = f * math.sqrt(np.sum(p * (y - np.sum(p * y)) ** 2))
    xe = np.linspace(xlim[0], xlim[1], grid + 1); ye = np.linspace(ylim[0], ylim[1], grid + 1)
    H, _, _ = np.histogram2d(y, x, bins=[ye, xe], weights=p)
    dx, dy = xe[1] - xe[0], ye[1] - ye[0]
    Z = gaussian_filter(H, sigma=(max(sy / dy, 0.5), max(sx / dx, 0.5)), mode="constant")
    X, Y = np.meshgrid((xe[:-1] + xe[1:]) / 2, (ye[:-1] + ye[1:]) / 2)
    return X, Y, Z


def cloud(ax, x, y, w, colour, lw=1.2, frac=0.8, label="", fill_alpha=0.0, grid=200):
    """Draw the contour holding `frac` of the weighted points (Gaussian kernel density). Returns True if drawn."""
    d = density(x, y, w, ax.get_xlim(), ax.get_ylim(), grid)
    if d is None:
        return False
    X, Y, Z = d
    lev = _hdr_level(Z, frac)
    if fill_alpha:
        ax.contourf(X, Y, Z, levels=[lev, Z.max() * 1.01], colors=[colour], alpha=fill_alpha)
    ax.contour(X, Y, Z, levels=[lev], colors=[colour], linewidths=lw)
    if label:
        ax.plot([], [], color=colour, lw=lw, label=label)
    return True


def ranges(vm: VolatilityMap) -> Tuple[Tuple[float, float], Tuple[float, float]]:
    """Axis ranges that hold every species of every method, with a margin for the clouds."""
    xs = np.concatenate([v[np.isfinite(v)] for v in vm.log_cstar.values()] or [np.array([0.0])])
    ys = vm.osc[np.isfinite(vm.osc)] if len(vm.osc) else np.array([0.0])
    x0, x1 = max(-16.0, math.floor(xs.min()) - 2), min(16.0, math.ceil(xs.max()) + 2)
    y0, y1 = min(-2.0, math.floor(ys.min() * 2) / 2 - 0.5), max(2.0, math.ceil(ys.max() * 2) / 2 + 0.5)
    return (x0, x1), (y0, y1)


def _frame(ax, vm: VolatilityMap, coa: Optional[float]):
    (x0, x1), (y0, y1) = ranges(vm)
    for lo, hi, name, col in CLASS_BANDS:
        lo, hi = max(lo, x0), min(hi, x1)
        if hi > lo:
            ax.axvspan(lo, hi, color=col, zorder=0)
            ax.text((lo + hi) / 2, y1 - 0.08 * (y1 - y0) / 4, name, ha="center", va="top", fontsize=8, color="#666")
    if coa:
        ax.axvline(math.log10(coa), color="k", ls="--", lw=0.8)
        ax.text(math.log10(coa), y0 + 0.1, f"  C_OA {coa:.3g} µg m⁻³", fontsize=8, rotation=90, va="bottom")
    ax.set_xlim(x0, x1); ax.set_ylim(y0, y1)
    ax.set_xlabel("log₁₀ C* (µg m⁻³, 298 K)"); ax.set_ylabel("OSc = 2 O/C − H/C − 5 N/C")


def _weight_text(vm: VolatilityMap, power: float) -> str:
    if vm.peak is None or power == 0:
        return "each species counted once"
    return {0.5: "weighted by √(peak concentration)", 1.0: "weighted by peak concentration"}.get(
        power, f"weighted by (peak concentration)^{power:g}") + f" in {os.path.basename(vm.run_file)}"


def draw_methods(fig, vm: VolatilityMap, chosen: str, power: float = 0.5, frac: float = 0.8):
    """One cloud per method; the chosen method bold, with its species as faint points."""
    fig.clear()
    ax = fig.add_subplot(111)
    _frame(ax, vm, vm.coa_max)
    w = vm.weights(power)
    order = [m for m in vm.log_cstar if m != chosen] + ([chosen] if chosen in vm.log_cstar else [])
    for mid in order:
        x = vm.log_cstar[mid]
        bold = mid == chosen
        lab = vol.METHODS[mid].label + (" (chosen)" if bold else "")
        if bold:
            ok = np.isfinite(x)
            ax.scatter(x[ok], vm.osc[ok], s=6 + 30 * (w[ok] / w.max()) ** 0.5 if vm.peak is not None and power else 8,
                       color=METHOD_COLOURS.get(mid, "k"), alpha=0.25, linewidths=0, zorder=2)
        if not cloud(ax, x, vm.osc, w, METHOD_COLOURS.get(mid, "k"), lw=3.0 if bold else 1.1, frac=frac, label=lab,
                     fill_alpha=0.10 if bold else 0.0):
            ax.plot([], [], color=METHOD_COLOURS.get(mid, "k"), lw=1, label=lab + " (too few species)")
    for mid, why in vm.unavailable.items():
        ax.plot([], [], color="#bbb", lw=1, ls=":", label=f"{vol.METHODS[mid].label} (not installed)")
    ax.legend(loc="lower left", fontsize=8, framealpha=0.9)
    ax.set_title(f"{len(vm.names)} species of the experiment; clouds hold {frac:.0%} of them ({_weight_text(vm, power)})",
                 fontsize=9)
    fig.tight_layout()
    return ax


def draw_groups(fig, vm: VolatilityMap, method_id: str, power: float = 0.5, frac: float = 0.8):
    """The cloud of one method split into structural groups; a group with few species as points."""
    fig.clear()
    ax = fig.add_subplot(111)
    _frame(ax, vm, vm.coa_max)
    if method_id not in vm.log_cstar:
        ax.set_title(f"{vol.METHODS[method_id].label}: not available ({vm.unavailable.get(method_id, '')[:80]})", fontsize=9)
        return ax
    x, w = vm.log_cstar[method_id], vm.weights(power)
    g = np.array(vm.groups)
    for grp in GROUPS:
        sel = (g == grp) & np.isfinite(x)
        n = int(sel.sum())
        if not n:
            continue
        col = GROUP_COLOURS[grp]
        size = (6 + 40 * (w[sel] / w.max()) ** 0.5) if vm.peak is not None and power else 10
        drawn = cloud(ax, x[sel], vm.osc[sel], w[sel], col, lw=2.0, frac=frac, label=f"{grp} ({n})", fill_alpha=0.08)
        ax.scatter(x[sel], vm.osc[sel], s=size, color=col, alpha=0.35 if drawn else 0.9, linewidths=0,
                   label=None if drawn else f"{grp} ({n}, points)", zorder=3)
    ax.legend(loc="lower left", fontsize=8, framealpha=0.9)
    ax.set_title(f"{vol.METHODS[method_id].label}: by group; clouds hold {frac:.0%} ({_weight_text(vm, power)})", fontsize=9)
    fig.tight_layout()
    return ax
