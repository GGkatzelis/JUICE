"""Where the organic carbon and oxygen of a run are: in the gas, in the particles, on the chamber wall.

For a run with gas–particle–wall partitioning (experiment_spec.partitioning or a GECKO-A base with its partitioning
include) every condensable species X has three companions in the results: SOA_X (in the particles), WALL_X (vapour
taken up by the wall, reversible) and WALL_PL_X (particles lost to the wall).  This module sums them, weighted by
the carbon or oxygen atoms of X, over time and, at one time, by the volatility of X (from the partitioning include
of the run).  No Qt here: the Phases tab of the results viewer draws it.

    mass of element E in a phase = sum_X n_E(X) * [X in that phase] * M_E / N_A      in ug m-3
    gas = every organic species of the run (n_C >= 1, not CO / CO2), partitioning or not.
    flushed out = the integral of DILUTE x (gas + particles) over time, when the run has a DILUTE column: gas and
    particles leave with the replenishment flow, the wall stays.  Gas + particles + wall + flushed out is then
    the organic carbon put into the chamber (less what became CO and CO2).
"""
from __future__ import annotations

import glob
import math
import os
from dataclasses import dataclass, field
from typing import Dict, List, Optional

import numpy as np

N_A = 6.02214076e23
ATOM_MASS = {"C": 12.011, "O": 15.999}
PHASES = ["gas", "particles", "wall (vapour)", "wall (particles)"]
DILUTED = "flushed out (cumulative)"
PHASE_COLOURS = {"gas": "#56B4E9", "particles": "#D55E00", "wall (vapour)": "#999999", "wall (particles)": "#555555",
                 DILUTED: "#E6E6E6"}
PREFIX = {"particles": "SOA_", "wall (vapour)": "WALL_", "wall (particles)": "WALL_PL_"}
NOT_ORGANIC = {"CO", "CO2"}
SUMS = ("SVOC", "LVOC", "ELVOC")                    # SOA_SVOC_C ... are diagnostic sums, not species


def ug_per_molec(element: str) -> float:
    """ug m-3 per (molecule cm-3 x atom)."""
    return ATOM_MASS[element] / N_A * 1e12


def soa_include_of_run(cs_path: str) -> str:
    """The partitioning include a run used: <run>/chemistry/*.SOA.include (MCM: <name>.SOA.include, GECKO-A:
    <base>.SOA.include).  '' if the run has none."""
    hits = sorted(glob.glob(os.path.join(os.path.dirname(os.path.abspath(cs_path)), "chemistry", "*.SOA.include")))
    return hits[0] if hits else ""


@dataclass
class PhaseBudget:
    t_hours: np.ndarray
    species: List[str]                             # the partitioning species X
    log_cstar: Dict[str, float]                    # X -> log10 C* (ug m-3), from the include
    mass: Dict[str, Dict[str, np.ndarray]]         # element -> phase -> ug m-3 over time
    per_species: Dict[str, Dict[str, np.ndarray]]  # element -> phase -> array (time x species) in ug m-3
    note: str = ""

    def at(self, element: str, i: int) -> Dict[str, float]:
        return {ph: float(self.mass[element][ph][i]) for ph in PHASES + [DILUTED] if ph in self.mass[element]}

    def by_volatility(self, element: str, i: int, bins: Optional[np.ndarray] = None) -> Dict[str, np.ndarray]:
        """Element mass per phase at output index i, in decade bins of log10 C* (partitioning species only)."""
        bins = np.arange(-12, 13) if bins is None else bins
        lc = np.array([self.log_cstar.get(x, np.nan) for x in self.species])
        idx = np.digitize(lc, bins) - 1
        out = {}
        for ph in PHASES:
            v = self.per_species[element][ph][i]
            h = np.zeros(len(bins) - 1)
            ok = np.isfinite(lc) & (idx >= 0) & (idx < len(h))
            np.add.at(h, idx[ok], v[ok])
            out[ph] = h
        return out


def build(df, db, cs_path: str = "", soa_include: str = "", t_seconds=None) -> PhaseBudget:
    """df: the run's results (columns = species in molecules cm-3); t_seconds: the output times (default the index
    of df, in seconds, as column_file_reader gives it; juice.Results: pass results.hour * 3600); db: a
    species_db.SpeciesDB with the run's species (formulas).  Raises ValueError for a run without partitioning."""
    from juice.core.particles import read_soa_include
    cols = {c.upper(): c for c in df.columns}
    inc = soa_include or (soa_include_of_run(cs_path) if cs_path else "")
    lc = {n.upper(): math.log10(c) for n, c, _cw in read_soa_include(inc)["species"] if c > 0} if inc else {}
    part = sorted(x for x in {c[4:] for c in cols if c.startswith("SOA_")} if x in cols and x.split("_")[0] not in SUMS)
    if not part:
        raise ValueError("this run has no gas-particle-wall partitioning (no SOA_ columns)")
    t = np.asarray(df.index if t_seconds is None else t_seconds, dtype=float)
    t_hours = (t - t[0]) / 3600.0
    counts = {}
    for n in cols:
        if n in NOT_ORGANIC or n.startswith(("SOA_", "WALL_")) or n not in db:
            continue
        sp = db[n]
        if sp.nC >= 1:
            counts[n] = (sp.nC, sp.nO)
    for x in part:                                  # a partitioning species the table lacks: from its own entry
        if x not in counts and x in db:
            counts[x] = (db[x].nC, db[x].nO)
    part = [x for x in part if x in counts]
    mass: Dict[str, Dict[str, np.ndarray]] = {}
    per: Dict[str, Dict[str, np.ndarray]] = {}
    for el, k in (("C", 0), ("O", 1)):
        f = ug_per_molec(el)
        gas_all = np.zeros(len(t))
        for n, cnt in counts.items():
            if cnt[k]:
                gas_all += cnt[k] * np.nan_to_num(df[cols[n]].to_numpy(dtype=float)) * f
        mass[el] = {"gas": gas_all}
        per[el] = {}
        n_el = np.array([counts[x][k] for x in part], dtype=float)
        per[el]["gas"] = np.column_stack([np.nan_to_num(df[cols[x]].to_numpy(dtype=float)) for x in part]) * n_el * f
        for ph, pre in PREFIX.items():
            arr = np.column_stack([np.nan_to_num(df[cols[pre + x]].to_numpy(dtype=float)) if pre + x in cols
                                   else np.zeros(len(t)) for x in part]) * n_el * f
            per[el][ph] = arr
            mass[el][ph] = arr.sum(axis=1)
        if "DILUTE" in cols and len(t) > 1:           # what the flow carried out of the chamber, cumulative
            rate = np.nan_to_num(df[cols["DILUTE"]].to_numpy(dtype=float)) * (mass[el]["gas"] + mass[el]["particles"])
            mass[el][DILUTED] = np.concatenate([[0.0], np.cumsum(0.5 * (rate[1:] + rate[:-1]) * np.diff(t))])
    missing = sum(1 for x in part if x not in lc)
    note = (f"{len(part)} partitioning species" + (f"; {missing} without a C* in {os.path.basename(inc)}" if inc and missing else "")
            + ("" if inc else "; no partitioning include found next to the run: no volatility axis"))
    return PhaseBudget(t_hours, part, {x: lc[x] for x in part if x in lc}, mass, per, note)


# ── drawing (matplotlib, no Qt) ─────────────────────────────────────────────────────────────────────────────

def draw(fig, pb: PhaseBudget, element: str = "C", i: Optional[int] = None):
    """(a) the element in each phase over the run, stacked; (b) at output index i (default the last), the
    partitioning species' element in each phase by decade of C*; the O:C of each phase in the legend."""
    fig.clear()
    i = len(pb.t_hours) - 1 if i is None else int(i)
    ax1 = fig.add_subplot(1, 2, 1)
    ax2 = fig.add_subplot(1, 2, 2)
    m = pb.mass[element]
    shown = PHASES + ([DILUTED] if DILUTED in m else [])
    ax1.stackplot(pb.t_hours, *[m[ph] for ph in shown], colors=[PHASE_COLOURS[p] for p in shown],
                  labels=[_label(pb, ph, i) for ph in shown], alpha=0.9)
    ax1.axvline(pb.t_hours[i], color="k", lw=0.8, ls="--")
    ax1.set_xlabel("hours since the start"); ax1.set_ylabel(f"organic {element} (µg m⁻³)")
    ax1.set_title(f"(a) organic {'carbon' if element == 'C' else 'oxygen'} by phase", loc="left", fontsize=9)
    ax1.legend(loc="lower left", fontsize=7, framealpha=0.9)
    ax1.set_xlim(pb.t_hours[0], pb.t_hours[-1])
    if pb.log_cstar:
        bins = np.arange(-12, 13)
        h = pb.by_volatility(element, i, bins)
        x = bins[:-1] + 0.5
        bottom = np.zeros(len(x))
        for ph in PHASES[1:] + PHASES[:1]:       # condensed phases first, the gas on top
            ax2.bar(x, h[ph], width=0.85, bottom=bottom, color=PHASE_COLOURS[ph], label=ph)
            bottom += h[ph]
        used = np.nonzero(bottom)[0]
        if used.size:
            ax2.set_xlim(x[used[0]] - 1, x[used[-1]] + 1)
        ax2.set_xlabel("log₁₀ C* (µg m⁻³)"); ax2.set_ylabel(f"{element} of the partitioning species (µg m⁻³)")
        ax2.legend(loc="upper left", fontsize=7)
    else:
        ax2.text(0.5, 0.5, "no C* for the species of this run", ha="center", va="center", transform=ax2.transAxes)
    ax2.set_title(f"(b) by volatility at {pb.t_hours[i]:.2f} h", loc="left", fontsize=9)
    fig.tight_layout()
    return ax1, ax2


def _label(pb: PhaseBudget, ph: str, i: int) -> str:
    c, o = pb.mass["C"][ph][i], pb.mass["O"][ph][i]
    oc = (o / ATOM_MASS["O"]) / (c / ATOM_MASS["C"]) if c > 0 else float("nan")
    return f"{ph}: {c:.3g} µg C m⁻³" + (f", O:C {oc:.2f}" if np.isfinite(oc) else "")
