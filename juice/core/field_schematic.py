"""
field_schematic.py — a field case as one picture: what goes into the box, and how each input is defined (no Qt).

The picture is drawn from the case itself (`FieldSpec` + its table), so it is right for any site or city:
the air (a parcel on a trajectory, a site, a plume) in the middle with its chemistry, and one card per
input around it.  A card says WHAT the case provides (species, ranges read from the table), HOW the model
uses it (the term of the equation, units, interpolation) and WHERE it comes from (`FieldSpec.sources`,
free text).  Inputs the case does not use are drawn grey, so it is also clear what is left out.  The
equation of the run stands underneath, each term in the colour of its card.

The middle also shows what the box does: it stays over a site while the air passes through, or it travels
with a plume or along a trajectory (cubes at successive positions).

    cards(spec, table)                 the content as data: [Card]; table None = what the kind of run needs
    draw(spec, table, path="x.png")    the figure (matplotlib); saved when a path is given
    draw(spec, table, fig=fig)         the same into an existing figure (the field window)
    hit(x, y)                          the card at a point of the picture ("centre" for the box, "" for none)
    python -m juice.core.field_schematic case.field.json [out.png]
    python -m juice.core.field_schematic --examples [folder]      one PNG per bundled example
"""
from __future__ import annotations

import os
import sys
import textwrap
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np

from juice.core import field_spec as fs

INK, MUTED, FAINT, SURFACE = "#0b0b0b", "#6b6a64", "#c3c2b7", "#ffffff"
COLOUR = {"conditions": "#52514e", "photolysis": "#b98300", "emissions": "#eb6834", "mixing": "#2a78d6",
          "constraints": "#6250d6", "deposition": "#8a5a2b", "observations": "#1baf7a", "aerosol": "#c2488f"}
TITLE = {"conditions": "Meteorology", "photolysis": "Sunlight", "emissions": "Emissions", "mixing": "Mixing with background air",
         "constraints": "Held to measurements", "deposition": "Deposition", "observations": "Compared with observations",
         "aerosol": "Aerosol: gas–particle partitioning"}
BOX = {"constrained": "Stationary box at one site", "plume": "Moving box: follows the plume",
       "lagrangian": "Moving box: follows the air parcel"}
SCENE_NOTE = {"constrained": "the same place at every time step; the air passes through",
              "plume": "time since emission  →   the box widens as the plume dilutes",
              "lagrangian": "along the trajectory  →   the box grows with the boundary layer"}
WRAP = 53                                    # characters per line of a card
MAX_LINES = 6
WIDE = 21                                    # more characters in the wider card at the top


@dataclass
class Card:
    key: str
    used: bool = True
    lines: List[str] = field(default_factory=list)     # what the case provides and how the model uses it
    source: str = ""                                    # where it comes from ("" = not stated)
    heading: str = ""

    @property
    def title(self) -> str:
        return self.heading or TITLE[self.key]


# ── content ──────────────────────────────────────────────────────────────────────────────────────
def _names(species, limit: int = 6) -> str:
    sp = list(species)
    return ", ".join(sp[:limit]) + (f" + {len(sp) - limit} more" if len(sp) > limit else "")


def _rng(v: np.ndarray, fmt: str = "{:.0f}", factor: float = 1.0) -> str:
    v = np.asarray(v, dtype=float) * factor
    v = v[np.isfinite(v)]
    if not len(v):
        return "?"
    a, b = fmt.format(float(v.min())), fmt.format(float(v.max()))
    return a if a == b else f"{a} to {b}"


def _step_text(table: fs.FieldTable) -> str:
    if table.n < 2:
        return "one row"
    dt = float(np.median(np.diff(table.times_s)))
    return f"one row every {dt / 60.0:g} min" if dt < 5400 else f"one row every {dt / 3600.0:g} h"


def _when(table: fs.FieldTable) -> str:
    import pandas as pd
    t = pd.to_datetime(np.asarray(table.times_s, dtype=float), unit="s", origin="2000-01-01")
    a, b = t[0], t[-1]
    return (f"{a:%d %b %Y %H:%M} to {b:%H:%M} UTC" if a.date() == b.date() else f"{a:%d %b %Y %H:%M} to {b:%d %b %Y %H:%M} UTC")


def needs(spec: fs.FieldSpec) -> List[Card]:
    """The cards before a table is chosen: what this kind of run takes, and the columns that define it."""
    rt, off = spec.run_type, set(spec.off or [])
    moving = rt != "constrained"
    text = {
        "conditions": (True, ["Columns T, p and H2O or RH, one row per time", "lat and lon per row along the path" if moving else
                              "Position of the site: lat, lon of the case", "Column h: boundary-layer height" if rt == "lagrangian" or spec.dilution == "pbl" else ""]),
        "photolysis": ("photolysis" not in off, ["Column jNO2 if it was measured: it scales the clear-sky photolysis of the mechanism",
                                                 "Without it: clear sky for time and position"]),
        "emissions": (rt == "lagrangian" and "emissions" not in off, ["One E_<species> column per emitted species, unit in brackets",
                                                                      f"Default unit of this case: {spec.default_emission_unit}"]),
        "mixing": (spec.dilution != "none" or rt == "constrained",
                   [f"Dilution: {fs.DILUTION_SOURCES.get(spec.dilution, spec.dilution)}", "One B_<species> column per species of the background air"
                    if "background" not in off else "Background switched off: diluted towards clean air",
                    f"Species that are not held are lost with a lifetime of {spec.generic_lifetime_h:g} h" if rt == "constrained" else ""]),
        "constraints": (rt == "constrained" and "constraints" not in off, ["One C_<species> column per measured species the model is held to"]),
        "deposition": (any(fs.deposition_of(spec)), ["Velocity or lifetime per species, set on the Mixing and deposition step"]),
        "observations": ("observations" not in off, ["One O_<species> column per measured species to compare the model with"]),
        "aerosol": (fs.partitioning_on(spec), ["Column CS: condensation sink of the ambient aerosol, 1/s", "Column OA: organic aerosol mass, µg/m³",
                                               "No chamber walls in a field case"] if fs.partitioning_on(spec) else
                    ["No partitioning: every species stays in the gas phase", "Switched on at the Chemistry step"]),
    }
    out = []
    for key in TITLE:
        used, lines = text[key]
        if key in off or (key == "mixing" and False):
            lines = ["Switched off in this case"]
        out.append(Card(key, used, [l for l in lines if l], (spec.sources or {}).get(key, "")))
    return out


def cards(spec: fs.FieldSpec, table: Optional[fs.FieldTable]) -> List[Card]:
    if table is None:
        return needs(spec)
    c, src = table.conditions, {k.lower(): v for k, v in (spec.sources or {}).items()}
    out: List[Card] = []
    gone = table.switched_off

    def switched(key: str) -> List[str]:
        return [f"Switched off in this case: {len(gone[key])} column(s) of the table are ignored ({_names(gone[key], 4)})"]

    def add(key: str, used: bool, lines: List[str]) -> None:
        out.append(Card(key, used, [l for l in lines if l], src.get(key, "")))

    # meteorology
    L = [f"T {_rng(c['t'])} K, p {_rng(c['p'])} hPa, H2O {_rng(c['H2O'], '{:.1f}')} %" if all(k in c for k in ("t", "p", "H2O")) else "",
         f"Boundary-layer height {_rng(c['h'])} m (column h)" if "h" in c else "",
         f"Table: {_step_text(table)}, interpolated to the {spec.step_s:g} s model step"]
    if "lat" in c and "lon" in c and spec.run_type != "constrained":
        L.append(f"Path: {c['lat'][0]:.2f}°, {c['lon'][0]:.2f}° to {c['lat'][-1]:.2f}°, {c['lon'][-1]:.2f}°" + (f" ({src['trajectory']})" if src.get("trajectory") else ""))
    add("conditions", True, L)

    # sunlight
    j = sorted(table.photolysis)
    sz = ("solar zenith angle from the table" if "sz" in c and np.isfinite(c["sz"]).all() else
          "solar zenith angle from time and position" + (" along the path" if "lat" in c and "lon" in c and spec.run_type != "constrained" else ""))
    if "NO2" in table.photolysis:
        others = [n for n in j if n != "NO2"]
        L = [f"Measured jNO2 ({_rng(table.photolysis['NO2'], '{:.1f}', 1e3)} ×10⁻³ s⁻¹) scales every clear-sky photolysis frequency of the mechanism (factor at most 1)",
             f"Also taken from the table: {_names('j' + n for n in others)}" if others else "", "Clear sky: " + sz]
    else:
        L = ["Clear-sky photolysis of the mechanism (no jNO2 column in the table)", "Clear sky: " + sz,
             f"In the table but not used: {_names('j' + n for n in j)}" if j else ""]
    if "photolysis" in gone:
        L = ["Clear-sky photolysis of the mechanism", "Clear sky: " + sz, "Measured photolysis switched off: " + _names("j" + n for n in gone["photolysis"])]
    add("photolysis", True, L)

    # emissions
    em = sorted(table.emissions)
    if em:
        groups = {g: [s for s in em if fs.emission_class(s) == g] for g in ("NOX", "VOC", "OTHER")}
        units = sorted({fs.resolve_emission_unit(table, spec, s) for s in em})
        how = []
        for u in units:
            needs_h, needs_a, _f = fs.EMISSION_UNITS.get(u, (False, False, 1.0))
            how.append(f"{u}: " + (f"E / (h · A), footprint A = {spec.area_km2:g} km²" if needs_a else "E / h, mixed over the boundary layer" if needs_h
                                   else "added directly to the air"))
        L = [f"{len(em)} species: " + " · ".join(f"{lab}: {_names(groups[g], 4)}" for g, lab in (("NOX", "NOx"), ("VOC", "VOC"), ("OTHER", "other")) if groups[g]),
             *how, "Held constant from one table row to the next",
             "Scaled in this case: " + ", ".join(f"{k} × {v:g}" for k, v in spec.emission_scale.items()) if spec.emission_scale else ""]
    else:
        L = switched("emissions") if "emissions" in gone else ["No emissions (no E_ columns in the table)"]
    add("emissions", bool(em), L)

    # mixing
    bg = sorted(table.background)
    if spec.dilution != "none":
        try:
            k = fs.dilution_rate(table, spec) * 3600.0
            rate = f"{_rng(k, '{:.2f}')} per hour"
        except Exception:
            rate = "rate not available"
        how = {"pbl": f"Rate from the growth of the boundary layer, d ln h / dt; {spec.k_collapse_per_h:g} per hour while it does not grow",
               "table": "Rate from the table (column kdil)", "constant": "One constant rate",
               "tracer": f"Rate from the decay of measured {spec.tracer} above its background"}[spec.dilution]
        L = [f"{how}: {rate}", f"Towards the background (B_ columns) of {_names(bg)}" if bg else ("Towards clean air: background switched off" if "background" in gone else "Towards clean air (no B_ columns)"),
             "Every other species is diluted towards zero" if bg else ""]
    else:
        L = ["No dilution in this case"]
    if spec.run_type == "constrained":
        L.append(f"Every species that is not held is lost with a lifetime of {spec.generic_lifetime_h:g} h (stands for transport and deposition)")
    add("mixing", spec.dilution != "none" or spec.run_type == "constrained", L)
    if spec.dilution == "none" and spec.run_type == "constrained":
        out[-1].heading = "Transport and other losses"

    # held species
    held = sorted(table.constraints)
    add("constraints", bool(held), switched("constraints") if "constraints" in gone else [f"{len(held)} species set to the measured value at every step: {_names(held)}" if held else
                                    "Nothing is held: every species is calculated", "Radicals and all products are calculated" if held else "",
                                    f"Held by your choice, from O_ columns: {_names(table.held_by_choice)}" if table.held_by_choice else ""])

    # deposition
    velocity, lifetime = fs.deposition_of(spec)
    dep = [f"{s} {v:g} cm/s" for s, v in velocity.items()]
    life = [f"{s} lifetime {t:g} h" for s, t in lifetime.items()]
    n_set = len(spec.deposition_velocity_cm_s) + len(spec.deposition_lifetime_h)
    add("deposition", bool(dep or life), [f"Switched off in this case: {n_set} species are set on the Mixing and deposition step"]
        if "deposition" in (spec.off or []) and n_set else ["Loss v / h over the boundary layer: " + ", ".join(dep) if dep else "",
                                          "First-order loss: " + ", ".join(life) if life else "",
                                          "" if dep or life else "No deposition in this case"])

    # aerosol
    if fs.partitioning_on(spec):
        cs, oa = fs.condensation_sink(table, spec), fs.organic_aerosol(table, spec)
        method = (spec.chemistry or {}).get("volatility_method", "") if (spec.chemistry or {}).get("partitioning") else ""
        L = [f"Condensation sink {_rng(cs, '{:.1f}', 1e3)} ×10⁻³ s⁻¹ ({'column CS' if 'cs' in c else 'value of the case'})"
             if (cs > 0).any() else "No condensation sink given: nothing condenses",
             f"Organic aerosol {_rng(oa, '{:.1f}')} µg/m³ ({'column OA' if 'oa' in c else 'value of the case'}) + what forms",
             " · ".join(x for x in (f"Volatility: {method}" if method else "",
                                    f"particle lifetime {spec.particle_lifetime_h:g} h" if spec.particle_lifetime_h > 0 else "") if x),
             "No chamber walls: the wall terms are zero"]
    else:
        L = ["No partitioning: every species stays in the gas phase"]
    add("aerosol", fs.partitioning_on(spec), L)

    # observations
    obs = sorted(table.observations)
    add("observations", bool(obs or held), switched("observations") if "observations" in gone and not held else [f"Measured, not used by the model (O_ columns): {_names(obs)}" if obs else "No O_ columns in the table",
                                            "Bias, NMB and R² are written after the run" if obs else "",
                                            "The held species are shown with the run as well" if held else "",
                                            f"Not held by your choice, compared only: {_names(table.released)}" if table.released else ""])
    return out


_NO_TABLE = fs.FieldTable(path="", times_s=np.array([]))


def centre_lines(spec: fs.FieldSpec, table: Optional[fs.FieldTable]) -> List[str]:
    table = table if table is not None else _NO_TABLE
    L = ["Chemistry: " + fs.chemistry_label(spec)]
    if spec.initials:
        L.append("Starts with (ppb): " + ", ".join(f"{k} {v:g}" for k, v in list(spec.initials.items())[:8])
                 + (f" + {len(spec.initials) - 8} more" if len(spec.initials) > 8 else ""))
    else:
        L.append("Starts from zero" + (" and the held species" if table.constraints else ""))
    if spec.run_type == "constrained" and spec.spin_up_days > 0:
        L.append(f"Spin-up: the first day repeated for {spec.spin_up_days:g} day(s) before the output")
    L.append(f"{_when(table)}, model step {spec.step_s:g} s" if table.n else "No input table chosen yet")
    return L


def equation(spec: fs.FieldSpec, table: fs.FieldTable) -> List[Tuple[str, str]]:
    """[(mathtext piece, colour)]: only the terms this case uses."""
    table = table if table is not None else _NO_TABLE
    parts = [(r"$\dfrac{dC}{dt}\;=\;$", INK), (r"  chemistry ", INK)]
    if table.emissions:
        units = {fs.resolve_emission_unit(table, spec, s) for s in table.emissions}
        area = any(fs.EMISSION_UNITS.get(u, (0, 0, 0))[1] for u in units)
        height = any(fs.EMISSION_UNITS.get(u, (0, 0, 0))[0] for u in units)
        parts.append((r"$\;+\;\dfrac{E}{h \cdot A}$" if area else r"$\;+\;\dfrac{E}{h}$" if height else r"$\;+\;E$", COLOUR["emissions"]))
    if spec.dilution != "none":
        parts.append((r"$\;-\;k_{dil}\,(C - C_{background})$", COLOUR["mixing"]))
    if spec.run_type == "constrained":
        parts.append((r"$\;-\;\dfrac{C}{\tau}$", COLOUR["mixing"]))
    if any(fs.deposition_of(spec)):
        parts.append((r"$\;-\;\dfrac{v_{d}}{h}\,C$" if fs.deposition_of(spec)[0] else r"$\;-\;\dfrac{C}{\tau_{dep}}$", COLOUR["deposition"]))
    if fs.partitioning_on(spec):
        parts.append((r"$\;-\;k_{CS}\,(C - C_{eq})$", COLOUR["aerosol"]))
    if table.constraints:
        parts.append((r"$\qquad C = C_{measured}$ for the held species", COLOUR["constraints"]))
    return parts


# ── drawing ──────────────────────────────────────────────────────────────────────────────────────
# card: (x0, y0, width, height); arrow: (from, to) in the same coordinates (x 0..100, y -9..67)
_SLOT = {"photolysis": (31, 47.5, 38, 12.5), "conditions": (0.5, 31.5, 29, 15), "constraints": (0.5, 14.5, 29, 15),
         "mixing": (70.5, 31.5, 29, 15), "observations": (70.5, 14.5, 29, 15), "emissions": (4.5, 0, 29, 13), "aerosol": (35.5, 0, 29, 13),
         "deposition": (66.5, 0, 29, 13)}
_CENTRE = (34.5, 18.5, 31, 25)
_ARROW = {"photolysis": ((50, 47.5), (50, 43.5), "->"), "conditions": ((29.5, 39), (34.5, 37), "->"), "constraints": ((29.5, 24), (34.5, 27), "->"),
          "mixing": ((70.5, 39), (65.5, 37), "<->"), "observations": ((65.5, 27), (70.5, 24), "->"), "emissions": ((30, 13), (39.5, 18.5), "->"),
          "aerosol": ((50, 13), (50, 18.5), "<->"), "deposition": ((60.5, 18.5), (70, 13), "->")}
_ASPECT = 1.16                               # y units per x unit of the same length on paper


FIGSIZE = (13.2, 8.6)
_AXES = (0.012, 0.01, 0.976, 0.98)
_XLIM, _YLIM = (0.0, 100.0), (-9.0, 67.0)


def to_data(fx: float, fy: float) -> Tuple[float, float]:
    """A point of the picture as fractions of its width and height (from the top left) → data coordinates."""
    x = (fx - _AXES[0]) / _AXES[2] * (_XLIM[1] - _XLIM[0]) + _XLIM[0]
    y = ((1.0 - fy) - _AXES[1]) / _AXES[3] * (_YLIM[1] - _YLIM[0]) + _YLIM[0]
    return x, y


def hit(x: float, y: float) -> str:
    """The card at a point of the picture (data coordinates): its key, 'centre' for the box, '' for none."""
    for key, (x0, y0, w, h) in {**_SLOT, "centre": _CENTRE}.items():
        if x0 <= x <= x0 + w and y0 <= y <= y0 + h:
            return key
    return ""


def _wrap(lines: List[str], width: int, limit: int) -> List[str]:
    out: List[str] = []
    for l in lines:
        out += textwrap.wrap(l, width, subsequent_indent="  ") or [""]
    if len(out) > limit:
        out = out[:limit]
        out[-1] = out[-1][:width - 2].rstrip() + " …"
    return out


def _cube(ax, x: float, y: float, w: float, h: float, lw: float = 1.5) -> None:
    """A box drawn as a cube: front face (x, y, w, h) and the edges to the back face."""
    from matplotlib.patches import Rectangle
    dx, dy = 0.32 * w, 0.32 * w * _ASPECT
    ax.add_patch(Rectangle((x, y), w, h, fc="#ffffff", ec=INK, lw=lw, alpha=0.75, zorder=6))
    X, Y = [x, x + dx, x + w + dx, x + w + dx, x + w], [y + h, y + h + dy, y + h + dy, y + dy, y]
    ax.plot(X, Y, color=INK, lw=lw, zorder=6, solid_joinstyle="round")
    ax.plot([x + w, x + w + dx], [y + h, y + h + dy], color=INK, lw=lw, zorder=6)


def _scene(ax, run_type: str, emissions: bool = True) -> None:
    """What the box does, drawn in the lower part of the centre box: it stays at a site, or travels."""
    from matplotlib.patches import FancyArrowPatch, Polygon, Rectangle
    x0, y0, w, _h = _CENTRE
    g = y0 + 2.6                                                        # the ground
    ax.plot([x0 + 1.2, x0 + w - 1.2], [g, g], color=MUTED, lw=1.2, zorder=5)
    ax.text(x0 + w / 2, y0 + 0.55, SCENE_NOTE.get(run_type, ""), fontsize=7.8, color=MUTED, ha="center", va="bottom", style="italic", zorder=6)

    def arrow(a, b, col=MUTED, lw=1.3, ls="-"):
        ax.add_patch(FancyArrowPatch(a, b, arrowstyle="->", mutation_scale=10, lw=lw, color=col, ls=ls, shrinkA=0, shrinkB=0, zorder=6))

    if run_type == "constrained":
        xm = x0 + w / 2
        ax.plot([xm, xm], [g, g + 2.2], color=INK, lw=1.6, zorder=6)                        # the mast of the site
        ax.add_patch(Polygon([(xm - 1.1, g), (xm + 1.1, g), (xm, g + 1.3)], fc=INK, ec=INK, zorder=6))
        _cube(ax, xm - 2.6, g + 2.2, 5.2, 5.2 * _ASPECT * 0.8)
        for xa, xb in ((x0 + 3, xm - 4.2), (xm + 5.6, x0 + w - 2.5)):
            for dy in (3.6, 5.6):
                arrow((xa, g + dy), (xb, g + dy))
        ax.text(x0 + 3, g + 6.5, "wind", fontsize=7.8, color=MUTED, va="bottom", zorder=6)
        ax.text(xm + 1.6, g + 0.35, "site", fontsize=7.8, color=INK, va="bottom", zorder=6)
    elif run_type == "plume":
        xs = x0 + 3.2
        ax.add_patch(Polygon([(xs, g + 0.9), (x0 + w - 1.6, g + 2.2), (x0 + w - 1.6, g + 9.2), (xs + 4, g + 3.4)], fc="#d9d6cb", ec="none", alpha=0.9, zorder=5))
        ax.add_patch(Polygon([(xs - 0.9, g), (xs + 0.9, g), (xs, g + 2.0)], fc=COLOUR["emissions"], ec="none", zorder=6))   # the source
        ax.text(xs - 1.2, g + 2.4, "source", fontsize=7.8, color=COLOUR["emissions"], va="bottom", zorder=6)
        for xc, yc, sz in ((x0 + 8.2, g + 2.0, 2.0), (x0 + 14.2, g + 2.9, 2.9), (x0 + 21.5, g + 3.6, 3.9)):
            _cube(ax, xc, yc, sz, sz * _ASPECT * 0.85)
        arrow((x0 + 11.2, g + 1.0), (x0 + 27.5, g + 1.0))
    else:
        for xb, hb in ((6.5, 0.9), (8.0, 1.5), (12.5, 1.1), (15.5, 1.7), (17.2, 1.0), (21.0, 1.4), (23.5, 0.8)):   # a city on the way
            ax.add_patch(Rectangle((x0 + xb, g), 1.2, hb, fc="#b9b6ab", ec="none", zorder=5))
        tops = []
        for xc, hc in ((x0 + 3.4, 3.4), (x0 + 11.6, 4.8), (x0 + 19.8, 6.2)):
            _cube(ax, xc, g + 1.9, 4.2, hc)
            tops.append((xc + 2.1, g + 1.9 + hc))
            if emissions:
                arrow((xc + 2.1, g + 0.2), (xc + 2.1, g + 1.8), col=COLOUR["emissions"], lw=1.6)
        arrow((x0 + 8.4, g + 3.3), (x0 + 11.2, g + 3.3)); arrow((x0 + 16.6, g + 3.3), (x0 + 19.4, g + 3.3))
        xp = x0 + w - 3.0                                                                    # the site the air arrives at
        ax.plot([xp, xp], [g, g + 2.2], color=INK, lw=1.6, zorder=6)
        ax.add_patch(Polygon([(xp - 0.9, g), (xp + 0.9, g), (xp, g + 1.2)], fc=INK, ec=INK, zorder=6))
        ax.text(xp, g + 2.4, "site", fontsize=7.8, color=INK, ha="center", va="bottom", zorder=6)


def draw(spec: fs.FieldSpec, table: Optional[fs.FieldTable], path: Optional[str] = None, dpi: int = 150, fig=None):
    """The schematic.  `table` None draws what the kind of run needs.  With `fig` the picture goes into that
    figure (cleared first); without it a new figure is made, and saved and closed when `path` is given."""
    from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
    own = fig is None
    if own:
        import matplotlib
        if path:
            matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig = plt.figure(figsize=FIGSIZE, facecolor=SURFACE)
    else:
        fig.clear(); fig.set_facecolor(SURFACE)
    ax = fig.add_axes(list(_AXES)); ax.set_xlim(*_XLIM); ax.set_ylim(*_YLIM); ax.axis("off")
    ax.text(0.5, 65.2, f"{spec.name or 'New field case'}: {fs.RUN_TYPES.get(spec.run_type, spec.run_type)}", fontsize=15, fontweight="bold", color=INK, va="top")
    sub = " · ".join(x for x in (spec.site, spec.description) if x)
    ax.text(0.5, 62.2, textwrap.shorten(sub, 175, placeholder=" …"), fontsize=9.5, color=MUTED, va="top")

    def box(x, y, w, h, edge, lw=1.4, ls="-", face=SURFACE):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=0.9", fc=face, ec=edge, lw=lw, ls=ls, zorder=2))

    for cd in cards(spec, table):
        x, y, w, h = _SLOT[cd.key]
        col = COLOUR[cd.key] if cd.used else FAINT
        box(x, y, w, h, col, ls="-" if cd.used else (0, (4, 3)))
        ax.text(x + 1, y + h - 1.1, cd.title, fontsize=10.5, fontweight="bold", color=col if cd.used else MUTED, va="top", zorder=3)
        body = _wrap(cd.lines, WRAP + (WIDE if cd.key == "photolysis" else 0), MAX_LINES)
        ax.text(x + 1, y + h - 3.6, "\n".join(body), fontsize=8.6, color=INK if cd.used else MUTED, va="top", linespacing=1.32, zorder=3)
        if cd.used:
            s = "Source: " + (cd.source or "not stated")
            ax.text(x + 1, y + 0.9, textwrap.shorten(s, WRAP + (WIDE if cd.key == "photolysis" else 0) + 4, placeholder=" …"), fontsize=8.2,
                    style="italic", color=MUTED if cd.source else FAINT, va="bottom", zorder=3)
        a, b, style = _ARROW[cd.key]
        ax.add_patch(FancyArrowPatch(a, b, arrowstyle=style, mutation_scale=16, lw=2.2 if cd.used else 1.2, color=col,
                                     ls="-" if cd.used else (0, (3, 3)), shrinkA=1.5, shrinkB=1.5, zorder=4))

    x, y, w, h = _CENTRE
    box(x, y, w, h, INK, lw=2.2, face="#f4f3ee")
    ax.text(x + w / 2, y + h - 1.3, BOX.get(spec.run_type, "Air"), fontsize=11.5, fontweight="bold", color=INK, va="top", ha="center", zorder=3)
    ax.text(x + 1.2, y + h - 4.3, "\n".join(_wrap(centre_lines(spec, table), 57, 5)), fontsize=8.6, color=INK, va="top", linespacing=1.3, zorder=3)
    _scene(ax, spec.run_type, emissions=table is None or bool(table.emissions))

    # the equation: pieces in the colour of their card, centred as one row
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    texts = [ax.text(0, -5, p, fontsize=13, color=col, va="center") for p, col in equation(spec, table)]
    inv = ax.transData.inverted()
    widths = [inv.transform((t.get_window_extent(rend).width, 0))[0] - inv.transform((0, 0))[0] for t in texts]
    x0 = 50 - sum(widths) / 2
    for t, wd in zip(texts, widths):
        t.set_x(x0); x0 += wd
    if path:
        fig.savefig(path, dpi=dpi, facecolor=SURFACE)
        if own:
            import matplotlib.pyplot as plt
            plt.close(fig)
    return fig


def draw_case(json_path: str, out: Optional[str] = None) -> str:
    """The schematic of a saved case; raises ValueError when the case or its table cannot be read."""
    spec = fs.load_spec(json_path)
    table, report = fs.check(spec)
    if table is None:
        raise ValueError(f"{os.path.basename(json_path)}: " + ("; ".join(report.errors) or "the table could not be read"))
    out = out or os.path.join(os.path.dirname(os.path.abspath(json_path)), spec.name + ".schematic.png")
    draw(spec, table, out)
    return out


def main(argv=None) -> int:
    a = list(sys.argv[1:] if argv is None else argv)
    if a and a[0] == "--examples":
        folder = a[1] if len(a) > 1 else os.getcwd()
        os.makedirs(folder, exist_ok=True)
        for name, (p, _d) in fs.list_examples().items():
            print(draw_case(p, os.path.join(folder, name + ".schematic.png")))
        return 0
    if not a:
        print(__doc__)
        return 1
    try:
        print(draw_case(a[0], a[1] if len(a) > 1 else None))
    except ValueError as e:
        print(f"ERROR {e}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
