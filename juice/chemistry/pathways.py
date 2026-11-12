"""Pathways to a product: which chains of reactions made a species during a run, and how much each contributed.

Backward flux attribution. Over a time window, every reaction that forms the target is weighted by what it produced
(its rate integrated over the window, times the product's yield). Its share of the target's production goes to its
*parent*: the carbon-carrying reactant (in RO2 + NO -> RO the parent is RO2; NO is the partner, shown on the arrow).
The parent's own production is split the same way, and the shares multiply along a path: a step that makes 60 % of
an intermediate which makes 50 % of the target contributes 30 %. Paths whose share falls below a threshold are
summed per species as "other routes".

Recycling (an equilibrium such as CH3CO3 + NO2 <-> PAN) is left out: a producer whose parent already lies on the
path towards the target is dropped and the rest renormalised, so the shares describe fresh production. A reaction
without a carbon-carrying reactant (an injection, a chamber source, a holding species of a delayed release) ends a
path in a source box.

The shares are an attribution: exact for chains, the usual proportional approximation where branches meet.

Inputs are plain: a RateSource gives, for a species, the reactions producing it with their rates (ppb/h) on the
run's time axis (results_viewer.BudgetEngine.pathway_source builds one from a run); a species lookup gives formula
and SMILES (species_db.SpeciesDB). Drawing uses matplotlib, and RDKit for the structures when it is installed.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np

_trapz = getattr(np, "trapezoid", None) or getattr(np, "trapz")
NOT_PARENTS = {"CO", "CO2"}
PARTNER_NAMES = {"HO": "OH", "HV": "hv", "O3P": "O"}


@dataclass
class Producer:
    reactants: List[Tuple[float, str]]       # (coefficient, name)
    yield_: float                            # stoichiometric coefficient of the species among the products
    rate_ppb_h: np.ndarray                   # reaction rate on the run's time axis
    label: str = ""
    photolysis: bool = False
    rate_expr: str = ""


@dataclass
class Edge:
    parent: str
    child: str
    share: float                             # fraction of the target's production that flows along this step
    partner: str                             # OH, O3, NO3, NO, HO2, hv, RO2, '' (unimolecular)


@dataclass
class PathwayGraph:
    target: str
    produced_ppb: float                      # target produced in the window (sum over its producers)
    window: Tuple[float, float]              # hours
    edges: List[Edge] = field(default_factory=list)
    other: Dict[str, float] = field(default_factory=dict)       # species -> share left in small routes
    sources: Dict[str, str] = field(default_factory=dict)       # source node -> what it is

    def node_share(self) -> Dict[str, float]:
        """Share of the target's production that passes through each species (the target: 1)."""
        out: Dict[str, float] = {self.target: 1.0}
        for e in self.edges:
            out[e.parent] = out.get(e.parent, 0.0) + e.share
        return out


class Tracer:
    def __init__(self, producers: Callable[[str], List[Producer]], t_hours: np.ndarray,
                 formula: Callable[[str], str], window: Optional[Tuple[float, float]] = None):
        self._producers = producers
        self.t = np.asarray(t_hours, float)
        lo, hi = window if window else (float(self.t[0]), float(self.t[-1]))
        self.window = (lo, hi)
        self._mask = (self.t >= lo) & (self.t <= hi)
        self._formula = formula
        self._cache: Dict[str, List[Tuple[str, str, float]]] = {}
        self._gross_cache: Dict[str, Dict[Tuple[str, str], float]] = {}

    # ── chemistry rules ──
    def is_organic(self, name: str) -> bool:
        if name in NOT_PARENTS:
            return False
        f = self._formula(name) or ""
        return bool(re.search(r"C(?![a-z])", f))

    def _parent_and_partner(self, p: Producer) -> Tuple[str, str]:
        names = [n for _c, n in p.reactants if n not in ("HV",)]
        organic = [n for n in names if self.is_organic(n)]
        if not organic:
            return "", ""
        parent = max(organic, key=lambda n: self._carbon(n))
        others = [PARTNER_NAMES.get(n, n) for n in names if n != parent]
        if p.photolysis or any(n == "HV" for _c, n in p.reactants):
            partner = "hv"
        elif others:
            partner = others[0]
        elif re.search(r"\bRO2\b", p.rate_expr or "", re.I):
            partner = "RO2"
        else:
            partner = ""
        return parent, partner

    def _carbon(self, name: str) -> int:
        m = re.search(r"C(\d*)(?![a-z])", self._formula(name) or "")
        return int(m.group(1) or 1) if m else 0

    def _integral(self, y: np.ndarray) -> float:
        t, y = self.t[self._mask], np.asarray(y, float)[self._mask]
        if len(t) < 2:
            return 0.0
        return float(_trapz(np.maximum(y, 0.0), t))

    def _gross(self, species: str) -> Dict[Tuple[str, str], float]:
        """{(parent or 'source:<what>', partner): ppb of `species` formed in the window} — gross, both directions."""
        if species not in self._gross_cache:
            out: Dict[Tuple[str, str], float] = {}
            for p in self._producers(species):
                amount = p.yield_ * self._integral(p.rate_ppb_h)
                if amount <= 0:
                    continue
                parent, partner = self._parent_and_partner(p)
                if not parent:
                    what = " + ".join(n for _c, n in p.reactants) or "emission / injection"
                    parent, partner = f"source:{what}", ""
                out[(parent, partner)] = out.get((parent, partner), 0.0) + amount
            self._gross_cache[species] = out
        return self._gross_cache[species]

    def contributions(self, species: str) -> List[Tuple[str, str, float]]:
        """[(parent or 'source:<what>', partner, ppb formed in the window)] for the reactions forming `species`, net of
        the way back: an equilibrium (CH3O2 + NO2 <-> CH3O2NO2) counts only with what it gave more than it took."""
        if species not in self._cache:
            out = []
            for (parent, partner), amount in self._gross(species).items():
                if not parent.startswith("source:"):
                    back = sum(v for (pp, _q), v in self._gross(parent).items() if pp == species)
                    amount -= back
                if amount > 0:
                    out.append((parent, partner, amount))
            self._cache[species] = out
        return self._cache[species]

    # ── the trace ──
    def trace(self, target: str, threshold: float = 0.05, max_depth: int = 8) -> PathwayGraph:
        target = target.upper()
        g = PathwayGraph(target, sum(v for _p, _q, v in self.contributions(target)), self.window)
        edges: Dict[Tuple[str, str, str], float] = {}

        def expand(sp: str, share: float, depth: int, path: frozenset):
            contrib = [(p, q, v) for p, q, v in self.contributions(sp) if p not in path and p != sp]    # no recycling
            total = sum(v for _p, _q, v in contrib)
            if total <= 0:
                return
            small = 0.0
            for parent, partner, v in sorted(contrib, key=lambda x: -x[2]):
                s = share * v / total
                if s < threshold:
                    small += s; continue
                edges[(parent, sp, partner)] = edges.get((parent, sp, partner), 0.0) + s
                if parent.startswith("source:"):
                    g.sources[parent] = parent[len("source:"):]
                elif depth < max_depth:
                    expand(parent, s, depth + 1, path | {sp})
            if small > 0:
                g.other[sp] = g.other.get(sp, 0.0) + small

        expand(target, 1.0, 1, frozenset())
        g.edges = [Edge(p, c, s, q) for (p, c, q), s in sorted(edges.items(), key=lambda kv: -kv[1])]
        return g


# ─────────────────────────────────────────────────────────────────────────────
# Drawing
# ─────────────────────────────────────────────────────────────────────────────

def layout(g: PathwayGraph) -> Dict[str, Tuple[float, float]]:
    """x = -steps from the target (longest path), y spread within each column, near the species it feeds."""
    depth = {g.target: 0}
    for _ in range(len(g.edges) + 1):                 # longest-path relaxation (the graph has no cycles by construction)
        changed = False
        for e in g.edges:
            d = depth.get(e.child, 0) + 1
            if e.child in depth and depth.get(e.parent, -1) < d:
                depth[e.parent] = d; changed = True
        if not changed:
            break
    share = g.node_share()
    pos: Dict[str, Tuple[float, float]] = {g.target: (0.0, 0.0)}
    for d in range(1, max(depth.values()) + 1):
        col = [n for n, k in depth.items() if k == d]
        want = {}
        for n in col:                                 # the share-weighted height of the species it feeds
            out = [(e.share, pos[e.child][1]) for e in g.edges if e.parent == n and e.child in pos]
            w = sum(s for s, _y in out)
            want[n] = sum(s * y for s, y in out) / w if w else 0.0
        col.sort(key=lambda n: (want[n], -share.get(n, 0)))
        ys = [want[n] for n in col]
        for i in range(1, len(ys)):                   # keep at least one unit apart
            ys[i] = max(ys[i], ys[i - 1] + 1.0)
        shift = (sum(want.values()) / len(col)) - (sum(ys) / len(ys)) if col else 0.0
        for n, y in zip(col, ys):
            pos[n] = (-float(d), y + shift)
    return pos


_IMG_CACHE: Dict[str, Optional[np.ndarray]] = {}


def structure_image(smiles: str, size: Tuple[int, int] = (350, 240)):
    """An RGBA array of the molecule (RDKit), or None without RDKit or for a SMILES it cannot read."""
    if not smiles:
        return None
    if smiles in _IMG_CACHE:
        return _IMG_CACHE[smiles]
    img = None
    try:
        import io
        from PIL import Image
        from rdkit import Chem, RDLogger
        from rdkit.Chem.Draw import rdMolDraw2D
        RDLogger.DisableLog("rdApp.*")
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            d = rdMolDraw2D.MolDraw2DCairo(*size)
            o = d.drawOptions()
            o.minFontSize = 26; o.bondLineWidth = 3; o.padding = 0.1
            d.DrawMolecule(mol); d.FinishDrawing()
            img = np.asarray(Image.open(io.BytesIO(d.GetDrawingText())).convert("RGBA"))
    except Exception:
        img = None
    _IMG_CACHE[smiles] = img
    return img


DX, DY = 2.7, 1.95          # column and row spacing (data units)
BW, BH = 1.75, 1.2          # structure box
PX_PER_UNIT = 85            # figure size per data unit, so a structure box is ~150 px wide whatever the chain length


def figure_size(pos: Dict[str, Tuple[float, float]], dpi: float = 100.0) -> Tuple[float, float]:
    """(width, height) in inches that draws the graph at PX_PER_UNIT (plus the colour bar and title)."""
    if not pos:
        return 8.0, 4.0
    xs = [p[0] for p in pos.values()]; ys = [p[1] for p in pos.values()]
    w = (max(xs) - min(xs) + 2.3 * BW) * PX_PER_UNIT / dpi + 1.4
    h = (max(ys) - min(ys) + 1.75 * BH + 0.5) * PX_PER_UNIT / dpi + 1.0
    return max(w, 6.0), max(h, 3.5)


def draw(fig, g: PathwayGraph, formula: Callable[[str], str] = lambda n: "", smiles: Callable[[str], str] = lambda n: "",
         structures: bool = True, cmap: str = "plasma_r"):
    """Draw the pathway graph into a matplotlib Figure; returns {node: (x, y)} (data coordinates) for click handling.
    Each species is a box with its structure (RDKit) or its name; arrows run from parent to product, their width and
    colour scaled with the share of the target formed along them, labelled with the partner and the share."""
    import matplotlib as mpl
    from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
    fig.clear()
    ax = fig.add_subplot(111)
    ax.set_axis_off()
    if not g.edges:
        ax.text(0.5, 0.5, f"No net formation of {g.target} found in the window "
                f"{g.window[0]:.2f}–{g.window[1]:.2f} h", ha="center", va="center", transform=ax.transAxes)
        return {}
    unit = layout(g)
    pos = {n: (x * DX, y * DY) for n, (x, y) in unit.items()}
    share = g.node_share()
    colours = mpl.colormaps[cmap]
    norm = mpl.colors.Normalize(vmin=0.0, vmax=1.0)
    xs = [p[0] for p in pos.values()]; ys = [p[1] for p in pos.values()]
    ax.set_xlim(min(xs) - BW, max(xs) + BW * 1.3); ax.set_ylim(min(ys) - BH - 0.5, max(ys) + BH * 0.75)
    pairs: Dict[Tuple[str, str], List[Edge]] = {}
    for e in g.edges:
        pairs.setdefault((e.parent, e.child), []).append(e)
    for (pa, ch), group in pairs.items():             # arrows from box edge to box edge
        a, b = pos[pa], pos[ch]
        skip = round((b[0] - a[0]) / DX) > 1
        base = 0.35 if skip else (0.0 if abs(a[1] - b[1]) < 1e-6 else 0.12)
        start, end = (a[0] + BW / 2, a[1]), (b[0] - BW / 2, b[1])
        if skip and abs(a[1] - b[1]) < 1e-6:          # a step over a column: leave from the top, arrive at the top
            start, end = (a[0], a[1] + BH / 2), (b[0], b[1] + BH / 2)
            base = -0.3
        m = len(group)
        for i, e in enumerate(sorted(group, key=lambda x: -x.share)):
            rad = base + (i - (m - 1) / 2) * 0.45      # several steps between the same pair (+OH, +NO3): fanned out
            ax.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=10 + 14 * e.share,
                                         linewidth=0.8 + 9 * e.share, color=colours(norm(e.share)), alpha=0.9,
                                         connectionstyle=f"arc3,rad={rad}", zorder=1, shrinkA=0, shrinkB=1))
            dx, dy = end[0] - start[0], end[1] - start[1]
            mx = (start[0] + end[0]) / 2 + rad * dy / 2      # the middle of the arc3 curve
            my = (start[1] + end[1]) / 2 - rad * dx / 2
            txt = ("hν " if e.partner == "hv" else (f"+{e.partner} " if e.partner else "")) + f"{100 * e.share:.0f} %"
            ax.text(mx, my + (0.06 if rad <= 0 else -0.06), txt, ha="center", va="bottom" if rad <= 0 else "top",
                    fontsize=8, zorder=3, bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.85))
    for n, (x, y) in pos.items():
        is_source = n.startswith("source:")
        edge = colours(norm(share.get(n, 0)))
        img = structure_image(smiles(n)) if (structures and not is_source) else None
        if is_source:
            label = "source: " + g.sources.get(n, n[len("source:"):])
            ax.text(x, y, label, ha="center", va="center", fontsize=8, style="italic", zorder=2,
                    bbox=dict(boxstyle="round,pad=0.35", fc="white", ec=edge, lw=1.3, ls="--"))
            continue
        ax.add_patch(FancyBboxPatch((x - BW / 2, y - BH / 2), BW, BH, boxstyle="round,pad=0.02,rounding_size=0.08",
                                    fc="white", ec=edge, lw=2.2 if n == g.target else 1.5, zorder=2))
        if img is not None:
            ax.imshow(img, extent=(x - BW / 2 + 0.05, x + BW / 2 - 0.05, y - BH / 2 + 0.05, y + BH / 2 - 0.05), zorder=2.5,
                      interpolation="antialiased")
        else:
            ax.text(x, y, n, ha="center", va="center", fontsize=10, zorder=3)
        f = formula(n)
        note = n + (f"\n{f}" if f else "")
        note += f" · {100 * share.get(n, 0):.0f} % of {g.target}" if n != g.target else ""
        if g.other.get(n, 0) >= 0.005:
            note += f"\nother routes {100 * g.other[n]:.0f} %"
        ax.text(x, y - BH / 2 - 0.05, note, ha="center", va="top", fontsize=7.5, zorder=3,
                fontweight="bold" if n == g.target else "normal")
    ax.set_aspect("equal", adjustable="box")      # boxes keep their shape; the figure fits around them
    fig.subplots_adjust(left=0.01, right=0.94, top=0.9, bottom=0.02)
    sm = mpl.cm.ScalarMappable(norm=norm, cmap=colours); sm.set_array([])
    cb = fig.colorbar(sm, ax=ax, fraction=0.022, pad=0.01)
    cb.set_label(f"share of the {g.target} formed")
    cb.ax.yaxis.set_major_formatter(mpl.ticker.PercentFormatter(1.0))
    ax.set_title(f"Pathways to {g.target}: {g.produced_ppb:.3g} ppb formed (net), {g.window[0]:.2f}–{g.window[1]:.2f} h. "
                 f"Arrow width and colour: share of the {g.target} formed along that step", fontsize=10)
    return pos
