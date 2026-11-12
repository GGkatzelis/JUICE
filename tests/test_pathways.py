"""Pathway tracing on a network small enough to check by hand, and on the SAPHIR ethanol run."""
import os

import numpy as np
import pytest

from juice.chemistry import pathways as pw

T = np.linspace(0.0, 10.0, 101)               # hours


def _rate(v, t0=0.0):
    return np.where(T >= t0, float(v), 0.0)


# injected A (via holding species AIN) -> B (+OH, 6 ppb/h) and C (+O3, 4 ppb/h); B + NO -> T, C + hv -> T;
# T <-> X (an equilibrium: 50 forward, 49 back -> net 1); AIN carries no carbon of its own
NET = {
    "A": [pw.Producer([(1.0, "AIN")], 1.0, _rate(10))],
    "B": [pw.Producer([(1.0, "A"), (1.0, "OH")], 1.0, _rate(6))],
    "C": [pw.Producer([(1.0, "A"), (1.0, "O3")], 1.0, _rate(4))],
    "T": [pw.Producer([(1.0, "B"), (1.0, "NO")], 1.0, _rate(6)),
          pw.Producer([(1.0, "C")], 1.0, _rate(4), photolysis=True, rate_expr="J41"),
          pw.Producer([(1.0, "X")], 1.0, _rate(49))],
    "X": [pw.Producer([(1.0, "T"), (1.0, "NO2")], 1.0, _rate(50))],
}
FORMULA = {"A": "C3H6", "B": "C3H6O", "C": "C2H4O", "T": "C2H4O2", "X": "C2H3NO5"}


def _tracer(window=None):
    return pw.Tracer(lambda sp: NET.get(sp, []), T, lambda n: FORMULA.get(n, ""), window)


def test_shares_multiply_along_paths():
    g = _tracer().trace("T", threshold=0.01)
    e = {(x.parent, x.child, x.partner): round(x.share, 4) for x in g.edges}
    assert e[("B", "T", "NO")] == 0.6 and e[("C", "T", "hv")] == 0.4      # the equilibrium with X is net zero here
    assert e[("A", "B", "OH")] == 0.6 and e[("A", "C", "O3")] == 0.4
    assert e[("source:AIN", "A", "")] == 1.0
    assert g.node_share()["A"] == pytest.approx(1.0)
    assert g.produced_ppb == pytest.approx(10 * 10.0, rel=1e-6)            # 6 + 4 ppb/h for 10 h


def test_equilibrium_counts_net():
    NET["T"][2] = pw.Producer([(1.0, "X")], 1.0, _rate(60))                # X now gives back more than it takes: +10 net
    try:
        g = _tracer().trace("T", threshold=0.01)
        e = {(x.parent, x.child): round(x.share, 3) for x in g.edges}
        assert e[("X", "T")] == 0.5 and e[("B", "T")] == 0.3
    finally:
        NET["T"][2] = pw.Producer([(1.0, "X")], 1.0, _rate(49))


def test_threshold_and_window():
    g = _tracer().trace("T", threshold=0.45)
    assert {(x.parent, x.child) for x in g.edges} == {("B", "T"), ("A", "B"), ("source:AIN", "A")}
    assert g.other["T"] == pytest.approx(0.4)
    NET["C"] = [pw.Producer([(1.0, "A"), (1.0, "O3")], 1.0, _rate(4, t0=5.0))]   # C only forms after 5 h
    try:
        early = _tracer((0.0, 4.0)).trace("C")
        assert not early.edges and early.produced_ppb == 0.0
    finally:
        NET["C"] = [pw.Producer([(1.0, "A"), (1.0, "O3")], 1.0, _rate(4))]


def test_draw_and_size():
    from matplotlib.figure import Figure
    g = _tracer().trace("T", threshold=0.01)
    fig = Figure()
    pos = pw.draw(fig, g, lambda n: FORMULA.get(n, ""), lambda n: {"A": "C=CC", "T": "OC(=O)C"}.get(n, ""))
    assert set(pos) == {"T", "B", "C", "A", "source:AIN"}
    assert pos["T"][0] > pos["B"][0] > pos["A"][0]
    w, h = pw.figure_size(pos)
    assert w > 6 and h > 3


def test_structure_image():
    img = pw.structure_image("CCO")
    assert img is not None and img.shape[2] == 4
    assert pw.structure_image("not a smiles") is None


CS = os.path.join(os.path.dirname(__file__), "..", "runs", "benchmark", "fast", "scene_05", "scene_05.cs")


@pytest.mark.skipif(not os.path.isfile(CS), reason="needs the ethanol benchmark run (runs/benchmark/fast/scene_05)")
def test_ethanol_run():
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PyQt6.QtWidgets import QApplication
    app = QApplication.instance() or QApplication([])                 # noqa: F841
    from juice.gui import results_viewer as rv
    v = rv.ResultsViewer(rv.load_column_file(CS), csv_path=CS)
    producers, t = v.budget_engine.pathway_source(v.df)
    formula, _smiles = v._species_info()
    g = pw.Tracer(producers, t, formula).trace("PAN", 0.05)
    e = {(x.parent, x.child, x.partner): x.share for x in g.edges}
    assert e[("CH3CO3", "PAN", "NO2")] == pytest.approx(1.0)
    assert e[("CH3CHO", "CH3CO3", "OH")] > 0.8 and e[("C2H5OH", "CH3CHO", "OH")] > 0.9
    h = pw.Tracer(producers, t, formula).trace("HCHO", 0.05)
    parents = {x.parent for x in h.edges}
    assert "CH3O2NO2" not in parents                                  # the CH3O2 <-> CH3O2NO2 equilibrium is net
    assert {"CH3O", "HOCH2CH2O", "CH3CO3"} <= parents
    # a double-click on a species in the list on the left draws it when the Pathways tab is open
    from PyQt6.QtWidgets import QTreeWidgetItemIterator
    v.right_tabs.setCurrentIndex(v._tab("Pathways"))
    it = QTreeWidgetItemIterator(v.var_tree)
    while it.value() and it.value().text(0) != "CH3CHO":
        it += 1
    assert it.value() is not None
    v._on_tree_double_click(it.value(), 0)
    assert v._pw_target == "CH3CHO" and v._pw_graph.target == "CH3CHO"
    assert ("C2H5OH", "CH3CHO") in {(x.parent, x.child) for x in v._pw_graph.edges}
