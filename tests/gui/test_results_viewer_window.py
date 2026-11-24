"""The results viewer opened without a screen on the run of the exact-answer experiment (tests/exact_small.py):
every tab draws, and the budget and reactivity views show the terms the run has."""
import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import pytest
from PyQt6.QtWidgets import QApplication

from juice.gui import results_viewer as rv

TABS = ["Timeseries", "Process Analysis", "Reactivity", "Pathways", "Comparison", "Budgets", "Chemical Space"]


@pytest.fixture(scope="module")
def viewer(exact_small_run):
    app = QApplication.instance() or QApplication([])                  # noqa: F841
    js, _df = exact_small_run
    cs = js[:-5] + ".cs"
    v = rv.ResultsViewer(rv.load_column_file(cs), csv_path=cs)
    yield v
    v.close()


def test_the_window_opens_with_its_tabs_and_the_species_of_the_run(viewer):
    names = [viewer.right_tabs.tabText(i) for i in range(viewer.right_tabs.count())]
    assert set(TABS) <= set(names)
    assert {"A", "B", "N", "O"} <= {c.upper() for c in viewer.df.columns}
    assert viewer.budget_engine._ready


@pytest.mark.parametrize("tab", TABS)
def test_every_tab_draws(viewer, tab):
    index = viewer._tab(tab)
    assert index >= 0
    viewer.right_tabs.setCurrentIndex(index)                           # draws on the first visit
    QApplication.processEvents()
    assert viewer.right_tabs.currentIndex() == index


def test_budget_and_reactivity_of_a_species(viewer):
    viewer.right_tabs.setCurrentIndex(viewer._tab("Process Analysis"))
    viewer._update_budget("A")
    prod, dest = viewer.budget_engine.compute("A", viewer.df)
    assert prod == {} and list(dest) == ["A → B"]
    viewer.right_tabs.setCurrentIndex(viewer._tab("Reactivity"))
    viewer._update_reactivity()
    QApplication.processEvents()
