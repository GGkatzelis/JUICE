"""Shared set-up for the JUICE test suite: the repository on sys.path, the reference scenes."""
import os
import sys

import pytest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO not in sys.path:
    sys.path.insert(0, REPO)

IDL_INPUT = os.path.join(REPO, "validation", "IDL Input")
SCENES = sorted(d for d in os.listdir(IDL_INPUT)
                if d.startswith("scene_") and os.path.isdir(os.path.join(IDL_INPUT, d)))


LAYERS = {
    "unit": "one function or class with a small hand-made input; no solver run, no window",
    "integration": "several modules in a chain (description -> EASY/ENZ -> scenario, importers, the API)",
    "regression": "results against a stored reference (FACSIMILE, an analytic solution, a tight solver run)",
    "gui": "windows opened without a screen",
}


def pytest_configure(config):
    config.addinivalue_line("markers", "slow: runs the solver (seconds to a minute); skip with -m 'not slow'")
    for name, what in LAYERS.items():
        config.addinivalue_line("markers", f"{name}: {what} (set from the folder, tests/{name}/)")


def pytest_collection_modifyitems(items):
    """Every test carries the marker of its folder, so `pytest -m unit` selects a layer."""
    here = os.path.dirname(os.path.abspath(__file__))
    for item in items:
        layer = os.path.relpath(str(item.fspath), here).split(os.sep)[0]
        if layer in LAYERS:
            item.add_marker(getattr(pytest.mark, layer))


@pytest.fixture(autouse=True)
def _repo_cwd(monkeypatch):
    """Several tools resolve paths relative to the repository; run every test from there."""
    monkeypatch.chdir(REPO)
