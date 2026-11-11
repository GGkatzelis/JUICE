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


def pytest_configure(config):
    config.addinivalue_line("markers", "slow: runs the solver (seconds to a minute); skip with -m 'not slow'")


@pytest.fixture(autouse=True)
def _repo_cwd(monkeypatch):
    """Several tools resolve paths relative to the repository; run every test from there."""
    monkeypatch.chdir(REPO)
