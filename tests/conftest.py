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


@pytest.fixture(scope="session")
def exact_small_run(tmp_path_factory):
    """The exact-answer experiment (tests/exact_small.py) run once with the fast profile:
    (scenario file, result table).  For tests of what is computed from a finished run."""
    import exact_small
    js = exact_small.prepare(tmp_path_factory.mktemp("exact_small"), "fast")
    return js, exact_small.run(js)


@pytest.fixture(scope="session")
def small_study(tmp_path_factory):
    """A real study, run: the UV-lamp chamber template with ethanol injected at 10:30 (four hours, lamps on
    09:00-13:00), swept over the initial NO2 (5 and 50 ppb) and the injected ethanol (x 0.5 and x 2), so a 2 x 2
    isopleth.  Four runs of a few seconds each.  Returns the path of the .study.json."""
    import contextlib
    import io
    import juice
    from juice.core import study as S
    from juice.core.experiment_spec import InjectionEvent
    root = tmp_path_factory.mktemp("small_study")
    spec = juice.new_experiment("uv_ethanol", template="uv_chamber_lamps")
    c = spec.conditions
    c.start, c.end = "08:00", "12:00"
    spec.precursors = ["C2H5OH"]
    spec.initials = list(spec.initials) + [("NO2", "M*1e-9*20")]
    spec.injections = [InjectionEvent("ethanol", "C2H5OH", "200", "10:30", 2.0, action="voc")]
    exp_dir = root / "uv_ethanol"
    with contextlib.redirect_stdout(io.StringIO()):
        files = juice.write(spec, str(exp_dir))
    toml = os.path.join(REPO, "settings.toml")
    st = S.new_study("iso", str(root / "studies" / "iso"), toml)
    st.experiments = [files["spec_path"]]
    st.mode = "sweep"
    st.workers = 2
    st.axes = [S.Axis("initial", "NO2", [5.0, 50.0]), S.Axis("injection_scale", "ethanol", [0.5, 2.0])]
    os.makedirs(st.folder, exist_ok=True)
    path = S.save_study(st, os.path.join(st.folder, "iso.study.json"))
    runner = S.StudyRunner(st)
    with contextlib.redirect_stdout(io.StringIO()):
        errors = runner.materialise_all()
        assert not errors, errors
        runner.run_all(workers=2)
    failed = [(p.name, p.message) for p in runner.plans if not S.is_done(p)]
    assert not failed, failed
    return path
