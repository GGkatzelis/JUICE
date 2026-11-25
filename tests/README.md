# The JUICE test suite

Install once with `pip install -e .[dev]`, then run from the repository folder.

| Command | What it runs | Time |
|---|---|---|
| `python -m pytest -n auto` | everything (about 900 tests), on all cores | about 10 min (18 min on one core) |
| `python -m pytest -n auto -m "not slow"` | everything except the long solver runs: what GitHub runs on every push | about 7 min on 4 cores |
| `python -m pytest -m unit` | the unit tests only (about 550) | about a minute |
| `python -m pytest tests/unit/test_particles.py` | one file | |
| `python -m pytest --cov --cov-report=html` | everything, and writes `htmlcov/index.html` showing every line no test reached | |

## The four layers

Each test file sits in the folder of its layer, and every test gets that folder's name as a marker.

| Folder | What belongs here | When it fails |
|---|---|---|
| `unit/` | One function or class, called with a small hand-made input, checked against a value worked out by hand. No solver run, no window. | The named function is wrong. |
| `integration/` | Several modules in a chain: experiment description → EASY/ENZ → scenario, the importers, the Python interface. | The pieces no longer fit together. |
| `regression/` | Results against a stored reference: FACSIMILE, an analytic solution, a tight solver run. | The numbers moved. |
| `gui/` | Windows opened without a screen. | A window no longer builds what it shows. |

On GitHub (`.github/workflows/tests.yml`) the quick job runs on every push, the full suite every night and on
request; a test that hangs fails after its time limit (`--timeout`) and is named in the log.

`fixtures/` holds input files for the tests; `conftest.py` holds the shared set-up. `exact_small.py` is a small
experiment whose every result has a closed form (decay, equilibrium, second order, light switch, injection,
dilution, temperature step); `regression/test_exact_scenarios.py` runs it in every solver profile, and the
`exact_small_run` fixture gives its finished run to tests of what is computed from a run (the viewer's budgets).

## Rules

1. **Every bug fix comes with a test** that fails before the fix and passes after it. Say in the test's
   docstring what went wrong and when.
2. **New code comes with unit tests** in `unit/test_<module>.py`, named after the source file it tests.
3. **File names are unique across the four folders** (pytest needs this, as the folders are not packages).
4. **A test that runs the solver for more than a few seconds** is marked `@pytest.mark.slow`.
5. **Coverage is a map, not a target.** It includes the solver, which runs as a program of its own (`patch = ["subprocess"]`
   in `pyproject.toml`). Use the report to find code no test reaches; a line being
   reached does not mean its result is checked.

## Decisions the tests record

| Date | Decision | Test |
|---|---|---|
| 24 Nov 2026 | In a study, the injection group `voc` means every organic injection: single compounds and injected solutions. `mixture` is the solutions alone (each scaled as a whole, composition unchanged), `compound` the single compounds alone. | `unit/test_study.py::test_voc_group_takes_single_compounds_and_injected_solutions` |
| 25 Nov 2026 | A scenario file without an expanded EASY file (the format of the early SAPHIR GUI) is refused with a message instead of run on the untested route; a failed run ends with exit code 1. | `regression/test_failed_runs.py` |
| 24 Nov 2026 | An injection timed before the run still enters in the first step and one after the end does not enter (the ENZ files stay as they were); the review of the experiment warns about both. | `unit/test_spec_grid_and_archive.py::test_pulse_timed_before_the_run_enters_in_the_first_step_and_the_review_says_so` |
