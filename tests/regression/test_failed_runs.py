"""What a run that goes wrong tells the user: the exact-answer experiment (tests/exact_small.py) with something
broken.  Until 25 Nov 2026 a failed run ended with exit code 0, so the windows, studies and juice.solve took the
.cs of an earlier run of the same experiment for its result; an unknown name in a reaction rate switched the
reaction off without a word; and a scenario file in the format of the early SAPHIR GUI took a route nothing tests."""
import json
import os
import subprocess
import sys
import time

import pytest

import juice
from exact_small import TOML, prepare, run

pytestmark = pytest.mark.slow


def solver(js):
    done = subprocess.run([sys.executable, "-m", "juice.core.scenario_runner", js, "--toml", TOML], capture_output=True,
                          env={**os.environ, "PYTHONIOENCODING": "utf-8", "PYTHONUTF8": "1"})
    return done.returncode, done.stdout.decode("utf-8", errors="replace")


def edit_elog(js, old, new):
    elog = json.load(open(js, encoding="utf-8"))["elog_path"]
    text = open(elog, encoding="latin-1").read()
    assert text.count(old) == 1
    open(elog, "w", encoding="latin-1").write(text.replace(old, new))


def test_a_failed_rerun_is_a_failure_and_not_the_earlier_result(tmp_path):
    js = prepare(tmp_path, "fast")
    first = juice.solve(js, toml=TOML)                                  # a good run leaves its .cs
    before = os.path.getmtime(first)
    os.remove(js[:-5] + ".enz")                                         # then the conditions file goes missing
    time.sleep(1.1)
    code, out = solver(js)
    assert code == 1 and "Simulation failed!" in out and "No .enz file found" in out
    with pytest.raises(RuntimeError, match=r"the solver stopped \(exit code 1\)"):
        juice.solve(js, toml=TOML)
    assert os.path.getmtime(first) == before                            # the earlier result is left as it was


def test_a_scenario_file_of_the_early_saphir_gui_is_refused(tmp_path):
    js = prepare(tmp_path, "fast")
    cfg = json.load(open(js, encoding="utf-8"))
    cfg.pop("elog_path")
    json.dump(cfg, open(js, "w", encoding="utf-8"))
    code, out = solver(js)
    assert code == 1 and "in the format of the early SAPHIR GUI" in out and not os.path.isfile(js[:-5] + ".cs")


@pytest.mark.parametrize("broken, name", [("CONST(1.E-4*NOSUCHNAME)", "NOSUCHNAME"), ("NOSUCHFN(1.E-4)", "NOSUCHFN")])
def test_an_unknown_name_in_a_rate_is_said(tmp_path, broken, name):
    js = prepare(tmp_path, "fast")
    edit_elog(js, "k[A --> B] = CONST(1.E-4)", f"k[A --> B] = {broken}")
    df, out = run(js, with_log=True)
    assert f"[WARNING] '{name}' in the rate of A --> B" in out and "this reaction does not run" in out
    assert float(df["A"].iloc[-1]) == pytest.approx(1e9)                # the value is still 0: A does not decay


def test_a_correct_experiment_gets_no_such_warning(tmp_path):
    _df, out = run(prepare(tmp_path, "fast"), with_log=True)
    assert "is defined nowhere" not in out and "Simulation failed" not in out


@pytest.mark.parametrize("profile", ["fast", "twostep"])
@pytest.mark.parametrize("rate, said", [
    ("CONST(1.E200)", "is 1e+200 /s, beyond any physical value for a reaction of order 1 (limit 1e+13 /s)"),
    ("CONST(1.E-4*EXP(800.))", "is inf (not a finite number) and is set to 0"),
])
def test_a_rate_beyond_physics_is_said(tmp_path, profile, rate, said):
    """Rate sanity check of 25 Nov 2026: before, a rate of 1e200 /s ran to nonsense (B = 6e210 cm-3) and an
    overflowing one was set to 0, both without a word.  The run itself is unchanged; the warning names the reaction
    and is repeated at the end of the run."""
    js = prepare(tmp_path, profile)
    edit_elog(js, "k[A --> B] = CONST(1.E-4)", f"k[A --> B] = {rate}")
    _df, out = run(js, with_log=True)
    assert f"[WARNING] the rate of A --> B" in out and said in out
    assert ">>> WARNING: 1 reaction rate(s) were not finite or beyond physical values" in out


@pytest.mark.parametrize("profile", ["fast", "precise", "twostep"])
def test_the_rates_of_a_correct_experiment_pass_the_check(tmp_path, profile):
    _df, out = run(prepare(tmp_path, profile), with_log=True)
    assert "[WARNING]" not in out and ">>> WARNING" not in out
