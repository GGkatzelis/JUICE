"""The small experiment with closed-form answers that several tests run (see regression/test_exact_scenarios.py for
what it holds): its EASY text, its ENZ table, the exact solution, and how to prepare and run it."""
import contextlib
import io
import math
import os
import subprocess
import sys

from juice._root import ROOT
from juice.core import scenario_prep as sp
from juice.core.column_file_reader import load_column_file

TOML = os.path.join(ROOT, "settings.toml")
END = 7200
M = 1013.0 / (300.0 * 1.379e-19)                 # air density of the EASY file's DENSITY at the start

EASY = """; exact_small - every result has a closed form (tests/regression/test_exact_scenarios.py)
GEOMETRY[NBOX]=1
CONST=FUNCTION[%1]
DENSITY=FUNCTION[%1/(%2*1.379E-19)]

A --> B
k[A --> B] = CONST(1.E-4)

C --> D
k[C --> D] = CONST(2.E-4)
D --> C
k[D --> C] = CONST(1.E-4)

E + F --> G
k[E + F --> G] = CONST(1.E-13)

H --> I
k[H --> I] = CONST(HV*5.E-4)

QJ --> J
QJ = INPUT(QJ)
k[QJ --> J] = CONST(M*1e-9/60*8.)

L + DIL -->
k[L + DIL --> ] = CONST(DILUTE)
DILUTE = CONST(FLOW/(276.*3600.))
DIL    = CONST(1.)

N --> O
k[N --> O] = CONST(1.E-2*EXP(-1000./T))

INITIAL[A] = CONST(1.E9)
INITIAL[C] = CONST(3.E9)
INITIAL[E] = CONST(2.E9)
INITIAL[F] = CONST(1.E9)
INITIAL[H] = CONST(1.E9)
INITIAL[L] = CONST(1.E9)
INITIAL[N] = CONST(1.E9)

HV   = INPUT(HV)
FLOW = INPUT(FLOW)
T    = INPUT(T)
P    = INPUT(P)
M    = DENSITY(P,T)
H2O  = CONST(0.)

FILES[ENZ]    = exact_small.enz
FACS[HMAX]    = 30
FACS[OUTSTEP] = 600
PLOT[SKIP] = CS,KV,PD
"""


def enz_text():
    lines = ["COLUMN 1=time", "COLUMN 2=HV", "COLUMN 3=QJ", "COLUMN 4=FLOW", "COLUMN 5=T", "COLUMN 6=P",
             "NUMBER OF COLUMNS=6", "&&&&&&&&&&&&&&&&&&"]
    for t in range(0, END + 1, 60):
        hv = 1.0 if 1800 <= t < 5400 else 0.0
        qj = 0.5 if t in (600, 660) else 0.0
        flow = 27.6 if t >= 3600 else 0.0
        temp, p = (250.0, 1013.0 * 250.0 / 300.0) if t < 3600 else (500.0, 1013.0 * 500.0 / 300.0)
        lines.append(f"{t:.1f},  {hv:.1f},  {qj:.1f},  {flow:.1f},  {temp:.1f},  {p:.6f}")
    return "\n".join(lines) + "\n"


def exact(t):
    """Every species at time t (seconds)."""
    lit = min(max(t - 1800, 0), 3600)
    injected = min(max(t - 600, 0), 120) / 120.0
    diluted = max(t - 3600, 0)
    k_cold, k_hot = 1e-2 * math.exp(-4.0), 1e-2 * math.exp(-2.0)
    n = 1e9 * math.exp(-k_cold * min(t, 3600) - k_hot * max(t - 3600, 0))
    c_eq = 3e9 * 1e-4 / 3e-4
    c = c_eq + (3e9 - c_eq) * math.exp(-3e-4 * t)
    f = 1e9 * 1e9 / (2e9 * math.exp(1e-13 * 1e9 * t) - 1e9)             # F0 (E0 - F0) / (E0 e^(k (E0 - F0) t) - F0)
    a = 1e9 * math.exp(-1e-4 * t)
    h = 1e9 * math.exp(-5e-4 * lit)
    return {"A": a, "B": 1e9 - a, "C": c, "D": 3e9 - c, "E": f + 1e9, "F": f, "G": 1e9 - f, "H": h, "I": 1e9 - h,
            "J": 8e-9 * M * injected, "L": 1e9 * math.exp(-27.6 / (276.0 * 3600.0) * diluted), "N": n, "O": 1e9 - n}


def prepare(folder, profile):
    """Write the two input files into `folder` and prepare the scenario; returns the scenario file."""
    folder = str(folder)
    with open(os.path.join(folder, "exact_small.easy"), "w", encoding="latin-1") as fh:
        fh.write(EASY)
    with open(os.path.join(folder, "exact_small.enz"), "w", encoding="latin-1") as fh:
        fh.write(enz_text())
    with contextlib.redirect_stdout(io.StringIO()):
        return sp.prepare_scenario({"easy": os.path.join(folder, "exact_small.easy"), "enz": os.path.join(folder, "exact_small.enz")},
                                   sp.load_paths(TOML)["include_dir"], export_kv_pd=False, toml_path=TOML,
                                   solver=sp.solver_from_toml(TOML, profile))


def run(js):
    """Run a prepared scenario; returns the result table with upper-case column names."""
    subprocess.run([sys.executable, "-m", "juice.core.scenario_runner", js, "--toml", TOML], capture_output=True, check=True)
    df = load_column_file(js[:-5] + ".cs")
    df.columns = [c.upper() for c in df.columns]
    return df
