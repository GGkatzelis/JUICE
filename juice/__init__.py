"""
JUICE from code
===============

    import juice

    exp = juice.new_experiment("isoprene_day", chamber="SAPHIR")
    ...
    run = juice.run(exp, "runs/isoprene_day")
    o3 = run.results().ppb("O3")

The functions (juice.new_experiment, juice.run, juice.results, ...) live in juice/api.py and are loaded on first
use, so that starting one part of JUICE (the solver, the results viewer) does not load all the others.

The package:

    juice/api.py       JUICE from code: everything the windows do, as plain calls
    juice/core/        experiment and field case as data, ENZ and EASY writers, parser, solver, runner, studies
    juice/chemistry/   mechanism registry, species databases, MCM converter, lamp photolysis
    juice/importers/   F0AM, GECKO-A, FACSIMILE and SAPHIR-wiki import; F0AM export
    juice/gui/         the windows: start menu, experiment, field case, study, results viewer
    juice/tools/       checks against the reference scenarios, the FACSIMILE benchmark, the desktop shortcut

The data (chemistry/, setups/, validation/) sits next to the package: juice/_root.py.
"""
from juice._root import ROOT, export_pythonpath as _export_pythonpath

_export_pythonpath()          # processes JUICE starts with "python -m juice...." find the package


def __getattr__(name):
    if name.startswith("__"):
        raise AttributeError(name)
    import importlib
    api = importlib.import_module("juice.api")
    try:
        return getattr(api, name)
    except AttributeError:
        raise AttributeError(f"module 'juice' has no attribute {name!r}") from None


def __dir__():
    import importlib
    return sorted(set(globals()) | set(n for n in dir(importlib.import_module("juice.api")) if not n.startswith("_")))
