"""Where the JUICE folder is: the one place that knows it.

The data (chemistry/, setups/, validation/, examples/), the settings file settings.toml and the logo sit next
to the 'juice' package, in ROOT.
Every module takes ROOT from here instead of working it out from its own location.
"""
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def module_args(name: str) -> list:
    """The arguments that start a JUICE module in its own process: [python] + module_args('juice.core.scenario_runner')."""
    return ["-m", name]


def export_pythonpath() -> None:
    """Processes JUICE starts ('python -m juice....') must find the package from any working directory."""
    parts = [p for p in os.environ.get("PYTHONPATH", "").split(os.pathsep) if p]
    if ROOT not in parts:
        os.environ["PYTHONPATH"] = os.pathsep.join([ROOT] + parts)
    if ROOT not in sys.path:
        sys.path.insert(0, ROOT)
