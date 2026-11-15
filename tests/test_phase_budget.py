"""Organic carbon and oxygen by phase (gas, particles, wall) for a run with partitioning; the Phases tab."""
import os

import numpy as np
import pytest

import juice
from juice.chemistry import phase_budget as pb
from juice.chemistry.species_db import load_species_db
from juice.core.column_file_reader import load_column_file


@pytest.fixture(scope="module")
def run(tmp_path_factory):
    e = juice.new_experiment("phases", chamber="SAPHIR")
    c = e.conditions
    c.source, c.date, c.start, c.end = "manual", "2024-07-09", "08:00", "10:00"
    c.temperature_mode, c.j_source, c.roof_open = "const", "clear_sky", [("08:30", "10:00")]
    e.precursors = ["APINENE"]
    e.injections = [juice.Injection(name="O3", products="O3", total_ppb="60", start="08:05"),
                    juice.Injection(name="apinene", products="APINENE", total_ppb="20", start="08:20")]
    e.partitioning, e.volatility_method, e.particles, e.particle_n_cm3 = True, "li2016", "mono", 2e4
    return juice.run(e, str(tmp_path_factory.mktemp("phases")), solver="fast", kv_pd=False, script=False)


def _df(r):
    d = load_column_file(r.cs)
    return d[0] if isinstance(d, tuple) else d


def test_carbon_closes(run):
    d = _df(run)
    b = pb.build(d, load_species_db(), run.cs)
    assert len(b.species) > 50 and len(b.log_cstar) == len(b.species)
    end = b.at("C", -1)
    assert set(end) == set(pb.PHASES) | {pb.DILUTED} and end["particles"] > 0 and end["wall (vapour)"] > 0
    injected = 20e-9 * float(d["M"].values[0]) * 10 * pb.ug_per_molec("C")          # 20 ppb of a C10
    co = (d["CO"].values[-1] + (d["CO2"].values[-1] - d["CO2"].values[0] if "CO2" in d else 0)) * pb.ug_per_molec("C")
    assert sum(end.values()) + co == pytest.approx(injected, rel=0.05)
    assert b.at("C", 0)["particles"] < 1e-6                      # every species starts at 1 molecule cm-3
    # the volatility bars hold exactly the partitioning species' carbon in the condensed phases
    h = b.by_volatility("C", -1)
    assert sum(h["particles"]) == pytest.approx(end["particles"], rel=1e-9)
    # with the results of juice.run (a plain row index) the times are passed in
    r = run.results()
    b2 = pb.build(r.df, load_species_db(), run.cs, t_seconds=np.asarray(r.hour) * 3600)
    assert b2.t_hours[-1] == pytest.approx(b.t_hours[-1]) and b2.at("C", -1)["particles"] == pytest.approx(end["particles"])


def test_a_run_without_partitioning_says_so(tmp_path):
    import pandas as pd
    with pytest.raises(ValueError, match="no gas-particle-wall partitioning"):
        pb.build(pd.DataFrame({"APINENE": [1.0, 2.0]}, index=[0.0, 60.0]), load_species_db())


def test_phases_tab(run):
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PyQt6.QtWidgets import QApplication
    app = QApplication.instance() or QApplication([])                # noqa: F841
    from juice.gui import results_viewer as rv
    v = rv.ResultsViewer(rv.load_column_file(run.cs), csv_path=run.cs)
    names = [v.right_tabs.tabText(i) for i in range(v.right_tabs.count())]
    assert "Phases" in names
    v.right_tabs.setCurrentIndex(names.index("Phases"))
    v._update_phases()
    assert len(v.ph_fig.axes) >= 2 and "partitioning species" in v.ph_info.text()
    v.ph_element.setCurrentIndex(1)                                     # oxygen
    assert "oxygen" in v.ph_fig.axes[0].get_title(loc="left")
    v.close()
