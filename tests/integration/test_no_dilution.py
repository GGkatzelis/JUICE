"""An experiment without dilution prepares and runs, and the EASY file reads back as 'no dilution' (bug of 10 Nov 2026:
'DILUTE is used in dilution reactions but not defined')."""
import juice
from juice.core import easy_builder as eb
from juice.core import experiment_spec as es


def _exp(name, dilution, partitioning=False):
    e = juice.new_experiment(name, chamber="SAPHIR")
    c = e.conditions
    c.source, c.date, c.start, c.end = "manual", "2024-07-09", "08:00", "09:00"
    c.temperature_mode, c.j_source, c.roof_open = "const", "clear_sky", [("08:10", "09:00")]
    e.precursors = ["C2H5OH"]
    e.injections = [juice.Injection(name="ethanol", products="C2H5OH", total_ppb="10", start="08:10")]
    e.chamber.dilution = dilution
    e.partitioning = partitioning
    return e


def test_easy_defines_dilute_as_zero_and_reads_back(tmp_path):
    text = es.render_easy(_exp("nodil", False))
    assert "DILUTE = CONST(0.)   ;no dilution in this experiment" in text and ".dil.mech" not in text
    p = tmp_path / "nodil.easy"
    p.write_text(text, encoding="latin-1")
    f = eb.load_form(str(p))
    statements = lambda t: [l.split(";")[0].strip() for l in t.splitlines() if l.split(";")[0].strip()]
    assert f.dilution is False and statements(eb.render_form(f)) == statements(text)
    on = es.render_easy(_exp("dil", True))
    assert "CONST(0.)   ;no dilution" not in on
    p2 = tmp_path / "dil.easy"; p2.write_text(on, encoding="latin-1")
    assert eb.load_form(str(p2)).dilution is True


def test_runs_without_dilution_also_with_partitioning(tmp_path):
    base = juice.run(_exp("dil", True), str(tmp_path / "dil"), solver="fast", kv_pd=False, script=False)
    nodil = juice.run(_exp("nodil", False, partitioning=True), str(tmp_path / "nodil"), solver="fast", kv_pd=False, script=False)
    # without the replenishment flow less ethanol is lost
    assert nodil.results().ppb("C2H5OH")[-1] > base.results().ppb("C2H5OH")[-1]
