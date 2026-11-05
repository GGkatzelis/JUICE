"""Plots say 'injection' instead of the Q column names; the files keep the Q names."""
import experiment_spec as es


def test_label_without_spec():
    assert es.injection_label("QNO") == "NO injection"
    assert es.injection_label("Q_ethanol") == "ethanol injection"
    assert es.injection_label("QH2O") == "H2O injection"


def test_label_uses_the_injection_name():
    s = es.ExperimentSpec(name="x")
    s.injections = [es.InjectionEvent(name="ethanol", products="C2H5OH", total_ppb="200", start="10:30", duration_min=10)]
    col = s.injections[0].enz_column()
    assert col.upper().startswith("Q")                 # the file column keeps its Q name
    assert es.injection_label(col, s) == "ethanol injection"
