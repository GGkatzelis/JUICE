"""The volatility map behind 'Compare methods': oxidation state, groups, the 80 % cloud, weights, the window."""
import os

import numpy as np
import pytest

import juice
from juice.chemistry import volatility as vol
from juice.chemistry import volatility_map as vmap


def test_oxidation_state():
    assert vmap.osc(1, 4, 0) == -4                       # methane
    assert vmap.osc(1, 0, 2) == 4                        # CO2
    assert vmap.osc(2, 6, 1) == pytest.approx(-2)        # ethanol
    assert vmap.osc(2, 3, 5, 1) == pytest.approx((10 - 3 - 5) / 2)   # PAN: the nitrogen counted as in a nitrate
    assert np.isnan(vmap.osc(0, 2, 1))


@pytest.mark.parametrize("smiles,expect", [
    ("CC(=O)OON(=O)=O", "PANs"), ("CC(O[N+](=O)[O-])C", "organic nitrates"), ("CC(ON(=O)=O)C", "organic nitrates"),
    ("Oc1ccccc1O", "aromatics"), ("CC(C)OO", "hydroperoxides and peracids"), ("CC(=O)OO", "hydroperoxides and peracids"),
    ("CC(=O)O", "acids"), ("CC(C)=O", "carbonyls"), ("CCO", "alcohols and other oxygenates"),
    ("CC1=CCC2CC1C2(C)C", "hydrocarbons")])
def test_groups_from_the_structure(smiles, expect):
    from rdkit.Chem import Descriptors
    from rdkit import Chem
    m = Chem.MolFromSmiles(smiles)
    n = {a: sum(1 for x in m.GetAtoms() if x.GetSymbol() == a) for a in "CON"}
    assert vmap.group_of(smiles, n["C"], n["O"], n["N"]) == expect


def test_cloud_holds_the_requested_fraction():
    rng = np.random.default_rng(1)
    x, y = rng.normal(0, 1, 4000), rng.normal(0, 1, 4000)
    from scipy.stats import gaussian_kde
    kde = gaussian_kde(np.vstack([x, y]))
    X, Y = np.meshgrid(np.linspace(-5, 5, 150), np.linspace(-5, 5, 150))
    Z = kde(np.vstack([X.ravel(), Y.ravel()])).reshape(X.shape)
    inside = kde(np.vstack([x, y])) >= vmap._hdr_level(Z, 0.8)
    assert inside.mean() == pytest.approx(0.8, abs=0.03)
    # the binned density used for drawing holds the same fraction
    Xb, Yb, Zb = vmap.density(x, y, np.ones_like(x), (-5, 5), (-5, 5))
    lev = vmap._hdr_level(Zb, 0.8)
    ix = np.clip(((x + 5) / 10 * 200).astype(int), 0, 199); iy = np.clip(((y + 5) / 10 * 200).astype(int), 0, 199)
    assert (Zb[iy, ix] >= lev).mean() == pytest.approx(0.8, abs=0.03)


def _map():
    from juice.chemistry.species_db import load_species_db
    names = ["PINONIC", "PINIC", "PAN", "APINENE", "HCHO", "PINONALDEHYDE", "C811OH", "LIMONIC", "C920PAN", "NLIMOOH",
             "MGLYOX", "GLYOX", "HCOOH", "CH3COOH", "CH3OOH", "PINAL", "C107OOH", "C108OOH", "C109OOH", "NORPINIC"]
    return vmap.from_species(names, load_species_db())


def test_map_and_weights(tmp_path):
    vm = _map()
    assert len(vm.names) >= 15 and set(vm.log_cstar) == {m for m in vol.METHODS if vol.METHODS[m].available()[0]}
    assert np.all(vm.weights(0.5) == 1)                                   # no run: each species once
    vm.peak = np.linspace(1.0, 100.0, len(vm.names))
    assert vm.weights(0.5) == pytest.approx(np.sqrt(vm.peak)) and vm.weights(0) == pytest.approx(np.ones(len(vm.names)))
    vm.peak = None
    from matplotlib.figure import Figure
    for draw in (lambda f: vmap.draw_methods(f, vm, "simpol", 0.5), lambda f: vmap.draw_groups(f, vm, "simpol", 0.5)):
        fig = Figure(); ax = draw(fig)
        (x0, x1), (y0, y1) = ax.get_xlim(), ax.get_ylim()
        assert x0 <= np.nanmin(vm.log_cstar["simpol"]) and y1 >= np.nanmax(vm.osc)


def test_window_from_a_run(tmp_path):
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PyQt6.QtWidgets import QApplication
    app = QApplication.instance() or QApplication([])               # noqa: F841
    from juice.gui.volatility_window import VolatilityWindow
    e = juice.new_experiment("vmap", chamber="SAPHIR")
    c = e.conditions
    c.source, c.date, c.start, c.end = "manual", "2024-07-09", "08:00", "09:00"
    c.temperature_mode, c.j_source, c.roof_open = "const", "clear_sky", [("08:10", "09:00")]
    e.precursors = ["APINENE"]
    e.injections = [juice.Injection(name="O3", products="O3", total_ppb="50", start="08:05"),
                    juice.Injection(name="apinene", products="APINENE", total_ppb="10", start="08:10")]
    run = juice.run(e, str(tmp_path / "run"), solver="fast", kv_pd=False, script=False)
    chosen = []
    w = VolatilityWindow(e, "simpol", os.path.join(run.out_dir, "vmap.cs"), on_choose=chosen.append)
    assert w.busy.isVisibleTo(w) and not w.tabs.isVisibleTo(w)          # opens at once, the work runs in a thread
    assert w._worker.wait(120000)
    app.processEvents()
    assert not w.busy.isVisibleTo(w) and w.tabs.isVisibleTo(w), w.error
    assert w.vm is not None and w.vm.peak is not None and w.weight.currentData() == 0.5
    assert np.nanmax(w.vm.peak[w.vm.names.index("PINONALDEHYDE")]) > 1e9 if "PINONALDEHYDE" in w.vm.names else True
    w.group_method.setCurrentIndex(w.group_method.findData("li2016")); w._use()
    assert chosen == ["li2016"] and w.chosen == "li2016"
    w.close()
