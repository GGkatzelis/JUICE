"""The study window opened without a screen, on a study that has really run (the `small_study` fixture, a 2 x 2
isopleth): what it reads back, the results table, the isopleth and its regimes, the video frames, the exports."""
import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import pytest
from PyQt6.QtWidgets import QApplication, QFileDialog

from juice.core import study as S
from juice.gui import study_window as sw

TOML = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "settings.toml")


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication([])


@pytest.fixture(scope="module")
def window(app, small_study):
    w = sw.StudyWindow(TOML, small_study)
    yield w
    w.close()


def test_the_study_is_read_back_into_the_window(window, small_study):
    st = S.load_study(small_study)
    assert window.mode_sweep.isChecked() and window.exp_path.text() == st.experiments[0]
    again = window._to_study()
    assert again.experiments == st.experiments and again.name == "iso" and again.folder == st.folder
    assert [(a.members(), a.values) for a in again.axes] == [(a.members(), a.values) for a in st.axes]
    # bug of 25 Nov 2026: an axis without a label (a study written from code) came back as VOC_x, named after the
    # window's "All VOC injections" row, and the finished runs, stored under ethanol_x, were not found
    assert [a.name() for a in again.axes] == ["NO2_t0_ppb", "ethanol_x"]
    assert window.ready.text().startswith("Ready: a sweep over 2 axis/axes of uv_ethanol")


def test_the_count_says_that_everything_has_run(window):
    window._count()
    assert "1 experiment(s) × 4 grid point(s) = 4 runs — all 4 already done" in window.count.text()


def test_the_check_finds_nothing_to_fix(window):
    assert window._check()
    assert "nothing to fix" in window.msgs.toPlainText()
    assert "axis 'NO2_t0_ppb': initial NO2 over 5, 50" in window.msgs.toPlainText()


def test_the_results_table_holds_every_run(window):
    window._refresh_results()
    assert window.results.rowCount() == 4
    status = [window.results.horizontalHeaderItem(c).text() for c in range(window.results.columnCount())].index("status")
    assert {window.results.item(r, status).text() for r in range(4)} == {"done"}


def test_the_window_opens_on_the_ozone_isopleth_with_its_regimes(window):
    assert window.plot_metric.currentData() == "o3_max"
    ax = window.fig.axes[0]
    assert ax.get_title(loc="left").startswith("O3 maximum — uv_ethanol") and ax.collections       # a filled contour
    assert "isopleth over ethanol_x × NO2_t0_ppb; 4 grid points" in window.plot_note.text()      # NOx on the y axis


@pytest.mark.parametrize("key", S.FOCUS["ozone"][1] + ["oh_mean_lit", "voc_consumed"])
def test_every_result_can_be_plotted(window, key):
    if window.plot_metric.findData(key) < 0:
        window.plot_metric.addItem(key, key)
    window.plot_metric.setCurrentIndex(window.plot_metric.findData(key))
    window._plot()
    assert window.fig.axes and window.fig.axes[0].get_title(loc="left")


def test_log_axes(window):
    window.plot_log.setChecked(True)
    try:
        assert "log10 scale" in window.fig.axes[0].get_xlabel()
    finally:
        window.plot_log.setChecked(False)


@pytest.mark.parametrize("quantity", list(S.TS_QUANTITIES))
def test_the_video_draws_every_quantity(window, quantity):
    if not window._v_data:
        window._v_load()
    assert window.v_slider.maximum() == 48                              # 08:00 to 12:00 in 5-minute frames
    window.v_qty.setCurrentIndex(window.v_qty.findData(quantity))
    window.v_slider.setValue(30)                                        # 10:30
    window._v_draw(rescale=True)
    assert window.v_time.text().startswith("10:30 UTC")
    assert window.v_fig.axes[0].get_title(loc="left").startswith(S.TS_QUANTITIES[quantity][0])


def test_the_video_can_be_written(window, tmp_path, monkeypatch):
    pytest.importorskip("imageio")
    out = tmp_path / "iso.gif"
    monkeypatch.setattr(QFileDialog, "getSaveFileName", lambda *a, **k: (str(out), ""))
    window.v_from.setText("10:00"); window.v_to.setText("10:10")
    try:
        window._v_export()
    finally:
        window.v_from.setText(""); window.v_to.setText("")
    assert out.stat().st_size > 0


def test_results_image_and_study_files_are_written(window, tmp_path, monkeypatch):
    csv_path, png = tmp_path / "res.csv", tmp_path / "plot.png"
    monkeypatch.setattr(QFileDialog, "getSaveFileName", lambda *a, **k: (str(csv_path), ""))
    window._export()
    lines = csv_path.read_text(encoding="utf-8").splitlines()
    assert lines[0].startswith("run,experiment,NO2_t0_ppb,ethanol_x,status") and len(lines) == 5
    monkeypatch.setattr(QFileDialog, "getSaveFileName", lambda *a, **k: (str(png), ""))
    window._save_image()
    assert png.stat().st_size > 0
    copy = window._save(str(tmp_path / "copy.study.json"))
    again = S.load_study(copy)
    assert [a.values for a in again.axes] == [[5.0, 50.0], [0.5, 2.0]] and again.experiments == window.study.experiments


def test_an_empty_window_says_what_is_missing(app):
    w = sw.StudyWindow(TOML)
    try:
        assert w.ready.text().startswith("Still needed:") and "the experiment to sweep" in w.ready.text()
        w.fig.clear(); w._plot()
        assert w.plot_note.text() == "no finished runs yet"
    finally:
        w.close()
