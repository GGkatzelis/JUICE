"""plot_interaction.py, the functions that work on a plot without a mouse: zooming an axis around a point, the
range of the data in view, switching an axis between linear and logarithmic, click or drag."""
from types import SimpleNamespace as NS

import matplotlib

matplotlib.use("Agg")
import numpy as np
import pytest
from matplotlib.figure import Figure

from juice.gui import plot_interaction as pi


@pytest.fixture
def ax():
    axes = Figure().add_subplot()
    axes.plot([0, 1, 2, 3, 4], [10.0, 20.0, 400.0, 20.0, 10.0])
    axes.set_xlim(0, 4); axes.set_ylim(0, 400)
    return axes


def test_zoom_keeps_the_point_under_the_mouse(ax):
    pi.zoom_axis(ax, "x", 1.0, 0.5)
    assert ax.get_xlim() == pytest.approx((0.5, 2.5))                  # the point 1.0 stays a quarter from the left
    pi.zoom_axis(ax, "y", 100.0, 2.0)
    assert ax.get_ylim() == pytest.approx((-100.0, 700.0))


def test_zoom_on_a_logarithmic_axis_works_in_decades(ax):
    ax.set_yscale("log"); ax.set_ylim(1.0, 1e4)
    pi.zoom_axis(ax, "y", 100.0, 0.5)
    assert ax.get_ylim() == pytest.approx((10.0, 1000.0))
    pi.zoom_axis(ax, "y", -5.0, 0.5)                                   # a point that a log axis does not have
    assert ax.get_ylim() == pytest.approx((10.0, 1000.0)) and pi.is_log(ax, "y") and not pi.is_log(ax, "x")


def test_range_of_the_data_counts_only_what_is_in_view(ax):
    assert pi._data_range(ax, "y", positive=False) == (10.0, 400.0)
    ax.set_xlim(2.5, 4)
    assert pi._data_range(ax, "y", positive=False) == (10.0, 20.0)     # the peak at x = 2 is out of view
    assert pi._data_range(ax, "x", positive=False) == (0.0, 4.0)


def test_range_leaves_out_markers_of_the_plot_hidden_lines_and_bad_numbers(ax):
    ax.axvline(1.5)                                                    # a marker line: not data
    ax.axhline(9999.0)
    hidden, = ax.plot([0, 4], [5000.0, 5000.0]); hidden.set_visible(False)
    ax.plot([0, 1, 2], [np.nan, -3.0, np.inf])
    ax.scatter([1.0, 3.0], [700.0, 0.5])
    assert pi._data_range(ax, "y", positive=False) == (-3.0, 700.0)
    assert pi._data_range(ax, "y", positive=True) == (0.5, 700.0)
    assert pi._data_range(Figure().add_subplot(), "y", positive=False) is None


def test_switching_to_log_fits_the_positive_data_and_cuts_the_floor(ax):
    ax.plot([0, 4], [1e-9, 1e-9])                                      # the floor of a flushed chamber
    pi.set_log(ax, "y", True)
    lo, hi = ax.get_ylim()
    assert pi.is_log(ax, "y") and hi > 400.0 and 400.0 * 1e-5 * 0.5 < lo < 400.0 * 1e-5
    pi.set_log(ax, "y", True)                                          # already log: nothing changes
    assert ax.get_ylim() == (lo, hi)


def test_switching_back_to_linear_fits_all_the_data(ax):
    pi.set_log(ax, "y", True)
    pi.set_log(ax, "y", False)
    lo, hi = ax.get_ylim()
    assert not pi.is_log(ax, "y") and lo == pytest.approx(10.0 - 0.05 * 390.0) and hi == pytest.approx(400.0 + 0.05 * 390.0)


def test_an_axis_without_lines_is_left_to_matplotlib():
    axes = Figure().add_subplot()
    axes.bar([1, 2, 3], [1.0, 10.0, 100.0])
    pi.set_log(axes, "y", True)
    assert pi.is_log(axes, "y") and axes.get_ylim()[1] >= 100.0


def test_click_or_drag():
    assert pi.is_click((100, 100), NS(x=101, y=102)) and not pi.is_click((100, 100), NS(x=110, y=100))
    assert not pi.is_click((100, 100), NS(x=None, y=None))             # the mouse left the plot
