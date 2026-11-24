"""report_basket.py: the figure list of a result — adding, ordering, captions, surviving a restart, and the
numbered set of files with the index and the slides."""
import os

import matplotlib

matplotlib.use("Agg")
import pytest
from matplotlib.figure import Figure

from juice.gui.report_basket import FigureList, _slug


def figure(value=1.0):
    fig = Figure(figsize=(4, 3))
    fig.add_subplot().plot([0, 1], [0, value])
    return fig


@pytest.fixture
def figures(tmp_path):
    cs = tmp_path / "run" / "day.cs"
    cs.parent.mkdir()
    cs.write_text("cs", encoding="utf-8")
    return FigureList.for_result(str(cs))


def test_the_list_lives_next_to_the_result(figures, tmp_path):
    assert figures.store_dir == str(tmp_path / "run" / "day.figures") and figures.run_name == "day" and figures.items == []


def test_file_names_from_titles():
    assert _slug("NOx & ozone (lit)") == "NOx_ozone_lit" and _slug("§§§") == "figure" and len(_slug("x" * 100)) == 48


def test_added_figures_are_kept_as_vector_and_bitmap_files(figures):
    item = figures.add(figure(), "NOx and ozone", source="Overview", caption="first try")
    assert (item.title, item.source, item.caption) == ("NOx and ozone", "Overview", "first try")
    assert item.svg.endswith(".svg") and item.png.endswith(".png")
    for name in (item.svg, item.png):
        assert os.path.getsize(os.path.join(figures.store_dir, name)) > 0
    assert (item.width_in, item.height_in) == (10.0, 6.5)             # a small panel is exported at presentation size
    fig = figure()
    figures.add(fig, "size check")
    assert tuple(fig.get_size_inches()) == (4.0, 3.0)                  # and the figure on screen keeps its own
    wide = Figure(figsize=(12, 5)); wide.add_subplot()
    assert (figures.add(wide, "wide").width_in, figures.items[-1].height_in) == (12.0, 5.0)
    assert figures.add(figure(), "no caption").caption == "day: no caption"
    second = figures.add(figure(), "NOx and ozone")
    assert second.key != item.key and second.svg != item.svg           # the same title twice: two figures


def test_order_captions_and_removal(figures):
    a, b, c = (figures.add(figure(), title) for title in ("a", "b", "c"))
    figures.move(c.key, -1)
    assert [i.title for i in figures.items] == ["a", "c", "b"]
    figures.move(a.key, -1); figures.move(b.key, +1)                   # already at the ends: nothing happens
    assert [i.title for i in figures.items] == ["a", "c", "b"]
    figures.set_caption(b.key, "the last one", title="B")
    assert (figures.get(b.key).caption, figures.get(b.key).title) == ("the last one", "B")
    figures.remove(c.key)
    assert [i.title for i in figures.items] == ["a", "B"] and figures.get(c.key) is None
    assert not os.path.exists(os.path.join(figures.store_dir, c.svg))
    figures.clear()
    assert figures.items == []


def test_the_list_survives_a_restart_without_the_figures_whose_files_are_gone(figures):
    a = figures.add(figure(), "kept")
    b = figures.add(figure(), "lost")
    for name in (b.svg, b.png):
        os.remove(os.path.join(figures.store_dir, name))
    again = FigureList(figures.store_dir, figures.run_name)
    assert [i.title for i in again.items] == ["kept"] and again.items[0].key == a.key
    assert again.add(figure(), "new").key not in (a.key, b.key)        # keys are not used twice


def test_an_image_from_a_file_is_added_like_a_figure(figures, tmp_path):
    png = tmp_path / "shot.png"
    figure().savefig(str(png))
    item = figures.add_image(str(png), "screenshot", source="Chemical Space")
    assert item.png and os.path.isfile(os.path.join(figures.store_dir, item.png)) and figures.items == [item]


def test_the_built_set_is_numbered_in_list_order_with_an_index(figures, tmp_path):
    figures.add(figure(), "NOx and ozone", source="Overview", caption="the day")
    figures.add(figure(), "OH budget", source="Process Analysis", caption="sources and sinks")
    out = figures.build(str(tmp_path / "talk"), title="SAPHIR day", subtitle="for the group meeting")
    names = sorted(os.path.basename(p) for p in out["figures"])
    assert names == ["01_NOx_and_ozone.png", "01_NOx_and_ozone.svg", "02_OH_budget.png", "02_OH_budget.svg"]
    index = open(out["index"], encoding="utf-8").read()
    assert index.startswith("# SAPHIR day") and "for the group meeting" in index
    assert "1. **NOx and ozone** (Overview) — the day" in index and "`figures/02_OH_budget.svg`" in index
    assert os.path.basename(out["index"]) == "day_figures.md"


def test_slides_are_written_when_python_pptx_is_there(figures, tmp_path):
    pptx = pytest.importorskip("pptx")
    figures.add(figure(), "NOx and ozone", caption="the day")
    figures.add(figure(), "OH budget", caption="sources and sinks")
    out = figures.build(str(tmp_path / "talk"), title="SAPHIR day")
    assert out["warning"] == "" and os.path.basename(out["pptx"]) == "day_figures.pptx"
    deck = pptx.Presentation(out["pptx"])
    assert len(deck.slides) == 3                                       # a title slide and one per figure
    text = " ".join(shape.text_frame.text for slide in deck.slides for shape in slide.shapes if shape.has_text_frame)
    assert "SAPHIR day" in text and "OH budget" in text and "sources and sinks" in text


def test_an_empty_list_cannot_be_built(figures, tmp_path):
    with pytest.raises(ValueError, match="empty"):
        figures.build(str(tmp_path / "talk"))
