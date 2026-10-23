"""
report_basket.py — the figure list of a result: collect figures while browsing, build a document.

No Qt here.  `FigureList` keeps every added figure as an SVG (editable) and a PNG (for the slides)
in `<run folder>/<stem>.figures/` together with `list.json`, so the list survives closing the
viewer and nothing is lost if it crashes.  `build(out_dir)` writes the final set:

    <out_dir>/figures/NN_<title>.svg + .png     numbered in list order
    <out_dir>/<run>_figures.pptx                title slide + one slide per figure (needs python-pptx)
    <out_dir>/<run>_figures.md                  index with the captions

The viewer's dock (report_dock.py) is the user interface; this module can also be used from a script:

    fl = FigureList.for_result(r"…/saphir_20240709.cs")
    fl.add(fig, "NOx and ozone", source="Overview")
    fl.build(r"D:/talks/2026-10 SAPHIR")
"""
from __future__ import annotations

import datetime as _dt
import json
import os
import re
import shutil
from dataclasses import dataclass, asdict, field
from typing import List, Optional


def _slug(text: str, n: int = 48) -> str:
    s = re.sub(r"[^A-Za-z0-9]+", "_", text).strip("_")
    return (s[:n] or "figure")


@dataclass
class FigureItem:
    key: str                    # unique within the list
    title: str
    caption: str
    source: str                 # the tab it came from
    svg: str                    # file names inside the store folder
    png: str
    added: str                  # ISO time
    width_in: float = 10.0
    height_in: float = 6.0
    notes: str = ""


class FigureList:
    def __init__(self, store_dir: str, run_name: str):
        self.store_dir = store_dir
        self.run_name = run_name
        self.items: List[FigureItem] = []
        self._next = 1
        self.load()

    # ── persistence ──
    @classmethod
    def for_result(cls, cs_path: str) -> "FigureList":
        stem = os.path.splitext(os.path.basename(cs_path))[0]
        return cls(os.path.join(os.path.dirname(os.path.abspath(cs_path)), stem + ".figures"), stem)

    @property
    def index_path(self) -> str:
        return os.path.join(self.store_dir, "list.json")

    def load(self) -> None:
        self.items = []
        try:
            with open(self.index_path, encoding="utf-8") as f:
                data = json.load(f)
            for d in data.get("items", []):
                it = FigureItem(**{k: d.get(k, "") for k in FigureItem.__dataclass_fields__})
                if os.path.isfile(os.path.join(self.store_dir, it.svg)) or os.path.isfile(os.path.join(self.store_dir, it.png)):
                    self.items.append(it)
            self._next = int(data.get("next", len(self.items) + 1))
        except Exception:
            self._next = len(self.items) + 1

    def save(self) -> None:
        os.makedirs(self.store_dir, exist_ok=True)
        with open(self.index_path, "w", encoding="utf-8") as f:
            json.dump({"run": self.run_name, "next": self._next, "items": [asdict(i) for i in self.items]}, f, indent=1)

    # ── editing ──
    def add(self, fig, title: str, source: str = "", caption: Optional[str] = None, dpi: int = 200) -> FigureItem:
        """Snapshot a matplotlib Figure now (figures are redrawn later; the files are what we keep)."""
        os.makedirs(self.store_dir, exist_ok=True)
        key = f"{self._next:03d}"; self._next += 1
        base = f"{key}_{_slug(title)}"
        svg, png = base + ".svg", base + ".png"
        orig = tuple(float(x) for x in fig.get_size_inches())
        w, h = orig
        if w < 8.0:                     # small on-screen panels: export at a presentation size, independent of the window
            w, h = 10.0, min(6.5, max(4.0, 10.0 * orig[1] / max(orig[0], 1e-6)))
            fig.set_size_inches(w, h, forward=False)
        try:
            fig.savefig(os.path.join(self.store_dir, svg), format="svg", bbox_inches="tight", facecolor="white")
            fig.savefig(os.path.join(self.store_dir, png), format="png", dpi=dpi, bbox_inches="tight", facecolor="white")
        finally:
            if (w, h) != orig:
                fig.set_size_inches(*orig, forward=False)
                try:
                    fig.canvas.draw_idle()
                except Exception:
                    pass
        it = FigureItem(key=key, title=title, caption=caption if caption is not None else f"{self.run_name}: {title}",
                        source=source, svg=svg, png=png, added=_dt.datetime.now().isoformat(timespec="seconds"), width_in=w, height_in=h)
        self.items.append(it); self.save()
        return it

    # ── several figures → one figure with sub-panels ──
    def add_composite(self, figs, title: str, source: str = "", caption: Optional[str] = None, grid=None,
                      panel_size=(6.4, 3.8), dpi: int = 200) -> FigureItem:
        """One SVG and one PNG with the given figures as sub-panels in a grid (rows, cols; chosen from
        the count when not given).  Every panel is exported at `panel_size` inches; in the SVG each panel
        stays editable vector content (nested <svg>, ids prefixed so they cannot clash)."""
        import math, io
        import xml.etree.ElementTree as ET
        figs = [f for f in figs if f is not None]
        if not figs:
            raise ValueError("no figures")
        n = len(figs)
        if grid:
            rows, cols = grid
        else:
            cols = 1 if n == 1 else (2 if n <= 4 else 3); rows = math.ceil(n / cols)
        os.makedirs(self.store_dir, exist_ok=True)
        key = f"{self._next:03d}"; self._next += 1
        base = f"{key}_{_slug(title)}"; svg, png = base + ".svg", base + ".png"
        pw, ph = panel_size
        svgs, pngs = [], []
        for fig in figs:
            orig = tuple(float(x) for x in fig.get_size_inches())
            fig.set_size_inches(pw, ph, forward=False)
            try:
                b = io.BytesIO(); fig.savefig(b, format="svg", facecolor="white"); svgs.append(b.getvalue().decode("utf-8"))
                b = io.BytesIO(); fig.savefig(b, format="png", dpi=dpi, facecolor="white"); pngs.append(b.getvalue())
            finally:
                fig.set_size_inches(*orig, forward=False)
                try:
                    fig.canvas.draw_idle()
                except Exception:
                    pass
        # PNG: paste the panels on a white sheet
        from PIL import Image
        ims = [Image.open(io.BytesIO(b)).convert("RGB") for b in pngs]
        cw, ch = max(i.size[0] for i in ims), max(i.size[1] for i in ims)
        sheet = Image.new("RGB", (cols * cw, rows * ch), "white")
        for k, im in enumerate(ims):
            r, c = divmod(k, cols); sheet.paste(im, (c * cw + (cw - im.size[0]) // 2, r * ch + (ch - im.size[1]) // 2))
        sheet.save(os.path.join(self.store_dir, png))
        # SVG: nest each panel's document, ids prefixed per panel
        NS = "http://www.w3.org/2000/svg"; ET.register_namespace("", NS); ET.register_namespace("xlink", "http://www.w3.org/1999/xlink")
        root = ET.Element(f"{{{NS}}}svg", {"width": f"{cols * pw:g}in", "height": f"{rows * ph:g}in",
                                          "viewBox": f"0 0 {cols * pw * 72:g} {rows * ph * 72:g}", "version": "1.1"})
        ET.SubElement(root, f"{{{NS}}}rect", {"width": "100%", "height": "100%", "fill": "white"})
        for k, text in enumerate(svgs):
            text = re.sub(r'id="([^"]+)"', lambda m: f'id="p{k}_{m.group(1)}"', text)
            text = re.sub(r'url\(#([^)]+)\)', lambda m: f'url(#p{k}_{m.group(1)})', text)
            text = re.sub(r'href="#([^"]+)"', lambda m: f'href="#p{k}_{m.group(1)}"', text)
            el = ET.fromstring(text)
            r, c = divmod(k, cols)
            nested = ET.SubElement(root, f"{{{NS}}}svg", {"x": f"{c * pw * 72:g}", "y": f"{r * ph * 72:g}", "width": f"{pw * 72:g}", "height": f"{ph * 72:g}",
                                                        "viewBox": el.get("viewBox", f"0 0 {pw * 72:g} {ph * 72:g}")})
            for child in list(el):
                nested.append(child)
        ET.ElementTree(root).write(os.path.join(self.store_dir, svg), encoding="utf-8", xml_declaration=True)
        it = FigureItem(key=key, title=title, caption=caption if caption is not None else f"{self.run_name}: {title}",
                        source=source, svg=svg, png=png, added=_dt.datetime.now().isoformat(timespec="seconds"),
                        width_in=cols * pw, height_in=rows * ph)
        self.items.append(it); self.save()
        return it

    def add_image(self, png_path: str, title: str, source: str = "", caption: Optional[str] = None) -> FigureItem:
        """A raster snapshot (e.g. a whole page grabbed from the screen); no SVG for these."""
        os.makedirs(self.store_dir, exist_ok=True)
        key = f"{self._next:03d}"; self._next += 1
        png = f"{key}_{_slug(title)}.png"
        shutil.copyfile(png_path, os.path.join(self.store_dir, png))
        w, h = 10.0, 6.0
        try:
            from PIL import Image
            with Image.open(os.path.join(self.store_dir, png)) as im:
                w, h = im.size[0] / 100.0, im.size[1] / 100.0
        except Exception:
            pass
        it = FigureItem(key=key, title=title, caption=caption if caption is not None else f"{self.run_name}: {title}",
                        source=source, svg="", png=png, added=_dt.datetime.now().isoformat(timespec="seconds"), width_in=w, height_in=h)
        self.items.append(it); self.save()
        return it

    def remove(self, key: str) -> None:
        it = self.get(key)
        if it is None:
            return
        self.items.remove(it)
        for fn in (it.svg, it.png):
            p = os.path.join(self.store_dir, fn) if fn else ""
            if p and os.path.isfile(p):
                try:
                    os.remove(p)
                except OSError:
                    pass
        self.save()

    def move(self, key: str, delta: int) -> None:
        it = self.get(key)
        if it is None:
            return
        i = self.items.index(it); j = max(0, min(len(self.items) - 1, i + delta))
        if i != j:
            self.items.insert(j, self.items.pop(i)); self.save()

    def get(self, key: str) -> Optional[FigureItem]:
        return next((i for i in self.items if i.key == key), None)

    def set_caption(self, key: str, caption: str, title: Optional[str] = None) -> None:
        it = self.get(key)
        if it is not None:
            it.caption = caption
            if title:
                it.title = title
            self.save()

    def clear(self) -> None:
        for it in list(self.items):
            self.remove(it.key)

    # ── building ──
    def build(self, out_dir: str, pptx_name: Optional[str] = None, title: Optional[str] = None,
              subtitle: Optional[str] = None) -> dict:
        """Write the numbered SVG/PNG set, the index and the PowerPoint into out_dir.
        Returns {'folder', 'figures': [paths], 'pptx': path or '', 'index': path, 'warning': str}."""
        if not self.items:
            raise ValueError("the figure list is empty")
        fig_dir = os.path.join(out_dir, "figures"); os.makedirs(fig_dir, exist_ok=True)
        written: List[str] = []; rows = []
        for n, it in enumerate(self.items, 1):
            base = f"{n:02d}_{_slug(it.title)}"
            for fn, ext in ((it.svg, ".svg"), (it.png, ".png")):
                src = os.path.join(self.store_dir, fn) if fn else ""
                if src and os.path.isfile(src):
                    dst = os.path.join(fig_dir, base + ext); shutil.copyfile(src, dst); written.append(dst)
            rows.append((n, it, base))
        index = os.path.join(out_dir, f"{self.run_name}_figures.md")
        with open(index, "w", encoding="utf-8") as f:
            f.write(f"# {title or self.run_name}\n\n{subtitle or ''}\n\n")
            for n, it, base in rows:
                f.write(f"{n}. **{it.title}** ({it.source}) — {it.caption}  \n   `figures/{base}.svg`" + (" / `.png`\n" if it.png else "\n"))
        out = {"folder": out_dir, "figures": written, "index": index, "pptx": "", "warning": ""}
        try:
            out["pptx"] = self._write_pptx(os.path.join(out_dir, pptx_name or f"{self.run_name}_figures.pptx"), rows, fig_dir,
                                           title or self.run_name, subtitle)
        except ImportError:
            out["warning"] = "python-pptx is not installed (pip install python-pptx): figures and index written, no slides"
        return out

    def _write_pptx(self, path: str, rows, fig_dir: str, title: str, subtitle: Optional[str]) -> str:
        from pptx import Presentation
        from pptx.util import Inches, Pt
        prs = Presentation()
        prs.slide_width, prs.slide_height = Inches(13.333), Inches(7.5)
        s = prs.slides.add_slide(prs.slide_layouts[0])
        s.shapes.title.text = title
        if len(s.placeholders) > 1:
            s.placeholders[1].text = subtitle or f"JUICE results · {_dt.date.today():%Y-%m-%d} · {len(rows)} figures"
        blank = prs.slide_layouts[6]
        box_w, box_h, top = 12.5, 5.7, 1.15
        for n, it, base in rows:
            sl = prs.slides.add_slide(blank)
            tb = sl.shapes.add_textbox(Inches(0.4), Inches(0.25), Inches(12.5), Inches(0.8))
            p = tb.text_frame.paragraphs[0]; p.text = f"{n}. {it.title}"; p.font.size = Pt(26); p.font.bold = True
            png = os.path.join(fig_dir, base + ".png")
            if os.path.isfile(png):
                ar = (it.width_in / it.height_in) if it.height_in else 1.6
                try:
                    from PIL import Image
                    with Image.open(png) as im:
                        ar = im.size[0] / im.size[1]
                except Exception:
                    pass
                w, h = box_w, box_w / ar
                if h > box_h:
                    h, w = box_h, box_h * ar
                left = (13.333 - w) / 2
                sl.shapes.add_picture(png, Inches(left), Inches(top), Inches(w), Inches(h))
            cb = sl.shapes.add_textbox(Inches(0.4), Inches(6.95), Inches(12.5), Inches(0.45))
            p = cb.text_frame.paragraphs[0]; p.text = it.caption; p.font.size = Pt(12); p.font.italic = True
            sl.notes_slide.notes_text_frame.text = f"{it.source} · added {it.added} · figures/{base}.svg"
        prs.save(path)
        return path
