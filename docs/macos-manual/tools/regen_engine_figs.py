#!/usr/bin/env python3
"""Regenerate the manual figures that come from the running engine.

Each figure has a journal in tools/figs/ that the CLI suite's driver (cli_tests/cli_cmd_gate.py)
runs on the /png device; the journal copies its giza pages into build/figs/<dir>/ and the SPEC
below assembles media/fig_engine_NN.png from them (one page, or pages side by side).  A journal
whose EXPECT fails stops the figure, so a figure cannot silently rot into a wrong picture
(2026-10-17 for figures 31, 56, 58; 2026-10-19 for the rest, Dave: "regenerate all the figures
using the current macos").  The figures that are not engine output (diagrams and schematics:
1, 2, 5, 7-9, 11-30, 32, 35, 38, 40-45, 50, 59, 65, 66) keep their rescued artwork; figure 46 is remade
on CassWithExitPupil.in (its 3.2 beam-expander deck is not in the tree).

Usage (from docs/macos-manual):  python3 tools/regen_engine_figs.py [<macos binary>] [--assemble-only]
(--assemble-only reuses the pages in build/figs/ -- for a label or layout change without a rerun.)
Default binary: ../../build_release/bin/macos.
"""
import csv, os, shutil, subprocess, sys, tempfile
from pathlib import Path
from PIL import Image

HERE = Path(__file__).resolve().parent
MAN = HERE.parent
ROOT = MAN.parent.parent
MEDIA = MAN / "src" / "media"
FIGS = MAN / "build" / "figs"

# figure number -> (page directory, [page indices side by side], [one label under each panel])
# The labels are drawn into the image under their panels, centered, the way the 3.2 manual set them
# (Dave 2026-10-19: "positioned under their respective panels"); the FIGURE line below the image is
# the caption, centered by tools/figure_captions.lua.
SPEC = {
    3: ("fig03", [0, 1], ["Spot diagram at entrance pupil", "Cross-section of ray-trace"]),
    4: ("fig04", [0, 1], ["Focal plane spot diagram\n(unaberrated system)",
                         "Focal plane spot diagram\n(primary mirror tilted 1e-4 radians)"]),
    6: ("fig06", [0, 1], ["Unaberrated (defocused)", "Aberrated OPD (coma)"]),
    10: ("fig10", [0, 1], ["Unaberrated (defocused)", "Aberrated (coma)"]),
    31: ("fig31", [0, 1], ["Elements 0-8", "The relay lens, elements 6-8"]),
    33: ("fig33", [0], []), 34: ("fig33", [1], []),
    36: ("fig36", [0, 1], ["Exit pupil (element 4)", "Detector (element 5)"]),
    37: ("fig36", [2, 3], ["OBS=ALL", "OBS=NEGATIVE"]),
    39: ("fig39", [1], []),
    46: ("fig46", [0, 1], ["Uniform beam", "Gaussian beam"]),
    47: ("fig47", [0], []), 48: ("fig47", [1], []),
    49: ("fig47", [2, 3, 4], ["Pupil (element 10)", "Lyot stop (element 11)", "Detector (element 16)"]),
    51: ("fig47", [5, 6], ["Unocculted", "Occulted"]),
    52: ("fig47", [7, 8], ["LOG10 stretch", "SQRT stretch"]),
    53: ("fig47", [9], []), 54: ("fig47", [10], []), 55: ("fig47", [11], []), 57: ("fig47", [12], []),
    56: ("fig56", [1, 2], ["0.1 \u00b5m kernel (below the pixel)", "20 \u00b5m kernel (14 pixels)"]),
    58: ("fig58", [2], []),
    60: ("fig60", [0, 1], ["Amplitude", "Phase"]),
    61: ("fig60", [2, 3], ["Real part", "Imaginary part"]),
    62: ("fig47", [13], []), 63: ("fig60", [5], []), 64: ("fig60", [6], []),
}
FONT = "/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf"

args = [a for a in sys.argv[1:] if not a.startswith("--")]
assemble_only = "--assemble-only" in sys.argv      # reuse build/figs/ pages (labels, layout)
binary = args[0] if args else str(ROOT / "build_release" / "bin" / "macos")
gate = ROOT / "cli_tests" / "cli_cmd_gate.py"
if FIGS.exists() and not assemble_only:
    shutil.rmtree(FIGS)
with tempfile.TemporaryDirectory() as td:
    if assemble_only:
        td = None
    if td is not None:
        out = Path(td) / "figs.csv"
        subprocess.run([sys.executable, str(gate), binary, str(HERE / "figs"), str(out)], check=True)
        bad = [r for r in csv.DictReader(open(out)) if r["status"] != "pass"]
        if bad:
            sys.exit("figure journals failed: " + ", ".join(f"{r['journal']} ({r['status']}: {r['missed']})" for r in bad))


def page(d, n):
    p = FIGS / d / f"giza_{n:04d}.png"
    if not p.exists():
        sys.exit(f"missing page {p}")
    return Image.open(p).convert("RGB")


def trim(im, pad=12):
    """Cut the white margin giza leaves around a page (the DRAW pages are strips in a 4:3 page)."""
    from PIL import ImageChops
    bg = Image.new("RGB", im.size, "white")
    bbox = ImageChops.difference(im, bg).getbbox()
    if not bbox:
        return im
    x0, y0, x1, y1 = bbox
    return im.crop((max(0, x0 - pad), max(0, y0 - pad), min(im.width, x1 + pad), min(im.height, y1 + pad)))


from PIL import ImageDraw, ImageFont


def labelled(im, text, font):
    """The panel with its label centered beneath (one or two lines)."""
    if not text:
        return im
    lines = text.split("\n")
    draw = ImageDraw.Draw(im)
    lh = int(font.size * 1.3)
    out = Image.new("RGB", (im.width, im.height + 10 + lh * len(lines) + 6), "white")
    out.paste(im, (0, 0))
    draw = ImageDraw.Draw(out)
    y = im.height + 8
    for ln in lines:
        tw = draw.textlength(ln, font=font)
        draw.text(((im.width - tw) / 2, y), ln, fill="black", font=font)
        y += lh
    return out


for fig, (d, pages, labels) in sorted(SPEC.items()):
    ims = [trim(page(d, n)) for n in pages]
    if labels:
        assert len(labels) == len(ims), fig
        font = ImageFont.truetype(FONT, max(18, int(min(im.height for im in ims) * 0.045)))
        ims = [labelled(im, lab, font) for im, lab in zip(ims, labels)]
    h = max(im.height for im in ims)
    w = sum(im.width for im in ims) + 20 * (len(ims) - 1)
    out = Image.new("RGB", (w, h), "white")
    x = 0
    for im in ims:
        out.paste(im, (x, 0))
        x += im.width + 20
    out.save(MEDIA / f"fig_engine_{fig:02d}.png")
    print(f"fig_engine_{fig:02d}.png  <- {d} pages {pages}  {out.size}")
