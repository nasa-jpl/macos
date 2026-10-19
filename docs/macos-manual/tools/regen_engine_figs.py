#!/usr/bin/env python3
"""Regenerate the manual figures that come from the running engine.

Each figure has a journal in tools/figs/ that the CLI suite's driver (cli_tests/cli_cmd_gate.py)
runs on the /png device; the journal copies its giza pages into build/figs/<dir>/ and the SPEC
below assembles media/fig_engine_NN.png from them (one page, or pages side by side).  A journal
whose EXPECT fails stops the figure, so a figure cannot silently rot into a wrong picture
(2026-10-17 for figures 31, 56, 58; 2026-10-19 for the rest, Dave: "regenerate all the figures
using the current macos").  The figures that are not engine output (diagrams and schematics:
1, 2, 5, 7-9, 11-30, 32, 35, 38, 40-45, 50, 59, 65, 66) and figure 46, whose 25-element beam
expander deck is not in the tree, keep their rescued artwork.

Usage (from docs/macos-manual):  python3 tools/regen_engine_figs.py [<macos binary>]
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

# figure number -> (page directory, [page indices side by side])
SPEC = {
    3: ("fig03", [0, 1]), 4: ("fig04", [0, 1]), 6: ("fig06", [0, 1]), 10: ("fig10", [0, 1]),
    31: ("fig31", [0, 1]), 33: ("fig33", [0]), 34: ("fig33", [1]),
    36: ("fig36", [0, 1]), 37: ("fig36", [2, 3]), 39: ("fig39", [1]),
    47: ("fig47", [0]), 48: ("fig47", [1]), 49: ("fig47", [2, 3, 4]), 51: ("fig47", [5, 6]),
    52: ("fig47", [7, 8]), 53: ("fig47", [9]), 54: ("fig47", [10]), 55: ("fig47", [11]),
    57: ("fig47", [12]), 56: ("fig56", [1, 2]), 58: ("fig58", [2]),
    60: ("fig60", [0, 1]), 61: ("fig60", [2, 3]), 62: ("fig60", [4]), 63: ("fig60", [5]), 64: ("fig60", [6]),
}

binary = sys.argv[1] if len(sys.argv) > 1 else str(ROOT / "build_release" / "bin" / "macos")
gate = ROOT / "cli_tests" / "cli_cmd_gate.py"
if FIGS.exists():
    shutil.rmtree(FIGS)
with tempfile.TemporaryDirectory() as td:
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


for fig, (d, pages) in sorted(SPEC.items()):
    ims = [trim(page(d, n)) for n in pages]
    h = max(im.height for im in ims)
    w = sum(im.width for im in ims) + 20 * (len(ims) - 1)
    out = Image.new("RGB", (w, h), "white")
    x = 0
    for im in ims:
        out.paste(im, (x, (h - im.height) // 2))
        x += im.width + 20
    out.save(MEDIA / f"fig_engine_{fig:02d}.png")
    print(f"fig_engine_{fig:02d}.png  <- {d} pages {pages}  {out.size}")
