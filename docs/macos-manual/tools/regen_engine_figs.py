#!/usr/bin/env python3
"""Regenerate the manual figures that come from the running engine (2026-10-17):
figure 31 (DRAW views of ImageDemo.in), 56 (GBLUR on CoroExample.in), 58 (Luneberg gain pattern).
Each figure has a journal in tools/figs/ that the CLI suite's driver (cli_tests/cli_cmd_gate.py) runs
on the /png device; the journal copies its giza pages into build/figs/figNN/, and this script
assembles media/fig_rescued_NN.png from them.  A journal whose EXPECT fails stops the figure, so a
figure cannot silently rot into a wrong picture.
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
def page(fig, n):
    return Image.open(FIGS / f"fig{fig}" / f"giza_{n:04d}.png")
def pair(a, b, dest):
    out = Image.new("RGB", (a.width + b.width, max(a.height, b.height)), "white")
    out.paste(a, (0, 0)); out.paste(b, (a.width, 0)); out.save(dest)
pair(page(31, 0), page(31, 1), MEDIA / "fig_rescued_31.png")
pair(page(56, 1), page(56, 2), MEDIA / "fig_rescued_56.png")
page(58, 2).save(MEDIA / "fig_rescued_58.png")
print("wrote fig_rescued_31.png (two DRAW views), fig_rescued_56.png (two GBLUR composites), fig_rescued_58.png (gain column)")
