#!/usr/bin/env python3
"""One-time migration helper: rasterize composite line-art figures.

The FrameMaker-era figures survive in the docx as absolutely-positioned
vector fragments + loose text labels; pandoc flattens those to garbage.
This script crops each such figure out of a correctly-rendered PDF of
the docx (figure region = between the caption line and the nearest
body-text line above it), saves it as src/media/fig_NN.png, inserts a
single image reference before the caption in the Markdown, and removes
the fragment debris (tiny images and short label lines) near the
caption.  Every removal is logged to FIGURE_RESCUE_LOG.md for review.

Usage: python3 tools/rescue_figures.py <rendered.pdf> <bbox.html> FIG [FIG...]
where FIG numbers are the composite figures to rescue.
"""
import re
import subprocess
import sys
from pathlib import Path
from PIL import Image

HERE = Path(__file__).resolve().parent
SRC = HERE.parent / "src"
MEDIA = SRC / "media"
LOG = HERE.parent / "FIGURE_RESCUE_LOG.md"

DPI = 150
PT2PX = DPI / 72.0
PAGE_W_PT = 612


def parse_bbox(bbox_html):
    """Return {page_no: [ (xmin,ymin,xmax,ymax,text), ... ]} lines."""
    pages = {}
    page_no = 0
    words = []
    for ln in Path(bbox_html).read_text().splitlines():
        if "<page " in ln:
            page_no += 1
            words = []
            pages[page_no] = words
        m = re.search(r'<word xMin="([\d.]+)" yMin="([\d.]+)" '
                      r'xMax="([\d.]+)" yMax="([\d.]+)">([^<]*)</word>', ln)
        if m:
            words.append((float(m.group(1)), float(m.group(2)),
                          float(m.group(3)), float(m.group(4)), m.group(5)))
    # group words into lines by y-center proximity
    lines = {}
    for p, ws in pages.items():
        ws = sorted(ws, key=lambda w: (round((w[1] + w[3]) / 2 / 4), w[0]))
        out = []
        cur, cur_yc = [], None
        for w in ws:
            yc = (w[1] + w[3]) / 2
            if cur and abs(yc - cur_yc) > 4:
                out.append(cur)
                cur = []
            cur.append(w)
            cur_yc = yc
        if cur:
            out.append(cur)
        lines[p] = [{
            "x0": min(w[0] for w in l), "x1": max(w[2] for w in l),
            "y0": min(w[1] for w in l), "y1": max(w[3] for w in l),
            "text": " ".join(w[4] for w in l), "nw": len(l),
        } for l in out]
    return lines


def find_caption(lines, fig):
    pat = re.compile(rf"^FIGURE {fig}(?![0-9])")
    for p, ls in lines.items():
        for i, l in enumerate(sorted(ls, key=lambda x: x["y0"])):
            if pat.match(l["text"]):
                return p, l
    return None, None


def is_body(l):
    return l["nw"] >= 8 and (l["x1"] - l["x0"]) > 250


def crop_figure(pdf, page, cap, lines, out_png, all_lines=None):
    """Crop from below the nearest body line above the caption down to
    the caption top.  If the caption sits near the top of its page the
    artwork usually finished the PREVIOUS page — crop that page's tail
    instead."""
    above = [l for l in lines if l["y1"] <= cap["y0"] - 2]
    body_above = [l for l in above if is_body(l)]
    top = max((l["y1"] for l in body_above), default=50) + 4
    bottom = cap["y0"] - 2
    if bottom - top < 80 and cap["y0"] < 160 and all_lines \
            and (page - 1) in all_lines:
        prev = sorted(all_lines[page - 1], key=lambda l: l["y0"])
        body_prev = [l for l in prev if is_body(l)]
        top = max((l["y1"] for l in body_prev), default=50) + 4
        bottom = 730            # above the footer/page number
        page = page - 1
    elif bottom - top < 20:    # nothing visible? take a taller slice
        top = max(40, bottom - 250)
    prefix = out_png.parent / (out_png.stem + "_page")
    subprocess.run(["pdftoppm", "-png", "-r", str(DPI), "-f", str(page),
                    "-l", str(page), "-singlefile", str(pdf),
                    str(prefix)], check=True)
    page_png = prefix.with_suffix(".png")
    img = Image.open(page_png)
    x0 = int(60 * PT2PX)
    x1 = int((PAGE_W_PT - 40) * PT2PX)
    img.crop((x0, int(top * PT2PX), x1, int(bottom * PT2PX))).save(out_png)
    page_png.unlink()
    return top, bottom


TINY_IMG = re.compile(
    r'!\[\]\(media/image\d+\.[a-z]+\)\{width="([0-9.e-]+)in" '
    r'height="([0-9.e-]+)in"\}')


def patch_markdown(fig, png_name, log):
    """Insert the rescued image before the caption; strip fragment
    debris in the 8 lines above (tiny images, short label-only lines)."""
    cap_re = re.compile(rf"^\*\*FIGURE {fig}\*\*")
    for f in sorted(SRC.glob("[0-9]*.md")):
        txt = f.read_text().splitlines()
        for i, ln in enumerate(txt):
            if not cap_re.match(ln):
                continue
            lo = max(0, i - 8)
            removed = []
            for j in range(lo, i):
                ln_j = txt[j]
                if ln_j is None or not ln_j.strip():
                    continue
                stripped = TINY_IMG.sub(
                    lambda m: "" if float(m.group(1)) < 0.6
                    and float(m.group(2)) < 0.6 else m.group(0), ln_j)
                if stripped != ln_j:
                    if not stripped.strip():
                        txt[j] = None
                        removed.append(f"dropped: {ln_j[:60]}")
                        continue
                    ln_j = stripped   # fall through to the label check
                    txt[j] = stripped
                    removed.append(f"trimmed imgs: {ln_j[:60]}...")
                # short label-only line: <=4 words, no sentence period
                words = ln_j.split()
                if (0 < len(words) <= 4 and not ln_j.rstrip().endswith(".")
                        and not ln_j.lstrip().startswith(
                            ("#", "-", "|", "+", "```", ">", "!"))):
                    txt[j] = None
                    removed.append(f"dropped label: {ln_j.strip()!r}")
            txt[i] = f"![](media/{png_name})\n\n" + ln
            f.write_text("\n".join(l for l in txt if l is not None) + "\n")
            log.append((fig, f.name, removed))
            return True
    return False


def main():
    argv = [a for a in sys.argv[1:] if a != "--redo"]
    redo = "--redo" in sys.argv   # regenerate PNGs only; md untouched
    pdf, bbox = argv[0], argv[1]
    figs = [int(a) for a in argv[2:]]
    lines = parse_bbox(bbox)
    log = []
    for fig in figs:
        page, cap = find_caption(lines, fig)
        if not page:
            print(f"FIGURE {fig}: caption not found in PDF — SKIPPED")
            continue
        png = MEDIA / f"fig_rescued_{fig:02d}.png"
        top, bot = crop_figure(pdf, page, cap,
                               sorted(lines[page], key=lambda l: l["y0"]),
                               png, all_lines=lines)
        ok = False if redo else patch_markdown(fig, png.name, log)
        print(f"FIGURE {fig}: page {page} crop y=({top:.0f},{bot:.0f})pt "
              f"-> {png.name}  md-patched={ok}")
    with open(LOG, "a") as fh:
        for fig, fname, removed in log:
            fh.write(f"\n## FIGURE {fig} ({fname})\n")
            for r in removed:
                fh.write(f"- {r}\n")


if __name__ == "__main__":
    main()
