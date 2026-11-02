#!/usr/bin/env python3
"""Tile two existing renders side by side (trim white margins, equal heights).
Composition only -- the panels are the tools' own output, never redrawn.
Usage: tile_pair.py <out basename> <png A> <png B>   (writes figs_dyson/<out>.png)"""
import os, sys
from PIL import Image, ImageChops
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'figs_dyson')
def trim(im, pad=12):
    bg = Image.new(im.mode, im.size, (255, 255, 255))
    bbox = ImageChops.difference(im, bg).getbbox()
    if not bbox: return im
    l, t, r, b = bbox
    return im.crop((max(0, l-pad), max(0, t-pad), min(im.width, r+pad), min(im.height, b+pad)))
out_name, srcs = sys.argv[1], sys.argv[2:]
panels = []
for p in srcs:
    if not os.path.exists(p): sys.exit(f'MISSING {p}')
    panels.append(trim(Image.open(os.path.expanduser(p)).convert('RGB')))
h = max(im.height for im in panels)
panels = [im.resize((round(im.width * h / im.height), h)) for im in panels]
gutter = 30
W = sum(im.width for im in panels) + gutter * (len(panels) - 1)
out = Image.new('RGB', (W, h), (255, 255, 255)); x = 0
for im in panels: out.paste(im, (x, 0)); x += im.width + gutter
os.makedirs(OUT, exist_ok=True); out.save(os.path.join(OUT, f'{out_name}.png'))
print('wrote', f'{out_name}.png', out.size)
