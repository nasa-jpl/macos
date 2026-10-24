#!/usr/bin/env python3
"""Tile a deck's engine renders (3-D view + dispersion-plane view) side by side
for deck_dyson, trimming white margins.  Composition only -- the panels are
macos.view_rx output (challenges/dyson5/dyson5_view_figs.m), never redrawn.
Usage: tile_views.py <tag> ...   (reads ../../MACOS_resources/mmacos/challenges/dyson5/<tag>_view{3d,yz}.png,
                                 writes figs_dyson/<tag>_views.png)"""
import os, sys
from PIL import Image, ImageChops
SRC = os.path.expanduser('~/dev/MACOS_resources/mmacos/challenges/dyson5')
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'figs_dyson')
def trim(im, pad=12):
    bg = Image.new(im.mode, im.size, (255, 255, 255))
    bbox = ImageChops.difference(im, bg).getbbox()
    if not bbox: return im
    l, t, r, b = bbox
    return im.crop((max(0, l-pad), max(0, t-pad), min(im.width, r+pad), min(im.height, b+pad)))
for tag in sys.argv[1:]:
    panels = []
    for suf in ('view3d', 'viewyz'):
        p = os.path.join(SRC, f'{tag}_{suf}.png')
        if not os.path.exists(p): sys.exit(f'MISSING {p}')
        panels.append(trim(Image.open(p).convert('RGB')))
    h = max(im.height for im in panels)
    panels = [im.resize((round(im.width * h / im.height), h)) for im in panels]
    gutter = 30
    W = sum(im.width for im in panels) + gutter * (len(panels) - 1)
    out = Image.new('RGB', (W, h), (255, 255, 255)); x = 0
    for im in panels: out.paste(im, (x, 0)); x += im.width + gutter
    os.makedirs(OUT, exist_ok=True); out.save(os.path.join(OUT, f'{tag}_views.png'))
    print('wrote', f'{tag}_views.png', out.size)
