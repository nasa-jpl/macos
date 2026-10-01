#!/usr/bin/env python3
"""Recompose a wide multi-panel figure for a 16:9 slide: trim white margins and,
for an N-panel strip, split it into R rows at column boundaries found from
white gutters.  Composition only (DECK_STYLE: recompose at the panel level,
never regenerate).  Usage: recompose.py <in.png> <out.png> [rows]"""
import sys
import numpy as np
from PIL import Image

def trim(im, pad=10):
    a = np.asarray(im.convert('L')); ink = np.where(a < 250)
    if ink[0].size == 0: return im
    t, b, l, r = ink[0].min(), ink[0].max(), ink[1].min(), ink[1].max()
    return im.crop((max(0, l-pad), max(0, t-pad), min(im.width, r+pad+1), min(im.height, b+pad+1)))

def split_cols(im, rows):
    """Cut the strip into `rows` row-strips at the widest white gutters."""
    a = np.asarray(im.convert('L')); col_ink = (a < 250).sum(axis=0)
    white = col_ink == 0
    # gutters = runs of white columns; pick the (rows-1) widest that are not at the edges
    runs, start = [], None
    for x, w in enumerate(white):
        if w and start is None: start = x
        if not w and start is not None: runs.append((start, x)); start = None
    runs = [(s, e) for s, e in runs if s > 0 and e < im.width]
    # choose cut points nearest to equal thirds among the gutters
    cuts = []
    for k in range(1, rows):
        target = im.width * k / rows
        s, e = min(runs, key=lambda se: abs((se[0]+se[1])/2 - target))
        cuts.append((s+e)//2)
    edges = [0] + sorted(cuts) + [im.width]
    return [im.crop((edges[i], 0, edges[i+1], im.height)) for i in range(rows)]

src, dst = sys.argv[1], sys.argv[2]; rows = int(sys.argv[3]) if len(sys.argv) > 3 else 1
im = trim(Image.open(src).convert('RGB'))
if rows > 1:
    parts = [trim(p) for p in split_cols(im, rows)]
    W = max(p.width for p in parts); gutter = 16
    H = sum(p.height for p in parts) + gutter*(len(parts)-1)
    out = Image.new('RGB', (W, H), (255,255,255)); y = 0
    for p in parts: out.paste(p, (0, y)); y += p.height + gutter
    im = out
im.save(dst); print('wrote', dst, im.size)
