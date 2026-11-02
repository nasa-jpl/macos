#!/usr/bin/env python3
"""Element table of a MACOS prescription for a deck slide (markdown).
Usage: rx_elt_table.py <file.in>  -> markdown table on stdout
Columns: E#, name, element type, surface, radius R (mm, |KrElt|; flat = -),
conic (KcElt), aperture radius (mm, ApVec(1) when ApType is Circular; '-' when none),
z of the vertex (mm).  Read from the .in as written; nothing is computed."""
import re, sys
txt = open(sys.argv[1]).read()
blocks = re.split(r'\n\s*iElt=\s*(\d+)\s*\n', txt)
rows = []
for k in range(1, len(blocks), 2):
    i, b = int(blocks[k]), blocks[k+1]
    g = lambda key, default='': (re.search(r'^\s*' + key + r'=\s*(.+?)\s*$', b, re.M) or [None, default])[1]
    name = g('EltName', f'E{i}'); typ = g('Element'); srf = g('Surface')
    kr = float(g('KrElt', 'nan').replace('D', 'E')); kc = float(g('KcElt', '0').replace('D', 'E'))
    R = '-' if abs(kr) > 1e10 else f'{abs(kr)*1e3:.1f}'
    ap = '-'
    if g('ApType').strip().lower().startswith('circ'):
        ap = f'{float(g("ApVec").split()[0].replace("D","E"))*1e3:.1f}'
    vpt = g('VptElt').split(); z = f'{float(vpt[2].replace("D","E"))*1e3:.1f}' if len(vpt) == 3 else '-'
    rows.append((i, name, typ, srf, R, f'{kc:g}' if R != '-' else '-', ap, z))
print('| E | name | type | surface | R (mm) | K | aperture r (mm) | z (mm) |')
print('|---|---|---|---|---|---|---|---|')
for r in rows: print('| ' + ' | '.join(str(x) for x in r) + ' |')
