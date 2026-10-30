#!/usr/bin/env python3
"""Appendix A against the 3.2 manual: a key-aware diff of every example deck, the 3.2 PDF listing
(pdftotext -layout of docs/macosMan3.2.pdf, a margin-split exponent rejoined) against
examples/<file> -- (element, key, occurrence) -> values; a value is 'same' when it rounds to the
digits 3.2 printed.  The editorial input of TRANSCRIPT_DIFF.md (BRIEF_to_manual_regen step 4).

Usage (from docs/macos-manual):  python3 tools/diff_appendix_32.py
"""
import re, subprocess, sys
from pathlib import Path
MAN = Path(__file__).resolve().parent.parent
pdf = subprocess.run(["pdftotext", "-layout", str(MAN.parent / "macosMan3.2.pdf"), "-"], capture_output=True,
                     text=True, check=True).stdout.splitlines()
ex = MAN / "examples"
start = [i for i, l in enumerate(pdf) if re.match(r'\s*APPENDIX A Examples', l)][0]
heads = [(i, re.match(r'\s*A\.\d\.\d (\S+)', l).group(1)) for i, l in enumerate(pdf) if i > start and re.match(r'\s*A\.\d\.\d \S+\s*$', l)] + [(len(pdf), None)]
def clean(lines):
    t = '\n'.join(l for l in lines if 'Modeling and Analysis for Controlled Optical Systems' not in l and not re.match(r'\s*A\.\d\s', l))
    return re.sub(r'([dDeE][-+])\n\s*(\d)', r'\1\2', t)
def parse(t):
    t = '\n'.join(l.split('%')[0] for l in t.splitlines())
    toks = t.split(); d = {}; elt = 0; key = None; cnt = {}
    for tok in toks:
        m = re.match(r'^(\w+)=(.*)$', tok)
        if m:
            key = m.group(1)
            if key == 'iElt': elt = None
            k = (elt, key); cnt[k] = cnt.get(k, 0) + 1
            cur = d.setdefault((elt, key, cnt[k]), [])
            if m.group(2): cur.append(m.group(2))
        elif key:
            cur.append(tok)
            if key == 'iElt' and elt is None:
                elt = int(float(tok)); d.pop((None, 'iElt', cnt[(None, 'iElt')]), None)
    return d
def num(s):
    try: return float(s.replace('d', 'e').replace('D', 'E'))
    except ValueError: return None
def sig(s):
    m = re.match(r'[-+]?(\d*)\.?(\d*)', s.replace('d', 'e').replace('D', 'E'))
    return max(1, len((m.group(1) + m.group(2)).lstrip('0')))
def same(a, b):
    x, y = num(a), num(b)
    if x is None or y is None: return a.lower() == b.lower()
    if x == y: return True
    n = sig(a); return float(f'%.{n-1}e' % y) == float(f'%.{n-1}e' % x)
for (i, name), (j, _) in zip(heads, heads[1:]):
    if not name.endswith('.in'): continue
    fn = {'coroExample.in': 'CoroExample.in'}.get(name, name)
    a = parse(clean(pdf[i+1:j])); b = parse(open(f'{ex}/{fn}').read())
    ch = [k for k in a if k in b and (len(a[k]) != len(b[k]) or not all(same(x, y) for x, y in zip(a[k], b[k])))]
    add = [k for k in b if k not in a]; rem = [k for k in a if k not in b]
    prec = [k for k in a if k in b and k not in ch and a[k] != b[k]]
    print(f'## {fn}: {len(prec)} keys re-printed at more digits (same value), {len(ch)} changed, {len(add)} added, {len(rem)} removed')
    for k in ch: print(f'   CHANGED elt {k[0]} {k[1]}: {" ".join(a[k])[:70]}  ->  {" ".join(b[k])[:70]}')
    from collections import Counter
    if add: print('   ADDED  ', ', '.join(f'{kk}x{n}' for kk, n in Counter(f'{k[1]}={" ".join(b[k])[:24]}' for k in add).most_common()))
    if rem: print('   REMOVED', ', '.join(f'{kk}x{n}' for kk, n in Counter(f'{k[1]}={" ".join(a[k])[:24]}' for k in rem).most_common()))
