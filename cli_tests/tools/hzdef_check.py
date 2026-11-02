#!/usr/bin/env python3
"""Validate macos_f90/hzdef.inc -- the NormHex (ZernType 7) polynomials ZerntoMon4 converts -- against an
independent construction (BRIEF_to_zernike_repairs.md (uu), Dave 2026-10-20).

The construction is the textbook one (Mahajan & Dai, "Orthonormal polynomials in wavefront analysis:
analytical solution", JOSA A 24, 2994 (2007), hexagonal pupil): the circle Zernikes in Noll's order,
Gram-Schmidt-orthonormalized over the unit regular hexagon (circumradius 1) with the area-mean inner
product, so each mode has RMS 1 over the hexagon.  Nothing here is read from the engine but the table under
test.  The moments are exact (Gauss-Legendre on the six centre triangles, degree 16 << 2 * 16 - 1).  Both
orientations are built: the table matches one of them, which is the orientation the manual states.
Independent of the construction, eight closed forms printed by Mahajan & Dai (modes 4-11) are checked.

    hzdef_check.py [hzdef.inc] [-v]        exit 1 unless every mode matches one orientation to 1e-9 rel

Also imported by docs/macos-manual/tools/gen_zernike_table.py (the analytic its NormHex rows are measured
against) -- the hexagon basis is written once.
"""
import math
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

NMODE = 45
ROOT = Path(__file__).resolve().parents[2]


# ---- polynomials as {(a, b): c}, the monomial x^a y^b -----------------------------------------------
def pmul(p, q):
    r = defaultdict(float)
    for (a, b), c in p.items():
        for (d, e), f in q.items():
            r[(a + d, b + e)] += c * f
    return dict(r)


def padd(p, q, s=1.0):
    r = defaultdict(float, p)
    for k, v in q.items():
        r[k] += s * v
    return dict(r)


def ppow(p, n):
    r = {(0, 0): 1.0}
    for _ in range(n):
        r = pmul(r, p)
    return r


def peval(p, x, y):
    return sum(c * x ** a * y ** b for (a, b), c in p.items())


RHO2 = {(2, 0): 1.0, (0, 2): 1.0}


def cis(m):
    """Re and Im of (x + i y)^m."""
    re_, im_ = defaultdict(float), defaultdict(float)
    for k in range(m + 1):
        c = math.comb(m, k)            # x^(m-k) (i y)^k
        [re_, im_, re_, im_][k % 4][(m - k, k)] += c if k % 4 < 2 else -c
    return dict(re_), dict(im_)


def noll_nm(j):
    """Noll index -> (n, m): within an order |m| increases; of a pair, even j is cos (m > 0)."""
    n = 0
    while (n + 1) * (n + 2) // 2 < j:
        n += 1
    k = j - n * (n + 1) // 2
    ms = [m for m in range(n % 2, n + 1, 2) for _ in ((0,) if m == 0 else (0, 1))]
    am = ms[k - 1]
    return (n, 0) if am == 0 else (n, am if j % 2 == 0 else -am)


def zernike(n, m):
    """Circle Zernike R_n^|m|(rho) cos(m A) (m >= 0) or sin(|m| A), unnormalized, in x, y."""
    am = abs(m)
    re_, im_ = cis(am)
    p = {}
    for s in range((n - am) // 2 + 1):
        c = (-1) ** s * math.factorial(n - s) / (
            math.factorial(s) * math.factorial((n + am) // 2 - s) * math.factorial((n - am) // 2 - s))
        p = padd(p, pmul(ppow(RHO2, (n - am) // 2 - s), re_ if m >= 0 else im_), c)
    return {k: v for k, v in p.items() if v != 0}


def moments(vertex_on_x=True, maxdeg=18, ng=16):
    """Area means of x^a y^b over the unit hexagon (circumradius 1)."""
    th0 = 0.0 if vertex_on_x else math.pi / 2
    V = [(math.cos(th0 + k * math.pi / 3), math.sin(th0 + k * math.pi / 3)) for k in range(6)]
    g, w = np.polynomial.legendre.leggauss(ng)
    u, wu = (g + 1) / 2, w / 2
    X, Y, W = [], [], []
    for k in range(6):
        (x1, y1), (x2, y2) = V[k], V[(k + 1) % 6]
        J = abs(x1 * y2 - x2 * y1)       # triangle (0, P1, P2) on the unit square: a((1-b)P1 + bP2)
        for a, wa in zip(u, wu):
            for b, wb in zip(u, wu):
                X.append(a * ((1 - b) * x1 + b * x2))
                Y.append(a * ((1 - b) * y1 + b * y2))
                W.append(wa * wb * a * J)
    X, Y, W = map(np.array, (X, Y, W))
    area = W.sum()
    assert abs(area - 3 * math.sqrt(3) / 2) < 1e-12, area
    return {(i, j): float((W * X ** i * Y ** j).sum() / area)
            for i in range(maxdeg + 1) for j in range(maxdeg + 1 - i)}


def inner(p, q, M):
    return sum(c * f * M[(a + d, b + e)] for (a, b), c in p.items() for (d, e), f in q.items())


def hexagon_basis(vertex_on_x=True):
    """The NMODE hexagon-orthonormal polynomials, Noll order (Gram-Schmidt, twice for round-off)."""
    M = moments(vertex_on_x)
    out = []
    for j in range(1, NMODE + 1):
        v = zernike(*noll_nm(j))
        for _ in range(2):
            for u in out:
                v = padd(v, u, -inner(v, u, M))
        v = {k: c / math.sqrt(inner(v, v, M)) for k, c in v.items()}
        out.append({k: c for k, c in v.items() if abs(c) > 1e-13})
    return out


# Mahajan & Dai (2007), hexagon: closed forms of modes 4-11, rho and A polar in the hexagon's frame
def _md(j, rho, A):
    s5 = math.sqrt
    return {4: s5(5 / 43) * (12 * rho ** 2 - 5),
            5: 2 * s5(15 / 7) * rho ** 2 * np.sin(2 * A),
            6: 2 * s5(15 / 7) * rho ** 2 * np.cos(2 * A),
            7: 4 * s5(42 / 3685) * (25 * rho ** 3 - 14 * rho) * np.sin(A),
            8: 4 * s5(42 / 3685) * (25 * rho ** 3 - 14 * rho) * np.cos(A),
            9: (4 / 3) * s5(10) * rho ** 3 * np.sin(3 * A),
            10: 4 * s5(70 / 103) * rho ** 3 * np.cos(3 * A),
            11: 3 / s5(1072205) * (6020 * rho ** 4 - 5140 * rho ** 2 + 737)}[j]


# ---- the table under test --------------------------------------------------------------------------
def read_table(hz_path):
    t = Path(hz_path).read_text()
    hzi = defaultdict(dict)
    for i, j, c, v in re.findall(r'hzi\((\d+),(\d+),(\d)\)\s*=\s*(\d+)', t):
        hzi[int(i)].setdefault(int(j), [None, None])[int(c) - 1] = int(v)
    hz = {}
    for blk in re.split(r'\n\s*iz\s*=\s*', t)[1:]:
        iz = int(re.match(r'(\d+)', blk).group(1))
        # one statement (continuations joined); blanks are not significant in fixed form, so the
        # compiler reads '38.71 8068' as one literal -- parse it the same way
        body = re.sub(r'\n\s*&', ' ', blk.split('=', 1)[1]).split('\n')[0]
        body = re.sub(r'[ \t]', '', body)
        lits = re.findall(r'[-+]?\d+\.\d*(?:[dDeE][-+]?\d+)?', body)
        assert ''.join(lits) == re.sub(r'[(/),;]', '', body), f'mode {iz}: unparsed text in {body!r}'
        hz[iz] = [float(x.replace('d', 'e').replace('D', 'e')) for x in lits]
        assert all(math.isfinite(v) and abs(v) < 1e6 for v in hz[iz]), f'mode {iz}: {hz[iz]}'
        # the literals must be double precision: without a d exponent Fortran reads REAL(4)
        assert all(re.search(r'[dD][-+]?\d+$', x) for x in lits), f'mode {iz}: a literal is REAL(4)'
    src = (Path(hz_path).parent / 'surfsub.F').read_text()
    g = re.search(r'Subroutine ZerntoMon4.*?nterm=\(/(.*?)/\)', src, re.S)
    nterm = [int(x) for x in re.sub(r'\n\s*&', ' ', g.group(1)).replace(' ', '').split(',')]
    polys = {}
    for i in range(1, NMODE + 1):
        ex = [tuple(hzi[i][j]) for j in sorted(hzi[i])]
        assert len(ex) == nterm[i - 1], f'mode {i}: {len(ex)} exponent pairs, nterm {nterm[i - 1]}'
        c = hz[i] * len(ex) if len(hz[i]) == 1 else hz[i]       # hz(1:n,iz) = scalar broadcasts
        assert len(c) == len(ex), f'mode {i}: {len(c)} coefficients for {len(ex)} terms'
        p = defaultdict(float)
        for e, v in zip(ex, c):
            p[e] += v
        polys[i] = dict(p)
    return polys


def main():
    args = [a for a in sys.argv[1:] if not a.startswith('-')]
    hz_path = Path(args[0]) if args else ROOT / 'macos_f90' / 'hzdef.inc'
    tab = read_table(hz_path)
    print(f'{hz_path}: {len(tab)} modes; exponent / coefficient counts = ZerntoMon4 nterm')
    rng = np.random.default_rng(0)
    rho, A = np.sqrt(rng.uniform(0, 1, 400)), rng.uniform(-np.pi, np.pi, 400)
    match = {}
    for orient, name in ((True, 'vertex on x'), (False, 'vertex on y')):
        B = hexagon_basis(orient)
        rows = []
        for j in range(1, NMODE + 1):
            p, q = tab[j], B[j - 1]
            d = max(abs(p.get(k, 0.0) - q.get(k, 0.0)) for k in set(p) | set(q))
            rows.append((j, noll_nm(j), d, d / max(abs(v) for v in q.values())))
        worst = max(rows, key=lambda r: r[3])
        nbad = sum(r[3] > 1e-9 for r in rows)
        match[name] = nbad == 0
        print(f'== {name}: worst mode {worst[0]} (n,m)={worst[1]}: largest coefficient difference '
              f'{worst[2]:.2e} ({worst[3]:.1e} of the mode\'s largest coefficient); {nbad} of {NMODE} modes '
              f'differ by more than 1e-9 of it')
        if '-v' in sys.argv:
            for j, nm, d, r in rows:
                print(f'   mode {j:2d}  Noll (n,m) = {str(nm):8s}  |table - construction| <= {d:.2e}  ({r:.1e})')
        if orient:                       # the closed forms, in the orientation the table uses
            x, y = rho * np.cos(A), rho * np.sin(A)
            md = max(float(np.max(np.abs(peval(tab[j], x, y) - _md(j, rho, A)))) for j in range(4, 12))
            print(f'   Mahajan & Dai closed forms, modes 4-11: table - printed form <= {md:.2e}')
            match['md'] = md < 1e-9
    ok = (match['vertex on x'] != match['vertex on y']) and match['md']
    print('PASS: hzdef.inc = the hexagon-orthonormal Noll set, ' +
          ('vertices on +-x' if match['vertex on x'] else 'vertices on +-y') if ok else 'FAIL')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
