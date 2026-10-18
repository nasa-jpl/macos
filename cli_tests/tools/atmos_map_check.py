#!/usr/bin/env python3
"""The ATMOS gate's arithmetic (PLAN_CONSOLIDATION 1.2 (p), 2026-10-18), pure Python so a journal can
run it with the system python3.  Reads the CLI's TEXT dumps and prints one ATMCHK line per check and
a verdict line `ATMCHK PASS` / `ATMCHK FAIL`, which the journal's EXPECT reads.

  atmos_map_check.py [--screen OPD.txt] [--ep OPD.txt] [--spot spot.txt] [--screens a,b,c]
                     --dx 0.08 --r0 0.25 --lam 1e-6 [--efl 22.4736] [--diam 4.0]

  --screen   the OPD map AT the screen element (= the screen when it sits at a flat pupil)
  --ep       the OPD map at the exit-pupil Reference downstream
  --spot     the focal-plane spot (x, y per ray, the rays in the order of the map's valid pixels)
  --screens  several seeds' screen maps, comma-separated, for the ensemble check
  --dir      the directory the file names are relative to (a journal's $TMP: the CLI's `!` line
             runs a leading `cd` itself and drops the rest of the line)

Checks (each runs when its inputs are given):
  kolm   (--screens) the structure function D(L) of each piston-removed screen (its tilt kept: the
         Kolmogorov D includes it) at lags 1..8 px, x and y pooled, over Kolmogorov 6.88 (L dx/r0)^(5/3),
         AVERAGED over the seeds: every lag within a factor 1.5.  One 51-px realization is not a test
         at the larger lags (24 seeds: 0.32-2.0 at L = 8, mean 0.98-1.06 at every lag)
  slope  the log-log slope of the tilt-removed screen's D between lags 1 and 4 within 5/3 +- 0.3 (24
         seeds: 1.43-1.71); a pixel checkerboard or white noise gives ~0
  noll   the tilt-removed rms over the pupil within a factor 1.5 of Noll's sqrt(0.134 (D/r0)^(5/3))
  carry  the exit-pupil map equals the screen (tilt removed): correlation >= 0.999, rms ratio 1 +- 1 %
  tilt   the spot against the engine's own 3-row (Prewitt) gradient of the exit-pupil map, interior
         pixels, per axis: spot = C x grad, C within 10 % of the EFL (sign included), R^2 >= 0.99
         (the TEXT map's lines are the grid's columns: spot x follows the line index's neighbour,
         spot y the line's)
The map's valid pixels are its nonzero ones (a random screen has no exact zeros inside the pupil).
"""
import math, os, sys


def arg(name, default=None):
    a = sys.argv
    return a[a.index(name) + 1] if name in a else default


def load(path):
    path = os.path.join(arg('--dir', ''), path)
    with open(path) as f:
        return [[float(t) for t in ln.split()] for ln in f if ln.strip()]


def solve3(A, b):
    """Gaussian elimination, 3x3."""
    M = [A[i][:] + [b[i]] for i in range(3)]
    for c in range(3):
        p = max(range(c, 3), key=lambda r: abs(M[r][c]))
        M[c], M[p] = M[p], M[c]
        for r in range(3):
            if r != c:
                f = M[r][c] / M[c][c]
                M[r] = [M[r][k] - f * M[c][k] for k in range(4)]
    return [M[i][3] / M[i][i] for i in range(3)]


def detilt(a, k):
    """Piston + tilt removed over the valid pixels; returns (map, mask), map scaled by k."""
    n = len(a); c = (n - 1) / 2.0
    pts = [(i, j) for i in range(n) for j in range(n) if a[i][j] != 0.0]
    A = [[0.0] * 3 for _ in range(3)]; b = [0.0] * 3
    for i, j in pts:
        v = (1.0, i - c, j - c); w = a[i][j]
        for r in range(3):
            b[r] += v[r] * w
            for s in range(3):
                A[r][s] += v[r] * v[s]
    x = solve3(A, b)
    out = [[0.0] * n for _ in range(n)]; m = [[False] * n for _ in range(n)]
    for i, j in pts:
        out[i][j] = k * (a[i][j] - x[0] - x[1] * (i - c) - x[2] * (j - c)); m[i][j] = True
    return out, m


def rms(a, m):
    v = [a[i][j] for i in range(len(a)) for j in range(len(a)) if m[i][j]]
    mu = sum(v) / len(v)
    return math.sqrt(sum((t - mu) ** 2 for t in v) / len(v)), v


dx = float(arg('--dx')); r0 = float(arg('--r0')); lam = float(arg('--lam'))
efl = arg('--efl'); diam = float(arg('--diam', '4.0'))
k = 2 * math.pi / lam
ok = True


def report(name, good, text):
    global ok
    ok = ok and good
    print('ATMCHK %-5s %s  %s' % (name, 'ok ' if good else 'BAD', text))


def dratio(S, ms, lags):
    n = len(S); rat = []
    for L in lags:
        d = [(S[i + L][j] - S[i][j]) ** 2 for i in range(n - L) for j in range(n) if ms[i + L][j] and ms[i][j]]
        d += [(S[i][j + L] - S[i][j]) ** 2 for i in range(n) for j in range(n - L) if ms[i][j + L] and ms[i][j]]
        rat.append(sum(d) / len(d) / (6.88 * (L * dx / r0) ** (5.0 / 3.0)))
    return rat


def depiston(a, k):
    n = len(a); v = [a[i][j] for i in range(n) for j in range(n) if a[i][j] != 0.0]; mu = sum(v) / len(v)
    return ([[k * (a[i][j] - mu) if a[i][j] != 0.0 else 0.0 for j in range(n)] for i in range(n)],
            [[a[i][j] != 0.0 for j in range(n)] for i in range(n)])


scr = arg('--screen'); ep = arg('--ep'); spot = arg('--spot'); scrs = arg('--screens')
if scrs:
    all_r = [dratio(*depiston(load(f), k), range(1, 9)) for f in scrs.split(',')]
    mean = [sum(r[L] for r in all_r) / len(all_r) for L in range(8)]
    report('kolm', all(1 / 1.5 <= r <= 1.5 for r in mean),
           'D(L)/Kolmogorov, L = 1..8, mean of %d seeds: ' % len(all_r) + ' '.join('%.3f' % r for r in mean))
if scr:
    S, ms = detilt(load(scr), k); n = len(S)
    r1, r4 = dratio(S, ms, (1, 4)); sl = math.log(r4 / r1 * 4 ** (5.0 / 3.0)) / math.log(4.0)
    report('slope', abs(sl - 5.0 / 3.0) <= 0.3, 'D(4)/D(1) log-slope %.3f (Kolmogorov 1.667); D(1)/Kolmogorov %.3f' % (sl, r1))
    rs, _ = rms(S, ms); noll = math.sqrt(0.134 * (diam / r0) ** (5.0 / 3.0))
    report('noll', 1 / 1.5 <= rs / noll <= 1.5, 'tilt-removed rms %.3f rad, Noll %.3f rad, ratio %.3f' % (rs, noll, rs / noll))
if scr and ep:
    E, me = detilt(load(ep), k)
    both = [(i, j) for i in range(n) for j in range(n) if ms[i][j] and me[i][j]]
    x = [S[i][j] for i, j in both]; y = [E[i][j] for i, j in both]
    mx = sum(x) / len(x); my = sum(y) / len(y)
    sxx = sum((t - mx) ** 2 for t in x); syy = sum((t - my) ** 2 for t in y)
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    cor = sxy / math.sqrt(sxx * syy); ratio = math.sqrt(syy / sxx)
    report('carry', cor >= 0.999 and abs(ratio - 1) <= 0.01, 'exit pupil vs screen: correlation %.6f, rms ratio %.5f' % (cor, ratio))
if ep and spot and efl:
    efl = float(efl)
    a = load(ep); n = len(a); sp = load(spot)
    valid = [(i, j) for i in range(n) for j in range(n) if a[i][j] != 0.0]
    if len(valid) != len(sp):
        report('tilt', False, 'map has %d valid pixels, spot %d rays' % (len(valid), len(sp)))
    else:
        pos = {p: sp[q] for q, p in enumerate(valid)}
        val = lambda i, j: a[i][j]
        rows = []
        for i, j in valid:
            if not all((i + di, j + dj) in pos for di in (-1, 0, 1) for dj in (-1, 0, 1)):
                continue
            g0 = sum(val(i - 1, j + t) - val(i + 1, j + t) for t in (-1, 0, 1)) / (6 * dx)
            g1 = sum(val(i + t, j - 1) - val(i + t, j + 1) for t in (-1, 0, 1)) / (6 * dx)
            rows.append((g0, g1, pos[(i, j)][0], pos[(i, j)][1]))
        res = []
        for gi, si in ((1, 2), (0, 3)):        # spot x ~ grad along the line, spot y ~ across lines
            g = [r[gi] for r in rows]; s = [r[si] for r in rows]
            mg = sum(g) / len(g); ms_ = sum(s) / len(s)
            g = [t - mg for t in g]; s = [t - ms_ for t in s]
            C = sum(p * q for p, q in zip(g, s)) / sum(p * p for p in g)
            r2 = 1 - sum((q - C * p) ** 2 for p, q in zip(g, s)) / sum(q * q for q in s)
            res.append((C, r2))
        good = all(0.9 <= C / efl <= 1.1 and r2 >= 0.99 for C, r2 in res)
        report('tilt', good, 'spot = C x grad W: C/EFL x %.4f y %.4f, R^2 %.5f %.5f (%d interior rays)'
               % (res[0][0] / efl, res[1][0] / efl, res[0][1], res[1][1], len(rows)))
print('ATMCHK ' + ('PASS' if ok else 'FAIL'))
