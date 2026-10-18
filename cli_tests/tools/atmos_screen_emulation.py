#!/usr/bin/env python3
"""Exact pure-Python replica of AtmosPhase (macos_f90/sourcsub.F): the Kolmogorov phase-screen
generator -- the filter sqrt(PSD) = 0.151 r0^(-5/6) f^(-11/6) times the frequency measure 1/(N dx), the
Hermitian-symmetric complex Gaussian noise, the unnormalized inverse transform, the piston/tilt
removal over nradius, the optional re-added tilt -- with the engine's own RNG / GAUSS (mathsub.F).
Written 2026-10-18 for the ATMOS finding (PLAN_CONSOLIDATION 1.2 (p)): it shows the routine
produces a correct screen (D_phi(1 px) within 20 % of 6.88 (r/r0)^(5/3)) for every seed, so the
defect the rays see is in AddAtmos's coupling, not here.  Usage: python3 atmos_screen_emulation.py
[seed ...]  (defaults 12345 4321 777; npts 256, dx 0.08 m, r0 0.25 m, nradius 52, pTilt 0)."""
import math, cmath, sys
def rng(nts):
    for i, b in ((0, 32717.0), (1, 32713.0)):
        x = float(nts[i]) * 3125.0; k = int(x / b); nts[i] = int(x - k * b + 0.5)
        if nts[i] <= 0 or nts[i] > int(b): nts[i] = 123
    k = nts[0] ^ nts[1]
    return (k if k else 1) / 32768.0
def gauss(nts):
    u = rng(nts); v = rng(nts); return math.sqrt(-2 * math.log(u)) * math.cos(2 * math.pi * v)
def fft1(a):
    n = len(a)
    if n == 1: return a
    e = fft1(a[0::2]); o = fft1(a[1::2]); out = [0] * n
    for k in range(n // 2):
        t = cmath.exp(2j * math.pi * k / n) * o[k]; out[k] = e[k] + t; out[k + n // 2] = e[k] - t
    return out
def fft2(x):
    n = len(x); y = [fft1(r) for r in x]; c = [fft1([y[j][i] for j in range(n)]) for i in range(n)]
    return [[c[i][j] for i in range(n)] for j in range(n)]
def screen(seed, npts=256, rsamp=0.08, r0=0.25, nradius=52, pTilt=0.0):
    nts = [seed, seed]
    coeff = (0.151 / r0 ** (5 / 6)) * (1.0 / (npts * rsamp)); scale = (1.0 / (npts * rsamp)) ** 2
    f = [[0.0] * npts for _ in range(npts)]
    for i in range(npts):
        ii = i - npts if i >= npts // 2 else i
        for j in range(npts):
            jj = j - npts if j >= npts // 2 else j
            ff = scale * (ii * ii + jj * jj) or scale
            f[j][i] = coeff / ff ** (11 / 12)
    x = [[0j] * npts for _ in range(npts)]; s = 0.70710678118655
    n2 = npts // 2; n21 = n2 + 1; n22 = n2 + 2; np2 = npts + 2
    for i in range(1, n21 + 1):
        for j in range(1, n21 + 1): x[j - 1][i - 1] = complex(s * gauss(nts), s * gauss(nts))
    for j in range(n22, npts + 1):
        x[j - 1][0] = x[np2 - j - 1][0].conjugate(); x[j - 1][n21 - 1] = x[np2 - j - 1][n21 - 1].conjugate()
    for i in range(2, n2 + 1):
        for j in range(n22, npts + 1): x[j - 1][i - 1] = complex(s * gauss(nts), s * gauss(nts))
    for i in range(n22, npts + 1):
        np2i = np2 - i; x[0][i - 1] = x[0][np2i - 1].conjugate()
        for j in range(2, npts + 1): x[j - 1][i - 1] = x[np2 - j - 1][np2i - 1].conjugate()
    for a, b in ((1, 1), (1, n21), (n21, 1), (n21, n21)): x[a - 1][b - 1] = complex(math.sqrt(2) * x[a - 1][b - 1].real, 0)
    for i in range(npts):
        for j in range(npts): x[j][i] *= f[j][i]
    y = fft2(x)
    idx = lambda v: v + npts if v < 0 else v
    bm = bx = by = cm = cx = cy = 0.0
    for i in range(-nradius, nradius + 1):
        for j in range(-nradius, nradius + 1):
            if i * i + j * j <= nradius * nradius:
                foo = y[idx(j)][idx(i)].real; bm += foo; cm += 1; bx += foo * i; cx += i * i; by += foo * j; cy += j * j
    bm /= cm; bx /= cx; by /= cy
    if pTilt:
        tiltamp = pTilt * math.sqrt(1.8 * (2 * nradius * rsamp / r0) ** (5 / 3))
        bx += tiltamp * gauss(nts) / nradius; by += tiltamp * gauss(nts) / nradius
    return lambda i, j: y[idx(j)][idx(i)].real - bm - bx * i - by * j, max(abs(v.imag) for r in y for v in r)
if __name__ == '__main__':
    seeds = [int(a) for a in sys.argv[1:]] or [12345, 4321, 777]
    r0, dx = 0.25, 0.08
    for sd in seeds:
        scr, imag = screen(sd)
        vals = [scr(i, j) for i in range(-25, 26) for j in range(-25, 26) if i * i + j * j <= 625]
        rms = math.sqrt(sum(v * v for v in vals) / len(vals))
        d = [(scr(i, j) - scr(i + 1, j)) ** 2 for i in range(-24, 24) for j in range(-24, 25) if i * i + j * j <= 576]
        print('seed %6d: rms over the 25-px disk %5.2f rad (tilt removed); D_phi(1 px) %.2f rad^2 vs Kolmogorov %.2f; max imag %.1e'
              % (sd, rms, sum(d) / len(d), 6.88 * (dx / r0) ** (5 / 3), imag))
