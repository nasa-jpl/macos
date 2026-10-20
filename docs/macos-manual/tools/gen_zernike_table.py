#!/usr/bin/env python3
"""Measure the engine's Zernike orderings and write src/transcripts/zernike_orderings.md (Appendix E;
`make zernike-table`; BRIEF_to_manual_zernike.md, Dave 2026-10-20: "make the mapping between Zernike types").

Nothing here is read off a listing.  For every ZernType the engine parses, a flat Zernike mirror is loaded
through pymacos, a unit coefficient is put in one mode at a time, the full ray grid is traced, and the sag
where each ray meets the mirror is fitted by least squares against every textbook (unnormalized) Zernike
term through radial order NMAX.  A mode that is one term -- one coefficient, nothing else above the noise --
is printed as that (n, m); the fitted coefficient is its normalization (1, or the RMS factor for Norm*).  A
mode that applies no surface, or is not one term, is reported, not printed as a row.

Runs under the pymacos venv (`make zernike-table` sets it up: oneAPI runtime, PYMACOS_SRC, MACOS_HOME).
The gate is MACOS_resources/pymacos/tests/test_zernike_orderings.py, which re-evaluates every printed
polynomial against the engine.  Never hand-edit the output: rerun this.
"""
import math, os, sys, tempfile
from math import factorial as fact
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE.parent / "src" / "transcripts" / "zernike_orderings.md"
OUT8 = HERE.parent / "src" / "transcripts" / "zernike_table8.md"      # Section 4's excerpt (orders 0-4)
NSHORT = 4
sys.path.insert(0, os.environ.get("PYMACOS_SRC", str(Path.home() / "dev/MACOS_resources/pymacos/src")))
import pymacos.macos as m  # noqa: E402
sys.path.insert(0, str(HERE.parents[2] / "cli_tests" / "tools"))
import hzdef_check as hexz  # noqa: E402  the hexagon-orthonormal construction (Mahajan & Dai), textbook

NHEX_PRINT = 15                # NormHex modes printed as polynomials (orders 0-4); the rest by name
HEX_TOL = 1e-9                 # NormHex: |sag - construction| / max|construction| (measured ~4e-12)
NMAX = 14                      # fit basis: every (n, m) through radial order 14 (ExtFringe's 66 reach 14)
COEF = 1e-4                    # the unit coefficient, scaled so the mirror stays flat to first order
TOL = 1e-8                     # a fitted term below TOL (relative to COEF) is noise; residuals run ~1e-11

# (column title, ZernType= value, normalized parent or None).  Malacara is ANSI's legacy name: measured,
# and refused if it differs.
TYPES = [("ANSI", "ANSI", None), ("Born & Wolf", "BornWolf", None), ("Noll", "Noll", None),
         ("Fringe", "Fringe", None), ("Ext. Fringe", "ExtFringe", None)]
NORM = [("NormANSI", "ANSI"), ("NormBornWolf", "BornWolf"), ("NormNoll", "Noll"), ("NormFringe", "Fringe")]
ANNULAR = ("NormAnnularNoll 0.0", "NormAnnularNoll 0.3")

DECK = """\
        ChfRayDir=  0.0E+00  0.0E+00  1.0E+00
        ChfRayPos=  0.0E+00  0.0E+00  -5.0E+00
          zSource=1.0E+22
        BaseUnits=  m
        WaveUnits=  m
           IndRef=1.0E+00
           Extinc=0.0E+00
          Wavelen=1.0E-06
             Flux=1.0E+00
         Aperture=2.0E+00
         Obscratn=0.0E+00
         GridType=  Circular
         nGridpts=  33
            xGrid=  1.0E+00  0.0E+00  0.0E+00
            yGrid=  0.0E+00  1.0E+00  0.0E+00
             nElt=  2
             iElt=  1
          EltName=  Zsurf
          Element=  Reflector
          Surface=  Zernike
         ZernType=  {ztype}
        nZernCoef=  {n}
         ZernCoef=  0.0E+00
            KrElt=-1.0E+22
            KcElt=0.0E+00
           psiElt=  0.0E+00  0.0E+00  -1.0E+00
           VptElt=  0.0E+00  0.0E+00  0.0E+00
           RptElt=  0.0E+00  0.0E+00  0.0E+00
           IndRef=1.0E+00
           Extinc=1.0E+22
             pMon=  0.0E+00  0.0E+00  0.0E+00
             xMon=  1.0E+00  0.0E+00  0.0E+00
             yMon=  0.0E+00  1.0E+00  0.0E+00
             zMon=  0.0E+00  0.0E+00  1.0E+00
             lMon=  1.0E+00
           ApType=  None
         PropType=  Geometric
             zElt=1.0E+22
             iElt=  2
          EltName=  FP
          Element=  FocalPlane
          Surface=  Flat
            KrElt=-1.0E+22
            KcElt=0.0E+00
           psiElt=  0.0E+00  0.0E+00  1.0E+00
           VptElt=  0.0E+00  0.0E+00  -3.0E+00
           RptElt=  0.0E+00  0.0E+00  -3.0E+00
           IndRef=1.0E+00
           Extinc=0.0E+00
           ApType=  None
         PropType=  Geometric
             zElt=1.0E+22
         nOutCord=  5
             Tout=  1.0E+00  0.0E+00  0.0E+00  0.0E+00  0.0E+00  0.0E+00  0.0E+00
                    0.0E+00  1.0E+00  0.0E+00  0.0E+00  0.0E+00  0.0E+00  0.0E+00
                    0.0E+00  0.0E+00  0.0E+00  1.0E+00  0.0E+00  0.0E+00  0.0E+00
                    0.0E+00  0.0E+00  0.0E+00  0.0E+00  1.0E+00  0.0E+00  0.0E+00
                    0.0E+00  0.0E+00  0.0E+00  0.0E+00  0.0E+00  0.0E+00  1.0E+00
"""


def radial_coefs(n, mm):
    """Textbook R_n^|m| as [(coefficient, power)], highest power first."""
    mm = abs(mm)
    return [((-1) ** k * fact(n - k) // (fact(k) * fact((n + mm) // 2 - k) * fact((n - mm) // 2 - k)), n - 2 * k)
            for k in range((n - mm) // 2 + 1)]


def zterm(n, mm, R, A):
    r = sum(c * R ** p for c, p in radial_coefs(n, mm))
    return r * (np.cos(mm * A) if mm >= 0 else np.sin(-mm * A))


def rms_factor(n, mm):
    return math.sqrt((2 if mm else 1) * (n + 1))


def poly_text(n, mm):
    """The polynomial the 3.2 way: (3R^3^ - 2R) sin(A); pandoc superscripts."""
    terms = []
    for i, (c, p) in enumerate(radial_coefs(n, mm)):
        mag = abs(c)
        body = ("" if mag == 1 and p else str(mag)) + ("R" if p else "") + (f"^{p}^" if p > 1 else "")
        terms.append(body if i == 0 else ("- " if c < 0 else "+ ") + body)
    r = " ".join(terms)
    if mm == 0:
        return r
    trig = ("cos" if mm > 0 else "sin") + (f"({abs(mm)}A)" if abs(mm) > 1 else "(A)")
    return (f"({r}) " if len(terms) > 1 else (r + " " if r != "1" else "")) + trig


BASIS = [(n, mm) for n in range(NMAX + 1) for mm in range(-n, n + 1, 2)]


def load(ztype, work):
    """Load the probe deck with the largest coefficient count the type accepts; None if it accepts none."""
    for n in (66, 45, 37):
        p = work / "zprobe.in"
        p.write_text(DECK.format(ztype=ztype, n=n))
        try:
            m.load(str(p))
            return n
        except Exception:
            continue
    return None


def trace_mode(k):
    """(x, y, sag / COEF) at every ray hit on the mirror, a unit coefficient in mode k alone."""
    m.setEltSrfZernCoef(1, [k], [COEF], reset=True)
    m.modify()
    _, nrays, _ = m.trace_rays(2)
    h = m.ray_pos_hist(int(nrays))
    ok = h["ok"][:, 1]
    return h["P"][0, ok, 1], h["P"][1, ok, 1], h["P"][2, ok, 1] / COEF


def measure_hex(work):
    """NormHex: [(k, err_vertex_on_x, err_vertex_on_y)], each the largest |sag - H_k| over max |H_k|."""
    nmodes = load("NormHex", work)
    if nmodes is None:
        return None
    m.ray_hist("on")
    B = {o: hexz.hexagon_basis(o) for o in (True, False)}
    rows = []
    for k in range(1, nmodes + 1):
        x, y, z = trace_mode(k)
        err = []
        for o in (True, False):
            want = hexz.peval(B[o][k - 1], x, y)
            err.append(float(np.max(np.abs(z - want)) / np.max(np.abs(want))))
        rows.append((k, *err))
    return rows


def hex_text(p):
    """A NormHex polynomial in x and y (xMon, yMon over lMon), 6 significant figures, highest degree first."""
    big = max(abs(c) for c in p.values())                # Gram-Schmidt round-off (~1e-13) is not a term
    terms = sorted(((e, c) for e, c in p.items() if abs(c) > 1e-9 * big),
                   key=lambda t: (-(t[0][0] + t[0][1]), -t[0][0]))
    out = []
    for i, ((a, b), c) in enumerate(terms):
        mono = "".join(v + (f"^{e}^" if e > 1 else "") for v, e in (("x", a), ("y", b)) if e)
        mag = f"{abs(c):.6g}"
        body = mag + mono if mono else mag
        out.append(("-" + body if c < 0 else body) if i == 0 else ("- " if c < 0 else "+ ") + body)
    return " ".join(out)


def measure(ztype, work):
    """[(k, kind, n, m, coef)] for every mode k; kind 'term', 'none' (no surface) or 'mixed'."""
    nmodes = load(ztype, work)
    if nmodes is None:
        return None
    m.ray_hist("on")
    rows = []
    for k in range(1, nmodes + 1):
        try:
            x, y, z = trace_mode(k)
        except Exception:
            rows.append((k, "refused", None, None, 0.0))   # the setter refuses the mode (an inert one)
            continue
        R, A = np.hypot(x, y), np.arctan2(y, x)
        M = np.stack([zterm(n, mm, R, A) for n, mm in BASIS], 1)
        c, *_ = np.linalg.lstsq(M, z, rcond=None)
        big = [i for i in range(len(BASIS)) if abs(c[i]) > TOL]
        resid = float(np.max(np.abs(M @ c - z)))
        if not big and np.max(np.abs(z)) < TOL:
            rows.append((k, "none", None, None, 0.0))
        elif len(big) == 1 and resid < TOL:
            n, mm = BASIS[big[0]]
            rows.append((k, "term", n, mm, float(c[big[0]])))
        else:
            rows.append((k, "mixed", None, None, 0.0))
    return rows


def main():
    m.init(128)
    work = Path(tempfile.mkdtemp(prefix="zern_"))
    meas = {}
    for zt in [t[1] for t in TYPES] + ["Malacara"] + [t[0] for t in NORM] + list(ANNULAR):
        meas[zt] = measure(zt, work)
        print(f"measured {zt}: {'refused' if meas[zt] is None else len(meas[zt])} modes", flush=True)
    hexr = measure_hex(work)
    print(f"measured NormHex: {'refused' if hexr is None else len(hexr)} modes", flush=True)

    errors = []
    def same(r, s):
        return r[1:4] == s[1:4] and abs(r[4] - s[4]) < 1e-6

    if len(meas["Malacara"]) != len(meas["ANSI"]) or not all(map(same, meas["Malacara"], meas["ANSI"])):
        errors.append("Malacara does not trace as ANSI")
    index = {}                                   # (n, m) -> {column: k}
    for title, zt, _ in TYPES:
        for k, kind, n, mm, c in meas[zt]:
            if kind == "term":
                if abs(c - 1) > 1e-6:
                    errors.append(f"{zt} mode {k}: coefficient {c:.9g}, not 1")
                index.setdefault((n, mm), {})[title] = k
    dead = {zt: [k for k, kind, *_ in meas[zt] if kind == "none"] for _, zt, _ in TYPES}
    refused = {zt: [k for k, kind, *_ in meas[zt] if kind == "refused"] for _, zt, _ in TYPES}
    for zt, ks_ in dead.items():
        if ks_:
            errors.append(f"{zt} modes {ks_}: accepted but apply no surface (refuse them, PLAN 1.2 (tt))")
    if hexr is not None:
        bad = [k for k, ex, ey in hexr if ex > HEX_TOL]
        if bad:
            errors.append(f"NormHex modes {bad}: not the hexagon-orthonormal polynomial (vertex on xMon)")
    mixed = {zt: [k for k, kind, *_ in meas[zt] if kind == "mixed"] for _, zt, _ in TYPES}
    for zt, ks in mixed.items():
        if ks:
            errors.append(f"{zt} modes {ks}: not a single Zernike term")
    norm_ok = {}
    for nz, parent in NORM:
        bad = []
        for (k, kind, n, mm, c), (_, pkind, pn, pmm, _) in zip(meas[nz], meas[parent]):
            if (kind, n, mm) != (pkind, pn, pmm) or (kind == "term" and abs(c - rms_factor(n, mm)) > 1e-6):
                bad.append(k)
        norm_ok[nz] = (len(meas[nz]), bad)
        if bad:
            errors.append(f"{nz} modes {bad}: not {parent} times the RMS factor")
    if errors:
        # a finding for PLAN_CONSOLIDATION 1.2, not something to print as the engine's behavior
        sys.exit("gen_zernike_table: the measurement does not fit the table's model:\n  " + "\n  ".join(errors))

    cols = [t[0] for t in TYPES]
    lines = ["| n | m | " + " | ".join(cols) + " | Polynomial |",
             # pandoc sizes a wrapped pipe table by its dashes: n and m narrow, the polynomial ~40%
             "|---:|---:|" + "-------:|" * len(cols) + ":" + "-" * 29 + "|"]
    order = sorted(index, key=lambda nm: (nm[0], nm[1]))
    for n, mm in order:
        cells = [str(index[(n, mm)].get(t, "–")) for t in cols]
        lines.append(f"| {n} | {mm} | " + " | ".join(cells) + f" | {poly_text(n, mm)} |")

    nterm = sum(1 for _, zt, _ in TYPES for r in meas[zt] if r[1] == "term")
    npairs = sum(len(v) for v in meas.values() if v) + len(hexr or [])
    nsettings = len([v for v in meas.values() if v]) + (hexr is not None)
    if meas[ANNULAR[0]] is not None or meas[ANNULAR[1]] is not None:
        sys.exit("gen_zernike_table: NormAnnularNoll loads again -- measure it before the appendix says what it does")
    hexB = hexz.hexagon_basis(True)
    hex_lines = ["| Mode | n | m | Polynomial (x, y in units of lMon) |",
                 "|---------:|----:|----:|:" + "-" * 60 + "|"]       # pandoc widths: "Mode" on one line
    for k in range(1, NHEX_PRINT + 1):
        n, mm = hexz.noll_nm(k)
        hex_lines.append(f"| {k} | {n} | {mm} | {hex_text(hexB[k - 1])} |")

    def ks(v):
        return ", ".join(str(k) for k in v)

    out = [
        "<!-- GENERATED by tools/gen_zernike_table.py (make zernike-table) -- do not edit; rerun. -->",
        "",
        f"**TABLE E.1.** Zernike orderings: the index of each term in every `ZernType=` ordering "
        f"(measured, {nterm} single-term modes; – where an ordering has no such mode)",
        "",
        *lines,
        "",
        f"Measured: {npairs} (type, mode) pairs over {nsettings} type settings.  "
        f"Every mode of ANSI, Born & Wolf, Noll and Fringe is one term of the table with coefficient 1; "
        f"`Malacara` traces exactly as `ANSI`.  "
        + "  ".join(f"`{nz}` = `{p}` with every term multiplied by its RMS factor "
                    f"({norm_ok[nz][0]} modes{', all' if not norm_ok[nz][1] else ''} measured equal to it)."
                    for nz, p in NORM),
        "",
        (f"`ExtFringe` modes {ks(refused['ExtFringe'])} (radial orders 11 to 14) are refused at load by mode, "
         f"and by the coefficient setters, when their coefficient is not zero: the engine converts Extended "
         f"Fringe through the 66-mode ANSI set, which stops at radial order 10, so they would apply no surface.  "
         f"They are not rows of the table." if refused["ExtFringe"] else
         "Every `ExtFringe` mode applies its term."),
        "",
        "`NormAnnularNoll` is refused at load (and by `MODIFY` and the bindings): its converter ignores the "
        "annulus ratio, so the surface traced as circular `NormNoll` at any ratio, and its modes 46–66 applied "
        "nothing.  Use `NormNoll`.",
        "",
    ]
    if hexr is None:
        out += ["`NormHex` is refused at load: the engine accepts no coefficients for it.", ""]
    else:
        worst = max(hexr, key=lambda r: r[1])
        out += [
            f"**TABLE E.2.** `NormHex`: hexagon-orthonormal modes 1–{NHEX_PRINT} (measured, all {len(hexr)} modes)",
            "",
            *hex_lines,
            "",
            f"Measured: each of the {len(hexr)} `NormHex` modes traces as the Noll term of the same index "
            f"orthonormalized over the unit hexagon with vertices on ±xMon (largest difference "
            f"{worst[1]:.0e} of the mode, mode {worst[0]}); the same set with vertices on ±yMon differs in "
            f"{sum(ey > HEX_TOL for _, _, ey in hexr)} modes.  Modes {NHEX_PRINT + 1}–{len(hexr)} continue the "
            f"construction through radial order 8 (Noll index = mode).",
            "",
        ]
    OUT.write_text("\n".join(out))
    short = [ln for ln in lines[2:] if int(ln.split("|")[1]) <= NSHORT]
    OUT8.write_text("\n".join([
        "<!-- GENERATED by tools/gen_zernike_table.py (make zernike-table) -- do not edit; rerun. -->",
        "",
        f"**TABLE 8.** Zernike terms through radial order {NSHORT}, by ordering (all {len(order)} terms, "
        "and the normalized types: Appendix E)",
        "",
        *lines[:2], *short, ""]))
    print(f"wrote {OUT.relative_to(HERE.parent)}: {len(order)} rows; {OUT8.relative_to(HERE.parent)}: {len(short)} rows")


if __name__ == "__main__":
    main()
