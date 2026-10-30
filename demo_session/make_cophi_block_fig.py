#!/usr/bin/env python3
"""COPHI block diagrams for deck_cophi, following Feng Zhao's slide 8
("Implementation Example using 4 Frame Method", COPHI_for_DM_gauge.pptx):
RF/timing chain top-left, optical chain across the middle, the bench at the
right, processing along the bottom -- with OUR components in every box (the
TG96 bench, dm_gauge_lib's stages, the surveyed cameras, the error terms the
model carries).  Four figures:

  figs/cophi_block_single.png   one color (632.8 nm HeNe; the 1.3 um option noted)
  figs/cophi_block_multi.png    two colors frequency-multiplexed on one camera
                                (beats 7.5 Hz and 3.75 Hz: 4 and 8 frames per beat
                                at 30 fps), the synthetic-wavelength unwrap; the
                                sequential (tunable-laser) form noted
  figs/cophi_block_lockin.png   one color on the lock-in camera (5 kHz beat)
  figs/cophi_block_lockin_multi.png  two colors on the lock-in camera, TIME-multiplexed:
                                the same 5 kHz beat for both, the trigger's /8 gates one
                                color's AOM pair on per 0.2 ms frame; the FPGA averages
                                odd and even frames apart (2026-10-29, Dave's question:
                                how are the colors separated, where is the processing)

Routing (re-laid 2026-10-29 after Dave: "some lines cross unrelated boxes, others
are not directly routed"): each RF chain runs straight across one row (synthesizer
-> trigger chain -> driver) and the driver drops straight to its AOM; the 10 MHz
reference runs down the left margin; the trigger goes over the top; the DM command
goes straight up from the driver; the camera data returns under the loop row.  The
two planar-unavoidable crossings (an RF drop over a fiber, the second color's fiber
over the first color's LO fiber) are drawn with a jumper hop.

Build:  ~/dev/MACOS_resources/pymacos/.venv/bin/python make_cophi_block_fig.py
Palette: make_brief_slides.py (accent 1F4E79, ink 282828); Feng's colors for
the signal classes (blue = RF / timing, red = fiber, green = free-space beam,
gray = data).  Dave 2026-10-11: "complete even if large".
"""
import math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Polygon, Rectangle, Arc

ACCENT = "#1F4E79"; INK = "#282828"; GRAY = "#6E6E6E"; LIGHT = "#ECF1F6"
RF = "#1F4E79"; FIBER = "#C0392B"; BEAM = "#2E8B57"; DATA = "#555555"; OPT = "#8E44AD"
GREEN = "#7DBB3C"; NOTE = "#8A6D1A"

W = 28.3


FS = 1.22   # global font scale, 2026-10-17


def box(ax, x, y, w, h, text, fc=ACCENT, tc="white", fs=10.5, bold=False, ec=None, lw=1.2, ls="-", r=0.25):
    fs = fs * FS
    p = FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={r}",
                       fc=fc, ec=ec or fc, lw=lw, linestyle=ls, zorder=3)
    ax.add_patch(p)
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", color=tc,
            fontsize=fs, fontweight="bold" if bold else "normal", zorder=4,
            linespacing=1.15)
    return (x, y, w, h)


def port(b, side, f=0.5):
    x, y, w, h = b
    return {"l": (x, y + h * f), "r": (x + w, y + h * f), "t": (x + w * f, y + h), "b": (x + w * f, y)}[side]


def wire(ax, pts, color=RF, lw=1.8, arrow=True, ls="-", z=2, hops=()):
    """Orthogonal polyline; `hops` = [(x, y), ...] points ON a vertical segment where the
    wire jumps over another line (a small arc, the schematic 'no connection' mark)."""
    for a, b in zip(pts[:-1], pts[1:]):
        segs = [(a, b)]
        for hx, hy in hops:
            if abs(a[0] - b[0]) < 1e-9 and abs(a[0] - hx) < 1e-9 and min(a[1], b[1]) < hy < max(a[1], b[1]):
                r = 0.16
                s = 1 if b[1] > a[1] else -1
                segs = [(a, (hx, hy - s * r)), ((hx, hy + s * r), b)]
                ax.add_patch(Arc((hx, hy), 2 * r, 2 * r, angle=0, theta1=-90, theta2=90,
                                 color=color, lw=lw, ls=ls, zorder=z))
        for p, q in segs:
            ax.plot([p[0], q[0]], [p[1], q[1]], color=color, lw=lw, ls=ls, zorder=z, solid_capstyle="round")
    if arrow:
        a, b = pts[-2], pts[-1]
        ax.annotate("", xy=b, xytext=a, arrowprops=dict(arrowstyle="-|>", color=color, lw=lw, ls=ls,
                                                        shrinkA=0, shrinkB=0, mutation_scale=14), zorder=z)


def label(ax, x, y, s, color=INK, fs=9, ha="left", va="center", style="normal", bold=False):
    fs = fs * FS
    ax.text(x, y, s, color=color, fontsize=fs, ha=ha, va=va, fontstyle=style,
            fontweight="bold" if bold else "normal", zorder=5, linespacing=1.15)


def note(ax, x, y, s, fs=8.2, ha="left"):
    fs = fs * FS
    ax.text(x, y, s, color=NOTE, fontsize=fs, ha=ha, va="center", fontstyle="italic", zorder=5,
            bbox=dict(boxstyle="round,pad=0.25", fc="#FFF8E1", ec="#E0C97A", lw=0.8))


def beam(ax, pts, w=0.28, color=BEAM, alpha=0.35):
    for a, b in zip(pts[:-1], pts[1:]):
        dx, dy = b[0] - a[0], b[1] - a[1]
        n = (dx * dx + dy * dy) ** 0.5
        ox, oy = -dy / n * w / 2, dx / n * w / 2
        ax.add_patch(Polygon([(a[0] + ox, a[1] + oy), (b[0] + ox, b[1] + oy),
                              (b[0] - ox, b[1] - oy), (a[0] - ox, a[1] - oy)],
                             closed=True, fc=color, ec="none", alpha=alpha, zorder=1))


def mirror(ax, x, y, ang_deg, size=1.0, color=INK, lw=3):
    a = math.radians(ang_deg)
    dx, dy = math.cos(a) * size / 2, math.sin(a) * size / 2
    ax.plot([x - dx, x + dx], [y - dy, y + dy], color=color, lw=lw, zorder=4, solid_capstyle="round")


def draw(multi, lockin=False):
    lm = multi and lockin                 # the lock-in two-color variant
    ro = 1.6 if multi else 0.0            # the RF rows ride up in the two-color figures
    Hh = 17.2 if multi else 14.9
    y0 = 10.4 if multi else 7.2           # the first color's band (laser / coupler level)
    yb = y0 - 4.0                         # the second color's band
    py = 1.0                              # the processing row
    fig = plt.figure(figsize=(26, 26 * Hh / W), dpi=200)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W); ax.set_ylim(0, Hh); ax.axis("off")
    beat2 = "80 MHz + 5 kHz" if lockin else "80 MHz + 7.5 Hz"
    beat4 = "80 MHz + 5 kHz" if lm else "80 MHz + 3.75 Hz"

    # ------------------------------------------------------------------ RF / timing (top left): one row per RF line
    R1, R2 = 12.2 + ro, 11.1 + ro
    rb = box(ax, 0.6, 13.5 + ro, 2.5, 0.9, "10 MHz Rb standard\n(SRS FS725 class)", fs=12.3)
    s1 = box(ax, 0.6, R1, 2.5, 0.9, "Synthesizer 1\n80 MHz", fs=12.3)
    s2 = box(ax, 0.6, R2, 2.5, 0.8, "Synthesizer 2\n" + beat2, fs=12.3)
    if lockin:
        cht = ("RF taps on both lines\n-> mixer, low-pass, x4\n-> 20 kHz quarter-period\ntrigger (0-5 V, Rb-locked)\n"
               "-> lock-in REFERENCE\ninput (DivideBy4)")
        if lm:
            cht += "\n/8 -> 2480 Hz color select\n(TTL) to the drivers' blanking"
    else:
        cht = ("RF taps on both lines\n-> x4 multipliers (320 MHz,\n320 MHz + 30 Hz) -> mixer,\nlow-pass -> 30 Hz frame clock\n"
               "-> camera trigger\n(30 fps = 4 frames / beat)")
    ch = box(ax, 6.6, R2, 2.5, R1 + 0.9 - R2, cht, fs=8.6 if lm else 9.2, fc=LIGHT, tc=INK, ec=ACCENT)
    d1 = box(ax, 9.6, R1, 2.5, 0.9, "AOM driver 1\n(RF amp, 2 W)", fs=11.7)
    d2 = box(ax, 9.6, R2, 2.5, 0.8, "AOM driver 2", fs=11.7)
    wire(ax, [port(rb, "b"), port(s1, "t")])
    y1m, y2m = port(s1, "r")[1], port(s2, "r")[1]
    wire(ax, [port(s1, "r"), (6.6, y1m)]); wire(ax, [(9.1, y1m), port(d1, "l")])
    wire(ax, [port(s2, "r"), (6.6, y2m)]); wire(ax, [(9.1, y2m), port(d2, "l")])
    label(ax, 3.3, 13.95 + ro, "one Rb reference, every synthesizer", fs=10.4, color=GRAY)
    # the 10 MHz reference down the left margin
    bus = [port(rb, "l"), (0.35, port(rb, "l")[1]), (0.35, y2m), port(s2, "l")]
    wire(ax, bus)
    if multi:
        RC1, RC2 = 4.4, 3.3
        s3 = box(ax, 0.6, RC1, 2.5, 0.8, "Synthesizer 3\n80 MHz", fs=11.7)
        s4 = box(ax, 0.6, RC2, 2.5, 0.8, "Synthesizer 4\n" + beat4, fs=11.7)
        d3 = box(ax, 3.6, RC1, 2.5, 0.8, "AOM driver 3", fs=11.7)
        d4 = box(ax, 3.6, RC2, 2.5, 0.8, "AOM driver 4", fs=11.7)
        wire(ax, [(0.35, y2m), (0.35, port(s3, "l")[1]), port(s3, "l")])
        wire(ax, [(0.35, port(s3, "l")[1]), (0.35, port(s4, "l")[1]), port(s4, "l")])
        wire(ax, [port(s3, "r"), port(d3, "l")]); wire(ax, [port(s4, "r"), port(d4, "l")])
        label(ax, 7.6, RC1 + 1.1, ("second color: the SAME two frequencies" if lm else "second beat 3.75 Hz: 8 frames per beat"), fs=10.0, color=GRAY)

    # ------------------------------------------------------------------ optical chain (middle)
    las = box(ax, 0.6, y0 - 0.5, 3.5, 1.3, "Single-frequency laser\nHeNe 632.8 nm, stabilized\n(+ isolator)", fs=12.0, fc=OPT)
    cp = box(ax, 4.5, y0 - 0.1, 1.8, 0.9, "1x2 PM\ncoupler", fs=11.7, fc=OPT)
    wire(ax, [port(las, "r"), port(cp, "l")], color=FIBER)
    yt, yl = y0 + 1.45, y0 - 0.75                    # the test and LO fiber levels
    a1 = box(ax, 12.8, yt - 0.45, 1.8, 0.9, "AOM 1\n80 MHz", fs=11.4, fc=OPT)
    a2 = box(ax, 12.0, yl - 0.45, 1.8, 0.9, "AOM 2\n" + beat2, fs=11.0, fc=OPT)
    wire(ax, [port(cp, "r", 0.7), (6.5, port(cp, "r", 0.7)[1]), (6.5, yt), port(a1, "l")], color=FIBER)
    wire(ax, [port(cp, "r", 0.3), (6.5, port(cp, "r", 0.3)[1]), (6.5, yl), port(a2, "l")], color=FIBER)
    label(ax, 6.8, yt + 0.3, "f_opt + 80 MHz", fs=11.1, color=GRAY)
    label(ax, 6.8, yl + 0.3, "f_opt + " + beat2, fs=11.1, color=GRAY)
    # the drivers drop straight to their AOMs (driver 2's drop hops the test fiber)
    xa1 = port(a1, "t")[0]
    wire(ax, [port(d1, "r"), (xa1, y1m), port(a1, "t")], lw=1.3)
    xa2 = port(a2, "t", 0.3)[0]
    wire(ax, [port(d2, "b"), (port(d2, "b")[0], R2 - 0.2), (xa2, R2 - 0.2), port(a2, "t", 0.3)], lw=1.3, hops=[(xa2, yt)])
    if lm:
        wire(ax, [(port(d1, "t")[0], R1 + 0.9 + 0.3), port(d1, "t")], ls=(0, (3, 2)), lw=1.2)
        label(ax, port(d1, "t")[0] + 0.15, R1 + 0.9 + 0.2, "blank on even frames", fs=8.6, color=RF)
        wire(ax, [(port(d2, "b", 0.15)[0], R2 - 0.35), port(d2, "b", 0.15)], ls=(0, (3, 2)), lw=1.2)
        label(ax, port(d2, "b", 0.15)[0] + 0.15, R2 - 0.45, "blank on even", fs=8.6, color=RF)
    note(ax, 0.6, y0 - (1.45 if multi else 1.95),
         ("electronics error terms: reference phase noise; AOM drive ripple -> contrast;\n"
          "ADC bits (2 x 10-bit I/Q, +-512 counts x 460 e-); one I/Q frame = one beat period (0.2 ms)") if lockin else
         ("electronics error terms the model carries: DDS phase noise -> P.cophi.jitter;\n"
          "AOM drive ripple -> contrast; trigger jitter; one beat = one cycle of loop latency (133 ms)"), fs=10.4)
    note(ax, 0.6, y0 - (2.05 if multi else 2.6),
         "fiber-phase walk WITHIN a scan (P.cophi.walk): the term the four-step cannot cancel\n"
         "-> fiber length / temperature stabilization" + ("; a per-frame piston on the lock-in (inert)" if lockin else ""), fs=10.4)

    if multi:
        las2 = box(ax, 0.6, yb - 0.5, 3.5, 1.3, "Laser 2 at lambda_2\n(tunable ECDL at 633 +- dl,\nor telecom tunable at 1.3 um)", fs=11.4, fc=OPT)
        cp2 = box(ax, 4.5, yb - 0.1, 1.8, 0.9, "1x2 PM\ncoupler", fs=11.7, fc=OPT)
        wire(ax, [port(las2, "r"), port(cp2, "l")], color=FIBER)
        yt2, yl2 = yb + 1.05, yb - 0.45
        a3 = box(ax, 12.8, yt2 - 0.45, 1.8, 0.9, "AOM 3\n80 MHz", fs=11.4, fc=OPT)
        a4 = box(ax, 12.0, yl2 - 0.45, 1.8, 0.9, "AOM 4\n" + beat4, fs=11.0, fc=OPT)
        wire(ax, [port(cp2, "r", 0.7), (6.5, port(cp2, "r", 0.7)[1]), (6.5, yt2), port(a3, "l")], color=FIBER)
        wire(ax, [port(cp2, "r", 0.3), (6.5, port(cp2, "r", 0.3)[1]), (6.5, yl2), port(a4, "l")], color=FIBER)
        label(ax, 6.8, yt2 + 0.3, "lambda_2 at 80 MHz", fs=11.1, color=GRAY)
        label(ax, 6.8, yl2 + 0.3, "lambda_2 at " + beat4, fs=11.1, color=GRAY)
        # drivers 3 / 4 come up from below (driver 3's rise hops the second LO fiber)
        wire(ax, [port(d3, "r"), (6.9, port(d3, "r")[1]), (6.9, yb + 0.3), (port(a3, "b")[0], yb + 0.3), port(a3, "b")], lw=1.3, hops=[(6.9, yl2)])
        wire(ax, [port(d4, "r"), (7.3, port(d4, "r")[1]), (7.3, yb - 1.2), (port(a4, "b")[0], yb - 1.2), port(a4, "b")], lw=1.3)
        if lm:
            wire(ax, [(port(d3, "t")[0], RC1 + 0.8 + 0.3), port(d3, "t")], ls=(0, (3, 2)), lw=1.2)
            label(ax, port(d3, "t")[0] + 0.15, RC1 + 0.8 + 0.2, "blank on odd frames", fs=8.6, color=RF)
            wire(ax, [(port(d4, "b")[0], RC2 - 0.3), port(d4, "b")], ls=(0, (3, 2)), lw=1.2)
            label(ax, port(d4, "b")[0] + 0.15, RC2 - 0.2, "blank on odd", fs=8.6, color=RF)
        # the WDM combiners on the first color's fibers; the second color rises to them
        w1 = box(ax, 14.9, yt - 0.45, 1.2, 0.9, "WDM\n(test)", fs=10.4, fc=OPT)
        w2 = box(ax, 14.9, yl - 0.45, 1.2, 0.9, "WDM\n(LO)", fs=10.4, fc=OPT)
        wire(ax, [port(a1, "r"), port(w1, "l")], color=FIBER)
        wire(ax, [port(a2, "r"), port(w2, "l")], color=FIBER)
        wire(ax, [port(a3, "r"), (14.65, yt2), (14.65, yt - 0.75), (port(w1, "b")[0], yt - 0.75), port(w1, "b")], color=FIBER, hops=[(14.65, yl)])
        wire(ax, [port(a4, "r"), (port(w2, "b", 0.7)[0], yl2), port(w2, "b", 0.7)], color=FIBER)
        xt, xl = port(w1, "r")[0], port(w2, "r")[0]
    else:
        xt, xl = port(a1, "r")[0], port(a2, "r")[0]

    # ------------------------------------------------------------------ the bench (right): a Michelson schematic
    bx0 = 17.1
    ax.add_patch(Rectangle((bx0, y0 - 3.4), 11.0, 8.3, fc="#F6F8FA", ec=ACCENT, lw=1.2, ls=(0, (5, 3)), zorder=0))
    label(ax, bx0 + 0.2, y0 + 4.65, "TG96 bench, mirror rig of record -- the reference flat REPLACED by the LO", fs=11.7, bold=True)
    ct = box(ax, bx0 + 0.2, yt - 0.45, 1.8, 0.9, "Test-arm\ncollimator", fs=10.4, fc=LIGHT, tc=INK, ec=OPT)
    cl = box(ax, bx0 + 0.2, yl - 0.45, 1.8, 0.9, "LO collimator\nat the flat's seat", fs=9.6, fc=LIGHT, tc=INK, ec=OPT)
    wire(ax, [(xt, yt), port(ct, "l")], color=FIBER)
    wire(ax, [(xl, yl), port(cl, "l")], color=FIBER)
    if not multi:
        label(ax, (xt + bx0 + 0.2) / 2, yt + 0.35, "PM fiber (test arm)", fs=10.4, color=FIBER, ha="center")
        label(ax, (xl + bx0 + 0.2) / 2, yl - 0.35, "PM fiber (LO arm)", fs=10.4, color=FIBER, ha="center")
    spx, spy = bx0 + 3.3, yt                         # the splitter node on the test beam's line
    mirror(ax, spx, spy, 45, 1.3, color=ACCENT, lw=4)
    label(ax, spx - 0.1, spy - 0.7, "splitter node (plate BS, 22.5 deg)", fs=9.6, ha="right", color=GRAY)
    oap = box(ax, bx0 + 4.9, spy - 0.45, 2.4, 0.9, "imaging tail\n(OAP pair / lens rig L1-L2)", fs=10.9, fc=LIGHT, tc=INK, ec=ACCENT)
    dm = box(ax, bx0 + 8.3, spy - 0.85, 2.3, 1.7, "DM 96 x 96\n(the STOP)\n9216 actuators", fs=12.0, bold=True)
    beam(ax, [port(ct, "r"), (spx, spy)]); beam(ax, [(spx, spy), port(oap, "l")]); beam(ax, [port(oap, "r"), port(dm, "l")])
    label(ax, bx0 + 6.1, spy - 0.75, "double pass: surface h = phi lambda / 4 pi", fs=10.8, color=GRAY, ha="center")
    beam(ax, [port(cl, "r"), (spx, yl), (spx, spy)], color="#2E8B57", alpha=0.22)
    label(ax, spx + 0.3, (spy + yl) / 2 - 0.1, "LO beam: fills the camera flat\n(overfills the pupil by 34 px)", fs=10.1, color=GRAY)
    if lockin:
        CH = "#B9770E"
        cam = box(ax, spx - 2.0, y0 + 2.3, 8.7, 1.9, "", fc="none", tc=INK, ec=ACCENT, lw=1.4)
        label(ax, spx + (1.0 if lm else 1.6), y0 + 4.02, ("IN THE PIXEL -- heliCam C4: one I/Q frame per beat (4960/s), alternate colors" if lm else "LOCK-IN camera (heliCam C4): 542 x 512 px, one I/Q frame per beat period (4960/s)"), fs=8.4, color=ACCENT, bold=True)
        pxw, pxh, pxy = 1.56, 1.15, y0 + 2.45
        p1 = box(ax, spx - 1.85, pxy, pxw, pxh, "photodiode\nI(t) = I0 [1 +\ng cos(phi - w_b t)]", fs=8.4, fc=ACCENT)
        p2 = box(ax, spx - 0.11, pxy, pxw, pxh, "transfer gates\nA|B by s_I (0/180)\nC|D by s_Q (90/270)", fs=8.2, fc=CH)
        p3 = box(ax, spx + 1.63, pxy, pxw, pxh, "storage nodes\nA  B  C  D\n(N cycles / sample)", fs=8.4, fc=LIGHT, tc=INK, ec=CH)
        p4 = box(ax, spx + 3.37, pxy, pxw, pxh, "A - B -> I\nC - D -> Q\n(the DC cancels)", fs=8.4, fc=DATA)
        p5 = box(ax, spx + 5.11, pxy, pxw, pxh, "ADCs on the sensor\n2 x 10 bit per pixel\n-> digital I, Q", fs=8.2, fc=DATA)
        wire(ax, [port(p1, "r"), port(p2, "l")], color=CH, lw=2.0); wire(ax, [port(p2, "r"), port(p3, "l")], color=CH, lw=2.0); wire(ax, [port(p3, "r"), port(p4, "l")], color=DATA); wire(ax, [port(p4, "r"), port(p5, "l")], color=DATA)
        beam(ax, [(spx, spy), (spx, port(p1, "b")[1] + 0.02)])
        label(ax, spx + 0.3, port(cam, "b")[1] - 0.2, "camera tail: the pupil image of the DM", fs=10.1, color=GRAY)
        label(ax, spx - 1.95, y0 + 4.28, "s_I, s_Q = the 5 kHz beat, from the 20 kHz trigger", fs=8.8, color=RF)
        lk = p5
    else:
        cam = box(ax, spx - 1.9, y0 + 2.55, 4.6, 1.45, "Camera, true global shutter\n1k x 1k: Feng baseline 60 ke / 30 fps\nPhotonfocus 200 ke / 150 fps; pco.edge GS", fs=10.4)
        beam(ax, [(spx, spy), (spx, port(cam, "b")[1] + 0.02)])
        label(ax, spx + 0.3, port(cam, "b")[1] - 0.35, "camera tail: the pupil image of the DM", fs=10.1, color=GRAY)
        lk = box(ax, spx + 2.95, y0 + 2.55, 2.6, 1.45, "or: lock-in camera\n(heliCam C4 class,\nbeat at kHz, I/Q out)", fs=10.4, fc=LIGHT, tc=INK, ec=ACCENT, ls=(0, (4, 2)))
    note(ax, bx0 + 0.3, y0 - 2.95, "LO wavefront drift (tilt / defocus per cycle, P.cophi.lo_drift)\n= the new drift row; frozen into the null on the flat DM", fs=10.1)

    # the trigger: over the top, down the bench's left edge, into the camera
    ytr = 13.4 + ro
    if lockin:
        entry = [(port(p1, "r")[0] + 0.12, y0 + 4.4), (port(p1, "r")[0] + 0.12, port(p2, "l")[1]), port(p2, "l")]
    else:
        entry = [(port(cam, "l")[0] - 0.02, y0 + 4.3), (port(cam, "l")[0] - 0.02, port(cam, "l")[1] + 0.4), (port(cam, "l")[0], port(cam, "l")[1] + 0.4)]
    wire(ax, [port(ch, "t"), (port(ch, "t")[0], ytr), (16.7, ytr), (16.7, entry[0][1])] + entry, lw=1.3)
    label(ax, 8.1, ytr + 0.35, ("the 20 kHz quarter-period trigger (Rb-locked) to the lock-in camera's reference input;\none 5 kHz beat period per 0.2 ms frame; its /8 (2480 Hz) selects the color, frame by frame" if lm else
                                "the 20 kHz quarter-period trigger (Rb-locked) to the lock-in camera's reference input:\nevery pixel demodulates one 5 kHz beat period per 0.2 ms frame; frames averaged on the camera; no bucket" if lockin else
                                "30 Hz camera trigger, phase-locked to the beat (4 buckets / beat);\na lock-in camera takes the beat itself as its reference input"), fs=9.6, color=RF)

    # ------------------------------------------------------------------ processing (bottom)
    fg = box(ax, 14.6, py - 0.1, 2.3, 0.9, ("GigE receiver\n(FPGA-averaged I/Q,\n2 x 16-bit, ~20/s)" if lockin else "Frame grabber\n(CXP / CLHS)"), fs=11.4, fc=DATA)
    src = lk if lockin else cam
    xr = bx0 + 10.75                                 # the data return down the right edge
    if lm:
        fp = box(ax, bx0 + 2.6, y0 + 5.3, 7.4, 1.1, "ON THE CAMERA FPGA: the digital I, Q of every frame -> two running averages,\n"
                                                    "odd frames = color 1, even = color 2 (interleaved: to confirm on the loaner; else blocks);\n"
                                                    "background suppression off; the two averaged I/Q images out over GigE (~20/s)", fs=8.0, fc=DATA)
        wire(ax, [port(src, "r"), (bx0 + 10.3, port(src, "r")[1]), (bx0 + 10.3, port(fp, "r", 0.27)[1]), port(fp, "r", 0.27)], color=DATA)
        wire(ax, [port(fp, "r", 0.75), (xr, port(fp, "r", 0.75)[1]), (xr, py - 0.5), (port(fg, "b")[0], py - 0.5), port(fg, "b")], color=DATA)
    else:
        wire(ax, [port(src, "r"), (port(src, "r")[0] + 0.3, port(src, "r")[1]), (port(src, "r")[0] + 0.3, y0 + 4.45), (xr, y0 + 4.45), (xr, py - 0.5), (port(fg, "b")[0], py - 0.5), port(fg, "b")], color=DATA)
    if multi:
        rc = box(ax, 10.5, py - 0.35, 3.9, 1.4, ("atan2(Q, I) per pixel, per color\n(the odd-frame and the even-frame\naveraged images) -> phi_1, phi_2" if lm else "N-bucket LSQ per pixel, TWO beats\n(dmg_cophi_recon at 4 and 8 frames)\n-> phi_1, phi_2"), fs=10.0, fc=DATA)
        sw = box(ax, 6.9, py - 0.35, 3.5, 1.4, "synthetic wavelength\nLambda = l1 l2 / |l1 - l2|\n-> fringe order -> phi unwrapped", fs=10.0, fc=DATA)
        wire(ax, [port(fg, "l"), port(rc, "r")], color=DATA); wire(ax, [port(rc, "l"), port(sw, "r")], color=DATA)
        nxt = sw
        if lm:
            note(ax, 7.0, py - 0.95, "range +-5.2 um (Lambda/4 at dl = 20 nm), coarse noise x33; record c4_capture_oap: zero wrong orders\n"
                                     "to 480 nm rms, 0.13 nm per pixel from ONE frame pair (0.4 ms); one rung: a 200 nm rms cold figure peaks ~1 um", fs=9.6)
        else:
            note(ax, 7.0, py - 0.95, "range +- Lambda/4, coarse noise x Lambda/lambda: a rung LADDER (dl ~ 20 nm then ~0.5 nm) captures\n"
                                     "200 nm rms from cold; sequential form = one tunable laser switched on alternate beats", fs=10.1)
    else:
        rc = box(ax, 10.3, py - 0.35, 3.95, 1.4, ("phase per pixel per averaged frame\natan2(Q, I); 4960 frames/s on the\ncamera -> ~20 images/s to the loop" if lockin else "four-bucket phase per pixel\natan2(I2 - I4, I1 - I3)\n(dmg_cophi_recon, N-bucket LSQ)"), fs=10.4, fc=DATA)
        wire(ax, [port(fg, "l"), port(rc, "r")], color=DATA)
        nxt = rc
        note(ax, 6.9, py - 0.95, ("temporal unwrap between averaged frames (free): the surface cannot move lambda/4 in 50 ms" if lockin else "temporal unwrap between beats (free): the surface cannot move lambda/4 in 133 ms"), fs=10.1)
    nl = box(ax, 3.85, py - 0.35, 2.85, 1.4, "null reference\n(frozen: flat DM;\nLO wavefront in it)", fs=10.9, fc=DATA)
    wire(ax, [port(nxt, "l"), port(nl, "r")], color=DATA)
    est = box(ax, 0.6, py - 0.35, 3.0, 1.4, "matrix estimator\n(dmg_act_fit, lambda_m)\n-> 9216 actuator heights", fs=10.9, fc=DATA)
    wire(ax, [port(nl, "l"), port(est, "r")], color=DATA)
    label(ax, 0.6, py + 1.55, "IN THE COMPUTER: " + ("the receiver -> the phases -> the synthetic wavelength -> the null -> the estimator -> the loop" if multi else "the receiver -> the phase -> the null -> the estimator -> the loop"), fs=10.4, color=DATA, bold=True)
    lp = box(ax, 21.3, py - 0.35, 3.0, 1.0, ("hold loop (dmg_loop)\ng = 0.5 at 8 Hz (0.125 s)" if lockin else "hold loop (dmg_loop)\ng = 0.5 at 7.5 Hz"), fs=11.2, fc=GREEN, tc=INK)
    wire(ax, [port(est, "t"), (port(est, "t")[0], py + 1.25), (port(lp, "t")[0], py + 1.25), port(lp, "t")], color=DATA)
    drv = box(ax, 24.7, py - 0.35, 3.0, 1.0, "DM driver\n(HV electronics, 9216 ch)", fs=11.2, fc=GREEN, tc=INK)
    wire(ax, [port(lp, "r"), port(drv, "l")], color=DATA)
    wire(ax, [port(drv, "t"), (port(drv, "t")[0], py + 1.3), (port(dm, "b")[0], py + 1.3), port(dm, "b")], color=DATA, lw=1.3)
    out = box(ax, 14.6, py + 1.9, 6.3, 0.9, "OUTPUT: DM surface figure vs time (pm);\nthe hold residual; the capture record", fs=10.6, fc="#DFF0D8", tc=INK, ec=GREEN)
    label(ax, 17.6, py - 0.8, ("per frame pair (0.4 ms): 0.13 nm per pixel, zero wrong orders to 480 nm rms (c4_capture_oap); the surface is static across the pair" if lm else
                               "per I/Q sample: sigma_phi^2 = (pi^2/8) 2 (I0 + r^2) / (I0^2 g^2), I0 = 2 x 236 ke (the output's full scale), r = 489 e- -- 0.77 s per pm (mirror rig)" if lockin else
                               "camera rows: photons / frame = fill x well; sigma_phi^2 = 2 (I0 + r^2) / (N I0^2 g'^2); g' = g x 0.900 (integrating bucket)"), fs=10.8, color=GRAY)

    name = "cophi_block_lockin_multi.png" if lm else ("cophi_block_lockin.png" if lockin else ("cophi_block_multi.png" if multi else "cophi_block_single.png"))
    fig.savefig(f"figs/{name}", dpi=200, facecolor="white", bbox_inches="tight", pad_inches=0.08)
    fig.savefig(f"figs/{name.replace('.png', '.svg')}", facecolor="white", bbox_inches="tight", pad_inches=0.08)
    print("wrote figs/" + name)


draw(False)
draw(True)
draw(False, lockin=True)
draw(True, lockin=True)
