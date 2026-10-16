#!/usr/bin/env python3
"""COPHI block diagrams for deck_cophi, following Feng Zhao's slide 8
("Implementation Example using 4 Frame Method", COPHI_for_DM_gauge.pptx):
RF/timing chain top-left, optical chain across the middle, the bench at the
right, processing along the bottom -- with OUR components in every box (the
TG96 bench, dm_gauge_lib's stages, the surveyed cameras, the error terms the
model carries).  Two figures:

  figs/cophi_block_single.png   one color (632.8 nm HeNe; the 1.3 um option noted)
  figs/cophi_block_multi.png    two colors frequency-multiplexed on one camera
                                (beats 7.5 Hz and 3.75 Hz: 4 and 8 frames per beat
                                at 30 fps), the synthetic-wavelength unwrap; the
                                sequential (tunable-laser) form noted

Build:  ~/dev/MACOS_resources/pymacos/.venv/bin/python make_cophi_block_fig.py
Palette: make_brief_slides.py (accent 1F4E79, ink 282828); Feng's colors for
the signal classes (blue = RF / timing, red = fiber, green = free-space beam,
gray = data).  Dave 2026-10-11: "complete even if large".
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Polygon, Rectangle

ACCENT = "#1F4E79"; INK = "#282828"; GRAY = "#6E6E6E"; LIGHT = "#ECF1F6"
RF = "#1F4E79"; FIBER = "#C0392B"; BEAM = "#2E8B57"; DATA = "#555555"; OPT = "#8E44AD"
GREEN = "#7DBB3C"; NOTE = "#8A6D1A"

W, H = 32.0, 20.0  # the two-color figure uses H + 2.4


def box(ax, x, y, w, h, text, fc=ACCENT, tc="white", fs=10.5, bold=False, ec=None, lw=1.2, ls="-", r=0.25):
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


def wire(ax, pts, color=RF, lw=1.8, arrow=True, ls="-", z=2):
    for a, b in zip(pts[:-1], pts[1:]):
        ax.plot([a[0], b[0]], [a[1], b[1]], color=color, lw=lw, ls=ls, zorder=z, solid_capstyle="round")
    if arrow:
        a, b = pts[-2], pts[-1]
        ax.annotate("", xy=b, xytext=a, arrowprops=dict(arrowstyle="-|>", color=color, lw=lw, ls=ls,
                                                        shrinkA=0, shrinkB=0, mutation_scale=14), zorder=z)


def label(ax, x, y, s, color=INK, fs=9, ha="left", va="center", style="normal", bold=False):
    ax.text(x, y, s, color=color, fontsize=fs, ha=ha, va=va, fontstyle=style,
            fontweight="bold" if bold else "normal", zorder=5, linespacing=1.15)


def note(ax, x, y, s, fs=8.2, ha="left"):
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
    import math
    a = math.radians(ang_deg)
    dx, dy = math.cos(a) * size / 2, math.sin(a) * size / 2
    ax.plot([x - dx, x + dx], [y - dy, y + dy], color=color, lw=lw, zorder=4, solid_capstyle="round")


def draw(multi, lockin=False):
    Hh = H + (2.4 if multi else 0) - 1.5
    fig = plt.figure(figsize=(26, 26 * Hh / W), dpi=200)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W); ax.set_ylim(0, Hh); ax.axis("off")

    # (no in-figure title: the slide carries it -- Dave 2026-10-12)
    # ------------------------------------------------------------------ RF / timing (top left)
    ro = 0.0 if not multi else 2.4                   # the RF column rides up in the two-color figure
    rb = box(ax, 3.2, 16.8 + ro, 2.6, 1.0, "10 MHz Rb standard\n(SRS FS725 class)", fs=12.3)
    s1 = box(ax, 0.6, 14.8 + ro, 2.5, 1.0, "Synthesizer 1\n80 MHz", fs=12.3)
    s2 = box(ax, 3.6, 14.8 + ro, 2.5, 1.0, ("Synthesizer 2\n80 MHz + 5 kHz" if lockin else "Synthesizer 2\n80 MHz + 7.5 Hz"), fs=12.3)
    wire(ax, [port(rb, "b"), (port(rb, "b")[0], 16.2 + ro), (1.85, 16.2 + ro), port(s1, "t")])
    wire(ax, [(port(rb, "b")[0], 16.2 + ro), (4.85, 16.2 + ro), port(s2, "t")])
    if multi:
        s3 = box(ax, 6.7, 14.8 + ro, 2.5, 1.0, "Synthesizer 3\n80 MHz", fs=12.3)
        s4 = box(ax, 9.7, 14.8 + ro, 2.5, 1.0, "Synthesizer 4\n80 MHz + 3.75 Hz", fs=12.3)
        wire(ax, [(4.85, 16.2 + ro), (7.95, 16.2 + ro), port(s3, "t")])
        wire(ax, [(7.95, 16.2 + ro), (10.95, 16.2 + ro), port(s4, "t")])
        label(ax, 12.5, 16.2 + ro, "one reference for all four synthesizers: the two beats are coherent", fs=11.1, color=GRAY)
    d1 = box(ax, 0.6, 13.0 + ro, 2.5, 0.9, "AOM driver 1\n(RF amp, 2 W)", fs=11.7)
    d2 = box(ax, 3.6, 13.0 + ro, 2.5, 0.9, "AOM driver 2", fs=11.7)
    wire(ax, [port(s1, "b"), port(d1, "t")]); wire(ax, [port(s2, "b"), port(d2, "t")])
    if lockin:
        m1 = box(ax, 0.6, 11.5 + ro, 2.5, 0.8, "RF tap\n80 MHz", fs=11.4, fc=LIGHT, tc=INK, ec=ACCENT)
        m2 = box(ax, 3.6, 11.5 + ro, 2.5, 0.8, "RF tap\n80 MHz + 5 kHz", fs=11.4, fc=LIGHT, tc=INK, ec=ACCENT)
        mx = box(ax, 1.6, 9.9 + ro, 3.5, 0.9, "Mixer + low-pass\n-> 5 kHz beat reference", fs=11.7)
        tg = box(ax, 1.6, 8.5 + ro, 3.5, 0.85, "Lock-in REFERENCE input\n(the beat itself; no frame trigger)", fs=10.8)
    else:
        m1 = box(ax, 0.6, 11.5 + ro, 2.5, 0.8, "x4 multiplier\n320 MHz", fs=11.4, fc=LIGHT, tc=INK, ec=ACCENT)
        m2 = box(ax, 3.6, 11.5 + ro, 2.5, 0.8, "x4 multiplier\n320 MHz + 30 Hz", fs=11.4, fc=LIGHT, tc=INK, ec=ACCENT)
        mx = box(ax, 1.6, 9.9 + ro, 3.5, 0.9, "Mixer + low-pass\n-> 30 Hz frame clock", fs=11.7)
        tg = box(ax, 1.6, 8.5 + ro, 3.5, 0.85, "Camera trigger\n(30 fps = 4 frames / beat)", fs=11.7)
    wire(ax, [port(s1, "l"), (0.35, port(s1, "l")[1]), (0.35, port(m1, "l")[1]), port(m1, "l")])
    wire(ax, [port(s2, "r"), (6.35, port(s2, "r")[1]), (6.35, port(m2, "r")[1]), port(m2, "r")])
    wire(ax, [port(m1, "b"), (port(m1, "b")[0], 11.2 + ro), (2.6, 11.2 + ro), (2.6, 10.8 + ro)])
    wire(ax, [port(m2, "b"), (port(m2, "b")[0], 11.2 + ro), (4.1, 11.2 + ro), (4.1, 10.8 + ro)])
    wire(ax, [port(mx, "b"), port(tg, "t")])
    if multi:
        d3 = box(ax, 6.7, 13.0 + ro, 2.5, 0.9, "AOM driver 3", fs=11.7)
        d4 = box(ax, 9.7, 13.0 + ro, 2.5, 0.9, "AOM driver 4", fs=11.7)
        wire(ax, [port(s3, "b"), port(d3, "t")]); wire(ax, [port(s4, "b"), port(d4, "t")])
        note(ax, 12.5, 13.45 + ro, "second beat 3.75 Hz = 8 frames per beat at the same 30 fps: both beats are\n"
                                   "integer frame counts, demodulated together by the N-bucket least squares", fs=10.4)

    # ------------------------------------------------------------------ optical chain (middle)
    y0 = 7.2 if not multi else 9.6
    py = 3.3 if not multi else 3.0
    las = box(ax, 0.6, y0 - 0.5, 2.9, 1.3, "Single-frequency laser\nHeNe 632.8 nm, stabilized\n(+ isolator)", fs=12.0, fc=OPT)
    cp = box(ax, 4.2, y0 - 0.1, 1.9, 0.9, "1x2 PM\ncoupler", fs=11.7, fc=OPT)
    wire(ax, [port(las, "r"), port(cp, "l")], color=FIBER)
    a1 = box(ax, 7.0, y0 + 1.0, 2.3, 0.9, "AOM 1\n80 MHz", fs=12.0, fc=OPT)
    a2 = box(ax, 7.0, y0 - 1.2, 2.3, 0.9, ("AOM 2\n80 MHz + 5 kHz" if lockin else "AOM 2\n80 MHz + 7.5 Hz"), fs=12.0, fc=OPT)
    wire(ax, [port(cp, "r", 0.7), (6.5, port(cp, "r", 0.7)[1]), (6.5, port(a1, "l")[1]), port(a1, "l")], color=FIBER)
    wire(ax, [port(cp, "r", 0.3), (6.5, port(cp, "r", 0.3)[1]), (6.5, port(a2, "l")[1]), port(a2, "l")], color=FIBER)
    wire(ax, [port(d1, "b"), (port(d1, "b")[0], 12.5 + ro), (7.55, 12.5 + ro), (7.55, port(a1, "t")[1] + 0.02)], lw=1.2)
    wire(ax, [port(d2, "b"), (port(d2, "b")[0], 12.2 + ro), (8.65, 12.2 + ro), (8.65, port(a2, "t")[1] + 0.02)], lw=1.2)
    label(ax, 8.15, y0 + 2.15, "f_opt + 80 MHz", fs=11.1, color=GRAY, ha="center")
    label(ax, 8.15, y0 - 1.55, ("f_opt + 80 MHz + 5 kHz" if lockin else "f_opt + 80 MHz + 7.5 Hz"), fs=11.1, color=GRAY, ha="center")
    note(ax, 0.6, y0 - 1.95, ("electronics error terms: reference phase noise; AOM drive ripple -> contrast;\n"
                              "ADC bits (2 x 10-bit I/Q); one I/Q frame = one cycle of loop latency (1 ms)") if lockin else
                             ("electronics error terms the model carries: DDS phase noise -> P.cophi.jitter;\n"
                             "AOM drive ripple -> contrast; trigger jitter; one beat = one cycle of loop latency (133 ms)"), fs=10.4)

    if multi:
        yb = y0 - 4.0                                  # the second color's band, below the first
        las2 = box(ax, 0.6, yb - 0.5, 2.9, 1.3, "Laser 2 at lambda_2\n(tunable ECDL at 633 +- dl,\nor telecom tunable at 1.3 um)", fs=11.4, fc=OPT)
        cp2 = box(ax, 4.2, yb - 0.1, 1.9, 0.9, "1x2 PM\ncoupler", fs=11.7, fc=OPT)
        wire(ax, [port(las2, "r"), port(cp2, "l")], color=FIBER)
        a3 = box(ax, 7.0, yb + 0.6, 2.3, 0.9, "AOM 3\n80 MHz", fs=12.0, fc=OPT)
        a4 = box(ax, 7.0, yb - 0.9, 2.3, 0.9, "AOM 4\n80 MHz + 3.75 Hz", fs=12.0, fc=OPT)
        wire(ax, [port(cp2, "r", 0.7), (6.5, port(cp2, "r", 0.7)[1]), (6.5, port(a3, "l")[1]), port(a3, "l")], color=FIBER)
        wire(ax, [port(cp2, "r", 0.3), (6.5, port(cp2, "r", 0.3)[1]), (6.5, port(a4, "l")[1]), port(a4, "l")], color=FIBER)
        # drivers 3 / 4 come down the gap right of AOM 1 / 2
        wire(ax, [port(d3, "b"), (port(d3, "b")[0], 12.2 + ro), (9.7, 12.2 + ro), (9.7, yb + 1.9), (7.55, yb + 1.9), (7.55, port(a3, "t")[1] + 0.02)], lw=1.2)
        wire(ax, [port(d4, "b"), (port(d4, "b")[0], 11.9 + ro), (9.95, 11.9 + ro), (9.95, yb + 0.3), (8.65, yb + 0.3), (8.65, port(a4, "t")[1] + 0.02)], lw=1.2)
        w1 = box(ax, 13.3, y0 + 1.0, 2.1, 0.9, "WDM combiner\n(test fiber)", fs=11.4, fc=OPT)
        w2 = box(ax, 13.3, y0 - 1.2, 2.1, 0.9, "WDM combiner\n(LO fiber)", fs=11.4, fc=OPT)
        wire(ax, [port(a1, "r"), port(w1, "l", 0.7)], color=FIBER)
        wire(ax, [port(a2, "r"), port(w2, "l", 0.7)], color=FIBER)
        wire(ax, [port(a3, "r"), (12.9, port(a3, "r")[1]), (12.9, port(w1, "l", 0.3)[1]), port(w1, "l", 0.3)], color=FIBER)
        wire(ax, [port(a4, "r"), (12.6, port(a4, "r")[1]), (12.6, port(w2, "l", 0.3)[1]), port(w2, "l", 0.3)], color=FIBER)
        label(ax, 9.45, yb + 1.05, "lambda_2 at 80 MHz", fs=11.1, color=GRAY)
        label(ax, 9.45, yb - 0.45, "lambda_2 at 80 MHz + 3.75 Hz", fs=11.1, color=GRAY)
        xt, yt = port(w1, "r"); xl, yl = port(w2, "r")
    else:
        xt, yt = port(a1, "r"); xl, yl = port(a2, "r")

    # the two PM fibers to the bench
    ct = box(ax, 17.0, y0 + 1.0, 2.5, 0.9, "Test-arm collimator\n(the bench's laser input)", fs=11.2, fc=LIGHT, tc=INK, ec=OPT)
    cl = box(ax, 17.0, y0 - 1.2, 2.5, 0.9, "LO collimator AT the\nreference flat's seat", fs=11.2, fc=LIGHT, tc=INK, ec=OPT)
    wire(ax, [(xt, yt), port(ct, "l")], color=FIBER)
    wire(ax, [(xl, yl), port(cl, "l")], color=FIBER)
    label(ax, (xt + 17.0) / 2, yt + 0.35, "PM fiber (test arm)", fs=11.1, color=FIBER, ha="center")
    label(ax, (xl + 17.0) / 2, yl - 0.35, "PM fiber (LO arm)", fs=11.1, color=FIBER, ha="center")
    note(ax, 13.4, y0 - 1.9, "fiber-phase walk WITHIN a scan (P.cophi.walk): the term\n"
                             "the four-step cannot cancel -> fiber length / temperature stabilization", fs=10.4)

    # ------------------------------------------------------------------ the bench (right): a Michelson schematic
    bx0 = 20.3
    ax.add_patch(Rectangle((bx0, y0 - 3.4), 11.4, 8.3, fc="#F6F8FA", ec=ACCENT, lw=1.2, ls=(0, (5, 3)), zorder=0))
    label(ax, bx0 + 0.2, y0 + 4.65, "TG96 bench, mirror rig of record -- the reference flat REPLACED by the LO (model 1024, grid 385)", fs=11.7, bold=True)
    spx, spy = bx0 + 3.6, port(ct, "r")[1]           # the splitter node on the test beam's line
    mirror(ax, spx, spy, 45, 1.3, color=ACCENT, lw=4)
    label(ax, spx - 1.0, spy + 0.75, "splitter node\n(plate BS, 22.5 deg)", fs=10.7, ha="right", color=GRAY)
    # test beam: collimator -> splitter -> tail -> DM, all on one line
    oap = box(ax, bx0 + 5.4, spy - 0.45, 2.4, 0.9, "imaging tail\n(OAP pair / lens rig L1-L2)", fs=10.9, fc=LIGHT, tc=INK, ec=ACCENT)
    dm = box(ax, bx0 + 8.8, spy - 0.85, 2.3, 1.7, "DM 96 x 96\n(the STOP)\n9216 actuators", fs=12.0, bold=True)
    beam(ax, [port(ct, "r"), (spx, spy)]); beam(ax, [(spx, spy), port(oap, "l")]); beam(ax, [port(oap, "r"), port(dm, "l")])
    label(ax, bx0 + 6.6, spy - 0.75, "double pass: surface h = phi lambda / 4 pi", fs=10.8, color=GRAY, ha="center")
    # the LO: fiber collimator BELOW the splitter, beam up into it (the flat's port)
    beam(ax, [port(cl, "r"), (spx, port(cl, "r")[1]), (spx, spy)], color="#2E8B57", alpha=0.22)
    label(ax, spx + 0.3, (spy + port(cl, "r")[1]) / 2 - 0.1, "LO beam: fills the camera flat\n(overfills the pupil by 34 px)", fs=10.1, color=GRAY)
    # the camera: ABOVE the splitter (the fourth port), its tail = the pupil image of the DM
    if lockin:
        # the lock-in camera drawn as the chain inside ONE pixel (the slide-6 signal flow), in series:
        # photodiode -> gates (the beat reference) -> nodes -> differences -> I, Q out (Dave 2026-10-16)
        CH = "#B9770E"
        cam = box(ax, spx - 2.0, y0 + 2.3, 9.1, 1.9, "", fc="none", tc=INK, ec=ACCENT, lw=1.4)
        label(ax, spx + 2.7, y0 + 4.0, "LOCK-IN camera (heliCam C4): 542 x 512 pixels, I/Q at 1 kHz", fs=8.4, color=ACCENT, bold=True)
        pxw, pxh, pxy = 1.62, 1.15, y0 + 2.45
        p1 = box(ax, spx - 1.85, pxy, pxw, pxh, "photodiode\nI(t) = I0 [1 +\ng cos(phi - w_b t)]", fs=8.4, fc=ACCENT)
        p2 = box(ax, spx - 0.03, pxy, pxw, pxh, "transfer gates\nA|B by s_I (0/180)\nC|D by s_Q (90/270)", fs=8.2, fc=CH)
        p3 = box(ax, spx + 1.79, pxy, pxw, pxh, "storage nodes\nA  B  C  D\n(N cycles / sample)", fs=8.4, fc=LIGHT, tc=INK, ec=CH)
        p4 = box(ax, spx + 3.61, pxy, pxw, pxh, "A - B -> I\nC - D -> Q\n(the DC cancels)", fs=8.4, fc=DATA)
        p5 = box(ax, spx + 5.43, pxy, pxw, pxh, "ADCs on the sensor\n2 x 10 bit per pixel\n-> digital I, Q", fs=8.2, fc=DATA)
        wire(ax, [port(p1, "r"), port(p2, "l")], color=CH, lw=2.0); wire(ax, [port(p2, "r"), port(p3, "l")], color=CH, lw=2.0); wire(ax, [port(p3, "r"), port(p4, "l")], color=DATA); wire(ax, [port(p4, "r"), port(p5, "l")], color=DATA)
        beam(ax, [(spx, spy), (spx, port(p1, "b")[1] + 0.02)])
        label(ax, spx + 0.3, port(cam, "b")[1] - 0.35, "camera tail: the pupil image of the DM", fs=10.1, color=GRAY)
        label(ax, spx - 1.95, y0 + 4.28, "s_I, s_Q = the 5 kHz beat (square wave)", fs=8.8, color=RF)
        lk = p5
    else:
        cam = box(ax, spx - 1.9, y0 + 2.55, 4.6, 1.45, "Camera, true global shutter\n1k x 1k: Feng baseline 60 ke / 30 fps\nPhotonfocus 200 ke / 150 fps; pco.edge GS", fs=10.4)
        beam(ax, [(spx, spy), (spx, port(cam, "b")[1] + 0.02)])
        label(ax, spx + 0.3, port(cam, "b")[1] - 0.35, "camera tail: the pupil image of the DM", fs=10.1, color=GRAY)
        lk = box(ax, spx + 2.95, y0 + 2.55, 2.6, 1.45, "or: lock-in camera\n(heliCam C4 class,\nbeat at kHz, I/Q out)", fs=10.4, fc=LIGHT, tc=INK, ec=ACCENT, ls=(0, (4, 2)))
    note(ax, bx0 + 4.9, y0 - 2.95, "LO wavefront drift (tilt / defocus per cycle, P.cophi.lo_drift)\n= the new drift row; frozen into the null on the flat DM", fs=10.1)
    # trigger to the camera (routed under the fibers' notes, above the processing row)
    ty = py + 1.4
    if lockin:
        wire(ax, [port(tg, "r"), (6.6, port(tg, "r")[1]), (6.6, ty), (bx0 - 0.4, ty), (bx0 - 0.4, y0 + 4.4), (port(p1, "r")[0] + 0.12, y0 + 4.4), (port(p1, "r")[0] + 0.12, port(p2, "l")[1]), port(p2, "l")], lw=1.3)
    else:
        wire(ax, [port(tg, "r"), (6.6, port(tg, "r")[1]), (6.6, ty), (bx0 - 0.4, ty), (bx0 - 0.4, y0 + 4.3), (port(cam, "l")[0] - 0.02, y0 + 4.3), (port(cam, "l")[0] - 0.02, port(cam, "l")[1] + 0.4), (port(cam, "l")[0], port(cam, "l")[1] + 0.4)], lw=1.3)
    label(ax, 6.8, ty + 0.22, ("the 5 kHz beat reference to the lock-in camera: every pixel demodulates ~5 cycles per 1 ms I/Q frame; no bucket, no frame trigger" if lockin else "30 Hz camera trigger, phase-locked to the beat (4 buckets / beat); a lock-in camera takes the beat itself as its reference input"), fs=10.8, color=RF)

    # ------------------------------------------------------------------ processing (bottom)
    fg = box(ax, 14.2, py - 0.1, 2.3, 0.9, ("GigE receiver\n(digital I/Q frames)" if lockin else "Frame grabber\n(CXP / CLHS)"), fs=11.4, fc=DATA)
    src = lk if lockin else cam
    wire(ax, [port(src, "r"), (port(src, "r")[0] + 0.3, port(src, "r")[1]), (port(src, "r")[0] + 0.3, y0 + 4.45), (bx0 + 11.2, y0 + 4.45), (bx0 + 11.2, py + 0.35), port(fg, "r")], color=DATA)
    if multi:
        rc = box(ax, 10.3, py - 0.35, 3.6, 1.4, "N-bucket LSQ per pixel, TWO beats\n(dmg_cophi_recon at 4 and 8 frames)\n-> phi_1, phi_2", fs=10.9, fc=DATA)
        sw = box(ax, 6.9, py - 0.35, 3.2, 1.4, "synthetic wavelength\nLambda = l1 l2 / |l1 - l2|\n-> fringe order -> phi unwrapped", fs=10.9, fc=DATA)
        wire(ax, [port(fg, "l"), port(rc, "r")], color=DATA); wire(ax, [port(rc, "l"), port(sw, "r")], color=DATA)
        nxt = sw
        note(ax, 7.0, py - 0.95, "range +- Lambda/4, coarse noise x Lambda/lambda: a rung LADDER (dl ~ 20 nm then ~0.5 nm) captures\n"
                                 "200 nm rms from cold; sequential form = one tunable laser switched on alternate beats", fs=10.1)
    else:
        rc = box(ax, 10.3, py - 0.35, 3.6, 1.4, ("phase per pixel per I/Q frame\natan2(Q, I) at 1 kHz\n-> decimate / average to the loop rate" if lockin else "four-bucket phase per pixel\natan2(I2 - I4, I1 - I3)\n(dmg_cophi_recon, N-bucket LSQ)"), fs=10.9, fc=DATA)
        wire(ax, [port(fg, "l"), port(rc, "r")], color=DATA)
        nxt = rc
        note(ax, 10.6, py - 0.95, ("temporal unwrap between I/Q frames (free): the surface cannot move lambda/4 in 1 ms" if lockin else "temporal unwrap between beats (free): the surface cannot move lambda/4 in 133 ms"), fs=10.1)
    nl = box(ax, 3.85, py - 0.35, 2.85, 1.4, "null reference\n(frozen: flat DM;\nLO wavefront in it)", fs=10.9, fc=DATA)
    wire(ax, [port(nxt, "l"), port(nl, "r")], color=DATA)
    est = box(ax, 0.6, py - 0.35, 3.0, 1.4, "matrix estimator\n(dmg_act_fit, lambda_m)\n-> 9216 actuator heights", fs=10.9, fc=DATA)
    wire(ax, [port(nl, "l"), port(est, "r")], color=DATA)
    lp = box(ax, 0.6, 0.55, 3.0, 1.0, ("hold loop (dmg_loop)\ng = 0.5 at up to 1 kHz" if lockin else "hold loop (dmg_loop)\ng = 0.5 at 7.5 Hz"), fs=11.2, fc=GREEN, tc=INK)
    wire(ax, [port(est, "b"), port(lp, "t")], color=DATA)
    drv = box(ax, 4.2, 0.55, 3.0, 1.0, "DM driver\n(HV electronics, 9216 ch)", fs=11.2, fc=GREEN, tc=INK)
    wire(ax, [port(lp, "r"), port(drv, "l")], color=DATA)
    wire(ax, [port(drv, "r", 0.25), (bx0 + 11.6, port(drv, "r", 0.25)[1]), (bx0 + 11.6, port(dm, "r")[1]), (port(dm, "r")[0], port(dm, "r")[1])], color=DATA, lw=1.3)
    out = box(ax, 7.6, 0.55, 8.4, 1.0, "OUTPUT: DM surface figure vs time (pm); the hold residual; the capture record", fs=11.2, fc="#DFF0D8", tc=INK, ec=GREEN)
    label(ax, 16.4, 1.05, ("per I/Q sample: sigma_phi^2 = (pi^2/8) 2 (I0 + r^2) / (I0^2 g^2), I0 = well x fill; well-capped -- 3.4 s per pm as surveyed" if lockin else "camera rows: photons / frame = fill x well; sigma_phi^2 = 2 (I0 + r^2) / (N I0^2 g'^2); g' = g x 0.900 (integrating bucket)"), fs=10.8, color=GRAY)

    name = "cophi_block_lockin.png" if lockin else ("cophi_block_multi.png" if multi else "cophi_block_single.png")
    fig.savefig(f"figs/{name}", dpi=200, facecolor="white")
    fig.savefig(f"figs/{name.replace('.png', '.svg')}", facecolor="white")
    print("wrote figs/" + name)


draw(False)
draw(True)
draw(False, lockin=True)
