#!/usr/bin/env python3
"""deck_cophi: one bench, several sensors (Dave 2026-10-17: the bench is a testbed,
multiple sensors are a virtue).  The TG96 mirror rig with the LO collimator at the
reference flat's seat, an INSERTABLE reference flat in the LO path (flip mount) that
gives the IFO back, and a pick-off behind the splitter feeding the Zernike sensor and
the pinhole / P-SRI arm from the test beam alone.
Build:  ~/dev/MACOS_resources/pymacos/.venv/bin/python make_cophi_multisensor_fig.py
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, Circle

ACCENT = "#1F4E79"; INK = "#282828"; GRAY = "#6E6E6E"; LIGHT = "#ECF1F6"
FIBER = "#C0392B"; BEAM = "#2E8B57"; LO = "#8E44AD"; GREEN = "#7DBB3C"; NOTE = "#8A6D1A"; ORANGE = "#B9770E"
W, H = 24.4, 14.0
fig = plt.figure(figsize=(26, 26 * H / W), dpi=200)
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W); ax.set_ylim(0, H); ax.axis("off")

FS = 1.3

def box(x, y, w, h, text, fc=ACCENT, tc="white", fs=13, bold=False, ec=None, ls="-"):
    fs = fs * FS
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.2", fc=fc,
                                ec=ec or fc, lw=1.4, ls=ls, zorder=3))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", color=tc, fontsize=fs,
            fontweight="bold" if bold else "normal", zorder=4)
    return (x, y, w, h)

def beam(pts, color=BEAM, w=0.26, alpha=0.35, z=1):
    xs, ys = zip(*pts); ax.plot(xs, ys, color=color, lw=w * 36, alpha=alpha, solid_capstyle="butt", zorder=z)

def mirror(x, y, ang, size=1.2, color=INK, lw=4):
    import math
    dx, dy = math.cos(math.radians(ang)) * size / 2, math.sin(math.radians(ang)) * size / 2
    ax.plot([x - dx, x + dx], [y - dy, y + dy], color=color, lw=lw, zorder=5)

def label(x, y, s, fs=12, color=INK, ha="left", bold=False, style="normal"):
    fs = fs * FS
    ax.text(x, y, s, fontsize=fs, color=color, ha=ha, va="center", fontweight="bold" if bold else "normal",
            fontstyle=style, zorder=6)

# ---- the test beam line (y = 8.0): input collimator -> splitter -> tail -> DM
ty = 8.0
ci = box(0.6, ty - 0.5, 2.8, 1.0, "test-arm collimator\n(PM fiber, AOM 1)", fc=LIGHT, tc=INK, ec=FIBER)
sx = 7.0
beam([(3.4, ty), (sx, ty)]); beam([(sx, ty), (17.4, ty)])
mirror(sx, ty, 45, 1.6, ACCENT, 5); label(sx - 0.2, ty + 1.05, "main splitter (plate BS)", fs=11.5, color=GRAY, ha="right")
tail = box(10.6, ty - 0.5, 3.2, 1.0, "imaging tail\n(OAP pair)", fc=LIGHT, tc=INK, ec=ACCENT)
dm = box(17.4, ty - 1.1, 2.6, 2.2, "DM 96 x 96\n(the STOP)", bold=True)
label(10.6, ty - 1.0, "double pass; the pupil imaged onto every sensor", fs=11, color=GRAY)

# ---- the LO arm BELOW the splitter: LO collimator at the flat's seat, insertable reference flat in front
beam([(sx, ty), (sx, 2.6)], color=LO, alpha=0.25)
lo = box(sx - 2.1, 1.3, 4.2, 1.1, "LO collimator\n(PM fiber, AOM 2) AT the flat's seat", fc=LIGHT, tc=INK, ec=LO, fs=11.5)
mirror(sx, 4.6, 0, 1.5, ORANGE, 6)
ax.annotate("", xy=(sx + 1.9, 4.6), xytext=(sx + 0.9, 4.6), arrowprops=dict(arrowstyle="<->", color=ORANGE, lw=2))
label(sx + 2.05, 4.6, "INSERTABLE reference flat (flip mount)\nin: the IFO's reference arm, LO shuttered\nout: COPHI, the LO fills the camera", fs=10.5, color=ORANGE)
label(sx - 0.3, 3.5, "LO beam, overfills the pupil", fs=10.5, color=LO, ha="right")

# ---- the camera port ABOVE the splitter with a pick-off cascade
cy = 12.3
beam([(sx, ty), (sx, cy + 0.2)])
px = 10.2                    # pick-off 1 on the way up
mirror(sx, px, 45, 1.1, ACCENT, 4); label(sx + 0.7, px - 0.45, "pick-off 1 (10-30 %)", fs=11, color=GRAY)
beam([(sx, px), (12.5, px)])
cam = box(sx - 2.2, cy, 4.4, 1.4, "pupil camera: COPHI / IFO\n(global shutter or lock-in)", fs=12.5)
label(sx - 2.2, cy - 0.35, "same pixels for COPHI and the IFO -> one registration", fs=10.5, color=GRAY)
# the pick-off line feeds a second splitter: ZWFS arm and the pinhole / P-SRI arm
qx = 12.5
mirror(qx, px, 45, 1.0, ACCENT, 4); label(qx + 0.6, px - 0.45, "pick-off 2 (50 %)", fs=11, color=GRAY)
beam([(qx, px), (qx, 12.9)]); beam([(qx, px), (17.4, px)])
zl = box(qx - 1.0, 10.9, 2.0, 0.7, "focusing lens", fc=LIGHT, tc=INK, ec=ACCENT, fs=10.5)
zd = box(qx - 1.4, 12.0, 2.8, 0.8, "Zernike dimple / metasurface\n(focal plane)", fc=ORANGE, fs=10.5)
zc = box(qx + 2.0, 12.6, 3.0, 1.0, "ZWFS pupil camera\n(stepped or vector reading)", fc=ACCENT, fs=11)
beam([(qx, 12.9), (qx, 13.3), (qx + 2.0, 13.3)], w=0.16)
ph = box(14.6, px - 0.45, 2.5, 0.9, "pinhole / P-SRI\n(focus; shutter)", fc=ORANGE, fs=10.5)
pc = box(17.4, px - 0.6, 3.4, 1.2, "pinhole / P-SRI camera\n(two arms, stepped)", fc=ACCENT, fs=11)

# ---- notes
ax.add_patch(Rectangle((0.4, 0.3), 27.2, 0.0, fc="none", ec="none"))
label(13.6, 4.9, "Which sensor sees what", fs=12.5, bold=True, color=ACCENT)
rows = [("COPHI", "LO on, flat out: test + LO beat on the pupil camera"),
        ("IFO (record)", "flat in, LO shuttered: piezo four-step on the same camera"),
        ("ZWFS", "pick-off, test beam alone, LO shuttered or at a far beat"),
        ("pinhole / P-SRI", "pick-off, test beam alone; its own two-arm stepping"),
        ("light", "pick-offs cost the pupil camera 10-30 %: s_pm x 1.1-1.4"),
        ("registration", "each pick-off camera has its own; COPHI and IFO share")]
for i, (a, b) in enumerate(rows):
    label(13.6, 4.2 - 0.66 * i, a, fs=11, bold=True); label(16.7, 4.2 - 0.66 * i, b, fs=10.5, color=GRAY)
label(0.6, 0.5, "Upstream of the splitter everything is the record bench: the LO collimator sits where the reference flat was;\nthe flip flat puts the flat back without moving the LO.", fs=10.5, color=NOTE, style="italic")
fig.savefig("figs/cophi_multisensor.png", dpi=200, facecolor="white", bbox_inches="tight", pad_inches=0.08)
print("wrote figs/cophi_multisensor.png")
