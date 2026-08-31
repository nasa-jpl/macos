#!/usr/bin/env python3
"""Slide-14 schematic: Twyman-Green vs Mach-Zehnder topology, plus the
plate-vs-cube splitter choice.

Build:  ~/dev/MACOS_resources/pymacos/.venv/bin/python make_ifo_schematic_fig.py
        ->  figs/fig_ifo_topologies.png
Palette matches make_brief_slides.py (accent 1F4E79, ink 282828).
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow, Rectangle, Polygon, Circle

ACCENT = "#1F4E79"
INK = "#282828"
GRAY = "#6E6E6E"
LIGHT = "#ECF1F6"
MID = "#D9E2EC"

fig, axes = plt.subplots(1, 3, figsize=(12.6, 2.75),
                         gridspec_kw={"width_ratios": [1, 1, 1.05]})
for ax in axes:
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 7.6)
    ax.set_aspect("auto")
    ax.axis("off")


def beam(ax, x0, y0, x1, y1, both=False, dashed=False, color=ACCENT, lw=1.6):
    ls = (0, (4, 2)) if dashed else "-"
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                arrowprops=dict(arrowstyle="<->" if both else "->",
                                color=color, lw=lw, linestyle=ls,
                                shrinkA=0, shrinkB=0))


def label(ax, x, y, s, size=8, color=GRAY, ha="center", bold=False,
          va="center"):
    ax.text(x, y, s, size=size, color=color, ha=ha, va=va,
            fontweight="bold" if bold else "normal")


# ------------------------------------------------------- Twyman-Green
ax = axes[0]
label(ax, 0.2, 7.15, "Twyman–Green — the demo rig", size=10,
      color=ACCENT, ha="left", bold=True)
# source
ax.add_patch(Circle((1.0, 3.5), 0.16, fc=INK, ec=INK))
label(ax, 1.0, 2.85, "laser")
beam(ax, 1.3, 3.5, 3.75, 3.5)
# PBS cube
ax.add_patch(Rectangle((3.8, 2.8), 1.4, 1.4, fc=LIGHT, ec=INK, lw=1.0))
ax.plot([3.8, 5.2], [2.8, 4.2], color=INK, lw=1.0)
label(ax, 3.0, 2.45, "PBS", size=8)
# reference arm (s), double pass
beam(ax, 4.5, 4.25, 4.5, 5.75, both=True)
ax.plot([3.85, 5.15], [5.95, 5.95], color=INK, lw=3.2)
label(ax, 5.45, 6.0, "reference flat", ha="left")
label(ax, 4.15, 5.15, "s", size=8, color=ACCENT)
# test arm (p), double pass
beam(ax, 5.25, 3.5, 8.1, 3.5, both=True)
th = np.linspace(0, 2 * np.pi, 60)
yy = np.linspace(2.7, 4.3, 60)
ax.plot(8.35 + 0.09 * np.sin(np.linspace(0, 2 * np.pi, 60)), yy,
        color=INK, lw=3.0)
label(ax, 9.35, 3.5, "DM", size=8.5, ha="left", color=INK)
label(ax, 6.6, 3.95, "p", size=8, color=ACCENT)
# output leg
beam(ax, 4.5, 2.75, 4.5, 2.35)
ax.add_patch(Rectangle((4.05, 2.05), 0.9, 0.22, fc=MID, ec=INK, lw=0.8))
label(ax, 5.25, 2.16, "λ/4", size=7.5, ha="left")
beam(ax, 4.5, 1.98, 4.5, 1.78)
ax.add_patch(Rectangle((4.05, 1.5), 0.9, 0.22, fc=MID, ec=INK, lw=0.8))
label(ax, 5.25, 1.61, "analyzer θ", size=7.5, ha="left")
ax.add_patch(Rectangle((4.12, 0.82), 0.76, 0.5, fc=INK, ec=INK))
label(ax, 5.25, 1.07, "camera", size=7.5, ha="left")
label(ax, 0.2, 0.22, "normal incidence, double pass: height counts twice · natural null",
      size=7.2, ha="left")

# -------------------------------------------------------- Mach-Zehnder
ax = axes[1]
label(ax, 0.2, 7.15, "Mach–Zehnder — the alternative", size=10,
      color=ACCENT, ha="left", bold=True)
# splitters and folds
for (cx, cy) in [(2.6, 2.2), (8.4, 5.4)]:      # BS1, BS2
    ax.add_patch(Rectangle((cx - 0.55, cy - 0.55), 1.1, 1.1, fc=LIGHT,
                           ec=INK, lw=1.0))
    ax.plot([cx - 0.55, cx + 0.55], [cy - 0.55, cy + 0.55], color=INK, lw=1.0)
ax.plot([2.05, 3.15], [4.85, 5.95], color=INK, lw=3.0)   # fold upper-left
ax.plot([7.85, 8.95], [1.65, 2.75], color=INK, lw=3.0)   # fold lower-right
label(ax, 2.6, 1.25, "BS1", size=8)
label(ax, 8.4, 4.45, "BS2", size=8)
label(ax, 1.35, 5.85, "fold", size=8)
label(ax, 9.3, 2.2, "fold", size=8, ha="left")
# source + beams (single pass)
ax.add_patch(Circle((0.8, 2.2), 0.16, fc=INK, ec=INK))
beam(ax, 1.1, 2.2, 2.0, 2.2)
beam(ax, 2.6, 2.8, 2.6, 4.8)          # up to fold
beam(ax, 3.2, 5.4, 7.8, 5.4)          # top arm to BS2
beam(ax, 3.2, 2.2, 7.8, 2.2)          # bottom arm to fold
beam(ax, 8.4, 2.8, 8.4, 4.8)          # up to BS2
# sample in one arm
ax.add_patch(Rectangle((5.0, 1.85), 0.9, 0.7, fc=MID, ec=INK, lw=0.8))
label(ax, 5.45, 1.2, "sample", size=7.5)
# two ports
beam(ax, 9.0, 5.4, 10.6, 5.4)
beam(ax, 8.4, 6.0, 8.4, 7.0)
label(ax, 10.8, 5.4, "port A", size=7.5, ha="left")
label(ax, 8.75, 6.9, "port B", size=7.5, ha="left")
label(ax, 0.2, 0.35, "single pass · isolated arms · two ports — dynamics and transmission, not figure",
      size=7.2, ha="left")

# --------------------------------------------- splitter: plate vs cube
ax = axes[2]
label(ax, 0.2, 7.15, "The splitter decides the error budget", size=10,
      color=ACCENT, ha="left", bold=True)
# plate (v1)
c, s = 2.7, 3.9
ax.add_patch(Polygon([(c - 1.0, s - 1.0 - 0.18), (c + 1.0, s + 1.0 - 0.18),
                      (c + 1.0, s + 1.0 + 0.18), (c - 1.0, s - 1.0 + 0.18)],
                     closed=True, fc=LIGHT, ec=INK, lw=1.0))
beam(ax, 0.5, s, c - 0.35, s)
beam(ax, c + 0.2, s, 4.9, s)                          # transmitted
beam(ax, c, s + 0.25, c, 6.1)                          # reflected
beam(ax, c + 0.5, s + 0.35, c + 0.5, 6.1, dashed=True, lw=1.1)  # ghost
label(ax, c + 1.15, 6.0, "ghost", size=7)
label(ax, 2.7, 1.95, "plate (v1)", size=8.5, color=INK, bold=True)
label(ax, 2.7, 0.9, "bare 45° coating: R$_p$ 2.1% → arm\nrotates 7.5°, gauge reads 11.7% high",
      size=7.2)
# cube (v2)
cx, cy, h = 8.9, 3.9, 1.25
ax.add_patch(Polygon([(cx - h, cy - h), (cx + h, cy - h), (cx + h, cy + h)],
                     closed=True, fc=LIGHT, ec=INK, lw=1.0))
ax.add_patch(Polygon([(cx - h, cy - h), (cx - h, cy + h), (cx + h, cy + h)],
                     closed=True, fc=MID, ec=INK, lw=1.0))
ax.plot([cx - h, cx + h], [cy - h, cy + h], color=INK, lw=1.6)
beam(ax, 6.6, cy, cx - h - 0.1, cy)
beam(ax, cx + h + 0.1, cy, 11.4, cy)
beam(ax, cx, cy + h + 0.1, cx, 6.1)
label(ax, 8.9, 1.95, "MacNeille cube (v2)", size=8.5, color=INK, bold=True)
label(ax, 8.9, 0.9, "cemented coated diagonal: the arms\nride its eigenaxes — rotation ≈ 0",
      size=7.2)

fig.subplots_adjust(left=0.01, right=0.99, top=0.98, bottom=0.02,
                    wspace=0.10)
fig.savefig("figs/fig_ifo_topologies.png", dpi=200, facecolor="white")
print("wrote figs/fig_ifo_topologies.png")
