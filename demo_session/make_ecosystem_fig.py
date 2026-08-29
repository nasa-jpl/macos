#!/usr/bin/env python3
"""Slide-3 ecosystem diagram: users -> surfaces -> engine, feeders,
validators, and engine outputs.

Build:  ~/dev/MACOS_resources/pymacos/.venv/bin/python make_ecosystem_fig.py
        ->  figs/fig_ecosystem.png
Palette matches make_brief_slides.py (accent 1F4E79, ink 282828).
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

ACCENT = "#1F4E79"
INK = "#282828"
GRAY = "#6E6E6E"
LIGHT = "#ECF1F6"
RULE = "#B8C6D4"

W, H = 15.2, 6.9
fig, ax = plt.subplots(figsize=(W, H), dpi=200)
ax.set_xlim(0, W); ax.set_ylim(0, H); ax.axis("off")


def box(x, y, w, h, head, sub=None, fill=LIGHT, edge=ACCENT, dashed=False,
        headcolor=None, headsize=13, subsize=10):
    fb = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.08",
                        linewidth=1.4, edgecolor=edge, facecolor=fill,
                        linestyle=(0, (4, 2)) if dashed else "solid")
    ax.add_patch(fb)
    hc = headcolor or (ACCENT if fill != ACCENT else "white")
    sc = GRAY if fill != ACCENT else "#D9E4EF"
    if sub:
        ax.text(x + w/2, y + h*0.62, head, ha="center", va="center",
                fontsize=headsize, fontweight="bold", color=hc)
        ax.text(x + w/2, y + h*0.28, sub, ha="center", va="center",
                fontsize=subsize, color=sc)
    else:
        ax.text(x + w/2, y + h/2, head, ha="center", va="center",
                fontsize=headsize, fontweight="bold", color=hc)


def arrow(x1, y1, x2, y2, label=None, lx=0.12, ly=0.0, both=False,
          color=INK, lw=1.3):
    ar = FancyArrowPatch((x1, y1), (x2, y2),
                         arrowstyle="<|-|>" if both else "-|>",
                         mutation_scale=14, linewidth=lw, color=color,
                         shrinkA=2, shrinkB=2)
    ax.add_patch(ar)
    if label:
        ax.text((x1+x2)/2 + lx, (y1+y2)/2 + ly, label, fontsize=9.5,
                color=GRAY, ha="left", va="center")


# ---- users band ----------------------------------------------------------
box(0.7, 5.95, 10.6, 0.72,
    "Users — optical engineers, and AI agents driving the same interfaces",
    fill="white", edge=RULE, headcolor=INK, headsize=13)

# ---- interface row -------------------------------------------------------
box(0.7, 4.45, 2.9, 0.95, "macos", "interactive command line")
box(4.35, 4.45, 3.25, 0.95, "mmacos", "analysis · design · simulation")
box(8.35, 4.45, 2.95, 0.95, "pymacos", "scripting · regression")

arrow(2.15, 5.95, 2.15, 5.42, "terminal")
arrow(5.97, 5.95, 5.97, 5.42, "MATLAB")
arrow(9.82, 5.95, 9.82, 5.42, "Python")

# ---- smacos library layer ------------------------------------------------
box(4.35, 3.15, 6.95, 0.78, "smacos — linkable engine library",
    "one shared wrapper layer serves every binding")
arrow(5.97, 4.45, 5.97, 3.95)
arrow(9.82, 4.45, 9.82, 3.95)

# ---- engine core ---------------------------------------------------------
box(0.7, 1.78, 10.6, 0.95,
    "MACOS ENGINE — Fortran ray trace + physical optics",
    "one prescription (Rx) language", fill=ACCENT, edge=ACCENT,
    headsize=15, subsize=10.5)
arrow(2.15, 4.45, 2.15, 2.75)          # CLI is the engine's own executable
arrow(7.8, 3.15, 7.8, 2.75)

# ---- feeders + validators (bottom row) -----------------------------------
box(0.7, 0.25, 2.4, 0.9, "SegMirMaker", "segmented-mirror\nprescriptions")
box(3.45, 0.25, 2.4, 0.9, "CODE V converter", "→ MACOS Rx", subsize=10)
box(6.2, 0.25, 2.4, 0.9, "CODE V", "reference\ndesign code",
    fill="white", edge=GRAY, dashed=True, headcolor=INK)
box(8.95, 0.25, 2.35, 0.9, "PROPER", "independent\ndiffraction code",
    fill="white", edge=GRAY, dashed=True, headcolor=INK)

arrow(1.90, 1.15, 1.90, 1.76, "Rx", both=True)
arrow(4.65, 1.15, 4.65, 1.76, "Rx")
arrow(6.20, 0.85, 5.85, 0.85)
ax.text(6.025, 0.48, ".seq", fontsize=9, color=GRAY,
        ha="center", va="center")
arrow(7.40, 1.15, 7.40, 1.76, "cross-validation\n6,601 tests", both=True)
arrow(10.12, 1.15, 10.12, 1.76, "cross-validation\n10⁻¹¹–10⁻¹³", both=True)

# ---- engine outputs (right column) ---------------------------------------
ax.text(13.55, 5.85, "outputs", fontsize=12, fontweight="bold",
        color=INK, ha="center", va="center")
box(12.35, 4.85, 2.6, 0.72, "files", "Rx · wavefront · image data",
    fill="white", edge=RULE, headcolor=INK, headsize=12, subsize=9.5)
box(12.35, 3.85, 2.6, 0.72, "screen", "plots · layouts · fringes",
    fill="white", edge=RULE, headcolor=INK, headsize=12, subsize=9.5)
box(12.35, 2.85, 2.6, 0.72, "MATLAB", "workspace arrays",
    fill="white", edge=RULE, headcolor=INK, headsize=12, subsize=9.5)
box(12.35, 1.85, 2.6, 0.72, "Python", "NumPy arrays",
    fill="white", edge=RULE, headcolor=INK, headsize=12, subsize=9.5)

for ytgt in (5.21, 4.21, 3.21, 2.21):
    arrow(11.32, 2.45, 12.33, ytgt, color=GRAY, lw=1.1)

fig.savefig("figs/fig_ecosystem.png", bbox_inches="tight",
            facecolor="white")
print("wrote figs/fig_ecosystem.png")
