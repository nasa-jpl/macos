#!/usr/bin/env python3
"""Slide-3 mmacos directory map: categories of code with major entries.

Build:  ~/dev/MACOS_resources/pymacos/.venv/bin/python make_mmacos_fig.py
        ->  figs/fig_mmacos.png
Palette matches make_brief_slides.py (accent 1F4E79, ink 282828).
Counts as of 2026-08 (168 +macos functions, 82 test classes).
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from matplotlib.lines import Line2D

ACCENT = "#1F4E79"
INK = "#282828"
GRAY = "#6E6E6E"
LIGHT = "#ECF1F6"
RULE = "#B8C6D4"

W, H = 17.6, 6.6
fig, ax = plt.subplots(figsize=(W, H), dpi=200)
ax.set_xlim(0, W); ax.set_ylim(0, H); ax.axis("off")

# ---- root band -----------------------------------------------------------
root = FancyBboxPatch((0.5, 5.8), 16.6, 0.62,
                      boxstyle="round,pad=0.02,rounding_size=0.08",
                      linewidth=1.6, edgecolor=ACCENT, facecolor=ACCENT)
ax.add_patch(root)
ax.text(8.8, 6.11, "mmacos/  —  the MATLAB toolbox over the MACOS engine",
        ha="center", va="center", fontsize=15, fontweight="bold",
        color="white")

COLS = [
    ("src/+macos/", "the engine veneer — 168 functions", [
        "init · load_rx · modify",
        "trace · opd · perturb",
        "spot · intensity · complex_field",
        "elt_* prescription get/set",
        "pupil_find · fex · xp",
        "jones_pupil · pol elements",
        "view_rx · view_std · draw_rays",
        "compose · segment_grid_basis",
    ]),
    ("+macos/ sensitivities", "the linear-model channels", [
        "dw_dx_multi — rigid body",
        "dw_dz_zernike_multi — figure",
        "dw_dsurf_multi — surface",
        "dw_dgrid_multi — grid pokes",
        "dmet_dx — metrology",
        "each per element / segment,",
        "stacked over fields and",
        "configs:  w = J·x + w₀",
    ]),
    ("src/+macos/+design/", "the design layer", [
        "System · Telescope · Bench",
        "tma_layout · offner_layout",
        "sz_tma · freeform_unobscured",
        "segment_rx · seg_apertures",
        "add_met · edge_sensors",
        "met_layout_opt · met_view",
        "twyman_green · pbs_macneille",
    ]),
    ("design/runners/", "one-call stage runners", [
        "run_sensitivities",
        "run_segmentation",
        "run_met",
        "run_compare",
        "run_simulator",
        "presets/ — saved winners",
        "checkpoint / resume built in",
    ]),
    ("templates/ · tests/", "worked examples · regression", [
        "10_telescopes …",
        "   … 90_polarization",
        "challenges/ — Rodgers studies",
        "every stage committed: the .in,",
        "   the report, the figures",
        "tests/: 82 classes,",
        "   600+ regression checks",
    ]),
]

cw, gap, x0 = 3.24, 0.10, 0.5
for i, (path, sub, items) in enumerate(COLS):
    x = x0 + i * (cw + gap)
    bx = FancyBboxPatch((x, 0.35), cw, 4.9,
                        boxstyle="round,pad=0.02,rounding_size=0.08",
                        linewidth=1.4, edgecolor=ACCENT, facecolor=LIGHT)
    ax.add_patch(bx)
    ax.plot([x + cw/2, x + cw/2], [5.25, 5.8], color=GRAY, lw=1.2)
    ax.text(x + cw/2, 4.93, path, ha="center", va="center",
            fontsize=10.5, fontweight="bold", color=ACCENT,
            family="monospace")
    ax.text(x + cw/2, 4.58, sub, ha="center", va="center",
            fontsize=9.5, color=GRAY, style="italic")
    ax.plot([x + 0.18, x + cw - 0.18], [4.38, 4.38], color=RULE, lw=1.0)
    for j, it in enumerate(items):
        ax.text(x + 0.16, 4.10 - j * 0.46, it, ha="left", va="center",
                fontsize=9.6, color=INK)

fig.savefig("figs/fig_mmacos.png", bbox_inches="tight", facecolor="white")
print("wrote figs/fig_mmacos.png")
