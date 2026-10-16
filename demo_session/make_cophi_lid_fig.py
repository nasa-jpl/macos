"""deck_cophi: 'How the lock-in detector works' -- one pixel's signal flow (left) and the
algebra of the signals (right, matplotlib mathtext).  Same palette and helpers as
make_cophi_block_fig.py (copied, so importing this does not redraw the block diagrams).
Writes figs/cophi_lid.png (+ .svg).  Dave 2026-10-16."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

ACCENT = "#1F4E79"; INK = "#282828"; GRAY = "#6E6E6E"; LIGHT = "#ECF1F6"
RF = "#1F4E79"; BEAM = "#2E8B57"; DATA = "#555555"; CHARGE = "#B9770E"; GREEN = "#7DBB3C"; NOTE = "#8A6D1A"

def box(ax, x, y, w, h, text, fc=ACCENT, tc="white", fs=11, bold=False, ec=None, lw=1.2, ls="-", r=0.25):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={r}", fc=fc, ec=ec or fc, lw=lw, linestyle=ls, zorder=3))
    ax.text(x + w/2, y + h/2, text, ha="center", va="center", color=tc, fontsize=fs, fontweight="bold" if bold else "normal", zorder=4, linespacing=1.15)
    return (x, y, w, h)
def port(b, side, f=0.5):
    x, y, w, h = b
    return {"l": (x, y + h*f), "r": (x + w, y + h*f), "t": (x + w*f, y + h), "b": (x + w*f, y)}[side]
def wire(ax, pts, color=RF, lw=1.8, arrow=True, ls="-", z=2):
    for a, b in zip(pts[:-1], pts[1:]):
        ax.plot([a[0], b[0]], [a[1], b[1]], color=color, lw=lw, ls=ls, zorder=z, solid_capstyle="round")
    if arrow:
        a, b = pts[-2], pts[-1]
        ax.annotate("", xy=b, xytext=a, arrowprops=dict(arrowstyle="-|>", color=color, lw=lw, ls=ls, shrinkA=0, shrinkB=0, mutation_scale=14), zorder=z)
def label(ax, x, y, s, color=INK, fs=10, ha="left", va="center", bold=False):
    ax.text(x, y, s, color=color, fontsize=fs, ha=ha, va=va, fontweight="bold" if bold else "normal", zorder=5, linespacing=1.15)
def note(ax, x, y, s, fs=9.5, ha="left"):
    ax.text(x, y, s, color=NOTE, fontsize=fs, ha=ha, va="center", fontstyle="italic", zorder=5,
            bbox=dict(boxstyle="round,pad=0.25", fc="#FFF8E1", ec="#E0C97A", lw=0.8))

W, H = 32.0, 15.0
fig = plt.figure(figsize=(W/1.6, H/1.6)); ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W); ax.set_ylim(0, H); ax.axis("off")
ax.text(0.5, H-0.55, "One lock-in pixel: the beat is demodulated by SORTING the photoelectrons, before any readout", color=ACCENT, fontsize=14, fontweight="bold", va="center")
ax.plot([17.3, 17.3], [0.6, H-1.1], color="#B8C6D4", lw=1)

# ---- left: the signal flow -------------------------------------------------
lt = box(ax, 0.6, 11.3, 4.2, 1.3, "test beam  $E_t\\,e^{i\\varphi}$", fc=BEAM, fs=11.5)
lo = box(ax, 0.6, 9.4, 4.2, 1.3, "LO  $E_{LO}\\,e^{i\\omega_b t}$\n(fiber, at the flat's seat)", fc="#8E44AD", fs=11)
pd = box(ax, 6.0, 9.9, 4.6, 2.2, "photodiode\n$I(t)=I_0\\,[1+\\gamma\\cos(\\varphi-\\omega_b t)]$\nelectrons at the rate $I(t)$", fc=ACCENT, fs=11)
wire(ax, [port(lt, "r"), (5.4, port(lt, "r")[1]), (5.4, port(pd, "l", 0.7)[1]), port(pd, "l", 0.7)], color=BEAM)
wire(ax, [port(lo, "r"), (5.4, port(lo, "r")[1]), (5.4, port(pd, "l", 0.3)[1]), port(pd, "l", 0.3)], color="#8E44AD")
# reference chain
clk = box(ax, 0.6, 6.6, 4.2, 1.3, "Rb clock → synthesizers\n$f_1$, $f_2 = f_1 + f_b$ (the AOMs)", fc=RF, fs=10.5)
mix = box(ax, 0.6, 4.6, 4.2, 1.3, "mixer / divider → $f_b$\nsquare wave $s_I(t)=\\mathrm{sgn}\\,\\cos\\omega_b t$", fc=RF, fs=10.5)
sh = box(ax, 0.6, 2.6, 4.2, 1.3, "quarter-period delay\n$s_Q(t)=\\mathrm{sgn}\\,\\sin\\omega_b t$", fc=RF, fs=10.5)
wire(ax, [port(clk, "b"), port(mix, "t")], color=RF); wire(ax, [port(mix, "b"), port(sh, "t")], color=RF)
# gates and storage
gI = box(ax, 6.0, 6.3, 4.6, 1.6, "transfer gates A | B\ndriven by $s_I$ (0° / 180°)", fc=CHARGE, fs=11)
gQ = box(ax, 6.0, 3.6, 4.6, 1.6, "transfer gates C | D\ndriven by $s_Q$ (90° / 270°)", fc=CHARGE, fs=11)
wire(ax, [port(pd, "b", 0.35), (port(pd, "b", 0.35)[0], 8.6), (port(gI, "t", 0.35)[0], 8.6), port(gI, "t", 0.35)], color=CHARGE, lw=2.4)
wire(ax, [port(pd, "b", 0.65), (port(pd, "b", 0.65)[0], 8.6), (5.7, 8.6), (5.7, port(gQ, "l", 0.75)[1]), port(gQ, "l", 0.75)], color=CHARGE, lw=2.4)
wire(ax, [port(mix, "r"), (5.5, port(mix, "r")[1]), (5.5, port(gI, "l")[1]), port(gI, "l")], color=RF)
wire(ax, [port(sh, "r"), (5.5, port(sh, "r")[1]), (5.5, port(gQ, "l")[1]), port(gQ, "l")], color=RF)
nA = box(ax, 11.4, 6.55, 1.6, 1.1, "node A", fc=LIGHT, tc=INK, fs=10.5, ec=CHARGE); nB = box(ax, 13.3, 6.55, 1.6, 1.1, "node B", fc=LIGHT, tc=INK, fs=10.5, ec=CHARGE)
nC = box(ax, 11.4, 3.85, 1.6, 1.1, "node C", fc=LIGHT, tc=INK, fs=10.5, ec=CHARGE); nD = box(ax, 13.3, 3.85, 1.6, 1.1, "node D", fc=LIGHT, tc=INK, fs=10.5, ec=CHARGE)
wire(ax, [port(gI, "r", 0.65), port(nA, "l")], color=CHARGE, lw=2.0); wire(ax, [port(gI, "r", 0.35), (11.1, port(gI, "r", 0.35)[1]), (11.1, 6.2), (14.1, 6.2), port(nB, "b")], color=CHARGE, lw=2.0)
wire(ax, [port(gQ, "r", 0.65), port(nC, "l")], color=CHARGE, lw=2.0); wire(ax, [port(gQ, "r", 0.35), (11.1, port(gQ, "r", 0.35)[1]), (11.1, 3.5), (14.1, 3.5), port(nD, "b")], color=CHARGE, lw=2.0)
label(ax, 11.4, 8.0, "nodes integrate N beat cycles\n(one I/Q sample, ~1 ms)", color=GRAY, fs=9.5)
dI = box(ax, 15.3, 6.55, 1.7, 1.1, "A − B\n→ I", fc=DATA, fs=11); dQ = box(ax, 15.3, 3.85, 1.7, 1.1, "C − D\n→ Q", fc=DATA, fs=11)
wire(ax, [port(nB, "r"), port(dI, "l")], color=DATA); wire(ax, [port(nD, "r"), port(dQ, "l")], color=DATA)
out = box(ax, 11.4, 1.0, 5.6, 1.5, "ADC → I, Q per pixel per sample\n$\\varphi=\\mathrm{atan2}(Q, I)$,   $h=\\varphi\\,\\lambda/4\\pi$", fc=GREEN, tc=INK, fs=11)
wire(ax, [port(dI, "b"), (port(dI, "b")[0], 5.6), (17.3-0.6, 5.6)], color=DATA, arrow=False); wire(ax, [port(dQ, "b"), (port(dQ, "b")[0], 2.75), (port(out, "t", 0.8)[0], 2.75), port(out, "t", 0.8)], color=DATA)
wire(ax, [(16.7, 5.6), (16.7, 2.75)], color=DATA, arrow=False)
note(ax, 6.0, 1.75, "the DC $I_0$ lands equally in A and B (and C, D):\nit cancels in the difference -- whether it still\nfills the nodes is the vendor question", fs=9.5)

# ---- right: the algebra ----------------------------------------------------
x0 = 18.0; y = H - 1.5
rows = [
 ("the light at the pixel", r"$I(t)=I_0\,[\,1+\gamma\cos(\varphi-\omega_b t)\,],\qquad \omega_b=2\pi f_b,\ \ f_b=|f_1-f_2|$"),
 ("the gate references (square waves at the beat)", r"$s_I(t)=\mathrm{sgn}\,\cos\omega_b t,\qquad s_Q(t)=\mathrm{sgn}\,\sin\omega_b t$"),
 ("what the nodes hold after N cycles ($\\tau = N/f_b$, $N_e=I_0\\tau$ electrons)", r"$A-B=\int_0^{\tau} I\,s_I\,dt=\frac{2}{\pi}\,N_e\,\gamma\cos\varphi,\qquad C-D=\frac{2}{\pi}\,N_e\,\gamma\sin\varphi$"),
 ("the DC is common to both nodes and cancels", r"$A+B=C+D=N_e\quad(\mathrm{the\ well\ sees\ } N_e \mathrm{\ unless\ the\ DC\ is\ drained})$"),
 ("the phase and the surface", r"$\varphi=\mathrm{atan2}(C-D,\ A-B)=\mathrm{atan2}(Q,I),\qquad h=\varphi\,\lambda/4\pi$"),
 ("shot + read noise per pixel per sample (the ×1.11 of a square-wave multiplier)", r"$\sigma_\varphi^2=\frac{\pi^2}{8}\cdot\frac{2\,(N_e+r^2)}{N_e^2\,\gamma^2}\ \ \Rightarrow\ \ \sigma_\varphi=\frac{\pi}{2\gamma\sqrt{N_e}}\ \ (r=0)$"),
 ("what a square wave also admits", r"$s_I=\frac{4}{\pi}\left[\cos\omega_b t-\frac{1}{3}\cos 3\omega_b t+\frac{1}{5}\cos 5\omega_b t-\dots\right]$: a clean two-frequency beat has no harmonics"),
]
for i, (cap, eq) in enumerate(rows):
    ax.text(x0, y, f"({i+1})  {cap}", color=ACCENT, fontsize=11, va="top", ha="left")
    ax.text(x0 + 0.3, y - 0.62, eq, color=INK, fontsize=13.5, va="top", ha="left")
    y -= 1.95
fig.savefig("figs/cophi_lid.png", dpi=200, facecolor="white"); fig.savefig("figs/cophi_lid.svg", facecolor="white"); print("wrote figs/cophi_lid.png")
