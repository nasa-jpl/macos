#!/usr/bin/env python3
"""deck_cophi slide 4: the reading's equations as one graphic (matplotlib mathtext; the
TeX source for a typeset version is deck_cophi_equations.tex).
Build: ~/dev/MACOS_resources/pymacos/.venv/bin/python make_cophi_equations_fig.py
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
ACCENT = "#1F4E79"; INK = "#282828"; GRAY = "#6E6E6E"
rows = [
 ("the irradiance during the beat",
  r"$I(t)=|E_t+E_{LO}\,e^{i\alpha(t)}|^2=I_0\,[\,1+\gamma\cos(\varphi-\alpha(t))\,],\quad \varphi=\arg(E_t\,\bar{E}_{LO}),\quad \gamma=\frac{2|E_t||E_{LO}|}{|E_t|^2+|E_{LO}|^2}$"),
 ("the integrating bucket: frame k averages the beat over a window D = 2$\\pi$/N",
  r"$I_k=\frac{1}{D}\int_{\alpha_k}^{\alpha_k+D} I\,d\alpha=I_0\,[\,1+\gamma\,\mathrm{sinc}(D/2)\cos(\varphi-\bar{\alpha}_k)\,],\quad \bar{\alpha}_k=\alpha_k+\frac{D}{2},\quad \mathrm{sinc}(\pi/4)=0.900\ (N=4)$"),
 ("the phase, by least squares over the N buckets (N = 4 is the four-step)",
  r"$I_k=a+b\cos\bar{\alpha}_k+c\sin\bar{\alpha}_k\ \Rightarrow\ \hat{\varphi}=\mathrm{atan2}(c,b)=\mathrm{atan2}(I_2-I_4,\,I_1-I_3)+\bar{\alpha}_1\quad (N=4)$"),
 ("the surface: double pass, about the null frozen on the flat DM",
  r"$h=\frac{\lambda}{4\pi}\,(\hat{\varphi}-\varphi_{\mathrm{null}})$"),
 ("a per-step error $\\epsilon$ on the four-bucket phase (the textbook ripple)",
  r"$\delta\hat{\varphi}\simeq\frac{\epsilon}{2}\sin 2\varphi$"),
 ("photon noise per pixel: $I_0$ e$^-$ per frame, r e$^-$ read noise, g$'$ = $\\gamma$ sinc(D/2)",
  r"$\sigma_\varphi^2=\frac{2\,(I_0+r^2)}{N\,I_0^2\,g'^2},\qquad \sigma_h=\frac{\lambda}{4\pi}\,\sigma_\varphi$"),
 ("two colors on one camera: the synthetic wavelength, its range, its noise",
  r"$\Lambda=\frac{\lambda_1\lambda_2}{|\lambda_1-\lambda_2|},\qquad \mathrm{range}\ \pm\frac{\Lambda}{4}\ \mathrm{of\ surface},\qquad \sigma_{h,\mathrm{coarse}}=\frac{\Lambda}{\lambda}\,\sigma_h$"),
]
fig = plt.figure(figsize=(22, 6.6), dpi=220)
ax = fig.add_axes([0, 0, 1, 1]); ax.axis("off"); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
cols = {0: (0.015, 0.03), 1: (0.515, 0.53)}      # caption x, equation x per column
plan = [(0, 0), (0, 1), (0, 2), (0, 3), (1, 4), (1, 5), (1, 6)]   # (column, row index)
y = {0: 0.96, 1: 0.96}
for c, i in plan:
    cap, eq = rows[i]
    ax.text(cols[c][0], y[c], f"({i+1})  {cap}", color=ACCENT, fontsize=12, va="top", ha="left")
    ax.text(cols[c][1], y[c] - 0.075, eq, color=INK, fontsize=16, va="top", ha="left")
    y[c] -= 0.245 if i not in (3, 4) else 0.21
fig.savefig("figs/cophi_equations.png", facecolor="white")
print("wrote figs/cophi_equations.png")
