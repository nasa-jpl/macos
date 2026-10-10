# NOTE to CCL — a plain-physics pupil-blur demo for the gauges story

From CCMac (Claude Code on Dave's Mac), 2026-10-09, after Dave's DM-gauge
talk (2026-10-08).  Dave asked me to put this into the resources tree as
part of the gauges story and leave you a note.  You own the gauge deck and
the lane reports, so the integration (deck slide, REPORT line) is yours if
you want it — this note hands you a finished, pushed tool and flags what is
open.

## What it is, and why

**Fang Shi is skeptical of the pupil-imaging analysis.**  `tg96_pupilsim`
answers his concern rigorously (its header even cites him) — but it is a
three-stage machine: engine-traced zone PSFs, the DM field through them, a
plane-to-plane Fourier cross-check.  Dave's instinct was that a skeptic
needs a version he can check by eye.  So this is the plain-physics
companion, **engine-free**, in four visible steps:

1. build a KNOWN DM surface from known commands (the `dm_influence_map`
   Gaussian-influence model),
2. BLUR it with a pupil-imaging kernel (a Gaussian PSF, 1/e radius sigma,
   swept) and add read noise — what the camera sees,
3. RECONSTRUCT the commands with the gauges' **own** reconstruction
   (`dmg_stencil` + `dmg_act_fit`, the Tikhonov lattice deconvolution — not
   a toy stand-in, so it cannot be waved off as disabling the real path),
4. plot actuator-command error vs blur, NAIVE (blur ignored) vs CALIBRATED
   (blur folded into the stencil), on the checkerboard (actuator Nyquist,
   the hardest pattern) and a random 30 nm working surface.

## Where it lives (resources, `dev-candidate`)

`mmacos/templates/40_benches/tg_psi_dm96_oap/pupil_blur_demo.m`, next to
`tg96_pupilsim`.  Writes `runs/pupil_blur_demo/pupil_blur_demo_{report.txt,
curve.png,maps.png}` (committed as the demo's results).  README updated in
two places: the "Fang Shi's question" run block and the Files table.
Run from that directory: `>> pupil_blur_demo` (~40 s, no mex, no engine).

## The result (so you can quote it)

- blur-free floor (sigma = 0): checker **8.1 %**, random **3.6 %** — this is
  intrinsic Nyquist regularization at lambda = 0.05, NOT blur.  State it as
  the reconstruction's own floor so no one reads it as a blur cost.
- a CALIBRATED read holds near that floor until the blur doubles it
  (checker 2x floor at sigma ~ 0.4 pitch); the NAIVE read doubles it by
  sigma ~ 0.2 pitch — i.e. calibration buys ~2x in tolerated blur.
- the structural knee is MTF(actuator Nyquist) = lambda at **sigma = 0.78
  pitch**; beyond it the checkerboard is lost even to calibration.
- the real tg96 leg (`tg96_pupilsim`) broadens a single poke by <~1 %
  (width/true ~ 1.0) — i.e. sigma ~ 0, the far-left flat part of the curve,
  consistent with the deck's 0.13 % / 0.29 % after-calibration numbers.

**The one line for a skeptic:** pupil blur is just an MTF roll-off of the
DM's high spatial frequencies; a calibrated read undoes it until that MTF at
the DM's highest frequency approaches the reconstruction's noise/regularization
floor.  The built leg is nowhere near that knee.

## Open items (yours to take or leave)

1. **Mark `tg96_pupilsim`'s measured PSF width on the curve** as a labelled
   point — the single strongest move, it ties the simple demo to the rigorous
   one on one axis.  I did not do it so the demo stays engine-free; you have
   the pupilsim run records to read the equivalent sigma from.
2. **Fold a line into `REPORT_gauge_ifo.md` / the deck's pupil discussion**
   if Dave wants the simple curve in the story (it is only in the README now).
3. The sigma = 0 floor depends on lambda (0.05) and the 20 pm read noise —
   both are `pupil_blur_demo` knobs; sweep if you want the SNR dependence.
4. At sigma > 0.8 pitch the CALIBRATED `pcg` does not fully converge (400-iter
   cap, ill-conditioned effective stencil) — immaterial (that regime is broken
   either way) but worth a word if anyone squints at the far right.

## Wider context (the other talk follow-ups)

Dave opened four threads after the talk; this closes the Fang one.  The others:
Kent has ZWFS-reconstruction experience to share (waiting on him); Feng
proposes **COPHI heterodyne** as a 5th candidate (1 Hz beat on a 1k x 1k
detector at 30 Hz -> full-field OPD; refs in `macos_sandbox/DM_gauge/`,
slides pending from him); and a **thorough per-candidate error budget** is the
frame they all feed — the blur row above is one cell of it.

Nothing here depends on an engine change.  Pushed on `dev-candidate`.

— CCMac, 2026-10-09, for CCL
