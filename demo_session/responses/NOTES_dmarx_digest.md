# David Marx's CTB-deck review — digest + answers from the record (2026-09-02)

Attachments: `deck_ctb_Redding_20260826_dmarx_notes.pptx` (his notes are in
the per-slide Notes panes) + `Jorge_VVC_SPIE_2024.pdf` (VVC non-idealities —
guidance for the field-estimation/knowledge-error phase).

## 1. His questions the committed record already answers

- **"Four-surface mask block — a 4-f relay?" (his slide 2)** No: it is
  flat-return at the mask plane + exit-pupil sphere (first half-propagation)
  + the mask plane (second half) + the same sphere again, both sphere
  distances identical to all digits — the condition that makes the block
  transparent with no mask applied.  **"Geometric bench == CODE V?"** No —
  it is the MACOS ray-traced bench model (CTB mech arm); chief-ray
  intersections match it to 1e-13 at every optic.
- **"Why 880 actuators?" (his slide 10 = deck slide 9)** Criterion is in
  `ctb_dm.m`: actuator CENTERS within beam_d/2 + ONE PITCH margin
  (pitch 0.666 mm = beam/32, Gaussian IF, 12% coupling) → 880 active per
  DM.  His guess (influence reaches the pupil) is the right reading; the
  margin is one pitch.  Worth a criterion note on the slide, as he asks.
- **"Control matrix measured, not modeled"** = engine-measured
  differential-poke Jacobian (`ctb_dm_jacobian`, cached + fingerprinted),
  not an analytic jk-propagation derivative.  And yes — his reading is
  exactly right: in those runs **control model == truth model**.
- **"Hard occulter: expected >1e-6 static, <<1e-9 after EFC"** Record
  (CTB_PROP_STATUS S10): static 2.93e-7 → 8.06e-9 fixed-G → 3.82e-9
  relinearized (Lyot 0.50, mono, 10 nm strokes).  Static is 3× under his
  1e-6 expectation (open Lyot, ideal plate, pre-aligned bench); the EFC
  floor is honest-but-shallow at fixed G — the vortex chain is the one
  that goes to the floor (6.8e-15 at half-nm strokes).
- **Mono vs broadband (his slides 4/10/11)** His eye is right: those
  figures are monochromatic; the band study is the later section, and he
  found it himself ("Ah, this answers my earlier comment about the
  magical FPM" — the scalar vortex's 5% chromatic floor).
- **Band averaging (his slide 13)** `ctb_bandpass`: nwf equally spaced
  wavelengths ACROSS the band INCLUDING the edges
  (`wvl0*(1+bf*linspace(-.5,.5,nwf))`), incoherent intensity sum, band-mean
  contrast normalized by the band-MEAN peak (not summed).  EFC band runs
  used 3 λ.  His CGI/Krist refinement — three 3.3% subbands, trapezoid
  within each, mean of subband images, control per subband — is noted for
  the replicate-reality phase (below).

## 2. The framing question: truth model vs control model

His diagram: only voltages go up, only intensities come down; every
top/bottom difference is a knowledge error.  **Answer: we are NOT using
FALCO as the control model — today MACOS is BOTH truth and control,
deliberately**: the campaign isolates the optical plant's controllable
floor from estimation/knowledge effects (his own note names this regime:
"optimal regularization... no model/testbed knowledge error" — exactly our
per-iteration α line search).  FALCO integration is on the post-demo
queue.  One experiment already crosses his line: **CF5b (e2e6m)** — truth
DRIFTS (10 nm/hr JWST-like) while control holds the FROZEN dug-state
Jacobian: a pure model-staleness knowledge error; holds ~2e-9 for 24 h.
Adopt his truth/control vocabulary in future decks and label each run's
regime — cheap and it answers his request directly.

## 3. New guidance to fold into the NEXT phase (aberrations/drifts/estimation)

1. **Jorge's VVC SPIE 2024 paper** — non-ideality budget guidance when
   field estimation + knowledge errors arrive.  Filed here.
2. **Subband structure** — replicate the NKT/Varia reality: 3×3.3%
   subbands, edge-including quadrature within each (Krist/CGISIM
   pattern), control per subband.
3. **HCIT control-strategy references** — add when control ≠ truth.
4. **Hybrid-Lyot coronagraph as the strongest EFC stress test** (large rms
   DM strokes) — connects directly to our stroke-economy machinery
   (push/bump + la stroke curves); a natural CTB family addition.
5. **His offer: PROPER reference-radius phase recipe** for fair
   phase comparison — ACCEPT; relevant to the phase-factor export
   (`BRIEF_ctb_phase_export`).
6. **Pure-PROPER slide (his slide 9)**: his two questions (second bullet
   wording; why the shipped-contrast curve is flat/no Airy rings) need a
   figure-side look — OPEN, not yet answered from the record.

## 4. For TODAY's Keysight talk

- Terminology bridge for the EFC ladder slides: HCIT calls our
  per-iteration regularization line search **"optimal regularization"** —
  usable verbatim in the room.
- If the truth-vs-control question comes up at slide 32: the honest line
  is "control == truth by construction in these runs — that isolates the
  plant's controllable floor; the knowledge-error program (estimation,
  gain maps, FALCO cross-model) is the next phase, and CF5b's frozen-G
  drift hold is its first data point."
