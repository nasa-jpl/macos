# MACOS for the CODE V team — session deck (TO independent draft)

**DRAFT — not for release.** Outward-facing; needs Dave's sign-off on
every claim about history, provenance and anyone else's work. Built to
the approved 25-slide structure in `BRIEF_demo_deck.md`; this file is the
*script*, and any `.pptx` must be generated from records, never
hand-typed.

_Written blind, before reading `deck_keysight.md`, for comparison._

---

## Before anything: who is in the room

The CODE V team are not a general audience. They know optical design
better than we do, they have seen every "our optimizer beat yours" deck,
and **Mike Rodgers is one of their own** — the offset-field problems in
Section V are his, and much of the unobscured-reflective design
literature is his too. Three consequences run through this whole deck:

1. **Reproduce before you compare.** Every Rodgers claim leads with "we
   reproduced his design at 0.99–1.00× his own metric", and only then
   says anything about improving on it. If we skip that, the room hears
   a vendor claim and stops listening.
2. **State the metric before the number.** This audience's first
   question is always "RMS of what, referenced to what surface, over
   which field?" Answering it *before* they ask is how the numbers get
   believed. Every quoted number in this deck carries its metric.
3. **Show the negatives.** A deck of wins reads as marketing. The
   retracted result, the unbuildable design, the 595 µm failure and the
   live refusal path are what make the wins credible.

The through-line, one sentence: **an open, scriptable physics engine,
anchored to references this room already trusts, driven by an AI that
shows its work and refuses when it cannot stand behind an answer.**

---

# I. Intro

### 1. Title / agenda

**Claim:** none yet — set expectations.
**On the slide:** title, date, presenters; four-line agenda (what it is
→ what it does → why you should believe it → watch it work).
**Say:** flag the live demo up front so the room knows something real is
coming, and flag that **we launch it now** — the audience spec is taken
in the next few minutes and the solve runs through the session.
**Note:** this is where the beat-22b spec is collected under Dave's
2026-08-28 ruling. See slide 22 for the mechanics.

### 2. The family: one engine, four surfaces

**Claim:** one physics core, four ways in, one prescription language.
**On the slide:** `macos` CLI · `libsmacos` · `mmacos` (MATLAB) ·
`pymacos` (Python), all over `macos_api_mod`; one `.in` Rx language
underneath.
**Say:** the point is not four products — it is that the MATLAB design
layer, the Python regression suite and the interactive CLI are provably
the *same* engine, so a result found in one is reproducible in the
others. That is what makes the test anchors on slide 9 mean anything.

### 3. History — ONE slide

**Claim:** this is not new code; it has been load-bearing on real
programs for three decades.
**On the slide:** timeline — late-1980s NASA technology development →
1991 Hubble recovery → JWST modeling and WFSC development → TMT →
validation against testbed and flight data.
**Say:** short. The credibility this buys is spent on slide 9, not here.
**⚠ DAVE'S SOURCES AND SIGN-OFF. Do not fabricate imagery, dates,
program names or roles.** If a claim is not sourced by Dave, it does not
go on the slide.

---

# II. Capability

### 4. Current applications

**Claim:** system-level modeling, not just lens design.
**On the slide:** error budgeting · integrated simulation · WFSC for
segmented telescopes · coronagraph testbeds and instruments.
**Say:** frames the difference from a design code. CODE V's centre of
gravity is designing the optic; ours is predicting what an *instrument*
does once the optic is real, perturbed, segmented, controlled and
diffracting. Sections IV and V then show we can do the design too.

### 5. mmacos capabilities (two columns)

**Claim:** the library is a toolkit, and the drivers are the product.
**On the slide:**
- *Drivers* — design layer (`System`/`Telescope`/`Bench`);
  `dw_dx`/`dw_dz`/`dw_dsurf`/`dw_dgrid` sensitivity supervisors;
  `run_sensitivities` / `run_met` / `run_compare` / `run_simulator` /
  `run_segmentation`; study orchestrators (`ctb_study`).
- *Utilities* — `trace`/`opd`/`perturb`/`fex`/`pupil_find`;
  `view_rx`/`view_std`/`draw_rays`;
  `jones_pupil`/`polarizer`/`waveplate`; `segment_grid_basis` + grid I/O;
  `compose`; `spot`/`intensity`/`complex_field`.
**Say:** one line only — "everything on the right is a call; everything
on the left is a study you can re-run from one command." Do not read the
list aloud.

### 6. The linear-model workflow

**Claim:** `w = J·x + w₀` is the spine, and the engine's job is an
honest `J`.
**On the slide:** the channels, multi-field × multi-configuration
stacking, groups.
**Say:** this is the shape every error budget, every WFSC loop and every
sensitivity study takes. What differentiates us is not the algebra — it
is that `J` comes from real rays through the real prescription, with
per-ray status, rather than from a paraxial surrogate.

### 7. Streamlining

**Claim:** the pipeline is one reconstructible study, not a folder of
scripts.
**On the slide:** design → segmentation → sensitivities → MET → compare
→ simulate; one-call studies; checkpoint/resume.
**Say:** the reproducibility point — a study can be re-run months later
and land on the same numbers. Bridges to slide 9.

### 8. The design layer, and why it exists

**Claim:** fully public-domain — **no proprietary prescription needed to
reproduce anything in this deck.**
**On the slide:** the design layer stack; the fact that every example
here is built from parameters, not from a licensed Rx.
**Say:** this matters to this room more than it looks. Everything they
are about to see they could rebuild themselves, without asking anyone
for a deck. It is also why the Rodgers work is a fair comparison —
we built his instrument from his published parameters.

---

# III. Trust

### 9. Test suites + anchors — the handshake slide

**Claim:** we do not ask you to trust the engine; we show what it is
pinned against, and we show the pins would catch a break.
**On the slide:** four anchors and one rule.

| anchor | against | agreement |
|---|---|---|
| geometric ray trace | **CODE V** | large comparison suite, all passing |
| physical optics | **PROPER** (Krist) | ~1e-11 … 1e-13 across propagation regimes |
| polarization | textbook closed forms (Born & Wolf; Abelès/Macleod) | ~1e-14 |
| polarization, absolute | **published lab data** (protected-Al Mueller) | ~1e-14 vs the paper's own ±0.01 |
| cross-compiler | ifx vs gfortran | bit-identical |

**THE RULE — say this one out loud:** *every gate is shown to FAIL the
pre-fix engine.* A test that passes against both the broken and the
fixed code is not evidence, and we treat it as a defect in the test.

**Say:** lead with CODE V — it is their reference and it is our
geometric anchor, which reframes the whole session from "competitor" to
"we check ourselves against you." Then note that PROPER and the
published Mueller data cover the regimes CODE V does not.
**Worked example of the rule, if there is time:** a polarization sign
error survived every existing gate for years because a reflection is an
involution (it cancels exactly on mirror pairs) and unitary (invisible
to a unitarity check) — and because the "analytic" reference had been
transcribed from the engine's own expression, which is circular in
exactly the sign it should test. The fix was found by demanding a gate
that fails the old engine. **Write analytics from the textbook, never
from the code.**

---

# IV. Illustrations

_Framing line for the section: each of these is a study you can re-run
with one call; the numbers on the slides are parsed from the committed
reports, not typed in._

### 10. e2e: telescope + relay

**Claim:** closed-form seed → diffraction-limited system, automatically.
**On the slide:** the telescope + Offner relay layout, the DL field
range, the seed-to-solved progression.
**Numbers, stated correctly:** D = 4 m, **f/18 system off an f/1.75
primary** (m2 = 8), 500 nm. Do NOT call it "an f/1.75 telescope" — that
is the primary's focal ratio and this room will hear the system's.
Telescope alone is DL over **±1′**; it is the Offner 1:1 relay that
takes it DL over the full **±2′** (−tilt ladder 0.015→0.043 waves,
Strehl 0.99→0.93, on pure spheres).
**Say:** the seed is *analytic* — no starting design was borrowed. That
is the design layer's actual claim: parameters in, buildable instrument
out.

### 11. e2e: segmentation

**Claim:** segments are physical apertures with per-segment channels,
not a pretty picture.
**On the slide:** the segmented pupil, physical apertures, per-segment
sensitivity channels.
**Say:** the honest gotcha, if asked — a segment's grid frame must match
its clocked monomial frame, or a per-segment poke images as a piston
instead of a localized figure. We found that as a real engine bug
(SrfType-9 was passing a null grid frame), and the *symptom* was exactly
that staircase.

### 12. e2e: MET + closed-loop simulator

**Claim:** design → metrology → control, with numbers at each hand-off.
**On the slide:** the MET layout, the closed-loop result.
**Say:** the two-pass ±u rule and the Tikhonov-regularized reconstructor
— raw pseudo-inverse diverges. Worth one sentence because it is the kind
of detail that tells this audience we have actually flown the loop.

### 13. Bench builder → the Twyman-Green pol-PSI gauge

**Claim:** you can build an *instrument*, not just an optical train —
and then discover the instrument is lying to you.
**On the slide:** the TG rig from visible `Bench` calls; fringes;
actuator poke; the 4-step PSI map beside truth.
**Say:** this is the setup for live beat (a). Two notes:
- **TG vs MZ trade** (Dave's ruling: TG for DM figure). MZ's case is
  isolation, two ports, transmission and dynamics.
- **Mechanism note, if anyone asks how the analyzer sweep is free:** it
  is the *exact bilinear analyzer basis* — three traces per arm span
  every analyzer angle. It is **not** post-trace Jones at the detector,
  which is wrong at the ~3% level because the field is ~2.8e-2
  non-transverse there. Get this right on the slide; it is the kind of
  claim this room will test.

### 14. CTB: coronagraph testbed

**Claim:** vortex EFC drives contrast to the floor, and the layers that
matter are modeled.
**On the slide:** the EFC contrast curve; band, polarization and
vector-vortex layers.
**Say:** the vector-vortex verdict is the interesting part — a scalar
model and a vector model do not agree, and we can say which regimes need
which.

### 15. e2e6m: telescope + coronagraph, closed loop under drift

**Claim:** the loop *holds* the dark zone while the system drifts.
**On the slide:** dark-zone contrast, open vs closed loop under the same
drift.
**Say:** the lesson worth stating: an optimizer will always find your
model error. The closed-loop result is only meaningful because the
model, the sensing and the control share one prescription.

---

# V. Rodgers — AI-driven design

_Section framing, say before slide 16: "These are Mike's problems. We
did not pick them, we did not tune them, and the first thing we did was
reproduce his answers."_

### 16. The challenges: the problems, the method, the ground rules

**Claim:** a fair test, run under rules stated in advance.
**On the slide:** the three challenge problems; the CC-driven method;
the ground rules (his parameters, his constraints, our metric stated
first).
**Say:** the method is the point of Section V — a human sets the
problem and the rules, the AI drives the tool, and every result is
gated. Also: **status glance at the live solve here** (Dave's ruling) —
one line, no dwelling.

### 17. Rodgers 1: the trust anchor

**Claim:** his three designs reproduced at **0.99–1.00×** his own
metric.
**On the slide:** the three designs, his number vs ours, ratio column.
**Say:** this slide is the price of admission for slide 20. Say plainly
that the residual difference is his edge-weighted sampling, not a
physics disagreement — we chased it down and can show it.

### 18. Rodgers 2: the buildability finding — an honest negative

**Claim:** we found a result we then had to retract, and the retraction
is the interesting part.
**On the slide:** the unconstrained design with the beam passing through
a mirror; the constrained design that actually builds; the price.
**Say:** **Dave caught this by reading the layout figure against the
gate** — the clearance model had three defects and the drawing showed
what the number did not. Two lessons for slide 24: *graphics are gates*,
and *a constraint's price is worth quoting* (paying real clearance costs
roughly a factor of two, and lands his own published number).

### 19. Rodgers 3: the family + solve discipline

**Claim:** the discipline is what makes the solve work, and it is
stateable.
**On the slide:** the five-stage ladder (S1 symmetric on-axis → S2 the
disaster map → S3 solve at the used field → S4 tilts/decenters + live
constraints → S5 freeform), plus the solve doctrine.
**Say:** the two rules worth naming: **solve at the field you actually
use** (S2→S3 is the whole lesson in one picture), and **pin power to the
radii and tilt to the pointing** — releasing them chokes the solve. Both
were learned the hard way.

### 20. Rodgers 3: the numbers — honest rows

**Claim:** **45.4 nm against his 53** — with every row shown, including
the ones that did not beat him.
**On the slide:** the full ladder table, ours vs his, with the metric
line at the top of the slide.
**Say:** state the metric first (strict RMS WFE, sphere on the spot
centroid, exit-pupil anchor, piston-only removal, headline = dense-map
maximum). Then the honest attribution: the gap was **solve-field count,
not iterations** — 9 fields under-determine 82 variables, so the solve
set converges while the dense map stalls. Give him the credit that his
53 was achieved without that diagnosis available.

### 21. The product: family + driver → adjacent problems

**Claim:** the deliverable is not a design, it is a driver that solves
the *next* problem — with or without CC.
**On the slide:** the walk frontier as compiled knowledge — a field-vs-
packaging frontier produced by one run, with the spec-compliant box
identified.
**Say:** the frontier is the hinge of the whole session. It is not a
table of results; it is *knowledge about the instrument* that the driver
can extend on demand. That is the claim slide 22 tests live.

### 22. LIVE DEMO

**Lead-in framing (Dave's ruling — say the contrast out loud):** the two
beats carry **different claims and different drivers.**
- **(a) is the PRODUCT in an engineer's hands** — Dave solo, desktop
  MATLAB. The script's pauses and animation gate on the desktop, and CC
  driving it would collapse the product beat into a second AI beat.
- **(b) is the designated CC-LIVE moment** — AI drives the design tool.

This is slide 21's "with or without CC" made concrete.

**Beat (a) — pol-PSI Twyman-Green DM gauge.** Visible `Bench` calls →
fringes → actuator poke → analyzer sweep with no re-trace → 4-step PSI
map beside truth.
*Open decision (Dave): keep beat 3 — the 11.7% finding, live.* The
beamsplitter's diattenuation rotates the test arm off orthogonal; the
gauge reads 11.7% high while fringe contrast moves 0.17% — a ~69×
blind spot — and the fix is solved in five traces. **TO and CC both
recommend KEEP:** it is the model-catches-what-the-instrument-cannot
beat, and it is the single most persuasive thing in the session for an
audience that builds real interferometers.

**Beat (b) — the adjacent problem, driven by CC.**
Spec taken at the top of the session; composed and launched **visibly on
the projector** (the driving is part of the show); status glance at
slide 16; reveal here.
Structure: spec relay → **frontier prediction stated before the answer
exists** → one warm-started continuation solve → reveal (map + layout +
gates).
**The claim being tested:** the walk table is compiled design knowledge,
and the driver extends it to a box nobody solved in advance — with the
answer *predicted before it is computed.*

**Backup snapshots behind BOTH beats. Rehearse twice.**
Full spoken script, scripted refusals and fallback ladder:
`templates/10_telescopes/offset_imager/demo_adjacent/REHEARSAL.md`.

### 23. Live results recap

**Claim:** what we just saw, in numbers, with the gates named.
**On the slide:** the verdict block from beat (b) and the PSI map from
beat (a) — the artifacts the runs actually produced, not a re-drawing.
**Say:** land on three things in order: the headline with its metric,
predicted vs measured, and the gates (PASS/FAIL, named). If a gate
failed, that is the *good* version of this slide — say why.
**Backups if needed:** the pre-generated bundle; then the committed walk
artifacts; then the frontier table.

---

# VI. Close

### 24. AI in design + analysis — the questions

**Claim:** the interesting question is not "can it design?" but "how do
you know it did?"
**On the slide:** four headings, with one concrete instance each.
- **Verification discipline.** Non-vacuity: every gate must fail the
  broken engine. Write analytics from the textbook, not from the code.
- **Graphics as gates.** Dave read the layout figure against the
  clearance number and found three defects the number hid. A drawing is
  a test.
- **Honest rows.** Retracted results, priced constraints, and a driver
  that refuses out-of-envelope asks in one accurate sentence rather than
  returning a plausible number.
- **Where the human rules.** Setting the problem, ruling on trades,
  signing off outward-facing claims — and noticing when a number and a
  picture disagree.
**Say:** this is the slide the room will actually argue with. Invite it.
**⚠ Dave + CC pass before this is final.**

### 25. Summary + where the public code lives

**Claim:** it is open, it is anchored, and you can run it.
**On the slide:** the repos, the branch model, how to build, where the
examples and the test suites are.
**Say:** close on the anchors, not the features — "check it against your
own CODE V models; that is exactly what our geometric suite does."

---

## Build notes for the `.pptx`

- Generator only: a `deck_keysight.py` in the pattern of
  `challenges/rodgers3/deck_rodgers3_walk.py`, with every number parsed
  from committed records via shared parsers. **Never hand-edit the
  `.pptx`** — edit the generator and rebuild.
- Figures: reuse committed PNGs (rodgers3 walk + endgame, CTB, e2e6m,
  tg_psi_dm, demo_adjacent). Do not regenerate figures for the deck;
  a figure that does not match its committed report is a defect.
- Run the `doc/STYLE_REPORTS.md` §5 pre-write gate before every build
  (it lives at the **repo root** `doc/`, not under `mmacos/`; the copy in
  `MACOS_resources/doc/` is stale — it predates the 2026-08-24
  `DECK_STYLE.md` ruling, which supersedes §3 on titles: plain
  descriptive titles, no aphorisms, kickers carry the headline number).
- Title carries **DRAFT** until Dave signs off. *(Practice, not
  doctrine — I checked, and there is no DRAFT-marking rule in
  `STYLE_REPORTS.md` or `DECK_STYLE.md`.)*

## Open items (mine, for Dave)

1. Slide 3 history — sources and sign-off.
2. Slide 24 — Dave + CC pass.
3. Slide 22 beat-3 keep/cut ruling (TO and CC both recommend KEEP).
4. Slide count vs time: 25 slides plus two live beats is a lot. If the
   session runs short, my cut order would be **7, 12, 14** (in that
   order) — 7 is process, 12 and 14 are the two illustrations whose
   claims are most nearly duplicated by 11 and 15. I would protect 9,
   17, 18 and 21 at all costs: they are the trust chain that makes 20
   and 22 mean anything.
