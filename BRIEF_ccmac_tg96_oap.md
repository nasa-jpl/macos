# BRIEF for CCMac: the 96x96 Twyman-Green DM gauge, all-reflective -- both lenses to OAPs

From CCL for Dave, 2026-09-10.  Dave: "brief CCMac on the 96x96 T-G IFO
and ask them to implement the reflective version, replacing both lenses
with OAPs."  Everything you need is on `dev-candidate` in both repos as
of today (resources `2db7e1d`, macos `d2d756a`) -- pull both first.

## What the instrument is (the record, do not re-derive)

`MACOS_resources/mmacos/templates/40_benches/tg_psi_dm96/` -- a
polarization phase-shifting Twyman-Green surface gauge for a Xinetics
96x96 DM (1.0 mm pitch, 96 mm beam), built by
`macos.design.twyman_green('polarizing',true, ...)` and driven by
`tg96.m` (Stage A clearance solve with real bodies and printed margins;
A2 sampling budget, asserted; B build; C battery; D actuator-lattice
transfer), then `tg96_tail.m` (tail retune), `tg96_eprime.m`
(actuator-space rescore), `tg96_s3/s4/s5noise/s6color.m`.  The README
is the full record; every stage has a `tg96_*_report.txt`.

The train (test arm; the reference arm shares the front end and the
tail): point source -> baffle -> **L1 (collimator, F1 = 857 mm, the
uniform 96/56 scale of the v1 rig)** -> input polarizer -> plate BS at
AOI 7 deg (the clearance solve's answer) -> compensator -> arm QWP
(double-passed) -> DM (GridData test optic, retro) -> BS -> Recomb ->
output QWP -> rotating analyzer -> **L2 (focusing, F2 = 429 mm)** ->
FocalMask -> field lens FL (pupil relay, tuned) -> detector at the DM
pupil image.  All lengths mm; conics on L1/L2 are the l2_trade values.

Numbers of record (run 10 + S3/S4/S6, model 1024, NGRID 385, 96x96):
null 0.1345 nm rms after the tail retune; single 10 nm actuator
differential 0.92 gain / 46 pm floor, base-independent (flat vs 30 nm
working state); dense random 10 nm rms read at 0.77 (instrument
roll-off); the roll-off is GEOMETRIC (identical at 480..780 nm, S6) --
the null-tuned tail's conjugate + 0.136 mm distortion.  The 1 pm target
is 46x below the single-actuator row.

## Why reflective (Dave, README "Next configurations" item 1)

Replace the lenses with OAPs.  It removes the transmitted-glass-path and
homogeneity rows from the cost budget, buys achromatic legs and no ghost
surfaces; the price is OAP alignment sensitivity -- and the battery
measures that unchanged.  The IFO's known lever is geometric (the tail
conjugate / distortion), so the reflective train is also the natural
place to re-attack the 0.77 dense-random gain.

## Rulings and constraints (standing law)

- **No engine work.**  Builder / example / runner level only.  If a
  capability is genuinely missing at the engine, STOP and write it up.
- **Both lenses, L1 and L2, become OAPs.**  The field lens FL in the tail
  is NOT in Dave's ask; keep it, but make its replacement a flag if it
  falls out naturally (see open question 1).
- **`'optics','lens'` (default) must emit BIT-IDENTICALLY** -- gate it in
  `tBench` the way `'polarizing'`, `'pbs'`, `'tail_arch'` and
  `'mask_prop'` are gated.  The lens rig is the record; the OAP rig is a
  variant beside it.
- **Keep the CONJUGATES, not the focal lengths.**  An OAP's pole-to-focus
  distance is `r = f_parent / cos^2(AOI)` (`Bench.add_oap` docstring;
  `e2e6m/s3_backend.m` paid 1.011x / 5.7 waves of pure focus for placing
  a marker at f instead of r).  Choose OAP1 with `'mode','collimate',
  'focus_dist', F1` (the source sits F1 from L1's pole today) and OAP2
  with `'mode','focus', 'focus_dist', F2`, so the collimated beam
  diameter, the mask-plane scale and the DM pupil conjugate stay what
  the lens rig has -- the sampling budget (Stage A2) and the tail
  bookkeeping then carry over.  The parent focal lengths follow from the
  fold angles.
- **Functional stops go on flat marker planes.**  `add_oap`'s `aprad` is
  metadata only: a circular ApVec on an off-axis section is applied
  about the parent VERTEX, far from the beam, and blocks everything.  Use
  `add_reference` / `add_baffle` (vertex == pole) for apertures.
- **Fold angles are free parameters; near-normal is the physics
  choice** (small AOI = small off-axis aberration; `e2e6m` used 5-6 deg)
  but the fold has to clear real bodies: redo Stage A for the folded
  layout with the same discipline -- hull bodies, field-box footprints,
  margins printed as NUMBERS (memory `feedback_clearance_gates`).  The
  BS at 7 deg and the polarizing elements (all in collimated legs, axes
  given in each leg's LOCAL transverse plane by `ax_local`) do not care
  which way the collimated beam arrives.
- **Perfect-conductor OAPs first** (the stock `IndRef=1, Extinc=1e22`
  idiom: RS=-1, RP=+1, polarization-neutral) so the lens/OAP comparison
  is geometric; a coated (protected-Al, `coat_set`) row is the STRETCH.
  Note where each OAP sits in the polarization train: OAP1 precedes the
  input polarizer (state is re-defined after it), OAP2 follows the
  analyzer (a diattenuator on a LINEAR state = common-mode amplitude,
  no fringe phase) -- so the reflective train adds no first-order
  polarization systematic.  Measure it anyway: print the arm-state
  departure from orthogonal, as v1's gate does.
- **The tail was tuned for L2's aberrations.**  After OAP2 the
  `tg96_tail.m` optimization (FL_F / FL_Kc / D_MASK_FL / DET_TRIM /
  MASK_TRIM) must be RE-RUN for the OAP rig; do not inherit the lens
  trims.  Report the null before and after.
- **The runner rule (Dave, standing, all build tasks):** the system
  ships with ONE parameterized runner a user edits and reruns without
  AI.  The ZWFS twin has it since today -- `zwfs_dm96/zwfs_params.m` +
  `zwfs_run.m` (+ `zwfs_run_figs`, `zwfs_run_batch`, `zwfs_batch.sh`;
  README "Run it yourself") -- copy the shape: `tg96_params.m` (every
  knob at its value of record, including `optics` = 'lens'|'oap' and
  the OAP fold angles) + `tg96_run.m` (stages bench / battery / figs;
  `<tag>_report.txt` + `.mat` + PNGs in `runs/<tag>/`).  Its
  equivalence gate is the lens rig reproducing the S3/S4 record numbers
  above; then the OAP rig runs through the SAME runner.  Machinery to
  reuse verbatim: `../dm_gauge_lib` (`dmg_ifo_gauge`, `dmg_frame`,
  `dmg_anchor`, `dmg_register`, `dmg_samp`, `dmg_act_fit`,
  `dmg_modal_corr`, `dmg_stencil`, `dmg_lit`, `dmg_say`).
- **Sampling / memory:** MODEL 1024, NGRID 385 (the 1 Mpix-class
  detector, tg96's rule) -- ~11 GB per MATLAB; one model size per
  MATLAB process; sequential runs.  If the box is memory-bound, the
  engine reads `macos_param.txt` from the CURRENT DIRECTORY first
  (`find_macos_file`), so a trimmed size table dropped into the run dir
  works with no engine change -- `zwfs_dm96/macos_param_2048.txt` is
  the worked example (mGridSrf 200->4 saves 1.7 GB at 1024; keep
  `mGridMat` >= the DM grid, 384 here -- the `nGridMat=` parse is
  unguarded).

## Deliverables, in priority order (stop at the timebox, never thin gates)

1. `twyman_green` option `'optics'` ('lens' default, bit-identical
   gate; 'oap' = OAP1 for L1 in `front_end`, OAP2 for L2 in `tail`),
   with `OAP1_AOI` / `OAP2_AOI` (deg) and the conjugate rule above; the
   emitted `.in`s beside the lens rig's.  Gate: lens default byte-equal;
   OAP rig's pole-to-focus distances == F1 / F2 (assert from the Bench
   objects); the collimated beam diameter at the BS == the lens rig's
   to 1e-6.
2. Stage A for the folded layout: the smallest OAP fold angles whose
   beams and bodies clear (margins printed) inside the 700 mm leg cap;
   layout PNG (`view_std`).
3. `tg96_params.m` + `tg96_run.m` (runner rule).  Lens rig through it
   = equivalence gate against `tg96_s3_report.txt` / `tg96_s4_report.txt`
   (state the tolerance you meet, not one you chose).  Then the OAP rig:
   tail retune, null, arm-state departure, the S3 modal transfer (12
   modes, 96x96 + 48x48), the S4 four rows + break ladder, all in
   actuator space, all pm -- side by side with the lens numbers.
4. Alignment sensitivity (the price): OAP1/OAP2 decenter (x, y) and tilt
   (about x, y) at 10 um / 10 urad, one at a time, on the null and on
   the single-actuator differential row; a table of pm per um and pm per
   urad.  The differential protocol should cancel most of it -- say how
   much survives.
5. STRETCH: the coated-OAP row (protected Al via `coat_set`, the
   published-stack pattern) -- polarization cost of the reflective train
   under 'polarizing' true.
6. README section in `tg_psi_dm96/` ("Reflective variant") + the
   memory-style summary at the top of your report; the runner's "Run it
   yourself" lines.

## Open questions for Dave (answer in your report; do not block on them)

1. FL (the tail's pupil-relay lens): keep, or make it a small OAP too
   (`'optics','oap'` would then be all-reflective end to end)?  Report
   the cost of keeping it (one transmitted singlet, sub-mm beam) so Dave
   can rule.
2. Fold plane: both OAPs folding in the same plane as the BS (compact,
   planar bench) vs. orthogonal planes (astigmatism partly cancels).
   Pick one by the Stage A margins and say why.

## Traps (paid for once already -- do not repay)

- matlab -batch: script files + `exit(0)` ONLY in the batch wrapper
  (`zwfs_run_batch` pattern), never in the user-facing runner; one
  model size per process; `MACOS_HOME` set.
- `add_oap` 'collimate' wants the incoming chief DIVERGING from a focus
  one conjugate back along the incoming chief -- the source-fed case.
- Registration is deck-dependent by doctrine (4 DOF classes; never
  register on a symmetric target; parity + sign from an off-center
  poke) -- the OAP rig gets its own registration, never the lens rig's.
- Read-tool PNGs go stale on overwritten paths: verify by printed
  numbers, unique filenames per run (`runs/<tag>/`).
- Report a designed NULL absolutely, never as a relative error.

## Standing rules

Work on a branch `tg96-oap` off `dev-candidate` in `MACOS_resources`;
commit as you go; push the branch when the lens-rig equivalence gate is
green (Dave orders merges).  Report:
`MACOS_resources/mmacos/templates/40_benches/tg_psi_dm96/REPORT_oap.md`
(numbers first; what was departed from in this brief and why).  If the
engine needs something, STOP and write it up -- no physics shortcuts.
