# BRIEF for TO: polarization phase-shifting Twyman-Green — the Keysight demo build

_Tasking: Terminal Opus.  Supervisor: CC.  Time-critical: the CodeV
demo is ~2026-09-01 — build + gates + rehearsal materials within ~2
days, then Dave's dry run.  Queue note: this may preempt
BRIEF_focal_surface Part 2 (Part 1, the engine fix, should land first
if in flight — Dave sequences).  Cold-start reads: root
`macos/CLAUDE.md` + `macos_f90/CLAUDE.md` (Phase-3 polarizing
elements + Tranche-1 sections) + `MACOS_res_dev/mmacos/CLAUDE.md`,
memory `project_bench_builder` + `project_dm_surface_models` +
`feedback_real_examples`, and the two substrate files:
`src/+macos/+design/twyman_green.m` and the
`templates/40_benches/bench_ifo/` example (its PSI processing)._

## What we're building

A **polarization phase-shifting Twyman-Green surface gauge with a DM
as the test optic**, as a shipped example plus a rehearsed live-demo
script.  Dave's topology ruling (2026-08-28): TG over Mach-Zehnder for
a DM — normal incidence, natural null, fewest reference surfaces;
polarization phase-shifting supplies the vibration-immunity story.
The demo beats: build the rig in visible Bench calls → layout view →
fringes → poke the DM → fringes show the actuator → sweep the analyzer
→ phase-shifting fringes walk with NO re-trace → 4-step PSI recovers
the DM map beside the truth map.

## Physics design (the standard pol-PSI TG)

- The splitter is a POLARIZING BS in concept: test arm carries p,
  reference arm s (orthogonal linear states).  Each arm double-passes
  a QWP at 45° so the return polarization is rotated 90° and the PBS
  routes ALL returned light to the output port (nothing back to the
  source) — state this in the README; in the model it is bookkeeping
  across the two arm decks (below).
- At the output the arms are orthogonal → no fringes until projected:
  output QWP at 45° makes them opposite CIRCULAR, then an analyzer at
  angle θ gives interference phase 2θ.  Four analyzer angles
  0/45/90/135° = phase steps 0/π/2/π/3π/2 → standard 4-step PSI.
- DM surface height h at normal incidence → OPD 2h (mind the factor 2
  in the recovery closure).

## Implementation map (all binding/example level — NO engine work)

1. **Two decks, two traces** — the engine does not split rays (a PBS
   is two traces, per the Phase-3 docs).  `twyman_green` already
   emits both arms as separate Bench objects sharing the detector
   plane; extend it with a `'pol'` option that inserts the REAL
   polarizing elements (`Bench.add_polarizer` / `add_waveplate`,
   both on this branch) into each arm: an arm-state polarizer + the
   double-passed QWP (double pass = the element appears TWICE in the
   sequential train, exactly like the rig's existing double-passed
   compensator).  `'pol', false` (default) must leave the existing
   rig BIT-IDENTICAL — gate it.
2. **Output QWP + analyzer are MATLAB-side Jones operations** on the
   per-arm vector fields at the detector: trace each arm with
   polarization + vector mode ON (`macos.polarization`,
   `macos.vector_diffraction`; `cfield_plane_get`/`complex_field(...,
   'plane',k)` REFUSES component planes unless ifVecDif3), read
   Ex/Ey(/Ez) per arm, apply QWP·analyzer(θ) Jones to each, coherent
   sum, square.  This is deliberate twice over: it respects the
   Tranche-1 seed rule, and it makes θ a FREE POST-TRACE KNOB — the
   live analyzer sweep costs zero traces.
3. **Tranche-1 rule (the trap that costs a rebuild if missed):** a
   polarizing element placed AFTER the first physical-optics leg
   transforms rays but never reaches the grid.  Each arm's deck must
   keep ALL polarizing elements before its single NFPlane→Geometric
   detector leg (the whole TG train is geometric until then, so this
   is natural — but ASSERT the emitted element order, don't assume).
   Reference pattern: `tests/Rx/Rx_PolElt.in` and its header comment.
4. **The DM test optic**: the rig's `'to_grid_file'` GridData path +
   `macos.write_grid_file`.  Build a small actuator-grid influence
   map (Gaussian influence functions, the project_dm_surface_models
   pattern); the live poke = rewrite the grid map (or grid setters) +
   **`macos.modify()`** — grid setters do NOT dirty the cached trace
   by themselves (the grid-setter-retrace memory).
5. **Detector registration**: both decks share the detector plane by
   the rig's construction — verify per-arm `dx_at` equality and pixel
   registration numerically (a gate), never by eye.
6. **PSI recovery**: reuse/extend the bench_ifo example's PSI
   processing for the 4 analyzer frames; closure = recovered surface
   vs the injected DM map (factor 2, mask edges excluded), report RMS
   residual.

## Deliverables

- `twyman_green` `'pol'` extension (+ help), example dir
  `templates/40_benches/bench_ifo_pol/` per the real-examples rules
  (self-contained driver, saves .in + .mat + figures, NO exit(0) in
  the example itself), README carrying the topology trade note
  (TG-vs-MZ for DM gauging — lift from the 2026-08-28 discussion:
  normal incidence/natural null/fewest references vs MZ's isolation/
  two ports/transmission testing; MZ's case is dynamics, not figure).
- `demo_tg_psi.m` — the LIVE script, in labeled beats with explicit
  pauses/sections: build (Bench calls visible) → `view_std` layout →
  null fringes → DM poke → analyzer sweep (4 frames + an animation
  loop) → PSI map beside truth.  Every beat seconds-fast at a modest
  model size.  **Backup PNGs of every beat pre-rendered** into the
  example dir — a live hang must cost ten seconds, not the demo.
- Gates (extend tBench or a new tTgPol, SUITE per its model size):
  (1) `'pol'` off ≡ existing rig bit-identical; (2) flat-DM null:
  analyzer sweep modulates total detector power as cos(2θ+φ₀) with
  pinned visibility; (3) 4-step PSI closure recovers the injected map
  (pin the residual); (4) emitted element order respects Tranche-1
  (polarizing elements precede the PO leg); (5) detector registration
  (dx + pixel alignment) between arms.
- A short numbers-first report section (fringe visibility, PSI
  residual, runtimes per beat) appended to the example README +
  memory update (extend `project_bench_builder`).

## Traps (paid for once already — do not repay)

- matlab -batch: script files + `exit(0)`; one model size per MATLAB
  process; MACOS_HOME=/home/dcr/dev/macos/macos_f90.
- Bench decks are BaseUnits mm — `macos.coating` takes BaseUnits (the
  tJonesPupil units trap); waveplate `Retardance=` is waves at
  parse-time λ, stored physically (a query after `set_src_wvl`
  legitimately changes).
- Polarizer axis is taken modulo its component along the element
  normal; axis ∥ normal extinguishes (prescription error).
- Masks in a deck must be `Element= Obscuring` to bite the grid.
- Model size ≥ the deck's nGridpts, or the engine silently resets and
  that path is NOT safe.
- Read PNGs are cache-stale on overwritten paths — verify figures by
  the numbers they print, use unique filenames per beat.
- Obscured rays carry RayE=0 → gate any per-ray probe on LRayPass
  before interpreting phases.

## Standing rules

Commits LOCAL (resources only; the engine is untouched by this task);
Dave orders pushes.  If a needed capability is genuinely missing at
the engine level, STOP and write it up for Dave's ruling — do not
work around with physics shortcuts.  The deck slides that present
this (demo deck slide 13 + 22) are assembled SEPARATELY under the
STYLE_REPORTS §5 gate — this task ships the example + figures + live
script, not the .pptx.

---

# DELIVERED 2026-08-28 — resolutions, and where the brief was departed from

Everything is LOCAL on `dev-candidate` in `MACOS_res_dev` (resources
only; the engine is untouched).  Full technical write-up:
`mmacos/templates/90_polarization/tg_psi_dm/README.md`.

## Shipped

| | |
|---|---|
| `templates/90_polarization/tg_psi_dm/example_tg_psi_dm.m` | the gated measurement (5 gates + closure) |
| `…/demo_tg_psi.m` | 7 narrated beats, a PNG each, live analyzer animation |
| `…/dm_influence_map.m` | actuator influence-function DM surfaces |
| `…/README.md` | topology trade + numbers-first report |
| `tests/tTgPol.m` | 9 gates, model 128, added to `SUITE_FAST` |

## Departures from the brief, and why

1. **Directory.**  The brief says `templates/40_benches/bench_ifo_pol/`.
   That name is already taken — the pol bench lives at
   `templates/90_polarization/bench_ifo_pol/` (slices 1–3, commit
   `1963408`).  A second directory of the same name in a different
   thread would be actively confusing, so this went to
   `templates/90_polarization/tg_psi_dm/` and `00_INDEX.md` was
   updated.

2. **Implementation-map item 1 was already done.**  `twyman_green`
   gained the `'polarizing'` option in slice 3 (`1963408`): real
   `TrPolarizer`/`WavePlate` elements, arm QWPs double-passed, default
   false emitting bit-identically.  **No builder change was needed**,
   and the bit-identity gate the brief asks for already exists as
   `tBench/test_twyman_green_polarizing` — not duplicated in `tTgPol`.

3. **Item 2 (output QWP + analyzer as MATLAB-side Jones) was NOT
   done, deliberately — it is wrong, and there is an exact alternative
   that delivers the same free knob.**  Applying the analyzer's Jones
   to the *detector* field is not the same operation as an analyzer
   *before* the tail: the tail transports each ray's vector first, so
   `P·T·E ≠ T·P·E`.  Measured at the detector, the field is 2.8e-2
   non-transverse (`|E_long|/|E_1|`), i.e. the shortcut is wrong at the
   several-percent level.  Instead the analyzer stays a real engine
   element and the sweep is made free EXACTLY: the detector field is
   bilinear in the analyzer axis, so
   `E(t) = c²·E(0) + c·s·[2E(45)−E(0)−E(90)] + s²·E(90)`.
   **Three traces per arm span every analyzer angle**, verified against
   direct traces at off-basis angles to 3.5e-10 relative, with the one
   term that breaks it measured and gated: `0.149·β²`, β = ray angle at
   the analyzer (a clean quadratic law over four decades on a Kr
   ladder).  Six traces now give the four PSI frames, a 64-angle
   least-squares fit and a 36-frame animation.

4. **Item 3 (Tranche-1) is moot on this rig but asserted anyway** — the
   whole train is `PropType= Geometric`, so no polarizing element can
   fall behind a physical-optics leg.  Gate 1 asserts the order from the
   emitted deck and separately checks that the grid responds to the
   analyzer.

## Resolutions found on the way (do not re-derive)

- **THE HEADLINE FINDING: the beamsplitter misaligns the rig, and the
  gauge reads 11.7% high.**  The design's arm-QWP azimuths (0 and 45)
  are supposed to leave the arms orthogonal-linear at ∓45.  Measured:
  ref exactly +45.0000, **test −37.5214** — 7.479° from orthogonal.
  Cause: every non-normal element between polarizer and recomb is a
  diattenuator (`t_s≠t_p`) and rotates a linear state; the ref arm's
  half-wave at 45° maps `θ→90−θ` and cancels the diattenuation before
  it against the diattenuation after it, while the test arm's at 0°
  maps `θ→−θ` and ADDS them.  A waveplate is unitary, so nothing
  downstream can repair non-orthogonality — the fix is the arm
  waveplate (+3.768°, solved in 5 traces).  **Every field phase in the
  model is exact throughout** (scalar, geometric OPD, polarized field
  at the detector, ray field at the analyzer all recover a 20 nm piston
  at gain 1.00000); the error is in the MEASUREMENT, which is why a
  model is the right place to find it.
- **Fringe visibility is not an alignment metric.**  That misalignment
  costs 11.7% of scale and 0.17% of contrast — a factor of ~69.  Both
  configurations are run and both numbers printed, so the alignment
  step cannot be silently deleted.
- **A `GridData` value IS the surface height** (the standing 1× vs 2×
  question in `project_bench_builder`): a uniform 20 nm grid piston
  recovers `4π·dz/λ` exactly, and agrees with a rigid 20 nm translation
  of the whole optic to 2.3e-10 rad.  The double pass supplies the 2.
- **`I(θ)` can contain only DC, 2θ and 4θ** — and the 6θ bin is a gate
  ONLY on directly traced frames (3.0e-14 measured over 12 traced
  angles): a frame synthesized from the quadratic basis is degree-2 in
  2θ by construction, so its 6θ bin is zero by algebra (5.1e-17) and
  would pass against any engine.  Both the example and `tTgPol` trace
  it.  The 4θ term (8.9e-4 of the fringe)
  is a real four-step systematic AT THE DETECTOR — zero at the analyzer,
  which is why slice 3's ray-level gate reads 1e-16 — and it is common
  to the poked and baseline runs, so the differential cancels it
  (1.7e-14 nm).
- **A circular state is analyzer-invariant in power**, so a single-arm
  "does the grid see the analyzer?" tripwire passes vacuously on the
  ALIGNED rig.  It must be run on the unaligned (linear) arm.
- **Diffraction-array parity is deck-dependent by documented policy**
  (`doc/opd_conventions.md` §2), so the gate calibrates the row/col↔x/y
  convention on one actuator and VERIFIES on a second rather than
  hard-coding one of the eight candidates.
- `set_elt_grid` already invalidates the cached trace — the brief's
  reminder to call `macos.modify()` after a grid poke is no longer
  needed.

## Numbers

visibility 0.998304 · analyzer basis exact to 3.5e-10 · piston gain
1.00000 (2e-5 nm error on 20 nm) · DM recovery corr 0.998434, residual
0.304 nm rms interior (0.363 whole pupil) on 6.35 nm rms of surface ·
9 traces ≈ 1.8 s for the whole analyzer basis · `tTgPol` 9/9 in 4.4 s.

## Not done

The `.pptx` slides (13 + 22), per the standing rule — this ships the
example, the figures and the live script.
