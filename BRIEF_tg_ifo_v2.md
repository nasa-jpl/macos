# BRIEF for TO: TG interferometer v2 — real splitter physics at real AOI

_Tasking: Terminal Opus (CC's pick: context continuity on a short
clock — TO built tg_psi_dm and owns its traps; the new physics is
fenced by existing engine gates).  Supervisor: CC.  **TIME-BOXED
SPRINT (Dave 2026-08-29): kicked off NOW, demo is ~2026-09-01.**
v1 (tg_psi_dm as delivered) is the DEMO DEFAULT and stays frozen —
v2 joins the demo ONLY if delivered, gated and rehearsed in time;
missing the demo is an acceptable outcome, breaking v1 is not.
Work the deliverables in the priority order below and STOP at the
timebox rather than thinning gates.  Origin: Dave, reviewing deck
slide 12 — "combining TO's work with the earlier pol-ifo work, using
more appropriate AOI at the surfaces."
Cold-start reads: `macos_f90/CLAUDE.md` polarization sections (Phase
2 coated-branch fixes, Phase 3 elements, the off-normal axis
convention, Tranche 1); `templates/90_polarization/tg_psi_dm/README.md`
+ the delivery log at the foot of `BRIEF_tg_demo.md`;
`templates/90_polarization/bench_ifo_pol/` (slices 1–3);
memory `project_tg_psi_dm`, `project_polarization_plan`._

## What v1 idealizes, and what v2 makes physical

The delivered gauge (v1) carries the polarization CONCEPTUALLY: the
PBS is bookkeeping across two arm decks — each arm opens with an
IDEAL `TrPolarizer` at normal incidence — and the splitter/compensator
plates enter only as uncoated Fresnel surfaces.  That was the right
call for the demo, and it already surfaced real physics (the 11.7%
diattenuation rotation).  v2 replaces the concept with the COMPONENT:

1. **A real polarizing beamsplitter at 45°, as a coated engine
   surface.**  The engine has carried trustworthy coated s/p physics
   at arbitrary AOI since 2026-07-27/28 (the r_p-sign fix, the
   incident-medium fix, the transmission radiometric factor — all
   gated, plus the published-data Mueller anchor at ~1e-14).  So the
   PBS needs NO new engine capability and NO ray splitting: per the
   standing two-decks doctrine, the TEST arm's deck carries the
   splitter in TRANSMISSION (coated `Refractor` at 45°) and the
   REFERENCE arm's in REFLECTION (coated `Reflector` at 45°), and on
   the return pass — after each arm's double-passed QWP has rotated
   its state 90° — the roles swap, which is exactly the pol-TG "all
   light to the output port" routing, one sequential deck per arm.
   `RfPolarizerElt` remains a stub and is NOT needed — a PBS is
   coating physics, not an ideal element.
2. **Real AOI everywhere:** splitter and compensator at 45° (or the
   cube geometry — open question 1), waveplates near-normal with a
   deliberate small tilt if desired (the projected-material-axis
   convention already handles tilted-element geometry exactly — the
   Korger-anchored rule), detector normal.
3. **The calibration story, re-derived.**  v1's 11.7% arm rotation
   came from ideal diattenuators + uncoated plates.  With a declared
   PBS stack the per-pass diattenuation/retardance are
   coating-dependent; re-measure the arm-state rotation, re-solve the
   compensating waveplate clock, and report v1-vs-v2 side by side.
   Both configurations keep printing (the cannot-silently-delete
   rule).
4. **Merge with the earlier pol-ifo machinery** (`bench_ifo_pol`
   slices 1–3): one PSI-processing source, one Jones-analysis path;
   port whichever pieces are better and delete the duplicate.
   (STRETCH under the timebox — see deliverable 6.)

## Constraints (standing law — do not relitigate)

- **No engine work, period** (ruling 3).  Everything is builder/
  example/coating-prescription level.  If a genuinely missing engine
  capability surfaces, STOP and write it up — it becomes a separate
  brief, never a side effect of this one.
- **Two decks, two traces** — the engine does not split rays.
- **Tranche 1:** every polarizing element stays ahead of each arm's
  single NFPlane→Geometric detector leg.  A 45° splitter mid-train
  does not violate this (the train is geometric until the detector
  leg), but ASSERT the emitted order as v1's gate does.
- **Coating inputs:** `coat_set` takes PHYSICAL thickness in
  BaseUnits (bench decks are mm — the tJonesPupil units trap);
  Rx `Coating=` thickness is waves at parse-time λ.  Incident-medium
  and per-interface conventions per the Phase-2 sections — and the
  radiometric factor is applied ONCE after the Airy recursion, never
  per interface.
- **Analytic gates are written from the textbook** (Macleod/Born &
  Wolf), never transcribed from the engine's own expressions — the
  r_p-sign lesson.
- v1's gates (`tTgPol` 9) stay green untouched; v2 is a variant
  (builder option or sibling template), not a replacement, until
  Dave rules otherwise.

## Rulings (Dave 2026-08-29 — settled, do not re-open)

1. **CUBE PBS** — a single cemented coated interface at 45°, no beam
   offset; substrate faces optional.
2. **Published stack** — a MacNeille-type design from the literature
   (declared with its source, indices and thicknesses in the README).
   Drive the model with the PUBLICATION'S OWN numbers, the
   protected-Al anchor pattern; the full curve-on-curve external
   anchor gate is the STRETCH deliverable, not the core.
3. **Keep it simple, done QUICK** — the ideal retarder stands; NO
   birefringent-plate engine work, no engine work of any kind.
   Axis geometry at tilt is already exact (the projection
   convention).
4. **v2 lives BESIDE v1** — own sibling template dir
   (`tg_psi_dm_v2/` or the builder-variant equivalent), linked
   separately for demo purposes: it must be RUNNABLE AS A DEMO —
   ship a `demo_tg_psi_v2.m` in the v1 seven-beat pattern (seconds-
   fast beats, a backup PNG per beat) so Dave can choose either rig
   on demo day.  v1 remains the rehearsed default.

## Deliverables, in priority order (stop at the timebox, never thin gates)

1. `twyman_green` v2 path (`'pbs','cube'` or equivalent): per-arm
   decks with the coated 45° interface (test = transmission leg,
   reference = reflection leg, roles swapped on return) + the
   declared published stack via `coat_set`; emitted `.in`s beside
   v1's.
2. The v2 measurement: alignment re-solve under real coating
   physics, 4-step PSI closure, v1-vs-v2 comparison table (arm
   rotation, gauge scale error, visibility, PSI residual,
   output-port efficiency).
3. Core gates (`tTgPol2` or `tTgPol` extension):
   - ideal-coating / normal-incidence LIMIT recovers v1 behavior
     (pick tilt/azimuth for the real assertions where the effect
     does NOT vanish identically — the degenerate-check lesson);
   - output-port energy accounting against the declared stack's R/T
     (the radiometric fix makes |T|² physical — use it);
   - PSI closure pinned at the v1 residual class;
   - non-vacuity: at least one gate fails against a v1-style
     ideal-polarizer rig.
4. `demo_tg_psi_v2.m` — the seven-beat demo cut + backup PNGs
   (ruling 4).  Measure the beat timings; if any beat is not
   seconds-fast, say so rather than shipping a slow demo.
5. Numbers-first README + memory update extending
   `project_tg_psi_dm`.
6. STRETCH (only if the clock allows): the external anchor gate
   against the published stack's own curves; the bench_ifo_pol
   machinery merge (single PSI/Jones source).  Neither blocks
   delivery — record what was deferred.
7. AFTER landing: deck slides 12–14 refresh + the separate demo link
   (CC, under the §5 gate).

## Traps (paid for once already)

- matlab -batch: script files + exit(0) in runners only; one model
  size per process; MACOS_HOME=/home/dcr/dev/macos/macos_f90.
- Obscured rays carry RayE=0 — gate per-ray probes on LRayPass.
- A circular state is analyzer-invariant in power — run
  grid-sees-the-analyzer tripwires on the unaligned (linear) arm.
- I(θ) contains only DC/2θ/4θ; the 6θ bin is a gate ONLY on directly
  traced frames (synthesized frames pass it by algebra).
- Read-tool PNGs cache-stale on overwritten paths — verify by printed
  numbers, unique filenames.

## Standing rules

Commits LOCAL on res dev-candidate; Dave orders pushes.  If a needed
capability is genuinely missing at the engine level, STOP and write
it up — no physics shortcuts.

---

# DELIVERED 2026-08-29 — resolutions, and where the brief was departed from

All LOCAL on `dev-candidate` in `MACOS_res_dev` (resources only; **the
engine is untouched**).  Full technical write-up:
`mmacos/templates/90_polarization/tg_psi_dm_v2/README.md`.

## Shipped

| | |
|---|---|
| `src/+macos/+design/pbs_macneille.m` | the MacNeille cube design: prism index from the design condition, symmetric stack, quarter-wave-AT-ANGLE thicknesses |
| `src/+macos/+design/thinfilm_rt.m` | the textbook reference — Macleod characteristic matrix, GENERAL incident medium, r **and** t |
| `src/+macos/+design/Bench.m` | `pbs_cube` token + `add_pbs_pass`; `Coating=` emission (Model A) |
| `src/+macos/+design/twyman_green.m` | `'pbs','cube'` path; `'plate'` still emits BIT-IDENTICALLY (tBench green) |
| `templates/90_polarization/tg_psi_dm_v2/example_tg_psi_dm_v2.m` | the gated measurement (6 gates + closure) |
| `…/demo_tg_psi_v2.m` | 7 self-timing beats, a PNG each |
| `…/README.md` | numbers-first report |
| `tests/tTgPol2.m` | 9 gates, model 128, added to `SUITE_FAST` (~6 s) |

`tTgPol` (9) and `tBench` (7) re-run green, 16/16 — v1 is untouched and
stays the demo default.

## Headline numbers

| | v1 plate | v2 cube |
|---|---|---|
| arm rotation from orthogonal | 7.4786° | **5.3e-06°** |
| PSI gain as designed | 1.11661 | **0.999999** |
| alignment step | +3.768° solve | **none** |
| compensator | required | **none** (arms balance 2.3e-13 mm) |
| delivered power | 0.1691 | **0.3836 (2.27×)** |
| sep. after a 10° waveplate error | 13.15° | 0.369° / **3.1e-05°** |
| DM closure | corr 0.998434, 0.304 nm | corr **0.999573, 0.183 nm** |

Coating, engine vs Macleod: `R_s` 9.6e-10, `T_p` 4.0e-12 relative;
`R_p` null at 4e-12; `R+T = 1.000000000000` **across the two arms' decks**;
extinction `T_p/T_s` = 2382:1.

## Findings the build produced (not asked for, and worth keeping)

1. **Brewster at the H/L interfaces is NOT enough.**  The condition
   equalizes the tilted p admittances (both 2.7101), so for p the stack is
   ONE HOMOGENEOUS SLAB — and its two boundaries with the PRISM are not
   Brewster.  `H(LH)^4` (9 quarter waves, ODD → quarter-wave layer) gives
   **R_p = 2.11e-02**; `(½H L ½H)^4` (8, EVEN → half-wave absentee) gives
   **R_p = 0**.  Both are textbook MacNeille designs; one is a polarizer.
   The odd form is now the non-vacuity counterexample for every p-null gate.
2. **The v1 arm rotation is STRUCTURALLY absent, not merely smaller.**  Each
   arm's state sits ON a coating eigenaxis (test on p, reference on s),
   where a diattenuator cannot rotate it.  v1 put the state at 45° to the
   splitter axes, where the rotation is FIRST order.
3. **The error budget INVERTS — the real instrument result.**  A waveplate
   azimuth error is *cleaned* on the return pass: to `r_p` (= 0) on the
   REFLECTING arm, to `t_s/t_p` = 0.0205 on the TRANSMITTING one.  Scale
   does not move at all; the error emerges as CONTRAST.  v1's lesson was the
   opposite (11.7% of scale for 0.17% of contrast).  **A PBS converts an
   invisible systematic into a visible one** — a better reason to buy one
   than throughput.
4. **A badly specified cube still reads only 0.21% high** (catalogue SF2
   prism instead of the design index), against v1's 11.7%.

## ENGINE FINDING — written up, NOT fixed (ruling 3 honoured)

**`Coating=` is unguarded against `mCoat = 10`.**  `msmacosio.inc:2753`
reads `EltCoat(iElt)` and loops `Do k=1,EltCoat(iElt)` writing
`IndRefArr(k,iElt)` / `ExtincArr(k,iElt)` / `EltCoatThk(k,iElt)`, then
`IndRefArr(i+1,iElt)`, with **no bound check** — while `IndRefArr` is
`(0:mCoat,mElt)` and `EltCoatThk` is `(mCoat,mElt)` with `mCoat = 10`
(`elt_mod.F:330`).  An 11-layer stack **loads without complaint** and writes
past the end; the only visible symptom is `coat_get` failing afterwards (it
DOES guard, `nLayer > mCoat` → FAIL).  Found the hard way — the first
MacNeille stack here was 11 layers.
*Suggested fix (separate brief):* reject `EltCoat > mCoat` in the parser
with a named message, the way the validator handles other bad values.
`pbs_macneille` asserts on it meanwhile, so nothing in this delivery can
trip it.

## Departures from the brief

1. **Ruling 2, "the publication's own numbers", is honoured for the DESIGN
   RULE and the MATERIAL INDICES, not for a published R/T curve.**  The
   condition is MacNeille (US 2,403,731, 1946) in Macleod's form, driven
   with the classic visible pair (ZnS 2.35 / cryolite 1.35), which fixes the
   prism at n = 1.655468 — a dense flint, i.e. the model reproduces *why*
   MacNeille cubes are made of flint.  The full curve-on-curve external
   anchor is the brief's own STRETCH item and was **not** attempted; what
   replaces it is the Macleod characteristic-matrix analytic, written from
   the textbook, agreeing at 1e-10.
2. **`nperiod` is capped at 4 (9 layers) by the engine ceiling above**, not
   by design taste.  At 2382:1 extinction it is if anything more
   representative of a real cube than a longer stack would be.
3. **The `v1`-limit gate in deliverable 3 was replaced.**  "Ideal-coating /
   normal-incidence limit recovers v1 behaviour" is not a true statement:
   an ideal cube recovers a *better* rig than v1, because v1's error comes
   from its 45° conductor plate and tilted transits, which the cube does not
   have.  The structural check that actually earns its place is a **bare
   cemented interface** — glass against the same glass, `R = 0` and `T = 1`
   exactly — which proves the COATING carries the split.  Note the
   consequence is stronger than it first reads: each arm reflects off the
   diagonal once, so with `R = 0` **both arms go dark**, and the gate scores
   delivered power (bare/coated < 1e-9) rather than an azimuth, because
   `arm_state` divides by a field that is zero.
4. **Deliverable 6 (both stretch items) deferred, deliberately**: no
   external curve anchor, and no merge of v1 / v2 / `bench_ifo_pol` onto a
   single PSI/Jones source.  The v2 example COPIES v1's helper block rather
   than sharing it — v1 is frozen as the rehearsed demo default, and a
   shared helper would make every v2 edit a v1 risk.
5. **Demo timing reported honestly.**  Beats 1 and 3–6 are seconds-fast
   (0.20–3.75 s); beat 2 is 8.06 s and beat 7 is 18.13 s, both dominated by
   figure rendering/export on a box with no graphics acceleration, both the
   same code paths v1 uses.  Whole-demo totals on the SAME box: v1 49.8 s,
   **v2 47.0 s** — v2 is marginally faster, because the coating beat it adds
   costs about what the alignment beat it removes used to.  Every beat
   writes its PNG, so a stall costs ten seconds.

## Traps paid for in this build

- **An empty matrix cannot mean both "default" and "none".**  `pbs_coat`
  defaults to `NaN(1,3)`; `zeros(0,3)` is an explicitly bare interface.  The
  first version silently ran the coated stack for the "bare" gate.
- **Report a designed NULL absolutely.**  `|engine − analytic| / analytic`
  on `R_p = 2.6e-30` is noise dressed as a number.
- **The cube subtracts its own half-side from `D_RC_L2`** so the
  Recomb→L2→mask→detector conjugate matches the plate rig and the
  `l2_trade` tail trims transfer verbatim.  Without that the closure would
  have silently degraded and looked like a v2 defect.

## Not done

Deck slides 12–14 and the separate demo link (deliverable 7, CC's, under the
STYLE_REPORTS §5 gate).  No pushes — standing rule, Dave orders those.
