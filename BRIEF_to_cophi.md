# BRIEF — COPHI campaign, the bench stages (TO), after the five-mode item

CC for TO, 2026-10-11, on Dave's word ("then COPHI").  The plan is
`macos/BRIEF_cophi_campaign.md` (read §1–§3 and the rulings in §7); the template
skeleton exists: `MACOS_resources/mmacos/templates/40_benches/cophi_dm96/`
(`cophi_params.m`, README) and the reading's algebra in `dm_gauge_lib/`
(`dmg_cophi_frames` / `dmg_cophi_recon` / `dmg_cophi_camera`, gate `tests/tCophi.m`,
8/8).  Read the README's "what the gate already settled" first — two of the
campaign's premises moved there (more buckets do NOT cure a step error; knowing
the steps does).  This brief is the bench half: the LO arm and the stages that
run on `tg96_run`'s bones.  Resources only; the engine is not touched.  The
consolidation's rules still apply (§0 of PLAN_CONSOLIDATION: gates red first,
records, one item in flight, nothing pushed without Dave).

## 1. The LO arm (ruling 1: traced, a drop-in at the reference flat's seat)

`tg_psi_dm96_oap`'s bench has the reference flat at the splitter's second port;
COPHI replaces it by a fiber-fed collimated LO entering the recombiner.  Build it
with the Bench builder as the P/SRI did its reference (`macos.design.psri_bench`,
the two-deck reference arm is the precedent): a second `ColSource`-class deck whose
source sits at the flat's seat, aimed into the splitter, with the SAME tail and
camera as the test arm, so `complex_field(iDET)` of the two decks are the two
fields `dmg_cophi_frames` wants.  Measure, before anything else, that the two
arms' fields at the camera register (the `dmg_frame` affine of each; they must
agree to a pixel — the LO fills the camera flat, the test arm images the DM).
Knobs: `P.cophi.lo_ratio` (LO power at the camera), `P.cophi.lo_drift` (tilt /
defocus of the LO's collimator per loop cycle: the new drift row).

Stage `bench`: layout + clearance with the LO arm placed (the three-leg render
with the LO leg a fourth color), both rigs (lens, OAP mirror).

## 2. Stages on `tg96_run`'s bones

Reading class 'H' in the runner = the two fields → `dmg_cophi_frames(Et, Elo,
P.cophi)` → `dmg_cophi_camera(Fr, cam)` → `dmg_cophi_recon` → `h = phi λ/4π`
about a frozen null (the flat DM), then the record's matrix estimator exactly as
the IFO reading — same `dmg_act_fit`, same `calib_mode 'matrix'`, same λ_m.

- `battery`: the readings on the flat and the 30 nm surface, the matrix on both;
  the accuracy gain and floor rows (the record's I column beside H for the same
  seeds).  EXPECT: gain 1.000 ± the matrix bias (no polarization gain by
  construction); floor at the photon level of the camera's well.
- `noise`: photons per pm at every camera in `P.cophi.cameras` — the frame mean
  at `fill × well`, shot + read per frame, N frames per reading; the formula
  `tCophi` pinned gives the per-pixel σ_φ, the matrix fit gives the per-actuator
  pm.  The table: camera, e⁻/frame, pm per reading at its cadence, seconds per
  pm.  Feng's baseline and the Photonfocus row are the deck's two.
- `loop`: hold at 3 pm with `dmg_loop` at the camera's cadence; drifts walk /
  thermal / cam (the record's) + 'lo' (the LO tilt/defocus per cycle: measure
  what a 1 µrad LO tilt per cycle does to the actuator estimate — it is a
  low-order term the matrix fit partly absorbs; say how much).  Intra-scan fiber
  walk and jitter enter the reading at `P.cophi.walk` / `.jitter` — the loop at
  walk = 0.05 rad/frame is the row that prices fiber stabilization.
- `capture`: the range; temporal unwrapping between consecutive readings as the
  free extra (a reading every 1/7.5 s cannot alias a drift slower than λ/4 per
  reading).
- `format` — RESCOPED by Dave's ruling 6 (2026-10-11): no coarse / low-order forms
  (EFC probes are commanded at full actuator resolution).  The row is a
  FULL-LATTICE, well-free reading at ≥ 1 element per actuator — a lock-in camera
  / ToF-demodulating sensor / 128 × 128 coherent array class, ASSUMED to exist at
  1:1 or better with the pitch — scored on cadence and photons per element per
  second: `P.format.lockin` (sites, sample rate, flux per element, phasemeter
  noise), the same reading algebra at a kHz beat, the pm-per-second table beside
  the cameras'.  Camera formats below 1:1 are not run.
- `figs`: the station figure in the deck's format (command, four frames, raw map
  vs engine, reconstruction, residual), mirror rig.

Each stage a commit with its record under `runs/<tag>/`; wall time in the log;
model 1024, one MATLAB at a time (the batch wrapper), never while another
lane's MATLAB runs on a relink.

## 3. Gates

`tCophi` (algebra) exists.  Add `tCophiBench` (SUITE_FAST if under a minute,
else its own class): the two arms register to a pixel; a flat DM reads 0 ± the
photon floor; a single-actuator poke of 1 nm reads 1.000 ± 0.01 nm through the
matrix (the differential protocol); the LO-tilt row moves the estimate by the
measured amount and the matrix fit's absorption of it is stated.  Red legs: the
reading with the LO field replaced by the test field (no interference: the
reconstruction must fail by name — amp/dc below the contrast floor), and the
camera at `well = 0`.

## 4. What is CC's

The design review of anything you change in the three `dmg_cophi_*` routines
(say before you change the model); `deck_cophi` as the records land; the
detector numbers checked against the vendor sheets before a slide; the
composite deck after.  Report per stage as in the consolidation: SHA, gate
counts, the record's path and its headline numbers.
