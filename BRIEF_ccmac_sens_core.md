> **SUPERSEDED 2026-09-08 by `BRIEF_ccmac_jpl_private_verification.md`** (cold-start version covering the engine changes of 2026-09-08 and the pupil-read ruling).  Kept for the protocol rationale.

# BRIEF for CCMac: sens-core verification on JPL-private prescriptions

Queued behind the EP-checks arc (start after that is done and
reviewed — Dave).  Goal: extend the sens-core acceptance evidence to
JPL-private Rx so Dave can merge `sens-core` → `dev-candidate`, then
file the merge request to `dev`.

## What sens-core is

MACOS_resources branch `sens-core` (tip `10dc0eb`): the four
`dw_d*_multi` sensitivity supervisors (~700 near-identical lines each)
rearchitected onto ONE shared core (`+macos/private/dw_multi_core.m`)
with thin per-family fronts.  Public signatures, options, outputs, and
error ids unchanged.  Accepted on this side with 10 A/B run-sets on
e5hex1 (all four families × defaults / configs+pupil_find+groups /
frozen-EP variants), every output struct BYTE-IDENTICAL to the
pre-refactor tree, plus the sensitivity suite green (87/87).  One
deliberate behavior change (Dave's ruling): a fully-vignetted field
WARNS once per run and contributes 0 rows — it no longer hard-errors
(dw_dx_multi's old behavior) nor passes silently (the other three's).

## Setup on the JPL Mac

1. Pull MACOS_resources `sens-core` (`10dc0eb`) AND note the PRE
   reference commit `eda6c5a` (the last pre-rearchitecture
   dev-candidate commit) — the harness runs a worktree there.
2. Pull the macos ENGINE at dev-candidate `6bab7af` or later (the NS
   flow-of-light fix; the sens-core suite's tNsFlowOfLight gate
   requires it).  Rebuild engine + relink the mmacos mex
   (FC=gfortran; `makems.sh release gfortran` then `make` in mmacos).
   ONE engine build serves both sides of every A/B.
3. The harness is committed on sens-core:
   `mmacos/tools/sens_core_ab/` (sens_core_ab.m + compare + README).

## Test 1 — the supervisor axis (the merge gate)

Per deck (IRIS family — start `iris_dp_v14.in` — OPTIIX, and whatever
Luis points at; prefer strongly curved NS optics, gratings,
historically-odd decks):

- `sens_core_ab(<pre_setup>, deck, out_pre/<deck>, fx, fy)` then the
  same against the sens-core tree, then `sens_core_ab_compare`.
- ONE matlab -batch per invocation (state-leak doctrine; drive from a
  shell loop).  Sequential MATLAB only.
- Same engine both sides → the engine axis cancels by construction.
- fx/fy INSIDE the deck's vignetting margin on the PRE tree (its
  dw_dx_multi still hard-errors on an empty field; nominal chief can
  be off-axis — e5hex1's is 3.5e-4 rad off in y).  Start 1e-5 rad.
- Error parity is a pass (deck refuses an option identically on both
  sides).  Model size: 256 default; raise for big-grid decks
  (mGridMat caps at 256 regardless — model 512 does NOT raise it).

**Expectation (Dave 2026-09-07): our NS decks are segment-class —
aperture-partitioned segments on a shallow common base, one crossing
per ray — so BYTE-IDENTICAL everywhere.  Any delta = STOP and report
(deck, set, field, max|diff|); it is a defect until proven otherwise,
not an adjudication case.**

## Test 2 — the engine axis (separate, cheap)

The engine NS fix (macos `6bab7af`: NS candidate probes pick the
smallest positive root > 1d-8 instead of the |L²−mpr| proximity
metric) is expected to be a NO-OP on segment-class decks.  Per NS
deck: trace to the last element on the OLD engine and the NEW engine
(two builds, or checkout+rebuild twice), compare ray positions +
ray-status counts bit-wise.  Segment decks: identical expected; a
difference is a finding to report with the die-map (launch-radius
histogram of differing rays).  The exceptional multi-crossing class
(Luneberg-like nested shells, deep NS conics) is the only place
legitimate differences live — Luneberg itself is already gated
(tNsFlowOfLight: pre-fix it skipped refractions on 20 of 100 fan
rays).

## Deliverable

One report back to Dave: per-deck table (deck | families run | sets |
supervisor verdict | engine verdict | notes), plus the compare logs.
Deck names and numbers only — no prescription content leaves the JPL
side.  On all-green: Dave merges sens-core → dev-candidate, deletes
the branch, and files the MR to dev.

## Context pointers (this repo)

`macos/BRIEF_luis_round3.md` (the whole arc + protocol rationale),
`mmacos/doc/SENSITIVITY_TOOLS.md` (the stack map),
`mmacos/tools/sens_core_ab/README.md` (harness rules).
