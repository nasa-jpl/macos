# BRIEF — PLAN_CONSOLIDATION phase 2, second list (TO) — the §1.2 items, then the bindings program

CC for TO, 2026-10-12, on Dave's word ("task them with their list").  You start from a
CLEARED session: this file is the whole hand-off.  Read, in this order: `macos/CLAUDE.md`
(root rules), `macos/PLAN_CONSOLIDATION.md` (§0 = the ten rules = the definition of done;
§1.2 and §1.4a/b = this list; §4a = Dave's rulings), `macos_f90/CLAUDE.md` (every engine
item below has its entry — the OPD-reference section, the per-ray status section, the
CALIB sections, the §1.4a paragraph), then `BRIEF_to_consolidation_p2.md` (what you just
finished, for the method).  Both repos are PUSHED and clean: macos `dev-candidate`
`e8c80ea`, MACOS_resources `dev-candidate` `2a07af9`, origin == local.  The modified
files in the resources tree (dyson5 PNGs, zwfs_dm96 PNGs/decks, an `.nohup`) are run
artifacts of earlier lanes — leave them.

Where we are: phase 1 exited 10-10; your phase-2 engine list (`BRIEF_to_consolidation_p2`:
1.1b A, B.1, 3.2, 3.3, 3.4, the five session modes, NSRefractor's grid frame) is DONE and
pushed; §3.3 is complete.  Dave's standing rulings: feature freeze PREFERRED (not
required) until the dev-candidate PR is approved; the `Opt*`/trace-state SAVE policy is
DEFERRED (do not touch the writer for it); OPL's two refusals (f3d7cea) are on Dave's desk
— leave them as they are until he rules.

## How we work (unchanged; the rules phase 1 and 2 paid for)

ONE item in flight.  Its gate FIRST, red on the pre-fix binary (`build_release_prefix*/`
side builds; CC gates from the clean worktree `~/dev/macos_cli_base`, never from the
shared tree).  Both compilers rebuilt.  The mex relinked only when no other MATLAB runs
(`ps -C MATLAB` first; say when you relink).  Fast suite after.  Corpus A/B before the
commit when the item touches surfsub / elemsub / tracesub / propsub / msmacosio / iosub /
macos_cmd_loop / zern_wf (`cli_tests/run_cli_tests.sh load` per commit; the full corpus
for anything that moves a trace), every moved deck named in a committed record
(`cli_tests/records/ab_*_moved_<sha>.txt`; `run_cli_tests.sh ab <pre> <post> load`).
Commits LOCAL, SHA + branch stated to CC; CC gates each; Dave says push.  **Measure before
you guard** (grep the value's later assignments first); **a guard warns and counts, a new
refusal needs Dave**; **user options at allocation, never in `*_init_vars`**; **a pymacos
red leg needs a pre-fix `.so` built in a side tree** (say so when you do it by hand).
Never `W ~= 0` for ray validity (`opdValid` / `opd_mask_get`).  One model-1024 MATLAB on
the box.  Cite SHAs and numbers only from output you have read.

## A. Engine items (plan §1.2 + the §1.1b leftovers), in this order

1. **`zern_fit` (zern_wf.F) finds the pupil as `OPDin /= 0`.**  Under the chief-ray
   reference (default since 10-10) a valid ray whose path equals the chief's reads
   EXACTLY 0 and is dropped from the fit — ZCO, CALIB's WFE_ZMODE target, every consumer
   of `zern_modes_coef`.  The bindings closed this with `opd_mask_get` (68ed8e2); the
   engine's own fit did not.  Fix shape: fit over `traceutil_mod%opdValid` (set at exactly
   the statements that write `OPDMat`; follows `iObsOpt`) instead of the value test.
   First MEASURE how many corpus decks have an exact-zero valid ray in their `opd nElt`
   map (the tDesignTelescope TMA had 4 of 276) — that count is the blast radius and goes
   in the commit.  Gate: a symmetric deck with exact-zero valid rays, coefficients with vs
   without those rays (pre-fix differs, post-fix the full set); plus `tCalibTargets`'s
   ZMODE leg unchanged where no zero ray exists.  `zern_fit` has other callers (ZRM, SZCO,
   the MAP path?) — list them in the commit; all take the same mask.
2. **`RayFailElt` stamps `nElt+1` for obscuration** (the end-of-trace stamp in CTRACE,
   tracesub.F ~4364, and propsub's twin) instead of the element that clipped the ray.
   The clipping element IS known where `ChkRayTrans` says no (`LRayTrans` false at that
   element) — carry it.  `SetRayFail` records the FIRST failure only, so the per-element
   stamp must be made at the clip, not at the end; the end-of-trace stamp then becomes the
   fallback for a ray clipped by a path that did not stamp.  Gate: `tGratingAperture`
   (already asserts `RayStat_Obscured` on clipped rays — add the element), `ray_status_get`
   on `Rx_Cass_FarField` (the central obscuration's element), the CLI `WARN` line per
   element.  Bindings: mmacos `get_ray_status`, pymacos `ray_status` read it unchanged.
3. **A ray that STARTS on a surface can be called a miss** (Grating_example_001's 16 / 8
   edge rays per compiler; the Return/Return/FocalPlane idiom puts vertices on
   predecessors' planes).  CC's ruling stands: **MEASURE, no surfsub change without Dave.**
   Engine-side, at each forward-branch miss, test `|L| <` a round-off band (say 16 ulp of
   the element scale, the `DeadBandUlp` convention) AND the vertex's coincidence with the
   predecessor's plane; count per deck, print ONCE per run; run the corpus on both
   compilers in a side build (the instrumented binary is NOT committed).  The count table
   (deck, misses within the band, of which on a coincident plane) goes to Dave with the
   proposal: nonzero on a legacy deck is the case for accepting `|L| < band` as a hit at
   `L = 0`; zero is the case for leaving it.  Record `cli_tests/records/start_on_surface_<sha>.md`.
4. **Two silent drops in the Opt* parse (rule 5, warn):** (a) the header chain matches
   `OptFOVWt` on 5 characters, so any `OptFo*` key is CAPTURED as FOV weights — a prefix
   the `HdrUnkNote` warning cannot see; audit every `LCMP(…, n)` in the header chain for
   a shorter prefix that swallows a longer sibling (a table in the commit), lengthen the
   matches that need it (the user-visible contract is the FULL key name; a shortened match
   that captures a different key is a bug, not a feature); (b) an `OptAsph=` above its
   element's `VarElt=` is ignored with no message (the parse is gated on `isVarElt`) —
   print one line naming the element and the fix ("OptAsph= must follow VarElt=").
   Gates: `cli_tests/cmd/Rx_header_unknown.jou` gains the look-alike key;
   a new `OptAsph_order.jou` (must-warn, and the asphere DOFs then absent from the solve).
5. **A trace still rewrites the prescription in two places** (§1.1b, your own 3b
   measurement: SAVE before vs after `opd nElt`, 403 decks): (a) `xObs` recomputed on 105
   decks (a different vector, not an ulp); (b) `Surface= Interpolated` SAVEd as `Conic`
   after a trace on four legacy decks (interpSurf, Int_Surf, MillsCross, MillsFilled) — the
   trace changes `SrfType`.  Each: which statement writes it, and whether the TRACE or the
   PARSER holds the truth.  For (a) the likely answer is that `ChkDf2`'s psiElt-based
   default and the trace's `define_local_csys` disagree on the obscuration frame — if so
   the parser should compute what the trace will use (one function, two callers), and the
   SAVE after a trace equals the SAVE before.  For (b) a trace must never change `SrfType`;
   find the assignment.  Gate: the 3b measurement itself, re-run — `xObs` 0 of 403 changed,
   `Surface` 0 changed; record `cli_tests/records/trace_rewrites_<sha>.md`.
6. **`ZernTypeL = 11` (ExtFringe) has no `ZerntoMon` converter** — it falls through the
   dispatch (warns in propsub, silent in tracesub/srtrace): a surface declared ExtFringe
   traces as the bare conic.  Decide by measurement: is ExtFringe a documented ordering
   the engine can convert (then write the converter from the textbook ordering, gate
   against Fringe on the modes they share), or not (then a NAMED load warning at
   `ZernType= ExtFringe` and the same fall-through, counted)?  No refusal without Dave.
7. **GRIN model 5** (§1.2 last item, four defects) is FEATURE-sized and frozen.  Do not
   start it; if an item above touches `GradCoef=`, say so and stop there.

Grating_example_001's compiler-sensitive edge rays (§1.1a) are resolved by item 3's
outcome — do not move the fixture before then.

## B. The bindings program (plan §1.4a, then §1.4b) — after A, or interleaved when an A item waits on Dave

1. **§1.4a — pymacos binds 136 of 146 api wrappers.**  Extend mmacos's
   `gen_mex_wrappers.py` (it already reads `macos_api_mod.F90` as the single source of
   truth — memory `reference_binding_codegen_source`) to emit the f2py layer and the Python
   veneer too: ONE generator, THREE targets.  Then the missing wrappers through it:
   `deadband_notes_get`, `fwd_root_flips_get`, `draw_rays_cmd/_get/3d_get`, `elt_info_get`,
   `elt_obs_get`, `elt_z`, `elt_srf_zrn_*`, `jmat_elt_get`, `polelt_get/_set`,
   `ray_hist_set`, `ray_pos_hist_get`, `ray_status_get`, `src_seg_get`,
   `met_calc/_geom_get/_get`, `translate`, `xps_cmd` — plus `modes_get`, `calib_stats_get`,
   `grid_idx_ovf_get` and whatever else landed since 10-10 (re-measure the gap first; the
   list above is the 10-10 count).  A pytest coverage gate: api routines listed FROM THE
   SOURCE, each must have a binding; must-fail = remove one.  A side-tree build recipe
   (`pymacos` against a `build_release_prefix*`) so a pre-fix red leg is routine — write
   it down in `pymacos/CLAUDE.md`.
2. **§1.4a defaults table + gate.**  Every session/user option the api exposes (OPD
   reference, FEX/SXP axis, `spcOption`, `ffcut`, obscuration option, the five modes,
   `remove_ptt` in the harvests, grid orientation, model size, `opd_mask` use): ONE
   default set in the engine at allocation, the three front ends inherit it.  Deliverables:
   a cmdref table (option, engine default, how each front end reads/sets it, the override)
   generated from the api source by the same generator; `tDefaults` (mmacos) +
   `test_defaults.py` (pymacos) reading the same deck and asserting the same answers; the
   rule in `macos_f90/CLAUDE.md`: a veneer never carries a default an api getter can supply.
3. **§1.4b — the cmdref BINDING column** (SMACOS command / direct setter / CLI-only /
   CLI-only legacy / missing), filled from the api source by the generator; the
   LINEAR-MODEL family and the differential commands are "CLI-only, legacy" by Dave's
   ruling (plan §1.4a text, 2026-10-10).  Then the missing group through the api in order
   of use: **single-ray tracing first** (RAY/SRAY/SRT — the gauge benches), the
   Zernike/analysis set second (ZCOEF/ZRM/SZCO/MAP/ELTS — item A.1 is its prerequisite),
   the image-plane post-processing third (veneer functions over the field the bindings
   already return; MULTISPEC/CADD through `compose`).  Each with generated mmacos +
   pymacos veneers and the cross-binding default gate.

## C. Small, yours when convenient

- `dyson_ladder` meniscus seeding (your 10-08 fix): confirm the gate exists and is in the
  suite (plan §1.4), tick it.
- `segment_rx` should emit `NSCount= 1` for single-hit segmented NS groups now that SAVE
  preserves it (plan §1.1b last item; resources) — CC will take this if you have not by
  the time A.2 lands; say which.

## Not yours

OPL's refusals, Mac-in-phase-4, the SAVE policy (Dave); the manual/cmdref words and the
§6.2 chapters (CCMac); the CLI journals, the weekly three-suite run, the gates on your
commits, COPHI's deck and the composite gauge deck (CC).  The `OpdReRayfLen=` misspelling
in the four IRIS/OPTIIX sandbox decks is CCMac's (he holds them).

## Reporting

One message to CC per landed item: the SHA (both repos if both moved), what the gate
showed red-then-green with the numbers, the corpus A/B count and the movers record if the
item touched a trace path, the fast-suite count, anything measured that overturned a
premise in this brief (the phase-1 record shows three of CC's did).  Dave hears the
summary from CC; anything that needs HIS ruling goes to him through CC with the
measurement attached, never a proposal alone.
