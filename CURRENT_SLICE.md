# Current Slice — in-flight working state

> **CC:** This file is the *ephemeral* working memory for the ONE sprint
> slice in progress.  It is the deliberate complement to the permanent
> record: `PLAN.md` / `PLAN_DESIGN_LAYER.md` hold *landed* state (sprint
> checkboxes, `CORE COMPLETE` blockquotes, the §10 Decisions ledger);
> the agent `MEMORY.md` holds durable learnings; `CLAUDE.md` (+ nested)
> holds mechanical gotchas.  THIS file holds only the half-done middle
> that compaction throws away — and it is **promoted then cleared** when
> the slice lands.  It is never a second source of truth.

> **CC (post-compaction / session resume):** read this file FIRST, then
> the plan section it anchors, then the CLAUDE.md set per the root
> directive.  If this file is at the empty template below, no slice is
> in flight — pick the next unchecked item from the active sprint.

> **CC (standing conventions, do not relax here):** work lands on
> `sls-dev` (macos) + companion `sls-dev` (MACOS_resources); every new
> function ships with a matlab.unittest test; `./run_mmacos_tests.sh
> fast` between edits, full suite pre-commit; every `matlab -batch`
> ends in `exit(0)`; each sprint tag ships a runnable worked-example.

---

## Active slice

> **LUIS SENS-TOOLS FIXES — IMPLEMENTED, LOCAL, UNPUSHED (2026-09-14).**
> Plan `~/.claude/plans/soft-hugging-peach.md` (Dave-approved + CCL/CCMac's
> six changes folded in). Luis reply draft `macos/DRAFT_email_luis_sens_groups.md`.
> Order run WS4→WS1→WS3→WS2. All targeted tests green; full `fast` = 469 pass,
> 3 fail ALL PRE-EXISTING (cross-platform bit-exact baselines: tPolElement
> material-flip = 1-ULP `isequaln` drift, err 2e-19; tSegMirMaker pie/hex2 =
> external SegMirMaker binary Hx bit-drift — neither is in my changed files).
> - **WS4** doc: `macos_f90/CLAUDE.md` Noll note fixed (Noll→ZerntoMon6 in all
>   3 chains; only ExtFringe(11) no-op).
> - **WS1** group dwdx speedup (`GroupedRigidBodyChannel` + builder + threaded
>   `group_smart_stop` default-on through dw_dx/dw_dx_multi/run_sensitivities):
>   Fix A skip-reaim-on-restore + one tail-settle/group; Fix B skip when all
>   members strictly downstream of the ENGINE-resolved stop, keep-when-ambiguous.
>   `tDwDxGroups` 15/15 (A/B gate-on==gate-off + StubGroupSession call-counts).
>   **Fix B engages only with an ELEMENT stop** (ApStop=/macos.stop); object-space
>   stop (default) → conservative, Fix A's ~1/3 still universal. Follow-up: resolve
>   object-space stop's effective element to widen Fix B.
> - **WS3** `zernike_grid_basis` +'noll'/'bornwolf' (verbatim ZerntoMon6/2 perm
>   tables → ANSI evaluator, same ndgrid orientation); `dw_dgrid` `zconv` recorded
>   in `out.zconv`. Fringe/NormHex/NormAnnularNoll deferred (error). `tZernikeGridBasis`
>   5/5 ((n,m) cross-check). Follow-up: non-segment engine round-trip per convention.
> - **WS2** engine: `PrtSingleEltInfo` emits `Link= <iElt>` when LnkElt>0;
>   `SAVE_KEYWORD_AUDIT.md` reclassified. Engine rebuilt (gfortran) + mex relinked.
>   `tLinkSave` 3/3 (emit-only-when-present; active+distinct-from-stock via OPD;
>   survives save→reload→save). `get_elt_csys` reads LCS not linked-reflector
>   motion → test uses OPD.
> **CCL REVIEW CLOSED (2026-09-14):** Linux gate green (fast 472/0, gfortran+ifx).
> WS1 follow-up landed: object-space stop (default) re-aim is a no-op (chief
> aimed at a fixed global point) -> now skipped UNCONDITIONALLY, so Luis gets the
> full gain with no knob (gated by test_smart_gate_preserves_jacobian on e5hex1).
> WS3 engine-exactness gate ADDED + PASSES both noll+bornwolf
> (tRunCompare/test_zern_grid_conventions_engine_equivalence; compare within the
> confined disk -- generator matches zern_seg_eval corr 1.0; norm applied
> post-reindex so bornwolf scale is exact). tLinkSave+tZernikeGridBasis added to
> SUITE_FAST. Also: 2 e2e6m_r2 template files were already modified in the
> resources tree (NOT mine).
> **CCL ROUND 2 (2026-09-14 eve) CLOSED:** Linux fast 481/0 on the push. WS3 gate
> STRENGTHENED to have teeth -- CCL showed Noll idx 8 and B&W idx 8 both land on
> ANSI slot 9 (coincide), so the mode-8-only gate couldn't see a convention swap.
> Now modes {4,7,8} per convention at ng 256 / sampling 31 (rays coarser than
> the grid pixel -> facets vanish, 1%/corr 1.0000), thresholds 1%/0.999, PLUS a
> negative control at mode 7 (Noll 7 coma vs B&W 7 trefoil; other-convention map
> vs this deck's poke must give |corr|<0.5). tRunCompare 4/4. Luis email UNBLOCKED
> (added the note: Noll & B&W coincide at 1-3 and 8). NEXT: Dave's new thread =
> predict as-built performance (surface errors, error budgets, drift, electronics
> noise) -- awaiting his steer on static-budget-first vs dynamic-loop-first.

> **TO / GAUGE-DECK PDI LANE, IN FLIGHT 2026-09-13 (`BRIEF_to_gauge_deck.md`;
> plan `BRIEF_gauge_deck.md`).  CODE COMMITTED, resources `5d99c20`, LOCAL.**
> New dir `mmacos/templates/40_benches/pdi_dm96/` (pdi_params / pdi_run /
> pdi_batch.sh / pdi_vfig_util + the psri decks and figures + README +
> REPORT_gauge_pdi.md + runs/); the P/PF README section MOVED there with a
> pointer left in zwfs_dm96; pre-split records stay in zwfs_dm96/runs.
> **Built:** (1) `dmg_pdi_gauge` REF_SHAPE 'deck' -- the P/SRI reference arm
> TRACED through psri_ref.in's physical pinhole, test through psri_test.in
> (+`pdi.ref_frozen`, the control that separates the arm's SHAPE from its
> MOTION); `zwfs_run` gains a self-contained `stage_bench_psri_`
> (`pdi.bench 'psri'`, readings must be {'PF'}).  ZWFS path BIT-IDENTICAL
> (pdi_dev3's G3/G4/G5/G6/G7 reproduce; G4 = 0.296 pm).  (2) `dmg_loop`
> knobs `start_rms`/`start_shape` (the DESCENT: the DM's initial figure),
> `recal_every` (+ the `ins.recal(cmd) -> struct('est','nstates')`
> contract CCMac mirrors), `intra` (drift WITHIN a scan, `aux.dstep`),
> `ref_walk` (non-common-path reference phase, `aux.ref_phase`); gates
> tDmgLoop G9-G12, 14/14.  A descent's STARTING matrix is measured on the
> STARTING surface.  Drift increments now drawn once ahead of the loop in
> the same order -> every earlier run reproduces bit-for-bit.
> **Measured so far (dev res, model 512/65, 48x48 DM):** the traced arm's
> reference under the 30 nm surface -- total change 0.1528, scale |0.847|,
> **SHAPE change 0.0057** (the FFT surrogate at 2 lam/D gives 0.1539 /
> 0.8463 / 0.0059: the real arm behaves like the surrogate); G5 absolute
> 31.06 pm = 0.257% of a 12 nm figure where the synthesized reference reads
> 0.000 by construction; rows single 0.9939 / 5 pm, grid 0.9936 / 2, dense
> 0.9949 / 164; N(1 pm) 3.15e13.
> **RECORD LANDED (resources `215a452`): deliverable 1 COMPLETE.**  At
> model 1024 / 193 rays / 96x96 with the matrix on the 30 nm surface:
> the traced arm's G5 absolute error **5.889 pm** on a 13 nm figure vs
> **0.000 pm** frozen and 0.000 synthesized -- so the error is ENTIRELY
> the reference MOVING, and the synthesized LP01 model IS a frozen
> reference; differentially the motion is a 10%-class noise-floor effect
> (dense 129 -> 144 pm) and nothing on gain or range; rows single
> 0.9935/2 pm/SNR 4895, grid 0.9924/1, dense 0.9926/144; **capture range
> 480 nm+ with gain inside 0.7%** (the synthesized model is +12.8%
> there; S and P fold at 120 nm); N(1 pm) 2.00e14 vs 1.9e14 -- the
> reference model does not set the photon cost, the 60/40 pickoff does.
> Deliverable 6 DONE (both layouts + psri_render redone in the recipe,
> `pdi_vfig_util`; parts lists in the README).  Gates: tDmgLoop 14/14.
> **UPDATE `BRIEF_to_capture.md` (Dave/CCL 2026-09-13): the wrap finding
> accepted; deliverables 8 and 9.**  D8 DONE: `dm_gauge_lib/dmg_unwrap.m`
> -- masked 2-D least-squares unwrapping (Ghiglia & Romero 1994: the
> unweighted Poisson solve by MIRRORED FFT, since `dct2` is a toolbox
> function and the release gate forbids external deps; then their sec-5
> PCG refinement on the WEIGHTED equations so the region outside the mask
> cannot pull on the answer inside it), with RESIDUES counted so a map
> past the pixel-gradient limit is reported not silently wrong.  Gates
> tDmgLoop G13, 15/15 (ramp 8.5e-14 rad; 1.5 / 3.0 waves PV on a disc
> 6.7e-13 / 1.3e-12, zero residues; 40 waves PV -> 644 residues).  Wired
> behind `battery.unwrap` (DEFAULT FALSE) for S/V/P/PF in BOTH the
> measurement differential and the calibration's class maps;
> `loop.unwrap 'auto'` turns it on exactly when `loop.start_rms` is set;
> L/F/I/I+ untouched.  Non-disturbance: `runs/uwoff_ref` bit-identical
> (v3dev G4 = 0.296 pm).  `loop.start_rms` now takes a VECTOR and the
> start matrix is measured ONCE per start for EVERY class (was once per
> reading) -- what makes the ladder affordable; the stage prints the
> OPENING differential's wrapped rms / residues / max gradient /
> unwrapped rms per start and reading.  **GATES GREEN: tDmgLoop 15/15,
> mmacos fast suite 469 pass / 0 fail.**  **PUSHED 2026-09-14** (Dave
> ran it: "I pushed them, keep going"; the session's auto-mode classifier
> refused `git push` throughout 2026-09-13).  Both repos in sync with
> origin; CCMac unblocked on the shared dmg_loop knobs.
> **DELIVERABLES 1, 2, 4, 5, 6, 7, 8, 9 DONE; 3 (the pinhole diameter of
> record) running as gseq2.**  Overnight, all exit 0:
> **CAPTURE (D9):** largest initial surface brought to 3 pm -- L 30,
> S 30, V 60, P 60, PF 60 nm, IDENTICAL with the unwrapper off and on.
> The wrap premise is right for PF ONLY: past ~60 nm S/V/P return 24-28
> nm whatever the truth is, with ZERO residues -- BLIND, not folded.
> Unwrapping AND on-surface recal TOGETHER raise PF and P-with-shutter
> to 100 nm = 200 nm WFE (0.245 / 0.207 pm, against 64095 with neither,
> 5383 unwrap-only, 63520 recal-only), at 1e13 as well as 1e15 --
> reference-limited, not light-limited.  Cadence is not the constraint.
> **V4 within-scan drift:** the DM's HELPS the stepped readings 26-32%
> (a scan reads the surface at its MIDPOINT = a free half-step of
> prediction); the camera's HURTS because it is additive BIAS.  Opposite
> signs, same knob; L and V unchanged to the digit.
> **D5 reference-arm walk:** 100x of walk costs 4% in photons -- a
> path-length change between arms is a PISTON, the mode the estimator
> nulls.  Benign for a DM servo, NOT for absolute E-field work.
> **Two bugs fixed:** the per-start calibrations held ~750 MB each and a
> six-rung ladder was OOM-killed (build one at a time, drop before the
> next); and `tr -d '-.'` reads the leading '-' as an option, collapsing
> three reference-walk runs onto one tag.
> **RUNNING:** gseq1's tail (pfdeck_loop, cap385p*, noise193p_b*), then
> `runs/gmaster3.sh` -> gcap (D9, the ladder BOTH WAYS: starts 30-300 nm,
> 193 rays and K 40, both departures stated), gseq3 (intra), gseq4 (ref
> walk), gseq2 (pinhole).  Superseded: gmaster.sh, gmaster2.sh.
> **[old plan below]**
> `runs/gmaster2.sh` -> gseq4 (D5 reference-arm walk, ~75 min), gseq2
> (D3 pinhole diameter, the model-2048 leg ~4 h), gseq3 (D4 descent +
> within-scan drift, ~7 h).  NEXT: fill REPORT_gauge_pdi.md sections
> 2-5 and 0 as they land; push only on Dave's review.

> **PDI READINGS (Dave 2026-09-12 "look at another sensor, using a
> point-diffraction IFO approach -- see Brandon Dube"), IN FLIGHT, LOCAL.**
> Plan + literature: `BRIEF_pdi_campaign.md` (the paper meant = Dube,
> Nejadriahi, Sidick, Jewell, Redding, Lou, Basinger, SPIE 13092-178
> (2024), "Absolute and differential complex E field reconstruction by
> phase shifting interferometry" -- non-common-path IFO with a photonic
> phase shifter; full text not retrievable here, ask Dave for the PDF).
> Built: `dm_gauge_lib/dmg_pdi_gauge.m` (readings P = stepped pinhole at
> the FocalMask, common path, surround t; PF = fiber reference with
> photonic steps, exact) threaded into `zwfs_run` as classes 5/6 (P.pdi
> knobs; opt-in, record defaults untouched), gates G5-G7 in the bench
> stage, figs palette.  Dev gates PASS (`runs/pdi_dev`, model 512):
> G5 P 1.9 pm of 11.8 nm, PF 0.000; G7 5e-15 (P at t=1 == S).  RECORD
> SEQUENCE RUNNING in background: `runs/pseq.sh` -> pdi193 (flat
> battery + noise), pdi193base, pdi193d1 (1 lam/D on the base), ploop193
> (loop P/PF).  **UPDATE 08:35: Dave pointed at the paper
> (`~/dev/MACOS_sandbox/pSRI/Dube_pSRI.pdf` + their MATLAB model); PF
> re-pinned to it (LP01 fiber-mode reference with per-state coupling
> kappa, pickoff 0.6, SH5 scheme option, step_err) -- resources 09c40da;
> the first sequence was stopped after pdi193 (record 1, PF = pinhole-
> shaped idealization, numbers in the brief); `runs/pseq2.sh` is the
> live sequence: pdi193f, pdi193fbase, pdi193state, pdi193d1,
> pdi193se_ls/sh5, ploop193 (~2.5 h from 08:31).**  V2 code committed (resources 7323d80); V2 runs v2g_*/
> v2e10a done (0.1 rad retardance error, uncalibrated: on-surface single
> 0.9934/4 pm -- the surface matrix absorbs it) -- README V2 still to
> write.  NEXT: fold numbers into README "P / PF" section + brief +
> memory; commit; push only on Dave's review.

> **2026-09-17 (TO): BRIEF_to_tg_redo PACKAGE A -- THE BENCH COLLIMATED FOR
> REAL, THE TAIL TUNER PUT ON THE READING.  Committed LOCAL, resources
> `23131f6` (A0 was CCL's `37193e5`).**  New tool `tg96_collimate.m`
> (runs/coll_lens): solves the collimator's RADIUS with its conic, the
> focuser's conic, and the FocalMask seat, against the engine's rays, over the
> rays that REACH the DM.  **THE FINDING: the lens rig's collimator was not
> merely fed 25 mm inside its conjugate -- its RADIUS was matched to that wrong
> conjugate.**  `L1_Kr` 236.866 = (n-1)*473.7 = (n-1)*(F1 - zsource), so with
> the source at F1 the lens has 5% of surplus focal length; a solve holding the
> radius runs the conic to -4.68, walks the source back 27 mm (undoing
> SRC_AT_FOCUS exactly) and still leaves 2.6e-4 rad rms.  Solved:
> **L1_Kr 249.246312, L1_Kc -0.583016, L2_Kc -0.581843, MASK_TRIM +1.231759**
> -> exit-ray spread **1.27e-3 -> 6.3e-9 rad rms** (47 waves over the beam ->
> 0.0), focal spot **0.17 um rms** (lambda F/D 2.8 um), marker ON the ray focus
> (was 10.4 mm off -- the same error the zwfs sheet carried as the constant
> MASK_TRIM -5.582 since S1; now `'scan'` there).  **The CONIC barely moved**
> (-0.583016 vs the record's -0.5829): a conic belongs to the shape and the
> plano orientation, not to the conjugate.  Beam now 58.3 mm at the 48 mm DM
> (the DM IS the stop, 68% of the grid through); `P.clear.beam_r` 56 -> 59; the
> DM leg (450) and BS angle (22.5) do not move.  **Builder trap closed:**
> `stage_B_` forwarded SRC_AT_FOCUS from INSIDE its `if optics=='oap'` block,
> so setting it on the lens rig was accepted by the sheet and silently dropped
> (the trap zwfs_params records for MASK_SUB); SRC_AT_FOCUS/SRC_TRIM/MASK_TRIM
> now forwarded on both rigs by tg96_run AND tg96_tail.  **Tail tuner (Dave's
> ruling):** `P.tail.objective 'reading'` -- the cost IS what the pupil stage
> measures the camera recovering off the DM (tg96_pupilsim stage 2, or its
> stage-1 band-edge-phase proxy, default, ~70 s/eval), the null computed and
> PRINTED every evaluation and never optimized; `P.tail.free` holds the field
> lens at the geometric seed station (lens: {FL_Kc,DET_TRIM}; mirror rig:
> {DET_TRIM} alone).  tg96_pupilsim gained 'stages'/'figs'; a failed evaluation
> now prints its error instead of a silent 1e6.  **Seed tail on the collimated
> bench, before any tuning:** band-edge phase 0.0122 rad rms / 0.0237 max (gate
> <0.06), phase gain 0.9997 worst (gate >=0.998), image surface flat to 0.7 mm,
> distortion 0.023 mm rms (the brief's 0.01 came from the UNcollimated bench;
> the same table already showed collimation RAISING distortion), null 59.3 nm.
> **IN FLIGHT:** `runs/redoseq.sh` (detached chain, ~3 h from 10:24): lens tail
> tune (tag lens96) -> mirror tail (oap96, DET_TRIM alone) -> both rigs
> re-emitted as `redo_lens` / `redo_oap` (stages clearance/bench/figs, on
> `redo_<optics>_tail.mat` copies so the record's lens_tail.mat/oap_tail.mat
> stay put until the gate passes) -> the pupil stage on each EMITTED deck =
> package A's gate record.  Watch: `runs/redoseq.nohup`, `runs/tail_*.log`.
> **NEXT:** harvest the chain into REPORT_bench_realism section 8 (8.3/8.4),
> copy the accepted tails onto `lens_tail.mat`/`oap_tail.mat`, then package B
> (tg96_pupil_s2s leg by leg on the improved bench: the collimated legs to
> <0.02 wave, the quartet's Airy spot and its far-side pupil radius -- 20%
> small today, suspects in BRIEF_to_tg_redo section 2 -- the exit step's zElt
> convention against a 10 mm known defocus, then the readout vs tg96_pupilsim).
> Then package C (the record's runs redone, the Mac carrying the mirror rig).
> Report: REPORT_bench_realism section 8.  Rules in force: one MATLAB at a
> time, kill by PID, no push.**

> **2026-09-17 (CCL): THE PUPIL IMAGE SIMULATED (Dave's ask), THE DM MADE THE
> STOP, RIY runners.**  `tg_psi_dm96_oap/tg96_pupilsim.m` (runs pupilsim_lens /
> pupilsim_oap; report section 7; deck slide 13 of 50; sheet block P.pupil;
> RIY `tg96_pupil_batch.{m,sh}`; README block).  Stage 1: the leg's coherent
> PSF per DM zone from the rays (41 tilts, 2-D over the band; the DM as the
> stop keeps each ray on its zone to 2 um; the zone wavefront = defocus +
> astig, higher orders 0.3 nm).  Stage 2: the DM field through the zone
> PSFs (overlap-add), four-step readout.  RESULTS: every mode observable on
> both rigs; lens rig as tuned Nyquist gain 0.99 center / 0.95 worst,
> half-Nyquist and below within 0.3%, pokes 0.99, 30 nm surface back to
> 1.2 nm; mirror rig 0.997 / 0.23 nm.  THE LENS RIG'S DETECTOR IS NOT AT
> THE PUPIL IMAGE: 2.6 mm ahead on axis, 4-6 at the edge (the null-tuned
> tail is blind to pupil defocus); +4.3 mm -> 0.993 worst, 0.4 nm.  THE BEAM
> OF RECORD WAS THE SOURCE CONE (77 mm lens / 82 mirrors on the 96 mm DM;
> nothing clips): Dave ruled the baffle must not constrain + an aperture on
> the DM -> sheets changed (tg96_params, zwfs_params: R_BAFFLE 12.5->18,
> D_LENS 60->66, R_TO_AP 30->28 = 48 mm); the simulation opens the baffle,
> widens the cone x1.06 and puts the 48 mm aperture on the DM.  ENGINE
> PLANE-TO-PLANE CHECK (Dave: reference surfaces in the .in, the CTB model):
> `tg96_pupil_engine.m` built (mask sandwich + a sphere concentric with the
> beam after the FL, NFS1surf to the detector); MEASURED that a
> geometric-to-physical hand-off is exact only at a pupil conjugate and this
> leg has none before the detector (the DM's image by L2 is virtual, 876 mm
> past the focus) -> the seed is the undiffracted DM pattern; the sandwich
> centered on the seat marker (5.5 mm off focus) seeds 54 waves of defocus.
> The station-to-station form was BUILT (`tg96_pupil_s2s.m`: 4 collimated
> NFPlane pairs, the Rx_Coro quartet with asymmetric radii, the lens per
> index, an exit sphere + NFS1surf, 4 zElt conventions x a 10 mm defocus
> discriminator).  Lens rig: cannot work -- its collimated space carries 41
> waves (collimator fed 25 mm inside its focus), the rays walk off the grid
> between NF legs (1.5 waves of spurious aberration at S1, a 10x focal
> spot).  Mirror rig (fed at the focus): runs but wrong -- the quartet's
> far-side pupil 20% small vs the rays; no convention gives a gain map.
> NEXT (half a day): the focal spot after the entrance sphere, the far
> sphere's pitch bookkeeping, the exit step alone with a known defocus.
> AFTERNOON: Dave asked how to IMPROVE the pupil images -> `tg96_pupil_options`
> (five lens-rig variants, REPORT_bench_realism 7.1, deck slide 15, brief
> section 8): THE TAIL GEOMETRY IS THE WHOLE STORY -- the seed field-lens
> station (10.8 mm past the focus) images the DM flat (0.9992 as built,
> 0.9999 with a 0.7 mm move, distortion 0.003 mm) where the null-tuned tail
> (39.8 mm, conic -2.59) gives 0.954 / 0.27; a flattener cannot fix the
> tuned bowl (5 mm-radius image surface); true collimation (the source is
> 14 mm inside the hyperbolic collimator's focus) halves the distortion and
> is needed for re-tunes + the physical-optics chain, not for the gain.
> Recommendation: hold the field lens at the seed station, re-tune conic +
> trim with the image surface in the objective, detector at the image mean;
> mirror rig the 0.6 mm move; fix the collimator first.  Deck 52 slides;
> Dave's edit copy was OPEN -> NOT re-synced (run sync_edit_deck.sh
> deck_gauges after closing it).  Committed: resources 02929cc 50358d4 +
> this; macos c13d4c5 fd56d45 + this.  Push on Dave's word.
> LATE AFTERNOON: Dave's rulings on BRIEF_to_tg_redo (substrates YES as
> proposed; tail objective = the IFO's own reading, null reported; keep the
> L2 re-solve; larger-pixel camera priced both ways; the Mac assigned to the
> redo's runs) written into the brief (sections 1b, 6) + a read-up section
> for TO (clear over compact).  The Mac's cycle 2 pulled (resources 12d03b0,
> rebased): oapsens385 (mirror rig at 385 px: capture S 36 / V 59 / P 52 nm,
> rows S 0.9911/4 pm) and oapcap22 (sensors' descent on the mirror rig:
> V/P 3 pm from 60 nm in 21/19 cycles, diverge from 100, S never) folded
> into the deck.  Deck 53 slides: parts slide 6 with the decided substrates
> + the seed field-lens station + the 9.8 mm pupil image + camera options;
> NEW slide 41 'Next: the interferometer redone on the improved bench';
> edit copy synced.  INTERIM DECK for Dave (~5 h from 15:00): this build.
> NEXT for TO: BRIEF_to_tg_redo packages A0, A, B, C; CCL writes the Mac's
> cycle-3 runbook once A lands.  Friday: the interim REPORT (current record,
> labeled), style gate, .docx.  Standalone paraxial chain (`fourier` true) works as
> a standby (flat pupil 4.81 vs 4.92 mm).  Interim report paragraph updated;
> BRIEF_pupil_quality_tg section 7 (three rulings for Dave: re-run the
> record on the 96 mm beam now or after Friday; DET_TRIM +4.3 / tuner
> objective; the s2s chain now or after Friday).  Final reruns of both rigs
> (corrected scoring) in flight at the time of writing; commit after.**

> **2026-09-16 evening (CCL): PUPIL IMAGE QUALITY for Fang Shi DONE (due
> 09-17).**  `tg_psi_dm96_oap/tg96_pupilq.m` (resources fd24ad9; runs
> pupilq_lens / pupilq_oap; REPORT_bench_realism section 6; deck slide 12
> after the IFO station slide, 49 slides; memory project_pupil_quality_tg).
> The DM is the STOP (Dave): bench as built, field = source shift at the
> collimator's focus; crossing cloud at the camera.  Distortion vs one
> affine 0.13/0.30 mm (lens), 0.45/0.85 mm (mirrors); blur 0.02-0.03 mm
> rms over the actuator band; pupil surface sag 3.8 mm (lens) / tilt 0.45
> (OAP); focal WFE lens 1.1-2.7 nm, OAP 0 on axis / 144 nm = 2.3 lam F/D at
> the band tilt (coma).  Traps: collimated source at the DM = a
> nonexistent bench (16 lam F/D); the seat marker is not the focus; the
> exit-pupil sphere idiom refuses to load via the mex (fitted focus removal
> instead); `range` is a toolbox call.  Mac cycle 2: oapnoise22 folded
> (photons for 1 pm: S 5.8e13, V 7.3e13, P 8.5e13); oapcap22 running;
> oapsens385 queued.  Interim report has the pupil paragraph.  Unpushed:
> macos since 303c209 (brief, deck, report), resources fd24ad9.**

> **2026-09-16 midday (CCL).**  Deck 48 slides (contents slide 2; Dave's
> editing copy synced via the gate; sidecar deck_gauges.geo.json).  TO's
> close-out landed on the slides: item 2(a) CLOSED (the raw four-step wraps
> between 60 and 120 nm on BOTH rigs: a base past lam/2 reads 91 nm =
> 316.4/sqrt(12); the lens "never flags" was the 7-deg record); item 4 rows
> on the tuned tail (10 mm plates: rows 1.00-1.01, reading attenuated 13% =
> floor/SNR not gain); item 6 station figures at 1800 px (stnoap on the
> deck; stnlens has a 62 nm misregistration, open); item 3 CLOSED (winner
> gate relative to the seed, enforcing).  MAC: job A died (no flock on
> macOS -> wrappers guarded, resources 3c88b46; oap_fold_batch.sh left for TO,
> it was executing); job B oaploop22 (Dave pushed 0cba0a6): sensors' servo
> on the redesigned mirror rig S 2.7e12/7.9e12, V 1.7e12/5.7e12, P
> 2.5e12/7.4e12, thermal floors 10 pm -- within 10% of the lens rig -> servo
> slide + measured slide + recommendation + report.  INTERIM REPORT full
> first draft (5200 words, sections 1-2 by subagent with a source-notes
> appendix; 3-5 CCL) at REPORT_interim_ctb_e2e6m_gauge.md.  Local unpushed:
> macos commits since a9eddc5; resources 3c88b46.**

> **2026-09-16 ~08:30 (CCL): DECK 47 SLIDES for Dave's JPL preliminary
> (macos this commit): Dave's two edit-copy changes folded ("The DM";
> slide 8 left block 15 pt via NEW sidecar deck_gauges.geo.json); the
> overnight results in: servo table (mirror rig redesigned 5.4e12 / 1.7e13
> / 10 pm; record 7-deg row kept), capture tables (redesigned captures from
> 200 nm unwrapped; raw fold 120 vs lens 240 = ladder arithmetic), NEW
> slide "Interferometer, station by station" (oapifol2_stations.png), NEW
> slide "The reflective front end, measured" (record / redesigned / lens
> table + vector amplitude-imbalance + open items), slide 3 for/against
> updated, recommendation front-end line = "either builds; ruling
> pending", camera 7.5 mm / 1155 px / bin 3 or 6, splitter 10 mm modeled
> (null 0.13 -> 20 nm, rows pending item4bseq), polarization at the built
> angle (1.56 deg, +0.9%), provenance tags.  Dave's edit copy is OPEN in
> Impress and NOT re-synced (gate refuses); he must close it, then
> `./sync_edit_deck.sh deck_gauges --folded`.  CTB doc pass DONE by
> subagent, committed resources d28fdeb (README 47% fill + Aperture
> convention, example_ctb label, deck_ctb slides 1/9, CTB_PROP_STATUS
> section).  Dave's CTB ruling: (1) docs now; (2) regenerate at 42.75 mm
> AFTER the gauge results clear.  Claude update pending restart; this
> slice is the handoff.**

> **2026-09-16 morning (CCL): TO's overnight results harvested (git log
> 09-15 18:00 -> 09-16 06:39, 20 commits; every run exit 0 except oapdesc2's
> figure bug after its descents landed).  HEADLINE: the CTB DM-model "slip"
> is RETRACTED by the traced footprint (beam 21.24 mm; the 32x32/0.67 mm
> lattice spans it; no EFC result affected); the real finding is the
> point-source Aperture convention (engine: FULL cone angle, sourcsub A =
> Aperture/2; example_ctb reads a half-angle -> 47% fill, README says 95%)
> -> Dave decides keep (CCL recommends) or regenerate; the 64x64/1 mm ruling
> withdrawn; brief 0b, note, memory, deck_ctb banner corrected (macos
> 1b1bdbe, resources ddaca96).  Vector pair on the reflective rig: G4 634
> bare / 319 lam/4 overcoat / 0.054 pm calibrated (PASS); the variable is
> channel AMPLITUDE (phase fixed at 0.062 rad); rows hold on every coating;
> fold lever unpulled.  Servo (oapifol2, 854 states, 5.5 h): 3 pm from
> 5.4e12 noise-only, 1.7e13 under the 2 pm walk; thermal floor 10.0 pm =
> rate/(gG), LOW-ORDER (9.2 of 10 pm below 4 cyc/ap), not photon-limited;
> both floors analytic to 1%.  Descent (oapdesc2, unwrap on): 100 and 200 nm
> starts reach 3 pm at cycle 17/19, rho 0.50/0.52.  Wrap: the record's lens
> ladder was a 7-deg bench; like-for-like the lens rig ALSO breaks (at 240,
> OAP at 120); wrap column identical; candidate = pixels per unknown in the
> solve (37 vs 41 px); the 'wrap' stage (dA vs dD, n_cross) staged;
> battery.unwrap never reached the ladder (vacuous) -- fixed as a stage.
> Realism: 10 mm plates cost the flat null 0.134 -> 20.09 nm retuned (the
> 1.39 mm shear), gate ate the winner -> item4bseq; the sensors' rows with
> substrates (sub22) and thickness (thk22) = gate22_193 to the digit
> (S 0.9888/5, V 0.9938/3, P 0.9937/3; capture 37/62/54); MASK_TRIM 'scan'
> found +0.86 mm; pupil image is 7.5-7.8 mm not 9.4 (binning 6 not 4 on the
> ZWFS rig).  Snapshot polarization at the built angles: arm rotation
> +/-1.56 deg, gain within 1% (record 7.48 deg / 11.7%); the 'corrected'
> column is a no-op by construction (rigid analyzer rotation cancels in a
> differential).  Tail-tuner gate FAILED its two-leg test (point sample
> tracks magnification) -> advisory; item 3 NOT done.  IFO station figure
> works (0.00 pm flat / 626 pm on 30 nm = 2% absolute); lens-rig figure +
> 1800 px re-run queued.  Field servo step 1 prescription written (dichroic
> at Apodizer_Pst, 300 mm lens F/9.4, 2 lam/D dimple 11.8 um).  Queue now
> chains on run artifacts (closefinal.sh).**

> **2026-09-15 ~18:00 (CCL).  TWO DEADLINES.**  (1) CTB DM-MODEL DEFECT (TO):
> `ctb_dm.m` beam_d_mm 21.3 used as the DIAMETER; the deck's beam is 21.4 mm
> RADIUS (Aperture 8.485e-3 x 2519 mm) -> the lattice covered the inner 28%
> of the pupil; Jacobians/EFC/README inherit; the deck_ctb slides 9-13 are
> self-consistent (aberration-free) but must be RE-SCORED.  RULED (Dave): beam stays 42.75, DM = a real 64x64 at 1 mm pitch, active =
> influence inside the beam, reach 21 cycles.  Ordered first
> in `BRIEF_to_gauge_close.md` section 0b (TO, by THU 09-17 NOON: probe
> committed, defaults fixed in ctb_dm/ctb_dm_jacobian/ctb_efc_physics +
> callers audited, tCtbDm gate, ctb_study re-derived jac/efc/relin/physics/
> bandwidth/vvc ONE MATLAB at a time, old/new table in CTB_PROP_STATUS).
> deck_ctb carries a DRAFT banner (resources, this session).  My note's
> separability numbers corrected to the real 1.34 mm pitch (12% at Nyquist:
> NOT separable in-band on the CTB either).  (2) INTERIM REPORT (CCL, due
> FRI 09-18): `REPORT_interim_ctb_e2e6m_gauge.md` skeleton + day plan (Wed
> draft 1-3 from the decks; Thu TO's numbers in + section 4; Fri style
> gate, .docx, Dave).  Deck_gauges: slide 3 both layouts + for/against,
> plain-language pass, slide 15 'magic' steps.  TO's ladder: the tg96
> ladder break is the wrapped-absolute subtraction (beyond-fold fraction is
> the meter; unwrapped capture is the number of record); control lensuw2
> running -- read it by the MEASURED fraction per rig (lens raw ladder held
> at 120 where the OAP broke).  Local macos commits f117c42..this; resources
> c9f92f6 + the banner + TO's.  Push on Dave's word.**

> **2026-09-15 ~16:00 (CCL).**  TO's reflective lane: items 1-3, 5-7 done;
> item 4 reads on the SEED tail (oapifo2: rows 0.9915/2.4 pm, 0.9905/161,
> cross-talk 0.006-0.02, null 28.9 pm; the break ladder WRAPS at 120 nm --
> smaller capture than the lens rig); the tail tuner's objective fix was
> DISPROVEN (objwin3 scores better and reads 0.0338) -- open why; item 5:
> pinhole RECOVERS (0.269 pm), vector pair raw G4 634 pm bare Al / 208 no
> coating = ~0.1 / 0.035 rad channel phase on the V3 scan where the
> on-surface rows still hold (my reply in BRIEF_to_reflective).  NEW BRIEF
> `BRIEF_to_gauge_close.md` (7 items: vector rows + lam/4 overcoat + verdict;
> item 4's loop/descent + the 120 nm wrap; tuner gated by a battery row;
> realism 3-5, 6, 8; the coronagraph FIELD SERVO on the CTB deck -- Dave's
> C+ concept, `NOTES_gauge_in_coronagraph.md` with the clearance scan tool
> `demo_session/gauge_in_coro_clearance.py`: face-on blocked in both
> packages, out-of-plane 15 deg clears the CTB, nothing clears the space
> relay; DM separability 33 cycles on the CTB, never in flight; 10 pm per
> 200 s on V=5).  Deck 45 slides (slide 13 pol-camera bullet, redesigned-rig
> slide).  Local macos commits since Dave's push: f117c42..this; resources
> c9f92f6 + TO's.  Push on Dave's word.**

> **2026-09-15 ~10:00 (CCL).**  TO's reflective items 1, 2, 3, 7 are
> COMMITTED (resources 5cb2651 layout + parts list, 5c170c0 P/SRI
> clearance PASS at 22.5 +36.9 mm / FAIL at 7, 69033b6 tail + seat; macos
> 62248d5).  THE STORY CHANGED: there is no fold coma -- the builder fed
> the collimator 25 mm inside its focus (zSource; 926 urad residual, blur
> linear in the fold angle, the 6.14 mm seat trim = that error refocused);
> fed at its focus the seat blur is 0.000 lam F/D at zero trim, and the
> designed rig's UNTUNED tail nulls 0.0289 nm (record's tuned OAP tail
> 12.9, lens rig tuned 0.134).  The blocker was the input polarizer 10 mm
> past OAP1's pole (in the cone + the sag envelope) -> moved to the source
> leg, ray loss 0.  DECK (45 slides, macos this commit): the
> lenses-vs-mirrors slide retitled "the record" with the cause corrected;
> NEW slide "The reflective front end, redesigned" (TO's oapdraw3 layout
> crop + record-vs-redesign table); recommendation's front-end bullet
> hedged; status line + provenance tags.  RUNNING in TO's lane:
> oap22d_tail (model 512) then ifoseq.sh -> oapifo (rows) + oapifol
> (loop); then the mask sensors on the rig (item 5; needs
> bench.MASK_TRIM 0).  When those land: replace the record's rows on the
> two slides and re-decide the recommendation's front-end line (Dave).**

> **2026-09-15 ~09:30 (CCL, post-compaction).**  The slice's TO DO is
> DONE: gate22_193 folded into the zwfs_dm96 README (Findings, first
> bullet: "Bench of record = the 22.5-degree splitter", rows vs the 7-deg
> record in parentheses, capture S 37 / V 62 / P 54 nm at 193 rays) and
> into the deck's side-by-side footnote + provenance tags (stations193,
> gate22_193); runs/gate22_193 COMMITTED (resources, this session's first
> commit); deck rebuilt (44 slides), render-checked.  TO IS ALIVE on
> `BRIEF_to_reflective.md`: items 1 and 2 committed (9a8d699 = item 1,
> pushed by me at 08:30 as part of Dave's push; faab9ce + 180258b = item
> 2, LOCAL: OAP1 20 deg / OAP2 25 deg, sides +1/-1, polarizer in the
> source leg, output optics 125 mm ahead of OAP2, worst +33.4 mm over 8
> parts) -- `tg_psi_dm96_oap/REPORT_reflective.md` has the status table
> (items 3-7 not started); a model-1024 tg96 run (oap22d, stages
> bench/figs/clearance) is LIVE in TO's lane, so NO model-1024 MATLAB from
> this lane until it exits; TO's tree edits (oap tg96_run.m, runs/oap22d,
> *.nohup) are THEIRS -- not staged, not touched.  NEXT (CCL): when TO's
> item 3 layout figure lands, redo the deck's lenses-vs-mirrors and OAP
> slides from it; the deck DRAFT awaits Dave's sign-off; push on Dave's
> word only (resources has TO's 2 + my 1 ahead of origin).**

> **HANDOFF 2026-09-15 ~08:30 (CCL, autocompact).**  PUSHED both repos at
> Dave's word (resources 1338920, macos c8a5472); local since: macos
> ccacaee (briefs re-addressed).  The 07:15-07:22 tg96 node-solve edits
> were MINE (lost from memory in the OOM crash), gated and committed as
> item 2 (1338920).  CCMac will NOT return: TO owns the realism brief's
> items 3-8 AFTER `BRIEF_to_reflective.md`; the queue is in that brief's
> section 6.  DONE: zwfs_dm96/runs/gate22_193 (exit 0): on the 22.5-deg bench the
> rows HOLD -- single 10 nm S 0.9888/5 pm, V 0.9938/3, P 0.9937/3 (7 deg:
> 0.9885/5, 0.9935/4, 0.9935/4); grid 1 nm S 0.9982/4, V 0.9980/3, P
> 0.9980/3 (0.9993/4, 0.9992/3, 0.9992/3); dense V 0.9986/278 pm, P
> 0.9972/287.  TO DO: fold into README (a 'bench of record 22.5' line)
> + the deck's side-by-side footnote; commit runs/gate22_193.  Deck 44 slides DRAFT at
> demo_session/deck_gauges.md (22.5 deg layouts, three jobs, parts x3,
> station figures S/V/P, analysis path, migration planned).  Rules learned
> today: one model-1024 MATLAB on this box, every wrapper waits (tg96
> fixed); a session's uncommitted edits survive its death in the shared
> tree -- read `git status` before assuming another lane's work.**

> **HANDOFF 2026-09-15 ~08:00 (CCL at 8% context).**  Bench of record =
> 22.5 deg (resources 9b2181c: twyman_green plates at the retro end,
> tg96_run forwards D_RECOMB/D_RC_L2, both sheets 22.5/150/55/comp 200;
> clearance every part >= +38 mm; gates green).  OOM 07:25: two
> model-1024 MATLABs (my gate22_193 + another lane's tg96 clear22) ->
> systemd-oomd killed VS Code + the gate run; tg96_batch.sh now waits +
> locks (fac4771); gate22_193 RE-QUEUED detached (zwfs_dm96/runs/
> gate22_193: read its S/V/P rows vs pdi193fbase; fold into README +
> deck).  UNCOMMITTED IN THE TREE, NOT MINE: tg96_run.m (+167: Stage A
> node-part solve + a 'clearance' stage), tg96_params.m (P.clear.MOUNT /
> node), dmg_bench_clearance.m (+63), edited 07:15-07:22 by a live
> session on this box AFTER and ON TOP OF 9b2181c (runs node22s/t,
> nodesolve, clear22) -- that session likely died in the crash; its
> owner (TO restarted? CCMac?) must commit them; do not overwrite.
> Dave: CCMac will NOT return; TO takes the realism brief's items 3-8 AFTER the reflective brief (item 2 done by CCL, resources 1338920).  NEXT for TO: `BRIEF_to_reflective.md` (the OAP front end from scratch:
> diagnose the OAP-body drawing defect first, then design / layout /
> IFO + sensors on it / fold lever / P-SRI clearance).  CCMac (budget
> out): `BRIEF_ccmac_bench_realism.md` items 2-8 (status table inside).
> Deck: 44 slides DRAFT, local macos commits cd97214 and later; push on
> Dave's word.  Memory: feedback_sequenced_batch_jobs (OOM rule),
> reference_gauge_flow_* (the path), project_gauge_deck_migration.**

> **2026-09-15 later (CCL): deck 44 slides -- three jobs (measure /
> capture / hold) on slide 2; parts slides (common, interferometer,
> sensors); the analysis-path slide + backup table; station-by-station
> walk-throughs for S / V / P from the runner's new `stations_fig_`
> (runs/stations193 at 22.5 deg; raw-map residuals on the 30 nm surface:
> S 4.75 nm = the self-reference moved, V 25 pm, P 0.18 nm) -- the
> interferometer's is CCMac's brief item 8.  Memory: the repeatable path
> (reference_gauge_flow_*) and the migration plan
> (project_gauge_deck_migration).  Local: resources + macos commits since
> the 736836a / 521cf9c push.**

> **2026-09-15 (CCL): splitter RULED 22.5 deg; the phantom "glass" on the
> node figure was view_rx joining the polarizer and the tilted splitter
> into one lens (fixed: faces join only when parallel and close; tBench
> 9/9, tPropLayout 3/3); the model's parts are the scaled 56 mm rig's
> (2.6 mm plates, 5/8 mm lens centers, zero-thickness plates) -> the
> realism round; the analyzer = a polarization camera (a rotating stage
> would make the snapshot sequential); two model-1024 jobs on 64 GB OK
> (measured 11.9 GB peak).  PUSHED both repos (resources 736836a, macos
> 521cf9c) at Dave's word for CCMac.  Deck 38 slides: the three-option
> splitter slide, the node at 22.5, "How every configuration was
> analyzed" (ten steps) + the backup table per configuration.  The
> repeatable path saved as reference memory (reference_gauge_flow_common
> / _ifo / _zwfs / _pdi); MIGRATION PLANNED at the cycle's end: the deck
> + tools from macos/demo_session to templates/40_benches/gauge_deck
> (plan section 12; memory project_gauge_deck_migration).**

> **BENCH NOT BUILDABLE (Dave 2026-09-15, deck review): the record's 7-deg
> splitter puts 8 of 9 node parts in another beam (L1, input polarizer,
> compensator, output QWP, analyzer, L2, both arm QWPs' node records);
> the Stage-A solve in tg96_run cleared only the three END bodies.**  Tool
> `dm_gauge_lib/dmg_bench_clearance.m` (parts grouped by name stem, beams
> tested where they cross a part's plane, the 103 mm DM aperture as the
> beam, +8 mm mounts; draws the bench from above + the node panel:
> `zwfs_dm96/bench_bs7/22/30.png`).  Scan with the output QWP + analyzer
> moved to 160 / 170 mm behind the splitter (`D_RECOMB` 150, `D_RC_L2` 55;
> L2 stays at 207 so the tuned tail is untouched): 15 deg still collides;
> 22.5 deg clears (compensator +16 mm); 30 deg clears every part by >= 57
> mm.  The sensors do not care (bs30_dev: rows 0.9943/5 vs 0.9942/5, G1/G3
> pass, arm channel phase 1.69 vs 1.63 mrad; the plate's diattenuation
> 10% but uniform).  Dave also wants: SUBSTRATES of real thickness on the
> polarizers / QWPs / analyzer / masks (mask plate in the F/4.2 beam: 0.03
> wave W040 + 0.7 mm focus shift for 2 mm; vector QWP 0.05 wave), and a
> REAL CAMERA pitch (9.4 mm pupil image: 6.5 um sCMOS = 1450 px across,
> bin 4 to the modeled 385; shrinking the image to 385 px of a 3-5 um
> camera needs an F/0.7 field lens -- keep the image and bin; small
> pixels help the well-depth time).  Brief `BRIEF_ccmac_bench_realism.md`
> (Stage-A extended to the node parts; forward D_RECOMB/D_RC_L2; arm QWPs
> at the retro for both passes; OAP folds re-solved; substrate option;
> camera pitch; re-runs: layouts, the snapshot form's polarization at the
> new angle, CCL's V3/V4 on the new bench + the mask-sensor gate; TO's
> P/SRI node through the tool).  DAVE RULES 22.5 vs 30 (CCL: 30).  Deck
> revised (35 slides): front-end + IFO layout slides show the proposed
> 30-deg bench, parts tables carry the substrates and the camera, a
> backup slide carries the clearance scan; DRAFT status line says the
> geometry is under revision.**

> **DECK ASSEMBLED 2026-09-14 (CCL): `demo_session/deck_gauges.md` ->
> `deck_gauges.pptx`, 34 slides, DRAFT pending Dave's sign-off.**  Main
> body 24 (title, jobs + scores, front end, IFO + layout, Zernike, vector
> + layout, pinhole + layout, side by side, capture range x2, photons,
> servo, capture (descent), lenses vs OAPs, systematics, recommendation,
> modes, complex amplitude, one bench, future work x3, run it yourself),
> backup 10 (phase-shift forms, OAP rig, lost Zernike readings + fold,
> PDI trades, drift, vector polarization terms, model + sampling + color,
> provenance).  Raw material with every run tag and figure size:
> `deck_gauges_material.md` (subagent extraction from the three reports;
> four source disagreements listed there, resolved in the deck by the
> lane reports' record values: IFO rows from lens_deck on the 30 nm
> surface, IFO aging range 322 / 480+, S's N(1 pm) by tag, 52 vs 47
> sites stated).  Figures: the tools' own PNGs; panel crops by
> `crop_panels.py` (lens/oap/psri/pdi/zwfs layouts); the flow diagram and
> the one-bench schematic drawn by MATLAB; the OAP layout figure NOT used
> in the main body (its OAP bodies sit off the drawn rays -- CCMac to
> check) -- the lenses-vs-OAPs slide is a table.  Render QA: soffice ->
> pdftoppm, no overflow.  Dave (evening): Luis email sent; CCMac's
> budget arc is a LATER add; deck first.  V4 + V5 records landed and
> committed (resources 983148a).  NEXT: Dave's review of the deck;
> TO's engine-drawn universal bench (follow-up); CCMac's WS3 gate fix.**

> **2026-09-14 LATE (CCL).**  Fast suite with CCMac's SUITE_FAST change:
> 481 / 0 on Linux.  TO's pin10_2048 re-run FINISHED 18:08 (exit 137 = the
> model-2048 exit segfault; report + .mat complete: P / PF capture 120 nm+,
> N(1 pm) 2.1e14 / 3.3e14) -> TO's lane is complete.  V4 RECORD RUNNING
> (v4watch.log): an193_ref G4 0.053 pm, rows single 0.9935/4, grid
> 0.9992/3; an193_cube 3.55 pm, rows 0.9934/4, 0.9992/3 (the cube is a
> non-term at record resolution too); q300 / q100 / az1 / bound follow
> (~8 min each).  V5 DONE at dev res (README V5, plan 11.2): the pair
> ALONE cannot give amplitude + phase (mirror ambiguity about the
> reference wave, square-root conditioning; ZW.solveVA kept for the
> record); the pair + the STATE's CLEAR FRAME does (mask.v_clear; gate G9
> mask.v_dip: 0.35 / 0.50 pm through 5% / 20% dips vs 241 / 967 with the
> flat's amplitude; rows unchanged) -- a second detached watcher
> (scratchpad/v5_watch.sh, log runs/v5watch.log) applies
> runs/v5_patch_clear.py after v4seq and runs runs/v5seq.sh (an193_clear).
> Modes flow diagram DONE (demo_session/gauge_modes_flow.m -> figs/
> gauge_modes_flow.png).  NEXT SESSION: read v4watch.log + v5watch.log;
> fold an193_* into README V4/V5 + plan 7.3/11.2 (replace the dev
> numbers); commit the V5 runner patch (v5_patch_clear.py applied by the
> watcher) with the run dirs; then items (4) QA of TO's universal-bench
> drawing and (5) the deck assembly.  Local commits awaiting Dave's push
> word: macos 4adf4e7 78af1ad + this; resources e0ec760 a9309bf + TO's 8
> (fe387ee..640db68).**

> **2026-09-14 EVENING (CCL, post-compaction).**  (a) CCMac's Luis
> follow-ups (macos c44d179, resources 4772c83) pulled and gated on Linux:
> tDwDxGroups 15/15, tZernikeGridBasis 5/5, tLinkSave 3/3, both tRunCompare
> zern_grid gates pass; WS1 object-space skip ACCEPTED; WS3 gate passes but
> is BLIND to a Noll/Born&Wolf swap (index 8 lands on ANSI slot 9 in both
> tables; the crossed pairings at mode 7 give corr -0.1, the right ones
> 1.039/0.969 = grid discretization at ng 128: at ng 256 every mode is 1% /
> 1.0000) -> fix = ng 256, modes 4/7/8, 1%/0.999, crossed negative control
> (`BRIEF_ccmac_luis_review.md` round 2); fast suite re-run in flight
> (scratchpad/fast_suite.log).  (b) CCMac's as-built question answered in
> `BRIEF_ccmac_asbuilt.md`: budget FIRST (static covariance roll-up, named
> allocation table), then the dynamic predictor consumes it; covariance
> primary, MC the check; one camera model shared with the gauge lane;
> contrast a follow-on; DAVE RULES which system (telescope sens-tools =
> CCMac's framing vs the gauge bench = plan 10.3) -- CCL assumed the
> telescope tools.  (c) CCL item (1) DONE at dev resolution: the vector
> sensor's ANALYZER leak (`dm_gauge_lib/dmg_analyzer_maps`, README V4,
> plan 7.3): cube extinction a non-term (3.6 pm on 100 nm pokes), plate
> errors coherent delta/2 (lambda/300 = 162 pm uncalibrated, 1.4% of the
> figure), absorbed by the on-surface matrix (gain within 0.3%); the
> runner patch (`zwfs_dm96/runs/v4_patch_analyzer.py`: mask.v_analyzer /
> v_qwp_err / v_qwp_az, frameV_ mixing, G4 priced) is applied by a
> DETACHED WATCHER (scratchpad/v4_watch.sh, log runs/v4watch.log) the
> moment TO's pin10_2048 MATLAB (PID 1728580, ~21:00) exits, then
> `runs/v4seq.sh` (an193_ref/cube/q300/q100/az1/bound, bench + battery
> rows) queues through zwfs_batch.sh.  NEXT SESSION: check v4watch.log,
> fold an193_* into README/plan 7.3 (replace the dev numbers), commit the
> patched runner; then items (2)-(5).  TO has 8 LOCAL commits on
> resources (fe387ee..640db68: reference-arm walk, README, report head,
> pin10 requeue) awaiting Dave's push word.**

> **HANDOFF 2026-09-14 15:30 (CCL at 8% context; compact next).**  State:
> all lanes' work is on origin (both repos level).  The gauge deck plan
> `BRIEF_gauge_deck.md` (sections 1-11) is the spec; slides = main body
> each approach once in its best configuration.  Lane status: CCMac IFO
> lanes COMPLETE (report REPORT_gauge_ifo.md; figures lens/oap_vlayout
> redone; OAP focus solved, trim 6.14: dimple ZWFS survives, vector
> degrades, pinhole breaks); TO PDI lanes complete except the pinhole
> trade (pin10_2048 must be RE-RUN: CCL's pkill killed it at 14:54;
> pin10_loop ran on); Luis sens-tools WS1-4 reviewed
> (`BRIEF_ccmac_luis_review.md`: WS2 accepted, WS1 accepted + object-space
> stop follow-up, WS3 needs the engine-equivalence gate / BornWolf
> normalization, add the 2 suites to SUITE_FAST; Linux gate 472/0).
> CCL's OWN remaining items: (1) vZWFS cube channel leakage (the coated
> diagonal in zwfs_v_camA/B.in; a crosstalk term between the two
> images), (2) solveV amplitude+phase from the pair + the amplitude-dip
> gate (plan 11.2), (3) the modes flow diagram (11.1), (4) QA of TO's
> universal-bench drawing (11.3), (5) assemble deck_gauges.md (~26
> slides) from the three reports via make_brief_slides.py, render QA,
> DRAFT for Dave.  Memory: project_tg96_gauge (capture, OAP, V4),
> feedback_capture_range, feedback_shell_self_kill (kill by PID only).**

> **CAPTURE RESULTS 2026-09-14:** self-referenced readings (S, V, P)
> capture to ~60 nm of surface -- their reference COLLAPSES (blind, not
> folded; unwrap changes nothing); PF captures from 100 nm only with
> unwrap + recalibration (0.245 pm; 3 pm at cycle 26); the IFO captures
> from 150 nm with unwrap alone (CCMac; 200/300 running).  Plan section
> 7.1; recommendation = hybrid / P/SRI / second color.  TO's 30 commits
> pushed by CCL (resources 0dc8767).  CCMac asked which next: zwfs-on-OAP
> (slide 15) before the OAP descent (backup).

> **CAPTURE = A WRAP PROBLEM (TO 2026-09-13, accepted): from a 100 nm rms
> start the differential to the set point is ~2 rad rms, wrapped for every
> reading (IFO included).  Fix = 2-D unwrapping of the differential before
> the estimator (`dmg_unwrap`, `battery.unwrap`), then the start-rms
> ladder both ways -> the deck's capture slide.  Briefs:
> `BRIEF_to_capture.md` (TO: push OK; deliverables 8-9), `BRIEF_ccmac_capture.md`
> (CCMac: accepted lens_deck; "re-measured collapses" needs its mechanism
> -- wrapped differentials everywhere?; mirror the unwrapper; the IFO's
> 47-site aging range is 45 nm, its single-site 120-158).  TO's PF result:
> reference moving with the state = 5.9 pm absolute on 13 nm, frozen 0.**

> **PUSHED 2026-09-13 (Dave: "Time to push!"): resources dev-candidate
> c24245f (27 commits: CCL V2/V3/capture range/layout, TO's PDI lane,
> the merge of CCMac's renders), macos dev-candidate b219135 (19
> commits: briefs, decks, slice, engine note).  Fast suite: no failures
> in a 25-min partial pass and the first 10 suites of the full pass; the
> full pass was still running at the push (Dave's call).**

> **GAUGE DECK PLAN RULED 2026-09-13 (Dave):** title "DM Surface Gauge
> Comparison", audience the JPL HWO WFS&C group; OAPs = the budget
> baseline, show why lenses (if true) and run the other gauges on the
> OAP front end; all three IFO phase-shift forms; TO's points ruled
> (paper's pickoff form; pinhole diameter = the better of 2.0@1024/193
> and 1.0@2048/385; own dir pdi_dm96/); the draft stance approved plus
> CAPTURING THE INITIAL FIGURE (100-200 nm WFE): the descent run
> (loop.start_rms / loop.recal_every) for every reading; CCMac's renders
> not deck quality -> redo in the zwfs_vlayout recipe; CCMac and TO
> tasked at the OPUS level (`BRIEF_ccmac_gauge_deck.md`,
> `BRIEF_to_gauge_deck.md`, standalone); style: succinct, more slides,
> one figure or table + <= 3 bullets per slide; MAIN BODY = each approach
> once in its best configuration, the rest backup.  CCL owns: the vZWFS
> cube leakage, figure QA, combined figures, the deck (deck_gauges.md,
> ~26 slides).  NEXT: CCL's cube-leakage run; then wait for the lanes;
> assemble.  The two briefs are the handoff text for Dave to relay.

> **GAUGE DECK PLAN 2026-09-13 (Dave: "bring all the DM surface gauge
> concepts together into one new deck ... Let's plan!"): `BRIEF_gauge_deck.md`
> -- six configurations on one front end (lens IFO, OAP IFO, ZWFS, vZWFS,
> PDI pinhole, P/SRI), the common-currency table with each cell's run tag
> or owner, the layout/parts rule (the zwfs_vlayout recipe), an 18-slide
> outline, the work split (CCMac: 30 nm rows + capture range + photons +
> layouts + parts + close the reflective design; TO: PF through the two
> decks + capture range + reference-arm drift + parts; CCL: cube leakage
> + assembly), five decisions for Dave.  Awaiting his ruling; nothing
> built.  CCMac's render commits (5e45851, b6ed9fe) merged locally
> (c24245f); their OAP table-plane panel is edge-on, labels E-numbers.**

> **DECK QUESTIONS ANSWERED 2026-09-12 (Dave: capture range to 10%
> recorded for every device -- "they will not be operating at null";
> the vZWFS layout drawn large; slide 2's mask figure swapped).  LOCAL
> (resources aa3bbb7 + this commit; macos this commit).**  Capture
> range: runner prints it after every ladder; aging calibration (matrix
> at 30 nm, 385 rays) L 44 / I+ 36 / S 42 / V 70 nm; re-measured on the
> surface, gain within 5% to 160 nm for all, photons for 1 pm x5 (L) to
> x40 (V); IFO 120-158 nm (its four-step wrap).  Layout: zwfs_vlayout.m
> (QWP + 12.7 mm MacNeille cube, two cameras; decks zwfs_v_camA/B.in;
> the cube's diagonal leakage = next item).  Mask figure: runs/mask385
> (2.0 lam F/D spot, 3.96 px per lam F/D at both record settings; the
> old figure was stage 1's spot 9).  Deck: 28 slides (11-12 capture
> range, 17 layout, 18 V2+V3), DRAFT.  Model-2048 MATLAB segfaults at
> EXIT after writing everything (mask385, exit 137) -- results intact.

> **ZWFS V3 DONE 2026-09-12 (LOCAL: resources 65788d9 + the record
> commit; macos this commit; push on Dave's review): the arm's
> polarization aberration per circular channel.**  `dmg_arm_maps` (two
> polarized vector traces -> J per grid pixel at the sandwich's entrance
> sphere, the mask's PROJECTED-axes basis per ray, J/sqrt(det J) --
> vector-mode field = scalar's exact conjugate, slope -2.0000 --, unit
> mean power per laser state); gauge V_ARM none|engine|synthetic,
> chained apodization (G8), cross-channel leak, per-pixel per-channel
> solver, v_cal ideal|amp|fit|map; runner knobs + prints.  Record (193):
> lens rig channel PHASE difference 1.6 mrad rms at laser 45 deg (nil on
> the fold plane's s axis), uncalibrated 9 pm on 100 nm pokes, oracle
> 0.053, on-surface rows = ideal to the digit; scan: diattenuation-type
> 6 nm/rad absolute, nothing through the matrix to 0.1 rad, 0.3 = grid
> row 2.6x + loop gain -10% (no fixed error); retardance-type 3.8x larger
> absolute but REMOVED by the per-channel unmasked frames ('amp').
> Record: README V3, brief V3, deck slide 18 (28 slides, DRAFT), memory.
> Queue: TO's pseq6 waits on my "v3seq2 done" (written 16:01).  NEXT: V4
> the stepped reading's between-frame drift; integral term (thermal);
> fold-aware I+; deck sign-off.**

> **ZWFS V2 DONE 2026-09-12 (LOCAL: resources this commit; macos this
> commit; push on Dave's review): metasurface retardance error = one
> complex constant kappa; uncalibrated absolute bias 1.3x the leaked
> amplitude in phase / quadratic in quadrature; 3-number flat
> calibration removes it (0.048 pm); differential rows + loop through
> the on-surface matrix IDENTICAL to the ideal mask at 0.1-0.2 rad ->
> not a servo-budget term.  Knobs mask.v_ret_err / v_leak_phase / v_cal.
> NEXT: V3 arm polarization aberrations per channel (Jones pupil); V4
> the stepped reading's between-frame drift in the loop.**
>
> **PUSHED 2026-09-12 (Dave: "Push yes. Merge yes. Brief CCMac!"):**
> tg96-oap MERGED into MACOS_resources dev-candidate (0fd6786; fast
> suite 463/0 incl. tBench) and pushed; macos dev-candidate pushed.
> CCMac's next brief: `BRIEF_ccmac_tg96_loop.md` (deliverable 7 on
> dev-candidate + the bare-Al / Jones-pupil runs).  Worktree
> `~/dev/MACOS_resources_wt_oap` can be removed.
>
> **ZWFS V1 -- THE VECTOR (POLARIZED-DIMPLE) READING, DONE 2026-09-12
> (LOCAL: resources 1cf5889 + e3e32e0; macos this commit; push on Dave's
> review).**  Reading V = +phi/-phi image pair, exact per-pixel solve, no
> fold (`dmg_zwfs_gauge` frameV/reconV/solveV/diffV); runner class 4 in
> every stage; G4 fold gate (100 nm sparse pokes: V 0.05 pm, single frame
> 9 nm).  Record (README V1; runs/v193flat, v193base, v193noise,
> vloop193): matrix on the 30 nm surface -> single 10 nm 0.9935 / 4 pm,
> grid 1 nm 0.9992 / 3 pm, dense 0.9999 / 0.33 nm (S 0.68); ladder gain
> within 2% to 60 nm rms (S 0.66); loop contraction 0.509, steps ->
> 0.000, 3 pm from 1.5e12 / 5.3e12 photons per cycle (S 2.6e12 /
> 7.5e12); N(1 pm) 4.7e13 (S 5.6e13).  Learned: a whole-pupil 60 nm
> figure collapses the core for every reading (not a fold test); b
> iteration exact modulo piston; V's differential must be the WRAPPED
> phase difference; matrix mode reports RAW (S10's ladders were
> Wiener-corrected); the noise stage prices the readings run; an eager
> ifelse_ killed any battery without I+ (fixed).  Deck slide 15 (DRAFT).
> NEXT: V2 metasurface retardance/leakage knob; V3 arm polarization
> aberrations per channel (Jones pupil); V4 the stepped reading's
> between-frame drift in the loop.**
>
> **CCMac tg96-oap f46585d REVIEWED 2026-09-11 (`BRIEF_ccmac_tg96_oap3.md`,
> LOCAL): accepted in substance; D2 lens numbers reproduced here to the
> digit (worktree of the branch + this box's mex; appendix in the brief).
> Six items before the merge: evidence dirs not committed; tuned tail
> keyed by TAG (any other tag = untuned bench, null 9.11 vs 0.134 nm);
> the OAP dense loss (0.75 / 4.8 nm, x-talk 0.42) not separated from
> window truncation / regularization (step 16 or Voronoi cells; lambda
> sweep); D1 "worst at centre" wants a map; non-vacuity stated on the
> lens (passes there, meaningful on the OAP); calib_surface 'base' +
> a wrap guard.  Then deliverable 7 (the loop).  Worktree
> `~/dev/MACOS_resources_wt_oap` holds the runs (remove when done).**
>
> **S11 LOOP METRIC BUILT + RUN (2026-09-11, LOCAL: resources 88fb9ac,
> 922c366 + the record commit; macos this commit; push only on Dave's
> review).**  `dm_gauge_lib/dmg_loop.m` (the ONE loop for both gauges) +
> `tests/tDmgLoop.m` (9 gates, SUITE_FAST) + `zwfs_run` stage 'loop' /
> `P.loop`.  Record `runs/loop193` (98 min, 2562 states): the loop
> propagates noise and a walk per theory; L and S cost the SAME light
> (3 pm from 2.1e12 / 2.6e12 photons per cycle noise-only, 7.3e12 /
> 7.5e12 under a 2 pm walk); S has NO fixed error (steps -> 0.000 pm,
> thermal = the lag 10 pm); L imprints high-frequency error under a
> persistent low-order residual (thermal 27.6 pm, creeping); I+
> DIVERGES (fold-flipped sites = negative gain).  Record: README S11,
> BRIEF_zwfs_campaign S11, deck_zwfs slide 14 (DRAFT), memory
> project_tg96_gauge.  IFO half handed to CCMac: BRIEF_ccmac_tg96_oap2
> addendum (deliverable 7).  `runs/loop385` (385-ray check, L + S) reproduces
> it (3 pm at 2.1e12 / 2.6e12 noise-only, 7.2e12 / 7.7e12 walk; thermal
> floors 24.7 / 10.1 pm).  Fast suite 462/0 with tDmgLoop.  OPEN: integral term for the thermal case;
> fold-aware I+ (re-take the branch prior); the IFO row.
>
> **ZWFS S10 LANDED (2026-09-10, LOCAL: resources 8181893, macos this
> commit; push only on Dave's review).**  Dave's measured response
> matrix (sparse multiplexed grids, dw/da) is the runner's default
> calibration (`battery.calib_mode` 'matrix'; the ZWFS piston null
> carried as a rank-one term -- without it the fit over-responds 2-4x
> at low order); flat single actuator 0.994 / 4 pm; matrix ON the
> working surface -> stepped reading 0.989 / 5 pm / SNR 2200 on a 30 nm
> surface, linear one-frame 1.05 / 23 pm; exact one-frame readings
> limited by the fold sensitivity.  Alternating +/- pokes: neutral in
> the model, kept for bench drift.  Decks re-stated (deck_zwfs 22
> slides, fold 3d: the matrix-calibration slide; deck_tg_fang head-to-
> head), both DRAFT.  Record: README S9/S10, BRIEF S9/S10, memory
> project_tg96_gauge.  OPEN: deck sign-off; the matrix for the IFO
> (CCMac's tg96 runner, brief addendum); vector ZWFS; the 2048-grid
> matrix run (not needed: 385/1024 == 385/2048 shown thrice).

> **dwd* PLOT SIZE -- DONE 2026-09-10, LOCAL, awaiting Dave's review
> (`BRIEF_dwd_plot_size.md`).**  Size is fixed FIRST and the pages follow:
> `mmacos/sensitivities/dw_page_layout.m` + `dw_page_fig` / `dw_page_axes`
> / `dw_draw_map` / `dw_block_keys` / `dw_canvas_tiles` / `plot_dw_index`
> / `write_page_index`; `plot_dw_channels` paginates on ELEMENT
> boundaries and returns a manifest; `plot_dw_per_element` gains `field`
> mode; `plot_opd_canvas` sized by tile count; `run_sensitivities` takes
> `panel_in` / `tile_in` / `page_in` / `page_max_in` / `max_per_page` and
> writes `<name>_pages_index.txt`.  Floors: 3.5 in per OPD map, 1.2 in per
> FIELD TILE of a canvas; a sparse page GROWS to fill 16:9, a page grows
> to [32 20] in to keep one element whole, only then splits (`_p02`).
> **Measured on the zoom fixture:** one dwdsurf channel's 567x567 canvas
> was drawn in a 14x15 px box with 40 non-white pixels (a ~9:1 subsample);
> now ~1500 px across, ~167 px per field point, 21 pages instead of one
> sheet.  Per-element centre page: 635 -> 999 px of drawn map on FEWER
> total page pixels (recovered subplot margin).  Two bugs closed en route
> (per-field cells are Nc x Nf -- index them 2-D; row-per-block slots must
> use the PAGE's column count).  10 gates in `tRunSensitivities`.
> Acceptance: the zoom drivers re-run (dwdsurf 200 s / 63 pages, dwdx
> 330 s / 92 pages) and the pages READ.  `<name>_<ch>_channels.png` KEEPS
> its name -- the single page when one suffices, else the page-labelled
> INDEX contact sheet, full-size pages in the (gitignored) `_pages/` --
> so no README or artifact reference breaks.  Deviation to flag:
> byte-identity with the old small-deck pages is NOT achievable together
> with the size floor (the old pages spend ~30% of every cell on subplot
> margin), so the gate asserts one page + the historical filename +
> panels no smaller than before.  Also regenerated: the zoom dwdsurf artifacts
> were STALE (4 channels; the powered set became 42 on 2026-09-05) --
> README corrected.  Bindings held: jet, zeros blank, no "leak" language,
> push only on Dave's review.  ZWFS S7 + Luis round 4 CLOSED in memory
> (`project_tg96_gauge`, `project_opd_conventions`); macos `ac8bf3b` LOCAL
> on top of pushed `6a0dc31`; resources was pushed `8bc8fcc`.**

> **LATE 2026-09-09 (Luis round 4, after S7): FIXED and PUSHED (Dave:
> "good results -- go ahead and push"; macos 097dd52, resources fe6c6bd;
> deck REBUILT on jwst_ote_designc: `demo_session/deck_dwdsurf_options.pptx`
> (macos 353e942 pushed); second runner bug fixed en route: 'elts' never
> reached the dwdsurf channel, resources 30b6627 pushed; THIRD bug =
> what Luis's pictures actually showed: the per-element centre-field
> page scrambled under orient xy, resources 1e0a956 pushed; DECK = the
> runner's own pages ONLY, unmodified -- Dave's standing rule, memory
> feedback_deck_plots_unmodified) --
> `run_sensitivities` never forwarded 'orient'/'sign'; no dw_d* driver
> set the OPD reference, so every harvest was mean-referenced and a
> single-segment poke pistoned every other segment by
> -(N_k/N)*mean(poked) (driver-path measurement e5hex1 seg 2 Kr/Kc:
> 14.7% / 12.0% of the poked rms; 'chief' -> exactly 0).  'opd_ref'
> {'mean','chief'} now in all 8 drivers + core (re-applied after every
> reload) + runner (+ 'orient'/'sign'); defaults 'mean'.  Record: PLAN
> 0.x dated block, `mmacos/doc/SENSITIVITY_TOOLS.md`, gate
> `tOpdRef/test_driver_single_segment_poke_is_local_under_chief`, reply
> `DRAFT_email_luis_round4.md`.  OPEN for Dave: the chief ray's OWN
> segment stays non-local under 'chief' -> nominal fixed-length
> reference (`opd_ref_len_set`, engine OPDRefRayLen branch) + the
> default flip.  Dave also asked for a few slides of the dwdsurf results
> under all options (figs: demo_session/figs/dwdsurf_*.png; deck md
> `demo_session/deck_dwdsurf_options.md`).**
>
> **CURRENT STATE (2026-09-09 evening).  S7 of the ZWFS campaign RAN:
> the iterated-reference exact reading landed AND a MODEL DEFECT was
> found and fixed on the way** -- the `twyman_green` 'nf' mask sandwich
> was ASYMMETRIC (exit sphere zElt/Kr 23.86 vs entrance 352.7 mm), so
> the engine's SPH2PL quadratic factor Fresnel-DEFOCUSED the reimaged
> pupil by an effective 4.86 m for EVERY ZWFS number in S1-S6 (ringed
> kernel, Talbot transfer null at ~30 cyc/ap, 29% amplitude modulation
> under the 30 nm state, the 744 pm base floors, the "undetected"
> grid-on-base).  Fixed: 'nf' symmetric (round trip 1.8e-15),
> 'nf_legacy' byte-identical (tBench gate); S1-S6 = legacy record.
> Record: `MACOS_resources/mmacos/templates/40_benches/zwfs_dm96/
> README.md` (banner + S7 bullet), `BRIEF_zwfs_campaign.md` S7 section,
> memory `project_tg96_gauge`.  Headline: on the corrected model the
> linear reading alone takes the single-on-base floor 744 -> 67 pm and
> grid-on-base SNR 1.46 -> 14; the one-frame iterated reading with the
> REFINED base prior (I+; `dmg_zwfs_gauge` priorS) reaches 13 pm / SNR
> 584 (96x96) and HOLDS to 60 nm rms working state (gain 0.85-0.91)
> where the four-frame stepped reading falls to 0.56; without the
> refined prior the iterated readings cliff between 30 and 40 nm (branch
> fold).  **NEXT (Dave's steer): (1) tell Dave -- the S1-S6 ZWFS record
> + deck_zwfs slides are on the defocused model; (2) NGRID 385 compliant
> run (the 0.90 raw hold-out at 96x96 is sampling); (3) S6 color RE-RUN
> on the corrected model; (4) S5 noise pricing of I+; (5) deck fold.**
> IFO untouched (all-geometric).  Everything below is the prior state.
>
> **PRIOR STATE (2026-09-09 midday).  NEXT was the ZWFS MODELING arc:
> the iterated-reference-wave reconstructor** -- read
> `REPORT_zwfs_lit_scan.md` FIRST (the ranked literature imports + the
> PZT verdict; "Suggested first ZWFS task" at the end is the spec), then
> `MACOS_resources/mmacos/templates/40_benches/zwfs_dm96/README.md` (S1-S6
> findings; S6 = multi-color) and `dm_gauge_lib/README.md` (the ONE scoring
> library; `dmg_zwfs_gauge` gains a third reading).  Memory
> `project_tg96_gauge` has the campaign state.  Everything below in this
> block is CLOSED and PUSHED: sens-core merged (resources `57a6ec0`,
> branch deleted), MR text `MR_dev_candidate_to_dev_2026-09-09.md`
> (Dave files it), CCMac round 2 closed (`BRIEF_ccmac_jpl_round2.md`),
> round 3 handed to CCMac (`BRIEF_ccmac_jpl_round3.md`: confirm the
> phantom-grid save fix `cda178e` on IRIS + run the instrumented
> supervisor for the reset_xp canvas fault), engine fixes `cdf8636`
> (Get_Values over-read) + `cda178e` (phantom-grid SAVE crash), decks
> carry the S6 color slide (`2b4821c`).  Open engine items in PLAN
> section 0: **IRIS save_rx -> reload SIGSEGV on the REAL ZrnGrData grids
> (iElt 17/19/21/35/37/39) -- STILL OPEN** (CCMac round 3; the phantom
> guard `cda178e` was only half; chase with a debug build on the IRIS
> deck, breadcrumbs in `REPORT_iris_save_crash.md`; off the merge path);
> the **dw_multi_core reset_xp=true per-field aggregation collapse to
> [12 12]/0** (CCMac localized it -- direct dw_dx is healthy [256 256];
> sens-core follow-up + fix my NaN-rays mis-route in the emptyOPD
> warning; needs the IRIS deck); lensarr TRACE-time overrun
> (tst_save_keys only); trace(26)->trace(27) stale OPD.  CCMac round-4
> note: `BRIEF_ccmac_jpl_round4.md`.  Tips: macos dev-candidate
> `6174b1f`+, resources `e36db7c`.
>
> **PRIOR STATE (2026-09-08 late) kept for the record.  NEXT was the Fang
> thread: IFO and ZWFS DM-gauge decks** -- `deck_tg_fang` (22 slides, macos `82ced2b` fold) and
> `deck_zwfs` (DRAFT, 9 slides); record in memory `project_tg96_gauge`
> (complete thru S5; deck folds DONE; next stages on Dave's steer) +
> `BRIEF_zwfs_campaign.md` (uncommitted edits present) +
> `templates/40_benches/{tg_psi_dm96,zwfs_dm96}` READMEs.  Read those before
> touching decks; STYLE_REPORTS section 5 gate before any rebuild.
>
> **Closed today, all PUSHED (SHAs in memory `project_fex_ep_radius_rework`,
> `project_ep_dome_review`, `project_zoom5x5_pupilfind`,
> `project_luis_sens_hardening`):** FEX Rx-order guard removed; STOP
> multi-value prompts + Segment stops; FEX four-probe medial pupil (five pins
> re-pinned); element-STOP handedness fix; EP-dome review -> pupil-read ruling
> (sens-core preflight + CCMac's single-DOF error variant); re-traces made
> idempotent (OrthoSrcFrame dead band -> dwdsurf speckle floor exactly 0).
> macos dev-candidate `48451d8`; resources dev-candidate `bf1b182`, sens-core
> `f02dc60` (merges after CCMac's JPL pass,
> `BRIEF_ccmac_jpl_private_verification.md`).  Open engine item: PLAN section 0
> `trace(26)->trace(27)` stale first OPD.  Luis draft ready:
> `DRAFT_email_luis_sens_noise.md`.
>
> **LATE ADD (2026-09-08 night, Dave's parting ask, run autonomously):
> S6 COLOR** -- both DM gauges at 480/532/632.8/700/780 nm + a multi-
> channel Wiener combination (`dm_gauge_lib/dmg_color_comb`).  ZWFS: the
> transfer nulls migrate as 1/lambda and the 5-color combination fills
> them (base-row floors 3-4x, dense random 2.3x, grid-on-base SNR 1.46 ->
> 3.43, still < 5).  IFO: identical at every color -> its roll-off is
> GEOMETRIC (tail conjugate + distortion), the joint tail objective is
> its lever.  Resources dev-candidate `b47dd5f`, macos `42cffde`
> (BRIEF_zwfs_campaign section); memory `project_tg96_gauge`; figure
> `zwfs_dm96/zwfs_s6color.png`.

> **PRIOR STATE (2026-08-19) kept below for the rodgers3 record.**

> **CURRENT STATE (2026-08-19 eve).  NEXT STEP = execute
> `BRIEF_rodgers3_s1.md`** (the offset_imager template + the
> challenges/rodgers3 instance; Stage 0 is DONE — all 5 rungs of
> Mike's ladder gate, artifacts in `~/dev/MACOS_sandbox/Design/
> Rodgers3/s0/`, conventions in memory `[[project_rodgers3]]` +
> `[[reference_codev_zrn_convention]]`).
>
> Board as of tonight (each item's record in the named file/memory):
> - **dev-candidate = the consolidated, validated pair, PUSHED**:
>   macos `136f353`, resources `052ef30` (includes #67 pol arc,
>   GridFile fixes, rodgers1 regen, pol-ifo + pol-core(res) merges,
>   the examples reorg templates/+challenges/ tree, 594/0 suite,
>   pymacos 6694/0 + PROPER 26/26).  `[[project_branch_model]]`,
>   `[[project_segmirmaker_audit]]`.
> - **Andy**: fast-forward both `dev` to dev-candidate tips, then
>   handoff steps 3–4 with fresh dry-runs (MERGE_HANDOFF_ANDY.md has
>   the status note).  pol-core/pol-ifo/develop deletable
>   (archive-tag develop first).
> - **`BRIEF_luis_round2.md` EXECUTED (2026-08-19 night), LOCAL,
>   unpushed.**  Report: `~/dev/MACOS_sandbox/notes_luis_08-19-2026/`
>   (`REPLY_luis_round2.md` = the three drafts for Dave to send).
>   Headline: **the brief's own correction of Luis was WRONG and his
>   diagnosis was right.**  `macos_cmd_loop.inc:366` does set
>   `LUseChfRayIfOK=.TRUE.` on every load, but ~40 lines BEFORE
>   `MBFile6`, whose first statement (macosio.F AND smacosio.F) is
>   `reinitialise_variables()` → `ray_mod_init_vars` → `.FALSE.`.
>   So the flag is FALSE on every path, the chief branch was
>   unreachable, and his missing-`'Y'` patch is the fix — adopted,
>   credited.  Nor is "chief dies on segmented decks" the mechanism
>   on e5hex1: `LRayOK(1)=1` at all five dw_dx_multi fields and the
>   map is mean-referenced anyway.  Measured piston on the 6 unpoked
>   segments: `+2.849e-06` (16.7% of peak) → `0.000e+00` exactly
>   under the chief reference; the poked segment's peak recovers by
>   the same constant.  Landed: `UseChfRay4OPD= Y` parses;
>   `opd_ref_set/get` + `macos.opd_ref` + Session method; `init`
>   resets `rxLoaded` (SystemCheck was passing on a wiped model);
>   `macos.init` screens model size (the engine's `stop` kills
>   MATLAB); `macos.unload()`; Session `opd(...)` now forwards
>   'orient'/'sign'; 25 new Session delegators + a coverage GATE;
>   tOpdRef (8) + tEngineMemory (4) + 2 tReadGridFile cases.
>   Item 2 verdict: the two grid paths are BIT-IDENTICAL on
>   dev-candidate — stale-version artifact (a bare readmatrix is the
>   transpose); now permanently gated.  Item 4: no steady-state
>   leak, but `clear mmacos` LEAKS ~720 MB/cycle — do not use it.
>   Item 5: API sketch only (`mmacos/design/PLAN_CONFIGURATIONS.md`),
>   awaiting Dave's review + a j18 deck (none exists in either repo).
>   OPEN for Dave: flip the global OPD default?  (PLAN.md §0.x has
>   the compatibility sweep + the regen list; nothing regenerated.)
> - **Keysight demo ~2026-09-01** (`[[project_keysight_demo]]`):
>   Week-1 reorg DONE; Week-2 rewalk of all design examples PENDING
>   (carry-overs listed in mmacos/PLAN_EXAMPLES_REORG.md); rodgers3
>   Stage 1 is demo material.
>
> **PRIOR STATE (2026-08-02) below — e2e2 TELESCOPE FLOW COMPLETE.**
> Report: `MACOS_res_dev/mmacos/design/examples/e2e2/E2E2_REPORT.md` (for
> Dave + CCL Fable).  Brief `macos/BRIEF_e2e2_implementation.md`; plan
> `MACOS_res_dev/mmacos/design/PLAN_TMA_E2E2.md`.  **8 commits, LOCAL on
> `MACOS_res_dev` `dev`, NOT PUSHED**, on top of the 4 pupil-fix ones:
>   `f7c6fd7` hoist the strict-WFE kernel rodgers1 → design/src
>   `5daaef7` fix the best-focus rung (it lost to the rung below it)
>   `f145d4e` e2e2 S0+S1
>   `87ad634` field widened to 0.6 deg (later narrowed)
>   `06c4b5f` fold before bias; price the two knobs together
>   `5f48f68` carry both frontier branches; why freeform cannot help
>   `0442d31` 0.4 deg box finishes the telescope; relay needs its own design
>   `d54e7f5` telescope closed at S2, relay parked, final score
>
> **ENGINE (2026-08-02, after CCLF's repo cleanup + PR #70 merge):**
> `~/dev/macos` is the PRIMARY tree and is now on **`dev` `ba23d93`**.
> The throwaway `macos_dev` pairing-check worktree is removed and the
> merged `colsource-pupil-fix` branch deleted; `macos_polfix` (pol-core)
> stays for polval.  NO Makefile/alias changes -- everything durable
> (shell aliases, both binding Makefiles' MACOS_BUILD_DIR, pymacos cmake,
> the CC project config) already anchors on `~/dev/macos`, so the
> physical build dir every e2e2 number came from is unchanged and is now
> correctly labelled.  Pairing is carried by the BUILD WIRING, not by
> directory names: `~/dev/macos` (dev) <-> `MACOS_res_dev` (dev),
> `macos_polfix` (pol-core) <-> `MACOS_resources` (pol-core).
> **`MACOS_res_polifo` still needs a macos counterpart decision when
> pol-ifo revives -- that question returns after PR #67.**
> RE-VERIFIED after the move: rebuilt release gfortran, mex relinked,
> `tStrictKernel` 4/0, `tE2E2Axial` 7/0, and **`s3_report.txt` re-scored
> BYTE-IDENTICAL** -- the 36a5f5a->ba23d93 delta is comment-only, checked
> rather than assumed.
>
> **FLOW (3 stages + scoring):** `e2e2_params` → `s1_axial` → `s2_fold`
> (the delivered telescope) → `s3_score` (designs nothing).  The relay is
> PARKED in `relay_followon/`.  Two stage reorders on Dave's call: fold
> BEFORE bias (geometry is free, bias costs ^1.80), then the off-axis
> stage folded into the relay stage, then the relay parked entirely.
>
> **DELIVERED**, uniform 13x13 over the 0.4 deg box, D=3 m f/20 500 nm,
> bias 13', M1 hole 0.2331 m (2.4% of the area, measured):
>   strict-chief      49.220 nm  Strehl 0.679  FAIL
>   strict-centroid   30.001     0.867         PASS   <- primary
>   + best focus      29.380     0.872         PASS
>   + LS tip/tilt     20.380     0.936         PASS
> Bar 35.71 nm / Strehl 0.80.  DL at the primary reference and everything
> more permissive; missed at the strictest.  Coma 0.11-5.33 um; mapping
> f.theta to -0.0075%, 630 um rms departure, 522 um nonlinear.
>
> **TESTS:** fast 236/0, tE2E2Axial 7/0, tStrictKernel 4/0.
>
> **FINDINGS, all recorded in-source, none fixed silently:**
> 1. Brief's hoist gate is STALE (baseline predates PR #70); replaced by a
>    same-engine A/B, 0.000e+00.  Baseline NOT regenerated.
> 2. rung-3 fix moves rodgers1 rung-3/4 artifacts by 2-3e-4; rungs 1-2
>    (incl. the centroid ruling) unchanged.  **Regen = a reviewed step,
>    NOT DONE** — `rodgers1_pupil_audit.mat`, `rodgers1_dense_field.mat`,
>    PACKET Addendum 10's rung-3/4 columns.
> 3. The offset does NOT scale with the field (RMS ~ bias^1.80 measured).
> 4. Freeform on pupil-conjugate mirrors cannot reach a field-VARYING
>    residual (astig reverses across the field, spread/mean 4.48).
> 5. A relay is sized by the IMAGE, not the aperture.
> 6. Bias beat extraction tilt ~100x here — the REVERSE of e2e VIS.
>
> **OPEN:** (a) `s3_score.m`'s `distortion_` had THREE successive errors
> (wrong quantity, forced parity, transposed scale) — all caught, the
> third now asserted, but the section wants fresh eyes before its numbers
> are quoted externally.  (b) The relay's field-corrector hypothesis is
> UNTESTED (see `relay_followon/README.md`).  (c) Nothing pushed.
>
> **RESUME:** `E2E2_REPORT.md`, then `examples/e2e2/README.md` (design
> point, FOV sweep, the numbered solve-order rules each with the failed
> run that earned it), then `s2_report.txt` / `s3_report.txt`.
>
> **TRAPS PAID FOR:** don't edit `run_mmacos_tests.sh` while a suite runs;
> `strict_ladder_deck`/`stage_score` need `macos.init` first; `pupil_gate`
> projects along the INCOMING chief (`get_src_fov`), not `get_ray_info`'s
> outgoing `.dir`; gate AOI **spread** on **powered** surfaces only
> (`aoi_report` includes the flat fold at its 45 deg); judge clearance
> BEFORE `add_pupil`; more DOFs must never make a reported design worse
> (take min per branch).

> **PRIOR STATE — ColSource pupil fix, LANDED.**
> macos **PR #70** open (`36a5f5a`, branch `colsource-pupil-fix`):
> collimated traces now use the DECLARED `Aperture=` (was
> `1+2/(nGridpts-1)` oversize — diagonal corner lobes only, which is
> why it hid; full diagnosis rodgers1 `PACKET.md` Addendum 10).
> Verification COMPLETE across three sessions (Terminal Opus pre-crash
> A/B + this session's completion runs): the ONLY test fallout is 3
> pinned-baseline tests + the uncommitted `tPupilAperture.m` gates.
> Wrap-up work list = **`macos/BRIEF_pupil_test_wrapup.md`** (re-pins
> with measured values, gate-class commit/registration, ordering:
> resources commits land AFTER #70 merges, no push until Dave).
> **WRAP-UP DONE 2026-08-01, committed LOCALLY on `MACOS_res_dev` `dev`,
> NOT pushed** — `851e33c` re-pins (tMacosPkg 12850→12454 / 1366→970
> launched+obscured, n_ok unchanged; CassFF peak/sum +1.6% = the pupil
> AREA change, mmacos + pymacos twin), `32f39d9` tDesignTelescope
> ray_bundle Y-slice made pitch-aware (the old `|x|<0.05` window passed
> only BECAUSE the oversize corner ray shrank the normalization; real
> cause = half-integer column lattice + the lone off-lattice chief ray),
> `8d3af7c` `tPupilAperture.m` gate class + its own 512 batch line in
> `run_mmacos_tests.sh` (all probes pinned to 512 so the class never
> transitions model size in-process).  Verified post-fix: tMacosPkg
> 25/0, tDesignTelescope 70/0, tProperCompareCassFF 4/0, tPupilAperture
> 5/0, pymacos `test_cass_ff.py` 4 passed; A/B against the preserved
> pre-fix mex gives tPupilAperture 1/4 (all of section A red).
> **Remaining: merge PR #70, then push the resources commits (Dave's
> call).**  Durable A/B evidence copied out of the volatile session
> scratchpads to `~/dev/MACOS_sandbox/pupil_fix_ab_20260801/`
> (pre/post mexes, probe `.mat`s, full-suite logs).
> Also committed: `26a1151` `rodgers1/dense_field_check.m` — run to
> completion (it never was pre-crash) and it ANSWERS the open average
> residual: max unchanged on all three designs (0.991/1.003/1.000× his)
> while the S3 average walks 1.082× → 1.000× on a uniform 9×9 grid.  The
> avg gap was the .seq quincunx's edge weighting, not the optics.
> **Untouched by me:** a CONCURRENT MATLAB run (not this session — my
> `-batch` log shows only the dense-field banner) rewrote a batch of
> committed `rodgers1_seq_*` decks/PNGs/`.mat`s plus
> `rodgers1_pupil_audit.mat` between 07:46–07:58, and wrote
> `rodgers1_dense_field.mat`.  Left MODIFIED and uncommitted — someone
> else's in-flight regen, and regens are a reviewed step anyway.
> **Second fallout surface — worked-example artifacts:** every
> committed e2e/e5 regen count and metric derived from a collimated
> trace (s3 "12106/12520, 376 gap clips" class counts, WEM tables,
> probe footprints) shifts slightly post-fix (~0.1–4% where corner
> rays reach the detector; point-source decks bit-identical).  Treat
> shifted counts on regen as EXPECTED, not regressions; regens are a
> separate reviewed step (same policy as the e5_pie stale-baseline
> flag).  Durable record: memory `project_pupil_oversize.md`.
>
> **QUEUED (Dave 2026-08-01): e2e2 — the improved TMA design-flow
> worked example** (params → Korsch axial → off-axis → fold → relay+FP,
> Rodgers doctrine folded in: joint solve, stated references, solve
> set ≠ scoring set, pupil gate, parameter provenance per stage).
> Full plan: `MACOS_res_dev/mmacos/design/PLAN_TMA_E2E2.md` — written
> for cold implementation by Opus/Sonnet/users.
>
> **PRIOR STATE (2026-07-23) below.**  No other half-done slice
> is in flight — everything below is LANDED + PUSHED.
>
> **§0 MODEL-TRANSITION HEAP CRASH — RESOLVED (macos `0b07046`).**  The
> reopened §0 crash (SIGSEGV at tViewRx setup after a 128→256 transition)
> is FIXED: 4 model-sized allocate-once buffers, never regrown on
> `macos_init_all(larger)` — `CumLStart`/`srcMap`/`ds1-ds2` (`54270af`)
> then the actual tViewRx culprit **`DrawRayVec_save`/`DrawEltVec_save`/
> `nDrawElt_save`** (the DRAW buffers `view_rx` harvests; `0b07046`).
> Found via `-fcheck=bounds` on a standalone `smacos_dvr` — **ASan
> structurally could not** (overflow clears its redzone into another live
> allocation).  CCMac converged on the diagnosis + split the no-arg
> runner per-model-size-group as a workaround (`bc5e8e1`); to VALIDATE
> the engine fix you must run ONE MATLAB process spanning model sizes
> (the split sidesteps it).  Full trail: `[[project_model_transition_crash]]`.
> Also this session: even-grid center fix (`[[project_even_grid_center_fix]]`),
> iris scrub, opt_example/test_calib/vsg adds, ASan/-fcheck infra
> (`build_asan`, `scratchpad/heartbeat.sh` liveness monitor).
>
> **⚠ FRIDAY (~2026-07-24) PUBLIC-RELEASE HISTORY REWRITE.**  `nasa-jpl/
> macos` goes public with history SCRUBBED of NPSOL / pgplot / etc.;
> `main` ← sls-dev functionality; a public `dev` branch keeps the
> developer files (stripped from `main`); users RE-CLONE.  Do NOT push
> from a stale pre-rewrite clone afterward (reintroduces scrubbed
> history).  Full-history safeguard taken: `~/macos-archive-20260722/`
> (bundles `git bundle verify`'d + clone-tested).  See root CLAUDE.md
> "Public-release strategy" + `[[project_branch_model]]`.
>
> **This session (2026-07-22) — all SHIPPED + COMMITTED + PUSHED:**
> (1) Luis's 3 mmacos gaps — SPOT `'beam'` obscured-chief-ray ENGINE fix
> (`tracesub.F` LocalCoord) + 25 `elt_srf_*` surface-inspection veneers +
> ~50 Session query methods (`[[project_spot_beam_veneer_sync]]`);
> (2) OPTIIX purge → `e2e_pie` fixture, GMI regression repointed 6/6,
> zero optiix refs (`[[project_optiix_removal]]`); (3) cmdref regenerated
> (`gen_cmdref.py`) + spot NOTES; (4) ifx `build_release` rebuilt (spot
> fix for pymacos).  Tips: macos `dae015e`, MACOS_resources `455fe26`.
>
> **This session, later (2026-07-22) — EVEN-GRID CENTER FIX SHIPPED**
> (`[[project_even_grid_center_fix]]`): captured upstream `ff85575`
> (half-pixel figure center for EVEN `nGridMat`) — bug lived at 4
> surfsub.F sites (SGSrf + FreeForm-refactor clones), fix
> `DBLE(nGridMat+1)/2d0`.  macos `708460f`; GMI/pymacos unaffected, grid
> test suites green.  Lightweight baselines regen'd + pushed
> (MACOS_res `78dc065` view_rx + `57a0965` e2e s4 dwdgrid).  **DEFERRED
> to after Friday:** e2e s6/s7 + e5_seg metopt regen.  **⚠ FLAG for
> Dave:** e5_pie committed baseline is STALE vs the e701f87 pie-frame
> convention (re-run reclocks xMon/yMon ~30°, unrelated to the grid
> fix) — regen as a separate reviewed step.
>
> **The e2e worked-example series (s1–s7) remains the active design
> thread; NEXT substantive item is still s7b** (unless Friday
> reprioritizes).  s1–s7 SHIPPED + COMMITTED on `sls-dev`; the s7 design
> + physics + numbers are in the "2026-07-20 / -21 (s7 SIMULATOR
> SESSION)" block further down this file.
>
> **NEXT = s7b:** upgrade the RBCS pose estimator from the static
> weighted-LS/BLUE form to the **steady-state Kalman filter** (Tesch
> *RBCS Algorithms* §2.3.3 eq 12-14, predict/update with the Riccati
> gain) and add **figure states** to the measurement model via the
> `dmdz`/`dmdgrid` blocks (macos.design.dmet_dfig) so the loop can
> SENSE and CORRECT the figure floor itself — not just the periodic
> image-based WFC.  The OSE single-step static estimator is this with
> converged gains.  Background PDFs in `MACOS_sandbox/Documents/`:
> `Tesch_RBCS_algorithms.pdf` (read), `OSE_Eqns_2019.pdf` (read),
> `2025_JATIS_HWO_Special_Issue-2.pdf` (PENDING — read before s7b).
>
> **Resume protocol:** read root+nested `CLAUDE.md`, then memories
> `[[project_recast_runners]]` (runners + RBCS loop + s7),
> `[[project_e2e_example]]` (s1–s6 heritage),
> `[[feedback_demo_plot_conventions]]` (movie plot rules), then the s7
> block below.  Runner = `design/runners/run_simulator.m`; driver =
> `examples/e2e/s7_simulate.m`; test =
> `tRunCompare/test_run_simulator_time_history` (SUITE_FAST).
> Everything OLDER than the s7 block (below, from 2026-07-16 on) is
> LANDED HISTORY — context, not in-flight work.

**2026-07-16: VISUALIZATION + MET-layout physicality/v2.** Session
record (chronological):
1. **Luis GridData transpose — DONE+PUSHED** (MACOS_res `17f1239`):
   `macos.read_grid_file`; engine GridInit reads file line=COLUMN;
   mmacos `elt_grid_add`=[x,y] vs pymacos=[y,x] OPPOSITE (memory
   `reference_grid_orientation_convention`).
2. **General visualizer — DONE, COMMITTED NOT PUSHED** (macos
   `032ddb8` + MACOS_res `d8011c6`): Dave: "work with any prescription
   — beam, optics, MET paths if present" (Lou's unfinished 3-D
   visualizer, modernized).  Engine: `Draw3DVec` 3-D DRAW capture
   (traceutil_mod; CTRACE fills at the 3 DrawRay sites via
   iDrawRay_global) + `draw_rays3d_get` + `met_geom_get` (endpoints ==
   met_get order; ride perturbations) + **OPD batch-hang fix: GO TO 6
   reprompt → abort — `trace(k)` on a Segment/NS elt spun FOREVER in
   SMACOS (same class as SPOT ed48d4f; sweep of remaining IACCEPT_S
   reprompt sites still OPEN)** + draw_rays_cmd vestigial stack args
   removed ("Unknown command" noise; load_rx still emits some —
   pre-existing).  mmacos: `view_rx` (DRAW-fan harvest, NOT per-elt
   trace(k)) / `draw_rays3d` / `met_geom`; design layer: `met_view`
   (3-D + face-on panels) / `hex_tile` / `seg_boundary` (hex + PIE
   wedges, arc-length `sample`); `segment_rx` exposes width/gap/grid.
   Examples: `examples/view_rx_demo` (Cass/Coro/e5mono+met) + e5_seg
   figures.  Tests: tMet 6/6, tMetView 4/4 (**findall can't reach
   sgtitle Text — title mirrored to fig.Name**), tSegmentRx 4/4,
   tReadGridFile 5/5.
3. **Dave's design-review constraints (figures drove these — the
   PURPOSE of the viz)**: boundary-true tiles (width/2 apothem, ONE
   global clocking per manual §segments — NOT per-seg face frames);
   launchers AT segment edges + edge_off 5 mm (add_met DEFAULT; edge
   placement alone: e5_seg edge+MET 5.98→4.34 nm); **MIN_SEP 50 mm
   between ANY two launchers (corner junctions!)**; **fiducials must
   MOUNT ON M2 ≤~25 mm inside its ApVec rim (615 mm here; r_fid≈590)**;
   **nf ≥3, likely 6**; **extra/M3 launcher ring hugs that element's
   physical radius + edge_off (was floating at r_fid)**; met_view hub
   disc now drawn at REAL hub radius (was fiducial-fit — masked the
   rim violation).  metopt v2: spread + CLUSTER-PAIR families,
   per-beam fiducial assignment enumeration (add_met `pair_map`),
   rim-zone RFID, MIN_SEP gate, hierarchical passes.  **DONE
   (3bdcda1+1103de0+369de7a)**: MIN_SEP killed the v1 corner ring
   (3 segments' launchers within 10 mm at every corner junction =
   INFEASIBLE); physical as-built prior 9577 → edge 3641 → MET 450 →
   edge+MET **232 nm** (MC 1.3%) — **aft/M3 ring at its physical
   100 mm = THE bottleneck**; optimizer (42-DOF SEGMENT sub-merit)
   → **3.358 nm**, winner spread [30 80 160]° + opposite-jump 6-fid
   map [1 4 2 5 3 6] at the rim (rfid=615), engine-FD 0.00%; the
   cluster family (64,800 layouts; kept FIRST-CLASS — deformation/
   rigid decoupling, Dave) lost to spread on rigid DOFs.
   PATTERN_FRAME 'radial'|'segment' knob (ring-uniform builder
   parts).  met_view readable views: face-on projected beams +
   color=segment association + **M2-M3 face-on inset** under the
   legend.
4. **e5_pie POLY-APERTURE ASSESSMENT DONE** (design/examples/e5_pie,
   6391b00; supersedes the sandbox copy): verdict YES — PolyApVec +
   explicit xObs emission loads clean; segments clip ZERO rays at
   nominal (fail_elt histogram); round-trip 6e-8 mm → 'rxpoly' reader
   trivial; pie seg_boundary/edge placement validated on the REAL
   fixture.  **Ray-loss question RESOLVED (Dave OPD review, MACOS_res
   f52ba4c): the RayFailElt=14 rays fall into the inter-segment GAPS —
   correct physical clipping, NOT an engine bug (elt-14 = cosmetic
   return-leg attribution).  Same commit: center segment (Elt 1) now
   emits a generated circumscribed 24-gon PolyApVec (Polygonal, no
   Circular special case); loss 456→442; all 7 polygons round-trip
   6e-8 mm.  rxpoly productization UNGATED.**
   Gotchas: macos.trace().nRays = SOURCE count (use ok_pass/fail_elt
   for parity); .presc segment blocks carry no ApType/nObs (append).
5. **Census for Dave's presentation**:
   `~/dev/MACOS_sandbox/improvements_census_draft.md` (~500 commits
   FreeForm→today, categorized, presentation arcs).  MEMORY.md index
   compacted 22→10 KB (no content loss).

**2026-07-16 (second session): the whole NEXT queue LANDED (all
local/unpushed on sls-dev — macos d09adc0+19d2c08, MACOS_resources
f52ba4c+10a955c+1dde42b+1dd4a64):**
1. **e5_pie manual example + aperture productization** (`10a955c` +
   `1dd4a64`): Dave's OPD review closed elt-14 (gap rays, correct
   clipping).  Center pie segment = HEXAGON from the traced footprint
   (poke-diff, NOT a disc/24-gon; corners at k·60°); ring-1 wedges
   abut it along straight CHORDS (flat (w−g)/2 + gap → chord (w+g)/2;
   obscuration = apex TRIANGLE, no inner arc — Dave's "circle around
   Seg 1").  `design.seg_apertures` (hex corners exact / pie hexagon +
   chorded sectors; xObs from psiElt NOT zMon) + `segment_rx
   emit_apertures/ap_pad/ap_obs` + `seg_boundary source=rxpoly`
   (auto when every segment declares PolyApVec; boundary = polygon
   minus obscuration via polyshape, LARGEST region — %.10E rounding
   leaves slivers, boundary #1 blind-take lost a tile).  Manual-grade
   runner `design/examples/e5_pie/e5_pie.m`, figure per step; pupil
   overlay needs affine centroid calibration (OPD grid transposed +
   mirrored; xGrid=(−1,0,0)).  met_view M2-M3 inset got axis labels.
   tSegmentRx 8 checks green.
2. **hub/aft DOFs in the analytic MET merit** (`1dde42b`):
   `design.met_bodies` = engine-truth frames (RptElt pivot + TElt
   rotation block via elt_rpt/elt_csys_get; segment triads reproduced
   EXACTLY).  metopt v3 merit = full 54 DOFs; fiducials RIDE the hub,
   aft ring rides the aft body.  Winner FLIPS to CLUSTER pairs at the
   corners (pmap [6 3 1 4 2 5], rim fiducials): rms 3.429 nm,
   engine-FD 0.00%; worst-mode 184 nm exposes the weakly-observed aft
   direction the 42-DOF merit hid.  tMet pins the 54-col identity.
3. **IACCEPT_S reprompt-loop sweep** (macos `19d2c08`): 21 sites in
   macos_cmd_loop.inc (BUILD/ORS/SRS/FEXIT/SXP/XPS/SPOT/MVAR/DVAR/
   GPERTURB/LPERTURB + AVAR's host-killing `stop`) now abort to the
   main loop; kept the 0=quit loops/MRESET/zernRange (self-
   terminating).  Both engines rebuilt, mmacos relinked, fast suite
   25 classes + tSegmentRx + tMet green; BUILD/SPOT at a Segment now
   abort in 0.01 s.

**2026-07-16/17 (third session): view_rx v2 — SOLID layout viewer
(Dave-planned + 6 review rounds; macos 30f4125+0c05406, MACOS_res
ac57bf4+b7c4c4b, all LOCAL):**
- Engine api: `ray_hist_set`/`ray_pos_hist_get` (traceutil RayPosHist
  — Lou's Vis3D substrate, now API-visible; slot 1 = source),
  `elt_info_get` (EltID/ApType/ApVec/xObs/lMon/PolyApVtx — IN-PLANE
  about VptElt, 5e-8 mm round-trip), `src_seg_get` (GridType/nSeg/
  width/gap).
- view_rx: rings-and-spokes FILLED bundle from the full traced grid
  (polylines connect the slots each ray REACHED — segmented ok is
  sparse), THIN flat two-tone sag-following shells (LightTools ref
  deck; no lighting), meridian profile curves both faces, sag sign
  calibrated vs crossings, consecutive-Refractor lens JOIN, EXACT hex
  Segment tiles (src_seg_get truth — no overlap, gaps read),
  source-plane ring (collimated-source cue), 'show'
  beam|beam+met|met, per-channel ray colors.
- macos.view_std (NEW): standard beam-aligned 4-panel figure (front
  from behind the SOURCE at the optic's face / back / iso / side),
  SOURCE AT LEFT, light → right, per-panel [az el] fine-tune;
  manual campos + camva('auto') + axis-off panels (hand-computed
  CameraViewAngle misframes; labels collide otherwise).
- Gotchas: ray_hist('on') must macos.modify() (grid-setter-retrace
  class); DRAW fan resamples its own rays (validate to grid pitch).
- Demo: 4 cases incl. e5hex1 via view_std; tViewRx 5/5 (SUITE_FREEFORM
  256), tMetView 4/4, tMet 6/6, quick 67/0.

**NEXT SESSION (Dave's objective, 2026-07-17): COMPLETE END-TO-END
worked example for users to hack, built from the parameterized
runners/utilities.  Dave 2026-07-17: EVERY stage runner must produce a
THOROUGH design report (saved text alongside the artifacts — identity,
first-order, field performance, and the stage's own metrics: Jacobian
condition/rank for stage 4, MET observability + post-control residual
for stage 5, time-history stats for stage 6) AS WELL AS graphics
(view_std/view_rx + the stage's metric figures).  Dave's SEQUENCE
(verbatim):
1. telescope design (TMA+FF feeding back end), with views and
   performance report;
2. add imaging instrument (3-4 mirrors to widen the field), new views
   and performance report;
3. segmentation, new views;
4. generate dwdx, dwdz, dwdgrid;
5. MET, MET-optimized performance report, dxdl, dxde, dwdl, dwde; new
   views with MET;
6. a SIMULATOR that generates PSFs using mmacos OR the linear model
   (user switch), driven by an x, z, grid TIME HISTORY.
Building blocks already shipped: Telescope/sz_tma builders +
design_report (1-2), segment_rx + emit_apertures (3), dw_d*_multi
harvests (4), add_met + metopt v3 + dmet_dx/dldx_analytic/met_bodies
+ edge_sensors (5; dxdl/dxde = the estimator gains from H, dwdl/dwde
= dwdx·gains), COMPOSE/psf + linear w = dwdx·x + dwdz·z + dwdgrid·g
(6).  view_std/view_rx figures at every stage.

Smaller queue: GMI mex not relinked against the new engine (next
makeall); pymacos export of the 4 new api routines deferred; rxpoly
for IMPORTED segmented Rx untested on an external fixture;
PLAN_DESIGN_LAYER promote-on-land.  Full record + gotchas: memory
`project_met_visualizer.md` + `project_e5pie_apertures.md`.**

**2026-07-17/18 (e2e sessions 2-3): STAGES 1+2 LANDED + VIS REDESIGN
(MACOS_res `e20d6b1` -> `5e94320` (s2 v1) -> `cbbb7c9` (VIS), all
PUSHED except cbbb7c9 LOCAL).  ARC: s1 v1 (f/1.25, m2=16, fold, bias
5') was NOT a VIS imager (Dave ran the .in) -- bias was the killer
(~bias^2).  VIS design point (Dave): f/1.75, m2=8, M3 EXTRACTION TILT
1.2 deg (return off the feed axis -> fold clears at bias 2').  s1 =
DL telescope (+-1' 0.0157 -tilt waves @500nm); s2 = 3-mirror bench
relay (M4 corrector/M5 collimator/M6 camera, radii DERIVED from the
collimator condition), +-2' at 0.26 -tilt, distortion 0.28" (M4
blur-guarded distortion stage stands down).  PROCEDURE (Dave: record
everything) = README 'design procedure' section, 11 rules each from a
failed run: WFE-based bias pick (small-bias conic basin K3->-3..-4 is
repeatable, continuation does NOT rescue), guarded M1 common-mode
null, SVD engine for degenerate bases (CALIB SIGSEGVs), joint-then-
null order, detector-plane re-fit last, distortion!=blur (M4 near
focus = the reflective distortion corrector; affine-projected chief
metric; blur-guarded).  DEAD ENDS on record: per-field Zernike patch
corrector (rank collapse), 4th near-pupil mirror (common-mode), full
BornWolf 3:25 (+13% only), relay tilts 6 deg (in/out beams overlap).
m2=12 A/B in ~/dev/MACOS_sandbox/e2e_m2_12 (bias 4', 0.25 waves).
OFFNER RELAY SHIPPED (MACOS_res 4ab5229): DL over the FULL +-2'
(Strehl 0.99->0.93, -tilt 0.015-0.043, pure spheres) via NEW
design/src/offner_layout.m (concentric chief solve -> Bauer chain;
spheres held, convex stop mirror out of solves, M4 = flat routing
fold).  Zigzag variant kept selectable (DL to 1.5', a7ed2f6).

**2026-07-18 (e2e session 4): FIELD SHIFT + HOLE CHAIN + s3 LANDED.**
(1) Field center -0.7' off the s1 bias (Dave read the s2 WFE map;
P.inst.field_dy_arcmin applied to set_field_bias -> artifact chief =
science center, s3-s6 inherit): worst +-2' 0.043 -> 0.0231 -tilt,
Strehl floor 0.965; [4f] center scan in runner+report; full -1.05'
re-solve = flatter interior / worse 2' edge, BOTH kept in
s2_variants/dy-0.70|dy-1.05 (Dave: keep both).  (2) M1-hole chain
(Dave: show the hole in all layout views): set_hole EMITS a real
ObsType=Circle obscuration (trace clips the 5 central rays); NEW
engine api elt_obs_get (macos_api_mod, codegen'd; pymacos export
deferred); view_rx renders circular obscurations (view_std inherits);
segment_rx carry_obs -> center segment.  (3) s3_segmentation.m: pie
(7, e2e_pie.in) + hex2 (19, e2e_hex2.in), physical apertures, 128-pt
grid (Dave: 41 too coarse), P.seg.variant="pie" feeds s4-s6; pie
12090/12520 pass (392 gap/rim), hex2 9838/9876 (hex tiling has no gap
rays); both VERIFIED.  TWO segment_rx product bugs fixed (tests
green): Surface=Zernike parent FIGURE dropped by SegMirMaker (15 um!)
-> carried into each segment's FF channel; SegMirMaker<->engine
tiling contract needs xGrid=(-1,0,0) AND SegXgrid=(-1,0,0) in the
merged Rx (design-layer (+1,0,0) left ray->segment 180 deg off the
frames; PSEG ignores SegXgrid = dead code, HSEG anchors to it).
**LATENT e5 FINDING for Dave: the e5 hex corpus has the SAME 180 deg
frames<->tiling offset (ray_hist-verified) — invisible bare, but
e5_seg's forward model pairs dwdx (tiling identity) with dedx/dldx
(frame identity) point-reflected; fix = SegXgrid -1 / regenerate,
touches committed Sprint-2D references — Dave's call.**
center_focal_plane now engine-truth (trace+get_ray_info; was
draw_rays plot-projection, U-sign follows grid handedness — why the
emitter KEEPS +1 and only segment_rx's output flips).  Fast suite
229/0; tSegmentRx 8/0 (+figure/obs-carry tests), tDesignTelescope
67/0 (+hole test).
**2026-07-18 (cont.): GENERAL FIX LANDED IN THE GENERATOR (Dave:
"fix it in e5 -- preferably in a general way").**  SegMirMaker.f, two
fixes (MACOS_res d54405f): (1) header SegXgrid now emits the in-plane
basis ACTUALLY USED — the back-facing-mirror 180-deg basis flip
(segment-numbering convention) negated xs/ys for frames+SegCoord but
the header kept the pre-negation vector → engine HSEG tiled rays
point-reflected from the frames (EVERY back-facing fixture, e5
included).  (2) Zernike-aware LoadParent: Surface=Zernike parent
figure merged into the FF channel (segments present the as-designed
surface; RptElt placement sees the figure sag).  Byte-identity refs
regenerated (delta = SegXgrid line only; Hx identical).  NEW
tSegmentRx gate test_state_consistency_rays_on_frames (ray_hist truth
vs frames, Pie+Hex, every off-center segment) = THE invariant: seg
k's DOFs move seg k's w/e/l only.  e5 corpus REGENERATED consistent:
e5_seg edge+MET 229.5 nm (was 232), metopt 3.421 nm (was 3.429),
engine-FD 0.00%; e5_pie unchanged behavior.  segment_rx text-carry
stands down when SMM already carried the figure.  Suites: tSegMirMaker
3/0, tSegmentRx 9/0, tMet 6/0, tEdgeSensors 3/0, tMetView 4/0, fast
230/0.

**s3 GRAPHICS + PIE-GEOMETRY FIXES (Dave review, 2026-07-18).**
Dave flagged the s3 views (EP dashed circles E15/E16, unclear pie
center segment) — pulling the thread found TWO real geometry bugs in
the pie aperture emission, not just rendering:
- **view_rx polish:** Return elements HIDDEN by default (`'returns',
  true` restores; tViewRx test); Segment plates draw with alternating
  face tints, center cell (Seg1, carries the hole) distinctly darker,
  and NO per-tile profile curves (7–19 overlapping meridian-curve
  sets were the spoke mush hiding the center segment).
- **Pie wedge apertures were apex SECTORS to r=0** (+ triangle
  obscuration faking the chord): every declared-polygon consumer
  (view_rx plates, step-3 figures) drew wedges COVERING the center
  hexagon with edges converging at the center.  A ring-1 wedge is
  CONVEX (disc ∩ 3 half-planes) → now emitted as the TRUE chorded
  polygon, no obscuration.  **Gaps are uniform-width slots** (side
  edges PARALLEL to the sector rays at ±g/2, not angular offsets that
  converge at the center — Dave).  New shared helper
  `macos.design.pie_wedge_geom` used by seg_apertures + seg_boundary
  + view_rx (ONE geometry source, no third desync).
- **pie_rings classification robustness:** e2e's Zernike-figured
  parent scatters ring-1 wedge radii ~5 µm (frame-tilt leak in the
  tiling-plane projection); the old 1e-6·max(rc) tolerance split the
  6 wedges into degenerate 1–2 member "rings" → 2π/nnz sector spans +
  go/sin(π) 1e14-vertex blowups in the emitted .in (the broken
  s3_footprints_pie Dave saw).  New `macos.design.pie_rings`
  (width-scaled tolerance, shared by all three sites) + jitter test.
  e5 never tripped it (exact symmetry) — the e2e example did its job.
- Regenerated: e5_pie (wedges 15 verts, 502 clips, no obscurations)
  + e2e s3 (pie 12106/12520, 376 gap/rim clips ≈ prior 392; hex2
  unchanged 9838/9876; both VERIFIED).  tSegmentRx 11/0 (+2 geometry
  tests), tViewRx 6/0 (+returns test).
- MET is NOT in s3 by design — metrology config + optimization is s5.

**s3 GRAPHICS REVIEW → PIE GEOMETRY FIXES + HELPER HOIST (Dave,
2026-07-18).**  Dave's s3-view flags exposed real emission bugs:
- **view_rx:** Return elements HIDDEN by default (E15/E16 EP dashed
  circles; `'returns',true` restores, tViewRx test); Segment plates
  alternate face tints (center cell darkest, carries the hole), NO
  per-tile profile curves (the spoke mush).
- **Pie wedges were apex sectors to r=0 + triangle obscuration** →
  declared-polygon consumers drew wedges covering the center hexagon.
  Ring-1 wedge is CONVEX → emitted as the TRUE chorded polygon, no
  obscuration; **gaps = uniform-width slots** (side edges PARALLEL to
  sector rays at ±g/2 — angular offsets converge at the center,
  Dave's rule).  Shared `macos.design.pie_wedge_geom` (seg_apertures
  + seg_boundary + view_rx = one geometry source).
- **`macos.design.pie_rings`:** ring classification needs a
  WIDTH-scaled tolerance — the e2e Zernike-figured parent scatters
  ring-1 wedge radii ~5 µm (frame-tilt leak in the tiling-plane
  projection); the old 1e-6·max(rc) split them into 1–2 member
  "rings" → 2π/nnz spans + go/sin(π) 1e14 vertices in the .in.  e5's
  exact symmetry never trips it — the second consumer found it.
- **Hoist (Dave: general-purpose runners):** `macos.design.
  seg_footprints` (poke-diff engine-truth footprint measurement) +
  `seg_footprint_view` (calibrated overlay figure) — e5_pie.m and
  s3_segmentation.m now thin narratives over them; duplicated
  pupil_axes_ copies deleted.  Regens BEHAVIOR-IDENTICAL pre/post.
- Regenerated: e5_pie (15-vert wedges, no obs, 502 clips) + e2e s3
  (pie 12106/12520, 376 gap/rim clips; hex2 9838/9876; both
  VERIFIED).  Suites: tSegmentRx 12/0 (3 new tests), tViewRx 6/0,
  tMet 6/0, tMetView 4/0, **fast 232/0**.
- **PLAN §0 model-transition heap crash REOPENED** (Dave asked; the
  full no-arg mmacos suite SIGSEGVs at tViewRx setup after ~26
  classes in one process — evidence + repro + suspects recorded in
  PLAN.md; `fast` subset is the green signal until refixed).
- MET is NOT in s3 by design — metrology config + optimization = s5.
- **s2 AUTO FIELD CENTER LANDED (Dave's ask):** sections [2]–[4e] now
  run in a 2-pass loop; [4g] maps ±3′ (13×13), finds the centroid of
  the raw-WFE<0.02 region, and if the chief is >0.15′ off it, pass 2
  re-solves there.  Result: chief moved another −0.71′ (bias 1.3′ →
  **0.59′ total**), worst ±2′ −tilt **0.0231 → 0.0215**, <0.02-region
  31→52/169 pts, centroid residual +0.03′; [4f] scan confirms the new
  center is the optimum (both ±0.35′ neighbors worse).  Clearance
  UNCHANGED by engine-truth ray-to-body margins (M2 −230 mm input-on-
  M2-back obscuration + FM mm-graze, both documented; M7 physical
  +126 mm — the report's new M7 flag is pupil-retrace bookkeeping,
  not light).  s2_wfe_field.png now ±3′ with the 0.02 contour +
  science patch + chief/centroid markers.  s3 rerun on the new
  parent: pie 12098/12520 (376 gap clips), hex2 9830/9876, both
  VERIFIED.  Old hand-picked variants preserved in s2_variants/.
- **s4 LANDED (MACOS_res 89f06f8): `s4_jacobians.m`** — dwdx (60412×90,
  rank 90; spectrum shows EXACTLY 21 strong modes = 7 segs ×
  piston/tip/tilt then a 2-decade cliff), dwdz (×56, FFZern 4..11 per
  seg, cond 4.4e2), dwdgrid (×42, 6 aperture-confined pokes/seg on
  grid-augmented `e2e_pie_grid.in` — flat 256 grids in each segment's
  CLOCKED Mon frame; cond 1.7), all via the production dw_d*_multi
  supervisors over C+4 corners at ±2′.  Per-segment column-norm table
  (piston dominant, clocking ~null) in s4_report.txt; s4_jacobians.mat
  carries the three outputs for s5/s6.
  POST-s4 FIXES (Dave review, e701f87): SMM pie frames = axis of symmetry (wedge bisector / center flat-normal; hex heritage unchanged; refs regen, tSegMirMaker 3/0); s4 dwdz = monzern (segment-LOCAL, cond 438→5.4; ffzern = parent-aperture basis, wrong channel); s4 dwdgrid = segment_grid_basis G-S (cond 1.26); segment_grid_basis ray-hist fallback for Segment pm_ref_elt; LESSON = copy the tested sensitivities runners CONFIGURATION (run_dwdz, run_dwdgrid_multi_multisegbasis), not just the API.
  NEXT: s5_met.m — MET config with the SHAPE-CLASS launcher constraint
  (below) + dedx/dldx joining s4's dwdx; then s6 simulator.

**2026-07-19 (s6 SESSION, post-compact resume point): run_compare
SHIPPED — the compare stage runner + s6 driver + tRunCompare.**

- `design/runners/run_compare.m` (Dave's spec): pokes each rigid DOF
  of the sensed bodies (100 nrad / 100 nm) on the optimized met Rx;
  per poke TWO graphics — mmacos ENGINE | LINEAR MODEL — each the
  center-field OPD change (parula, shared clim) above stacked bars of
  l / e_piston / e_gap / e_shear (shared ylims); dwell 1.6 s (Dave:
  0.25 read too fast — runner default + driver + GIF DelayTime all
  1.6); frames in <name>_frames/ + <name>_compare.gif + agreement
  report + <name>_compare.mat (exports dwdu = segments+SM columns for
  s7; ~7 MB, committable).
- Engine truth per channel: w = RigidBodyChannel apply → trace(nElt−1)
  → OPD delta (same convention as the s4 harvest; plain bodies need no
  FP tracking); l = macos.met() (met points ride the element under
  programmatic perturb); e = FINITE-ROTATION kinematics at the Hx
  SensorPos points — axis recovered from the row's own translation
  block (a_w = T_s·blkᵀ), arm = SensorPos − rpt_s, full Rodrigues;
  first order == the dedx row, difference = true linearization error.
- **STALE-dedx catch (the debugging story):** first run showed e
  disagreeing O(1) on in-plane DOFs while translations are
  algebraically forced to agree → the s5 e2e_pie_met.mat dedx was ONE
  Hx GENERATION BEHIND (s5 ran before the gap/shear axis-relabel
  regen; piston rows agreed, in-plane pair rotated).  run_compare now
  builds dedx FROM THE Hx SIDECAR (pad + rot cols ×cbm, same as
  run_met) and only cross-checks the met .mat, warning
  'run_compare:stale_dedx'.  s5 rerun queued/in-flight to refresh the
  .mat estimator products (dxde/dwde) for s7 — WEM numbers are
  rowspace-invariant, layout unchanged.
- **s6 RESULT (e2e pie, 54 DOFs): engine == linear to ≤9.1e-5 (w),
  ≤1.2e-6 (l), ≤6e-8 (e rot; translations to machine 1e-16).**
  Null-response DOFs get a floor (opts.w_floor 1e-12 m): segment Rz
  (clocking through the near-flat parent leaves ~0.01 nm — REAL, above
  floor; true nulls are the flat fold's in-plane DOFs + FPA-class
  bodies at ~5e-15 m FD noise) → reported 'null', not a noise ratio.
- Bodies without dwdx channels (aft ring riding the FPA in the e5
  fixture) get ZERO linear columns — matching the engine's null plain-
  trace response; control bodies are asserted actuatable.
- tRunCompare (SUITE_FAST): e5 pie fixture, REAL single-field coarse
  harvest ('grid','1x1', ngridpts 15, fp_mode none — no ApStop in the
  SMM corpus) + run_met-with-no-jac met struct; gates = the physics
  (w<5%, l<1%, e<1e-3, FPA null, dwdu 48 cols).  36 s.  PASSES.
- **s7 NOTES (Dave, this session): single-step STATIC estimator** —
  the OSE form x̂ = x̄ + K(m−m̄), converged steady-state gains (= the
  dxde/dxdl MMSE partitions run_met exports), m = [l; e]; OSC
  controller u = −[c_wu·I + (DwDu)ᵀDwDu]⁻¹(DwDu)ᵀ·dwdx·x̂.
  Background papers: MACOS_sandbox/Documents/OSE_Eqns_2019.pdf +
  2025_JATIS_HWO_Special_Issue-2.pdf (read OSE; JATIS pending).
- s5 rerun DONE: e2e_pie_met.mat refreshed on the current Hx — WEMs
  reproduce EXACTLY (as-built 15.57/3.918, optimized 3.766/1.777, FD
  0.00%, MC 1.9%) confirming rowspace invariance of the relabel.

**z + grid EXTENSION + FIGURE SENSING (Dave, same session):** "add
dwdz and dwdgrid visualization to s6 — keep going after dwdx, z then
grid" and "evaluate the effect of grid and z DOFs at each sensor
location (l and episton) → generate dmdz and dmdgrid":
- run_compare now phases x → z (MonZernChannel pokes on the met Rx) →
  grid (GridChannel pokes on og.rx_path, basis REBUILT via the same
  segment_grid_basis call — mismatch fails the agreement gate).  One
  GIF, frames p%03d, per-channel worst summary.
- **macos.design.dmet_dfig** (+ macos.design.zern_seg_eval): dmdz /
  dmdgrid = mode shape at each Hx SensorPos + met_geom launcher point
  (src_elt gives the launcher→segment map engine-truth), projected on
  the measurement axis: dl = −û·n̂·f(p_src), de = (â·n̂)·f(p_q).
  Piston + l rows carry it; gap/shear are SLOPE-ORDER small (~2-3% on
  the e5 f/1.75 parent), NOT zero — in-plane axes are ⊥ the LOCAL
  normal, Mon sag displaces along the FACE normal zMon (real model
  response; flagged for Dave's model review).  Exported in
  s6_compare.mat as [l;e]-ordered dmdz/dmdgrid for the s7 H.
- **zern_seg_eval convention PINNED BY ENGINE GATE** (tRunCompare
  test_zern_grid_engine_equivalence): the same MonZern mode sampled
  onto a grid channel + poked via elt_grid_add reproduces the
  MonZernCoef poke's OPD (scale within 2%, corr 0.998 = grid
  discretization) → lMon-normalized, UN-normalized ANSI (MonZernType=
  ANSI; NORM_RMS gates to Norm* types only — the E1 fix), Mon frame =
  clocked face frame.
- **UNIT GOTCHA (found via 999×≈1/cbm mismatch on the mm fixture):**
  the dwdz/dwdgrid supervisors (dwdz_for_current_source) return
  OPD-BaseUnits per coef-BaseUnits — NO cbm scale, UNLIKE dwdx's
  OPD-metres convention.  run_compare ×cbm on the figure-poke linear
  maps; invisible on e2e (m), 1000× on e5 (mm).
- The engine METcalc/Hx hold met + sensor points RIGID → z/grid
  engine bars read zero vs the dmdz/dmdgrid linear bars — the movie
  shows the model-vs-engine gap explicitly (labeled 'l mdl'/'e_pist
  mdl' in the report; natural follow-up = engine met points riding
  figure, Dave's call).
- tRunCompare EXTENDED (z leg with a real monzern mini-harvest — the
  e5 FreeForm REFRACTOR elt 9 legitimately sweeps in, counts dynamic
  — + the engine equivalence gate): 2/2 green, ~60 s.
- **.doc REPORTS (Dave):** MACOS_sandbox/e2e_reports/make_reports.py
  (pandoc) → s1–s6_report.docx: front matter + stage summary + stage
  figures + verbatim report; s5 carries the MET-OPTIMIZER DEEP DIVE
  (merit/WEM algorithms, shape-class search flow, code flow, results
  tables — s5_met_optimizer_deepdive.md).  s1–s5 built; s6 after its
  rerun.
- **GRID-BASIS PERSISTENCE (the deep find):** run_compare's grid
  agreement exposed that the s4 og dwdgrid block could NOT be
  reproduced by a rebuilt basis — modes 1–5 bit-exact, but each
  segment's LAST G-S mode came out ORTHOGONAL across sessions (|Δ|=√2;
  in-session double-builds bit-identical; raw detrended modes well-
  conditioned sv≥0.19 — a discrete tie-break below the G-S waterline,
  root cause unchased).  Fix = ARCHITECTURAL: the influence basis is
  PART OF THE JACOBIAN'S DEFINITION → run_sensitivities persists it
  (og.sgb) and run_compare/dmet_dfig consume it VERBATIM (rebuild
  path kept with a loud warning).  Also: run_compare must RELOAD the
  grid Rx after any segment_grid_basis call (it traces/ray-hist-pokes;
  leftover state biased grid FDs ~20% — the harvest supervisor
  reloads, manual loops must too).
- **FINAL s6 (persisted basis): x ≤2.8e-3, z ≤1.8e-3, grid ≤0.235.**
  The grid residual is an **OPEN ENGINE QUESTION for Dave**: the
  engine's grid-surface response is NOT proportional across sub-µm
  amplitudes on METRE-based Rx — central FD at the harvest delta
  (1e-6 BU) reproduces the Jacobian to 1.7e-5, but forward pokes give
  dW(h)/h drifting 0.15→0.24→0.42→0.98 over h=1e-8..1e-6 and BLOWING
  UP 84× at h=1e-5 (10 µm figure → ~570 µm OPD response?!); fex
  referencing ruled out; the mm fixture at the same physical poke
  agrees <5% → units-dependent (absolute-EPS class?) behavior in the
  FreeForm grid intersection path (GSZPSolve/SFFZPSolve suspects).
  x/z channels are clean at the same amplitudes (z fwd 1e-7 vs 1e-6
  agree to 1.5e-3).  run_compare report carries the note.
- Fast suite FINAL: **242/0** (tRunCompare = 2 tests, x/z/grid legs +
  engine convention gate).
- .doc reports (Dave): s4 doc gained the SENSITIVITIES GENERATION
  STORY deep dive (channels/referencing/grid substrate/basis
  persistence/unit conventions/validation); standalone CONCISE
  met_optimizer_concise.docx (with the optimizers-used note: discrete
  block coordinate descent + exhaustive per-class enumeration +
  top-K refinement, closed-form MMSE inside, no NLP solver).
  All in ~/dev/MACOS_sandbox/e2e_reports/ (make_reports.py).

**2026-07-20 (s7 SIMULATOR SESSION): run_simulator SHIPPED (Dave's
two-part spec) + the grid non-proportionality question ANSWERED
(numerical, not small-angle).**
- **Q2 answered (Dave asked if the grid rel 0.235 is just too-large-
  poke nonlinearity): NO — three discriminators.** (1) 100 nm on the
  8 m aperture is ~12 ppb sag and the z path is proportional to
  1.5e-3 at the same amplitudes; (2) the FD ratio does NOT converge
  as h→0 (0.15 at 1e-8 → 0.98 at 1e-6) — real nonlinearity improves
  as h shrinks, a granularity/EPS floor doesn't; (3) the mm fixture
  at the SAME PHYSICAL poke is clean <5% — physics can't see the Rx
  units, an absolute solver tolerance can.  Engine investigation
  still open (GSZPSolve/SFFZPSolve).
- **run_simulator (design/runners/, SHIPPED): time-series x/z/grid →
  movie, UNCORRECTED + CORRECTED legs** (Dave 2026-07-20: history
  opens with µm-to-mm misalignments; initial image-based WFC
  u = −pinv(dwdu)·w(frame 1) solved from the ENGINE wavefront and
  HELD; then nm-to-µm random-walk steps; 1000 s at 1–10 s steps).
  Frames: OPD unc | OPD corr | PIX psf | COMPOSE broadband psf +
  m=[l;e] bars (piston/gap/shear colors) + ACCUMULATING rms-WFE
  (log, both legs) and Strehl curves (unc λ0/corr λ0/corr bb, peak
  ratio to nominal).  m bars = the validated linear model
  [dldx;dedx]·(x+u) + dmdz·z + dmdgrid·g (dedx rebuilt from Hx;
  dmet_dfig blocks; grid Rx has no metrology — engine-l cross-check
  runs when the history has no grid states).  Per-frame corrected-leg
  engine-vs-linear w_rel = the s6 gate extended to MIXED states,
  computed on the DRIFT INCREMENT (frame t − frame 1): the absolute
  frame-1 state is deliberately nonlinear at µm-mm amplitudes and the
  control compensates it, so an absolute comparison would only
  re-measure that compensation.
- **WFC ITERATES (wfc_iters, default 3):** one linear solve at 202 µm
  left an ENGINE residual of 1.3 µm vs 40 nm predicted — the ~0.5%
  per-column nonlinearity of µm-mm states, not a solve error.  Real
  image-based WFC iterates; each Gauss–Newton refinement re-measures
  the engine wavefront at the corrected state and ridge-solves the
  update (monotone state path — never toggles back, respecting the
  two-pass non-closure rule).
- **TWO-PASS ENGINE SCHEDULE (correctness find):** toggling ±u per
  frame through the incremental perturb path does NOT close — fixed-
  order single-axis rotation increments leave a SYSTEMATIC ~|u_rot|²
  non-closure per cycle that accumulates LINEARLY (~µrad phantom over
  100 frames at 50 µrad u).  run_simulator plays the whole
  uncorrected history (storing per-frame OPD + psf peak), reloads the
  Rx, then plays the corrected history — within a pass the large
  state applies once and increments are nm-scale.
- **WFC SOLVE = TIKHONOV RIDGE (both failure modes hit, fixture +
  e2e):** plain pinv (tol 1e-6) noise-amplified near-degenerate dwdu
  combos (piston-like Tz directions) into huge canceling commands —
  engine honors them only to first order → WORSE psf (0.72→0.21)
  despite lower fitted rms; a hard SVD cutoff at 1e-2 over-truncated
  (624 nm of the 202 µm e2e initial error uncorrected → corr Strehl
  0.25 at 500 nm).  Ridge u = −V·(sv/(sv²+λ²))·Uᵀw, λ = wfc_tol·s1
  (default 1e-3; e2e driver tunes 3e-4 → predicted residual ~30 nm)
  corrects every direction sv ≫ λ and bounds each command by
  |w|/(2λ).  This IS the OSC controller form — carries straight into
  the estimator/controller loop.
- tRunCompare/test_run_simulator_time_history (SUITE_FAST, fixture
  reuse): µm-scale opening state + drift + z + g; gates = control
  collapse (corr < 0.2×unc), Strehl improvement, mixed-state w_rel
  <0.15, artifacts/m_hist/dmdz/dmdgrid dims.  3/3 green.
- **METROLOGY LOOP = RBCS estimator/controller (Tesch "RBCS
  Algorithms" ch 2.3 + 3.3; several live-review corrections):**
  `met_loop` (default on): the post-WFC state IS the control TARGET;
  each frame the sensed drift δm = m − m_ref drives a POSE ESTIMATOR
  then a CONTROLLER (Dave: "estimate the state without weighting the
  WF impact").
  - **Estimator = weighted-LS / BLUE (Tesch §2.3.2 eq 11):**
    δx̂ = R_meas·δm, R_meas = (Hᵀ N⁻¹ H + R_x⁻¹)⁻¹ Hᵀ N⁻¹, H = dmdx,
    N = sensor-noise cov (default [1 pm laser truss, 1 nm edge]),
    R_x = state/disturbance prior (default "auto" from the drift
    std, PTT ≫ lateral ≫ pinned-aft).
  - **Controller = min pose error (Tesch §3.3.1 eq 16-17):**
    u_t = u_{t−1} − k_p·δx̂(control DOFs), k_p 0.5 (<1 for margin;
    the TCE integrator makes the loop robust to gain error).
  - **THE BUG (Dave: "MET control is not correctly implemented" →
    pointed at Tesch):** my first loop used a RAW pinv(dmdx) = the
    basic-LS estimator (Tesch §2.3.1 eq 10).  It amplifies un-modelled
    δm content (figure drift, linearization residual) by 1/σ_min in
    the weakly-observed rigid directions; the integrating loop RAN
    AWAY — the engine e2e went 0.02 nm at t=10 → 2e5 nm at t=20 →
    5.3e6 nm.  Standalone linear reproduction: basic-LS max 5.6e6 nm
    WF residual / 34 mm commands vs BLUE max 6.98 nm / 74 nm commands
    on the identical drift (open-loop drift was 81.6 nm).  The R_x
    prior is the fix: weak/pinned DOFs fall back to prior-0 instead of
    being inverted.  This is STATE weighting (noise + disturbance
    stats), NOT wavefront-impact weighting (Tesch §3.3.2 eq 19,
    deliberately unused per Dave).
  - Bars show the SENSED DRIFT δm, not absolute m (µm-scale post-WFC
    offsets swamped the nm drift — Dave's "MET results not changing").
    Figure drift ALIASES into x̂ via the l/e_piston rows (no figure
    states in the simple estimator) — negligible at realistic figure
    drift; s7b's H adds dmdz/dmdgrid.
  - Loop-STABILITY gate added to tRunCompare: corrected rms ≤
    uncorrected rms every frame (the old pinv went 100× above).
- **DOF-class statistics (Dave):** support structure puts 10× more
  error in PTT (local Tz/Rx/Ry — the class MET+edge find and correct)
  than lateral (Tx/Ty/Rz); equal allocation had loaded the ridge-
  dropped weak directions and left a pessimistic ~40 nm one-shot
  residual.  **Figure drift = few nm TOTAL** (1 nm/step had the 98
  coefs each walking to ~10 nm → >100 nm rigid-uncontrollable figure
  dominating the corrected leg's drift-away).
- e2e driver s7_simulate.m (final): ~1 µm PTT initial (0.25 µrad
  tip/tilt, 1 µm piston; lateral /10), drift 4 nrad-nm/step (5×
  reduced; lateral /10), figure 0.004 nm/step, T=100 @ 10 s; M3→FPA
  + aft ring pinned (truly stable structure, Dave).  Running WFE
  printed in the OPD panels' xlabels (Dave).  Stale-frame cleanup
  built into the runner (s6 lesson).
- **WF-MAINTENANCE RECONTROL + TWO SCENARIOS (Dave 2026-07-21):**
  `wfc_reset_times` re-runs the image-based WFC mid-history (Tesch's
  WF Maintenance Activity) + `wfc_on_frame` delays the initial WFC so
  the movie opens uninitialised ("no system starts perfect" — first
  two data points at the as-deployed ~100 µm, then control turns on).
  Two 500s runs, reset@400s, one GIF each (s7A/s7B):
  - **s7A metrology-bias:** the metrology zero-point drifts; the loop
    holds the BIASED reading so the true wavefront walks off UNSEEN
    (`meas_bias = −dmdx·p(t)`); the 400s image-based reset
    re-references and knocks it back.  ENGINE: corr holds ~1 nm →
    **44 nm by 390s → reset → 1.0 nm** → 12 nm by 500s.
  - **s7B focus/astig figure:** per-segment focus(5)+astig(4,6) trend
    (60 nm, 2× for visibility per Dave); the truss reads RIGID POSE
    only (`loop_senses_figure=false`) so figure accumulates unseen; the
    tight-ridge reset (`wfc_reset_tol=1e-5`) engages the LATERAL DOFs.
    ENGINE: 20 nm → **45.6 nm by 390s → reset → 24.4 nm** → 31 nm.
  - **PLOT conventions (Dave, [[feedback_demo_plot_conventions]]):**
    delayed init (first 2 pts at as-deployed 100 µm); Strehl EXACT from
    OPD `|<exp(i2πW/λ)>|²` (≤1, not psf-peak >1); autoscale the
    corrected OPD panel (shared 100 µm scale hid the nm structure);
    broadband Strehl trace dropped; legends non-obscuring (hide the
    reset-marker auto-legend entry); running WFE in the OPD xlabels.
  - **KEY PHYSICS (Dave, corrects my earlier error):** RB control CAN
    counter SEGMENT focus/astig on a parabolic parent — a segment x
    move changes its local best-fit radius (focus + a bit astig),
    y/twist add astig.  Verified on the s4 dwdx: **focus RB-residual
    0.017, astig 0.31–0.49 (via lateral DOFs), higher order 0.6–0.98
    (uncorrectable)**; the reset needs tol 1e-5 to reach the weak
    lateral DOFs (3e-4 leaves astig 0.78).  My earlier per-segment
    reset failed 24→24 because the loop sensed+pre-consumed the figure
    AND the loose reset tol missed the lateral DOFs.
  - There is NO wavefront error BOTH MET-invisible AND rigid-
    correctable in a well-instrumented segmented system (rigid-
    correctable ⟹ segment pose ⟹ MET-visible); the WFC reset earns its
    keep against (A) MET bias drift and (B) figure the MET is DEFINED
    not to sense.
- NEXT (s7b): upgrade the estimator to the STEADY-STATE KALMAN form
  (Tesch §2.3.3 eq 12-14, predict/update with the Riccati gain) +
  figure states via dmdz/dmdgrid in H + sensor noise; the OSE
  single-step static estimator is this with converged gains.

**2026-07-19 (RECAST SESSION, earlier): the
sensitivity diagnostic + runner recast + SMM EDGE-SENSOR REWORK all
LANDED (pushed: MACOS_res 20d7f03+b41503d, macos 3a89432).**

DIAGNOSTIC (Dave's four questions, all closed):
- s4 == the runners BIT-IDENTICAL (same dw_d*_multi supervisors).
- Segment coordinates CORRECT (engine-truth poke footprints match
  frames on e2e AND e5pie; apparent reflection = canvas convention).
- Nominal subtraction CORRECT (reference frozen per field; median
  col-vs-w0 corr 1-3%; top col Elt 17 Rx 0.63 = physical FP signature).
- "e5 doesn't compare": the E5 CORPUS is the trap, not s4 --
  SegMirMaker replicates the PARENT's grid channel (pData=0, full-
  aperture span) into every segment block; segment-frame bases poked
  against those paint a CENTRAL DOT and rank-collapse (e5pie dwdgrid
  rank 15/42 cond 1e7 vs healthy 42/1.26).  ALSO: SMM corpus ships NO
  ApStop (exit-pupil machinery fails).  Last-key-wins parsing means
  appending correct grid lines cannot fix stale ones.

RECAST LANDED (MACOS_resources, all local):
- macos.design.grid_augment_rx: per-segment grid channels in the
  CLOCKED Mon frames, REPLACING stale lines.  Span default = the
  dxGrid convention (Dave): GridSrfdx = Aperture/(ng-1) (span = the
  beam; e5mono heritage 31.25 = 8000/256).  NEVER size from lMon:
  pie-wedge lMon is the hex 'length', not a circumscribing radius --
  a 2.2*lMon span clipped wedge corners and made s4 dwdgrid
  non-physical (Dave caught it; segment_grid_basis now WARNS when
  footprint rays fall outside the grid span).
- design/runners/run_sensitivities.m: general stage runner (dwdx/
  dwdz/dwdgrid + dwdsurf opt-in; grid_basis multi|single; influence
  passthrough for DM maps; ApStop injection preflight 'stop' opt;
  zkinds; segment-only conditioning; per-segment column norms).
- design/runners/run_segmentation.m: parent .in -> verified segmented
  .in (parity, engine-truth footprints, apertures, reload gate).
- s3/s4 = thin drivers; regen VALIDATED (e2e_pie.in/e2e_hex2.in
  byte-identical; s4 dwdx/dwdz numerics unchanged, dwdgrid healthy).
- PLOT CONVENTIONS (Dave): piston removed in the shared
  plot_dw_per_element/plot_dw_channels (parula stays); per-element
  pages center+multi in <name>_pages/ subfolder (page flood).
- sensitivities/run_dwd*_multi.m x6: PRESERVED as same-name,
  same-CONFIG thin wrappers over run_sensitivities (Dave: people use
  them -- do NOT delete); examples/run_* dirs = thin drivers too.
- Tests: tRunSensitivities (full-rank + poke-localization regression
  gates on the SMM trap corpus) + tRunSegmentation; both in
  SUITE_FAST.

SMM EDGE-SENSOR REWORK (Dave's spec, landed + validated):
- Per SHARED EDGE: 2 sensor locations at +/-SensorOff (new prompt,
  default 0.25*width; new stdin answer BEFORE the final Y --
  fixtures + segmirmaker_run 'sensor_off' updated) x 3 axes:
  1 piston = surface normal, 2 gap = in-plane PERP to edge, 3 shear =
  in-plane ALONG edge.  NO absolute-piston anchor row (not a
  measurement -- Dave).  Pie = 72 rows = 24 locations x 3 axes; hex2
  config-2 = 252 (42 edges x 6).  Prescriptions BYTE-IDENTICAL.
- QUIRKS FIXED (flagged for Dave): legacy rows reused rhoi for BOTH
  segments + projected triad-before-cross; new rows use each
  segment's OWN arm: dm/d(del_s)=+/-a'T_s, dm/d(th_s)=+/-(rho_s x
  a)'T_s.
- Hx.m now carries MeasAxis/MeasLoc/SensorPos; edge_sensors ingests
  (axis/loc/sensor_pos/has_anchor; legacy files still parse).
- tEdgeSensors REWRITTEN (axis recovery both sides, per-segment
  sensor-point coincidence, on-edge at +/-SensorOff; rows determine
  the point only PERP to the axis -- all checks project the axis
  out; in-plane axes are point-local).  tSegMirMaker refs
  regenerated.  3/0 + 3/0.
- s5 RERUN with 72-row dedx: as-built edge+MET WEM 10.81 -> 3.918,
  optimized 4.59 -> 1.777 (in-plane DOFs now edge-observable);
  MET-only optimized 3.766 unchanged (truss stands alone).  dldx FD
  1.1e-7 PASS, merit FD 0.00%.  (Hx radhat/tanhat->gap/shear relabel
  afterward leaves WEM invariant -- same rowspace, isotropic noise;
  e2e sidecar regenerated with final semantics.)

DAVE'S s6/s7 SPEC (VERBATIM INTENT, next to build):
- Linear model: w = dwdx*x + dwdz*z + dwdgrid*grid + dwdu*u + w0;
  u = control = SEGMENT + SM rigid-body DOFs (dwdu = those dwdx
  columns).  Measurement m = dmdx*x + m0, m = [l; e], dmdx = [dldx;
  dedx].
- s6 = run_compare: poke each DOF in turn (default 100 nm / 100
  nrad), display 2 graphics -- mmacos vs linear -- EACH = OPD map
  above stacked bar charts of l, e_piston, e_gap, e_shear; settable
  dwell (default 0.25 s); + saved frames + per-poke agreement report.
  Engine side: reuse dw_dx/met_calc paths for w/l; e engine-truth
  from finite rotations at es.sensor_pos.  s7 = run_simulator
  (estimator/controller, uncontrolled + controlled).

COMMIT GATE IN FLIGHT: s3 final Hx regen + fast suite running; on
green commit BOTH repos (macos: SegMirMaker untouched -- SMM lives in
MACOS_resources).  NO PUSH until Dave asks.  DEFERRED: hole->SM
chain rerun; run_design (s1/s2 recast); deletions list (inventory in
hand: freeform_unobscured retired-candidate, examples/design/oatma
orphan, e5hex1 old runners -- Dave must sign off; NEVER
coro_planet_demo); GMI mex relink; pymacos api exports.

**2026-07-19: RUNNERS DOCTRINE (Dave) + s5 BUILT AS THE GENERAL
run_met RUNNER.**  Dave: the PRODUCT is a small set of reusable stage
runners handing the .in from stage to stage — design → segmentation →
sensitivities → met → **compare** (NEW stage: step x/z/grid states,
emit w/e/l from mmacos AND the linear model) → simulate (needs an
ESTIMATOR/CONTROLLER; uncontrolled + controlled outputs).  Decision:
build s5/s6 as general runners from day one; recast s1–s4 into
run_design/run_segmentation/run_sensitivities in ONE pass after s6
(interfaces settle when the last consumer exists), then rethink the
mmacos file structure for new users + DELETE obsolete code/examples.
Also queued: run_dwdgrid_multi should adopt the ray-history footprint
approach.  Landed this session (mmacos, LOCAL until suites green):
- `design/runners/` (on the mmacos_setup path) + README (pipeline
  table, handoff contract = .in + declared sidecars) + **run_met.m**:
  as-built add_met → trace parity + aft-beam HOLE-margin table → dedx
  (Hx→SI: rot cols ×cbm ONLY) → dldx engine-FD vs analytic GATE →
  merit table → gains dxde/dxdl/dwde/dwdl → met_layout_opt →
  realized winner + engine-FD validation → MC acceptance → report +
  met_view/view_rx/metric figures + <name>_met.mat.
- `macos.design.met_layout_opt` = e5_seg_metopt v3 hoisted, SHAPE-
  CLASS aware: classes by boundary CONGRUENCE (vert count + edge
  lengths of seg_boundary polygons — pie → hexagon+wedge, hex → one,
  rxpoly imports classify free), one pattern per class about each
  member's own frame-x (pattern_frame 'segment'; 'radial' = e5
  heritage), coordinate-descent sweeps + top-K cross refinement.
  e5_seg_metopt.m is now a thin consumer — regen reproduces the v3
  winner EXACTLY (3.421 nm, FD 0.00%, same cluster angles; pmap/fclock
  differ by the fiducial-ring rotation symmetry; v3 baseline was
  ALREADY Inf — infeasible 10 mm corner sep, not a regression).
- `macos.design.seg_from_rx`: rehydrates the seg struct from the
  segmented .in ALONE (elt-type scan + met_bodies triads + lMon +
  src_seg_get tiling) so runners take files, not in-memory structs.
- `dldx_analytic` grew optional `unit_to_m` (default 1e-3 mm
  heritage): **e2e BaseUnits = METRES — the hard-coded mm arm scale
  would have squashed rotation rows 1000×.**  e2e Hx needs NO scaling
  (readings already metres); general rule in run_met: rot cols ×cbm.
- e2e s5 driver `s5_met.m` (thin): hub=M2 (elt 8), **aft ring on the
  FOLD FM (elt 9), NOT M3 — M3 is 2.1 m off-axis, its beams cross the
  M1 plane at ~2 m radius; FM is the on-axis bench body whose 0.10 m
  ring reaches the M2 rim fiducials THROUGH the 0.1225 m hole (run
  verified: crossing radii 0.079–0.10, margin 22 mm CLEAR)**.  54
  DOFs = 7 segs + M2 + FM.  Verified before the OOM (below): trace
  parity exact, FD gate 1.1e-7 PASS, 54 cols found.
- **OOM lesson (the two "VS Code crashed" reports = THIS)**:
  run_met's first merit form trace(Dk·P·Dk') with Dk 60412×54 built a
  29 GB dense → killed the 30 GB box twice.  Fixed to trace(P·G),
  G=D'D (memory feedback_trace_gram_not_outer); long MATLAB runs now
  caged via systemd-run MemoryMax.
- Tests: tMet +test_dldx_analytic_unit_to_m; NEW tRunMet (seg_from_rx
  rehydration identity; shape-class discovery [1 6] + WEDGE-PATTERN
  CONGRUENCE in segment frames; run_met end-to-end on the e5 pie
  fixture with synthetic jac + trimmed grids).  tRunMet added to
  SUITE_FAST.
- IN FLIGHT at session edge: caged s5 e2e run + tMet/tRunMet batch;
  then fast suite → commit both repos.

**DAVE'S s5 DESIGN REVIEW (2026-07-19, mid-session — supersedes the
first s5 numbers):** the first optimized layout (tight ±5° clusters,
WEM 4.9) "scores far worse than others we have evaluated" — it
maximizes neither beam angles nor outer-edge coverage.  Directives,
ALL IMPLEMENTED in run_met/met_layout_opt (rerun pending):
1. HEADLINE METRICS DIMENSIONLESS (WEM = wavefront per unit gauge
   noise, suite ratios held; floor as fraction of prior) — absolute
   nm at arbitrary priors "invites panic" (HWO context).  MET tracks
   the post-WFSC DRIFT state (≪1 nm class), gauge roadmap ~1 pm →
   e2e scenario sigmas now 1e-10 prior / 1e-11 edge / 1e-12 met.
2. Report WEM full / TILT / NON-TILT for MET-only AND edge+MET (tilt
   control absorbs tilt; orthogonal split via per-field piston+tilt
   projection of the s4 per-field blocks — Gram-only, no row order).
3. OPTIMIZE THE TRUSS WITHOUT EDGE SENSORS (must stand alone) on the
   NON-TILT Gram; edge sensors stay in the reported cases + estimator
   products.  Edge-sensor PLACEMENT is never optimized — SMM Hx as-is
   (the figures' open gray circles were the as-built LAUNCHER overlay;
   real sensor midpoints now drawn as gray dots via met_view
   'sensor_pts').
4. MANUAL BENCHMARK (auto 'corner_pairs' extra in met_layout_opt):
   launcher pairs at each class's two outermost boundary corners + a
   pair on the inside edge; evaluated GATE-BYPASSED with its true
   min-sep reported (adjacent-wedge corner pairs sit ~gap apart vs
   the 50 mm rule — a design datum, not hidden).
5. AFT STRUCTURE TO THE SM RADIUS (0.232 m ring on the fold/bench;
   'hole_r' override = SM shadow radius for the clearance check;
   "the hole will be that large anyway").  QUEUED (Dave's call):
   enlarge the REAL Rx hole via P.tel → full s1→s5 chain rerun.
6. AFT LEG SOLVED FIRST (Dave: "solve the simpler M3-SM truss first,
   fix it"): the aft ring's clocking + ITS OWN fiducial map = the
   first coordinate block in the descent, frozen for the class
   sweeps, revisited once.  add_met grew extra_clock/extra_pair_map.
7. PRESET LAYER (Dave: "save the optimized configuration as a new
   as-built for future PIE builds"): met_layout_opt exports every
   winner SCALE-FREE (class pattern angles, fiducial RIM INSET, aft
   block) + 'apply' mode realizes a preset on any same-tiling build;
   run_met 'preset' option makes it the as-built; s5 driver promotes
   to design/runners/presets/pie_met.mat; tRunMet round-trip test.
8. SYMMETRIC (ROTATIONAL) FIDUCIAL ASSIGNMENT + nf=6 (Dave: "with 6
   [fiducials] the segments can be ~interchangeable"): sym_assign
   (default ON) shifts each member's fiducial map by its clocking →
   congruent beam geometry per member; add_met accepts nseg×6
   pair_map; e2e forces nf_grid=6 (the nf=3 raw-merit winner broke
   the 60° symmetry).
**SESSION CLOSE-OUT (Dave 2026-07-19): pushed (MACOS_res 51c2f77 /
macos 5669ffe).  Decisions: (a) hole→SM-radius CHAIN APPROVED but
DEFERRED — P.hole_min_r_m floor is COMMITTED in s1/s2 [5], artifacts
NOT yet regenerated; run the chain (one MATLAB per stage — the §0
model-transition bug punishes mixed-size single processes; script
pattern in the transcript) AFTER the s3/s4 cleanup.  (b) NO corner
hardware sharing (corner_pairs stays a benchmark).  NEXT SESSION =
RECAST: s1–s4 → run_design/run_segmentation/run_sensitivities +
sensitivities/examples/* onto robust runners + file-structure rethink
+ obsolete deletion + run_dwdgrid_multi ray-history; AND investigate
Dave's flag: "improve the pie sensitivity results — they don't
compare well to the e5 examples" (pie dwdx/dwdz/dwdgrid quality vs
the e5 fixtures — start by diffing s4's spectra/conditioning and maps
against run_dwdx/dwdz/dwdgrid_multi outputs on e5).  Then s6
(run_compare first, then run_simulator w/ estimator/controller).**

FINAL RESULTS (s5 v3, nf=6 + sym assignment + solved aft block):
as-built WEM 15.6 MET-only / 10.8 edge+MET → **optimized 3.77 /
4.59, floor 0.28% of prior** — the SYMMETRIC nf=6 config BEATS the
asymmetric nf=3 raw-merit winner (4.86) outright: interchangeability
AND merit.  Dave's corner_pairs: 3.97 (ties the machine) but min-sep
15 mm violates the 50 mm rule at wedge-corner junctions.  WEM_tilt ≡
0 by CONVENTION: s4 dw_dx channels are fp_mode='track' harvested →
global tilt out at the source (prior tilt split 8e-13 nm);
per-segment tilt fully in the merit.  MC 0.7%, FD 0.00%.  Preset
promoted: design/runners/presets/pie_met.mat.  e5 regen under the
new optimizer: **3.421 → 3.255 nm, worst-mode 184 → 171 nm, FD
0.00%** (aft block + sym assignment — improvement, documented).

**s5 DESIGN CONSTRAINT (Dave 2026-07-18): MET-configuration
optimization solves ONE launcher pattern per segment SHAPE CLASS,
expressed in the segment frame, replicated to all same-shape segments
(pie: hexagon class + wedge class; hex2: one hexagon class) — the
tier-3 symmetry-first collapse, now shape-class-aware.  add_met
'launch_pts' + seg_boundary/rxpoly polygons are the realization
hooks.**

NEXT: s4_jacobians.m — dwdx/dwdz/dwdgrid on e2e_pie.in (dw_d*_multi
harvests; per-segment 6-DOF x ordering per Sprint-2D), then s5 MET
(shape-class patterns per above), s6 simulator.  Views + THOROUGH
report every stage.
RESUME RECIPE (post-compaction): (1) re-read root+nested CLAUDE.md,
MEMORY.md (project_e2e_example), this file; (2) artifacts live in
mmacos/design/examples/e2e/ (s1_telescope.in/.mat = DL telescope,
s2_instrument.in/.mat = Offner system, reports+PNGs beside them;
e2e_params.m = every knob; README = 12-rule design procedure);
(3) s3_segmentation.m = NEW runner: segment M1 of the e2e system
(segment_rx parent = s2_instrument.in, elt 1, P.seg block: Hex,
rings 1, gap 25 mm, emit_apertures, model 512 -- the e5_seg pattern,
watch the FF-parent + ApStop notes in project_sprint2d_segmentation);
then s4 dwdx/dwdz/dwdgrid (dw_d*_multi), s5 MET (add_met + metopt v3
+ met_bodies), s6 simulator (COMPOSE/psf + linear model switch).
Views + THOROUGH report at every stage (Dave's standing rule).**

**(superseded record of session 2 below)**
**2026-07-17 (e2e session): STAGE 1 LANDED (MACOS_res `e20d6b1`,
LOCAL).  `design/examples/e2e/` = the 6-stage worked example; all
knobs in `e2e_params.m`.  Dave's decisions this session: (a) user
specifies BASIC telescope parameters, f/# a FREE input (tma_layout,
NOT D-scaling); (b) the case = D=4 m, primary f/1.25, system f/18,
ON-AXIS Korsch taken SLIGHTLY OFF-AXIS; (c) 90-deg FOLD after M2
moves M3 + image + FP BEHIND M1; (d) every stage runner emits a
THOROUGH design report + graphics; (e) stage 2 = JOINT refinement:
as instrument optics are added, keep refining M1-M3 (M2/M3 apertures
may GROW, Zernike terms deepen) to improve field performance.
Stage-1 as-built: m2=16 -> f/20 int focus in front of M1 (met
injection), near-unit M3 relay, ~6% M2 obscuration, fold at z=0.3,
bench M3 [2.1 0 0.3] / FP [-0.5 -0.12 0.3], M1 hole r=0.182, bias
sweep -> 5' (least that fully clears; only M2 obstructs, by design);
solve ladder (each step was a debugged lesson): conics AT the bias
point with FP align FIRST (align between two solves; else -2 mm
field-curvature defocus = 1.6 waves poisons the conic solve) -> joint
FF field solve (1.07 worst) -> M1 stop common-mode null (bias point
0.0024 waves, Strehl 0.97; corners pay: worst +-1' = 1.58 raw/0.97
-tilt = the pure field DIFFERENTIAL stage 2 corrects).  ORDER
MATTERS: static-null-then-joint lands WORSE (1.48) than
joint-then-null -- LM basin path dependence.  NEW
`design/src/field_zone_lmon.m` (doctrine field-zone lMon; tested in
tDesignTelescope).  Reload 1305/1305; fast suite 225/0.
NEXT: s2_instrument (3-4 mirror relay widening toward +-2', joint
M1-M3 re-solve per (e)), then s3..s6 per the sequence below.**

**2026-07-12: PLAN_DESIGN_LAYER Sprint 2D — SEGMENTATION + SENSING
(Dave's objective): segment the design flow incl. per-segment local
coordinate systems, edge sensors, laser metrology; MET truss from M2
to points around M3; future = optimize the MET configuration.**

### The frame (Dave, 2026-07-12 — supersedes §6.6 tier-3 draft merit)
Deliverable = the **linear forward model of a segmented system on one
shared per-segment DOF vector x**, for closed-loop control analysis:
- optics `w = dwdx·x + dwdz·z + dwddm·dm + w0` (channels exist)
- edge sensors `e = dedx·x + e0` (dedx = SegMirMaker `Hx`)
- laser MET `l = dldx·x + l0` (poke → METcalc → FD; NEW channel)
- simple control `dx = −pinv(dwdx)·(w − wtarg)`, w from the ESTIMATE.
**Configuration merit = post-control wavefront residual.**  With
full actuation the correctable part cancels exactly and
`w_post = dwdx·(x − x̂)`, so
`merit = E‖w_post‖² = trace(dwdx · P_δx · dwdxᵀ)`,
`P_δx = X − X·Hᵀ(H·X·Hᵀ + R)⁻¹·H·X`, `H = [dedx; dldx]`, X = prior
cov of x (deploy tolerances), R = sensor noise.  Prior wavefront cov
`W = dwdx·X·dwdxᵀ` is the baseline to beat (report the ratio).  When
actuation ⊄ states, report the uncorrectable projection
`(I − dwdx·pinv(dwdx))·w` separately.  Unobservable directions
saturate at their prior wavefront cost (no singular inverse).  MET
placement optimization (tier 3, future) minimizes this merit over
launcher/fiducial placement + beam topology — pure MATLAB on the S4
outputs.

### Decisions (Dave, 2026-07-12)
1. Dev/test parent = **e5mono import** (SegMirMaker's canonical parent,
   committed references); showcase on the design-layer 3M after.
2. **Hx normal-height edge-sensor model accepted** (relative
   surface-normal displacement at edge midpoints — piston+dihedral, no
   in-plane gap/shear); richer sensor models later = pluggable dedx
   backends behind the same Jacobian contract.
3. **MET default geometry = Stewart-platform trusses**: 6 launchers on
   EACH segment illuminating 3–6 fiducials on M2; 6 launchers around
   M3 illuminating 3–6 fiducials on M2; one measurement per
   launcher/fiducial pair = change in straight-line distance;
   ≥ as many measurements as DOFs.

### Slices
- [x] **S0 — SegMirMaker refresh + batch driver (Q7).  DONE 2026-07-12,
  commit pending fast-suite gate.**  CMakeLists repointed
  build_release_giza→build_release, npsol/lapack/blas→slsqplib;
  **`-fp-model strict` added → a rebuilt binary reproduces the
  committed test_in references BYTE-IDENTICALLY** (first non-strict
  build differed at 1 ulp; strict closed it; verified 2× fresh-dir
  sha256).  Batch mode = scripted stdin (NO control-file mode exists;
  README's "nine questions" is really ~15–17 prompts; DOF default is
  3 not the documented 6 — README fixed).  Committed answer fixtures
  `test_in/e5pie.stdin`/`e5hex2.stdin`.  MATLAB driver
  `macos.design.segmirmaker_run` (fresh scratch dir; copies parent +
  macos_param.txt + GridFile= refs; positional answers; 'Done.' gate)
  + `tSegMirMaker` (3/3: pie + hex2 byte-identity, size-args error)
  wired into SUITE_FAST.  GOTCHA: rerunning in a used dir trips the
  overwrite prompt and shifts the answer stream — always fresh dir.
- [x] **S1 — `segment_rx` splice.  DONE 2026-07-12** (`macos.design.
  segment_rx` + `tSegmentRx` 4/4).  As-built: engine numbers elements
  by READ ORDER (`iElt=`/stale `nElt=` are cosmetic — e5hex2.in loads
  25 blocks despite nElt=24); splice replaces the parent's source
  GridType with the segment tiling + inserts nSeg/width/gap/SegXgrid/
  SegCoord after yGrid; nElt from BLOCK COUNT; downstream iElt
  renumbered cosmetically.  **ApStop never renumbers** (header form =
  3-vector StopPos POSITION msmacosio.inc:90; element form = 2-vector
  offset inside the stop element's own block :2887, last-wins) — but
  if the SEGMENTED element carried the element-form stop it drops
  with the block → warning + out.dropped_apstop.  OptTgtElt/RefElt/
  tMetElt refused pre-run.  Validation = load + trace parity (seg
  0 lost rays, rms 3.564e-5 vs parent 3.920e-5 mm = sampling) +
  frames orthonormal.  **FF-parent gotcha: ring centers are NOT
  3-D-equidistant (mm-scale astig figure) — tiling-plane invariants
  only.**  Original plan line: [ ] Telescope.segment() builder method
  — deferred to the 3M showcase.
  splice .presc blocks in place of the parent elt + downstream
  renumbering (emitter owns numbering; watch ApStop/element-index
  refs) → reload → validate (ValidatePrescription + segment-center
  ray check ≈ SEGRAYTRACE's: center-ray endpoint vs RptElt, tol
  ~1e-11).  Spec carries per-segment frames; `segment_frames()`
  readback.  Segment local csys = RptElt + TElt/pMon,xMon,yMon,zMon
  triad (SEGMENT = full parent surface: same conic/FF/psi/Vpt every
  segment; Mon slot zeroed = reserved per-segment figure channel;
  DOF columns = [rot·x̂ rot·ŷ rot·ẑ | trans·x̂ ŷ ẑ] in the triad).
- [x] **S2 — edge sensors (dedx).  DONE 2026-07-12** (`macos.design.
  edge_sensors` + `tEdgeSensors` 3/3).  Ingests Hx.m via isolated-
  workspace run(); pads row-sparse assignments to nMeas×nState; per-
  segment columns [rot_xyz | trans_xyz] IN THAT SEGMENT'S TRIAD;
  row 1 = master piston.  Validation WITHOUT surface re-eval, from
  the generator's algebra (SegMirMaker.f:588-645): translation
  triplet = normhatᵀ·T → unit norm + T_i·del_i' == −T_j·del_j'
  (same world normal both ways); rotation rows = shared-ρ cross
  terms (the generator uses rhoi for BOTH segments) → [del]ₓρ=th'
  solves ρ = pSeg_i − pr, pr lands on the shared-edge midpoint
  (lateral < 0.05 width; residual limited by Hx TEXT precision
  ~5e-10 — tolerances 1e-8/1e-6, not 1e-12).  Frame-mixing quirk
  noted: Fortran multiplies triad-component rows by world [ρ]ₓ —
  mirror the code, don't "fix" it silently (model review with Dave
  if it matters downstream).
- [x] **S3 — DONE 2026-07-12 (engine leg + add_met Stewart emitter, tMet 4/4; 48-beam truss == geometry to 1e-12)** (`met_calc`/`met_get` in
  macos_api_mod after ray_status_get; capacity mMetSrf 20→64 /
  mSysMetBeam 128→512 in elt_mod; codegen Path A → `macos.met()`
  veneer with SI/native units; `tMet` 3/3 = Q8 closed-form gate:
  exact baseline lengths on a hand-inserted m2→fpa 2-beam fixture,
  gauge Δ == −û·d LOS projection under global perturb (5e-6), ⊥-null
  < 1e-7).  Both engine trees + mex rebuilt.  **Build gotcha: makems
  must run from ~/dev/macos (repo root) — a bad cwd 'succeeds' via
  stale logs; mmacos make needs explicit
  MACOS_BUILD_DIR=~/dev/macos/build_release_gfortran with
  FC=gfortran (its default is still build_release_giza — fix
  sometime).**  REMAINING: builder `add_met(...)`:
  Stewart trusses per decision 3 (fiducial/launcher placement in the
  segment triad → global; emits nMetPos/SrfMetPos rows + tMetElt +
  metBeamFlg).  Engine (small): `met_calc` + `met_get` wrappers in
  macos_api_mod (codegen Path A → mmacos veneer `macos.met()`), and
  CAPACITY BUMP `elt_mod.F:332-336` mMetSrf 20→64, mSysMetBeam
  128→512 (19 segs + M2 + M3 = 21 surfaces; 19·6+6 = 120 beams —
  2-ring barely fits today, 3-ring doesn't).  Ship with Q8-style
  closed-form tests (gauge Δ == LOS projection of relative fiducial
  motion; orthogonal-motion null; reciprocity) = the §6.6 tier-2 gate
  for SrfMetCalc.  Document limitation: straight-line lengths, NO
  LOS/obscuration check (PLAN §4.5 ShowMetObscur deferred).
- [ ] **S4 — SHAPE PER DAVE 2026-07-12: `design/examples/e5_seg/` — a MODIFIABLE runner: e5mono.in -> segmented .in INCLUDING the MET points (segment_rx + add_met) -> dedx (edge_sensors) + dldx (dmet_dx FD channel) -> MET metric performance = trace(dwdx*P_dx*dwdx') vs prior W=dwdx*X*dwdx' ratio.  Example rules: save .in+.mat, figures in dir, README, NO exit(0).  Underlying pieces: `dmet_dx` channel + `forward_model()` + demo.**  FD
  Jacobian via perturb→METcalc→read (dw_dx supervisor pattern; the
  API perturb path CPERTURB_PROG DOES move SrfMetPos).
  `forward_model()` = {dwdx, dedx, dldx, w0,e0,l0} on ONE x ordering
  (per-segment 6-DOF local + M2/M3 rigid).  Worked example: segmented
  e5-class primary → forward model → draw x~X, measure with noise,
  estimate, control dx=−pinv(dwdx)·ŵ → simulated RMS w_post ==
  analytic trace(dwdx·P_δx·dwdxᵀ).  That demo = sprint acceptance.

### Engine facts (2026-07-12 surveys; file:line current that day)
- MET: `nMetPos`/`tMetElt`/`metBeamFlg` parsed msmacosio.inc:1937-1975
  (nMetPos MUST precede tMetElt or parser STOPs); global-frame points
  `SrfMetPos(3,48,20)`; ONE command `METcalc` (3-char min-match,
  macos_cmd_loop.inc:2511) → `SrfMetCalc` utilsub.F:1810 = **pure
  straight-line Euclidean distance** source-point→target-point;
  flat output `metMeasBuf(1:nMetMeas)` (GMI output #9, pflg(27)).
  No LoadStack entry, no api_mod wrapper, no mmacos surface yet.
- SrfMetPos moves under CPERTURB + CPERTURB_GRP + **CPERTURB_PROG
  (funcsub.F:372-380 — the API/mmacos/pymacos path, so FD works)**;
  CPRead / CPERTURB_2 / LnkEltCPERTURB do NOT (PLAN §4.5 gaps — avoid
  those paths for met work).
- Engine `EdgeSensors` keyword = parsed+stored+SAVEd but **dead** (no
  consumer; PLAN §4.5 recommends not building on it) — dedx comes
  from SegMirMaker Hx per §6.6 tier 1.  SAVE round-trips all met
  keys since 662e86e (nMetPos precedes tMetElt in the writer).
- SegMirMaker: standalone fixed-form Fortran, own CMake vs
  build_release; loads the parent through the REAL engine loader
  (needs macos_param.txt + GridFile data in cwd);
  `SEGRAYTRACE` = engine-side QA command (macos_cmd_loop.inc:1391).

### Tier-3 MET-layout optimization plan (Dave's Q 2026-07-12, answered)
Merit evaluation is ANALYTIC once frames are known: dldx rows are
closed-form LOS projections + moment arms (validated FD==analytic in
tMet), dedx fixed, dwdx fixed per optical design -> thousands of
layouts/sec of pure linear algebra, NO engine in the loop (engine FD
only validates the winner).  Approach = HIERARCHICAL, not raw
combinatorics: (1) SYMMETRY first -- one launcher pattern per segment,
mirror-symmetric about the segment centerline, replicated by the hex
rotation group -> collapses the discrete space to one segment's
pattern x global ring params; (2) COMBINATORIC enumeration of that
small symmetric set (affordable analytically); (3) STEP-AND-EVALUATE /
patternsearch polish of the continuous knobs (radii, angles, standoff,
nf) on the shortlist; (4) worst-mode check (max eig of P_w, not just
trace) + a symmetry-BREAKING perturbation stage if a symmetric blind
mode saturates; (5) engine-FD validation of the final layout.  BUILT 2026-07-12: dldx_analytic + e5_seg_metopt.m -- 20160 layouts/12.8 s, as-built 5.156 -> 3.732 nm, engine validation 0.00% off analytic.  Queued: hub/extra bodies in the analytic rows need engine TElt/RptElt frame parse; fid ring on M2 BACK side; beam-clearance check.

## In-session state NOT yet committed
- MACOS_resources: segmirmaker CMakeLists/makesegmirmaker.sh/README/
  CLAUDE.md edits + test_in/*.stdin + mmacos segmirmaker_run.m +
  tSegMirMaker.m + run_mmacos_tests.sh SUITE_FAST entry.  Fast suite
  running as the pre-commit gate → commit on green.
- macos: this CURRENT_SLICE rewrite only (no engine change yet;
  engine work starts in S3).

## Just tried / ruled out (with why)
- Byte-parity chase vs the retired giza-era binary — RESOLVED, not
  ruled out: `-fp-model strict` reproduces the committed references
  exactly; no reference refresh needed.
- Reusing a scratch dir for repeat runs — ruled out (overwrite prompt
  shifts the positional answer stream; driver always makes fresh).

## Parked threads (carry, do not lose)
- **eac5 thread (Dave 2026-07-07):** reflective e5 back end =
  `~/dev/tst_dir/eac5mono.in` m3 (R=0.328 m, lMon 20 mm) + m4
  (R=0.407 m, lMon 30 mm) compact mini-relay ~0.2 m past the f/21
  intermediate focus — the reference for the freeform_unobscured
  "+2" (supersedes single-M4).  File breakage: CRLF + tabs-in-TElt
  (l.146) + `/* */` comment blocks (l.151+); Phase-1 validator
  rejects ("blank line inside TElt block").  NEXT: does the parser
  know `/*`? clean a copy, load/trace, extract m3/m4 conjugates.
- **±1′ design fork (Dave's call pending):** (a) ANSI-45 mode-depth
  test via one Jacobian assembly + column-space projection (~15 min);
  (b) pupil-aware re-layout (Sprint 5); (c) re-scope field/λ.  The
  shipped freeform_unobscured 3M numbers rest on pathological
  surfaces (see [[project_zern_solve_doctrine]]); re-solve under the
  doctrine needs Dave's OK.  Probe scripts:
  `~/dev/MACOS_sandbox/freeform_4m_probes/`.
- Others' actions on the AppScan response: Luis closes #61 + deletes
  ft-dev (both repos); IT rescans + dispositions the 2 conftest FPs.

### Tier-3 AS LANDED (e06c254) + NEXT SESSION (Dave 2026-07-12)
Constrained optimizer SHIPPED: launchers ON the segment hex boundary
offset OUTWARD by EDGE_OFF=5 mm (clearance off the reflecting surface,
by construction; sign-flip for inside-the-rim), 3 MIRROR PAIRS about
each segment's radial centerline, free pair angles (center seg = xhat
line).  add_met 'launch_pts' override realizes any winner.  Hex run:
54400 layouts/34.7 s; edge ring 3.857 -> 2.973 nm rms (worst 215.8 ->
155.9); angular spread [30 90 150] CONFIRMED optimal by search; the
gain came from the HUB FIDUCIAL geometry (rfid 300->150, fclock 105);
engine-FD validation 0.00%.  e5_seg now GRID='Hex'.
**NEXT SESSION (Dave): (1) annotation + graphics for the e5_seg
example (layout views of the truss/launchers/fiducials, metric
figures); (2) more examples.  Queued: hub/extra bodies in analytic
rows (TElt/RptElt parse); fiducials on M2 back side; beam clearance;
full fast suite ran only per-class today (Dave deferred); PLAN
promote-on-land pass.**

### RESUME RECIPE (written for a fresh session/model — Opus next)
1. Read THIS file + memory `project_sprint2d_segmentation.md` +
   `mmacos/CLAUDE.md` (+ root CLAUDE.md per directive).  All Sprint-2D
   code is PUSHED (macos f740cf6 / MACOS_res e06c254 tips).
2. Regenerate working state (~8 min):
   `cd ~/dev/MACOS_resources/mmacos && matlab -batch "run('mmacos_setup.m'); run('design/examples/e5_seg/e5_seg.m'); run('design/examples/e5_seg/e5_seg_metopt.m'); exit(0)"`
   Expect: edge+MET 5.978 nm, MC ~2.3%; optimizer 3.857→2.973 nm,
   engine validation 0.00%.  Tests: `./run_mmacos_tests.sh tMet` (and
   tSegMirMaker/tSegmentRx/tEdgeSensors), 19 green.
3. NEXT TASK (Dave): annotation + graphics for e5_seg, then more
   examples.  Graphics hooks that already exist: `macos.draw_rays` /
   `Telescope.view_layout` (real ray bundle), `am.src_pts/tgt_pts`
   (3xN global truss endpoints — plot3 beams launcher→fiducial),
   `seg.frames` (segment centers/triads for hex outlines + labels),
   `e5_seg.mat`/`e5_seg_metopt.mat` (all Jacobians + metric tables).
   Figures land in the example dir (house rule; no exit(0) in
   examples — the -batch wrapper above supplies it).
4. TRAPS a fresh model WILL hit: macos.modify() after every poke
   (cached OPD); macos.opd() = no args, N×N, WaveUnits; run makems
   from ~/dev/macos ROOT; mmacos mex: make FC=gfortran
   MACOS_BUILD_DIR=~/dev/macos/build_release_gfortran; SegMirMaker
   scratch dirs must be FRESH (overwrite prompt shifts stdin answers);
   engine numbers elements by read order; ad-hoc triads for hub/fpa
   DON'T match engine TElt/RptElt (why analytic rows are seg-only).

## Next concrete step
Fast suite green → commit S0 (MACOS_resources + this file on macos)
→ S1 `segment()` splice: decide splice host (Telescope method vs
System.from_rx path — likely both share one private splice helper),
renumbering rules for downstream element-index keywords, and the
segment-center ray check via mmacos trace.

## Open micro-questions (slice-local)
- S1: which element-index-bearing keywords must renumber on splice
  besides nElt (ApStop index, tMetElt targets, OptTgtElt, RefElt,
  LnkElt/group refs?) — enumerate from msmacosio.inc when building.
- S3: met-point placement API — accept explicit point lists AND the
  Stewart preset; where do M2-fiducial defaults sit (radius on M2
  face)?
- S4: x ordering — segments first (6/seg, triad order rot|trans) then
  M2/M3 rigid?  Must match dw_dx channel DOF naming.

## Promote-on-land  →  then CLEAR this file
> Same commit as the `design-sprint-N` tag: move each item to its
> permanent home, then reset this file to the empty template below.
- [ ] PLAN_DESIGN_LAYER Sprint 2D checkboxes ticked (Q7, segment(),
      post-splice validation) + NEW: S2/S3/S4 sensing items recorded
- [ ] `CORE COMPLETE <date>` blockquote added to Sprint 2D
- [ ] §10 Decisions: Stewart MET default; post-control-residual merit
      (supersedes tier-3 min-σ(H)); Hx sensor scope — all dated
- [ ] PLAN.md §4.5: met wrapper/capacity items ticked when S3 lands
- [ ] CLAUDE.md / nested gotcha captured (segmirmaker batch gotchas →
      segmirmaker/CLAUDE.md — DONE in S0)
- [ ] agent MEMORY.md: sprint-2d project memory
- [ ] worked-example committed + named (S4 demo)
- [ ] **reset CURRENT_SLICE.md to empty template**

---

## Empty template (reset state — copy over the above on land)

```
## Active slice
- Sprint / item: —
- Plan anchor: —
- Branch / worktree: sls-dev + sls-dev @ —
- Definition of done (honest): —

### Tier-3 MET-layout optimization plan (Dave's Q 2026-07-12, answered)
Merit evaluation is ANALYTIC once frames are known: dldx rows are
closed-form LOS projections + moment arms (validated FD==analytic in
tMet), dedx fixed, dwdx fixed per optical design -> thousands of
layouts/sec of pure linear algebra, NO engine in the loop (engine FD
only validates the winner).  Approach = HIERARCHICAL, not raw
combinatorics: (1) SYMMETRY first -- one launcher pattern per segment,
mirror-symmetric about the segment centerline, replicated by the hex
rotation group -> collapses the discrete space to one segment's
pattern x global ring params; (2) COMBINATORIC enumeration of that
small symmetric set (affordable analytically); (3) STEP-AND-EVALUATE /
patternsearch polish of the continuous knobs (radii, angles, standoff,
nf) on the shortlist; (4) worst-mode check (max eig of P_w, not just
trace) + a symmetry-BREAKING perturbation stage if a symmetric blind
mode saturates; (5) engine-FD validation of the final layout.  BUILT 2026-07-12: dldx_analytic + e5_seg_metopt.m -- 20160 layouts/12.8 s, as-built 5.156 -> 3.732 nm, engine validation 0.00% off analytic.  Queued: hub/extra bodies in the analytic rows need engine TElt/RptElt frame parse; fid ring on M2 BACK side; beam-clearance check.

## In-session state NOT yet committed
—

## Just tried / ruled out (with why)
—

### Tier-3 AS LANDED (e06c254) + NEXT SESSION (Dave 2026-07-12)
Constrained optimizer SHIPPED: launchers ON the segment hex boundary
offset OUTWARD by EDGE_OFF=5 mm (clearance off the reflecting surface,
by construction; sign-flip for inside-the-rim), 3 MIRROR PAIRS about
each segment's radial centerline, free pair angles (center seg = xhat
line).  add_met 'launch_pts' override realizes any winner.  Hex run:
54400 layouts/34.7 s; edge ring 3.857 -> 2.973 nm rms (worst 215.8 ->
155.9); angular spread [30 90 150] CONFIRMED optimal by search; the
gain came from the HUB FIDUCIAL geometry (rfid 300->150, fclock 105);
engine-FD validation 0.00%.  e5_seg now GRID='Hex'.
**NEXT SESSION (Dave): (1) annotation + graphics for the e5_seg
example (layout views of the truss/launchers/fiducials, metric
figures); (2) more examples.  Queued: hub/extra bodies in analytic
rows (TElt/RptElt parse); fiducials on M2 back side; beam clearance;
full fast suite ran only per-class today (Dave deferred); PLAN
promote-on-land pass.**

### RESUME RECIPE (written for a fresh session/model — Opus next)
1. Read THIS file + memory `project_sprint2d_segmentation.md` +
   `mmacos/CLAUDE.md` (+ root CLAUDE.md per directive).  All Sprint-2D
   code is PUSHED (macos f740cf6 / MACOS_res e06c254 tips).
2. Regenerate working state (~8 min):
   `cd ~/dev/MACOS_resources/mmacos && matlab -batch "run('mmacos_setup.m'); run('design/examples/e5_seg/e5_seg.m'); run('design/examples/e5_seg/e5_seg_metopt.m'); exit(0)"`
   Expect: edge+MET 5.978 nm, MC ~2.3%; optimizer 3.857→2.973 nm,
   engine validation 0.00%.  Tests: `./run_mmacos_tests.sh tMet` (and
   tSegMirMaker/tSegmentRx/tEdgeSensors), 19 green.
3. NEXT TASK (Dave): annotation + graphics for e5_seg, then more
   examples.  Graphics hooks that already exist: `macos.draw_rays` /
   `Telescope.view_layout` (real ray bundle), `am.src_pts/tgt_pts`
   (3xN global truss endpoints — plot3 beams launcher→fiducial),
   `seg.frames` (segment centers/triads for hex outlines + labels),
   `e5_seg.mat`/`e5_seg_metopt.mat` (all Jacobians + metric tables).
   Figures land in the example dir (house rule; no exit(0) in
   examples — the -batch wrapper above supplies it).
4. TRAPS a fresh model WILL hit: macos.modify() after every poke
   (cached OPD); macos.opd() = no args, N×N, WaveUnits; run makems
   from ~/dev/macos ROOT; mmacos mex: make FC=gfortran
   MACOS_BUILD_DIR=~/dev/macos/build_release_gfortran; SegMirMaker
   scratch dirs must be FRESH (overwrite prompt shifts stdin answers);
   engine numbers elements by read order; ad-hoc triads for hub/fpa
   DON'T match engine TElt/RptElt (why analytic rows are seg-only).

## Next concrete step
—

## Open micro-questions (slice-local)
—

## Promote-on-land → then CLEAR this file
- [ ] PLAN checkbox(es)
- [ ] CORE COMPLETE blockquote
- [ ] §10 Decisions entry
- [ ] CLAUDE.md / nested gotcha (if any)
- [ ] agent MEMORY.md learning (if any)
- [ ] worked-example committed
- [ ] reset this file
```
