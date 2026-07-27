# Expose & extend MACOS polarization physics — IFO metrology + coronagraph modeling

> Revised from the Opus draft 2026-07-25 (CC review): engine claims re-verified against
> source, two factual corrections folded in (coating-wavelength semantics, Ex0Ey0
> keyword state), the **interferometer/bench application added as a co-equal track**
> (Dave: "add polarization to this model, AND to the coronagraph models"), binding
> mechanics aligned with the actual mmacos codegen path, and the validation ladder
> extended with an IFO self-measurement cross-check.

## Context — two applications, one core

**A. Interferometer / bench metrology** (`examples/design/bench_ifo*`, the VSG-class
testbeds). The rigs run 45° plate beamsplitters and folds. The 2026-07-25 AOI sweep
found scalar metrology AOI-neutral in the compensated design — **the entire case for
smaller AOI lives in polarization** (R_s/R_p diattenuation, coating retardance), and
today none of it is modeled. Polarization sets fringe-visibility loss, s/p phase
splitting (a PSI systematic), and enables the classic instrument architectures we
cannot currently build: PBS + quarter-wave-plate Twyman-Green, polarization
phase-shifting.

**B. Coronagraph contrast floor.** `CORONAGRAPH_DESIGN_RULES.md` §3/§9 flag
polarization aberration as *"a hard floor for 1e-10 contrast … does not respond to
scalar DM control."* The budget needs a **number** — the polarization-limited contrast
floor and its sensitivity to coating choice — not just maps.

MACOS already contains substantial polarization physics in the engine, but **none of it
is reachable from pymacos/mmacos**, so neither application can use it. This plan
exposes the existing physics first (Dave's directive), builds the shared Jones-pupil
machinery, then delivers per-application analyses, then (optionally) extends the engine
with polarizing elements and spatially-variable coatings.

**Deliverable framing.** Phases 1+2 are the committed, shippable core: *given a
prescription, return the Jones pupil, the polarization-aberration decomposition, the
IFO visibility/PSI-systematic numbers (track A) and the contrast floor (track B).*
Phases 3 and 4 add modeling reach, carry most of the schedule risk, and are
independent of each other — sequence them by need (PBS/waveplate IFO and VVC trades →
Phase 3 first; segmented-primary coating non-uniformity → Phase 4 first).

## Branching & execution sequencing (agreed with Dave 2026-07-25)

**Branches — split at the consumer boundary** (single-issue rules; stacked-branch
pattern per bench-builder/expose-beam):

1. **`pol-core`** — single-issue branch in **both repos** (macos: engine API +
   `Ex0Ey0` fix; MACOS_resources: mmacos/pymacos bindings + tests). Scope: Phases
   0–1 plus the application-neutral part of Phase 2 (Jones pupil, `pol_maps`,
   unitarity/2θ tests). **CUT 2026-07-25 and pushed:** macos `pol-core` off
   `expose-beam` (742ce29 — carries the `beam_set` template the plan cites; same
   commit as macos bench-builder), MACOS_resources `pol-core` off **`bench-builder`**
   (cea3d60 — Dave's call: Phase 1's track-A smoke test drives
   `macos.design.twyman_green`, which lives there, and bench-builder ⊇ expose-beam).
   Stacked PRs reduce as parents merge. Later phases (elements, spatial coatings)
   become their own follow-on single-issue branches/PRs; the phase boundaries are
   PR-sized on purpose.
2. **`pol-ifo`** — thin consumer branch **stacked on `bench-builder`**: Phase 2d, the
   Bench coating emission, `add_polarizer`/`add_waveplate` wiring, the
   `bench_ifo_pol` example. Needs both parents; lands after bench-builder and
   pol-core merge (stacked PR auto-reduces as parents land).

⚠ The public history rewrite landed ~2026-07-24: confirm the current base branch
name (`dev` vs `sls-dev`) before cutting branches, and cut only from a fresh
post-rewrite clone (never push from a stale pre-rewrite clone).

**Execution — split at the judgment boundary, not the phase boundary.
REVISED 2026-07-26** (Dave): CCMac is out of tokens for the month, and Fable
time is the scarce resource — **anything Opus-executable is handed to Opus**.
The Opus lane is defined by the *model*, not the machine: it runs on CCMac
when its tokens reset, or as an Opus session on this Linux box in the
meantime.  (The standing ifx-smoke gate can only run on the Linux box either
way — CCMac is gfortran-only.)

- **Opus lane — template-following with objective, machine-checkable gates.**
  Delegation is safe where passing the gate *is* the definition of done.
  Worklist, in priority order (each item's spec + gates are already written
  in this document; the executor should need no new decisions):
  1. **Phase 3a Tranche 1** (§3a.2 — all call sites line-mapped): K=1,3
     kernel loops, 3-plane `FFObscure`, validation-ladder tests 1/2/4/5.
     Includes implementing the two §3a.1 fixes **to spec** (assembly
     seed-once-then-update; `PFFPROP`→`FFPROP`×3 unification) — revised
     from "written here": the spec is now detailed enough that a Fable
     line-review of the finished diff is cheaper than Fable authorship.
  2. ~~**Phase 2b low-order expansion**~~ — **DONE 2026-07-26.**
     `macos.pol_zernike` / `pymacos.pol_zernike` (pure binding-layer,
     no engine change) + 3 gates each side + report §2.5 + cmdref
     NOTES + manual §5.  Result: the Al Cassegrain reduces to pure
     **polarization astigmatism** exactly as the literature form
     predicts — astig0 in s1, astig45 in s2, equal to 1.9e-7
     (discretization, verified by grid scaling), every forbidden term
     (piston/tilt/defocus/coma/trefoil/spherical AND the whole
     circular component) at ~1e-15 of it, the rho^4 companion
     sub-dominant at 2.6e-3, and the rho^2 radial law's extrapolated
     on-axis diattenuation vanishing to 1e-4 — which nothing in the
     fit arranges.  Ladder rung 5 moves to PARTLY DONE: form matches,
     a numeric regression against a specific published system still
     wants that system set up here.
     Two side-findings, both recorded: (a) the shared ANSI evaluator
     now lives in `+macos/private/ansi_zernike_eval.m` so influence-
     basis and aberration-report mode indices cannot drift; (b) a
     **units trap** in `tJonesPupil` — `coating` takes BaseUnits and
     the class's two fixtures differ (`Rx_Cass_FarField`=m, Bench fold
     rig=mm), so one shared Al thickness silently meant 200 um on the
     Cassegrain.  Split; the bindings now agree to 11 digits.
  3. ~~**Phase 2c contrast floor**~~ — **DONE 2026-07-27.**  Landed on the
     analyzer-at-detector design: `macos.pol_contrast_floor` + `tPolContrast`
     (256) + `tPolContrastCoro` (512) + report §5 + cmdref/manual, and §2c's
     plan text amended to the shipped signature.  All four gates pass —
     x-pol reduces to the scalar contrast curve at round-off (cross exactly
     0, curve 2.4e-15), Parseval on the split 1.3e-24, energy closure
     1.8e-16, floor reported by component.  Deliverable numbers and the
     Tranche-1 scope finding (the grid carries 0.84 / 0.57 of the ray-level
     cross-pol on `Rx_Coro`, and reports the coating sensitivity with the
     WRONG SIGN) are in §2c.  Two side-effects: the polval driver is now
     split into per-model-size batches (128 / 256 / 512) with a
     `merge_numbers.py` stage, since 2c cannot run at 128; and the coherency
     matrix's conjugation order needed a circular-input gate to catch.
     *History of how it got unblocked, kept for the record:*
     **UNBLOCKED 2026-07-27 pm — the s/p sign fix is LANDED (Fable lane):
     standard `+r_p` restored in BOTH `Reflector` branches (uncoated RP;
     coated innermost RP + per-layer RP1 — the Airy recursion propagates
     the flip exactly, so the coated fix is the same clean restoration).
     Decision + verification + the one packet correction (the scratch
     suite prediction came from a HALF-patched engine; the local-sp
     artifact assertions are KEPT — measured 247× on the corrected
     engine) are appended to `REVIEW_POL_SP_SIGN_2026-07-27.md`.
     **TAIL CLOSED 2026-07-27 pm (Opus)** — all three follow-on items
     landed: (1) the odd-mirror gates, which turned out to support more
     than the bound the spec asked for — on the perfect-conductor idiom
     the WHOLE single-reflection Jones has a closed form, and the engine
     matches it at median 2.1e-15 / max 5.9e-14, with AOI and azimuth
     taken from ray DIRECTIONS so no pupil-grid mapping is assumed
     (`tPolarization/test_odd_mirror_crosspol_{pec_analytic,rho2_law}`,
     mmacos `fc2e22e`; non-vacuity measured by rebuilding the flipped
     engine — 7 of 8 assertions fail, residual median 1.14e+02, radial
     slope 0.033); (2) `Refractor` normalized (macos `25c4386`) with the
     transmitted field verified BIT-IDENTICAL and the inconsistent flip
     built to prove that A/B non-vacuous (−3.2% power); (3) polval
     regenerated with a new §4 evidence section, GMI 6/6 bit-identical,
     full mmacos suite re-run.  **New engine finding, open:** coated and
     uncoated `Refractor` transmission use different amplitude
     normalizations (the coated branch omits the radiometric
     `sqrt(n2 cos02/(n1 cos01))` factor — measured coated/uncoated |Ex| =
     0.816442 = 1/sqrt(1.5) at normal incidence on an index-matched
     layer, i.e. a coated lens under-transmits ~18% in amplitude), and
     the coated Refractor branch has no analytic gate at all.  Its own
     Opus-lane slice; see `macos_f90/CLAUDE.md`.
     2c may proceed on the analyzer-at-detector design below.**
     ~~BLOCKED 2026-07-27 on an ENGINE DEFECT, not on the 2c design.~~
     Second attempt (2026-07-27) got as far as validating the machinery and
     then found that the pupil polarization state itself is untrustworthy on
     any train with an ODD number of mirrors: `Reflector` assembles the
     reflected field in a p̂-follows-the-outgoing-ray basis but multiplies by
     `−r_p`, so near normal incidence the transverse field is REFLECTED about
     the local p̂ instead of negated.  One on-axis mirror turns x-polarized
     light into a 50/50 x/y mixture (`Py/Px = 1.0163` on `Rx_Cass_FarField`
     at <2° AOI, where physics allows ~1e-3); a mirror PAIR cancels it
     exactly, which is why every Phase-1/2/3a gate passed.  Diagnosed,
     verified by scratch build (restoring `+r_p` gives `2.07e-4` for one
     mirror and a BIT-IDENTICAL number for two), reverted, nothing landed.
     **Decision packet: `REVIEW_POL_SP_SIGN_2026-07-27.md` (Fable lane — it
     is a convention change in the Fresnel core).  Reproducer:
     `MACOS_resources/mmacos/tools/pol_sp_sign_probe/`.  Engine-side summary
     in the polarization section of `macos_f90/CLAUDE.md`.**  2c cannot
     produce an honest floor until this is settled: the co/cross split would
     be dominated by the artifact.
     **The 2c DESIGN is, however, now settled and simpler than either
     previous sketch** (both superseded — do not resume from them).  This
     is the design that shipped; §2c below is the authority on what the
     code actually does, and these bullets are the reasoning that led to
     it:
     - The chain is **linear in the input Jones state** — measured
       4.2e-16 (45°) / 5.8e-16 (circular) on `Rx_Coro` at model 512.
     - Therefore a spatially uniform analyzer **commutes with propagation**
       (Tranche 1 propagates every component with the identical kernel), so
       the co/cross split is taken **at the detector** on the engine's own
       component planes: `complex_field(det,'plane',1..2)`, projected on the
       analyzer and its complement.  Parseval on that split is exact
       (3.0e-16 … 4.5e-16), and the incoherent sum is automatic.
     - **No pupil multiplier is needed at all**, which retires BOTH the
       Jones-pupil-multiplier design (Finding 2 below) and its "divide by the
       scalar run" successor.  `pupil_propagator` / `apodize_complex` are not
       in the loop.
     - Analyzer choice: the dominant eigenvector of the 2×2 **pupil
       coherency matrix** `C_ij = Σ E_i E_j*`.  Phase-insensitive (the common
       wavefront cancels), so unlike a plain pupil mean it does not collapse
       on an aberrated pupil, and it is by construction the analyzer that
       minimizes cross-polarized power — which is the operational definition
       of "co-polarized" that Finding 1 demands.
     - Incidental: `complex_field(..., 'reset_trace', false)` returns
       bit-identical planes ~100× faster (0.01 s vs 0.83 s at model 512), so
       three component planes cost ONE propagation.
     **Fixture correction:** `Rx_Coro.in` declares `nGridpts=511` and must be
     run at model **≥ 512** (Dave: grid size ≤ model size, or `MREset`).
     Below that the engine prints `Too many grid points. Resetting npts to N`
     and then intermittently SIGSEGVs in `intensity` (~30–50% of runs at 128;
     deterministically under `trace`).  This supersedes the "runs at model
     128" note below.  Separately, `macos.trace(e)` on `Rx_Coro` loses every
     ray past element 7 while `intensity`/`complex_field` propagate fine —
     unrelated to polarization, not chased.
     **First-attempt findings (2026-07-26), retained for the record.  WIP
     parked at `~/dev/MACOS_sandbox/pol_2c_wip/` (moved OFF the MATLAB path
     deliberately: a half-working `+macos/` function is reachable by users);
     it is now superseded by the design above.**
     - `macos.pupil_propagator(pupilElt,detElt)` WORKS and its contract
       gate is exact: `p(ones)` is BIT-IDENTICAL to the plain scalar run
       (intensity→apodize_complex→intensity with `reset_trace=false`;
       passing the default `true` silently discards the imprint and
       returns the nominal PSF, which reads as "polarization has no
       effect" rather than as a bug).
     - **FINDING 1 — "co-polarized" must be defined against the MEAN
       OUTPUT state, not the input state.**  A real train rotates
       polarization geometrically with zero diattenuation and zero
       retardance; projecting on the input bills that uniform rotation as
       cross-polarized light.  It is neither an aberration (a perfectly
       unitary system produces it) nor a floor (you would orient the
       analyzer to it).  On the stock conductor `Rx_Coro.in` chain that
       error reports a ~50% cross fraction on a system whose measured
       diattenuation is 1.2e-15.  Same mean/variation discipline
       `pol_maps` already enforces.
     - **FINDING 2 (the blocker) — the Jones pupil CANNOT be used
       directly as a pupil multiplier.**  `jones_pupil` is built from
       `RayE`, so it carries the accumulated OPL phase: measured
       `|mean J11| / median|J11| = 3e-4` at the Rx_Coro pupil, i.e. the
       phase is fully scrambled across the pupil (many waves of
       wavefront).  Multiplying `WFElt` by it DOUBLE-COUNTS the
       wavefront the engine already applied.  Worse, the RayE↔WFElt
       phase relation is TRAIN-DEPENDENT (the Tranche-1 finding), so
       there is no fixed conjugation to divide out.
       *Refuted along the way:* it is NOT a ray-grid vs diffraction-grid
       mismatch — the masks match exactly (3176/3176, same centroid and
       extent), so do not spend time there.
     - ~~**Next attempt should build the multiplier entirely inside
       `WFElt`**, using the NEW plane getter: per-component pupil fields
       `complex_field(pupilElt,'plane',k)` divided by the scalar-run
       `complex_field(pupilElt)`.~~  **SUPERSEDED 2026-07-27** — no pupil
       multiplier of any kind is needed, because a uniform analyzer
       commutes with propagation (see the settled design above).  The
       division would also have re-imported the train-dependent RayE↔WFElt
       phase relation through the vector run's seed, so it was not the
       clean sidestep it looked like.
     - Fixture note (**CORRECTED 2026-07-27, see above — `Rx_Coro.in`
       needs model ≥ 512; it only appears to run at 128**): the original
       note read "`Rx_Coro.in` and `Rx_Coro_noLyot.in` both run at
       model 128"; `Rx_Coro_FPM.in` returns an ALL-ZERO intensity there
       (it SIGSEGVs at 256 and is only usable at 1024 — pre-existing, see
       `mmacos/CLAUDE.md`).  So a 2c evidence section either uses the
       weaker chain at 128 or the driver must split into per-model-size
       batch invocations.
  4. ~~**Phase 3 polarizer + waveplate**~~ — **LANDED 2026-07-27; the one
     convention decision it raised is now SETTLED AND GATED (see the
     material-axis block at the end of this item).**  `TrPolarizerElt(15)` finished from its
     name-table-only stub and `WavePlateElt(18)` added (`mEltTypes`
     17→18), both served by a new `PolElt` in `elemsub.F`; `PolAxis=` /
     `Retardance=` keywords in both parse chains + SAVE; `polelt_set`/
     `polelt_get`/`jmat_elt_get` (codegen Path A) → `macos.polarizer` /
     `macos.waveplate` / `macos.elt_jones`.  `JmatElt(2,2,mElt)` is no
     longer dead.  Gates: `tPolElement`, **23 pass**, all closed-form and
     written from the textbook — Malus 1e-12, crossed extinction EXACTLY
     0, QWP `S3/S0 = −1` (signed, so it pins the retardance convention),
     HWP 2θ slope 2 to 1e-10, two-QWP≡HWP 1e-15, unitarity 1e-14 for
     linear AND circular states, and pol-off BIT-IDENTICAL to a
     Reference-surface twin fixture.  Each mechanism has an in-suite
     non-vacuity A/B (R=0 collapses the QWP and the 2θ law; the polarizer
     is checked to FAIL unitarity).
     **The off-normal axis convention — SETTLED AND GATED 2026-07-27
     (Fable decision, Opus landing).**  Off normal an ideal polarizer
     carries a real O(sin²AOI) question: declaring the PASS axis and
     projecting it is not equivalent to declaring the BLOCK axis and
     transmitting the complement.  Size: **3.56° of transmitted-axis
     orientation at 20° AOI**, closed form `acos(2cos a/(1+cos²a))` at 45°
     azimuth, **identically zero at normal incidence and when the axis
     lies in or perpendicular to the plane of incidence** (so the obvious
     test tilt is degenerate and reports a meaningless zero).  **The rule
     is PROJECT THE MATERIAL AXIS** — the absorbing direction for a
     polarizer, the crystal fast axis for a waveplate.  Not a taste call:
     both constructions are in the literature (pass-axis = Fainman &
     Shamir, *Appl. Opt.* **23**, 3188 (1984); material-axis = Korger et
     al., *Opt. Express* **21**, 27032 (2013) Eq. (5)–(6)) and Korger et
     al. decided between them by measuring the Mueller matrix of a tilted
     dichroic polarizer.  `PolElt`'s polarizer branch was flipped
     accordingly (keyword semantics unchanged, `PolAxis=` still declares
     the pass axis); the waveplate branch was already the settled rule.
     **Normal incidence is bit-identical across the flip** — verified
     against pre-flip captures, and `PolElt` deliberately omits a
     redundant `DUNITIZE` to keep it so.  New fixture
     `Rx_PolElt_Tilt.in` tilts the BEAM (tilting the element does nothing
     to a collimated on-axis bundle); gates in `tPolElement` section F,
     both dispatch chains: transmitted axis on the material rule to
     6e-17 rad and missing the pass-axis one by the closed form to 1e-13°,
     and — grid side — the crossed-analyzer null moves 7.11° between the
     rules, giving 9.1e-33 vs 1.53e-2 of relative detector power.
     Non-vacuity measured by rebuilding the pass-axis engine: both new
     engine gates fail (3.5616°, and the null/leak values swap).  The
     degenerate-azimuth gate passes on BOTH engines, which is exactly its
     job.  Decision packet: `REVIEW_POL_ELEMENTS_2026-07-27.md`;
     evidence: polval §6.7.  Landed as macos `216c56c` +
     MACOS_resources `a3417ce`.
     **Still NOT landed: `RfPolarizerElt(14)`.**  The convention no longer
     blocks it; a reflective wire grid simply carries more physics (grid
     reflection efficiency, the substrate's own s/p response) and nothing
     needs it yet.
     **Finding worth propagating: there are TWO element dispatch chains.**
     `propsub.F`'s `CPROPAGATE` re-traces the rays that seed the
     diffraction grid through its own `EltID` chain, so an element wired
     only into `tracesub.F` works perfectly in `ray_field` and is INVISIBLE
     to `intensity`/`complex_field` — measured at 3.6e-33 ray power against
     a full 9.69e-01 detector plane, which reads as "polarization has no
     effect on the image" rather than as a bug.  Both chains are now wired;
     `test_grid_carries_the_polarizing_train` is the tripwire.
     Not attempted here: VVC, and the polarization phase-shifting
     Twyman-Green example (Phase 2d, needs the Bench emitters on
     `pol-ifo`).  pymacos shims for the three new api routines are
     deferred — mmacos-only, as with the `view_rx` engine leg.
  5. **Phase 2d mechanical half** (on `pol-ifo` once parents merge): Bench
     coating emission (`add_mirror`/`add_bs_reflect` coating options,
     `Coating=` block ordering handled), `twyman_green('coating',...)`
     pass-through, visibility-budget driver.
  6. **Phase 4 spatial coatings** (§4): zone-map model per the AmplMat
     template + the two stated departures.
  7. ~~**Validation-report skeleton + Phase 0–2b evidence**~~ — **DONE
     2026-07-26** (see the STATUS block in the Validation document
     section below).  Shipped the `polval/` report, the `make polval` /
     `polval-pdf` / `polval-regen` / `polval-check` targets, the
     regeneration driver, AND — beyond the stated scope — the Phase 3a
     Tranche 1 evidence sections, which back-fills the gap Tranche 1
     left open.
  Per the new standing rule, every item above ships its cmdref + manual
  entries **and its validation-report evidence section** as part of its
  definition of done.
- **Fable lane — scarce; keep to short, judgment-critical sessions:**
  phase-gate reviews of the Opus items (line-review the diff, run the ifx
  smoke + both-compiler suites, check the gates weren't satisfied
  vacuously); the **VVC** (retardance/handedness conventions, broadband
  leakage); **Phase 3a Tranche 2** (per-ray `J_run` engine surgery) when
  IFO `DO_NEARFIELD` needs it; the **Phase 2d interpretation half** (PSI
  polarization systematic, the AOI-trade conclusions — error-budget
  numbers, where "tests pass, plausible number" is not enough).  The
  2026-07-25 bench_ifo_dm diagnosis (a plausible "9 nm algorithm floor"
  that decomposed into a branch-cut artifact + a sign convention + real
  retrace physics) remains the reference case for what stays here — as
  does Phase 2a/2b itself, where building the Jones-pupil layer exposed
  three latent engine bugs the templates would have built on top of.

**Sequence:** Phases 0–2b are DONE (gate-passed, pushed, documented).  The
Opus lane works its list top-down on `pol-core` (Tranche 1 first — it
unblocks the Phase-3 VVC acceptance test and the polarized PROPER re-runs);
Fable gates each landing.  `pol-ifo` (2d) follows after `bench-builder` and
`pol-core` merge.

**Opus-lane session discipline (added 2026-07-27, after the end-of-session
error cluster on 2026-07-26).**  The 26th's session landed items 7 + 2 +
the plane getter + the review packet + started 2c in ONE continuous run;
the errors (a 740 MiB `git add -A` sweep, a provenance stamp citing a
commit its own reset had destroyed, a watcher script that matched its own
process) all clustered at the END.  That is a session-length failure mode,
not a competence one — the same session's early work was clean.  Rules:

1. **One worklist item per session.**  Land it (code + gates + packet
   entry + push), then END the session.  Start the next item with fresh
   context.  Do not "start the next item while the suite runs."
2. **Stage by explicit path, never `git add -A`/`git add <dir>`.**  Run
   `git diff --cached --stat | tail -3` before every commit; anything
   over ~2 MB staged needs a stated reason in the commit message.
   (The gitignore now blocks the known artifact trees — MACOS_resources
   `cfe2fbb` — but the rule stands for trees it can't foresee.)
3. **After any commit rewrite (reset/amend/rebase), grep the worktree
   for the old SHAs** before pushing: `git grep -I <oldsha>` plus a
   plain `grep -rn` over untracked docs.  A stamp, packet, or doc
   written minutes earlier may cite a commit that no longer exists
   (the f10b234 dangling-stamp bug, fixed in macos `f8c0f34`).
4. **Watchers/pollers must watch external state** (a file mtime, a PID
   captured at spawn, a log line count) — never a `pgrep -f`/`grep`
   pattern that appears in their own command line.  Prefer foreground
   waits; if polling, print what changed each tick so a self-watching
   loop is visible in the first two lines of its output.
5. **When a gate needs more than two fix attempts, or a push takes
   minutes instead of seconds, STOP and write the state down** for the
   next session instead of pushing through — both 26th incidents would
   have been caught at this tripwire.

---

### What already exists in the engine (verified 2026-07-25, line-checked)

- **Per-ray vector E-field** `RayE(3,mRay)`, `Complex*16` end-to-end
  (`elt_mod.F:450,846` — double precision confirmed; updated per-ray in `CPROPAGATE`
  at `propsub.F:507-509` and via the `Reflector`/`Refractor` `Evec→Eout` argument pair,
  `elemsub.F:206-216`).
- **Fresnel + recursive multilayer thin-film coatings**, complex index stored as
  `DCMPLX(n, −κ)`, in `Reflector` (`elemsub.F:432-547`) and `Refractor`
  (uncoated Fresnel + coated recursion with **live TP/TS transmission**).
  All **gated on `ifPol`** — inert when polarization is off.
  ⚠ `Reflector`'s transmittance sub-blocks are `if (.false.)` dead code — reflective
  elements return R only (fine for mirrors; relevant if anyone expects a coated
  reflector to report transmitted leakage).
- **Vector diffraction** (3 independent FFTs of Ex/Ey/Ez), **far-field FFT leg only**
  (`propsub.F:1380-1426, 1604`, gated `ifVecDif3`).
- **Source polarization** `Ex0/Ey0` (`src_mod`); CLI `POLARIZATION`/`NOPOLARIZATION`/
  `VECTOR`/`SCALAR` (`macos_cmd_loop.inc:1072-1114`). `VECTOR` requires `ifPol` first;
  the `POLARIZATION` path enables `ifVecDif3` only when `mWF≥3` — **all stock models in
  `macos_param.txt` have `mWF=3`**, so this is available at every model size, but the
  API must still assert it.
- **Coating wavelength semantics (verified — corrects the earlier draft).** The parser
  converts the Rx thickness spec (units of waves) to **physical thickness** at load
  time via `Thk·Wavelen/IndRef` (`msmacosio.inc:2648-2675`); the trace-time recursion
  then applies phase `2·(2π/λ_now)·d·N` with the **current** wavelength
  (`elemsub.F:512-518`). So a loaded stack is physically fixed and **broadband sweeps
  are already correct**. The SAVE writer already inverts the normalization
  (`iosub.inc:1747-1766`), so save/load round-trips preserve the physical stack. The
  real remaining wrinkles: (a) the Rx spec means "waves at whatever `Wavelen` was
  current at parse" — a file-order / cross-λ-interpretation coupling to document, and
  (b) `Coating=` must follow `IndRef=` in the file (the parser snapshots the boundary
  indices). Phase 1's `coat_set` API takes **physical thickness** and sidesteps all of
  it.

### Gaps this plan closes

- Zero binding exposure (can't even turn `ifPol` on from Python/MATLAB).
- `RfPolarizerElt(14)`/`TrPolarizerElt(15)` are **name-table-only stubs** (verified: no
  dispatch anywhere in `tracesub.F`/`propsub.F`/`elemsub.F`); `JmatElt(2,2,mElt)` is
  allocated/zeroed/deallocated and never otherwise touched (dead — confirmed).
- **`Ex0Ey0` prescription keyword is silently swallowed** (worse than dead): a live
  `ELSE IF (LCMP(VAR_NAM,'Ex0Ey0',6)) THEN GO TO 50` no-op at `msmacosio.inc:176`,
  next to the commented-out `PolBeam`/`PolSrc` block (150-174). A user setting source
  polarization in a `.in` gets no error and no effect.
- No Jones/Stokes/Mueller math, no polarization metrics, no contrast floor, no IFO
  polarization diagnostics.
- No waveplate/retarder, no polarizer, no vector-vortex mask.
- Coatings are per-element uniform; the one spatial map (`AmplMat`) is real-amplitude
  only.
- **`srtrace.F` hardcodes `ifPol` — as a `Logical,parameter :: ifPol=.false.`**
  (`srtrace.F:67`), so lifting it is a declaration + argument-threading change, not a
  flag flip. Whether it matters is a Phase-0 scoping question (srtrace serves the
  single-ray/chief-ray paths — SPOT beam frames, LocalCoord, sensitivity probes — not
  the pupil-grid trace the Jones pupil uses).
- **Vector propagation is not known to close the coronagraph chain.** Only the
  far-field FFT leg is vectorized. A VVC test needs pupil → FPM → Lyot → focal with the
  vector field preserved at every leg. Audited in Phase 0; if it doesn't close,
  "vectorize the chain" becomes its own Phase 3a.

---

## Conventions to pin down before any code

Cross-code polarization disagreements are almost always convention mismatches. Fix
once, document in `macos_f90/CLAUDE.md` (new engine-polarization section), assert in
tests:

| Convention | Decision |
|---|---|
| Time-harmonic phase | Derive from the engine, don't legislate: the coating recursion uses `N = n − iκ` with propagation factor `exp(−i·2kd·N)` (`elemsub.F:512-518`), and the 2026-07-25 IFO work established empirically that **engine field phase advances as OPL shortens** (frame-ratio diagnostic, exact to 1e-4 pm; same convention the pymacos↔PROPER `opd_sign_flip` reconciles). Phase 0 writes down the single consistent (ω,k) convention these imply and pins it with a test. |
| Absorbing index | n − iκ as implemented; κ > 0 is loss under the engine's phase convention (verify sign once in Phase 0 with a lossy-metal reflectance check against tabulated Al). |
| Jones storage basis | **Linear (x,y)**, with a basis-change utility for circular. Circular is a fixed unitary transform; storing linear and converting is strictly simpler than the reverse. |
| Reference frame | **Double-pole (dipole)** by default; local ŝ/p̂ and global as options. See Phase 2a. |
| Retardance sign | Fast axis leads. State it; VVC handedness depends on it. |
| Diattenuation | Intensity-based, D = (T_max − T_min)/(T_max + T_min) ∈ [0,1]. |

---

## Phase 0 — Audit & scoping (no new features; gates everything downstream)

Short, cheap, can reshape the plan. Deliverable is a findings note, not code. Items 3
and 6 were pre-answered during this review and only need the recorded confirmation.

1. **Vector propagation leg coverage.** Enumerate every propagation leg (`PFFPROP`,
   near-field/Fresnel, angular-spectrum, pupil-plane products, Lyot application) and
   record vector-capable vs scalar-collapse. Output: legs × capability table. **If the
   coronagraph chain does not close, Phase 3's VVC acceptance test is rescoped and
   "vectorize the chain" becomes Phase 3a.** Note for track A: the IFO detector-mode
   PSI needs only trace-to-detector fields — it closes with Phase 1 alone, so track A
   does not wait on this finding.
2. **`srtrace.F` `ifPol` scoping.** Which analyses route through srtrace (chief-ray/
   SPOT/LocalCoord/single-ray probes), and do any polarization use cases need them?
   Expected answer: not blocking for Jones-pupil work; record it.
3. **Coating/wavelength semantics — CONFIRMED during review** (see above). Remaining
   Phase-0 work: one empirical round-trip test (load → save → reload at a different
   `Wavelen`, compare physical stacks) to pin the SAVE-writer inversion, plus document
   the parse-order couplings.
4. **Precision audit.** `RayE` confirmed `Complex*16`; sweep the remaining chain (the
   vector FFT arrays, coating recursion locals — spot-checked `COMPLEX*16` in
   `Refractor`) and record. At 1e-10 contrast we track field amplitudes near 1e-5; a
   single-precision stage is fatal and invisible.
5. **Coating subsystem inventory — two models confirmed:** the polarization-path stack
   `EltCoat/IndRefArr/ExtincArr/EltCoatThk` and the separate `nCoatElt/CoatIndxElt/
   CoatThkElt` (+`AirGap`) model (both appear in the SAVE writer, `iosub.inc:1738-1766`
   — note the pre-existing `nCoatElt` size ToDo there). Document both; decide with
   Dave whether to unify or extend only the polarization-path one.
6. **Phase-convention consistency** (from the conventions table): one note deriving
   (ω,k) sign from the recursion + the IFO empirical finding, plus a lossy-reflectance
   sanity check.

---

## Phase 1 — Expose existing polarization physics (FIRST deliverable)

**Binding mechanics (corrected to match the actual toolchain).** New API routines go
in `macos_api_mod.F90` following the `beam_set`/`beam_get` template
(`macos_api_mod.F90:4482-4541`): `SystemCheck()` guard → marshal → drive `SMACOS` or
read `use`d module state → `OK=PASS`. Then:

- **mmacos:** most routines ride the **codegen Path A** — add the routine, re-run
  `python3 gen_mex_wrappers.py` from `mmacos/src/`, done (no `mmacos_mex.F` edits).
  Only `rayfield_get` needs a **hand-written Path B helper** (complex split + many
  outputs; follow `do_complex_field`, add to `HAND_WRITTEN`). Veneers in `+macos/*.m`
  + `Session.m` per the naming/validation conventions in `mmacos/CLAUDE.md`.
- **pymacos:** bind in `pymacos_f2py.f90` (NOT "pymacosf90.pyf" — the binding is an
  f2py Fortran shim) + `macos.py` methods.
- Every routine: `matlab.unittest` test + pymacos pytest. Work lands on `sls-dev`
  (both repos); full `run_mmacos_tests.sh` + pymacos pytest green pre-commit; engine
  rebuild chain per `reference_build_test_workflow` (rm mex to force relink).

**Engine API (`macos_api_mod.F90`):**

- `pol_set(OK, on, ExRe,ExIm, EyRe,EyIm)` / `pol_get(OK, on, ...)` — drive the
  `POLARIZATION`/`NOPOLARIZATION` logic and set `Ex0/Ey0`; read back `ifPol`,
  `ifVecDif3`, `Ex0/Ey0`. Mirror `macos_cmd_loop.inc:1072-1097`. Assert `mWF≥3`
  before enabling vector diffraction and return a distinct failure code if not (all
  stock models pass, but a custom param file may not).
- `vecdif_set(OK, on)` — `VECTOR`/`SCALAR`. Enforce the CLI's ordering rule
  (polarization must be on) with a clear error rather than the CLI's silent reset.
- `coat_set(OK, iElt, nLayer, IndRef(:), Extinc(:), Thk(:))` / `coat_get(...)` —
  read/write `EltCoat/IndRefArr/ExtincArr/EltCoatThk` **in physical thickness units**
  (element length units); the waves-normalization is an internal storage detail applied
  on set and inverted on get. `coat_get∘coat_set == identity` to round-off is a
  required test. Must set `modified_rx` (grid-setter dirty-the-trace convention). Must
  also populate the boundary-layer entries (`IndRefArr(0,·)`, `(nCoat+1,·)`) the way
  the parser does — the recursion reads them.
- `rayfield_get(OK, ExRe,ExIm, EyRe,EyIm, EzRe,EzIm, kx,ky,kz, nx,ny,nz, status, N,
  iElt)` — harvest `RayE(3,:)` per ray **plus ray direction cosines, surface normal
  at `iElt`, and the ray-status mask** (buffer-getter pattern; `ray_status_get` is the
  precedent for per-ray `elt_mod` harvests). Phase 2 cannot build a correct Jones
  basis without the geometry, nor honest statistics without the mask; adding them later
  changes a published signature.

**Fix `Ex0Ey0`** — replace the silent swallow at `msmacosio.inc:176` with a real parse
(complex pair ×2), and add the SAVE-writer line (cross-ref `SAVE_KEYWORD_AUDIT.md`) so
source polarization round-trips through a `.in`. Decide with Dave whether to also
revive `PolSrc=` (the commented block) as the on-switch keyword or leave state control
to the API/CLI only.

**Bindings:** `macos.polarization(...)`, `macos.vector_diffraction(...)`,
`macos.coating(...)`, `macos.ray_field(...)` + Session methods + pymacos equivalents.

**Tests/examples:**
1. **Geometry invariance (prerequisite for Phase 2).** On a coating-free system,
   polarization-off and polarization-on traces produce **bit-identical `ray_status`
   and OPD**. Phase 2 assumes ray-by-ray correspondence between traces with different
   input states; this test makes that assumption safe.
2. **Round-trip + wavelength stability.** `coat_get∘coat_set` identity; a stack set at
   λ₁ and evaluated at λ₂ reproduces the analytic Fresnel result for that physical
   stack at λ₂ (should pass already — the semantics are correct; the test pins them).
3. Turn `ifPol` on for `Rx_Cass_FarField.in`, Al coating (tabulated n,κ) on an oblique
   fold; vector-diffraction PSF differs from scalar; coating changes throughput;
   reconcile with the PROPER far-field baseline where polarization-neutral.
4. **Track A smoke:** `ifPol` on the bench_ifo Twyman-Green with a bare-metal Al BS
   coating — fringe visibility drops below the scalar 0.868 by the s/p reflectance
   split; number checked against the single-surface Fresnel analytic at 45°.
5. **`ifPol` off is a no-op.** GMI regression stays 6/6 (its Rx are polarization-off).

---

## Phase 2 — Jones pupil, polarization metrics, and the two application deliverables

> **STATUS 2026-07-26: 2a + 2b core LANDED on pol-core** (mmacos
> `jones_pupil`/`pol_maps` + pymacos ports + `tJonesPupil`/`test_jones_pupil.py`).
> Gates green at round-off: unitarity (stock conductor Cass, D/ret < 3e-15),
> Fresnel-analytic fold (Bench rig + thick Al, RS/RP and D vs closed form at
> 1e-14), 2θ symmetry (orientation locks azimuth to 5e-14), D basis-invariance,
> synthetic polar-decomposition identity.  En route, THREE engine fixes (see
> `macos_f90/CLAUDE.md` polarization section): `pol_set`/`vecdif_set` now dirty
> the trace (stale-`RayE` on state change); coated-branch incident medium
> (conductor bleed via `IndRefArr(0)`); coated-branch signed incident cosine
> (reciprocal 1/r coefficients, |R|>1 — |D| survived, which is why intensity
> tests never saw it).  The "no engine change beyond Phase 1" premise below
> did not survive contact with the coated branches — the exposure work was
> exactly how these were found.  Remaining in Phase 2: the 2b low-order
> polarization-aberration expansion (optional), 2c contrast floor, 2d IFO
> deliverable.

The engine has no Jones matrix. Build it the standard PRT way **in the binding layer**
(no engine change beyond Phase 1). This is where the physics can go quietly wrong, so
it is specified in more detail.

### 2a. Jones pupil — and the reference-frame problem

`macos.jones_pupil(srf, 'basis', ...)` (mmacos + pymacos): trace twice with orthogonal
input states ((1,0), (0,1)) via `pol_set`, harvest `RayE` + geometry at the pupil via
`rayfield_get`, assemble the 2×2 Jones per pupil point.

**The basis is not a detail.** The 3×3 PRT matrix is basis-independent; the 2×2 Jones
pupil is not. The naive local ŝ/p̂ basis carries a coordinate singularity that shows up
as **spurious retardance with a tilt/astigmatism-like signature** — at 1e-10-relevant
levels the artifact can exceed the physics. Standard remedy: the **double-pole
(dipole) coordinate system** (Chipman/Lam/Young).

- `basis`: `'double-pole'` (**default**), `'local-sp'`, `'global'`; docstring states
  `'local-sp'` is diagnostic-only, never for budget numbers.
- Return N×N×2×2 complex **plus the vignetting mask**; vignetted points are `NaN`, not
  zero (zero-filling seeds RMS statistics with a ring of perfect nulls).

**Gating acceptance test — unitarity.** Reflectivity forced to unity, no coating: the
Jones pupil is **unitary at every unvignetted point** and reduces to identity (within
the system's geometric rotation) in the correct basis. Cheap; catches basis,
normalization, and sign errors in one check. Regression on every commit touching this
path.

### 2b. Polarization-aberration metrics

`macos.pol_maps(jones)` — pure MATLAB/NumPy:

- **Polar/SVD decomposition** J = H·U per point; diattenuation from singular values,
  retardance from the eigenphases of U.
- **Retardance branch handling:** eigenphase differences are mod 2π; unwrap across the
  pupil and flag ambiguous points instead of silently choosing.
- **Separate the pupil-mean from the spatially-varying part and report both.** Uniform
  retardance is a state change, not a contrast term; only the variation drives the
  floor (track B) or the PSI systematic (track A). One conflated RMS is the likeliest
  way to produce a plausible-looking wrong number.
- **Optional expansion into the standard polarization-aberration terms** (piston/tilt/
  defocus/astig in diattenuation and retardance) — makes two-mirror results directly
  comparable to the published literature; worth it for validation alone.

### 2c. Coronagraph contrast floor (track B deliverable)

> **LANDED 2026-07-27** (Opus lane item 3) — `macos.pol_contrast_floor(pupil,
> det, ...)` + `tPolContrast` (model 256, 14 tests) + `tPolContrastCoro`
> (model 512, 6 tests) + report §5 + cmdref/manual.  The signature below is
> the SHIPPED one; the original `pol_contrast_floor(jones, stokes_in,
> coronagraph_fn)` is retired — the Jones pupil cannot be a pupil multiplier
> (Finding 2, endorsed by the Fable lane 2026-07-26) and no multiplier of any
> kind is needed once the split is taken at the detector.  Headline numbers
> and the scope caveat are in the closing block of this section.

`macos.pol_contrast_floor(pupil, det, 'input', ..., 'coatings', ..., 'dark_zone', ...)`:

1. Trace with vector diffraction on and read the detector's **own component
   planes** `complex_field(det,'plane',1..3)`.  No pupil multiplier: the chain
   is linear in the input Jones state (4.2e-16) and Tranche 1 propagates all
   three planes with the identical kernel, so a spatially uniform analyzer
   **commutes with propagation** and the split may be taken at the end.
2. Project on an **analyzer derived from the pupil coherency matrix**
   `C_ij = Σ E_i conj(E_j)` (dominant eigenvector) and on its orthogonal
   complement → co-polarized / cross-polarized; `|Ez|²` is the third,
   longitudinal channel.  Referencing to the MEAN OUTPUT state rather than the
   input is Finding 1; coherency is used because it is phase-insensitive, so
   it does not collapse on an aberrated pupil.
3. Sum intensities **incoherently** across the run set — an unpolarized source
   is TWO traces (x-in, y-in), each with its own analyzer; the second state is
   never synthesized from the first.
4. Return the floor **broken out by component** (`.floor.co/.cross/.long`,
   `.contrast_cross`, `.floor.dark_zone`), plus `.sweep` — the coating
   sensitivity, one full recomputation per coating set.
5. **Report the scope, measured.**  `.scope` compares the pupil cross-polarized
   fraction seen by the GRID against the same quantity from `RayE`, per input
   state; `macos:pol_contrast_floor:tranche1` warns when they disagree.  Ratio
   maps NaN-mask small denominators; they are never zero-filled.

**Measured (2026-07-27).**  `Rx_VecChain` — a polarization-neutral train — gives
cross **exactly 0**, and the co-polarized channel reproduces the scalar run to
1.33e-15 of peak and its contrast curve to 2.36e-15: the decomposition invents
no floor.  `Rx_Cass_FarField` (both mirrors before its single far-field leg, so
the grid carries the whole train): Parseval on the split 1.3e-24, energy
closure 1.8e-16, uncoated floor 7.0612e-07 of the co-polarized power — the same
number §4's ray-level probe reports, reached through the grid planes instead.
Coating sensitivity, annulus-mean cross contrast: bare 9.59e-12 → bare Al
3.34e-10 (27.9× the cross power) → MgF₂/Al 1.99e-09 (151.3×).  **The coating,
not the geometry, sets the floor on that train.**

**Scope — and a Tranche-2 finding this produced.**  `Rx_Coro` puts six of its
seven mirrors AFTER the first physical-optics leg, which is outside Tranche 1's
validity condition, so its floor (peak cross contrast 1.27e-09; 20–80 px annulus
mean 5.79e-13) is a **lower bound**.  The cost is now quantified rather than
qualitative: the grid carries **0.8412** of the ray-level cross-polarized
fraction bare and **0.5653** with the mirrors coated, and the coating
sensitivity it reports there comes out at **−3.2%** while the ray-level
fraction RISES **+59%** — the wrong *sign*, not merely the wrong size, because
only the first mirror's coating precedes the seed leg.  That is the sharpest
argument yet for §3a.3 Tranche 2, and `tPolContrastCoro`'s last two tests pin
these numbers so Tranche 2 has to come back and change them.

**One bug worth remembering.**  The coherency matrix was first written with
MATLAB's `'` on the wrong operand, building `conj(C)` — whose dominant
eigenvector is the CONJUGATE analyzer.  Identical for any linear input state
(real eigenvector), and exactly ORTHOGONAL to the truth for a circular one:
reported cross/co jumped from 1.4e-6 to 7.1e+05.  x / y / 45° all passed
vacuously.  The circular input state is load-bearing in `tPolContrast`, for the
same reason the 45°/circular states are load-bearing in `tVecChain`.

### 2d. IFO polarization metrology (track A deliverable — NEW)

`bench_ifo_pol` example (+ small `+macos/+design` helpers), built on the existing
bench_ifo/bench_ifo_dm machinery (differential PSI, `PSI_MODE`, per-ray pupil map):

1. **Visibility budget.** Coated 45° BS + folds: predicted fringe visibility from the
   s/p amplitude split vs the measured visibility in the simulated interferograms.
2. **PSI polarization systematic.** Run the PSI pipeline with the source polarization
   at 0°, 45°, 90°: the recovered surface shifts by the s/p retardance difference of
   the train. Report the apparent surface error vs input state — this is the number a
   metrology error budget wants (and it is measured with machinery already proven to
   3e-5 pm, so anything it shows is real polarization physics, not pipeline noise).
3. **The AOI trade, closed.** Repeat at BS AOI 45°/30°/22.5° (the `ifo_aoi_sweep`
   scaffold exists): scalar metrology was AOI-neutral; this quantifies the actual
   polarization payoff of smaller AOI — the open question from 2026-07-25.
4. **Bench builder wiring:** coating options on `add_mirror`/`add_bs_reflect`/`plate`
   (emit `Coating=` blocks + `IndRef` ordering handled), `twyman_green('coating',...)`
   pass-through, so polarization studies run through the same one-call rig builder.

**Defer:** populating dead `JmatElt(2,2,mElt)`. Per-element Jones is not what either
application needs; the per-pupil-point export is.

**Tests/examples (beyond the two application deliverables):**
- Unitarity (2a) — gating.
- Folded 2-mirror at ~12° AOI, Al: diattenuation scales with AOI, **vanishes at normal
  incidence**; RMS retardance vs analytic single-surface Fresnel.
- **Rotational-symmetry invariant.** On-axis rotationally-symmetric system → the
  characteristic 2θ symmetry in diattenuation/retardance maps. Strong, cheap,
  basis-sensitive; a broken frame breaks the symmetry.
- Basis comparison `'double-pole'` vs `'local-sp'` on the same system, artifact
  magnitude documented in the worked example.

---

## Phase 3 — New polarizing elements (engine)

**Gated on Phase 0 finding 1 for the VVC only** — the polarizer and waveplate serve
track A through the ordinary trace path and do not wait on vector-propagation
closure.

Each element needs: EltID dispatch in `tracesub.F`/`propsub.F`/(as needed) `srtrace.F`,
a surface routine in `elemsub.F` transforming `Eout` from `Evec`, keywords in
`msmacosio.inc` + defaults in `iosub.inc`, SAVE writer, Phase-1-style API/binding
exposure — and, for track A, a Bench `add_*` emitter.

> **STATUS 2026-07-27: the polarizer and the waveplate LANDED** (Opus
> worklist item 4), and the off-normal axis convention they raised is
> **settled and gated** — project the MATERIAL axis (the second landing,
> same day).  `PolElt` in `elemsub.F`, dispatched from BOTH `tracesub.F`
> and `propsub.F`; gates in `tPolElement` (27), including section F on the
> tilted-beam fixture `Rx_PolElt_Tilt.in`.  The VVC is untouched and
> remains Fable-lane.  `RfPolarizer` is still a stub, now only because a
> reflective wire grid carries grid-efficiency and substrate s/p physics
> beyond the axis rule — see `REVIEW_POL_ELEMENTS_2026-07-27.md` and the
> worklist entry above.  Text below is the original spec; the deltas from
> it are recorded in the worklist entry.

- **Ideal linear polarizer** — ~~finish `RfPolarizerElt(14)`/~~`TrPolarizerElt(15)`.
  Keyword for the transmission-axis vector; `Eout` = projection onto the axis, reflect
  variant sends the rejected component. ~~Model on the s/p projection in `Reflector`
  (`elemsub.F:385-428`).~~  **Modelled on `RefSrf` instead** — the s/p basis is
  DEGENERATE at normal incidence (`Reflector` falls back to an arbitrary
  `svec` when `|ihat x Nhat| < 1e-10`), which is precisely the regime a
  polarizer or waveplate is used in, so the transverse basis is built from
  the element's declared axis rather than from the plane of incidence.
  **Document:** sequential tracing yields one output port per
  run — a **PBS requires two traces** (the Twyman-Green driver already runs per-arm
  traces, so this fits the existing example structure exactly).
- **Waveplate / retarder** — new element type (extend `mEltTypes`): settable retardance
  + fast-axis orientation; applies the retarder Jones in the local transverse basis.
  **Documented as a thin, non-ray-splitting idealization** — no o/e walk-off. Also the
  primitive for bounding **stress birefringence** in transmissive elements (e.g. the
  wedged Fresnel plate upstream of DM1).
  Landed as `WavePlateElt = 18`.  Retardance is stored PHYSICALLY
  ((n_slow−n_fast)·d) with the Rx keyword in waves at parse-time
  `Wavelen` — the same treatment `Coating=` thickness gets, so a plate is
  fixed glass and a wavelength sweep is chromatic.
- **Vector vortex (VVC)** — charge-l geometric-phase mask coupling circular states
  with exp(±i·l·θ). **Retardance is a first-class parameter, not idealized away:**

  ```
  J(θ) = cos(δ/2)·I − i·sin(δ/2)·(n̂(θ)·σ),  fast axis at l·θ/2
  ```

  The identity term is the non-vortex leakage, amplitude cos(δ/2); with δ = π + ε,
  leakage intensity = sin²(ε/2) ≈ (ε/2)² falls out automatically. Make δ spatially
  variable and wavelength-dependent so vendor retardance-vs-λ curves load directly and
  broadband null depth is computed honestly.

**Tests/examples:** crossed polarizers → ~0 transmission; QWP converts linear →
circular (Stokes check); **polarization phase-shifting Twyman-Green** (PBS + QWPs —
the classic instrument, now buildable; recovered surface must match the PZT-stepping
result on the same rig); VVC null depth vs ε against the analytic (ε/2)², with and
without a circular analyzer built from polarizer + waveplate.

---

## Phase 3a — Vector near-field propagation (chain closure)

> **STATUS 2026-07-26: TRANCHE 1 LANDED on `pol-core`** (Opus lane item 1;
> awaiting the Fable line review).  All of §3a.2 plus both §3a.1 fixes are
> implemented in `propsub.F`; ladder tests 1/2/4/5 ship as `tVecChain`
> (mmacos) + `test_vec_chain.py` (pymacos), and ladder test 3 was run.
> Engine-side detail lives in the new Phase-3a section of
> `macos_f90/CLAUDE.md`.  Summary of what a reviewer should check:
>
> - **Coverage** went beyond the enumerated call sites by two items, both
>   deliberate: `SPH2PL`/`PL2SPH` (PropType 10/11) got the same K-loop —
>   they are the same scalar-only class and leaving them out is a silent
>   hole; and the 13 ray-side clip sites + 2 taper sites now extinguish /
>   scale all three planes (via new `WFZeroPt`/`WFScalePt` helpers), which
>   the seed-once fix makes load-bearing.
> - **A THIRD defect of the §3a.1 class was found and fixed**: `RayE`
>   carries no aperture clipping, so seeding from it resurrected vignetted
>   rays.  Both polarized branches now gate the seed on `LRayPass`.  With
>   that, **polarization-ON/vector-OFF is bit-identical to
>   polarization-OFF** — it was wrong by 21% after one leg, 38% after two.
> - **The seed applies NO phase-convention bridge** — CONFIRMED CORRECT by
>   the Fable review (2026-07-26), with the mechanism now understood: the
>   Return-leg `RayL=-RayL` flip makes `CumRayL` (subtracts return legs)
>   and `RayE`'s C1 phases (add them) diverge structurally, so the
>   RayE-vs-scalar phase relation is TRAIN-DEPENDENT — same-convention on
>   the Return-terminated far-field train (measured: EP-field circular
>   concentration 0.994 same / 0.002 conjugate; signed tilt response
>   equal), conjugate on the plain trace-to-detector flow (measured slope
>   −0.9995 vs the OPD map).  No universal bridge exists; the behavioral
>   gates carry the correctness claim, and **Tranche 2's `J_run` must
>   track phases explicitly against the `CumRayL` bookkeeping** rather
>   than assume a sign.  Full statement + debug gotchas in the Phase-3a
>   section of `macos_f90/CLAUDE.md`.
> - **Non-vacuity was checked explicitly** (the pre-fix engine was rebuilt
>   and re-run): it fails the new gates at 0.21…0.38 relative error and
>   mis-states total power by 4–7%.  The 45°/circular input states are
>   load-bearing — an x-pol-only gate passes vacuously because all the
>   energy sits in plane 1, which the old single-plane propagator carried
>   correctly.
>
> **Gate results** (both compilers built; gfortran for mmacos/GMI, ifx for
> pymacos = the standing ifx smoke):
>
> | Gate | Result |
> |---|---|
> | Ladder 1 — energy per leg | vector total == scalar total, 0…2.2e-16 |
> | Ladder 2 — x-pol ≡ scalar | 4.5e-16 … 6.8e-16 on `Rx_VecChain` (also 45°, circular) |
> | Ladder 3 — polarized PROPER re-run | scalar + pol-scalar reproduce the committed 4.836e-13 macos↔PROPER residual exactly; vector differs 1.3e-2 from scalar at identical total power (attribution **closed 2026-07-26** — see below) |
> | Ladder 4 — chain closure | two-leg mask chain, vector ≡ scalar at round-off; mask throughput identical to 1e-14 |
> | Ladder 5 — single-hop A/B | vector far-field total 8.9377e-01 → 1.8155e+06 == scalar total (2.03e6 in intensity) |
> | mmacos full suite | 412 pass, 0 fail (fast 281 / masks 62 / freeform 46 / proper-512 10 / proper-1024 13) |
> | pymacos suite | 6645 pass; PROPER-compare 26/26 |
> | GMI regression | 6/6, bit-identical (`vs-ref = 0.000e+00`) |
>
> Docs shipped with it: manual §6 (VECTOR/SCALAR scope rewritten), cmdref
> `10_engine_commands.md` (VECtor/SCAlar NOTES) + `20_bindings.md`
> (`vector_diffraction` NOTES, regenerated), both binding docstrings, and
> the `macos_f90/CLAUDE.md` Phase-3a section.
>
> **ATTRIBUTION — CLOSED 2026-07-26** (was carried as an open item).  A
> plane-selectable getter now exists: `cfield_plane_get` in
> `macos_api_mod` → `macos.complex_field(srf,'plane',k)` /
> `pymacos.complex_field(srf, plane=k)`, k=1..3 = Ex/Ey/Ez, k=0 = the
> element's own wavefront (the historical behaviour, bit-identical).
> Requesting a component plane with vector diffraction OFF is REFUSED —
> in scalar mode plane k is an unrelated wavefront, not a component.
>
> With it, the vector/scalar difference on the off-normal Cass-FF train
> is measured rather than guessed, and **the original one-line
> attribution was half wrong**.  Two mechanisms, both driven by the
> out-of-plane content: (1) POWER REDISTRIBUTION, dominant — the scalar
> run seeds from `|RayE|`, so ALL the power including the out-of-plane
> part propagates in ONE plane, while the vector run leaves only
> f = 0.997890 in Ex; that is a near-pure rescale (1−corr = 4e-8).
> (2) Ey/Ez diffract to their own pattern.  So `Iv ≈ f·Is + Iy + Iz`,
> which drops the difference from 2.5638e-3 to 2.8983e-4.  The naive
> expectation (difference == out-of-plane intensity, 2.11e-3) is wrong
> by ~2×.  The 2.9e-4 that remains is a shape difference between the
> scalar field and Ex, consistent with their different seeds, and is NOT
> further verified.  Gates: `tVecChain` (3 new) + `test_vec_chain.py`
> (3 new); report §3.5.
>
> **Not done here:** Tranche 2 (§3a.3) is untouched, as scoped.  ~~Nor the
> validation-report evidence section~~ — **back-filled 2026-07-26** with
> worklist item 7: Tranche 1's evidence is now §3.1–§3.6 of the `polval/`
> report (per-leg energy, x-pol≡scalar residual maps, mask throughput,
> the single-hop A/B including the pre-fix number, and the polarized
> PROPER cross-check).  The attribution that section carried as an open
> item was CLOSED 2026-07-26 by the plane-selectable field getter — see
> the Tranche-1 status block above.

**Goal:** promote the near-field propagators (sphere→plane, sphere→sphere,
plane→plane, and the DFT legs) from scalar-only to vector, propagating each
Cartesian field component separately on the model of the far-field `PFFPROP`
(K=1,3 loop of the identical scalar kernel), so a physical-optics **chain**
(pupil→FPM→Lyot→focal; IFO recomb→detector under `DO_NEARFIELD`) preserves the
vector field.  Planned 2026-07-25 from a code+math review on pol-core; the
Phase-0 audit's leg table is the coverage map.

### 3a.0 Mathematical basis (checked — cite, don't re-derive)

- **Per-component propagation is rigorous.**  In a homogeneous isotropic
  medium each Cartesian component of a monochromatic E-field independently
  satisfies the scalar Helmholtz equation (∇·E-coupling enters only through
  boundary data).  The angular-spectrum propagator is a solution operator for
  scalar Helmholtz with planar data, so applying the SAME scalar kernel to
  Ex, Ey, Ez separately introduces **no vector-specific approximation**.  The
  implemented kernels use the Fresnel (quadratic) transfer function
  `exp(-iπλΔz(fx²+fy²))` (`NFPROP`/`PPPROP`, propsub.F:1955/2013) — paraxial;
  the vector legs inherit exactly the scalar legs' F/# validity envelope,
  no better and no worse.
- **Component coupling lives at surfaces, not in free space.**  Maxwell's
  inter-component coupling is (i) surface boundary conditions — already
  handled per-ray by the s/p machinery in `elemsub.F` (RayE transport), and
  (ii) the divergence constraint, spectrally
  `Ez(fx,fy) = −(fx·Ex + fy·Ey)/fz`.  The engine samples Ez from rays
  (`RayE(3,:)`) rather than enforcing it spectrally — consistent to paraxial
  order (fz ≈ 1/λ).  Optional diagnostic (cheap, worth adding): report the
  spectral-divergence residual of the seeded field as a fidelity meter.
- **Sziklas–Siegman reference-sphere legs** (`NFPROP`'s `(Z1/Z2)` factors,
  PropType 2/5): the S-S transform is a scalar conformal rescale + quadratic
  phase, polarization-independent, applied identically per component — valid
  per component to the same paraxial order.  True high-NA vector focusing
  (Richards–Wolf basis rotation across a fast sphere) is carried only through
  the ray-sampled boundary data; at high NA the error is O(NA²) relative —
  the SAME order at which the scalar Fresnel kernel itself degrades.  Record
  this envelope in the VVC acceptance tolerances; a Richards–Wolf integrator
  is out of scope.
- **Amplitude masks are scalar.**  `FFObscure` is a 0/1 transmittance =
  diagonal Jones t·I ⇒ apply the identical mask to all 3 planes.  (A future
  Jones-valued mask — VVC — is a per-pixel 2×2 in the transverse basis;
  Phase 3 applies it at the mask element via the same 3-plane hook.)

### 3a.1 Two defects found in the existing vector path (fix first)

1. **`PFFPROP` omits the Fresnel-integral output factors.**  Scalar `FFPROP`
   applies `1/(iλΔz)·dx1²` and the output-plane quadratic phase via
   `applyfac2` (propsub.F:2054-2056); `PFFPROP` (2273-2297) applies only the
   bare per-component FFT.  Harmless for a terminal intensity hop (global
   scale + curvature nobody reads), **fatal for a chain** — the output
   curvature is the input to the next leg.  Fix: replace `PFFPROP` with a
   K=1,3 loop over scalar `FFPROP` at the call site (propsub.F:1611), so
   vector and scalar far-field legs share one kernel by construction.
   **This changes existing vector-run intensity normalization** — intended;
   A/B the single-hop case and note it.
2. **The `ifPol` field assembly discards diffracted fields.**  At every
   physical-leg assembly the vector branch RELOADS the grid from rays
   (`WFElt(i,j,k)=RayE(k,iRay)`, propsub.F:1389-1395; the scalar-pol branch
   1397-1414 likewise rebuilds from |RayE| + CumL phase), whereas the
   non-pol branch (1456-1474) MULTIPLIES the existing `WFElt` by the
   incremental geometric phase `exp(i·TPL·(CumRayL−CumLStart))`, preserving
   prior diffraction.  So polarized multi-leg diffraction is broken by
   *bookkeeping*, independent of kernel coverage: leg 2 erases leg 1's
   diffraction.  Fix (vector mode): **seed once, then update** — on the
   FIRST physical-leg assembly of a trace, seed the 3 planes from `RayE`;
   on subsequent assemblies multiply all 3 planes by the per-ray
   incremental geometric phase (one factor shared by k=1..3) and maintain
   `CumLStart` exactly as the non-pol branch does.  A `LOGICAL`
   seeded-flag reset at trace start.
   Post-Phase-2 amendments (2026-07-26 review):
   - Correct seeded phases/amplitudes now depend on THREE engine fixes:
     `NaCmplx` (e2f680a) AND the coated-branch incident-medium +
     signed-cosine fixes (d137a97) — without the latter two, any coated
     surface upstream of the seed leg poisons the seeded field.  All are
     on pol-core.
   - `pol_set` now calls `modified_rx` (bef4598), so a pol-state change
     forces a full re-trace — key the seeded-flag reset off the same
     trace-start path and the lifecycle is automatic.
   - **Scope decision for the ifPol-SCALAR branch** (propsub.F:1397-1414,
     same erase-diffraction class): Tranche 1 fixes BOTH branches with the
     same incremental-multiply update.  Note the scalar-pol branch's
     |RayE| amplitude reload is only load-bearing when surfaces BETWEEN
     legs change ray amplitudes — excluded by the Tranche-1 mask-type
     validity condition, and handled by `J_run` in Tranche 2 — so the
     incremental-multiply form is exact within each tranche's stated
     scope.

### 3a.2 Tranche 1 — mask-type chains (the coronagraph case)

Valid when the elements BETWEEN physical legs are non-polarizing (Obscuring /
Reference / FocalPlane — true for pupil→FPM→Lyot→focal): between legs the
per-ray transfer is a scalar phase, so per-plane scalar update is exact.

- Call-site K=1,3 loops (kernels already single-plane; compute the leg's
  dx/z bookkeeping ONCE outside the loop): `NFPROP` at propsub.F:1592
  (PropType 2) and :1669 (5); `PPPROP` at :1650 (4) and :1691 (6); `SFPROP`
  at :1722/:1742 (7/8); `FRPROP` at :1804 (12); `NFPropDFT` at :1820/:1841
  (13/14); `FFPropDFT` at :1852 (15).  `DWF` scratch is per-call, reusable
  across k.
- `FFObscure` (:1635-1639): loop the call over the 3 planes under
  `ifVecDif3`.
- Assembly fix + far-field unification per 3a.1.
- `Ca2Int`/`Ca2Log`/`Ca2Gain` (utilsub.F) already sum 3 planes ✓; `ReGrid`
  is ray-side only, no change ✓.
- **Constraint to document:** vector mode repurposes the `mWF=3` planes as
  Ex/Ey/Ez — one wavefront only (no multi-WF/COMPOSE concurrently).

### 3a.3 Tranche 2 — Jones-through-chain (deferred until needed)

Chains with COATED/reflective surfaces between physical legs (IFO
recomb→detector with folds under `DO_NEARFIELD`; VVC layouts with an OAP
between masks) need the surface Jones applied to the GRID field mid-chain,
not only to rays.  Design: maintain a per-ray running 3×3 complex transfer
`J_run(3,3,mRay)` — reset to identity at each assembly, left-multiplied in
`elemsub.F` at every point where `RayE` is transformed (the C1 propagation
factor, s/p projections, RP/RS) so `RayE = J_run·E_prev_assembly` by
construction; the assembly update becomes
`WFElt(i,j,:) ← J_run(:,:,iRay)·WFElt(i,j,:)`.  ~7 MB at model 256.  This
subsumes the Tranche-1 scalar-phase update (J_run diagonal for mask-type
chains) — Tranche 1 remains as the fast path and the regression anchor.

### 3a.4 Validation ladder

1. **Energy conservation per leg**: Σ|E|² preserved by each vectorized
   kernel (unitary FFT core) to round-off.
2. **x-pol ≡ scalar equivalence**: near-normal-incidence train (Cass FF Rx),
   no coatings, Ex-only source: vector-chain Σk|Ek|² must equal the scalar
   run's intensity to round-off at every leg (Ey/Ez ≈ 0 ride along).
   Frame note (2026-07-26): "Ex-only" means the ENGINE launch frame —
   uniform (xGrid,yGrid) for collimated sources, the per-ray
   yray=unit(RayDir×xGrid) frame for point sources (ssrcray.inc; see the
   engine-pol section of macos_f90/CLAUDE.md).  Energy sums are
   frame-independent, so the gate is unaffected — but don't be surprised
   when per-component maps on point-source rigs aren't globally aligned.
3. **PROPER cross-checks re-run polarized**: pymacos proper_compare Phases
   2/3 (NF plane-to-plane + sphere legs) with `pol_set` on, x-pol — must
   reproduce the committed scalar↔PROPER tolerances (2.4e-14 … 4e-8 class).
   The strongest gate: vectorization must not perturb validated scalar
   physics.
4. **Chain closure**: CoroExample-class pupil→FPM→Lyot→focal, vector x-pol
   vs scalar baseline — contrast curves equal to round-off (masks are
   scalar; any difference is a Tranche-1 bug).
5. **Single-hop A/B** for the `PFFPROP`→`FFPROP`×3 unification (intended
   normalization change, documented magnitude).

### 3a.5 Execution split & sequencing

- Tranche 1 — kernel loops + FFObscure + tests 1/2/4/5, AND the two §3a.1
  fixes implemented to spec — is **Opus-lane item 1** (see the revised
  execution split above; CCMac when tokens allow, or an Opus session on the
  Linux box).  This section is the spec.  The §3a.1 diffs (assembly
  seed-once; far-field unification) get a **Fable line review** at the
  gate before merge; test 3 (polarized PROPER re-runs) and the ifx smoke
  run on the Linux box regardless.
- Slots after Phase 2 exposure work is stable (Jones-pupil machinery is the
  main consumer of correct chain phases) — Phase 2a/2b landed 2026-07-26,
  so Tranche 1 is now unblocked; must land before the Phase-3 VVC
  acceptance test.  Track A phase 1 (single far-field hop) is unaffected
  either way; Track A `DO_NEARFIELD` additionally wants Tranche 2
  (Fable lane, deferred until needed).

---

## Phase 4 — Spatially-variable coatings

Extend per-element uniform coatings to spatially varying, modeled on the `AmplMat`
ray→pixel projection (`elemsub.F:3348-3375`, `AmplSrfdx`/`nAmplMat`) — **with two
deliberate departures**:

**Departure 1 — representation.** Per-pixel full layer stacks run to gigabytes and
force a recursion per ray per surface. The physically-motivated variation is nearly
always smooth-and-low-order (radial grade) or **piecewise constant** (segment-to-
segment witness-sample variation across a segmented primary — the actual coronagraph
case). So: **primary form = per-pixel zone-index map + small coating table**; escape
hatch = per-pixel stacks, off by default; optionally precompute per-zone r_p/r_s at
the working wavelength to skip the recursion in the piecewise case.

**Departure 2 — interpolation.** Do **not** inherit `AmplMat`'s nearest-neighbor
lookup: quantization steps in a smooth phase map generate speckles exactly where we
can't afford them. Bilinear minimum; grid-resolution sensitivity an explicit test.

Mechanics otherwise parallel AmplMat: `CoatMapdx`, `ifCoatMap`/`iEltToCoatSrf`,
`CoatFile` keyword parallel to `AmplFile`/`AmplInit` (`surfsub.F:3278-3335`), SAVE
writer, API/bindings. In `Reflector`/`Refractor`, when `ifCoatMap(iElt)` is set, build
`nb_arr/kxb_arr` per-ray from the map — **the recursion is unchanged, only its inputs
become spatially varying.**

**Tests/examples:** radial grade → Jones pupil varies as designed; per-segment coating
→ Jones pupil shows the segment pattern; grid-resolution sweep shows no speckle-floor
dependence on map sampling.

---

## Validation document (deliverable — Dave 2026-07-26)

> **STATUS 2026-07-26: SKELETON + PHASE 0–2b + PHASE 3a EVIDENCE LANDED**
> (Opus lane item 7).  `macos/docs/macos-manual/polval/` — six sections
> (frontmatter/provenance, conventions + engine-fix provenance, Phase 1
> exposure, Phase 2a/2b Jones pupil, Phase 3a Tranche 1, gate index +
> coverage/gaps), six generated figures, `make polval` / `polval-pdf` /
> `polval-regen` / `polval-check`.  Driver:
> `MACOS_resources/mmacos/tools/pol_validation_report/`.
>
> **The no-hand-copied-numbers requirement is enforced, not aspirational.**
> Prose lives in `polval/*.md.in` and contains no numeric literals — only
> `@@TOKEN@@` placeholders.  `render_polval.py` resolves every token
> before writing anything (so a failed render leaves no half-updated
> `.md`); `tools/check_polval.py` runs as a prerequisite of `make polval`
> and refuses to build if a template was edited without re-rendering, a
> figure is newer than `numbers.json`, or a placeholder survives.  The
> report stamps engine + binding SHA, branch, model size, MATLAB and host.
> **The driver also asserts 19 gate thresholds mirroring the CI tests and
> ABORTS on a regression** — the report cannot document a broken gate as
> round-off (guard verified non-vacuous against degraded values for all
> three comparison operators and the missing-measurement path).
>
> Numbers the driver cannot produce (pymacos/ifx suite, PROPER compare,
> GMI, and the HISTORICAL pre-fix engine A/Bs) are in `external.json`
> with their producing command and capture date, and are labelled
> *(external, captured DATE)* in the report — not silently omitted, not
> presented as regenerated.
>
> **Finding, made while writing the unitarity section.** The
> transmission-uniformity gate reports `std/mean` ≈ 5.1e-14, but the true
> spread of that map is 6.1e-15 p-v (only 30 distinct doubles): `mean()`
> over 11484 points accumulates 5.1e-14 of summation error, *larger than
> the quantity being measured*.  The gate is a valid upper bound and is
> unchanged; the report now publishes the honest median-referenced spread
> and RMS alongside it, plus the summation floor itself, and the figure
> panel is referenced to the median (which was otherwise painting a
> spurious uniform offset across the whole pupil).  Left `tJonesPupil` as
> Fable reviewed it — flagging rather than silently editing a landed gate
> — but a tighter median-referenced assertion would be a small
> improvement if anyone touches it.
>
> **Still to append** (each phase's definition of done): Phase 2b
> low-order expansion, 2c, 2d, Phase 3, Tranche 2, Phase 4.  Ladder rungs
> 4/5/6 are enumerated as not-yet-climbed in the report's own
> coverage-and-gaps section rather than left implicit.

When the polarization work completes, we deliver a **validation report** with
PNG evidence from the validation suite — the reviewer-facing companion to the
test code.  Requirements:

- **Home & build:** `docs/macos-manual/polval/POLARIZATION_VALIDATION.md` +
  `polval/media/*.png`, built by the existing pandoc toolchain (`make polval`
  → docx/HTML/PDF alongside the manual and cmdref).  Committed media, like
  the manual's.
- **Regenerable by one command:** every figure and every quoted number comes
  from a driver (`mmacos/tools/pol_validation_report/`) that re-runs the
  validation cases and rewrites `media/` + a generated numbers include.  No
  hand-copied numbers — a stale figure must be impossible to ship silently.
  (Follow the demo-plot conventions: exact Strehl from OPD where used,
  autoscaled panels, non-obscuring legends.)
- **Structure mirrors the validation ladders.**  One evidence section per
  gate, each with: the claim, the figure(s), the measured number vs the
  analytic/reference truth, and the test that pins it in CI.
  - *Phase 0–2b (evidence exists TODAY):* unitarity-gate D/retardance maps
    on the conductor Cass (round-off); Fresnel-analytic fold — measured vs
    closed-form RS/RP and D across the AOI spread, residual panel (1e-14);
    2θ-symmetry diattenuation orientation map (quiver over the pupil) +
    azimuth-lock residual; double-pole vs local-sp retardance maps (the
    basis-artifact figure); conventions table + the three engine-fix
    provenance notes (NaCmplx, incident medium, signed cosine) with the
    tests that pin them.
  - *Phase 3a:* per-leg energy-conservation table; x-pol ≡ scalar residual
    maps per leg; polarized PROPER cross-check residuals at the committed
    tolerances (reuse the proper_compare artifact machinery); CoroExample
    chain-closure contrast curves (vector x-pol vs scalar overlay); the
    single-hop A/B for the normalization change.
  - *Phase 3:* crossed-polarizer extinction, QWP linear→circular Stokes
    check, VVC null depth vs retardance-error ε against the analytic
    (ε/2)² curve.
  - *Track A / pol-ifo:* visibility budget (predicted vs simulated fringe
    visibility), PSI polarization systematic vs input state, the AOI trade
    figures.
  - *Phase 4:* per-segment/radial-grade Jones-pupil maps + the
    grid-resolution speckle-floor sweep.
- **Lanes:** the report skeleton + the Phase 0–2b evidence sections and
  their figure driver are **Opus-lane work, addable now** (the gates and
  numbers are landed; specs above).  Each later phase appends its evidence
  section as part of its definition-of-done — this extends the standing
  docs rule (cmdref + manual + **validation-report section**).  Fable
  reviews the claims/interpretation before the document ships.

## Validation ladder (polarization ground truth)

Reconciling with PROPER only proves the scalar path is intact. In increasing effort:

1. **Analytic single-surface Fresnel Jones** — exact, available now; the unit-level
   check (both tracks).
2. **Unitarity** of a lossless Jones pupil (2a) — the single most diagnostic test.
3. **Rotational-symmetry 2θ invariant** — cheap, basis-sensitive.
4. **IFO self-measurement cross-check (NEW).** Put a known retarder/coating in the
   test arm of the simulated Twyman-Green and *measure* its retardance with the PSI
   pipeline (proven to 3e-5 pm); compare against the same coating's Fresnel analytic
   and against the Phase-2 Jones-pupil decomposition. Three independent paths to one
   number, entirely in-house.
5. **Published two-mirror polarization-aberration results** — compare form and
   magnitude of the diattenuation/retardance patterns (the 2b expansion exists to make
   this direct).
6. **One cross-check against a commercial PRT code (Polaris-M) if a JPL license is
   reachable** — the only item here that satisfies an outside reviewer on its own.

**Per-phase mechanics:** `run_mmacos_tests.sh fast` between edits; full suite +
pymacos pytest pre-commit; every `matlab -batch` ends `exit(0)`. **Every phase's
definition-of-done includes its user documentation**: cmdref entries (run
`make cmdref-regen` so new binding functions are cataloged, then fill their
NOTES blocks), the matching manual section (`docs/macos-manual/src/04` for
new Rx keywords, `src/05`–`06` for trace/diffraction behavior), **and its
evidence section in the validation report** (see the Validation document
deliverable below). The Phase 0–2b docs landed 2026-07-26; do not let the
gap re-open at Phase 3. Build both compilers
(`makems.sh release` + `makems.sh release gfortran`) — new fixed-form code passes
gfortran's stricter checks (LOGICAL `.eqv.`, ≤72-col, cpp `//` gotcha per
`mmacos/CLAUDE.md`). GMI regression stays 6/6 (its Rx are polarization-off; Phase-1
gating must be a no-op there).

---

## Files touched (representative, per phase)

- **Engine API:** `macos_f90/macos_api_mod.F90` (all phases).
- **Engine physics:** `macos_f90/elemsub.F` (Ph3 surface routines, Ph4 map lookup),
  `macos_f90/tracesub.F` + `propsub.F` (Ph3 dispatch, Ph4 per-ray inputs; possibly
  Ph3a vector legs), `macos_f90/elt_mod.F` (element constants + arrays),
  `macos_f90/surfsub.F` (Ph4 CoatFile loader), `macos_f90/srtrace.F` (Ph0/1 `ifPol`
  parameter, if in scope).
- **Parse/IO:** `macos_f90/msmacosio.inc` (Ex0Ey0 + Ph3/4 keywords), `macos_f90/
  iosub.inc` (defaults + SAVE writers; note the existing coating-SAVE audit comments),
  `macos_f90/macosio.F` (interactive MOD), `macos_f90/SAVE_KEYWORD_AUDIT.md`.
- **mmacos:** `src/gen_mex_wrappers.py` rerun (+ one Path-B helper in `mmacos_mex.F`),
  `src/+macos/*.m` + `Session.m`, `tests/t*.m`; `src/+macos/+design/Bench.m` +
  `twyman_green.m` (coating emission, Ph3 `add_polarizer`/`add_waveplate`).
- **pymacos:** `src/cmake/source/pymacos_f2py.f90` + `src/pymacos/macos.py` + tests.
- **Docs:** `macos_f90/CLAUDE.md` (new engine polarization/coating section — home of
  the conventions table), `mmacos/CLAUDE.md` / `pymacos/CLAUDE.md` cross-refs, a
  worked polarization example under the coronagraph design docs, and the bench_ifo_pol
  example README.

## Reuse (do not reinvent)

- `beam_set`/`beam_get` — API setter/getter template; `gen_mex_wrappers.py` Path A for
  everything mechanically regular.
- `cfield_get`/`cfield_apodize` — split-real/imag buffer + WFElt-multiply patterns.
- `ray_status_get` — per-ray `elt_mod` harvest precedent (for `rayfield_get`).
- `AmplMat`/`AmplInit`/`ChkRayTrans` — Phase 4 model (with the two departures).
- Existing s/p projection in `Reflector` — Phase 3 Jones-apply model.
- `modified_rx` dirty-the-trace convention for all setters.
- bench_ifo/bench_ifo_dm PSI machinery (proven to 3e-5 pm) — track A test vehicle.

---

## Open items to confirm with Dave

- **Reference-frame default** — double-pole proposed; confirm, and whether any JPL
  coronagraph polarization results we'll be compared against used a different frame.
- **Coating subsystem unification** — unify the two coating models or extend only the
  polarization-path one? (Phase 0 finding 5 informs.)
- **`PolSrc=`/`Ex0Ey0=` Rx semantics** — should a `.in` file be able to turn
  polarization ON (revive `PolSrc`), or only carry the source state (`Ex0Ey0`) with
  on/off left to API/CLI?
- **Phase 3 vs Phase 4 priority after 1+2** — PBS/QWP + VVC trades (Ph3) vs segmented
  coating non-uniformity (Ph4). Note track A gives Phase 3 a customer immediately
  (polarization phase-shifting IFO).
- ~~**Scope of Phase 3a**~~ — RESOLVED 2026-07-25: scoped in the Phase 3a
  section above (two tranches; two pre-existing defects found: PFFPROP's
  missing Fresnel output factors, and the ifPol assembly reloading RayE at
  every physical leg, which erases prior diffraction).  Tranche 1 closes
  mask-type chains (coronagraph); Tranche 2 (running per-ray Jones) covers
  chains with coated/reflective surfaces between legs.
