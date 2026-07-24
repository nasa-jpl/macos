# Constrained optimization (SLSQP / NPSOL) — subsystem cheatsheet

> **Post-compaction / resume:** if you are resuming any constrained-optimization work (`design_slsqp_optim.F`, SLSQP/NPSOL dispatch, the pre-scaling fix), re-read THIS
> file — nested CLAUDE.md files are NOT auto-injected after compaction;
> they reload only when CC next reads a file in this directory. Root rules,
> conventions, and the subsystem index live in `../../CLAUDE.md`.

*Sections below are lifted verbatim from the former monolithic root
CLAUDE.md. Move text, don't paraphrase — engine gotchas are exact.*

---

## Constrained optimization (SLSQP default, NPSOL opt-in)
- Default back end: `design_slsqp_optim.F` + vendored Kraft SLSQP under
  `macos_f90/slsqp/` (BSD; see `slsqp/LICENSE.txt`).  Same external
  signature as `np_optim_dvr` so call sites in `macos_cmd_loop.inc`
  (start_optim ~3661, restore ~3964) don't change.
- Dispatch (in `macos_cmd_loop.inc`): `#ifdef USE_NPSOL` calls
  `np_optim_dvr`; `#else` calls `slsqp_optim_dvr`.  Both drivers exist
  under USE_NPSOL=ON so the same binary can be A/B'd via prescription.
- MACOS's "constrained" optimization is bound-constrained only
  (`nclin=ncnln=0` historically in `design_cons_optim.F`); SLSQP runs
  in its simplest regime — smooth objective + box bounds.
- **Numerical gotcha: variable pre-scaling.**  SLSQP's QP solver
  underflows when `|grad|/|bound|` is large (e.g. macos's native
  `dtt=1e-9` FD step gives ~1e5 gradient against ~5e-4 tip/tilt
  bounds → 1e9 ratio → zero step on first call, mode=0 silent
  exit).  `slsqp_optim_dvr` rescales each DOF by
  `s_i = 1/max(|bl_i|, |bu_i|, eps)` before every `slsqp()` call and
  un-rescales `aparams` inside the funcobj evaluation — SLSQP is
  mathematically scale-invariant so this is a pure numerical-health
  fix.  Comment block above the slsqp call has the derivation.
  Do not bump `dtt`: ray-trace linearity at large FD steps was the
  reason the alternate fix was rejected.
- A/B verification reference: `ZGD_test_files/opt_example_constrained.in`
  (Elt 7 TIP+TILT, ±0.5 mrad bounds, WFE target, OptMaxItrs=200).
  NPSOL and SLSQP agree to 11+ digits on final WFE.
- Open work: setbeam.inc still uses `npoptn` + `npsol` directly
  (Phase 3.1 — port deferred); `m.calib_set_algorithm(...)` setter
  in pymacos/mmacos bindings (Phase 5).
- See `macos_f90/slsqp/README.md` for the full NPSOL→SLSQP mapping,
  reverse-comm protocol, and convergence story.

