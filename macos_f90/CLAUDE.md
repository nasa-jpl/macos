# MACOS engine core — subsystem index (the cheatsheet is split into `notes/`)

> **Post-compaction / resume:** nested CLAUDE.md files and the notes are NOT auto-injected
> after compaction.  Before resuming engine work, re-read THIS index, then the `notes/<topic>.md`
> for the subsystem you touch (table below).  Root rules, conventions and the tree-wide
> index live in `../../CLAUDE.md`.  Nothing in `notes/` is paraphrased — engine gotchas are
> exact; add a new entry to the topic's notes file, and a one-line rule here only if it is
> every-session.

## Where the gotchas live

| Touching… | Read | Holds |
|---|---|---|
| any engine file | `notes/hygiene.md` | §0 closed latent bugs (macos_realloc, DFOURN, psiElt renorm, Get_Values over-read…), per-ray status (`SetRayFail`, the clip stamp, restarts), `DMPROD`, name collisions |
| `surfsub.F`, `elemsub.F`, `tracesub.F`, `propsub.F`, `srtrace.F`, `elt_mod.F` | `notes/surfaces.md` | FreeForm (14) + user input, GridSrf null frame + NS threading, grid index guard, `NSCount`, ZernType parsing + the ZernTypeL≠0 dispatch (ExtFringe wired), Conforming Reference is PASSIVE, GRIN model 5 only, root picks (ifLNsrf = nearest `RptElt`; forward = the VERTEX SHEET, `nFwdRootFlip`), chord-ruled gratings + OPL mod λ, grating apertures, medium-aware kernels, far-field evanescent cut |
| `FEX`/`SXP`/`XPS`/`STOP`, `LocalCoord`, `sourcsub.F` | `notes/pupil_fex_stop.md` | SPOT beam frame with an obscured chief, FEX axis = chief ray, FEX radius = the SURFACE of iElt+1, the four-probe medial pupil (five re-pinned tests), EP-radius rework + SXP, idempotent re-traces (`OrthoSrcFrame`, dead bands, `DeadBandNote`), STOP keeps handedness, STOP on a Segment + multi-value prompts |
| `SUBROUTINE OPD`, `opd_ref_*`, `opd_mask`, `spcOption`, any `OPDMat` consumer | `notes/opd_reference.md` | chief-ray default (session option at allocation; per-deck `UseChfRay4OPD`), `opdValid` — **never test `W ~= 0` for validity**, `zern_fit` on `opdValid`, the e5hex1 segment-piston measurement, SPOT/dcdx centring (`spcOption` at allocation) |
| `msmacosio.inc`, `iosub.inc`, `macosio.F`, `smacosio.F`, `validate_prescription.F90`, SAVE | `notes/rx_io.md` | the two parser chains, header keys match their FULL name + `HdrUnkNote`, `OptAsph=` after `VarElt=`, comments (block / whole-line / in-line), `GlassElt=` catalog at allocation, PolyObsVec deferral, `ReadRealsPad`, the validator, MOD empty-value guard, a trace never writes the Rx (`SrfTypeTr`, `OrthoObsFrameLoad`), test prescriptions |
| `design_optim.F`, `nls.F`, `macos_ops.F`, `calib_*`, `Opt*` keywords | `notes/calib.md` | CALIB wrappers, targets (WFE_ZMODE works, OPL removed), the LM stopping rule and its derivation (K=10 floor, OptTol abs AND 1e-6 rel), beam rows on any target, asphere DOFs + step, the failure path and the leaked SAVE arrays, `calib_setup_vars` |
| `macosio.F` prompts, `mhist.c`, `macos_cmd_loop.inc` | `notes/cli.md` | readline sub-prompt cache (`set_sub_prompt`; `LEN_TRIM` not `ICLEN`), non-readline builds, IACCEPT_S reprompt sweep (abort to label 1, never re-enter), LOG cleanup |
| `macos_api_mod.F90`, pymacos/mmacos, coronagraph fixtures | `notes/api_bindings.md` | the wrapper layer (`ray_status_get`, buffer-wrapper template), pymacos CodeV/PROPER regression phases, masks must be `Element= Obscuring` |
| `ifPol`, `RayE`, coating branches, `PolElt`, `ifVecDif3`, `CPROPAGATE` legs | `notes/polarization.md` | conventions table (exp(+iωt), N = n − iκ, coating thickness in waves), Model A vs B coatings, the coated-branch fixes, the r_p sign fix and its gate-blindness lesson, radiometry (ONE factor, never per interface), the external anchor, overcoat quarter-wave reversal, Jones basis, TrPolarizer/WavePlate (project the MATERIAL axis), **two dispatch chains**, Tranche-1 vector propagation and the seed phase convention |
| plots / display | `giza/CLAUDE.md` | (unchanged) |
| constrained optimisation | `slsqp/CLAUDE.md` | (unchanged) |

## Every-session rules (each is argued in its notes file)

- **A trace reads the Rx and writes TRACE state.**  A modified copy of an Rx value is a `*Tr`
  array or a local, never the Rx array (`notes/rx_io.md`).
- **User options live at allocation, never in `*_init_vars`** (glass catalog, `spcOption`,
  `LOpdRefChief`, the five session modes).  The test: can a user set it with no deck loaded and
  expect it to stay?  (`notes/opd_reference.md`, `notes/surfaces.md`.)
- **Zero is not "no ray".**  Pupil validity is `opdValid` / `opd_mask`, never `W ~= 0`
  (`notes/opd_reference.md`).
- **Two element dispatch chains** — `tracesub.F` CTRACE and `propsub.F` CPROPAGATE.  A new
  element or field effect goes in both, or it is invisible to `intensity`/`complex_field`
  (`notes/polarization.md`).
- **Guards warn and count; a new refusal needs Dave.**  Precedent for a refusal: a crash or a
  corrupted trace becomes a clean NAMED refusal (`notes/rx_io.md`).
- **Measure before you guard** — grep the value's later assignments (the `nGridpts` clamp
  already existed; a parse-time refusal broke 195 tests).
- **Header keys match their FULL name**; a shorter `LCMP` captures look-alikes
  (`notes/rx_io.md`).  Element-chain matches are not audited.
- **A Zernike surface traces only with `ZernTypeL ≠ 0`** (ANSI default on every
  Zernike-carrying surface; ExtFringe wired 2026-10-12) (`notes/surfaces.md`).
- **Masks for diffraction are `Element= Obscuring`**; a Reference with obscurations clips rays
  only (`notes/api_bindings.md`).
- **Never declare a dummy `A(1)`; write `A(*)`** — the gfortran debug build bounds-checks it
  (`notes/surfaces.md`, grid section).
- **Any surface / element change that can move a trace gets the corpus A/B** before its
  commit, gated from the clean worktree `~/dev/macos_cli_base`, every mover named in
  `cli_tests/records/ab_*_moved_<sha>.txt`.
- **A gate's analytic is written from the textbook, never transcribed from the engine**
  (the r_p sign lesson), and every vignetting element's gate puts one ray OUTSIDE its aperture.
