# BRIEF — OPD reference: make the engine match the manual (chief ray), keep the mean as an option

For **TO**, 2026-10-10, Dave's ruling (PLAN_CONSOLIDATION §4 item 2): "Match the manual —
but highlight the option to switch to the mean in the manual and cmd ref."  Same
contract as the phase-1 items: gate RED first on the pre-fix binary, both compilers,
corpus A/B mandatory (this moves every per-ray OPD map), fast suite, CC gates from
the clean worktree, one commit pair, nothing pushed.

## The facts

- The manual (§5 Ray-Trace Analysis) and the cmdref's OPD entry define the OPD as
  "the difference between the optical pathlength of that ray and the chief ray."
- The engine has ALWAYS subtracted the aperture MEAN (`DAvgl`): `macos_cmd_loop.inc:378`
  sets `LUseChfRayIfOK=.TRUE.` on every LOAD, but `MBFile6`'s first act is
  `reinitialise_variables` → `ray_mod_init_vars` → `LUseChfRayIfOK=.FALSE.`
  (`traceutil_mod.F:278`), so the chief branch of `OPD` (tracesub.F) and its `propsub`
  twin never ran unless the Rx keyword `UseChfRay4OPD= Y` (msmacosio.inc:474, made to
  work 2026-08) or the API `opd_ref_set(OK,.TRUE.)` (macos_api_mod.F90:4139) asked
  for it after the load.  Full history: `macos_f90/CLAUDE.md` "OPD reference:
  `LUseChfRayIfOK` is FALSE on every path".
- Why it matters: on a segmented pupil the mean reference prints a spurious piston on
  every untouched segment when one segment moves (e5hex1: 2.85e-6, 17 % of the poked
  segment's peak); the chief reference gives exactly 0.  Monolithic pupils differ by a
  piston only.  An OBSCURED chief ray still serves (the branch gates on `LRayOK(1)`,
  geometric, not on vignetting) — no deck in the corpus has a geometrically dead chief.
- `design_optim.F:659-664/880` deliberately forces the flag FALSE inside optimisation
  loops and restores it after: leave that alone, and make sure "restore" now restores
  the NEW default.

## The change

1. **Default = chief ray.**  The user option lives where user options live (rule 6:
   NOT in `ray_mod_init_vars` / `elt_mod_init_vars`): set the default once at
   allocation (`ray_mod_init` or the start-up path, the `spcOption` / glass-catalog
   pattern), remove the per-load reset, delete the dead `.TRUE.` at
   `macos_cmd_loop.inc:378` and its "dead in place" note.  The Rx keyword keeps both
   values: `UseChfRay4OPD= Y` (now the default) and `UseChfRay4OPD= N` (the mean);
   `opd_ref_set` unchanged.  Add a CLI command pair — `OPDREF CHIEF | MEAN` (or extend
   an existing toggle in the Behavior Toggles group; your call, say which) — so the
   interactive user can switch without editing the deck, matching how `CHIEFRAY` /
   `CENTROID` select the FEX axis.
2. **Every OPD run PRINTS its reference** — one line, `OPD reference: CHIEF RAY
   (default)` / `OPD reference: APERTURE MEAN (UseChfRay4OPD= N)` — the FEX-axis
   precedent.  The chief-ray-dead fallback (the branch already falls back to the mean
   when `LRayOK(1)` is false) prints `OPD reference: APERTURE MEAN (chief ray lost)`.
3. **Fallback counted, not silent:** if the chief ray is geometrically lost the mean is
   used and the line above says so; no error (rule 5).

## Gates (red first on the pre-fix binary)

- `tOpdRef` (mmacos, exists, 8 tests): extend — the DEFAULT after a plain `load_rx`
  is the chief reference (e5hex1 single-segment poke: untouched segments at 0, not
  2.85e-6); `UseChfRay4OPD= N` in a deck twin restores the mean; `opd_ref_set` still
  overrides either way; the obscured-chief case (`Rx_Cass_FarField`) uses the chief.
- CLI journal `cli_tests/cmd/OPD_reference.jou`: `opd` prints `OPD reference: CHIEF
  RAY (default)`; after `opdref mean` (or the chosen command) it prints the MEAN line
  and the `Average OPD` becomes ~0; the `UseChfRay4OPD= N` twin deck prints MEAN.
- pymacos twin of the default test.
- **Corpus A/B (rule 3), the real gate:** `cli_tests/run_cli_tests.sh corpus` from the
  clean worktree at pre and post.  EXPECTED: `nPass` and `lost` unchanged on every
  deck; `RMS OPD error` and `P-V` — say in the record whether the CLI's RMS is
  mean-removed (then unchanged except round-off) or raw (then it moves by the piston
  on every deck with a non-zero chief-vs-mean offset); `Average OPD` moves everywhere.
  Every moved deck in a committed record (`cli_tests/records/ab_opdref_moved_<sha>.txt`,
  the `ab_6a` form); the CC gate will PASS on status/pass/lost and list the movers.
- Re-pins: every mmacos / pymacos test that pins a per-ray OPD MAP (not an rms) moves
  by a piston — re-pin each with the mechanism and the legacy value in its comment
  (the FEX-rework convention), never a tolerance bump; name them in the commit.
- **The bindings' drivers default to the MEAN today** (`opts.opd_ref = 'mean'` in
  `dw_dx`, `dw_dx_multi`, `dw_dz_zernike(_multi)`, `dw_dsurf(_multi)`, `dw_dgrid_multi`,
  set on the session after every load): change the default to `'chief'` in all of
  them (resources, same commit pair) so the bindings, the engine and the manual agree.
  EXPECTED and to be measured: the segmented harvests (e5hex1 / e5pie / e2e6m
  families) move by exactly the piston the mean reference was attributing to the
  untouched segments; the monolithic ones by round-off.  Re-pin `tDwDx`, `tDwDgrid`
  and every committed harvest the fast suite compares, with the mechanism and the
  legacy value in the comment.  The per-load `session.opd_ref(...)` call stays only
  for an explicit `'mean'` once the default lives at allocation.

## Documentation (CC/CCMac after the engine lands, but the words are the engine's)

- Manual §5: keep the chief-ray sentence; add one paragraph: "By default the OPD is
  referenced to the chief ray.  The alternative — the mean pathlength over every valid
  ray in the aperture — is selected per prescription with `UseChfRay4OPD= N`, per
  session with `OPDREF MEAN`, or from the bindings with `opd_ref_set`; every OPD run
  prints which reference it used.  On a segmented pupil the mean reference attributes
  one segment's motion as a piston to every other segment; use the chief ray.  If the
  chief ray is geometrically lost the engine falls back to the mean and says so."
- cmdref OPD entry: the same, plus the printed line; a new entry for the command;
  the `UseChfRay4OPD` keyword in the prescription-keyword list with both values.
- `macos_f90/CLAUDE.md`: the existing section rewritten to the new state (one
  paragraph: default chief at allocation, the mean opt-in, the print, the fallback;
  keep the segmented-pupil measurement and the "obscured chief still serves" fact).
- PLAN.md §0.x and PLAN_CONSOLIDATION §1.1 / §4 closed with the SHAs.

## Rules that bite

Rule 6 is the whole trap here: the last three "defaults" that went wrong (`spcOption`,
the glass catalog, this flag) all lived in a per-load reset.  Set once at allocation.
≤ 72 columns.  Guards warn.  No `--amend`; push on Dave's word.
