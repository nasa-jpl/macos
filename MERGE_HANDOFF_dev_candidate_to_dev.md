# Push request: `dev-candidate` -> `dev`, both repos (prepared 2026-09-29)

For Dave to review and push, then open the two PRs.  Engine first, resources
second (the branch model's merge-ordering rule: a resources veneer promoted
ahead of its engine fails at prescription-load time, not build time).

## Preflight -- done, all local, nothing pushed

| | macos | MACOS_resources |
|---|---|---|
| dev-candidate vs origin/dev | **336 commits ahead** (2026-08-06 .. 09-29) | **675 commits ahead** (2026-07-25 .. 09-29) |
| dev's stragglers merged in | `5b60927` -- Andy's 09-05 NS flow-of-light commit; **content no-op** (identical patch already here as `6bab7af`) | `7f9490e` -- Andy's four 09-05 sensitivity commits; **content no-op** except `SENSITIVITY_TOOLS.md` add/add, resolved to ours (dev's text + 117 lines inserted at l.43) |
| verified by | `git diff --stat <pre-merge> HEAD` empty | same, empty |
| local commits awaiting push | 8 (`223a6ff` .. `5b60927`, incl. this file) | 7 (incl. `49364b3` TO's descent tool, `0d0edd8` tRxBlockComment, `7f9490e` merge) |
| gates | block-comment fixture, both compilers, gfortran==ifx byte-identical; eac5mono.in loads | **fast suite: RESULT PENDING (fill in)** |

Also in this batch: `DEV_FILES.md` + `release-exclude.txt` re-audited
(`de0f3c8`).  That list governs the LATER promotion `dev -> main`; the PR to
`dev` carries every file, developer-facing ones included, by design.

Push (Dave):
```
cd ~/dev/macos           && git push origin dev-candidate
cd ~/dev/MACOS_resources && git push origin dev-candidate
```
Open (either of us, after the push; engine first):
```
cd ~/dev/macos           && gh pr create --base dev --head dev-candidate --title "<title below>" --body-file <this file's macos section>
cd ~/dev/MACOS_resources && gh pr create --base dev --head dev-candidate --title "<title below>" --body-file <resources section>
```

---

## PR 1 -- nasa-jpl/macos: `dev-candidate` -> `dev`

**Title:** Engine: FEX/SXP exit-pupil rework, OPD reference, NS flow-of-light, Rx comment fixes, CLI/SAVE hardening (Aug-Sep 2026)

**Body:**

336 commits; 24 touch the engine (`macos_f90/`, build).  The rest are
developer-facing files (briefs, runbooks, reports, the gauge decks under
`demo_session/`, slice/plan state) that `dev` keeps and the `dev -> main`
strip list removes -- see `DEV_FILES.md`, re-audited in this batch.

### Ray trace / geometry
- **FEX/SXP exit-pupil rework** -- EP radius is the chief-ray distance to
  `iElt+1`'s SURFACE, not its tangent plane (`240b59b`; curved focal
  surfaces carried pure defocus in the off-axis OPD, zero on axis, which is
  why no gate saw it for a year); the probe is frame-independent
  (`82d8148`: four differential chief rays, medial pupil -- the legacy
  single probe's answer depended on the sign and azimuth of xGrid); the
  pupil-sphere axis defaults to the CHIEF RAY on every platform
  (`448f5a8`, Dave's ruling; centroid is opt-in); the Rx-order guard that
  fired on most legacy decks is gone (`eb84095`).
- **Re-traces are idempotent** (`81d3308`): `OrthoSrcFrame` keeps the
  incoming source frame when the re-orthogonalised one differs by round-off
  only, and the STOP-OBJ aim gets a dead band -- the 2-cycle ulp
  oscillation that put a 1e-6 speckle floor under every finite-difference
  sensitivity column.
- **Element STOP preserves the source frame's handedness** (`44fc362`):
  an element stop on a left-handed deck mirrored the ray grid against
  `EltToSegMap` and obscured 732/985 rays on segmented sources.
- **NS trace: flow-of-light root selection** for non-sequential candidates
  (`6bab7af`, = Andy's `a610358` on dev).
- STOP accepts `Segment` elements; multi-value prompts gather across
  tokens (`eb84095`); stop state no longer outlives LOAD (`22b8d4b`).

### OPD reference
- `UseChfRay4OPD= Y` was unreachable -- the keyword had only an `N`
  branch, and the cmd-loop's `.TRUE.` was reset by `MBFile6`'s
  reinitialise before it could matter (`9b4932a`, Luis's diagnosis);
  `opd_ref_set/get` in `macos_api_mod` (`4e538c5`); an OBSCURED chief
  still serves as the reference -- the gate is `LRayOK`, geometric, not
  `LRayPass` (`39c7015`).

### Prescription I/O
- **Block comments** `/* ... */` / `CommentBegin..End` (`223a6ff`, Scott's
  report): the parser always handled them; the Phase-1 validator ran first
  and refused them, and SAVE lost them.  Both fixed; validator also treats
  CR as whitespace, judges `Key= % note` as empty, and refuses an
  unterminated block.  Whole-line `%` comments already round-tripped;
  in-line ones are deliberately not preserved.
- `Get_Values` bounded by its buffer, not the command length (`cdf8636`):
  `ArrWaveLen=`/`ArrIndRef=` parsed nondeterministically (a fixture loaded
  3 of 5 times).
- `GridFile=`/`AmplFile=` name buffer 24 -> 256 chars (`ab9fc45`); the
  grid read deferred until `nGridMat=` is known (`d7e64ff`, silent 0x0
  grid); phantom grids (`nGridMat>0`, no data) no longer crash
  save->reload (`cda178e`).
- `Link=` survives SAVE (`d6b17c8`); PSEG honours `SegXgrid` (`1ce1f5c`).

### CLI / build
- Non-readline builds print sub-prompts; cmake auto-builds the bundled
  readline when its `.a` is missing (`443140a`).
- cmdref: CHIefray/CENTRoid convention table, FEXit/SXP entries
  (`af749d5`, `f33d1e2`).

### Gates
mmacos fast suite on the exact tree (see PR 2) is the gate for this pair.
pymacos: `pymacosf90.so` rebuilt against this tree (2026-09-29) and the
block-comment fixture smoke-tested through it (load, KrElt inert, SAVE ->
load -> SAVE identical, blocks kept) -- so the fix is proven on all three
surfaces.  Its full suite (6601 + PROPER-compare) was last green on
2026-09-14's tree and has NOT been re-run on this one.  The FEX re-pins carry the mechanism, not a tolerance bump
(`9948617`).  Block
comments: `ZGD_test_files/tst_block_comment.in` byte-identical SAVE ->
load -> SAVE on gfortran and ifx, gfortran == ifx.

### Not in this PR (open, tracked in CURRENT_SLICE)
The redo-bench descent fix (aperture pulled in ~3 actuators + per-column
regularisation) -- diagnosed, not yet implemented.

---

## PR 2 -- nasa-jpl/MACOS_resources: `dev-candidate` -> `dev`

**Title:** mmacos: template reorganisation, sensitivity core, pupil_find, the DM-gauge bench campaign (TG96 / ZWFS / PDI), design-layer benches (Jul-Sep 2026)

**Body:**

675 commits: 93 library (`mmacos/src`), 60 design layer, 422 templates,
120 tests, 14 pymacos, 14 docs.  Requires the engine PR (macos
`dev-candidate`) merged first -- the sensitivity supervisors and the bench
runners call `opd_ref_set`, the FEX probe and the STOP fixes it carries.

### Library (`+macos`)
- **Sensitivity core** -- ONE multi supervisor behind the four `dw_d*_multi`
  fronts (`481fc24`); engine-truth powered-element discovery and loud
  explicit `elts` (`8e7363d`, Luis round 3); `wf_elt_auto` reads at the
  EXIT PUPIL and errors on a pupil-less powered element (`8c64323`,
  `3034fff`); the empty-OPD guard warns once per run, never errors
  (`10dc0eb`, Dave's ruling); the OPD reference is carried by every driver
  (`fe6c6bd`); element GROUPS (`acd55cb`, `de8a8cc`) with the grouped-dwdx
  speedup and the non-ANSI Zernike grid basis (`25a5a25`, engine-exact,
  `4772c83`); the `configs` axis on all four multi drivers (`4a6659c`..);
  `dw_dgrid` honours `elts` (`3346d04`); `dw_dx` emits the Jacobian's OPD
  in BaseUnits with baselines regenerated (`f5d648f`) and the units adapter
  in `run_compare`/`run_simulator` (`c0c6bad`).
- **pupil_find / fit_chief / focal_surface** -- the cone-fit sphere as a
  supervisor XP method (`e3d08ea`), `fit_chief` written vertex
  (`7006f85`), per-(config, field) mini-cone placement (`be2c505`),
  honours a deck-declared object-space ApStop (`d51c70a`), the
  stop-enforced chief (`2dd62ae`); `focal_surface` measures the best-focus
  surface (`07b2edd`).
- `macos.opd_ref` + Session `opd()` forwarding (`b95bbfc`); `macos.unload()`
  (`5960add`); `macos.noll_mode` replaces the JPL-internal
  `zernike_mode.m` -- self-containment (`98d3320`); `orient`/`sign` on
  `opd` and every driver (`0e556e3`); `append_rx` (`03bd70f`);
  `prop_layout` (`30bc7bb`); `view_rx` multi-field bundles and glass
  solids (`fe4a5d5`, `6545942`); dwd* figures sized first, paginated
  second (`eac0fc8`).
- Polarization bindings (Phase 1-3a: Jones pupil, pol maps, vector chain,
  plane-selectable complex field) -- `880c217` .. `9f2eed4`.

### Design layer
- `twyman_green` gains `optics` `'lens'|'oap'` (`2e56522`); `add_oap`
  geometry fixed and bodies drawn at the OAP pole, not the parent vertex
  (`02c291c`, `9a8d699`); `psri_bench` (`19d9514`); the real MacNeille PBS
  cube (`5cd2112`); transmitting Refractors default transparent
  (`cafc53c`); `aperture_full_field` in the element frame (`6703a38`).

### Templates -- the reorganisation and the arcs
Stage 1 moved the corpus into a ladder (`8cd731a` .. `8feaee0`):
`10_telescopes`, `20_segmentation`, `30_instruments`, `40_benches`,
`50_sensitivities`, `80_end_to_end`, `90_polarization`, plus
`challenges/`.  On it: the **DM-gauge bench campaign** (`40_benches`:
`tg_psi_dm96_oap` 123 commits, `zwfs_dm96` 89, `pdi_dm96` 37,
`dm_gauge_lib` 37 -- interferometer, Zernike, vector-Zernike and
point-diffraction gauges scored on one bench, the redo bench collimated for
real, the descent-stall diagnosis tool `tg96_ring_analysis`); the
**coronagraph test bench** (`bench_ctb`, 40); **e2e6m** rounds 1-2 (53);
**offset_imager** + `oi_clear`/`oi_score` (46); **zoom_5x5** (21);
**afocal4** / rodgers challenges.

### Tests
120 commits; new classes include `tRxBlockComment`, `tLinkSave`,
`tZernikeGridBasis`, `tDwDxGroups`, `tPupilFindMethod`, `tFocalSurface`,
`tStopReload`, `tDmgLoop`, `tPolElement`, `tPolRadiometric`,
`tPolExternal`, `tVecChain`, `tJonesPupil`, `tEltTypeCoverage`,
`tNsFlowOfLight`, `tBench`, `tOpdRef`.  Gate for this PR: the fast suite
on the exact tree (result above).

### pymacos
Polarization bindings + tests (`80ab9ad`), `spot()` return-code check
(`8c172ff`), plane-selectable `complex_field` (`9f2eed4`), `fex_axis`
help (`25eae0b`).  `pymacosf90.so` rebuilt against this tree 2026-09-29
and smoke-tested (block-comment fixture); the full pytest suite not re-run
here.

### Not in this PR
The descent fixes (aperture / per-column lambda) and the three open bench
decisions (D_BS_CMP, field-lens conic, overcoat quarter-wave).
