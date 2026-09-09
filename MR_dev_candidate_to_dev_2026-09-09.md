# Merge request: `dev-candidate` -> `dev`, both repos (engine first)

Prepared 2026-09-09 by CCL for Dave.  Order of operations (the standing
rule): merge **macos** first, then **MACOS_resources** -- the resources
veneers emit prescription keywords and call API routines that only the
engine at this tip parses.  Both merges dry-run CLEAN against
`origin/dev` (no conflicts; the four resources cherry-picks and the one
macos cherry-pick already on `dev` are patch-identical).

Tips: resources `dev-candidate` = `57a6ec0` (the sens-core merge; branch
deleted after it); macos `dev-candidate` = the commit carrying this file
(engine tip `cdf8636`, the Get_Values fix).

---

## 1. macos (engine): `dev-candidate` -> `dev`

**Title:** Engine: exit-pupil finder made frame-independent and surface-true, idempotent re-traces, OPD reference selectable, STOP hygiene, NS flow-of-light, parser hardening

**Body:**

Engine work since `dev` `a610358` (140 commits incl. records; the engine
commits are listed below), verified on the public corpus (mmacos full
suite, pymacos main + PROPER-compare) and on the JPL-private IRIS and
OPTIIX decks by CCMac (report: `BRIEF_ccmac_jpl_round2.md`).

Exit pupil (FEX/SXP):
- Frame-independent four-probe crossing = the MEDIAL pupil (`82d8148`);
  off-axis decks move by half the tangential/sagittal split (e5hex1 1%,
  jwst 3e-4); symmetric decks bit-identical.  Five fixture pins re-pinned
  with value + mechanism (`9948617`).
- EP radius intersects a CURVED next surface, not its tangent plane
  (`240b59b`): removes a pure off-axis defocus term (a4 = r0^2 h^2 /
  8 R_ep^2 R_next); 11 of 208 corpus decks move, all JWST-OTE focal spheres.
- Chief-ray axis is the default on all platforms (`448f5a8`); centroid opt-in.
- Rx-order guard dropped (it flagged the manual's own recipe) (`eb84095`).

Traces and stops:
- Re-traces are idempotent (`81d3308`): the source-frame
  re-orthogonalisation and the object-space STOP re-aim keep the incoming
  frame within round-off.  Removes a strict 2-cycle of a few ulp that set
  a finite-difference noise floor on the jwst zoom deck (elt-4 dw/dsurf
  column 1.19e-6 -> exactly 0).
- Element STOP preserves the source frame's handedness (`44fc362`): on
  left-handed decks an element stop mirrored the ray grid against the
  segment map (732 of 985 rays obscured on a segmented deck; now 2).
- STOP accepts Segment elements and multi-value one-line prompts
  (`eb84095`); stop state no longer outlives LOAD (`22b8d4b`).
- Non-sequential trace: flow-of-light root selection (`6bab7af`; already
  on dev as `a610358`).

OPD reference:
- `UseChfRay4OPD= Y` honoured (`9b4932a`, Luis Marchen's diagnosis): the
  chief-ray reference was unreachable on every path; `opd_ref_set/get`
  API (`4e538c5`); obscured chief is not a dead chief (`39c7015`).

Parser / IO:
- `GridFile=` read deferred until `nGridMat=` is known (silent 0x0 grid,
  `d7e64ff`); GridFile/AmplFile name buffer 24 -> 256 chars (`ab9fc45`);
  PSEG honours SegXgrid (`1ce1f5c`); `Get_Values` bounded by its buffer
  (36-byte over-read made `ArrWaveLen=`/`ArrIndRef=` lines parse
  nondeterministically) (`cdf8636`).
- CLI prints sub-prompts in non-readline builds; bundled readline
  auto-builds (`443140a`).

Verification:
- mmacos full suite green on the paired resources tip; pymacos main
  6694/0 + PROPER 26/26 (ifx) on the consolidated pair; GMI 6/6.
- JPL-private decks (CCMac, 2026-09-08/09): supervisor axis byte-identical
  or error-parity; engine axis classified (IRIS FEX -0.35 mm = half the
  T/S split; OPTIIX reset_xp=true differences = the `81d3308` floor).

Open, tracked in PLAN section 0 (not blocking): `trace(26)`->`trace(27)`
stale first OPD on the jwst zoom deck; IRIS `save_rx`->reload SIGSEGV
(not reproducible on the public corpus; two measurements requested);
IRIS `reset_xp=true` producing no rows under a stop (two numbers requested).

Records: `macos_f90/CLAUDE.md` (engine cheatsheet sections dated
2026-08-27 .. 2026-09-08), `REPORT_fex_probe_frame_independent.md`,
`REPORT_ep_dome_review.md`, `REPORT_sens_noise_center.md`.

---

## 2. MACOS_resources: `dev-candidate` -> `dev`

**Title:** mmacos: one sensitivity core behind the four multi supervisors, pupil-read ruling, pupil_find, configuration axis, benches + design-layer campaigns

**Body (depends on the macos merge above):**

Sensitivities (Luis rounds 2-3 + sens-core, verified on IRIS/OPTIIX):
- ONE multi supervisor (`private/dw_multi_core.m`) behind
  `dw_dx/dw_dz_zernike/dw_dsurf/dw_dgrid _multi` (`481fc24`); empty-OPD
  guard warns once per run, never errors (`10dc0eb`, Dave's ruling);
  read-surface preflight: the wavefront is read at the PUPIL, a
  pupil-less powered nElt-1 is refused with the add_pupil recipe
  (`a76a487`, `8c64323`, `3034fff`); `sens_core_ab` deck-agnostic A/B
  harness (`f144b3a`).
- Engine-truth powered-element discovery + loud explicit `elts`
  (`cdf7bc7`, on dev); `tEltTypeCoverage` (`a31fbdb`, on dev);
  `tNsFlowOfLight` (`9dbc19e`, on dev); `SENSITIVITY_TOOLS.md`
  (`5b6beb0`, on dev).
- Stop-enforced chief in supervisors + pupil_find + focal_surface;
  `fex_axis` option; `reset_xp_method='pupil_find'` with `pf_scope`
  (field default), per-config pristine state; configuration axis
  (`configs` option on the four multi drivers, tiled row order);
  element groups in `dw_dx_multi`; `dw_dx` Jacobian OPD in BaseUnits
  (baselines regenerated); `dw_dgrid` honours `elts`; dwdsurf removes
  piston/tip/tilt per Kr/Kc response.
- OPD conventions: `orient`/`sign` on `opd` + all drivers; `macos.opd_ref`;
  `macos.unload()`; `macos.noll_mode` replaces the JPL-internal
  `zernike_mode.m` (self-containment, Luis 2026-09-03).

Benches and campaigns (templates/40_benches, 80_end_to_end, challenges):
- Twyman-Green polarization-PSI DM gauge (v1 plate, v2 MacNeille cube,
  option-3 shallow plate at 96x96) and the Zernike-sensor twin, through
  S6 (multi-color); `dm_gauge_lib` one scoring implementation.
- e2e6m 6 m unobscured telescope + coronagraph back end; rodgers1/2/3 +
  afocal4 challenges; offset_imager; examples reorg (templates/ +
  challenges/) with the self-containment sweep.
- Polarization IFO slices (pol-ifo merged): rotating-analyzer PSI,
  BS-AOI trade, transmitting Refractors Extinc=0, material-axis gates.

Verification: mmacos full suite green at the tip (per-model-size
batches); pymacos 6694/0 + PROPER 26/26; JPL-private IRIS + OPTIIX per
`macos/BRIEF_ccmac_jpl_round2.md`.

Open (post-merge follow-ups, not blocking): emptyOPD guard keys on
nonzero samples rather than the ray-pass mask (pending CCMac's IRIS
numbers); remaining IRIS decks on the reduced protocol (round-2 brief
section 4).
