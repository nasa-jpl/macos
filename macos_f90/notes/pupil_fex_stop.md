# MACOS engine notes — Pupils, FEX/SXP, STOP, source frames, idempotent re-traces

> **Read this file before touching:** `FEX`/`SXP`/`XPS`/`STOP`, `LocalCoord`, `sourcsub.F`, `OrthoSrcFrame`, `UpdSrcGrid`.  These notes are the engine record
> (what moved, the measurements, the gates); they are NOT auto-injected — the index is
> `macos_f90/CLAUDE.md`, the root rules `../../CLAUDE.md`.  Text is moved verbatim from the
> former monolithic cheatsheet (2026-10-15 split); add new entries here, not to the index.

## SPOT BEAM reference frame: obscured chief ray (tracesub.F LocalCoord)
*(User rule -- an obscured but geometrically defined chief ray still defines the Beam frame -- is Manual §5 SPOT.  Mechanism below.)*
- `SPOT ... beam` (and `WINDOW/PLOCATE ... beam`) build the beam coordinate
  frame from a **chief-ray trace** in `LocalCoord` (tracesub.F ~2612): a
  nominal chief ray + two differential rays (x/y translate/rotate), each via
  `CRTrace` (single ray, `SetSourceRayGrid(nMinPts=3)`).  Each `CRTrace`
  result was gated on `IF (nBadRays.NE.0)` → "chief ray becomes undefined" →
  abort.  **The on-axis chief ray of any centrally-obscured telescope is
  VIGNETTED**, so `nBadRays≥1` from obscuration alone, and beam-frame SPOT
  failed at **every** element of every obscured Rx (reproduced identically in
  the CLI; TOUT/TELT work because they use the stored `Tout`/`TElt` frame, no
  chief-ray trace).  Obscuration sets only the flux flag (`L1`/`LRayPass`),
  NOT the geometric intersection `RayPos` (CTRACE comment ~4189), so the frame
  is well-defined.  Fix: the 5 `LocalCoord` failure checks now read
  `nBadRays.NE.0 .AND. RayStatus(1).NE.RayStat_Obscured` — tolerate a
  vignetted-but-geometrically-valid chief ray, fail only on a genuine
  geometric miss/bracket.  Do **not** try to fix this by toggling `iObsOpt`:
  the obscuration `nBadRays` bump (tracesub.F ~4364) is unconditional
  (`obs_set 0` does not suppress it), and it MUST stay that way so WARN's
  per-category loss accounting fires.  mmacos regression: `tSpot` (obscured
  Cass, beam ray count == tout).

## FEX/SXP axis convention: CHIEF-RAY DEFAULT on all platforms (Dave 2026-08-27)
*(Documented: Manual §5 FEX (axis paragraph), §3 "Session options and LOAD"; cmdref CENTRoid / CHIefray / FEXit / SXP.  Below: the sites, the history and the gate pattern.)*
`ifCentroid` governs ONLY the FEX/SXP pupil-sphere AXIS (nothing else
consumes it): chief (`psi=-cr1dir`) vs centroid (`psi=unit(centroid@
iElt+1 - CrossPt)`).  On obscured/segmented/comatic beams the centroid
walks off the chief -> tilted sphere -> PURE TIP/TILT (+piston vs a
chief OPD reference) in the OPD -- frame terms, zero aberration (Luis's
2026-08-27 report; same physics as the pupil_find fit_chief ruling).
Historic trap: the CLI defaulted to CENTROID (`ifCentroid=.TRUE.` in 5
init sites) while the api (`xp_fnd(mode)`, save/set/restore per call)
defaulted to chief -- t/t/p appeared interactively but never in
supervisor runs.  FLIPPED to `.FALSE.` at all 5 sites (macos_mod,
macos_init.inc, smacos_glass.inc, macos_ops.F, macos_cmd_loop.inc).
Opt-in: `CENTROID` cmd (match-5; `CHIEFRAY` restores), api mode 0, and
the Rx keyword `FEXCentroid= N` remains a per-deck VETO.  Every FEX/SXP
run now PRINTS its axis.  SXP radius = plane intersection with
iElt+1's vertex/normal -- TYPE-AGNOSTIC (mask/Reference/Obscuring at
iElt+1 handled right despite the "FP" name); FEX and SXP are
geometrically identical, FEX adds the guards.  Gate: pty-driven CLI
(scratchpad fex_cli_gate.py pattern -- readline needs a tty; the
model-size prompt comes FIRST).  cmdref FEXit/SXP/CHIefray/CENTRoid
entries carry the convention table.

## FEX/SXP radius on a CURVED iElt+1: intersect the SURFACE (2026-08-28)
*(Documented: Manual §5 FEX, the KrElt bullet.  Below: the engine record -- the four load-bearing facts, the measurements, the blast radius.)*

The 2026-07-03 rework below made the EP radius the chief-ray distance
from the EP crossing to element iElt+1's **tangent PLANE**
(`VptElt`/`psiElt`).  Its note says "type-agnostic" -- right about the
element TYPE, but it drops CURVATURE.  When iElt+1 is a real focal
SURFACE and the chief lands at image height `h` off its vertex, the
plane sits `h**2/(2*Kr)` beyond the surface, and the EP sphere carries
that as **pure DEFOCUS in the OFF-AXIS OPD**: `a4 = r0^2 h^2 /
(8 R_ep^2 R_next)`.  **Identically zero on axis**, which is why every
gate missed it for a year (found running `BRIEF_wnom_cli_ab.md`;
`macos/REPORT_wnom_cli_ab.md` + `REPORT_focal_surface.md`).

`FEXConicLeg` (tracesub.F, module procedure, called by BOTH `FEX` and
`SXP`) refines the plane leg to the conic intersection.  Four things
about it are load-bearing:

1. **It does NOT go through ConSrf.**  ConSrf picks its quadratic root
   by `|L**2-mpr|` PROXIMITY with no flow-of-light sense -- the trap
   that produced the CENTROID branch's far-root garbage.  Here the root
   NEAREST the plane leg is taken, which is well posed: the sag
   correction is `h**2/2R` while the roots are ~`2R` apart.
2. **It gates on SrfType, NOT on Kr alone.**  A `Surface= Flat` element
   may legally declare a small `KrElt` -- `docs/.../SegDemo.in` declares
   `KrElt= 0d0` -- so a Kr-only test mis-reads a flat detector as a
   curved one.  Measured: with Kr-only gating, `set_elt_kr(13,-100)` on
   e5hex1's Flat FocalPlane moves the FEX radius `-2548.0018521245 ->
   -2570.2492488727`; with the SrfType gate it does not move at all.
   Note the **Rx path cannot reach this** -- `msmacosio.inc:1338` AND
   `ChkDf2` (`iosub.inc:1105`) both force `KrElt=-1d22` whenever
   `Surface= Flat`, so a hand-edited .in is normalised before FEX sees
   it.  The API's `elt_kr` setter writes `KrElt` directly and leaves
   `SrfType` alone -- and that is exactly the path the design layer
   uses, so the gate earns its keep.
3. **Only the conic-BASE family is refined** (Conic/Aspheric/Monomial/
   Zernike/the grid composites).  Toric, Anamorphic, Interpolated and
   UserDefined give `Kr` other meanings; they keep the plane leg and say
   so.  A grid/FreeForm FIGURE on iElt+1 is ignored -- um-class against
   a sub-mm sag.  `Kr = 0` is a degenerate POINT, not a surface (ConSrf
   kills every ray on such an element), so it keeps the plane leg too.
4. **Flat-next decks are bit-identical BY CONSTRUCTION** (the plane form
   IS the `|Kr| -> inf` limit) and that is ASSERTED, not assumed:
   e5hex1 / SegDemo3conic / iris_dp / keckFF / eac2_7seg all give
   `|drad| = |dvpt| = |dpsi| = max|dW| = 0.000e+00` pre vs post.

**Measured on the zoom fixture** (`templates/50_sensitivities/zoom_5x5/
jwst_ote_designc`, model 128, ng 63, stop 25, OPD at 27, +-2.90888e-4
rad; elt 28 = curved focal surface `Kr=-3017.560611`, h = 54.2 mm,
sag 0.487 mm): corner nominals `4.26/4.37/3.09/3.23e-5` ->
`7.33/7.96/11.47/9.92e-6` mm (5.8x/5.5x/2.7x/3.3x), center unchanged to
1.4e-9 relative.  FEX == SXP to 0.000e+00 on every deck checked.

**Blast radius, from an ENGINE-truth scan** (`fs_fix/scan_engine.csv`;
396 decks loaded, 208 with a FEX-writable EP at nElt-1): **11 move, all
of them the JWST OTE focal sphere `Kr=-3017.56`** -- the zoom deck + its
grid sibling, `j18dcWithStop`, `j18mono`, `j18sc`, plus the A/B run's
generated copies.  197 are flat there and are bit-unchanged.  **Do NOT
text-parse .in files to reproduce this count**: several decks declare an
`nElt` that disagrees with their `Element=` block count (e5hex2 declares
24, the engine reports 25; eac2_7seg declares 47 with 48 blocks), which
shifts the indexing and attributes the WRONG element's Kr.  That error
is what put keckFF (`Kr=-0.0200`), iris_dp (`-2516.07`) and eac2_7seg
(`-306`) on the earlier affected list -- the engine reports `Kr=-1e22`
(flat) for element nElt on all three, and none of them moves.

**Still open, deliberately:** the FEX radius doubles as the far-field
PROPAGATION distance, where a tangent PLANE target is genuinely the
right thing.  This slice moved the OPD reference; if a physical-optics
leg ever needs the plane distance separately, the two uses have to be
split.  And the legacy `zp_iEm1` leg is NOT a general substitute -- it
equals the surface leg on the jwst deck only because elt 26 is
coincident with elt 28.

## Re-traces are IDEMPOTENT: OrthoSrcFrame + the aim dead bands (2026-09-08)
*(Documented: Manual §5 Background, last paragraph (bit-identical re-traces, the 16/2 ulp bands, the NOTE).  Below: the engine record.)*
Luis's 'centre-channel speckle' in dwdsurf_9zoom_5fov (Terminal Opus,
`REPORT_sens_noise_center.md`): a strict 2-CYCLE of the OPD across
identical traces on the jwst zoom deck -- max 2.2e-11 mm = 6 ulp of the
24 459 mm path on ~2200 of 2207 rays, W3 == W1 exactly -- which every
finite-difference column divides by 2*delta (1.19e-6 at delta 1e-6),
fixed per (field, zoom) because the pattern is deterministic.  CAUSE
(measured, not the stop -- it is there with NO stop set, and
ChfRayPos never moved): the source-frame re-orthogonalisation done at
every grid setup (`z = +/-ChfRayDir`, `y = unit(z x x)`, `x = y x z`;
ColSource, PtSource, ssrcray.inc, the STOP-OBJ cleanup at cmd_loop
label 205) has NO floating-point fixed point on a general frame: fed
the frame it produced, it returns a 1-ulp neighbour and the two
alternate (xGrid(1) -0.99999991481031303 <-> ...292 every trace).  A
full-precision orthonormal input alternates too; e5hex1's exact
(-1,0,0) is its own fixed point, which is why it was idempotent.  FIX:
`OrthoSrcFrame` (math_mod) computes the candidate frame and KEEPS the
incoming one when they differ by round-off only (every component within
1e-14); all four sites call it.  Same disease, second place: the
object-space STOP translation `ChfRayPos += x*xGrid + y*yGrid` (run at
every load with an `ApStop=` header) alternated eac2_7seg's SAVE
round-trip by 2 ulp of ChfRayPos(3) -- the 'pre-existing 2-ulp
ChfRayPos re-aim oscillation' the SAVE work met in July; dead band
1e-13 x max(1, |StopPos|, |ChfRayPos|) on |x|,|y|, and on the
point-source direction re-aim (1e-14 on the unit vector).  Measured
after: ten traces bit-identical on jwst (with the header ApStop, with
`stop 25`, and with no stop) and e5hex1; the elt-4 (virtual
CenterSegment) dw/dsurf column is EXACTLY 0 (was 1.19e-6 rms of
salt-and-pepper), elt 5's columns unchanged and smooth (roughness
0.01-0.02); a single load+SAVE is bit-identical pre/post.  First-trace
results never change (the band only suppresses the ulp re-application),
so single-trace gates are unaffected by construction.  `srcaim.inc` is a
commented-out include (dead) and still carries the old recipe.
**Bands tightened and made visible (2026-10-02, Dave: "how do we work
around the 'didn't change' consequence of the LSB dead band?").**  Each
band compares the FRESHLY computed target with the current state and skips
the update only when the residual is below it -- a residual, not an
increment, so a sub-band change is never accumulated and the error is
bounded by one band.  MEASURED on the corpus (instrumented scratch build,
`scratchpad/measure_dband.m`: jwst zoom with header / element stop,
e5hex1, Cass, dyson5 R4 and telescope decks, eac2_7seg; element, object
and header stops): round-off residuals <= 0.5 ulp in the frame and <= 0.7
ulp (of the position scale) in the aim; the smallest GENUINE update 7e5
ulp.  So `math_mod` now carries `DeadBandUlp = 16` (was 1e-14 = 45 ulp on
the frame and 1e-13 = 450 ulp on the positions) and `DeadBandQuietUlp = 2`:
a residual in (2, 16] ulp is suppressed AND reported ONCE per run by
`DeadBandNote` ("** Note: a <site> update of N ulp ... was NOT applied"),
every such event counted in `nDeadBandNote` (api `deadband_notes_get`,
mex command of the same name; zero on every deck measured).  A nonzero
count means round-off larger than measured or a caller stepping the
source direction / stop position at round-off level -- both worth
knowing, neither flooded.  Found while gating: `set_src_fov` wrote
`ChfRayDir`/`ChfRayPos` WITHOUT invalidating the cached trace, so the next
trace re-used the old chief ray's frame; it now calls `modified_rx`.  And a
hand-written deck's `ChfRayDir` is unitised once at load (1 ulp), after
which SAVE -> load -> SAVE is byte-identical: assert the SECOND round
trip.  Gate `tRetraceIdempotent` (mmacos, SUITE_FAST): ten traces
bit-identical (header stop, element stop, no stop, segment stop), the
round-trip fixed point, and a direction change of 4500 ulp applied / 4 ulp
kept-and-counted / 1 ulp kept-and-silent.  The 5-iteration smoke test
`tOffsetImager/test_s3_resolve_recovers` went red at 81d3308 because its
LM path, not its property, depended on the old Jacobian noise; re-pinned
at 7 iterations (ratio 0.18, converged from 7 up) with the mechanism in
its header.

## Element STOP preserves the source frame's HANDEDNESS (2026-09-08)
*(Documented: Manual §5 STOP, the segment / handedness paragraph.  Below: the mechanism and the measurements.)*
`UpdSrcGrid` (sourcsub.F, the only caller of `define_local_csys` on the
SOURCE frame; reached from ChiefRayAiming, i.e. every ELEMENT stop)
rebuilt (xGrid,yGrid,ChfRayDir) right-handed.  The corpus is full of
LEFT-handed deck frames (`xGrid= -1 0 0` with ChfRayDir +z: e5hex1,
6MST, every Telescope-emitted deck), so an element stop flipped xGrid to
+1 0 0 while an object-space stop (cmd_loop label 205 re-orthogonalises
the deck frame) did not.  Two measured consequences: (1) FEX's legacy
single probe `5d-6*xGrid` flipped sign (the 1.58 mm e5hex1 finding;
moot since the four-probe FEX); (2) **on a SEGMENTED source the mirrored
xGrid mirrors the ray grid while `EltToSegMap` does not move, so rays
hit segments they are not mapped to and are OBSCURED** -- e2e6m
`s3_imager_full` + the add_pupil pair: 732 of 985 rays obscured after
`macos.stop(1)` (Segment) and 962 after `stop(25)` (Reflector), 2 after
`stop_obj(0,0,0)`; the 253 survivors are the 5 segments on the x=0
column, the mirror-symmetric ones.  This is what made `add_pupil` on
segmented decks look like it "kills" rays, and what every supervisor
harvest with an explicit `stop_elt` on such a deck was silently
suffering; accepting Segment stops widened its reach to add_pupil's
default `stop_elt=1`.  Fix: UpdSrcGrid records the incoming triad's
handedness and negates xGrid after define_local_csys when it was
left-handed (y and z are unchanged, so this is exactly the OBJ path's
triad).  Right-handed decks: bit-identical.  Gate: after `stop elt 1`
on e5hex1 the SAVEd xGrid is still `-1 0 0`; s3 obscured count 2.
Engine-side reporting gap seen en route (CLOSED 2026-10-12, see the
per-ray status section): the end-of-trace obscuration
stamp writes `RayFailElt = nElt+1` (29 on a 28-element deck), not the
clipping element.

## FEX probe is FRAME-INDEPENDENT (2026-09-08, Dave)
*(Documented: Manual §5 FEX, the five-ray paragraph and the `EP crossing` line.  Below: the magnitudes and the five re-pinned tests.)*
`FEXProbeCross` (tracesub_mod, called by BOTH `FEX` and `SXP`) traces
FOUR differential chief rays, `+/-5d-6` about two orthonormal axes
perpendicular to the source chief ray (azimuth seeded by xGrid, yGrid
fallback), crosses each with the chief (`FindCrossPt`, point ON the
chief) and returns the MEAN crossing distance; `CrossPt = cr1pos +
zpLeg*cr1dir`.  Why: the legacy single probe `th = 5d-6*xGrid` gave a
crossing with a term LINEAR in the probe angle on off-axis decks, so the
answer depended on the SIGN of xGrid (e5hex1: 2548.0019 vs 2549.5813,
1.58 mm along the chief ray, between the deck's left-handed frame and
the right-handed one an element STOP rebuilds) and on its azimuth
(tangential vs sagittal pupil).  The +/- pair is a central difference
(sign term cancels EXACTLY); the two azimuths average to the medial
pupil (rotation-invariant at first order).  Symmetric decks: all four
crossings coincide -> unchanged to round-off.  Telecentric test now uses
the LARGEST chief/probe sine; `nGood=0` (all probes lost) takes the same
station fallback with its own message; `nGood<4` warns and averages the
survivors.  Each run prints `EP crossing = mean of N probes; spread S`
-- S is the pupil-astigmatism / probe-asymmetry measure.
**MAGNITUDE, measured (the medial pupil is NOT the legacy tangential
one on off-axis decks):** e5hex1 axis1 2548.79 / axis2 2498.68 mm ->
medial 2523.74 (legacy 2548.00; T/S split 50 mm, 1%); j18sc stop-elt-4
3037.064 -> 3037.968 (split 1.8 mm); jwst zoom fixture null radii
3017.5x -> 3018.3-3018.7 (+0.7..1.2 mm, 2.4e-4..4.1e-4).  Independently
confirmed with the LEGACY engine on 90-degree-rotated source frames
(e5hex1 sagittal 2498.68, j18sc 3000.29 vs tangential 2997.89), so the
split is pupil astigmatism, not a probe defect.  **tFocalSurface pins
3 tests to the legacy radii** (`test_null_radii_are_pinned` 1e-7,
`..._match_the_ab_report_...` 5e-3 abs vs REPORT_wnom_cli_ab's V4,
`..._fex_radius_follows_the_fit` 3017.5444) -- they failed by exactly the
shift above and were **RE-PINNED to the medial values 2026-09-08 (Dave's
ruling)**, each with the legacy value and the mechanism in a comment.
Same re-pin: `tPupilFindMethod` -- the zoom cross-config vertex
separation is pinned at 5.536e-4 mm (the sagittal half of the medial
crossing sees the FSM deflection; the tangential one did not) and
e5hex1's FEX-vs-cone-fit gap at 1.302 mm (was ~23 mm: the medial FEX
AGREES with pupil_find's cone station).  Rule for the next FEX-definition
change: these five are the pins that move; re-pin with the value AND the
mechanism, never a tolerance bump.
NOT changed:
XPS (per-ray crossing cloud, still one probe about xGrid -- its vertex
can now differ from FEX by the probe-sign term on asymmetric decks),
the STOP-ELT entrance-pupil crossing, FPP/PFP and the WINDOW/PLOCATE
beam-frame probes (those DEFINE a frame from xGrid by design).

## FEX EP-radius rework (2026-07-03) + SXP command (Set eXit Pupil)
*(Documented: Manual §5 FEX, the KrElt bullet and the guards paragraph (telecentric, footprint, lost probes, SXP, XPS); cmdref FEXit / SXP / XPS.  Below: the engine record and the SMACOS dispatch rule.)*
**FEX now defaults to the EP→next-element radius** (Dave's spec): the
EP Return radius is ALWAYS the chief-ray distance from the EP
(CrossPt) to the **iElt+1 plane — whatever iElt+1 is** (FP,
coronagraph mask, ZWFS…) — because that is the far-field propagation
distance for physical optics.  No element-type scan.  Sign per EP
EltID (Return reverses the beam → −cr1dir; Reference passes →
+cr1dir).  Legacy `zp_iEm1` = fallback (no iElt+1 / degenerate plane)
and the footprint-autoswitch alternative.  Guards, all noisy:
(1) **telecentric** — parallel probe chief rays make FindCrossPt
divide 0/0 (Cp2→0, NO guard there); FEX detects sin<1d-12, keeps the
element station, radius = station→iElt+1 plane, FLAT 1d22 fallback;
(2) **beam-footprint sanity** — a reference sphere smaller than the
beam footprint at the EP guarantees k2<0 "surface miss" for marginal
rays (the SegDemo3 failure); autoswitches to the other leg if usable;
A third guard, the **Rx-order flag** (warn whenever the element
before the EP return is a Return), was **REMOVED 2026-09-08** (Dave,
j18sc `fex 27`): the pattern it prescribed — Reference@FP, Return@EP,
Return@FP — exists nowhere.  The manual's FEX setup is two Returns
then a FocalPlane/Reference, the FIRST Return AT the focus; a
514-deck census finds 119 ending Return/Return/FocalPlane and 0 ending
Reference/Return/Return; and the legacy `zp_iEm1` leg assumes iEm1 IS
the focal Return.  "Fires on most legacy Rx" was the guard being
wrong, not the corpus.  The SegDemo3 defect it was meant to catch is a
POSITION error (focal Return not at the focus), visible as the two
printed legs disagreeing — do not re-add a TYPE check.
**Compatibility (fex_sweep 2026-07-03):** conforming double-pass Rx
have the pre-EP Return AT the focus, so both legs are equal by
construction → round-off-level no-op (e5hex1, 6MST, iris, j18*,
dmt6mono, CoroExample, HOEExample, lst3zern, CassWithExitPupil).
Divergent: manual `Cassegrain.in` (EP Reference right after SM)
legacy −1.2295 → new +6.7907 == exactly CassWithExitPupil's value —
the rework computes the right answer where legacy was degenerate;
`eac2_7seg` −289.7 → +52597.6.  GMI regression all-pass.  Sweep tool:
`MACOS_resources/mmacos/tools/fex_sweep/`.  Negative-L background
(why small spheres fail SILENTLY): ordinary train rejects L<0 (GO TO
98 → miss); near Reference/Obscuring/Return/LensArray
`ifLNsrf=.TRUE.` allows it and ConSrf picks the quadratic root by
|L²−mpr| PROXIMITY to the ray's current vertex distance
(surfsub.F:~148) — no flow-of-light sense.

SXP (retained, now largely redundant with FEX): FEX clone whose EP
radius is the chief-ray distance EP→FP (assumed iElt+1) — the
geometry FEX now uses by default.  Captures FP-Tz perturbations
through the EP radius and catches EP shifts when upstream optics are
perturbed (use case: sensitivity / dw/dx loops).
- `SUBROUTINE SXP` in `tracesub.F` — cloned from FEX; after
  `FindCrossPt` replaces zp with `t_FP = -((FP_vpt - CrossPt) · psi_FP)
  / (cr1dir · psi_FP)`, i.e. the chief-ray distance from CrossPt to
  the FP plane (chief ray reverses at EP, so post-EP direction is
  `-cr1dir`).  Falls back to `zp_iEm1` if iElt+1 is out of range.
- Dispatch added in `macos_cmd_loop.inc` alongside FEXIT (uses label
  8073 to avoid the existing 73 collision) — reachable from both
  interactive macos AND SMACOS, since SMACOS `#include`s
  `macos_cmd_loop.inc` from `smacos.F:210`.
- `LoadStack` entry in `smacosutil.F` — one IARG (iElt) + one CARG
  (YES/NO), same signature as FEXIT.
- pymacos wrapper `m.sxp(mode=1)` (in MACOS_resources/pymacos/dr-dev2).

Pre-existing dispatch puzzle worth knowing: `macos_ops.F:MACOS_OPS`
is NOT the SMACOS dispatcher — it's only called from the inner
optimization loop (`smacos_compute.inc`).  SMACOS proper dispatches
via `#include "macos_cmd_loop.inc"` from `smacos.F`.  When adding
new commands callable via SMACOS, add them to `macos_cmd_loop.inc`
and to `smacosutil.F`'s LoadStack (so the args get on the STACK
that IACCEPT_S reads in SMACOS mode), NOT to `macos_ops.F`.


## STOP on a Segment element + multi-value prompt gather (2026-09-08)
*(Documented: Manual §5 STOP -- one-line entry, Segment allowed, NS refused, the parent-vertex trap.  Below: the mechanisms and the measured traps.)*
Two CLI fixes from Dave's j18sc `fex 27` session (no `ApStop=` in that
deck, so a STOP had to be set by hand):
- **`stop obj 0 0 0` on ONE line crashed** (`forrtl severe (24): end-of-
  file during read, unit -5`).  `DACCEPT`/`RACCEPT`/`IACCEPT` (macosio.F)
  get ONE blank-delimited token per `READ_LOH`, then list-directed READ
  ABSN values from it -- EOF on a short token, and `ERR=` does not catch
  EOF.  Now they gather: `READ(acc,*,IOSTAT=ios)`; on `ios<0` pull the
  next token (rest of the line, then further lines / journal lines, with
  a `(N values expected -- enter the rest)` note when the line is used
  up) and retry; `ios>0` re-prompts as before.  The comma form `0,0,0`
  was always one token and still is.  Journal contract unchanged: READ_LOH
  journals per token, and the gather reads across lines on replay.
  smacosio.F's copies (SMACOS stack path) are NOT touched -- the stack
  packs complete values.
- **STOP/CENTER rejected `Segment` (EltID 11) with `Invalid element type`**
  -- and indexed `EltID(iStop)` BEFORE the range check.  Now: range check
  first; the veto keeps only `NSRefractorElt`/`NSReflectorElt` (named in
  the message); Segment is allowed.  The one Segment-specific mechanism:
  CTRACE visits a Segment element only for rays whose `RayToSegMap`
  contains `EltToSegMap(iElt)`, so the chief ray is mapped to the stop
  segment right before `CALL STOP` (AFTER label 210, so the Rx `ApStop=`
  entry path gets it too); `SetSourceRayGrid` rebuilds the maps at the
  next grid setup.  `stop_info_set` (macos_api_mod) dropped the same
  Segment veto, so `macos.stop(4)` / `m.stop(4)` work on segmented decks.
  `COORD` still carries the old 9/11/12 veto -- not asked, not touched.
  Gate: j18sc `stop elt 4 0,0` (CenterSegment) then `fex 27` reproduces
  `stop obj 0,0,0` to 1e-13 (aiming convergence).  MEASURED TRAP: `stop
  elt 5 0,0` (an off-axis segment) gives the IDENTICAL stop -- not a
  bug: in segment-class decks every segment's `VptElt` IS the parent
  vertex (j18sc: all 19 at 0,0,0; only `RptElt`/apertures differ), and
  the offset is relative to VptElt.  A segment CENTRE needs an
  RptElt-based offset.
  **SECOND MEASURED TRAP (FEX side FIXED the same day -- see "FEX probe
  is FRAME-INDEPENDENT"; the ELT-stop frame flip itself remains):** on
  e5hex1, `stop elt 1 0,0` then `fex 12` gives EP radius 2549.581280 vs
  2548.001852 from `stop obj 0,0,0` (6.2e-4; vertex moved 1.58 mm ALONG
  the chief ray; chief-ray direction and StopPos identical to 1e-10).
  Cause, pinned by SAVE-diff: the OBJ path keeps the deck's source frame
  (`xGrid = -1 0 0`, LEFT-handed with yGrid/ChfRayDir), while every
  ELEMENT stop rebuilds it right-handed via ChiefRayAiming -> UpdSrcGrid
  -> `define_local_csys` (`xGrid = +1 0 0`); FEX's differential probe is
  `th = 5d-6*xGrid`, so the probe sign flips and the crossing walks
  1.58 mm on this off-axis deck.  Carried by the FRAME, not the stop
  type: `stop elt 1` then `stop obj 0,0,0` (which preserves xGrid)
  reproduces 2549.581280 exactly.  j18sc (right-handed deck frame) shows
  no difference (1e-13).  Gate `tStopReload/
  test_stop_accepts_a_segment_element` asserts direction + on-line
  vertex and tolerates the radius to 2e-3.

