# MACOS engine notes — Polarization: conventions, coatings, s/p sign, radiometry, TrPolarizer/WavePlate, vector propagation

> **Read this file before touching:** `ifPol`, `RayE`, `Reflector`/`Refractor` coating branches, `PolElt`, `ifVecDif3`, `CPROPAGATE` legs.  These notes are the engine record
> (what moved, the measurements, the gates); they are NOT auto-injected — the index is
> `macos_f90/CLAUDE.md`, the root rules `../../CLAUDE.md`.  Text is moved verbatim from the
> former monolithic cheatsheet (2026-10-15 split); add new entries here, not to the index.

## Polarization physics (PLAN_POLARIZATION, Phase 0/1)
See `PLAN_POLARIZATION.md` (root) + `POLARIZATION_PHASE0_AUDIT.md` for the
full audit.  Quick map of what the engine has and the conventions to hold.

**What exists (all gated on `ifPol`):**
- **Per-ray vector E-field** `RayE(3,mRay)` (`Complex*16`, `elt_mod.F`) --
  the polarization state carrier, s/p-decomposed at each surface in
  `Reflector`/`Refractor` (`elemsub.F`).  Source state is `Ex0/Ey0`
  (`src_mod`), mapped onto each ray's local x̂/ŷ in `sourcsub.F`/`ssrcray.inc`.
- **Fresnel + recursive multilayer thin-film coatings** with complex index
  in `Reflector`/`Refractor` (`elemsub.F:432-547`).  `Reflector`'s TP/TS
  transmittance sub-blocks are `if(.false.)` dead code (mirrors return R only).
- **Vector diffraction** = the 3 Cartesian components Ex/Ey/Ez propagated as
  independent scalar fields, gated `ifVecDif3`.  Since **Phase 3a Tranche 1**
  this covers EVERY leg (see the Phase 3a section below); before it, only the
  far-field sphere→plane leg (PropType 3, via `PFFPROP`) was vectorized.
- CLI: `POLARIZATION`/`NOPOLARIZATION` (sets `ifPol`; enables `ifVecDif3`
  when `mWF≥3` -- all stock model sizes have `mWF=3`), `VECTOR`/`SCALAR`.
  `POLARIZATION` is SMACOS-dispatchable (LoadStack packs `Ex0/Ey0` as
  `DARG(1..4)`, `smacosutil.F:160`).

**Stubs / gaps (do not assume these work):**
- ~~`RfPolarizerElt(14)`/`TrPolarizerElt(15)` are **name-table-only**~~ --
  **PARTLY CLOSED 2026-07-27, see the Phase-3 elements section below.**
  `TrPolarizerElt(15)` and the NEW `WavePlateElt(18)` are implemented and
  gated; `JmatElt(2,2,mElt)` is no longer dead (PolElt fills it).
  `RfPolarizerElt(14)` is STILL a name-table-only stub, deliberately --
  a reflective wire grid carries grid reflection efficiency and the
  substrate's own s/p response beyond the (now settled) axis rule, and
  nothing needs it yet.
  No VVC; no Mueller math in the engine.
- `srtrace.F`'s `ifPol=.false.` is a **local `parameter` in the dead
  `SRTRACE_Test` driver only** (its caller is under `#if 0`).  The production
  single-ray paths (`SRTrace`/`CRTrace`, tracesub.F) thread `ifPol` as an
  argument -- so there is nothing to "lift" here.

**Conventions (pinned in Phase 0 -- assert in tests, do not relegislate):**
| Convention | Value |
|---|---|
| Time-harmonic | `exp(+iωt)`; spatial propagator `exp(−ikz)`.  Derived from `elemsub.F:387` (`C1=exp(−i·2π·L·N/λ)`, phase decreases as `L` grows) + the coating recursion `elemsub.F:512-516`, consistent with the 2026-07-25 IFO finding (field phase advances as OPL shortens) and the pymacos↔PROPER `opd_sign_flip=True`. |
| Absorbing index | `N = n − iκ`, κ>0 = loss (as stored in `IndRefArr`/`ExtincArr`, applied `DCMPLX(n,−κ)`). |
| Jones storage basis | Linear (x,y); circular via a unitary change of basis (Phase 2 decision). |
| Coating thickness | Rx `Coating=` layer thickness is **waves at parse-time `Wavelen`**, converted to physical (`·Wavelen/IndRef`) at load (`msmacosio.inc:2660`); the trace applies phase at the *current* λ so broadband sweeps are already correct.  `Coating=` must follow `IndRef=` (boundary media snapshot).  **`coat_set` takes PHYSICAL thickness** and sidesteps all of this. |

**Two coating subsystems** (do not conflate): Model A = `Coating`/`EltCoat`/
`IndRefArr`/`ExtincArr`/`EltCoatThk` (complex index, drives the polarization
path, `coat_set`/`coat_get` target this).  Model B = `nCoatElt`/`CoatIndxElt`/
`CoatThkElt` + `AirGap` (real index, non-sequential refractive path only).
Extend A, leave B (Dave 2026-07-25).

**Phase-1 API surface** (`macos_api_mod.F90`, `beam_set`/`beam_get` template):
`pol_set`/`pol_get`, `vecdif_set`, `coat_set`/`coat_get`, `rayfield_get`.
`Ex0Ey0=` Rx keyword now parses (4 reals `ExRe ExIm EyRe EyIm`) + round-trips
through SAVE (STATE only; on/off is API/CLI).  See `SAVE_KEYWORD_AUDIT.md`.
`pol_set`/`vecdif_set` call `modified_rx` (Phase-2 fix): the POLARIZATION
command changes trace-relevant state (`Ex0/Ey0` seed `RayE` at source-grid
setup) but does NOT reset the cached trace -- without the dirty, a pol-state
change + re-trace harvests the PREVIOUS state's `RayE` (verified stale).

**Coated-branch fixes (Phase 2, this box).**  Two latent bugs in the
`ifPol`+`nCoat/=0` recursion, both in `Reflector` AND `Refractor`
(`elemsub.F`), both invisible to intensity-level tests and fatal to Jones
work; found by the Fresnel-analytic gate (`tJonesPupil`), which now pins
them at 1e-12:
1. **Incident medium**: the recursion read `nb_arr(0)` = the parser's
   `IndRefArr(0,iElt) = IndRef(iElt-1)` -- for a coated mirror FOLLOWING
   another mirror that slot holds the previous mirror's conductor-idiom
   substrate (`Extinc=1e22`), i.e. light modeled as arriving from inside a
   perfect conductor.  Fixed: the coated branch now uses `na,kxa`
   (= `CurIndRef/CurExtinc`, the medium the ray actually travels in --
   what the uncoated branch always used).  The stored slot 0 is now unused
   by the trace; the parser convention is unchanged.
2. **Signed incident cosine**: `ccfb_arr(0)=DDOTC(ihat,Nhat)` is NEGATIVE
   when the normal faces the beam (psi convention), which turns every
   interface coefficient into its RECIPROCAL (1/r): |R|>1, s/p roles
   swapped, retardance sign flipped -- while |D| survives (the sneaky
   part).  Fixed with `DABS` (mirrors the uncoated branch's `ccfa`).
   Diagnostic signature if it ever returns: measured/analytic RS/RP ratio
   = (RP/RS)^2 exactly.

**CLOSED 2026-07-27 (macos cb29ea5 + 25c4386) -- reflected-p-hat vs
Fresnel-r_p sign conflict.**  The diagnosis is kept below because it is
the reference case for how a whole gate set can be structurally blind;
the FIX and its evidence are in the next block.  Odd-mirror polarized
results are trustworthy again as of cb29ea5.  Original text:
`Reflector` assembles `Eout = prhat*Epr + shat*Esr` (elemsub.F:592) with
`prhat = shat x rhat` (:431) -- p-hat follows the OUTGOING ray -- but
line :455 uses `RP = (Na*ccfb-Nb*ccfa)/(Nb*ccfa+Na*ccfb)`, which is
**-r_p** in that convention (at normal incidence it gives `RP = RS`, the
signature of the p-hat-held-fixed convention).  The standard `+r_p` sits
commented out one line above as `! dcr's original`; the flip arrived in
the 2022 bulk import `e1fa721` with no recorded rationale.  Net effect:
near normal incidence the transverse field is REFLECTED about the local
p-hat instead of negated.  Measured on `Rx_Cass_FarField`, x-polarized,
perfect conductors: after ONE mirror `Py/Px = 1.0163` (a 50/50 x/y
mixture, at <2 deg AOI, where physics allows `O(sin^2 beta)` ~ 1e-3);
with `+r_p` restored it is `2.07e-4`.  Two things make it near-invisible:
a reflection is an INVOLUTION, so a mirror PAIR cancels it exactly (the
two-mirror number is BIT-IDENTICAL either way, 7.0612e-07), and it is
UNITARY, so the unitarity gate cannot see it.  `Rx_Cass_FarField` has
exactly two mirrors, `Rx_VecChain` has none, and the Fresnel gate builds
its "analytic" RPa from the engine's own expression (`tJonesPupil.m:166`)
and compares a RATIO -- circular in this sign.  Second, fixture-free
tell: the effect is FLAT IN PUPIL RADIUS (median |Ey/Ex| = 1.014 / 1.016
/ 1.010 / 1.005 at rho = 32/64/96/128 px), and any real isotropic-surface
effect must vanish on axis and grow as rho^2.  NOT the point-source
launch frame -- that contributes `~NA^2/2` (measured 1.02e-5 at NA
0.0045).  The same `-r_p` form is in the COATED branch (innermost `RP`,
per-layer `RP1`) and `Refractor` (:1102-1103, :1226-1228) -- unaudited;
transmission (`TP`/`TS`) is a separate question since p-hat does not flip
there.  Assembling with `pihat` instead is NOT an admissible alternative
fix: `pihat` is perpendicular to the INCIDENT ray, so it emits a
longitudinal component along the outgoing ray.  Full analysis, suite
impact and the decision: `macos/REVIEW_POL_SP_SIGN_2026-07-27.md`.
Reproducer: `MACOS_resources/mmacos/tools/pol_sp_sign_probe/`.

**Reflected-p̂ / r_p sign fix (2026-07-27).**  `Reflector` assembles
`Eout` on `prhat = shat x rhat` (p̂ follows the OUTGOING ray) but the
2022 bulk import ("fixed by jzlou") negated the standard Born&Wolf r_p
in BOTH branches (uncoated `RP`; coated innermost `RP` + per-layer
`RP1` -- the Airy recursion propagates the flip exactly, |RP|
unchanged).  Effect: the transverse field was REFLECTED ABOUT the local
p̂ instead of negated -- an involution (cancels EXACTLY on mirror
pairs) and unitary, so every two-mirror/unitarity/ratio gate passed
while a single near-normal mirror turned x-pol into Py/Px ≈ 1 (physical
bound O(sin²β) ≈ 1e-3).  Post-fix the single-mirror cross-pol is 2.07e-4
in power AND grows as ρ² (slope-driven, as physics requires; pre-fix it
was FLAT in radius).  Two-mirror results bit-identical -- all published
2a/2b/polval two-mirror evidence survives.  The commented `! dcr's
original` line was correct all along; the same file's (dead)
transmittance code was already in the standard convention.  Diagnosis:
Opus, `REVIEW_POL_SP_SIGN_2026-07-27.md` + reproducer
`mmacos/tools/pol_sp_sign_probe/`.  Refractor's copy of the coated
recursion is deliberately UNTOUCHED (its transmitted output contains
only r·r products -- invariant under the consistent flip); Refractor
audit + odd-mirror ρ² gate = scoped follow-on.  Gate-blindness lesson:
an "analytic" reference transcribed from the engine's own expression is
circular in exactly the sign it should check -- write gates' analytics
from the textbook form.
**Tail landed 2026-07-27 (Opus):** (1) the odd-mirror gates
`tPolarization/test_odd_mirror_crosspol_{pec_analytic,rho2_law}` --
after ONE mirror the whole PEC single-reflection Jones is pinned to a
closed form written from Born&Wolf (`Ey/Ex = -sin2phi sin^2 a / den`,
`Ez/Ex = -sin2a cos phi / den`, `den = 1-2 sin^2 a cos^2 phi`), AOI and
azimuth taken from RAY DIRECTIONS so no pupil-grid mapping is assumed;
measured median 2.1e-15 / max 5.9e-14, radial slope 1.87, and the
pre-fix engine (rebuilt) fails 7 of 8 assertions (residual median
1.1e+02, slope 0.033).  (2) **`Refractor` normalized too** (25c4386,
innermost `RP` + per-layer `RP1`): there the flip was INTERNAL -- the
element emits TP/TS, and RP reaches the transmitted field only via the
`RP1*RP` products in the Airy denominators, invariant under a
CONSISTENT flip -- verified bit-identical on a coated singlet, with an
INCONSISTENT flip built and measured (-3.2% transmitted power) so the
invariance is not an untested path.  (3) Report section
`polval/50_sp_sign.md.in`.

**Coated-`Refractor` transmission radiometry -- FIXED 2026-07-28 (macos
`a5e4288`).**  Coated and uncoated `Refractor` transmission used
DIFFERENT amplitude normalizations: the uncoated branch multiplies by the
radiometric factor `sqrt(n2 cos02/(n1 cos01))` (the `S1` at elemsub.F
~:1146), so `|TP|**2` IS the POWER transmittance; the coated branch
composed plain Fresnel FIELD coefficients through the Airy recursion and
omitted it.  Measured with an index-matched single layer (optically a
bare interface): coated/uncoated amplitude **0.8164965809 = 1/sqrt(1.5)
exactly** at normal incidence, **0.7311104457** for BOTH s and p at 45
deg (`= 1/sqrt(1.5*cos28.13/cos45)` -- the equality across s and p is the
tell that a COMMON scalar was missing, since a Fresnel error would
separate them), and **0.6666666667 = 1/1.5 in INTENSITY at the detector
plane** -- propsub's `CPROPAGATE` chain inherited it, so the image
under-reported flux by a third.  A coated lens under-transmitted by ~18%
in amplitude.
**Decision (`REVIEW_POL_RADIOMETRIC_2026-07-28.md`, Fable lane):** keep
the incumbent power-amplitude convention; bring the coated branch to it
with ONE factor `sqrt(Re(n_sub)*cos_sub/(na*cos_inc))` applied ONCE
after the Airy recursion, to `TP` and `TS`, using `na` (the medium the
ray is actually in -- the Phase-2 incident-medium fix) and not the
stale `nb_arr(0)` slot.
**NEVER per interface.**  The conversion is boundary-to-boundary and the
interior-layer factors cancel identically.  Do not be talked into the
per-interface form by noticing that a plain chain's per-interface factors
telescope to the same number: that holds only while every interface is
present, in order, and paired with the RIGHT media -- and this branch has
already been bitten there once, by `nb_arr(0)`.  The site comment carries
the argument.
**What did NOT move**, and is gated: `RP`/`RS` and the `RP1*RP` Airy
denominators are untouched, so nothing REFLECTED changes; and the factor
is a common real scalar, so it cancels in `t_p/t_s` -- the transmitted
polarization STATE (D, retardance, every Jones-pupil quantity) is
unchanged, which is what lets the Phase 2a/2b coated results stand.
Only the absolute transmitted amplitude moved.
**Scope, deliberately ungated:** transmission into an ABSORBING (metal)
substrate has no meaningful power convention -- the transmitted wave is
evanescent and carries no propagating power.  Coated LENSES (dielectric
substrate) are the use case.  `Reflector`'s transmittance blocks stay
`if(.false.)` dead code -- untouched.
**Gates:** `tPolRadiometric` (mmacos, SUITE_FAST, 13 tests) on two new
fixtures, `Rx_Refract.in` (normal-incidence air/glass/air plate with a
physical-optics leg, so BOTH dispatch chains are covered) and
`Rx_Refract45.in` (45 deg AOI -- load-bearing: at normal incidence both
cosines are 1, so a normal-only suite gates the INDEX half of the factor
and is blind to the cosine half).  Every analytic is the Abeles
characteristic matrix written from Macleod ch. 2, NOT transcribed from
`elemsub.F` -- the r_p-sign lesson.  Non-vacuity: 7 of the 13 fail
against the rebuilt pre-fix engine.  Report section `polval/80_radiometric`.
**Analytic gotcha that cost this slice a cycle:** Macleod's
`2*eta0/(eta0*B+C)` is the TANGENTIAL amplitude coefficient; for p the
tangential component is `E*cos(theta)`, so it exceeds the ordinary
Fresnel `t_p` by `cos_sub/cos_inc` (1.2472 at 45 deg into n=1.5) --
exactly the size of a plausible radiometric error, so it reads as a
failed gate rather than a units slip.  `T` is unaffected either way,
which is why every transmittance gate passed while the amplitude-ratio
one did not.
**Still open (predates this landing, deliberately uncoupled):** the
`cos` term in the radiometric factor overlaps conceptually with the
beam-area bookkeeping ray DENSITY also carries when the grid is seeded
for diffraction.  Whether that double-counts in oblique refractive
seeding is an open audit question -- the uncoated branch has always
carried the factor, and the PROPER comparisons are mirror trains that
never exercised it hard.

**Coated REFLECTOR anchored to publication (2026-07-28) -- no engine
change.**  The coated `Reflector` recursion is now checked against an
OUTSIDE source on a DIELECTRIC-ON-METAL stack (a protected mirror), the
one configuration class our own gates never covered: `tJonesPupil`'s
Fresnel gate is an optically THICK SINGLE layer (algebraically a bare
interface) and `tPolRadiometric` is TRANSMISSION.  Method: reproduce the
publication's own configuration via `coat_set` -- their indices, their
4.12 nm Al2O3 on 220 nm Al, their 500/550/600/650 nm, their 6-70 deg --
and compare curve-on-curve with their model (van Harten, Snik & Keller,
PASP 121, 377 (2009), arXiv:0903.2740v1).  **Result: 2.828e-14
(diattenuation) / 4.937e-14 (retardance) against their stated +-0.01 per
normalized Mueller element.**  Driving it with THEIR indices is the
point: aluminium index tables genuinely disagree (the paper says so and
FITS k), and an index-table difference is NOT a machinery error.
Their index convention `N = n - ik`, k>=0, is IDENTICAL to ours, so
nothing translates -- and it is load-bearing, since they report the
opposite sign made their fitted oxide come out ~50 nm instead of ~4 nm.
Tool `mmacos/tools/pol_external_anchor/`, gate `tPolExternal`, polval
section 8, packet `macos/REVIEW_POL_EXTERNAL_2026-07-28.md`.
**Two measurement traps this exposed, both of which are self-cancelling
at exactly 45 deg -- the one angle the pre-existing Fresnel gate runs
at, so the suite could not have caught either:**
(1) a fold mirror's DEVIATION is `180 - 2*AOI`, NOT `2*AOI`; building a
rig the other way sweeps the COMPLEMENT of the intended angles.
(2) the single-trace s/p extraction divides out the input state's
projection on the s/p axes and therefore needs the engine's per-ray
LAUNCH FRAME, in practice a hard-coded `xGrid` -- wrong elsewhere, and
the symptom is diattenuation FLAT in AOI where an isotropic surface must
give `D ~ theta^2` (the same "flat where physics demands a power law"
tell that exposed the r_p sign defect).  Use the frame-free TWO-trace
form instead: `M = diag(r_s,r_p) R(phi)` so `r_s/r_p = M11/M22 =
-M12/M21`, with `phi` cancelling and the two estimates cross-checking.
**p-hat bridge, MEASURED:** on the perfect-conductor idiom the engine
returns `r_s/r_p = +1` exactly (imag part 0) at every AOI, so the ratio
measured in the `(s, p_r)` basis is already in the publication's
fixed-transverse frame -- bridge ZERO.  Assuming the `pi` the
ray-following doctrine suggests is wrong by exactly 180 deg everywhere;
`tPolExternal/test_phat_convention_bridge_is_zero` pins it.

**Overcoat quarter-wave reversal, engine-measured on BOTH sides
(2026-07-28) -- no engine change, no fixture change.**  The external
anchor above also found that the 110 nm MgF2 overcoat the Phase-2c ladder
applies is **0.607 quarter-waves at `Rx_Cass_FarField`'s own 1 um**, not
the 0.96 its `% ~quarter wave at 632.8 nm` comment describes -- and that
the overcoat polarization trade **REVERSES** across the quarter-wave
condition, inverting the §5.5 design rule.  That correction rested on an
independent analytic; it now has engine numbers on both sides.
**Ruling that shaped the fix: gated fixtures do not move.**
`Rx_Cass_FarField` stays at `Wavelen=1.0E-06` and `tPolContrast`'s
27.898/151.31 stay asserted; the companion wavelength is applied at
RUNTIME with `macos.set_src_wvl` AFTER `load_rx`.  Nothing on disk moves.
**Why that is a real measurement and not a unit conversion:**
`macos.coating`/`coat_set` take PHYSICAL thickness and the trace divides
by the CURRENT `Wavelen` when applying the layer phase, so a film is
fixed glass under a wavelength change and its optical thickness moves --
the same treatment `Retardance=` gets.  Measured (model 256, x-pol, both
mirrors coated, cross-polarized TOTAL power): MgF2/bare = **5.2707x at
1 um (costs)**, **0.2035x at 632.8 nm (suppresses)**, reversal
**25.899x**; a TRUE quarter wave gives **0.0532x** at either wavelength.
The 1 um leg reproduces the recorded 27.8977/151.311 exactly, in the same
harness.
**Report TWO ratios, they answer different questions:** the TOTAL power
ratio (what a designer sees -- carries the irreducible GEOMETRIC cross
term) and the coating-EXCESS ratio over the uncoated baseline (== the
ratio of the two `d_cross_rel`; the coating-only quantity, and the only
one comparable with a pure-Fresnel analytic).  They diverge most at the
true quarter wave, where the coating term is nearly extinguished and the
total sits at 1.5376x the UNCOATED floor -- which is why the engine's
total (0.0532) cannot follow the analytic (0.0157) down and is not a
disagreement.
**NON-VACUITY, and the pattern worth reusing:** the counterfactual is an
engine that pinned thickness in WAVES (achromatic) -- reachable from the
REAL engine by asking for 110 nm x (632.8/1000) = 69.6 nm at 632.8 nm,
the film with the same optical thickness in waves.  It lands back on the
1 um answer to 2.1e-8 and shows NO reversal; 3 of 3 of the gate's
632.8 nm assertions fail against it.  Three further identities localise
the effect IN THE FILM: the metal-only leg is lambda-invariant (1.3e-8),
the quarter-wave condition is lambda-invariant (1.9e-8, and 181.2 nm at
1 um vs 114.6 nm at 632.8 nm is DIFFERENT GLASS giving the same answer),
and the uncoated geometric floor is lambda-invariant (1.6e-15).
**Measure TOTAL POWER, never a fixed-pixel annulus mean, when comparing
across wavelengths** -- the annulus subtends a different lambda/D range
at each lambda and would mix the coating effect with the diffraction
scale.  Tool `mmacos/tools/pol_overcoat_chromatic/` (`oc_ladder` +
`oc_nonvacuity`), gate `tPolContrast/
test_overcoat_trade_reverses_across_the_quarter_wave_condition`, polval
section 8.3.1, packet `macos/REVIEW_POL_OVERCOAT_CHROMATIC_2026-07-28.md`.
**Design rule, chromatic form:** a protective overcoat is neither
inherently a polarization cost nor a benefit -- the sign is its optical
thickness AT THE WORKING WAVELENGTH.  Specify overcoats at lambda/4 of
the working wavelength.

**Jones input basis (engine launch frames, `ssrcray.inc`).**  Collimated
sources launch every ray with `E = S*(Ex0*xGrid + Ey0*yGrid)` -- the
source-frame pair, UNIFORM over the grid.  Point sources use a PER-RAY
frame: `yray = unit(RayDir x xGrid)`, `xray = yray x RayDir` (reduces to
the global frame on the chief).  `S` = flux normalization (~1/sqrt(nRays))
-- a common real scalar carried by the Jones pupil.  ColSource re-
orthogonalizes the Rx frame as `z=+-Chf; y=unit(z x x); x=y x z`.  The
perfect-conductor mirror idiom (`IndRef=1, Extinc=1e22`) gives RS=-1,
RP=+1 +O(1e-22) in the engine's ray-following basis (prhat = shat x
rhat): polarization-neutral (E_out = -E_in at normal incidence), which
makes any stock conductor Rx a unitarity gate for free.  **RP=RS=-1 at
normal incidence is the signature of the 2022-import sign DEFECT**
(elemsub.F "fixed by jzlou" negated the standard r_p; fixed 2026-07-27,
see the s/p sign section below) -- neutral only in a fixed-transverse
basis, which is NOT the basis the assembly uses.

**Phase-2 binding layer** (mmacos `+macos/jones_pupil.m` + `pol_maps.m`,
pymacos `macos.py` same names): two-trace Jones pupil (double-pole default
basis / local-sp / global) + closed-form 2x2 Pauli polar decomposition
(D/Dvec, ret/retvec in [0,pi] w/ near-pi ambiguity FLAG, T, phase; pupil
mean and variation reported SEPARATELY -- the mean absorbs the system's
geometric rotation and is a state change, not an aberration).  D/T are
singular-value invariants (basis-independent); ret is exactly what
double-pole makes artifact-free (local-sp inflates ret var ~10-250x --
asserted in tests as the documented artifact).  Gates: `tJonesPupil`
(mmacos, incl. the Bench fold Fresnel gate) + `test_jones_pupil.py`
(pymacos, ifx-linked = the standing ifx smoke).

**Phase 2b low-order expansion** (`pol_zernike`, 2026-07-26, binding-side
only -- no engine change).  Expands the Dvec/retvec maps onto ANSI
Zernikes so results are comparable with the polarization-aberration
literature term-by-term.  The engine-relevant fact: an on-axis
rotationally symmetric two-mirror system MUST reduce to pure
**polarization astigmatism** -- astig0 in Pauli s1, astig45 in s2, equal
magnitude, no circular component, no defocus, and a rho^2 radial law
whose on-axis extrapolation vanishes (measured: forbidden terms at 1e-15
of the astigmatic term, D(rho=0)/D(rho=1) ~ 1e-4).  This is now the
sharpest cheap check that the coating/s-p machinery and the pupil
reference frame are both right: a broken frame breaks the mode pattern,
not just the map.  Residual asymmetries are DISCRETIZATION, verified by
scaling (astig-pair mismatch 1.9e-7 at model 128 -> 5.8e-8 at 256; the
symmetry-breaking magnitude term is quadrafoil-X, aligned with the pixel
grid's own axes, while quadrafoil-Y stays at 1e-17).

## Phase 3 polarizing elements -- TrPolarizer + WavePlate (2026-07-27)

`PolElt` in `elemsub.F` serves `TrPolarizerElt(15)` (ideal linear polarizer,
finished from a name-table-only stub) and the NEW `WavePlateElt(18)` (linear
retarder; `mEltTypes` 17 -> 18).  Rx keywords `PolAxis=` (3-vector) and
`Retardance=` (waves at parse-time Wavelen), `ChkDf2` requires both on the
types that use them, SAVE writer inverts the Wavelen scaling.  API
`polelt_set`/`polelt_get`/`jmat_elt_get`; mmacos `macos.polarizer` /
`macos.waveplate` / `macos.elt_jones`.  Gates: `tPolElement` (27, SUITE_FAST).

**Geometry is RefSrf's, verbatim** -- conic intersection, `rout=ihat`, the
same `C1=exp(-i*2pi*L*N/lambda)` propagation phase, the same ChkRayTrans
block.  So with `ifPol` off the elements ARE Reference surfaces; that is
gated against a twin fixture (`Rx_PolElt_Ref.in`) rather than argued from
the source.

**`JmatElt` is now live and its per-ELEMENT shape is exact, not a
shortcut.**  The element Jones in its OWN eigenbasis -- `diag(1,0)` for a
polarizer, `diag(1,exp(-i*delta))` for a retarder -- carries no ray
dependence; ALL of it lives in the basis (`ahat` = transmitted / fast axis,
`bhat` = blocked / slow axis, built per ray by convention 2 below).

**Conventions -- all five EXTEND the pinned set, none is new law:**
1. *Axis as a 3-VECTOR*, not an angle in some element frame, so no
   "which in-plane direction is zero degrees" convention has to be invented.
2. *Project the MATERIAL axis* into the ray's transverse plane -- see the
   off-normal section below for why this is a physical statement and not
   bookkeeping.  For a WavePlate the material axis is the declared (fast)
   axis; for a TrPolarizer it is the ABSORBING direction, the in-element
   complement `psi x PolAxis` of the declared pass axis, which is projected
   and extinguished while its partner is transmitted.
3. *Orthonormalize the partner axis* rather than projecting a second
   declared axis.  FORCED, not chosen: a lossless retarder must be unitary,
   and only an orthonormal eigenbasis makes `diag(1,exp(-i*delta))` unitary.
4. *Retardance sign read off the engine.*  `C1=exp(-i*2pi*L*N/lambda)` means
   the slow axis takes the more negative phase, so fast-axis-declared gives
   `diag(1,exp(-i*delta))`.  Pinned by the SIGNED circular Stokes parameter
   (`S3/S0 = -1` for linear-in / QWP at 45 deg), which flips if the
   convention flips -- a `|S3|` gate would accept either.
5. *Retardance is stored PHYSICALLY* ((n_s-n_f)*d), Rx value scaled by
   Wavelen at load, divided by the current lambda at trace -- the same
   treatment `Coating=` thickness gets, so a plate is fixed glass and sweeps
   are chromatic.

**THE TRAP THAT COST THIS SLICE A REBUILD: there are TWO element dispatch
chains.**  `tracesub.F` traces rays for ray-level queries; `propsub.F`
(`CPROPAGATE`) RE-TRACES the rays that seed the diffraction grid, through
its own `ELSE IF (EltID(iElt).EQ.n)` chain.  Wiring only `tracesub.F` gives
an element that works perfectly in `ray_field` and is INVISIBLE to
`intensity`/`complex_field`.  Measured before the fix: crossed polarizers
took the ray power to 3.6e-33 while the detector plane sat at the full
9.69e-01 -- i.e. the grid reported that polarization had no effect, which
reads as "no polarization aberration" rather than as a bug.  Any future
element that touches the field needs BOTH chains.  (`srtrace.F` has a third
chain but it is inside `SRTRACE_Test`, under `#if 0` -- dead, nothing to
add.)  `tPolElement/test_grid_carries_the_polarizing_train` is the tripwire.

**SETTLED 2026-07-27 -- the off-normal axis convention: PROJECT THE
MATERIAL AXIS.**  Away from normal incidence the declared axis has to be
brought into the plane transverse to the ray, and orthographic projection
does not preserve orthogonality, so projecting the declared PASS axis is
NOT equivalent to projecting the BLOCK axis and transmitting the
complement.  The invariant object is the direction fixed in the SUBSTANCE
of the element -- the absorbing direction (a wire grid's wires, a dichroic
sheet's dipole chains) for a polarizer, the crystal fast axis for a
waveplate -- so that is what gets projected.  Not a taste call: the same
ambiguity was adjudicated experimentally for a tilted dichroic polarizer by
**Korger et al., Opt. Express 21, 27032 (2013)**, in the dipole model's
favour.

The waveplate branch was already the settled rule (its declared axis IS the
material axis); the polarizer branch was flipped, keyword semantics
unchanged -- `PolAxis=` still declares the TRANSMISSION axis, and the
absorbing direction `psi x PolAxis` is derived from it.  Size of the
change: **3.56 deg of transmitted-axis orientation at 20 deg AOI**, closed
form `acos(2 cos a/(1+cos^2 a))` at 45 deg azimuth, bounded by `sin^2 a`,
and **identically zero at normal incidence AND when the axis lies in or
perpendicular to the plane of incidence** -- so the obvious test tilt (axis
along x, tilt in x-z) is degenerate and reports a reassuring zero that
means nothing.  Two consequences worth knowing:

* **Normal incidence is bit-identical across the flip**, by construction
  and verified: the partner axis is completed by a cross product with the
  ray direction, which returns the declared axis to the last bit up to an
  overall sign a projector cannot see.  `PolElt` deliberately does NOT
  `DUNITIZE` the polarizer's transmitted axis -- normalizing an already-unit
  vector would perturb it by an ULP and break exactly that.
* **A polarizer's declared axis is now taken modulo its component along the
  element normal** (only `psi x PolAxis` enters), and a pass axis PARALLEL
  to the normal is a prescription error that extinguishes, alongside the
  pre-existing material-axis-along-the-ray case.

Gated engine-side at 20 deg AOI / 45 deg azimuth on BOTH dispatch chains by
`tPolElement` section F, on the tilted-BEAM fixture `Rx_PolElt_Tilt.in`
(tilting the ELEMENT does nothing to a collimated on-axis bundle).
`RfPolarizerElt(14)` is still a stub -- the convention no longer blocks it,
but a reflective wire grid carries grid reflection efficiency and the
substrate's own s/p response beyond the axis rule, and nothing needs it.
Packet: `macos/REVIEW_POL_ELEMENTS_2026-07-27.md`.

**Idealizations (documented, not worked around):** no ray splitting (a PBS
is two traces, or better a coated `Reflector` at 45 deg where s/p is
physical), no walk-off, no face Fresnel loss, no substrate; the output is
purely transverse, so any longitudinal component at the surface is dropped
(exactly zero for a collimated normal-incidence beam, O(NA) otherwise).
The waveplate is an IDEAL retarder: its retardance does not vary with
incidence angle, where a real crystal plate's does (the field-of-view
effect behind compound and Pancharatnam designs).  Anyone needing that
needs a birefringent-plate model with o/e indices and thickness.

**Tranche-1 interaction to respect when building fixtures:** a polarizing
element placed AFTER the first physical-optics leg transforms rays but never
reaches the grid (the seed already happened).  `Rx_PolElt.in` puts all four
polarizing elements BEFORE its single `NFPlane`->`Geometric` leg for exactly
that reason, and says so in its header.

## Phase 3a Tranche 1 -- vector propagation across the whole chain

All in `propsub.F` (`CPROPAGATE` + two new module helpers).  Vector mode
now closes a multi-leg chain; the coronagraph case (pupil→FPM→Lyot→focal)
works, the VVC acceptance test is unblocked.

**What changed.**
- **Every leg loops the SAME scalar kernel over the 3 component planes.**
  One `kWF1..kWF2` range is computed once before the leg dispatch
  (`1..3` under `ifVecDif3`, else `iWF..iWF`), so every `DO kWF=kWF1,kWF2`
  degenerates to the original single call in scalar mode -- those paths
  stay **bit-identical**.  Covered: `NFPROP` (PropType 2, 5), `PPPROP`
  (4, 6), `SFPROP` (7, 8), `SPH2PL`/`PL2SPH` (10, 11 -- NOT in the plan's
  list, added for coverage completeness), `FRPROP` (12), `NFPropDFT`
  (13, 14), `FFPropDFT` (15).  Leg dx/z bookkeeping stays OUTSIDE the loop;
  the `DWF` scratch is reusable across planes.
- **`PFFPROP` retired** (3a.1(1)).  It was a bare per-component FFT that
  omitted the Fresnel output factors `FFPROP` applies via `applyfac2`
  (`1/(i·λ·dz)·dx1²` + the output quadratic phase).  The far-field leg now
  runs `FFPROP` per plane.  The routine survives in the file, uncalled,
  with a deprecation banner.  **A/B (Rx_Cass_FarField, model 128):** vector
  total power was `8.9377e-01`, is now `1.8155e+06` == the scalar total
  exactly (Parseval).  Factor 2.03e6 in intensity, 1.43e3 in amplitude.
- **Assembly: seed once, then update** (3a.1(2)).  Both `ifPol` branches
  used to RELOAD the grid from `RayE` at every physical-leg assembly,
  erasing earlier legs' diffraction (the non-pol branch always MULTIPLIED
  by the incremental geometric phase).  A SAVE'd `LWFSeeded`, reset with
  `CumLStart` at trace start (`iStartElt==0`), now seeds on the first
  assembly and multiplies thereafter, advancing `CumLStart` identically in
  both branches.
- **Vignetting at the seed** (found while implementing, not in the plan).
  `RayE` carries NO aperture clipping -- the surface routines report it via
  `LRayTrans` and never zero `Evec` -- so a seed straight from `RayE`
  resurrected rays the ray-side zeroing had already extinguished.  Both
  polarized branches now gate the seed on `LRayPass(iRay).OR.(iRay.EQ.1)`
  (the `.OR.` mirrors the `iRay.GT.1` guard the zero sites carry for the
  chief ray).  Consequence worth knowing: **polarization-ON / vector-OFF is
  now BIT-IDENTICAL to polarization-OFF** -- it was wrong by 21% after one
  leg and 38% after two.
- **Masks are 3-plane.**  A 0/1 transmittance is a diagonal Jones `t·I`.
  `FFObscure` is called per plane under `ifVecDif3`, each call handed a
  FRESH copy of `xObs` (it re-orthogonalizes its `xGrid` argument IN PLACE,
  so reusing the mutated vector could move an edge pixel by a ULP between
  components); the scalar call is untouched.  The 13 ray-side clip sites
  and the 2 taper sites go through new module helpers `WFZeroPt` /
  `WFScalePt`, which are the original single-plane statement in scalar
  mode.

**PHASE CONVENTION AT THE SEED -- read before "fixing" it.**  The seed is
`WFElt(i,j,k)=RayE(k,iRay)` with **no** phase correction, and that is a
MEASURED result.  Reading `elemsub.F:395` alone (`S1=-TWOPI*L/lambda`,
`C1=exp(i*S1*Na)` applied to `Evec` at every surface) suggests `RayE`
advances phase OPPOSITE to `WFElt`, whose non-pol assembly uses
`exp(+i*TPL*Cphi)` and whose kernels carry `exp(-i*pi*lambda*dz*f**2)`
(paraxial forward `exp(+ikz)`).  Two independent experiments say `RayE` is
already in `WFElt`'s convention: (a) seeding `RayE*exp(i*c*TPL*CumRayL)`
and scanning `c` over `{-2,-1,0,1,2}` on a tilted Cassegrain moves the
far-field centroid by exactly `(c+1)x` the scalar tilt shift, so `c=0`
reproduces the scalar pupil phase; (b) on `Rx_VecChain.in` the `c=0` seed
reproduces the scalar intensity to 4.5e-16 at every leg for x-polarized,
45° and circular input.  **RESOLVED by the Fable review (2026-07-26):**
no bridge at the seed is CORRECT, but the reason is structural, not "RayE
happens to share WFElt's convention".  The Return-leg sign flip
(`IF (ifReturnElt(iElt)) RayL(iRay)=-RayL(iRay)`, applied AFTER the
surface call) makes the two phase bookkeepings diverge by construction:
`CumRayL` SUBTRACTS return legs (physical OPD bookkeeping, equalized at a
proper reference sphere) while `RayE`'s per-surface
`C1=exp(-i*2pi*L*N/lambda)` phases ADD them (`L` pre-negation).  So the
RayE-vs-CumRayL phase relation is **TRAIN-DEPENDENT**: on the
Return-terminated Cass FF the assembled vector EP field matches the
scalar field's convention (measured: circular-concentration same=0.994 /
conjugate=0.002; signed far-field tilt response equal at +0.94x), while
plain trace-to-detector `RayE` measures CONJUGATE to the OPD map (slope
-0.9995, corr -0.9995).  Confirmed on the GRID 2026-09-12 (zwfs_dm96
Twyman-Green test arm: retro DM, two lenses, PEC plate BS, NF sandwich;
model 512): the vector-mode field at the detector has the scalar trace's
amplitude structure but its common phase fits **-2.0000 x the scalar
pupil phase, residual 4e-9** -- an exact conjugate, plus the Fresnel
losses the scalar trace never applies (T 0.66).  `dmg_arm_maps` (the
vector-ZWFS arm maps) does not bridge it: it takes J/sqrt(det J), which
drops the common scalar and keeps only the polarization part.  Anyone
comparing a vector-mode field to a scalar one phase-for-phase on such a
train must expect the conjugate.  There is NO universal convention bridge to
apply -- correctness rests on the behavioral gates (exact on the
mask-type/normal-incidence class Tranche 1 claims), and **Tranche 2's
J_run must carry phases explicitly relative to the CumRayL bookkeeping**
rather than assuming either sign.  The tests are the authority; re-run
`tVecChain` before changing anything here.  Note the scalar-pol branch
cannot serve as the vector seed regardless: it rebuilds from `|RayE|` and
throws away the per-component phases.
Debug gotcha that cost this review an hour: obscured rays carry
`RayE=0`, and `ATAN2(0,0)=0` reads as "zero phase, zero Ez" -- gate ANY
per-ray probe on `LRayPass` before interpreting phases.

**Gate prescription.**  `mmacos/tests/Rx/Rx_VecChain.in` (copied to
`pymacos/tests/Rx/`) -- collimated on-axis source, flat normal-incidence
uncoated planes, TWO bracketed near-field legs (`NFPlane`→`Geometric`, the
`Rx_Coro.in` idiom; consecutive same-PropType elements MERGE into one leg,
which silently defeats a multi-leg test), aperture + central obscuration.
The ray E-field direction is a CONSTANT unit vector, so `E_k = e_k·u(x,y)`
and the vector sum `Σ|E_k|²` must equal the scalar intensity EXACTLY.  On a
real off-normal train it cannot: `Rx_Cass_FarField` carries `|Ez|/|Ex| ≈
8.8e-2` at the exit pupil and vector/scalar legitimately differ ~2.6e-3.
**Run the gate at 45° and circular input, not just x-pol** -- with an x-only
source all the energy sits in plane 1, which the OLD single-plane
propagator carried correctly, so an x-pol-only gate passes VACUOUSLY.

**Constraint (document, don't work around):** vector mode repurposes the
`mWF=3` planes as Ex/Ey/Ez of ONE wavefront -- no multi-WF / COMPOSE
concurrently.

**`cfield_plane_get` (2026-07-26) -- reaching the component planes.**
`cfield_get` only ever returned `iEltToiWF(iElt)`, so in vector mode the
per-component field was unreachable from the bindings: callers could see
the summed intensity and never how Ex/Ey/Ez made it.  New sibling
`cfield_plane_get(OK, RE, IM, N, iElt, iPlane)`: `iPlane=0` is
`cfield_get` exactly (bit-identical, and that is what the bindings pass by
default); `iPlane=1..3` are Ex/Ey/Ez and are REFUSED unless `ifVecDif3` --
in scalar mode plane k is an unrelated wavefront, not a component, and
handing it back would look plausible and be wrong.  Surfaced as
`macos.complex_field(srf,'plane',k)` (mex cmd gained an optional 4th arg)
and `pymacos.complex_field(srf, plane=k)`.

This is what CLOSED the Tranche-1 attribution, and it corrected it: the
vector/scalar difference on an off-normal train is NOT simply the
out-of-plane intensity (that guess is off by ~2x).  Two mechanisms --
(1) the scalar seed `|RayE|` puts ALL the power, including the
out-of-plane part, into ONE propagating plane while the vector run leaves
only f=0.9979 in Ex (a near-pure rescale, 1-corr = 4e-8), and (2) Ey/Ez
diffract to their own pattern.  `Iv ~ f*Is + Iy + Iz` takes 2.56e-3 down
to 2.90e-4.  Pinned by `tVecChain/test_vector_scalar_difference_decomposition`.

**Still open -- Tranche 2.** Between physical legs the grid field is
advanced by a scalar phase, exact only when the intervening elements are
non-polarizing (Obscuring / Reference / FocalPlane).  A COATED or
reflecting surface BETWEEN legs (IFO recomb→detector with folds under
`DO_NEARFIELD`; VVC layouts with an OAP between masks) needs the per-ray
running Jones `J_run(3,3,mRay)` design in PLAN_POLARIZATION §3a.3.

**Gates:** `tVecChain` (mmacos, SUITE_FAST) + `test_vec_chain.py`
(pymacos, ifx-linked).  Both are non-vacuous -- checked 2026-07-26 against
the pre-fix engine, which fails them at 0.21..0.38 relative error and
mis-states total power by 4-7%.

