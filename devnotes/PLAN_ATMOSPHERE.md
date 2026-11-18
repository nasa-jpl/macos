# PLAN_ATMOSPHERE — atmospheric turbulence in MACOS: the single screen, the layered chain, and the regrid at the entrance pupil

**Status: DESIGN (Dave + CC, 2026-10-18).  Build is a future activity, after the ATMOS fix of
PLAN_CONSOLIDATION 1.2 (p) lands.**  Nothing here is built; the engine facts cited are measured.

## 1. What exists

- `ATMOS` (CLI; no api entry yet) puts ONE Kolmogorov phase screen on the ray state at an element:
  `AtmosPhase` (sourcsub.F) generates the screen (the generator and the RNG are right —
  `cli_tests/records/atmos_kolmogorov_check_2026-10-18.md`), `AddAtmos` adds the screen's OPD to each
  ray's `CumRayL` and tilts each ray by the screen's gradient.  Its coupling to the rays is defective
  (checkerboard, seed-dependent amplitude, ifx non-determinism, the 3x stencil constant, the disabled
  field multiply) — the fix is TO's, `BRIEF_to_findings_rulings.md`.
- A grid-surface Reference element between two physical-optics legs imprints its sag on the propagating
  field as phase (the DM phase-imprint validation); the field between legs is advanced by the rays'
  geometric phase.
- `REGRID <elt>` (`ReGrid`, sourcsub.F, called after a near-field leg) rebuilds a REGULAR ray grid at an
  element from the chief ray and the grid pitch, on a plane or a sphere — geometry only; it does not read
  the field.
- With polarization on, each ray carries a complex field `RayE` and the diffraction seed takes it,
  amplitude included; the scalar seed uses a uniform amplitude and `CumRayL`.

## 2. The three workflows

### 2.1 Single screen at the entrance pupil (ray-geometric) — exists, to be fixed

The whole column collapsed to one phase screen at a flat Reference / Obscuring at the pupil.  Right for
the phase structure function, Strehl, long- and short-exposure PSF (the tilt participation), seeing-
limited image quality in the ray trace; exact to the degree that scintillation is negligible
(sqrt(lambda h) small against r0, the aperture averaging the rest).  No amplitude, no anisoplanatism,
no wind.  Placement rule: a flat element at a pupil with the full square ray grid — the screen is a
square grid at the ray-grid pitch read by ray index (the routine's own comment: "this should be at a
pupil!").  Under SEED once (p) is fixed.

### 2.2 Layered atmosphere, rays all the way

Layers as grid-surface Reference elements at altitudes h_i with Kolmogorov / von Karman screens as
their sag, weighted so the r0 compose (r0_i^(-5/3) proportional to Cn2(h_i) dh).  Each ray takes its
OPL from the sag and its tilt from the local surface normal — exact, unwrapped — and marches to the
next layer; the ray grid becomes irregular, which is beam wander and is physical.  Every ray output
(OPD, spot, sensitivities, metrology) sees the atmosphere.  No diffraction amplitude.  Needs no engine
change beyond a screen GENERATOR (2.5).

### 2.3 Layered atmosphere, field all the way

The same layers carried on the grid: plane-to-plane near-field legs between layers and down to the
pupil, then the telescope as today's physical-optics chain to the detector.  Scintillation comes from
the Fresnel propagation between layers; anisoplanatism from the screens being FIXED in global
coordinates (an off-axis field passes each layer displaced by theta h — generate the screens larger
than the pupil by that displacement); time from shifting each layer's origin by v_i dt (frozen flow).
The ray outputs are blind to it; the results are `intensity` / `complex_field`, and the pupil wavefront
is recoverable from the field by unwrapping in weak turbulence.  Sampling: the near-field legs are
within their condition at kilometre distances for metre apertures (dx^2 >> lambda z / N); the pitch
sets the smallest scale carried.  Constraint: vector mode and COMPOSE exclusive, as today.

### 2.4 Rays AND field to the entrance pupil, then REGRID at the pupil — the design (Dave 2026-10-18)

Trace both carriers through the layers (2.2 and 2.3 on the same screens).  At the pupil:

1. Bring the ray OPL onto the regular grid: the rays have wandered by centimetres, so interpolate
   their OPL from the displaced positions to the nominal grid points (linear is enough at that
   displacement).
2. Unwrap the field's phase AGAINST the ray phase, pixel by pixel: unwrapped phase = arg(E) + the
   multiple of 2 pi that brings it nearest the ray's phase.  Local, one pass, no path dependence, no
   error propagation, immune to dark pixels.  The residual (field minus ray phase) IS the Fresnel
   correction — small and smooth in the valid regime, and the diagnostic.
3. Re-seed a regular ray grid at the pupil: L from the unwrapped field phase, tilt from its gradient
   (the mean slope over the pixel, the same approximation the grid makes), amplitude from |E| into
   `RayE` (so the next diffraction seed carries it; the scalar seed would need a per-ray weight).

The hook is `ReGrid` (it already rebuilds the regular ray grid at an element after a near-field leg);
the extension is to take L, tilt and amplitude from the field as above instead of the chief-ray
geometry alone.  Downstream, the ray OPD / spot / sensitivities see the atmosphere AND the next
physical-optics leg sees the scintillated amplitude: what neither pure workflow gives.

**Refusal, built in from the start.**  Where rays cross the geometric phase is multivalued and cannot
guide; where |E| has a zero the field phase has a branch point and no ray exists.  Both arrive with
strong scintillation (Rytov variance above about 0.3).  The regrid counts the pixels whose residual
exceeds pi/2 and REFUSES, naming the count, rather than guessing (rule 5 of PLAN_CONSOLIDATION: a named
refusal, never a silent one).

### 2.5 The generator (a resources runner, Dave's parameterized-runner rule)

The engine's generator is pure Kolmogorov, no outer scale, no subharmonics (the tilt participation is
its patch).  The layered model wants von Karman screens with L0, subharmonics for the low frequencies,
Cn2(h) weighting from a profile (Hufnagel-Valley or a measured one), wind per layer, the field-of-view
margin, written as grid files for the Reference layers and as the pupil screen file for 2.1.  One
runner `<atm>_params` + `<atm>_run` under `templates/`, with the screens' structure function and the
scintillation index of a standard profile as its own gates (Hardy ch. 3 / Roddier numbers).

## 3. Gates for the build

- (a) same screens, workflow 2.4's downstream PSF matches 2.3 to the Fresnel correction; its pupil OPD
  matches 2.2 to the same.
- (b) a run with the residual check tripped refuses and says how many pixels.
- (c) the regrid is bit-identical to today's `ReGrid` when no field is present (scalar mode, no legs).
- (d) the generator's screens: D_phi(r) within 10 % of 6.88 (r/r0)^(5/3) over 1-16 px including the
  subharmonic correction at the largest lags; the scintillation index of the standard profile within
  the published range.
- (e) anisoplanatism: the isoplanatic angle of the standard profile from two field points agrees with
  theta_0 = 0.314 r0 / h_bar.

## 4. Sequence

1. TO: fix (p) — the coupling, the seed, the stencil, the placement refusal, the model memo
   (`BRIEF_to_findings_rulings.md`).  The single screen (2.1) is then usable.
2. Generator runner (2.5), gated by (d).  Enables 2.2 and 2.3 with no engine change.
3. The pupil regrid (2.4): `ReGrid` extended, the refusal, an api entry (`regrid_field`) so the
   bindings drive it; gated by (a)-(c).  Then (e).

Dates: after phase 4 of PLAN_CONSOLIDATION (the release gate); a future activity.
