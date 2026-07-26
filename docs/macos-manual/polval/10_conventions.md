<!-- GENERATED FILE -- do not edit.
     Source: 10_conventions.md.in
     Numbers: generated/numbers.json (MATLAB driver)
     Regenerate: make polval-regen  (docs/macos-manual)
-->
# Conventions

Cross-code polarization disagreements are almost always convention
mismatches. These were pinned during the Phase-0 audit by *deriving them from
the engine source*, not by legislating them, and they are asserted in the
tests. They are reproduced here because no measured number below is
interpretable without them.

| Convention | Value | Where it comes from |
|---|---|---|
| Time-harmonic phase | `exp(+iωt)`; spatial propagator `exp(−ikz)` | `elemsub.F:387` (`C1 = exp(−i·2π·L·N/λ)`, phase decreases as path grows) and the coating recursion at `elemsub.F:512-516`. Consistent with the independently established interferometer result that field phase advances as optical path shortens, and with the `opd_sign_flip` the pymacos↔PROPER comparison already carries. |
| Absorbing index | `N = n − iκ`, with κ > 0 = loss | As stored in `IndRefArr`/`ExtincArr` and applied as `DCMPLX(n,−κ)`. |
| Jones storage basis | Linear (x, y); circular by unitary change of basis | Phase-2 decision — storing linear and converting is strictly simpler than the reverse. |
| Jones reference frame | Double-pole (Chipman) by default; local s/p and global as diagnostics | The naive local s/p basis carries a coordinate singularity that appears as spurious retardance. Quantified in the basis-artifact gate below. |
| Retardance sign | Fast axis leads | Stated so that vector-vortex handedness (Phase 3) has a fixed reference. |
| Diattenuation | Intensity-based, `D = (T_max − T_min)/(T_max + T_min)` ∈ [0,1] | |
| Coating thickness | Rx `Coating=` layer thickness is **waves at parse-time `Wavelen`**, converted to physical at load; the trace applies phase at the *current* wavelength | `msmacosio.inc:2660`. Broadband sweeps are therefore already correct. `Coating=` must follow `IndRef=` in the file (the parser snapshots the boundary media). The `coat_set` API takes **physical** thickness and sidesteps all of this. |

Two coating subsystems exist in the engine and must not be conflated. Model A
(`Coating`/`EltCoat`/`IndRefArr`/`ExtincArr`/`EltCoatThk`, complex index) is
the one that drives the polarization path and the one `coat_set`/`coat_get`
target. Model B (`nCoatElt`/`CoatIndxElt`/`CoatThkElt` + `AirGap`, real index)
serves the non-sequential refractive path only. All results in this document
are Model A.

## Engine defects found by this validation work

Building the validation layer is what exposed these. Each was invisible to
intensity-level testing, each is now pinned by a named test, and each is
recorded here because a reviewer comparing MACOS against another polarization
ray-trace code needs to know which engine revision the comparison is valid
for.

**1. Stale `RayE` on a polarization-state change.** The `POLARIZATION`
command changes trace-relevant state — `Ex0`/`Ey0` seed `RayE` at source-grid
setup — but did not reset the cached trace. A polarization-state change
followed by a re-trace harvested the *previous* state's field. Fixed by
having `pol_set`/`vecdif_set` call `modified_rx`, the standard
dirty-the-trace convention for setters. This is a prerequisite for the
Jones-pupil layer, which by construction traces twice with different input
states.

**2. Coated-branch incident medium.** The `ifPol` + coated recursion read the
incident index from the parser's `IndRefArr(0,iElt)` slot, which for a coated
mirror *following* another mirror holds the previous mirror's
perfect-conductor substrate — i.e. light modeled as arriving from inside a
perfect conductor. Fixed to use the medium the ray actually travels in, which
is what the uncoated branch always used.

**3. Coated-branch signed incident cosine.** The incident cosine
`ccfb_arr(0)` is negative when the surface normal faces the beam, which turns
every interface coefficient into its reciprocal: |R| > 1, s and p roles
swapped, retardance sign flipped — while |D| *survives*, which is precisely
why intensity-level tests never caught it. Fixed with an absolute value,
mirroring the uncoated branch. Diagnostic signature if it ever returns: the
measured-to-analytic `RS/RP` ratio becomes exactly `(RP/RS)²`.

Defects 2 and 3 are both pinned by the Fresnel-analytic fold gate below, at
the 1.20e-14 level. Defect 1 is pinned by the Phase-1 state tests.
