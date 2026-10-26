# BRIEF: an IODC lens design problem as an mmacos challenge -- the verified statements and the build (CC, 2026-10-23)

Source: the SPIE papers Dave pulled to `~/dev/MACOS_sandbox/IODC/` (R. C. Juergens; 2006 with P. K. Manhart).  Read
from the text layer 2026-10-23; the figures are layouts at 1/4 scale, NO prescriptions are printed in any of them, so
every comparison is against the published MERIT values and the top-ten tables, never against a copied lens.
Entries were CODE V, Zemax, OSLO and in-house programs (2014: 17 CODE V, 15 Zemax, 4 OSLO).

## The statements, verbatim where it matters

| | 2010 "the green lens" (76522I) | 2014 "the Cinderella lens" (929303) | 2021 "the down under lens" (1207820) | 2023 "the Yin-Yang lens" (1279837) | 2006 "the lens shuffler" (63421U) |
|---|---|---|---|---|---|
| Focal length | magnitude 100 mm (intermediate images allowed) | magnitude 100.00 +- 0.01 mm | magnitude 100.00 +- 0.1 mm, both versions | 100 +- 0.1 mm at 550 nm | 100 +- 0.1 mm paraxial, both lenses |
| Wavelength | monochromatic 532 nm | monochromatic 587.56 nm | monochromatic 587.56 nm | 540 and 550 nm | 587.56 nm |
| Glass | Schott N-BK7 | Schott N-BK7 (n = 1.5168) | Schott N-BK7 | Schott N-BK7 and N-SF6 | N-BK7 |
| Surfaces | spherical or plano; every spherical surface has the SAME radius magnitude, either sign | spherical or plano | spherical or plano | spherical or plano; no cemented surfaces; no intermediate images | spherical or plano |
| The twist | one radius | all elements cut from ten N-BK7 blanks 100 mm dia x 30 mm thick (multiple lenses per blank allowed, any dissection, no material lost) | the lens must meet spec flipped end to end with the radii NOT reversed | diffraction-limited over the field at 550 nm AND bad on axis at 540 nm at the SAME image location; length <= 500 mm | one 100 mm N-BK7 cylinder divided into TWO lenses that shuffle back into the cylinder |
| Vignetting | none: the stop fully filled over the field | stop fully filled over the field | stop fully filled, both versions | stop fully filled at 550 nm over the field, on axis at 540 nm | none |
| Image quality | RMS WFE <= 0.070 wave over the field, "chief ray reference; tilt and focus not removed" | RMS WFE <= 0.070 wave over the field, "piston and tilt removed; focus not removed" | the same as 2014, both versions | diffraction-limited at 550 nm over the field (the limit in the paper's text) | RMS WFE <= 0.07 wave |
| Distortion | <= 5.0 % everywhere | no requirement (entries ran to -100 %) | no requirement | not specified | <= 5.0 % |
| Object / image | infinity / flat, real, in air | infinity / flat, real, in air | infinity / flat, real, may differ between versions | infinity / flat in air, same location for both wavelengths; paraxial image distances equal within 0.0001 mm | infinity / flat, outside the cylinder |
| Stop | anywhere | anywhere | anywhere, may differ between versions | anywhere | floating |
| Merit | D[mm] x semi-field[deg] | D x semi-field; the semi-field is where RMS exceeds 0.070 OR vignetting starts, whichever first | D x semi-field, common to both versions | D x semi-field x on-axis RMS WFE at 540 nm (higher is better) | D x theta, same for both lenses |
| Entries / winner | 37 from 8 countries; EPD 81.3 mm, 43.5 deg, 3537 | 45 from 13; Jon Ehrmann 25 elements, 33.90 mm, 62.50 deg, 2119 (the only entry with an intermediate image; 870.7 mm long; f/2.95 on axis; -53 % distortion); 2nd Hoffman 64.80 mm, 28.46 deg, 1844 | 43 from 13; Damien Gawron 198.4 mm, 90 deg, 17857 | 24 from 7; Arnaud Davenel; merit values 13 to 152,013 | -- |

Known loopholes the 2014 paper names (so our rules must state them): zero axial / edge thicknesses allowed; no
mounting material beyond the clear aperture; no distortion requirement; intermediate images allowed.

## What the statements change in the ranking (CC)
- **2010 first.**  Its wavefront convention, "chief ray reference; tilt and focus not removed", is MACOS's native OPD
  (the chief-ray reference of notes/opd_reference.md, evaluated on the fixed image plane) -- no evaluator argument to
  settle.  The single-radius rule makes the design a small search: one radius magnitude R, a sign per surface, the
  spacings, the stop -- CALIB on the spacings (and R) inside a sign / element-count search a MATLAB driver owns.
  Distortion <= 5 % is a stated, checkable limit.  37 entries, winner 3537, a published top table.  It is also the
  cheapest refractive stress of the engine: N-BK7 at 532 nm through the glass table, 43.5 deg semi-field with every
  stop ray kept, |EFL| = 100 with the intermediate-image option.
- **2014 second.**  "Piston and tilt removed" is NOT MACOS's native reference (a best-fit tilt must be removed from
  the OPD map per field -- the Zernike fit, modes 2-3, exists); the blank-packing rule is a result scorer we write;
  the winner is an intermediate-image form the others did not find.  Harder, richer, 45 entries.
- **2023 as the stress test**: two glasses, two wavelengths, a merit that REWARDS on-axis error at 540 nm while
  demanding diffraction limit at 550 nm -- an outer loop over D, field and the chromatic balance around CALIB.
- 2021 (the flip symmetry at 90 deg) and 2006 (the shuffler) as later rounds; 2017 excluded.

## The build (the chosen year becomes `mmacos/challenges/iodc<year>/`)
1. **State the metric before the number** (challenges/README.md): the wavefront reference used, the field sampling,
   how the RMS limit sets the semi-field (RMS or vignetting, whichever first), the distortion check, the paraxial EFL
   check.  Cross-check MACOS's RMS WFE against CODE V on ONE lens (Luis's CODE V, or the rx_converter round trip) so
   the evaluator is argued from the definition, not transcribed.
2. **The lens as a design-layer System / Bench** emitting the .in (refractors with the glass from the table at
   the problem's wavelength; apertures = the blank diameter; the stop as an Obscuring/aperture element so "fully
   filled" is a ray count, not a hope).  Glass: `macos_f90/macos_glass_list.txt` (201 Sellmeier lines) carries
   Schott's BK7 and UBK7 but NOT N-BK7 -- add the N-BK7 (and, for 2023, N-SF6) Sellmeier line from Schott's data
   sheet with its provenance in the comment; BK7 and N-BK7 agree to 1e-4 in n_d but the entry must be the named glass.
3. **CALIB** on spacings (and curvatures where the year allows) for the WFE target across the field.  CALIB's
   `VarElt=` mask has 8 DOFs per element: 1-6 rigid body (axial position = the spacing), 7 = radius, 8 = conic
   (notes/calib.md), so curvature and spacing are available without an engine change; `OptAsph` is not needed
   (spherical only).  For 2010 the shared |R| is ONE variable: drive it from the MATLAB loop and let CALIB solve the
   spacings, or add a linked-radius DOF later if the scan is slow.
4. **The outer loop** (MATLAB, `iodc<year>_params.m` + `iodc<year>_run.m`): for each (D, theta) candidate, the
   CALIB solve; the merit boundary where RMS = 0.070 or vignetting starts; the published winner as the target line
   on the (D, theta) plot (the paper's Figure 29 form).
5. **Gates**: the evaluator against the CODE V cross-check (1e-3 wave); a must-fail leg (a lens with a vignetted
   stop is refused by the scorer); the challenge's numbers reproduced by a re-run.
6. **Report** the way the README says: merit, D, theta, element count, length, distortion; beside the published
   top ten.  Any engine finding on the way -> PLAN_CONSOLIDATION 1.2.

**RULED 2026-10-26 (Dave): the revised order -- 2010 first, 2014 second, 2023 the stress test.**  Owner: TO's lane, after
round 2 of the mmacos manual; the chosen year's directory is `mmacos/challenges/iodc2010/`.  First deliverable, before any
optimization: the EVALUATOR (step 1), cross-checked on one lens against CODE V; its gate is the challenge's own proof.
