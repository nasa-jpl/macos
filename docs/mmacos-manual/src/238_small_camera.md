### `small_camera`: the small dedicated imager: a Cassegrain with a refractive field corrector

A single-purpose camera whose telescope and instrument are one design: two mirrors, a corrector behind the primary's hole, a filter and one detector.  The relay's power is part of the camera's focal length.  That is what sets this family apart from an observatory front end feeding a `Bench` (the mmacos manual's chapter 2 thread).

The worked case is the camera whose deck the MACOS manual carries as `ImageDemo.in` (Appendix A): a 350 mm, f/10 Ritchey–Chrétien imager with a 2048-pixel line array of 13 µm pixels (3.7 µrad per pixel, a 0.44° field) and a 500–900 nm band.  These design values are the runner's defaults, and the stripped deck is checked against them.

The design layer under it is `macos.design.cass_camera` and `macos.design.camera_emit`.  The mirrors come from `Telescope`'s closed form, or are given outright.  Relay elements are `macos.design.Lens` objects.  The y-nu trace runs through the folded train, and the deck is scored with `macos.design.lens_eval`.

The runner has six stages:

1. **The deck as found, stripped.**  It removes the primary's Zernike figure and puts the secondary, the exit pupil and the image plane back on the axis.  The builder then rebuilds the stripped deck from the deck's own numbers, and the two decks trace alike.
2. **Refocus.**  The figure carried the focus, so without it the paraxial focus falls inside the filter.  The secondary's spacing brings the image back to the deck's detector plane (CALIB on that one DOF).  The stage then reports the field at the design wavelength, the clear apertures the beam needs, and whether the primary's hole clears.
3. **The glasses.**  The deck's fixed indices are matched to the engine's catalog (`macos_glass_list.txt`), and the band is traced with the matched glasses' dispersion.
4. **The bare telescope.**  This is the closed-form Ritchey–Chrétien with the same m and β and no relay, which shows what the corrector is for.
5. **The camera from its spec** (aperture, focal length, pixel, band):
   - the RC closed form;
   - the found corrector, scaled and in catalog glass;
   - the first-order solve that sets the focal length and puts the focus behind the filter;
   - CALIB rounds over the band and the field.
6. **The sampling**: λN/p at each wavelength.  The camera is undersampled by design, and the stage says by how much.

**Why rounds, and why "best kept".**  CALIB has no focal-length row.  Its corrector radii move the relay's magnification, so each round ends with the first-order re-solve: the two-mirror f-number and back focal distance, by two nested 1-D solves, put the focal length and the focus back on the spec.

The first re-solves shift m and undo part of the balance.  The rounds settle as the re-solve's move shrinks.  The design kept is the best one seen by an independent merit: the worst `lens_eval` RMS over the band at the CALIB fields and the full field, among designs that lose no ray beyond the shadow.  The LM merit does not charge lost rays, as the lens templates found.  A failed CALIB, a failed re-solve, or a round that loses rays ends the rounds.

**Run it.**

```matlab
cd('<path-to>/mmacos/templates/36_cameras/small_camera');
out = small_camera();
```

**The knobs** (`small_camera.m`):

<!-- transcript: p2_small_camera_knobs -->

**What it prints:**

<!-- transcript: p2_small_camera_run -->

![](media/p2_small_camera_fig.png)

**FIGURE 64** RMS wavefront error (log scale) and distortion against field for the small dedicated imager of ImageDemo.in: the found design, stripped of its primary's figure and refocused by the secondary, in its matched catalog glasses at 500, 700 and 900 nm; the camera re-designed from its spec; and the bare Ritchey-Chretien without its corrector (`small_camera_field.png`).

**The number it reproduces.** Stripped of the primary's figure, the deck's focus falls %pin:p2_small_camera_zfocus% mm short, inside the filter; a %pin:p2_small_camera_dsep% um despace of the secondary restores it.  Re-designed from its spec (D 350 mm, f 3500 mm, 500-900 nm), the camera's worst RMS over the band and the field is %pin:p2_small_camera_merit% waves, and at 700 nm the pixel samples the image at Q = lambda N / p = %pin:p2_small_camera_q700%, under Nyquist by design.

**Adapt it.** *"Using templates/36_cameras/small_camera, design the camera for a 250 mm aperture at f = 2500 mm with 10 um pixels over 450-850 nm, and tell me the field curve at each wavelength, the distortion at the edge, which CALIB rounds were kept and how undersampled the detector is."*

