    ==================== e2e6m S1 -- the 6 m unobscured telescope
    D = 6.00 m | lambda = 500 nm | model 256 | grid 41
    metric: RMS WFE @ 500 nm, exit-pupil referenced (add_pupil FarField),
            on the stage's fitted focal plane; raw keeps tilt, -tilt removes
            piston+tip/tilt.  TELESCOPE-ONLY config, anchored at the
            TELESCOPE best-focus exit pupil.
    gates : f/# in [12 20] | shroud <= 8.0 m (deployed, diameter-only)
            | UNOBSCURED | DL <= 0.071 waves | AOI spread < 15 deg (preference)
    apertures OFF (rule 1); clearance from check_clipping footprints

    [1] all-sphere layout  R = [38.400 3.800 2.000] m, T = [17.500 14.000] m,
        tilt = [-5.60 5.90 8.00] deg -> uncorrected 8219 waves

    [2] staged Zernike (freeform) correction, waves RMS @ 500 nm:
        full-field footprint radius (m): 3.0164 0.2972 0.4162
        lmon mode "auto" -> 3.4688 0.3418 0.4786
        modes (BornWolf): [3 4 5 9 10 11 12 13 19 20 21 22 23 24 25]
        S0 centre        : 8345 -> 0.8594
        S1 inner 2-D     : worst 166.503 -> 0.3512
        S2 full 2-D +-0.35' : worst 0.0702, area-weighted 0.0430 -> DIFFRACTION-LIMITED
        EFL after correction 70.573 m -> f/11.76 (spheres alone: see [4])

    [3] fitted FP from a 5x5 field-foci grid: tilt 0.835 deg,
        defocus removed +1.332 mm, field-curvature sag -6.1 to +6.5 um

    [4] first order + packaging
        EFL 70.579 m -> f/11.76  [FAIL]
        shroud 7.332 m diameter (1.222 x D), train 17.37 m  [PASS]
        clearance 4/4 bodies clear  [PASS]
        dense 7x7 map over +-0.35': raw max 0.1059, -tilt max 0.0930 waves  [FAIL]
        AOI spread per mirror: 8.9 9.8 3.3  deg (preference < 15)
        per-element radial extent (m): M1 2.970   M2 3.666   M3 3.442   FP 3.218   

    [5] saved <work>/templates/80_end_to_end/e2e6m/s1_telescope.in (+ .mat)
        RELOAD GATE (rule 3): 6 elements, 1185/1185 rays pass  [PASS]

    [6] figures: layout, WFE field map, FP curvature, shroud fit
        shroud figure measures 7.361 m against the 8.0 m gate  [PASS]
    check_clipping  (family=tma, 6 elements)
      name       kind            ap_r    foot_r    margin  clearnce obstruct  status
      M1         Reflector        3.3     2.999    0.3006        --        0  OK
      M2         Reflector          3    0.2738     2.726    0.1649        0  OK
      M3         Reflector          3   0.03424     2.966    0.2821        0  OK
      FP_return  Return           3.3 5.698e-07       3.3        --        0  OK
      ExitPupil  Return           3.3   0.02268     3.277        --        0  OK
      FP         FocalPlane       1.8 6.045e-07       1.8 -5.976e-07      159  ** CLIP
      => layout has CONFLICTS: margin<0 = own beam overfills aperture; clearance<0 = body cuts a foreign beam

    aoi_report  (limit 15 deg, 1 field(s))
      mirror   chief(deg)        min        max  SPREAD(deg)   ok
      M1             5.66       1.15      10.03         8.88   yes
      M2             5.68       0.69      10.49         9.79   yes
      M3             8.55       6.92      10.24         3.33   yes

    [7] design report -> s1_design_report.txt

    S1 DONE in <t> min
