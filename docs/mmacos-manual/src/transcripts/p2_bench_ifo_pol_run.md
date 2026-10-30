    Building Twyman-Green rig...
      test arm: coated BSrefl = elt 4, Recomb = elt 12
      ref  arm: coated BScrefr = elt 8, Recomb = elt 10

    === GATE 1: single-surface 45-deg Fresnel analytic (BS reflect) ===
      3001 rays good at BS
      BS mean AOI            : 45.000 deg
      RS/RP magnitude resid  : < 1e-12
      RS/RP phase resid (rad): 3.318e-14
    GATE 1 PASS: coated BS reflection matches Fresnel to < 1e-12

    Building arm Jones pupils at recombination...
      test leak=5.47e-23  ref leak=4.82e-23 (longitudinal residual)
      3001 rays good in both arms

    Arm-differential polarization at recombination (3001 rays):
      Diattenuation D : mean=7.2133e-02  var(rms)=2.0099e-08  max=7.2133e-02
      Retardance (rad): mean=8.3521e-02  var(rms)=3.5997e-08  max=8.3521e-02

    PSI phase-error contribution (co-pol fringe, x-pol input):
      RMS phase : 2.2848e-08 rad = 3.6364e-09 waves = 2.3011e-06 nm @ 632.8 nm

    === GATE 2: pol-off bit-identity (coated vs uncoated) ===
      test-arm OPD max-diff: 0.000e+00 mm
      ref-arm  OPD max-diff: 0.000e+00 mm
    GATE 2 PASS: coating is inert when polarization is off (<1e-12 mm)

    === ALL GATES PASS ===
