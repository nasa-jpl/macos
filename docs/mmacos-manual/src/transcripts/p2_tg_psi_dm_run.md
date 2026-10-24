    === GATE 1: Tranche-1 element order + the grid carries the pol train ===
      tg_test.in   polarizing elements at [4 8 10 16 17]; 0 physical-optics legs
      tg_ref.in    polarizing elements at [4 7 9 14 15]; 0 physical-optics legs
      grid power at analyzer 90 / at 0 = 0.5897  (1.0 would mean the grid is blind)
    GATE 1 PASS

    === GATE 2: detector registration between the two arms ===
      pixel pitch  test 0.068842384 mm | ref 0.068842386 mm | relative diff 3.11e-08
      flux-weighted centroid offset 1.11e-06 px | second-moment radius diff 1.94e-07 px
    GATE 2 PASS: pitch matches to 3.1e-08, beams registered to 1.1e-06 px

    === ALIGNMENT: measure the arms, then null the beamsplitter rotation ===
      reference arm azimuth at recomb      +45.0000 deg (design +45.0000)
      test arm azimuth at recomb, QWP @0   -37.5214 deg (design -45.0000)
      --> the arms are 7.4786 deg from orthogonal; the beamsplitter has
          rotated the TEST arm by 7.4786 deg
      solved test-arm QWP azimuth  +3.7680 deg -> arm azimuth -45.0000 (want -45.0000), 5 traces
      output QWP azimuth set to +90.0000 deg (bisector + 45)
      after the output QWP: test |b/a| 1.00000 arg  -90.00 | ref |b/a| 1.00000 arg  +90.00
      (circular is |b/a| = 1 and arg = -+90);  |<test|ref>| = 2.98e-07
    ALIGNED

    === GATE 3: analyzer-basis synthesis is exact (and why) ===
      9 traces (3 per arm x 3 decks) in 3.55 s -- every analyzer angle now free
      theta   23.70 deg : synthesized vs direct trace, rel 7.541e-11
      theta   71.30 deg : synthesized vs direct trace, rel 7.337e-11
      theta  137.20 deg : synthesized vs direct trace, rel 3.507e-10
      theta  -18.00 deg : synthesized vs direct trace, rel 1.422e-10
      non-vacuity: the O(beta^2) law (beta = ray angle at the analyzer)
        test optic Kr= 4000 : beta 1.1233e-02 rad, rel 1.8397e-05, rel/beta^2 0.1458
        test optic Kr= 2000 : beta 2.2467e-02 rad, rel 7.3446e-05, rel/beta^2 0.1455
        DM deck          : beta 2.2728e-05 rad (the DM slope error itself)
      quadratic-law check: (rel ratio)/(beta ratio)^2 = 1.0020  (1 = exactly 2nd order)
    GATE 3 PASS: exact to 3.5e-10 collimated; degrades as 0.149*beta^2 with power

    === GATE 4: analyzer sweep modulates as cos(2 theta) ===
      from 12 TRACED angles (9.29 s): visibility 0.998310 | 4t/2t 8.8951e-04 | 6t/2t 2.9871e-14
      from the synthesized basis        : visibility 0.998310 | 4t/2t 8.8951e-04 | 6t/2t 5.0647e-17
      (the synthesized 6t is zero BY CONSTRUCTION -- the traced one is the gate)
    GATE 4 PASS

    === GATE 5: a grid value is a SURFACE height (sign + scale) ===
      uniform grid piston 20 nm -> +0.397167 rad (expect +0.397167), gain 1.00000
      rigid optic shift   20 nm -> +0.397167 rad
      grid vs rigid shift  8.730e-11 rad  => a grid value IS the surface height
      scale error 0.00002 nm of surface (0.0001%), pupil rms 5.74e-04 rad
      COUNTERFACTUAL, same rig with the waveplates at their DESIGN angles:
        gain 1.11661 (11.7% high) with visibility 0.996612
        aligned: gain 1.00000 with visibility 0.998304
        -> the misalignment costs 11.7% of SCALE and 0.17% of CONTRAST; a fringe
           monitor is ~69x less sensitive to it than the measurement is
    GATE 5 PASS: h = +1 * psi * lambda/(4 pi), gain 1.00000

    === MEASUREMENT: four-step polarization PSI (differential) ===
      fringe visibility 0.9983 | recovered surface 6.271 nm rms in the pupil
      four-step vs 64-angle least squares: < 1e-10 nm rms of surface
      (the 4-theta systematic is COMMON to the poked and baseline runs, so
       the differential cancels it -- it survives only in a single-shot map)

    === CLOSURE: recovered surface vs the injected DM map ===
      pupil image: mag 10.5252 (det->DM), anamorphic 0.000%, rotation -180.000 deg,
                   nonlinear distortion 0.0068 mm rms (0.03% of the beam)
      pixel axes x=-col,y=-row, correlation 0.998434
      truth 6.35 nm rms | recovered 6.27 nm rms
      RESIDUAL 0.363 nm rms whole pupil | 0.304 nm rms interior (rim excluded)
    CLOSURE PASS: 95.3% of the injected map recovered (interior)
