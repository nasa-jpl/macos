    === the coating ===
      prism index (MacNeille condition)  n_g = 1.655468
      internal angles: theta_H 29.8760  theta_L 60.1240  (sum 90.000000, Brewster)
      9 layers, symmetric, 8 quarter waves
      layer thickness (nm): 38.82 235.25 77.64 235.25 77.64 235.25 77.64 235.25 38.82 
      ANALYTIC (Macleod characteristic matrix, glass in / glass out):
        symmetric (1/2H L 1/2H)^4 : R_s 0.99958011  R_p 2.601e-30  T_p 1.00000000
        odd-QW    H(LH)^4         : R_s 0.99981430  R_p 2.112e-02  T_p 0.97888307
      -> Brewster at the H/L interfaces equalizes the tilted p admittances
         (both 2.7101), so for p the whole stack is ONE HOMOGENEOUS SLAB.
         Its two boundaries with the prism are NOT Brewster, so what is
         left depends only on the slab's total p phase: 8 quarter waves
         is a half-wave ABSENTEE (R_p = 0); 9 is a quarter-wave layer
         (R_p = 2.11%).  Both satisfy the textbook condition.

    DM map: 5.598 nm rms, 65.717 nm PtV, 89.2 mm span, 2.98 mm influence radius

    v2 test arm 22 elements | cube faces+diagonal at [5 6 7] (out) / [11 12 13] (back)
    v2 ref  arm 22 elements | [5 6 7] (out) / [11 12 13] (back)
    v1 test arm 23 | v1 ref arm 21 elements (the plate rig, for comparison)

    arm geometric path difference: 2.274e-13 mm (no compensator plate in the rig)

    === GATE 1: Coating= round-trips, and Tranche-1 order holds ===
      v2_test.in   pol elements at [4 8 10 15 16]; 0 physical-optics legs
      v2_ref.in    pol elements at [4 8 10 15 16]; 0 physical-optics legs
      diagonal (elt 6): 9 layers | thickness round-trip 1.75e-16 | index 0.00e+00
      face     (elt 5): 1 layer  n 1.380, 114.64 nm (MgF2 quarter wave = 114.64 nm)
      grid power at analyzer 90 / at 0 = 2.4855e+11  (1.0 would mean the grid is blind)
      detector pitch test 0.068838097 mm | ref 0.068838097 mm | relative diff 0.00e+00
      centroid offset 1.75e-06 px | second-moment radius diff 2.72e-11 px
    GATE 1 PASS

    === GATE 2: engine coated diagonal vs the Macleod analytic ===
                           engine             analytic     rel err
         T_s   4.198875192923e-04   4.198884783816e-04    2.28e-06
         T_p   9.999999999960e-01   1.000000000000e+00    3.97e-12
         R_s   9.995801124807e-01   9.995801115216e-01    9.59e-10
         R_p   3.976191109163e-12   2.600962280092e-30    3.98e-12 (abs)
      ENERGY (the two arms' decks measure R and T of the SAME coating):
        R_s + T_s = 1.000000000000      R_p + T_p = 1.000000000000
      extinction  T_p/T_s = 2382 : 1   (a real MacNeille cube is ~1e3)
      NON-VACUITY, the odd-QW stack: engine R_p = 2.1116e-02 (analytic 2.1117e-02), T_p = 0.978884
    GATE 2 PASS

    === GATE 3: arm states at the recombination plane ===
      v2 cube : test  -89.999997 deg | ref   -0.000002 deg | from orthogonal -5.264e-06 deg
      v1 plate: test  -37.521390 deg | ref  +45.000000 deg | from orthogonal -7.4786 deg
      -> the cube removes the beamsplitter rotation by a factor 1.42e+06
      after the output QWP (fixed at 45 deg): test |b/a| 1.000000 arg  -90.00
                                              ref  |b/a| 1.000000 arg  +90.00
      circular is |b/a| = 1, arg = -+90;  |<test|ref>| = 8.17e-08
    GATE 3 PASS

    === GATE 4: the analyzer basis is still exact on the cube rig ===
      9 traces (3 per arm x 3 decks) in 3.71 s
      theta   23.70 deg : synthesized vs direct trace, rel 7.553e-11
      theta   71.30 deg : synthesized vs direct trace, rel 7.348e-11
      theta  137.20 deg : synthesized vs direct trace, rel 3.512e-10
      theta  -18.00 deg : synthesized vs direct trace, rel 1.424e-10
      12 TRACED angles (9.77 s): visibility 1.000000 | 4t/2t 8.895e-04 | 6t/2t 2.987e-14
      synthesized             : visibility 1.000000 | 4t/2t 8.895e-04 | 6t/2t 4.753e-17 (0 by construction)
    GATE 4 PASS

    === GATE 5: PSI scale on a known 20 nm surface piston ===
      grid piston 20 nm -> +0.397167 rad (expect +0.397167)  gain 0.999999, err 0.0000 nm
      rigid optic shift    -> +0.397167 rad  (grid vs rigid 4.95e-10 rad)
      visibility 1.000000

      --- arm-waveplate azimuth error: what each rig does with it ---
        eps |  v2 |sep|-90 (ref-arm err) |  v1 |sep|-90 |    v2 gain     v2 vis
        0.0 |   -5.264e-06   -5.264e-06 |      -7.4786 |   0.999999   1.000000
        2.0 |   -2.337e-06   -7.093e-02 |      -3.5520 |   0.999999   1.000000
        3.8 |   -8.459e-06   -1.342e-01 |      -0.0001 |   0.999999   0.999980
        5.0 |   -1.275e-05   -1.789e-01 |      -2.5215 |   0.999999   0.999926
       10.0 |   -3.101e-05   -3.692e-01 |     -13.1490 |   1.000000   0.998447
      -> the cube's REFLECTING arm is cleaned to r_p (= 0): 3.1e-05 deg at 10 deg of error.
         Its TRANSMITTING arm is cleaned only to t_s/t_p = 0.0205: 0.369 deg.
         v1, with no cleanup at all: 13.15 deg.
         (the v1 column is NON-MONOTONIC on purpose: v1 STARTS 7.479 deg off, so
          its curve dips through a minimum near eps = 3.768 -- which IS the
          waveplate clock v1 has to solve for.  v2 has no such minimum to find.)
      -> and in v2 the error costs CONTRAST (1.000000 -> 0.998447), not SCALE
         (gain 0.999999 -> 1.000000).  That is the INVERSE of the v1 finding:
         the PBS turns an invisible systematic into a visible one.
    GATE 5 PASS: h = +1 * psi * lambda/(4 pi), gain 0.999999

    === GATE 6: output-port efficiency ===
      detected power  v2 cube 0.383616 | v1 plate 0.169079 | ratio 2.2689
      per-arm coating budget T_p*R_s = 0.999580, four AR faces at T = 0.995110 each
        -> predicted per-arm throughput 0.980172 (bare glass would be 0.777342)
    GATE 6 PASS

    === MEASUREMENT: four-step polarization PSI (differential) ===
      fringe visibility 1.000000 | recovered surface 6.259 nm rms in the pupil
      four-step vs 64-angle least squares: < 1e-10 nm rms of surface

    === CLOSURE: recovered surface vs the injected DM map ===
      detector->DM magnification 10.5258, rotation +0.00 deg, anamorphic 0.0000%,
      nonlinear distortion 0.0068 mm rms; axis convention x=-col,y=-row
      truth 6.35 nm rms | recovered 6.26 nm rms | correlation 0.999573
      residual 0.183 nm rms interior (0.207 whole pupil)

    =========================================================
     SUMMARY -- v1 plate (ideal polarizers) vs v2 MacNeille cube
    =========================================================
                                              v1 plate        v2 cube
     arm rotation from orthogonal (deg)        -7.4786      -5.26e-06
     PSI scale gain, as designed               1.11661       0.999999
     fringe visibility                        0.996612       1.000000
     delivered power (arb, same source)         0.1691         0.3836
     sep after 10 deg WP error (deg)           13.1490         0.3692
     alignment step needed              yes (+3.768 deg)           none
     compensator plate                        required      2.274e-13

     coating: engine vs Macleod  R_s 9.59e-10 | T_p 3.97e-12 | R_p 3.98e-12 (null)
     DM closure: corr 0.999573, residual 0.183 nm rms interior on 6.35 nm
