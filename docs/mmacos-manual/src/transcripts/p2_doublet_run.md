    ====================================================================
     Cemented achromat | f = 100 mm f/5 | N-BK7 (V 64.17) + F2 (V 36.37) | t = [6 2.5] mm
    ====================================================================

    [1] thin lens: phi1 = 0.023081, phi2 = -0.013081 /mm (f1 = 43.326, f2 = -76.446 mm)
          partial dispersions P = (n_d - n_C)/(n_F - n_C): 0.3076, 0.2938
          textbook secondary spectrum f (P1 - P2)/(V1 - V2) = 0.0499 mm

    [2] thick-lens achromat, equi-convex crown (y-nu, EFL 100.000000 at d, BFL_F - BFL_C = 1.4e-14 mm):
          R = [44.7813 -41.5473 -701.5949] mm
          RMS at F / d / C: 0.2672 / 0.0239 / 0.0497 waves

    [3] CALIB on R1..R3 at F, d, C (158 iterations, converged):
          R = [46.4388 -42.3535 -476.9006] mm, EFL 99.7644 at d
          RMS at F / d / C: 0.0403 / 0.0759 / 0.0509 waves (paraxial achromat: 0.2672 / 0.0239 / 0.0497)
          contrast, CALIB at d alone (3000 iterations, cap): RMS F/d/C 0.0922 / 0.0056 / 0.1270, BFL_F - BFL_C -0.0877 mm

    [4] focal shift, back focal distance minus its d value (mm), one engine trace per wavelength:
          paraxial achromat:  F +0.04895   C +0.04895   (y-nu: F +0.04895  C +0.04895)
          CALIB lens:         F -0.00347   C +0.06869   (y-nu: F -0.00347  C +0.06869)
          secondary spectrum (paraxial achromat, d focus - F/C focus): engine -0.04895 mm, textbook -0.04993 mm (thin lens)
          worst RMS across the band: CALIB 0.0759, paraxial achromat 0.2672 waves

    [5] PAC095AR.14 (Newport): R = [140.890 -101.720 -454.867 ] mm, t = [14.50 4.00 ] mm, N-BK7 / N-F2
          at 589 nm: EFL paraxial 250.322 / engine 250.322 / vendor 250 mm; BFL paraxial 239.883 / vendor 240.1 mm

    DONE.  Report <work>/templates/35_lenses/doublet/out/doublet_report.txt
