    ====================================================================
     Small dedicated imager | the found deck ImageDemo.in
     spec: D 350 mm, f 3500 mm (f/10), 2048 x 13 um pixels, band [500 700 900] nm
    ====================================================================

    [1] the deck as found, stripped (36 Zernike terms removed from PRI_MIR; SEC_MIR, EXIT_RET, IMG_PL put back on the axis):
          mirrors R1 933.7886 K1 -1.009724, R2 156.5935 K2 -1.932284, separation 400.2602 mm; D 350, obscuration 0.2669
          paraxial at 700 nm: f 3489.94 mm (f/9.971), chief height 3506.89 mm per tan(theta) on the deck's plane; two-mirror f 3134.47, relay x1.1134
          the paraxial focus is 0.3323 mm in front of the deck's image plane (the figure carried the focus)
          builder vs deck (spiders and the overfilled beam dropped from both), RMS waves at [0 0.109 0.217] deg:
            deck [0.145078 0.157236 0.174872] | builder [0.145078 0.157236 0.174872] | max |diff| 1.7e-10

    [2] refocus by the secondary (CALIB, 24 iterations, floor): separation -6.08 um
          f 3507.66 mm, chief height 3506.80 mm per tan(theta); the paraxial focus now 0.0167 mm from the detector plane
            field (deg)     0.0000  0.0544  0.1087  0.1631  0.2175
            RMS  700 nm     0.0133  0.0100  0.0024  0.0121  0.0295
            distortion %     0.000   0.094   0.378   0.866   1.579
            rays lost           53      53      53      53      53 of 741 (the shadow)
          clear semi-diameters the beam needs at 0.217 deg (paraxial): M1 175.00, M2 26.50, CORS1 12.15, CORS2 12.28, CORS3 12.59, FLTS1 12.92, FLTS2 13.29
          the primary's hole: 13.45 mm radius needed at 0.217 deg, inside the 46.70 mm shadow (clears)

    [3] the glasses: the deck's indices at 700 nm against the engine's catalog:
          n = 1.617830: BaF8 1.617833 (+2.7e-06), then SK10 1.617834 (+3.8e-06)
          n = 1.657863: SF19 1.658044 (+1.8e-04), then BaSF2 1.656392 (-1.5e-03)
          n = 1.535018: LLF2 1.535595 (+5.8e-04), then BaK2 1.535830 (+8.1e-04)
          the refocused design in BaF8 / SF19 / LLF2 over the band:
            field (deg)     0.0000  0.0544  0.1087  0.1631  0.2175
            RMS  500 nm     0.0638  0.0597  0.0494  0.0402  0.0535
            RMS  700 nm     0.0137  0.0105  0.0027  0.0118  0.0293
            RMS  900 nm     0.0010  0.0034  0.0105  0.0199  0.0308
            distortion %     0.000   0.094   0.378   0.866   1.579
            rays lost           53      53      53      53      53 of 741 (the shadow)
          paraxial focus [51.84335 51.77524 51.74508] mm at [500 700 900] nm: spread 0.0983 mm against the +/-2 lambda N^2 = 0.1406 mm depth of focus at 700 nm

    [4] the bare Ritchey-Chretien (closed form at the found m = 6.7134, beta = 0.1008: K1 -1.007387, K2 -1.906616; the deck's -1.009724, -1.932284)
          f 3134.47 mm at its own paraxial focus, the same image heights:
            field (deg)     0.0000  0.0608  0.1217  0.1825  0.2433
            RMS  700 nm     0.0001  0.0489  0.1958  0.4407  0.7836
            distortion %     0.000   0.001   0.004   0.009   0.017
            rays lost           53      53      53      53      53 of 741 (the shadow)

    [5] the camera from its spec: D 350 mm, f 3500 mm, primary f/1.334, BaF8 / SF19 / LLF2, the found corrector x1
          merit = the worst RMS over [500 700 900] nm at [0 0.109 0.218] deg; seed (RC conics, the found corrector) 0.1395 waves
          round 1: CALIB (corrector radii, both conics) 22 iterations, converged; after the re-solve 0.1466 waves -- not the best yet
          round 2: CALIB (corrector radii, both conics) 22 iterations, converged; after the re-solve 0.0447 waves -- KEPT (best so far)
          round 3: CALIB (corrector radii, both conics) 38 iterations, floor; after the re-solve 0.0485 waves -- not the best yet
          round 4: CALIB (corrector radii, both conics) 42 iterations, floor; after the re-solve 0.0479 waves -- not the best yet
          result: two-mirror f/8.8865, BFD 47.3378 mm -> separation 399.7766, R1 933.7886, R2 157.9451, K1 -1.010134, K2 -1.946829
          corrector R [-43.42805 22.94005 174.3436] mm; f 3500.0000 mm (asked 3500), f/10.0000, image 0.1270 mm behind the filter
            field (deg)     0.0000  0.0545  0.1090  0.1634  0.2179
            RMS  500 nm     0.0380  0.0356  0.0310  0.0293  0.0447
            RMS  700 nm     0.0040  0.0071  0.0139  0.0180  0.0099
            RMS  900 nm     0.0139  0.0161  0.0216  0.0259  0.0214
            distortion %     0.000   0.076   0.309   0.710   1.310
            rays lost           53      53      53      53      53 of 741 (the shadow)

    [6] sampling: pixel 13 um = 3.714 urad (f 3500.0 mm); field 2048 pixels = 0.4358 deg full
          Q = lambda N / p at [500 700 900] nm: [0.385 0.538 0.692]; Nyquist sampling is Q = 2, so the camera is undersampled x[5.2 3.71 2.89]
          (lambda/D at 700 nm is 2.000 urad: the pixel is 1.86 lambda/D)
