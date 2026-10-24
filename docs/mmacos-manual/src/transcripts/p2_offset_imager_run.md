    =================================================================
     S1  coaxial symmetric-asphere solve, on-axis box
    =================================================================
      oi_solve S1  : 15 vars, 9 fields, start  13461.974 nm
        it  1:    575.770 nm   (lam 3.7e-05)
        it  2:    345.616 nm   (lam 1.4e-06)
        it  3:     47.374 nm   (lam 1.0e-07)
        it  4:     25.755 nm   (lam 1.0e-07)
        it  5:     23.812 nm   (lam 1.0e-07)
        it  6:     23.694 nm   (lam 1.0e-07)
        it  7:     23.691 nm   (lam 8.0e-07)
      s1: map max 37.8 nm  avg 25.8  (gates: exit PASS, clear FAIL)

    =================================================================
     S2  the SAME mirrors at the offset box -- FPA refit only
    =================================================================
      oi_solve fpa :  2 vars, 9 fields, start  37196.566 nm
        it  1:  24341.528 nm   (lam 3.7e-05)
        it  2:  24340.765 nm   (lam 9.9e-05)
        it  3:  24340.765 nm   (lam 2.6e-04)
      s2: map max 40824.4 nm  avg 17815.8  (gates: exit PASS, clear FAIL)

    =================================================================
     S3  re-solve symmetric aspheres/conics AT the offset field
    =================================================================
      S3 start candidates: carry 37196.6 nm, fresh seed 40247.3 nm -> carry
      oi_solve S3  : 16 vars, 9 fields, start  37196.566 nm
        it  1:  36682.801 nm   (lam 3.7e-05)
        it  2:  36400.947 nm   (lam 1.4e-06)
        it  3:  35805.323 nm   (lam 1.0e-07)
        it  4:  35591.201 nm   (lam 1.0e-07)
        it  5:  35292.592 nm   (lam 1.0e-07)
        it  6:  35001.959 nm   (lam 1.0e-07)
        it  7:  34002.207 nm   (lam 1.0e-07)
        it  8:  33425.758 nm   (lam 1.0e-07)
        it  9:  33261.885 nm   (lam 1.0e-07)
        it 10:  33217.751 nm   (lam 1.0e-07)
        it 11:  32685.632 nm   (lam 1.0e-07)
        it 12:  32452.454 nm   (lam 1.0e-07)
      s3: map max 41849.9 nm  avg 34306.9  (gates: exit PASS, clear FAIL)

    =================================================================
     S4  open tilt/decenter + radii under the constraint set
    =================================================================
      oi_solve S4  : 23 vars, 9 fields, start  32452.454 nm
        it  1:  27806.439 nm   (lam 3.7e-05)
        it  2:  26061.602 nm   (lam 1.4e-06)
        it  3:  25812.950 nm   (lam 1.0e-07)
        it  4:  25708.455 nm   (lam 1.0e-07)
        it  5:  25064.465 nm   (lam 1.0e-07)
        it  6:  24513.854 nm   (lam 1.0e-07)
        it  7:  24168.321 nm   (lam 1.0e-07)
        it  8:  23783.168 nm   (lam 1.0e-07)
        it  9:  23581.144 nm   (lam 1.0e-07)
        it 10:  23144.110 nm   (lam 1.0e-07)
        it 11:  22759.569 nm   (lam 1.0e-07)
        it 12:  22393.959 nm   (lam 1.0e-07)
      s4: map max 28907.4 nm  avg 23649.9  (gates: exit PASS, clear FAIL)

    =================================================================
     S5  open Zernike departures (BornWolf, lMon frozen)
    =================================================================
      oi_solve S5  : 83 vars, 9 fields, start  11810.889 nm
        it  1:  10305.818 nm   (lam 3.7e-05)
        it  2:   9789.750 nm   (lam 1.4e-06)
        it  3:   9431.743 nm   (lam 1.0e-07)
        it  4:   9237.638 nm   (lam 1.0e-07)
        it  5:   9086.290 nm   (lam 1.0e-07)
        it  6:   8892.939 nm   (lam 1.0e-07)
        it  7:   8734.527 nm   (lam 1.0e-07)
        it  8:   8506.421 nm   (lam 1.0e-07)
        it  9:   8258.537 nm   (lam 1.0e-07)
        it 10:   8051.323 nm   (lam 1.0e-07)
        it 11:   7835.824 nm   (lam 1.0e-07)
        it 12:   7604.716 nm   (lam 1.0e-07)
      s5: map max 12415.1 nm  avg 5993.6  (gates: exit PASS, clear FAIL)

      LADDER (dense 11x11 map max / avg / std, strict centroid nm)
       s1         37.8       25.8        9.2
       s2      40824.4    17815.8    10945.5
       s3      41849.9    34306.9     6248.0
       s4      28907.4    23649.9     4368.8
       s5      12415.1     5993.6     3421.9

    saved <work>/templates/10_telescopes/offset_imager/oi_run.mat + <work>/templates/10_telescopes/offset_imager/oi_REPORT.md
