    ==================== DESIGN REPORT ====================
     family TMA | D = 4.000 m | lambda = 0.5 um | 11 elements
       M1         R =   14.0000 m  K =  -0.98972  Vpt=[  0.000   0.000   0.000]
       M2         R =    1.6508 m  K =  -1.51153  Vpt=[  0.000   0.000  -6.278]
       FM         Reflector  (flat)        Vpt=[  0.000   0.000   0.300]
       M3         R =    3.2625 m  K =  -0.97231  Vpt=[  2.100   0.000   0.300]
       M4         R = 1000000.0000 m  K =  +0.00000  Vpt=[ -2.079   0.000   0.300]
       M5         R =    2.0000 m  K =  +0.00000  Vpt=[ -0.569  -0.266   0.300]
       M6         R =    1.0000 m  K =  +0.00000  Vpt=[ -1.494   0.153   0.300]
       M7         R =    2.0000 m  K =  +0.00000  Vpt=[ -0.481   0.234   0.300]
       FP_return  Return     (flat)        Vpt=[ -2.436   0.553   0.300]
       ExitPupil  R =    1.8489 m  K =  +0.00000  Vpt=[ -0.606   0.291   0.299]
       FP         FocalPlane (flat)        Vpt=[ -2.436   0.553   0.300]
     -- first order --
       plate-scale EFL at the field point 73.88 m (system EFL 72.71 m) | f/# at M1 1.75 | f/# at FP 18.47
       plate scale 358.2 um/arcsec | lambda/D 25.78 mas
     -- field performance (RMS WFE about the field center) --
       bias +0.3' |    field       raw     -tilt    Strehl
              center     0.008     0.007     0.998
          0.50' ring     0.009     0.007     0.998
          1.00' ring     0.018     0.010     0.996
          1.50' ring     0.043     0.014     0.994
          2.00' ring     0.086     0.022     0.982
       (Strehl = far-field peak / unaberrated-aperture peak, exact coherent sum over the
        EP-referenced OPD, tilt removed; worst ring point.  PSF products: INT at the FP.)
     -- focal plane --
       FP tilt 0.309 deg wrt the arriving chief ray
       aligned from 25 field foci: defocus removed +0.046 mm; field-curvature
       sag -5.6 to +5.6 um (rms 3.3 um); best-focus blur 2.64e-07 m (center)
     -- exit pupil (add_pupil / FEX) --
       1.849 m from the image, at [ -0.606   0.291   0.299] (elt 10)
     -- packaging / clearance --
       shroud 1.25 x D | train length 6.69 m
       clearance: body-in-beam at M2, FM (FP vs its own pupil-retrace legs excluded)
       M1         AOI   0.1 deg (spread 8.0 deg across the beam)
       M2         AOI   0.5 deg (spread 9.3 deg across the beam)
       FM         AOI  45.0 deg (spread 4.1 deg across the beam)
       M3         AOI   1.0 deg (spread 0.5 deg across the beam)
       M4         AOI   2.3 deg (spread 3.1 deg across the beam)
       M5         AOI   5.4 deg (spread 0.1 deg across the beam)
       M6         AOI  11.0 deg (spread 0.1 deg across the beam)
       M7         AOI   6.1 deg (spread 0.0 deg across the beam)
    =======================================================
     -- stage-2 addendum: the imaging instrument --
       relay: M4 corrector 0.45 m past the telescope focus / M5 collimator / M6 corrector / M7 camera
       relay R = [1e+06 2 1 2] m, K = [0 0 0 0] (M5/M7 solved at the bias point)
       joint solve: CALIB M2/M3/M5/M6 (11 pts) + SVD M4 + polish M2..M6 + M1 common mode (5x5 grid)
       field-zone lMon (m): [2.1 0.221 0.166 0.0767 0.167 0.0564 0.0629]  (s1: [2.1 0.219 0.137])
       baseline +-2': 0.101 -tilt -> after joint solve 0.0561 raw / 0.0232 -tilt waves (DL)
       bias point 0.0178 -tilt waves | max|coef| [0.00064 0.00084 0.00061 5.3e-08 0.0004 0 0.00065] m
       field center (auto): s1 bias 2' start -0.7' -> adopted 0.261538' (chief at the <0.02-waves-region centroid; residual [+0.00 -0.03]'); center scan [dy' raw -tilt]:
         +0.70/0.063/0.038  +0.35/0.036/0.029  +0.00/0.056/0.023  -0.35/0.096/0.031  -0.70/0.146/0.056
       distortion (M4 reflective field-corrector solve): 207.6 -> 207.6 um rms at the detector
    =======================================================
