    ==================== DESIGN REPORT ====================
     family TMA | D = 4.000 m | lambda = 0.5 um | 7 elements
       M1         R =   14.0000 m  K =  -0.98972  Vpt=[  0.000   0.000   0.000]
       M2         R =    1.6508 m  K =  -1.51153  Vpt=[  0.000   0.000  -6.278]
       FM         Reflector  (flat)        Vpt=[  0.000   0.000   0.300]
       M3         R =    3.2625 m  K =  -0.97231  Vpt=[  2.100   0.000   0.300]
       FP_return  Return     (flat)        Vpt=[ -1.628   0.091   0.300]
       ExitPupil  R =    1.7569 m  K =  +0.00000  Vpt=[  0.127   0.014   0.300]
       FP         FocalPlane (flat)        Vpt=[ -1.628   0.091   0.300]
     -- first order --
       plate-scale EFL at the field point 72.32 m (system EFL 72.14 m) | f/# at M1 1.75 | f/# at FP 18.08
       plate scale 350.6 um/arcsec | lambda/D 25.78 mas
     -- field performance (RMS WFE about the field center) --
       bias +2.0' |    field       raw     -tilt    Strehl
              center     0.005     0.005     0.999
          0.25' ring     0.005     0.004     0.999
          0.50' ring     0.006     0.005     0.999
          1.00' ring     0.015     0.014     0.992
          2.00' ring     0.069     0.065     0.846
       (Strehl = far-field peak / unaberrated-aperture peak, exact coherent sum over the
        EP-referenced OPD, tilt removed; worst ring point.  PSF products: INT at the FP.)
     -- focal plane --
       FP tilt 0.002 deg wrt the arriving chief ray
       aligned from 25 field foci: defocus removed +0.025 mm; field-curvature
       sag -1.8 to +1.8 um (rms 1.1 um); best-focus blur 1.40e-07 m (center)
     -- exit pupil (add_pupil / FEX) --
       1.757 m from the image, at [  0.127   0.014   0.300] (elt 6)
     -- packaging / clearance --
       shroud 1.05 x D | train length 6.69 m
       clearance: body-in-beam at M2 (FP vs its own pupil-retrace legs excluded)
       M1         AOI   0.1 deg (spread 8.0 deg across the beam)
       M2         AOI   0.3 deg (spread 9.2 deg across the beam)
       FM         AOI  45.0 deg (spread 4.1 deg across the beam)
       M3         AOI   1.0 deg (spread 0.5 deg across the beam)
    =======================================================
     -- stage-1 addendum: requested vs achieved --
       requested: D 4.00 m | primary f/1.75 | system f/18 | m2 8
       first-order: EFL 72.00 m (f/18.00); achieved live EFL 72.32 m (f/18.08)
       fold FM at z=+0.300 m -> bench behind M1: M3 [2.10 0.00 0.30], FP [-1.63 0.09 0.30]
       bias sweep [1 2 3 4 6]' -> recommended 2' (least whose folded design clears)
       conic solve (bias point): K = [-0.98972 -1.51153 -0.97231], 0.0470 waves
       freeform:  lMon = [2.097 0.219 0.137] m (field-zone), max|coef| = [6.4e-04 9.1e-04 2.4e-03] m
       bias point 0.0003 waves (DL); worst-field @ +-1': 0.0138 raw / 0.0137 -tilt waves @ 500 nm
       the +-1' spread is the telescope's field differential -- stage 2's instrument widens toward +-2'
       M1 hole r 0.232 m | intermediate focus f/14.0 at z=-0.500 m
    =======================================================
