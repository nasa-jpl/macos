    ====================================================================
     Best-form singlet | f = 100 mm | f/5 (EPD 20 mm) | N-BK7 n = 1.516800 at 587.562 nm | t = 4 mm
     object at infinity (magnification 0, position factor p = -1)
    ====================================================================

    [1] the plano singlet (R = 51.6800 mm, the benches' add_lens form with a sphere):
          curved face to the object (X = +1): RMS 1.4006 waves
          flat face to the object   (X = -1): RMS 5.4548 waves  (3.9x worse)

    [2] textbook minimum-SA bending: X = -2(n^2-1)p/(n+2) = 0.7397
          R1 = 59.2348  R2 = -395.8859 mm (R2/R1 = -6.683), scaled to EFL 100.000000 mm
          RMS 1.2941 waves at the paraxial focus (1.08x better than convex-plano)

    [3] CALIB on R1, R2 from the equi-convex start (RMS 2.0097 waves at X = 0):
          721 iterations (floor), image plane held at the start's paraxial focus
          solved R1 = 60.1205  R2 = -401.0366 mm, EFL 101.4660 mm -> X = 0.7393
          re-scaled to EFL 100: RMS 1.2941 waves at the paraxial focus
          shape factor: engine 0.7393 vs textbook 0.7397 (difference -0.0004)

    [4] bending scan, X = -1.50 .. 2.50 (41 lenses at EFL 100):
          minimum on the grid: X = 0.70 at the paraxial focus, X = 0.70 at best focus
          best focus: textbook 0.3269, convex-plano 0.3536, plano-convex 1.4074 waves

    DONE.  Report <work>/templates/35_lenses/singlet/out/singlet_report.txt
