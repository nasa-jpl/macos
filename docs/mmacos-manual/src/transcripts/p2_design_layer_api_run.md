    Stage 3 — single-point evaluation (on-axis)
    ------------------------------------------------------------------
      on-axis RMS WFE (nominal) = < 1e-09 m   (RC: spherical + coma corrected)

     Stage 4 — same layout, different conics (RC vs classical Cassegrain)
    ------------------------------------------------------------------
      shared layout : R1 = 24.0000 m   R2 = 2.6263 m   M1-M2 sep = 10.8182 m
      RC          : K1 = -1.002185   K2 = -1.524259   (primary slightly hyperbolic)
      Cassegrain  : K1 = -1.000000   K2 = -1.493827   (parabolic primary)
      => RC trades a slightly-hyperbolic primary to also null coma (aplanat),
         widening the usable field; the classical Cass is simpler but coma-limited.

     Stage 5 — align a perturbed secondary (M2 despace)
    ------------------------------------------------------------------
      perturbed (+2.000 mm) : RMS WFE = 1.752e-05 m
      recovered             : RMS WFE = < 1e-09 m   at despace =  0.0000 mm
      (optimizer returned the secondary to within 0.0 um of nominal)

     Stage 6 — add an accessible exit pupil (deliverable)
    ------------------------------------------------------------------
      exit pupil @ z = -11.9892 m,  reference radius = 12.9892 m
      elements now: image-return #3, exit-pupil #4, detector #5
