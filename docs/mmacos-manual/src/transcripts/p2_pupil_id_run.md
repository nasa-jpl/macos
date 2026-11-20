    [2] EP at element 1 via anchor='rim' (EP as the beam rim -- classic telescope)
    [3] cone field grid: 9 cones (3x3, +/-2 arcmin); blur.rms=7.659e-06 m  blur_ratio=0.010
    [4] XP sphere -> Rx: vertex=[-0.00000 -0.00000 0.26358] m (fit_chief), rad=-0.146931 m (follows vertex)
        pupil-quality: convergence curvature R=0.07644 m, departure-from-sphere RMS=1.853e-07 m
        (bundle diag vs FEX: 5.16e-06 m; written vs FEX along chief: 0.00e+00 m)
    [5] XPS cross-check (pupil_map anchor='stop' vs engine pupil_quality):
        defocus 2.738e-04 vs 2.737e-04 (0.02%)  astig -3.071e-04 vs -3.071e-04 (0.00%)  [rescale s=1.0511]
    [6] EP->XP imaging: 25 zones, median spot=1.737e-03 m, worst=3.002e-03 m
    [7] pupil walk over -2..+2 arcmin: XP z 0.263557..0.263582 m, rad -0.147057..-0.146915 m
    [8] revised Rx written: <work>/templates/10_telescopes/pupil_id/results/tma_onaxis/tma_onaxis_xp.in (fit_chief vertex + vertex-following radius)
        wrote <work>/templates/10_telescopes/pupil_id/results/tma_onaxis/pupil_id_cloud.png
        wrote <work>/templates/10_telescopes/pupil_id/results/tma_onaxis/pupil_id_zernikes.png
        wrote <work>/templates/10_telescopes/pupil_id/results/tma_onaxis/pupil_id_walk.png

    =====================================================================
     PUPIL ID -- tma_onaxis, beyond FEX
       FEX (1 chief ray):  XP z=0.26358 m   propagation rad=-0.146931 m
       cone-convergence :  blur=7.66e-06 m   departure-from-sphere=1.85e-07 m
       pupil aberration :  defocus=2.66e-04  astig=-2.19e-08  spherical=8.72e-08 (rim anchor)
       XPS cross-check  :  defocus 0.02% / astig 0.00% vs engine pupil_quality
       EP->XP imaging   :  median zone spot=1.74e-03 m (25 zones)
       revised Rx       :  <work>/templates/10_telescopes/pupil_id/results/tma_onaxis/tma_onaxis_xp.in
    =====================================================================

    =====================================================================
     PUPIL-ID RUNNER SUMMARY (1 case)
    =====================================================================
     case        | FEX rad   | conv rad   | dep RMS    | XPS df%  | XPS as% | revised Rx
     --------------------------------------------------------------------------------------------
     tma_onaxis  |  -0.14693 |    0.07644 | 1.85e-07 |    0.02 |    0.00 | tma_onaxis_xp.in
     --------------------------------------------------------------------------------------------
     (FEX rad = XP->detector propagation radius, written to the Rx; conv rad =
      pupil-imaging convergence-surface curvature, a QUALITY diagnostic, NOT the Rx
      radius. XPS df/as% = pupil_map anchor=stop vs engine pupil_quality (% rel).)
