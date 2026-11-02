    pupil beam DIAMETERS (mm) -- the source Aperture is the full cone angle, so R_DM*FILL lands as the diameter: DM 21.38  apod 16.03  Lyot 8.02  backend 8.02
    ctb: 17 elements, chief path 13494.603 mm
       1  OAP1         Reflector   s= 2519.136  vpt=[ 2462.019   434.120     0.000]
       2  DM1          Reflector   s= 3119.136  vpt=[ 1928.251  -104.189     0.000]
       3  DM2          Reflector   s= 3619.136  vpt=[ 2428.251  -104.189     0.000]
       4  OAP2         Reflector   s= 4319.136  vpt=[ 1773.704  -490.382     0.000]
       5  Focus23      Reference   s= 5854.801  vpt=[ 3274.551  -225.743     0.000]
       6  OAP3         Reflector   s= 7006.550  vpt=[ 4400.186   -27.263     0.000]
       7  Apodizer     Reference   s= 7356.550  vpt=[ 4081.617  -286.519     0.000]
       8  OAP4         Reflector   s= 7856.550  vpt=[ 3620.056  -607.769     0.000]
       9  FPM          Reference   s= 9216.883  vpt=[ 4949.546  -373.344     0.000]
      10  OAP5         Reflector   s= 9897.049  vpt=[ 5614.291  -256.131     0.000]
      11  Lyot         Reference   s=10247.049  vpt=[ 5285.030  -434.120     0.000]
      12  OAP6         Reflector   s=10697.049  vpt=[ 4856.374  -622.529     0.000]
      13  FieldStop    Reference   s=11336.910  vpt=[ 5481.727  -512.262     0.000]
      14  OAP7         Reflector   s=11976.770  vpt=[ 6107.080  -401.996     0.000]
      15  Backend      Reference   s=12276.770  vpt=[ 5826.145  -564.357     0.000]
      16  OAP8         Reflector   s=12726.770  vpt=[ 5400.391  -774.818     0.000]
      17  FPA          FocalPlane  s=13494.603  vpt=[ 6150.814  -642.498     0.000]
      ok/elt [3002  3002  3002  3002  3002  3002  3002  3002  3002  3002  3002  3002  3002  3002  3002  3002  3002]  chief max < 1e-09 mm
      RMS WFE 0.001423 waves | DM1->FPA zone spot med 0.02713 um worst 0.05121 um
      pol @ FPA: ret mean 0.02109 (var 4.90e-05), diatten mean 0.00362 (var 8.41e-06)

    === POLARIZATION at FPA (Al-coated mirrors) ===
      retardance : mean 0.02109 rad   var_rms 4.896e-05 rad
      diattenuat.: mean 0.00362       var_rms 8.414e-06
      NOTE: fold polarization compensation is non-trivial -- in-plane fold
      alternation does NOT cancel retardance, and naive crossed-plane folding
      cuts diattenuation but not retardance-variation and costs WFE.  Proper
      compensation pairs equal-AOI folds with exactly-crossed planes AND
      balances the induced astigmatism (flagged follow-up).

    === PIL DESIGN STEP (pupil-imaging lens) ===
      PIL f=150mm: pupil image diam 3.29 mm, camera +0.0 mm vs f=150mm  (<work>/templates/30_instruments/bench_ctb/ctb_pil150.in)
      PIL f= 75mm: pupil image diam 1.23 mm, camera -100.7 mm vs f=150mm  (<work>/templates/30_instruments/bench_ctb/ctb_pil75.in)
      (size ratio 2.68x and camera move ~101 mm mirror Marx's 150mm/75mm 1000/500-px, +106mm trade)

    DONE.  Optimized source-imaging Rx: <work>/templates/30_instruments/bench_ctb/ctb_planar_opt.in
           PIL pupil-imaging Rx: <work>/templates/30_instruments/bench_ctb/ctb_pil150.in, <work>/templates/30_instruments/bench_ctb/ctb_pil75.in
