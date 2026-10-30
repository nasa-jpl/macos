    [1] chain (eccentric-pupil sections, decenter 5.245 m):
        M1   Reflector  Vpt=[   0.000    0.000    0.000]
        M2   Reflector  Vpt=[   0.000    0.000  -13.375]
        M3   Reflector  Vpt=[   0.000    0.000    3.963]
        FP   FocalPlane Vpt=[   < 0.0001    0.000   -5.614]

    [2] correction ladder (waves RMS @ 500 nm, field center):
        [2a] as found:           raw   0.604 | -tilt   0.306 | -focus   0.214 | -astig   0.098
        [2b] off-axis refigure:  raw   0.001 | -tilt   0.001 | -focus   0.001 | -astig   0.001
        [2c] field balance:      raw   0.004 | -tilt   0.004 | -focus   0.003 | -astig   0.003
        [2d] + freeform M1:      raw   0.006 | -tilt   0.003 | -focus   0.003 | -astig   0.002
        field center: 0.006 waves @ 500 nm -> DIFFRACTION-LIMITED
    [2e] true FP from 5x5 field foci (+center): tilt 1.208 deg, defocus removed -1.439 mm,
         field-curvature sag -0.5 to +0.5 um

    [3] FOV ladder (worst on ring, waves RMS @ 500 nm; blur on -tilt):
             ring       raw     -tilt
            0.10'     0.009     0.004  <- DL
            0.25'     0.010     0.005  <- DL
            0.50'     0.018     0.008  <- DL
            1.00'     0.047     0.022  <- DL
    check_clipping  (family=tma, 4 elements)
      name       kind            ap_r    foot_r    margin  clearnce obstruct  status
      M1         Reflector      3.633     3.302    0.3309     1.894        0  OK
      M2         Reflector      3.303    0.6354     2.667    0.2972        0  OK
      M3         Reflector      3.303    0.2383     3.064        --        0  OK
      FP         FocalPlane     1.982 < 0.0001     1.981    0.2804        0  OK
      => layout is CLEAR (no body-in-beam, beams fit apertures)

    [4] clearance: UNOBSCURED -- every body clear (minimal decenter = M2 as close to the
        source->M1 beam as the 5%-of-D margin allows)

    [5] saved: <work>/templates/10_telescopes/tma_unobscured/tma_unobscured.in
               + <work>/templates/10_telescopes/tma_unobscured/tma_unobscured.mat
        figures: tma_unobscured_layout.png + tma_unobscured_fpmap.png

    [6] design report:
