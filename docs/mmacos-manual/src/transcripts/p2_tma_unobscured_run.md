    [1] chain (eccentric-pupil sections, decenter 10.989 m):
        M1   Reflector  Vpt=[   0.000    0.000    0.000]
        M2   Reflector  Vpt=[   0.000    0.000  -13.375]
        M3   Reflector  Vpt=[   0.000    0.000    3.963]
        FP   FocalPlane Vpt=[   < 0.0001    0.000   -5.614]

    [2] correction ladder (waves RMS @ 500 nm, field center):
        [2a] as found:           raw  12.278 | -tilt   3.443 | -focus   2.277 | -astig   0.487
        [2b] off-axis refigure:  raw   0.010 | -tilt   0.010 | -focus   0.010 | -astig   0.010
        [2c] field balance:      raw   0.005 | -tilt   0.005 | -focus   0.005 | -astig   0.004
        [2d] + freeform M1:      raw   0.011 | -tilt   0.005 | -focus   0.004 | -astig   0.003
        field center: 0.011 waves @ 500 nm -> DIFFRACTION-LIMITED
    [2e] true FP from 5x5 field foci (+center): tilt 2.411 deg, defocus removed -7.531 mm,
         field-curvature sag -0.8 to +0.8 um

    [3] FOV ladder (worst on ring, waves RMS @ 500 nm; blur on -tilt):
             ring       raw     -tilt
            0.10'     0.030     0.005  <- DL
            0.25'     0.033     0.007  <- DL
            0.50'     0.046     0.013  <- DL
            1.00'     0.097     0.032  <- DL
    check_clipping  (family=tma, 4 elements)
      name       kind            ap_r    foot_r    margin  clearnce obstruct  status
      M1         Reflector      3.633     3.302    0.3309     7.323        0  OK
      M2         Reflector      3.303     0.676     2.626    0.3777        0  OK
      M3         Reflector      3.303    0.2351     3.067        --        0  OK
      FP         FocalPlane     1.982 < 0.0001     1.981    0.8575        0  OK
      => layout is CLEAR (no body-in-beam, beams fit apertures)

    [4] clearance: UNOBSCURED -- every body clear (minimal decenter = M2 as close to the
        source->M1 beam as the 5%-of-D margin allows)

    [5] saved: <work>/templates/10_telescopes/tma_unobscured/tma_unobscured.in
               + <work>/templates/10_telescopes/tma_unobscured/tma_unobscured.mat
        figures: tma_unobscured_layout.png + tma_unobscured_fpmap.png

    [6] design report:
