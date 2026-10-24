    [1] CENTERED: axial 0.0037 | worst ring field 0.0575 waves -> DIFFRACTION-LIMITED

     field(arcmin x,y)     raw    -tilt   -focus  -astig   Z4(nm)  |Z5,6|(nm)
     ( +2.50, +0.00)     0.0558  0.0558  0.0372  0.0147    +72.2      85.0
     ( +1.77, +1.77)     0.0555  0.0555  0.0375  0.0146    +71.0      83.9
     ( +0.00, +2.50)     0.0546  0.0546  0.0365  0.0146    +70.6      83.2
     ( -1.77, +1.77)     0.0555  0.0555  0.0375  0.0146    +71.2      83.9
     ( -2.50, +0.00)     0.0558  0.0558  0.0372  0.0147    +72.4      85.1
     ( -1.77, -1.77)     0.0572  0.0572  0.0385  0.0147    +73.6      86.3
     ( -0.00, -2.50)     0.0570  0.0570  0.0378  0.0148    +74.1      86.7
     ( +1.77, -1.77)     0.0572  0.0572  0.0385  0.0148    +73.5      86.3
     ( +1.25, +0.00)     0.0477  0.0260  0.0128  0.0094    -39.3      21.5
     ( -1.25, +0.00)     0.0476  0.0259  0.0128  0.0094    -39.2      21.5
     worst:               0.0572  0.0572  0.0385  0.0148  (waves @ 1 um)
    [2] ladder: the centered family's residual survives the astig fit
        only at high order -- symmetric anastigmat behaviour.
    [3] clearance: 2/4 optics clear -- the CENTERED price is the
        central obscuration (M2 shadows the pupil).

    [4] SECTION (unobscured, decenter 0.97 D): worst ring field 0.1530 waves

    ================== A/B: 5'-dia circular field @ 1 um ==================
                           worst    -tilt   -focus  -astig    (waves)
     CENTERED (obscured)   0.0575   0.0572  0.0385  0.0148
     SECTION (unobscured)  0.1530   0.1180  0.0924  0.0354
     Read: the section's wavefront collapses once ASTIG is removed --
     the unobscured cut buys its clear pupil with induced field
     astigmatism.  j18's remaining edge is its 2.3 um yardstick
     (divide the waves by 2.3).
    ==========================================================================
    [5] saved: <work>/templates/10_telescopes/tma_centered/tma_centered.in
               + <work>/templates/10_telescopes/tma_centered/tma_centered.mat
        figures: tma_centered_wfe.png + tma_centered_layout.png
