    ====================================================================
     3+1 coronagraph front end | D=6.61 m f/20 | 30" patch | j18 parent
    ====================================================================

    [1] TMA: axial 0.0037 waves; K=[-0.9950 -1.6350 -0.8514]; M3->focus 4.862 m
    [2] 3+1 coaxial, conics carried: 0.0034 waves on-axis at FP2
    [3] compact section (0.84 D) + axial refigure: 0.0143 waves
    [4] 30"-dia patch balance: worst 0.0337 waves
    [5] clearance: 5/5 clear -> UNOBSCURED
    exit-pupil surface @ elt 6  vertex [0.002 0.068 -4.361] m  diam 32.4 mm
      defocus -0.133  astig [+0.199 +0.000]  coma [+0.000 +0.001]  sph +0.000  (mm)
      sag 1.574 mm RMS (4.133 P-V)  fit resid 0.000 mm  chief-vs-FEX 1.0e-03 m
    [6] relayed pupil: dia 32.4 mm, defocus -0.133 mm, astig 0.199 mm   (bare TMA: +0.7 / 1.1)

    aoi_report  (limit 15 deg, 11 field(s))
      mirror   chief(deg)        min        max  SPREAD(deg)   ok
      M1            19.30       8.09      29.02        20.93   ** >limit
      M2            21.97       9.15      33.23        24.07   ** >limit
      M3             1.35       0.41       2.25         1.84   yes
      M4             1.11       0.49       1.73         1.24   yes
    [7] worst AOI spread 24.1 deg -- OUT of the 15-deg preference on
        M1/M2 (f/1.2 parent; see tma_3plus1_aoi_search.m).
    [8] saved: <work>/templates/10_telescopes/tma_3plus1/tma_3plus1.in
               + <work>/templates/10_telescopes/tma_3plus1/tma_3plus1.mat
        layout: tma_3plus1_layout.png
