    [1] layout: R=[3.000 0.349 0.580] m  t=[1.347 1.947] m  EFL=8.000 m (f/8.00)
        intermediate focus z=-0.125 m (between M1 at 0 and M2 at -1.347); M3 z=0.600 m behind M1

    [2] bias sweep (least bias that clears the FP wins -- WFE grows ~bias^2):
          bias(arcmin)   RMS WFE(waves)   FP out of beam?
               1.0           0.0059       YES
               2.0           0.0098       YES
               3.0           0.0164       YES
               4.0           0.0296       YES
               6.0           0.0716       YES

    [3] RECOMMENDED bias = 1' (least bias that clears the FP and is diffraction-limited)
    [4] design @ bias 1': 0 -> 0.0059 waves -> DIFFRACTION-LIMITED  (K=[-0.9828 -1.4158 -0.6067])
        clearance (obscured baseline -- M2 central obscuration expected):
          M1         : OBSTRUCTS beam
          M2         : OBSTRUCTS beam
          M3         : clear
          FP         : clear
        exit pupil after M3 (FEX, assessed): z=0.264 m (0.26*D from M1), radius=0.147 m

    [5] saved: <work>/templates/10_telescopes/tma_onaxis/tma_onaxis.in
               + <work>/templates/10_telescopes/tma_onaxis/tma_onaxis.mat
        layout: tma_onaxis_layout.png

    --------------------------------------------------------------------
     D=1.00 m | f/8 | EFL 8.00 m | Cass focus between M1 & M2 | bias 1' | WFE 0.0059 wv
     obscured baseline (M2 central obscuration); FP biased out of beam; EP after M3 assessed.
     Unobscured eccentric-pupil cousin: templates/10_telescopes/tma_offaxis.
    ====================================================================
