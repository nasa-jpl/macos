    MACOS error budget: e5hex1
    metric=opd  base_units=  cbm=0.001  nw=51225  Nz=66
    control: 48 DOFs projected (post-WFC)

    group                    class               sigma        alpha      WFE(nm)
    PM Seg Tilt              angle               5e-08    < 1e-06    < 1e-06
    PM Seg Clocking          angle               1e-07    < 1e-06    < 1e-06
    PM Seg Decenter          length_si           1e-07     < 1e-06     < 1e-06
    PM Seg Piston            length_si       1.155e-07    < 1e-06    < 1e-06
    SM Tilt                  angle               5e-08    < 1e-06     < 1e-06
    SM Clocking              angle               1e-07    < 1e-06    < 1e-06
    SM Decenter              length_si           1e-07    < 1e-06    < 1e-06
    SM Despace               length_si       1.155e-07    < 1e-06    < 1e-06
    Lens Tilt                angle               5e-08      0.00402     0.000201
    Lens Clocking            angle               1e-07    0.0006775    6.775e-05
    Lens Decenter            length_si           1e-07     0.001926    0.0001926
    Lens Despace             length_si       1.155e-07      0.01309     0.001512
    FPA Tilt                 angle               5e-08     0.004771    0.0002385
    FPA Clocking             angle               1e-07       0.1229      0.01229
    FPA Decenter             length_si           1e-07      0.00204     0.000204
    FPA Despace              length_si       1.155e-07    0.0003093    3.571e-05
    WFS Sensing Error        wfs                 3e-06      0.03061      0.09183
    PM Seg Tilt Actuator     act             6.619e-10         1936        1.281
    PM Seg Clocking Actuator act             6.619e-10       0.7881    0.0005216
    PM Seg Decenter Actuator act                 1e-09        37.47      0.03747
    PM Seg Piston Actuator   act                 1e-09         2595        2.595
    SM Tilt Actuator         act             1.626e-09        769.2        1.251
    SM Clocking Actuator     act             1.626e-09       0.6792     0.001104
    SM Decenter Actuator     act                 1e-09        86.55      0.08655
    SM Despace Actuator      act                 1e-09        13.75      0.01375
    pre-WFC TOTAL                                                          317.4
    TOTAL (RSS)                                                            3.156

    Monte-Carlo (1000 draws) total: 3.179 nm  (cov 3.156 nm)

    summary                      preWFC    postWFC     preMET    postMET      Total
    Primary Mirror WFE            314.9   < 1e-06          0          0   < 1e-06
    Secondary Mirror WFE          39.46  < 1e-06          0          0  < 1e-06
    Lens WFE                    0.08615   0.001539          0          0   0.001539
    Detector WFE                  1.485     0.0123          0          0     0.0123
    WFC WFE                           0      3.156          0          0      3.156
    Everything Else               39.48      3.156          0          0      3.156
    TOTAL                         317.4      3.156          0          0      3.156
    Design WFE                      NaN        NaN        NaN        NaN          0
    Payload WFE                     NaN        NaN        NaN        NaN      3.156
    WFE Requirement                 NaN        NaN        NaN        NaN        NaN
    WFE Margin                      NaN        NaN        NaN        NaN        NaN

    method: Redding 2014, linearized ray-trace covariance budget.
