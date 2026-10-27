    MACOS error budget: e5hex1
    metric=opd  base_units=  cbm=0.001  nw=51225  Nz=66
    piston removed (harvest opd_ref chief)
    control: 48 DOFs projected (post-WFC)

    group                    class               sigma        alpha      WFE(nm)
    PM Seg Tilt              angle               5e-08    < 1e-06     < 1e-06
    PM Seg Clocking          angle               1e-07    < 1e-06    < 1e-06
    PM Seg Decenter          length_si           1e-07    < 1e-06    < 1e-06
    PM Seg Piston            length_si       1.155e-07    < 1e-06    < 1e-06
    SM Tilt                  angle               5e-08    < 1e-06    < 1e-06
    SM Clocking              angle               1e-07    < 1e-06    < 1e-06
    SM Decenter              length_si           1e-07    < 1e-06    < 1e-06
    SM Despace               length_si       1.155e-07     < 1e-06     < 1e-06
    Lens Tilt                angle               5e-08     0.003225    0.0001613
    Lens Clocking            angle               1e-07    0.0006174    6.174e-05
    Lens Decenter            length_si           1e-07     0.001631    0.0001631
    Lens Despace             length_si       1.155e-07      0.01162     0.001342
    FPA Tilt                 angle               5e-08      0.00345    0.0001725
    FPA Clocking             angle               1e-07       0.1152      0.01152
    FPA Decenter             length_si           1e-07     0.001635    0.0001635
    FPA Despace              length_si       1.155e-07    0.0003083     3.56e-05
    WFS Sensing Error        wfs                 3e-06      0.03061      0.09183
    PM Seg Tilt Actuator     act             6.619e-10         1936        1.281
    PM Seg Clocking Actuator act             6.619e-10       0.7881    0.0005216
    PM Seg Decenter Actuator act                 1e-09        37.47      0.03747
    PM Seg Piston Actuator   act                 1e-09         1834        1.834
    SM Tilt Actuator         act             1.626e-09        769.2        1.251
    SM Clocking Actuator     act             1.626e-09       0.6792     0.001104
    SM Decenter Actuator     act                 1e-09        86.54      0.08654
    SM Despace Actuator      act                 1e-09        12.67      0.01267
    pre-WFC TOTAL                                                          236.2
    TOTAL (RSS)                                                            2.567

    Monte-Carlo (1000 draws) total: 2.573 nm  (cov 2.567 nm)

    summary                      preWFC    postWFC     preMET    postMET      Total
    Primary Mirror WFE            232.9   < 1e-06          0          0   < 1e-06
    Secondary Mirror WFE          39.45  < 1e-06          0          0  < 1e-06
    Lens WFE                     0.0859   0.001363          0          0   0.001363
    Detector WFE                  1.485    0.01152          0          0    0.01152
    WFC WFE                           0      2.567          0          0      2.567
    Everything Else               39.48      2.567          0          0      2.567
    TOTAL                         236.2      2.567          0          0      2.567
    Design WFE                      NaN        NaN        NaN        NaN          0
    Payload WFE                     NaN        NaN        NaN        NaN      2.567
    WFE Requirement                 NaN        NaN        NaN        NaN        NaN
    WFE Margin                      NaN        NaN        NaN        NaN        NaN

    method: Redding 2014, linearized ray-trace covariance budget.
