    ==== run_sensitivities: <thread>/e2e_pie.in ====
    field set: center + 4 corners at +-0.0005818 rad
    OPD conventions: orient=raw sign=opl opd_ref=chief
    segments: 7

    [dwdx] rigid-body 6-DOF per element, LOCAL triads...
        dwdxall 59750 x 90 over 5 fields

    [dwdz] segment-LOCAL MonZernike modes [4 5 6 7 8 9 10 11]...
        dwdzall 59750 x 56

    [dwdgrid] grid-augmenting in the clocked Mon frames...
        s4_grid.in: 17 elts, 11950/12128 rays, gdx 0.01569
        dwdgall 59750 x 42 (influence basis saved in og.sgb)

    [conditioning] finite rows only:
        channel          size   rank     sv_max    sv_min+      cond+
        dwdx     59750x90          <90  4.984e+02  < 5e-05  > 1e+06
        dwdx     segment-only: 59750x42  cond+ > 1e+06
        dwdz     59750x56         56  4.903e+02  3.385e+01  1.448e+01
        dwdgrid  59750x42         42  8.570e+02  1.566e+02  5.471e+00

        dwdx per-segment column norms (rms native-unit per unit DOF):
          seg     rx        ry        rz        tx        ty        tz
            1   2.39e-01  2.40e-01  < 5e-05  1.71e-02  1.70e-02  1.89e+00 
            2   3.32e-01  2.71e-01  1.01e-04  1.90e-02  2.35e-02  7.66e-01 
            3   3.32e-01  2.71e-01  1.01e-04  1.90e-02  2.35e-02  7.66e-01 
            4   3.33e-01  2.70e-01  1.02e-04  1.90e-02  2.36e-02  7.66e-01 
            5   3.32e-01  2.71e-01  1.03e-04  1.90e-02  2.35e-02  7.66e-01 
            6   3.32e-01  2.71e-01  1.03e-04  1.90e-02  2.35e-02  7.66e-01 
            7   3.33e-01  2.70e-01  1.02e-04  1.90e-02  2.36e-02  7.66e-01 

    figures: s4_opdall/svspec + 90 dW pages in s4_pages/ (>= 3.5 in per OPD map, 1.2 in per field tile) -- index: s4_pages_index.txt

    Done: <thread>/s4_sens.mat + s4_sens_report.txt
