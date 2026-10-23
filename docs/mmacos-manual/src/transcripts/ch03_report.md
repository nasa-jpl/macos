    ==== e2e stage 3: segmentation of M1 (parent s2_instrument.in) ====
    parent: D = 4 m, BaseUnits=m; gap 25 mm; per-variant detail in
    e2e_<variant>_seg_report.txt (run_segmentation runner reports)

    ---- variant "pie": rings=1, grid=Pie ----
        -> e2e_pie.in VERIFIED (11950/12128 pass, 0 clipped)

    ---- variant "hex2": rings=2, grid=Hex ----
        -> e2e_hex2.in VERIFIED (9698/9876 pass, 0 clipped)

    ==== stage-3 summary ====
    name    nseg   width_m   pass/src   clipped    rmsWFE d%
    pie        7     1.333  11950/12128        0       +3.44
    hex2      19       0.8   9698/9876         0      -25.48
    (rmsWFE deltas vs parent = pupil-sampling coverage: the segment tilings
     sample the pupil differently from the parent Circular grid)
    downstream variant (P.seg.variant) = "pie" -> stages 4-6 consume e2e_pie.in
