    ==== run_met: <thread>/e2e_pie.in ====
    [0] e2e_pie: 7 segments (Pie tiling), 17 elements; 11950/12128 rays pass

    [1] as-built truss: 48 beams (6 launchers/segment + aft ring on elt 9), 6 fiducials at r=0.1963 on elt 8
        trace parity: 11950/12128 rays pass (bare 11950/12128), rmsWFE 4.197e-09 (bare 4.197e-09)
        aft-beam hole clearance (effective hole r=0.232 (design override)):
          crossing radii [0.1798 0.226 0.1807 0.182 0.226 0.1763] -> min margin 0.005971 (CLEAR)

    [2] dedx 72x54 (SI);  dldx 48x54 engine-FD vs analytic: < 1e-05 max rel  PASS

    [3] dwdx: 54 columns for bodies [1  2  3  4  5  6  7  8  9] from the sensitivities .mat
        post-control wavefront residual (fraction of prior; nm at the
        configured scenario: prior 1e-10 rad / 1e-10 m, noise 1e-11 / 1e-12 m):
          prior (no sensing)       1.0000   (0.2173 nm)
          edge sensors only        0.4348   (0.09447 nm)
          MET only                 0.1357   (0.02949 nm)
          edge + MET               0.0944   (0.0205 nm)
        as-built WEM (per unit gauge noise, suite ratios held;
        prior split: < 1e-12 nm tilt / 0.2173 nm non-tilt):
          config        WEM_full  WEM_tilt WEM_ntilt  floor_ntilt
          MET only        15.322     0.000    15.322       0.1159
          edge+MET         3.943     0.000     3.943       0.0926
        observability: H 120x54, sv(H*sqrt(X)) 7.185e-10 .. < 1e-15 (cond 3.24e+06)
        products: dxde 54x72, dxdl 54x48, dwde 59750x72, dwdl 59750x48

    [4] met_layout_opt (one pattern per shape class; MET-only, NON-TILT merit):
        class 1 (1 members, cluster): angs [170 190 50 70 -70 -50] deg
        class 2 (6 members, cluster): angs [175 185 55 65 -65 -55] deg
        engine-FD validation (non-tilt merit): 0.0204 nm (analytic 0.0204, 0.00% off)

    [4b] configuration comparison (WEM per unit gauge noise):
          config         | MET_full MET_tilt   MET_nt | e+M_full e+M_tilt   e+M_nt |    flr_nt  min_sep
          as-built       |   15.322    0.000   15.322 |    3.943    0.000    3.943 |    0.0926 
          optimized      |    3.877    0.000    3.877 |    1.782    0.000    1.782 |    0.0922 
          corner_pairs   |    4.102    0.000    4.102 |    1.799    0.000    1.799 |    0.0922  0.01502

    [5] Monte-Carlo (200 draws): 0.020 nm rms  (analytic 0.020) [2.0%]

    [6] figures: <thread>/e2e_pie_met_layout.png, <thread>/e2e_pie_metopt_layout.png, <thread>/e2e_pie_met_view_rx.png, <thread>/e2e_pie_met_metric.png
    preset saved: <thread>/e2e_pie_met_preset.mat (reusable as-built for future Pie builds)

    artifacts: e2e_pie_met.in, e2e_pie_metopt.in, e2e_pie_met.mat, this report + figures
    next stage: run_compare (linear model vs mmacos on x/z/grid states).
