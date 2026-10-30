    ==== run_simulator: <thread>/e2e_pie.in ====
    [0] 9 bodies (7 segments), T = 50 frames, dt = 10 s
        states: x 54x50;  sim Rx <thread>/e2e_pie_metopt.in
        OPD reference chief (the harvest's), piston removed on both sides
    [1] baseline on <thread>/e2e_pie_metopt.in: 11950 rays in pupil; psf @ elt 17, peak 6.53e+03;
        COMPOSE band [4.85e-07 5e-07 5.15e-07] WaveUnits, 128 px @ 2.18 um; engine l cross-checked (corrected leg)
    [2] channels in motion: 48 rigid, 0 zernike, 0 grid
    [3] pass 1 (uncorrected): 10 20 30 40 50 done
        initial wavefront control @ frame 2 (t=20 s): uncorrected w rms 100.8 um -> predicted residual 6.556 nm
        (u on 8 bodies, ridge 0.0003, 21/48 directions with sv > lam, |u| max 0.0001353)
        wfc iteration 2: engine residual 2159 nm -> update |du| max 2.158e-06
        wfc iteration 3: engine residual 1.155 nm -> update |du| max 1.313e-06
    [3b] metrology loop ON (RBCS weighted-LS estimator + pose-error controller):
         R_meas = (H'N^-1 H + Rx^-1)^-1 H'N^-1, kp 0.5, rho 0; sensor sigma [1e-12 1e-09] m;
         state prior sigma 2.1e-11..2.1e-08 m (auto from drift stats); loop senses figure: 1; meas bias: 1
         WF-maintenance recontrol at t = 400 s (reset ridge 0.0003)

    [5] final: rms WFE 1.008e+05 nm uncorrected vs 11.81 nm corrected;
        Strehl 0.0000 uncorrected vs 0.9782 corrected (1.0140 broadband); worst w_rel 2.14e-02
        engine-l cross-check (corrected leg): worst l_rel 3.47e-05

    artifacts: s7A_sim.mat, frames in s7A_sim_frames/, s7A_sim.gif, this report
