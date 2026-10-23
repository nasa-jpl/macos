    ==== run_simulator: <thread>/e2e_pie.in ====
    [0] 9 bodies (7 segments), T = 50 frames, dt = 10 s
        states: x 54x50, z 56x50, g 42x50;  sim Rx <thread>/s4_grid.in
    [1] baseline on <thread>/s4_grid.in: 11950 rays in pupil; psf @ elt 17, peak 6.53e+03;
        COMPOSE band [4.85e-07 5e-07 5.15e-07] WaveUnits, 128 px @ 2.18 um; (grid Rx: no metrology) engine l model-only
    [2] channels in motion: 48 rigid, 21 zernike, 42 grid
    [3] pass 1 (uncorrected): 10 20 30 40 50 done
        initial wavefront control @ frame 2 (t=20 s): uncorrected w rms 100.8 um -> predicted residual 24.33 nm
        (u on 8 bodies, ridge 0.0003, 21/48 directions with sv > lam, |u| max 0.0001354)
        wfc iteration 2: engine residual 2159 nm -> update |du| max 6.671e-06
        wfc iteration 3: engine residual 19.78 nm -> update |du| max 6.011e-06
    [3b] metrology loop ON (RBCS weighted-LS estimator + pose-error controller):
         R_meas = (H'N^-1 H + Rx^-1)^-1 H'N^-1, kp 0.5, rho 0; sensor sigma [1e-12 1e-09] m;
         state prior sigma 2.1e-11..2.1e-08 m (auto from drift stats); loop senses figure: 0; meas bias: 0
         WF-maintenance recontrol at t = 400 s (reset ridge 1e-05)

    [5] final: rms WFE 1.008e+05 nm uncorrected vs 30.15 nm corrected;
        Strehl 0.0000 uncorrected vs 0.8695 corrected (0.9830 broadband); worst w_rel 2.14e-02
        m bars are the LINEAR measurement model (grid Rx declares no metrology;
        the engine METcalc/Hx hold met/sensor points rigid for figure states)

    artifacts: s7B_sim.mat, frames in s7B_sim_frames/, s7B_sim.gif, this report
