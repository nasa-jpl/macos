    ==== run_compare: <thread>/e2e_pie.in ====
    [0] 9 bodies (7 segments), channels [x z grid]; met Rx <thread>/e2e_pie_metopt.in
        jac: center field 11950 x 54 x-cols; dedx 72x54; dldx 48x54
        baseline: 11950 OPD pixels (jac mask mismatch 0 px), 48 met beams, 72 edge rows

    [3] worst engine-vs-linear agreement per channel:
        x     ( 54 pokes): w 9.45e-01 (Seg 1 Tz +100 nm), l < 1e-05, e < 1e-06
        z     ( 56 pokes): w 9.90e-01 (Seg 1 Z5 +100 nm), l     null, e     null
        grid  ( 42 pokes): w 9.87e-01 (Seg 1 G2 +100 nm), l     null, e     null
        (4 poke(s) below the 0.001 nm response floor -- no physical
         response, engine delta is FD noise; rel = null)
        (rel = engine-vs-linear residual over the engine response;
         at these poke sizes this measures the linearization error.
         z/grid l/e: the LINEAR side shows the dmdz/dmdgrid figure-
         sensing model (mode shape at each mount point); the ENGINE
         METcalc/Hx hold met and sensor points RIGID, so its bars
         read zero -- a known engine-model gap, not a failure.)
        dmdz: |l| model up to 409.4 nm, |e_piston| up to 92.97 nm per 100 nm figure
        dmdgrid: |l| model up to 463.3 nm, |e_piston| up to 432.5 nm per 100 nm figure
        NOTE (grid w rel): the truth is a plain forward poke at 100 nm;
        the harvest Jacobian is central-FD at 1e-06 (BaseUnits).  The
        engine's grid-surface response is NOT proportional across
        sub-um amplitudes on metre-based prescriptions (central FD at
        the harvest delta reproduces the Jacobian to ~2e-5; forward
        100 nm deviates ~15-25% on e2e, <5% on the mm fixture) --
        flagged for engine review (FreeForm grid intersection path,
        2026-07-19).  x and z channels do not show this.
        control u = bodies [1  2  3  4  5  6  7  8] -> dwdu 59750 x 48 (multi-field stack) for the simulator

    artifacts: s6_compare.mat, frames in s6_frames/, s6_compare.gif, this report
    next stage: run_simulator (single-step static estimator + controller over a time history).
