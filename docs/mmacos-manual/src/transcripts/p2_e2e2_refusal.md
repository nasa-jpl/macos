    ====================================================================
     E2E2 stage 1: axial Korsch TMA | D=3.00 m | f/20 | 500 nm | box 0.4 deg
    ====================================================================

    [1] conic-solver fixture gate (tma_fixture.json)
        K fixture  = [  -0.9747480   -1.0345250   -4.5527288]
        K solved   = [  -1.0004068   -1.3371575   -4.4894914]   max|dK| = 3.03e-01 (bar 1e-06)
        EFL 49.0000 m (fixture 5.0342), t_focus -1.866667 m (fixture 1.064384)

    REFUSED: CONIC SOLVER GATE FAILED: max|dK| = 3.026e-01 > 1e-06.  STOP AND FIX -- do NOT widen this bar and do NOT hand-edit the fixture; regenerate it with optical_design/make_tma_fixture.py and explain what moved.
