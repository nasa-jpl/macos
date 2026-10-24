    sampling budget (P.samp, enforce = warn):
      dimple 13.08 px across at the mask plane (min 6)  -> PASS
      pinhole 13.08 px across at the mask plane (min 6; the dimple's rule)  -> PASS
      DM 96x96: 1.66 detector px per actuator (min 2)  -> NOT MET
      DM 48x48: 3.31 detector px per actuator (min 2)  -> PASS
    G3 30 nm rms working state (DM 96x96, lit 5048): detector |E|/|E_flat| on msk std 3.708e-16   (gate < 1e-12: A = |E0| is exact)
    G4 vector pair on 100 nm single-actuator pokes every 8th actuator (12117 pm rms on msk, peak 1.92 rad; 3.03% of msk beyond the one-frame fold): V rms error 2.574 pm (gate < 0.1% of the figure), one-frame exact I 9132 pm (non-vacuity: must exceed 10x)   -> PASS
    stations figure pdi_stations_S.png: residual vs the engine's field on msk 0.00 pm (flat), 7063.27 pm (30 nm rms working surface), gain 0.7829
    stations figure pdi_stations_V.png: residual vs the engine's field on msk 0.00 pm (flat), 272.12 pm (30 nm rms working surface), gain 0.9954
    stations figure pdi_stations_P.png: residual vs the engine's field on msk 0.00 pm (flat), 744.20 pm (30 nm rms working surface), gain 0.9808
    P: pinhole 2.00 lam/D, surround t 0.8626, 4 steps, |b|^2 flat, NITER 3; frames per measurement 4; eta_pin 0.8228; throughput (detected/incident, flat) 0.9493; visibility on the flat 0.5021; surrogate vs engine reference 2.11e-15; the flat reads 1.86e-15 rad rms, amplitude 2.81e-16 rel
    G5 P on 100 nm single-actuator pokes every 8th actuator (12117 pm rms on msk, 3.03% of msk beyond the one-frame fold): rms error 4.508 pm (gate < 0.1% of the figure)   -> PASS
    PF: P/SRI, reference = LP01 mode V 2.30 b 0.50 core radius 0.50 lam/D, coupling eta_c 0.7414 on the flat; pickoff 0.60, a 0.003506 (match 0.003719, budget 0.003506), 4 steps (ls), |kappa| flat; frames per measurement 4; eta_pin 0.8228; throughput (detected/incident, flat) 0.8448; visibility on the flat 0.8608; surrogate vs engine reference 2.11e-15; the flat reads 3.16e-16 rad rms, amplitude 1.59e-16 rel
    G5 PF on 100 nm single-actuator pokes every 8th actuator (12117 pm rms on msk, 3.03% of msk beyond the one-frame fold): rms error 0.000 pm (gate < 0.1% of the figure)   -> PASS
    G6 P/SRI fiber reference under the 30 nm rms working state: shape FIXED by construction (the mode's); coupling relative to the flat |kappa| 0.86262, arg +0.0200 rad (an amplitude scale + a piston; the solver takes kappa = 1 in |kappa| 'flat' mode)
    G6 reference motion under the 30 nm rms working state (pinhole-diffracted reference, flat -> state, on msk), by pinhole diameter:
       lam/D    total  |scale|    SHAPE   (total = |b1-b0|/|b0|; scale = best complex factor, a Strehl-class amplitude drop the |b|^2 frame / iteration absorbs; shape = the rest, what the solve must iterate out)
        0.50   0.1240   0.8771   0.0002
        1.00   0.1284   0.8728   0.0007
        1.50   0.1350   0.8662   0.0024
        2.00   0.1429   0.8586   0.0076
        3.00   0.1519   0.8534   0.0342
        2.00   0.1429   0.8586   0.0076   <- the ZWFS dimple
    G7 the pinhole reading at t = 1, 2.00 lam/D, frozen reference vs the stepped Zernike reading S on 20 nm sparse pokes: rel rms diff 3.02e-15   (gate < 1e-2: same instrument, same frames)
    registration: center-poke anchor at detector px (512.19, 512.19) <-> DM (-0.42, -0.42) mm
      parity search on the off-center poke (60,40): winner [1  2  1  1] sign +1, |corr| 0.980, runner-up 0.009 (gates >= 0.40, sep >= 0.30); all 8: 0.98 -0.00 -0.00 -0.01 -0.01 -0.00 -0.00 -0.00 
      center poke in the DM frame: raw peak gain 0.6948, corr(map, truth) 0.9832
    bench stage <t> min
