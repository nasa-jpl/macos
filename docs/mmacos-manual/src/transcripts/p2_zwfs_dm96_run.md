    ---- bench ----
    MASK_TRIM re-scanned: 1.2342 mm (mask-plane peak/sum 1.868e-02); the sheet asked for a scan.  <s> s
    deck run_test.in: TO elt 7, FocalMask 18, Detector 22; mask sandwich spheres zElt 355.720 / 355.720 mm (SYMMETRIC)
    mask: phase 1.5710 rad (|c| 1.414), dimple 2.000 lam/D = 5.6500e-03 mm = 13.08 px at the mask plane; msk 19601 px
    focal spot: 0.432 um per px at the mask plane, 6.54 px per lam F/D; the dimple encloses 82.6% of the flat DM's focal-plane light
    G1 mask-sandwich round trip (unmasked, entrance -> exit sphere): 1.849e-15   (gate < 1e-12)
    G2 reference-wave surrogate T(D Ti(E0)) vs the engine's Eb on msk: 2.107e-15   (gate < 1e-10)
    reference wave |Eb|/|E0| vs pupil radius (20 bins, 0..1): 1.3075 1.3015 1.2895 1.2717 1.2488 1.2200 1.1859 1.1469 1.1038 1.0577 1.0075 0.9543 0.8989 0.8421 0.7834 0.7233 0.6629 0.6029 0.5440 0.4857 
    frame: ray magnification 9.9752 DM-mm per detector-mm, detector px 6.0492e-02 mm -> 0.6034 DM-mm per px
    camera: pupil image 9.56 mm across (158 modeled px at 60.5 um); sCMOS 2048x2048 at 6.50 um -> 1470 raw px across the pupil; binned 4 = 368, and binning 9 = 163 lands nearest the modeled 158
      snapshot analyzer: polarization sCMOS (micro-polarizer array, 0/45/90/135) at 3.45 um -> 2770 px across the pupil, 1385 per orientation (every 2nd pixel in each direction on the 2x2 superpixel -- N/2 across a line, not N/4); four simultaneous frames, no rotating stage
    sampling budget (P.samp, enforce = warn):
      dimple 13.08 px across at the mask plane (min 6)  -> PASS
      DM 96x96: 1.66 detector px per actuator (min 2)  -> NOT MET
      DM 48x48: 3.31 detector px per actuator (min 2)  -> PASS
    G3 30 nm rms working state (DM 96x96, lit 5048): detector |E|/|E_flat| on msk std 3.708e-16   (gate < 1e-12: A = |E0| is exact)
    G4 vector pair on 100 nm single-actuator pokes every 8th actuator (12117 pm rms on msk, peak 1.92 rad; 3.03% of msk beyond the one-frame fold): V rms error 2.574 pm (gate < 0.1% of the figure), one-frame exact I 9132 pm (non-vacuity: must exceed 10x)   -> PASS
    stations figure run_stations_S.png: residual vs the engine's field on msk 0.00 pm (flat), 7063.27 pm (30 nm rms working surface), gain 0.7829
    stations figure run_stations_V.png: residual vs the engine's field on msk 0.00 pm (flat), 272.12 pm (30 nm rms working surface), gain 0.9954
    registration: center-poke anchor at detector px (512.19, 512.19) <-> DM (-0.42, -0.42) mm
      parity search on the off-center poke (60,40): winner [1  2  1  1] sign +1, |corr| 0.980, runner-up 0.009 (gates >= 0.40, sep >= 0.30); all 8: 0.98 -0.00 -0.00 -0.01 -0.01 -0.00 -0.00 -0.00 
      center poke in the DM frame: raw peak gain 0.6948, corr(map, truth) 0.9832
    bench stage <t> min
