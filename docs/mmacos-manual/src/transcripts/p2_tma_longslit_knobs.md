    % ---- spec (Joe's sheet, deck_dyson spec slide)
    P.f_m            = 0.330;
    P.D_m            = 0.183;
    P.strip_half_deg = 4.7;
    P.slit_m         = 0.054;
    P.pixel_m        = 18e-6;
    P.telecentric_deg = 0.5;
    P.cone_fnum      = [1.7 1.8];
    P.work_dist_m    = 0.250;
    P.spec_joe   = struct('smile_px', 0.1, 'keystone_px', 0.1, 'crf_px', 1.5, 'srf_px', [1.5 2.0], 'eip', 0.75);
    % Bradley 2024 Table 1, in MICROMETRES (its own units: SRF in CO-ADDED pixels of 2 x 18 um -- "64.8 um (1.8 co-added
    % pixels)"; smile "1.8 um (5 % of co-added pixel)" = 0.10 of OUR 18-um pixel, not 0.05; keystone 1.8 um (10 % of a
    % pixel); ARF / CRF 50.4 um (2.8 px)).  Until 2026-10-07 these were compared in our 18-um pixels: every "fail (1.8)" on
    % SRF and the "< 0.05" smile were units errors (addendum 49 step 0).
    P.spec_paper = struct('smile_um', 1.8, 'keystone_um', 1.8, 'srf_um', 64.8, 'crf_um', 50.4, 'arf_um', 50.4);
    % ---- seed: SBG VSWIR Fig. 4b (f 345), digitised +-5 mm / +-2 deg
    P.seed_f_m       = 0.345;
    P.seed_legs_m    = [0.270 0.257 0.313];
    P.aoi_deg        = [30 37 15];
    P.turn           = [+1 -1 +1];
    P.stop           = 'M2';
    P.axis_side      = [+1 +1 +1];
    P.axis_theta_deg = [];
