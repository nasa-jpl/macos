    function P = e2e_params()
    %E2E_PARAMS  Every user knob of the end-to-end example, in one place.
    %   The e2e stages (s1_telescope ... s6_simulator) all read THIS file --
    %   hack it for your system, then re-run the stages in order.  Each
    %   stage consumes the previous stage's saved artifacts, so a knob
    %   change re-runs from the first stage it affects.
    %
    %   Stage 1 derives the telescope from the BASIC parameters below via
    %   the first-order solver (macos.design.tma_layout): aperture, both
    %   f/#s, the Cassegrain-feed magnification and the packaging fractions
    %   are free inputs.  The telescope is an ON-AXIS Korsch TMA taken
    %   SLIGHTLY OFF-AXIS: the science field is biased just far enough off
    %   the axis to carry the focal plane out of the beam (the runner
    %   sweeps the bias and picks the least that clears -- off-axis
    %   aberration grows ~quadratically with bias, so least is best).

    % ================= telescope basics (stage 1) =======================
    P.D_m          = 4.0;       % aperture diameter (m); packaging scales with D
    P.lambda_m     = 0.5e-6;    % design wavelength (m) -- visible telescope
    P.fov_arcmin   = 1.0;       % HALF-field (about the bias) the telescope
                                % holds; stage 2's imaging instrument pushes
                                % the corrected field toward P.inst.fov_arcmin
    P.primary_fnum = 1.75;      % M1 f/#  (f1 = fno1*D, R1 = 2*f1)
    P.system_fnum  = 18;        % final f/# at the FP  (EFL = fno_sys*D)
    P.secondary_mag  = 8;       % Cassegrain feed magnification m2 (>1).
                                % Lower mag relaxes M2's curvature (a
                                % field-aberration source) at the cost of
                                % a somewhat longer bench; the f/1.25+m2=16
                                % starting point could not reach VIS blur
                                % (Dave 2026-07-17: f/1.75, m2=8)
    P.int_focus_frac = -0.125;  % intermediate-focus z, fraction of D
                                % (NEGATIVE = in front of M1, between M1
                                % and M2 -- the field-stop / metrology-
                                % injection plane stays accessible)
    P.m3_behind_frac = 0.6;     % M3 vertex z, fraction of D behind M1
    P.fold_frac    = 0.075;     % 90-deg fold (FM) station z, fraction of D
                                % behind M1: the M2->M3 feed is folded into
                                % +x there, so M3 and everything after live
                                % on a flat bench BEHIND M1 (Dave: move the
                                % back end behind the primary)
    P.m3_tilt_deg  = 1.2;       % EXTRACTION TILT on M3, about the bench
                                % normal (the builder's tilt axis maps
                                % there through the fold isometry): the
                                % M3->FP return leaves the feed axis
                                % geometrically, so the fold clears WITHOUT
                                % field bias -- the bias sweep then settles
                                % near its minimum.  (Bias was the VIS
                                % killer: aberration ~ bias^2.)
    P.fold_margin  = 1.15;      % fold body radius = margin * local feed
                                % radius.  A slim mount: the M3->FP return
                                % passes the fold body with mm-scale
                                % daylight at the chosen tilt/bias -- 1.4x
                                % grazed the clearance judge
    P.fp_body_r    = 0.075;     % focal-plane BODY radius (m) for the
                                % clearance judge (detector + mount)
    P.bias_sweep_arcmin = [1 2 3 4 6];  % off-axis biases to explore; the
                                % least whose folded design fully CLEARS
                                % (only M2's central obscuration allowed)
                                % wins -- aberration grows ~bias^2
    P.hole_margin  = 1.3;       % M1 central hole radius = margin * the
                                % measured beam radius at the M1 plane
    P.hole_min_r_m = 0.232;     % FLOOR on the hole radius = the SM radius
                                % (M2 lMon, see the s2 report): the SM
                                % shadow obscures this much anyway, and the
                                % aft MET beams route through the hole to
                                % the M2 rim fiducials (Dave 2026-07-19).

    % ---- freeform refinement on top of the conic solve (the "+FF") -----
    P.modes     = [3 4 5 9 10 11 12 13 19 20 21 22 23 24 25];
    P.ztype     = 'BornWolf';   % Zernike ordering for the departure basis
    P.max_iters = 150;          % CALIB iterations per solve
    P.max_iters_ff = 300;       % ... for the joint freeform field solve
    P.dl_waves  = 0.071;        % Marechal diffraction-limit bar (waves RMS)

    % ---- engine sampling ----------------------------------------------
    P.model_size = 256;         % stage 1-2 engine model
    P.grid_npts  = 41;          % circular ray-grid points (~1300 rays)

    % ================= stage 2: imaging instrument ======================
    % A THREE-mirror bench relay behind the telescope focus (M4 weak
    % corrector at an intermediate conjugate / M5 collimator / M6 camera,
    % unit magnification -> final f/# preserved).  Field correction =
    % freeform on surfaces at staggered intermediate conjugates; no active
    % control in this part of the telescope (Dave 2026-07-17).  A 4-mirror
    % variant (extra weak corrector near the relayed pupil) was probed and
    % conditioned WORSE -- near-pupil weak correctors act common-mode and
    % collapse to near-zero field-solve rank.
    % The solve is JOINT: M1-M3 keep refining with the instrument.
