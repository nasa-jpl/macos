    % ====================  USER DESIGN CHOICES  ==========================
    D          = 6.605;             % aperture (m)
    LAM        = 0.5e-6;            % center wavelength: 500 nm (visible)
    SYS_FNUM   = 20;                % system f/# (EFL = 132 m; lambda*F = 10 um
                                    % -> Nyquist at ~5 um detector pixels)
    F1_STEPS   = [1.5 2.0 2.5 3.0]; % the slower-M1 ladder (primary f/#)
    F1_PICK    = 2.5;               % the design point (Dave 2026-07-06:
                                    % "go with f/2.5" -- comfortably under
                                    % the AOI rule with margin for the
                                    % coronagraph, at acceptable tube)
    FEED_FNUM  = 10;                % Cassegrain FEED f/# (= f1 * m2), held
                                    % constant across the ladder (j18: 1.2 x
                                    % 8 = 9.6).  Holding m2 instead drives
                                    % the feed toward the system f/# as M1
                                    % slows -- the M3 relay degenerates to
                                    % 1:1 and the clearing decenter blows up
                                    % (seen at m2=8: f/2.5 needs > 2.5*D).
    FIELD_RAD  = 0.5;               % science-field RADIUS (arcmin) for the
                                    % conic solve + the demo's FOV ladder
    AOI_MAX    = 15;                % per-mirror AOI SPREAD preference (deg,
                                    % coronagraph polarization rule)
    MARGIN     = 0.05;              % section clearance margin (x D)
    % =====================================================================
