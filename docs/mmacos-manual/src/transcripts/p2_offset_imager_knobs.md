        % ---- the rodgers3 challenge instance (defaults) ---------------------
        P.name       = 'rodgers3';
        P.tag        = 'oi';
        P.EPD_m      = 0.075;
        P.Fno        = 4.0;
        P.lambda_m   = 1.0e-6;

        P.box_deg    = [20 20];
        P.offset_deg = 22;
        P.nsolve     = 3;
        P.map_n      = 11;

        P.z_m1_m     = 0.6649568;
        P.spacings_m = [-0.7228968  0  0.7408280];
        P.bfd_m      = [];

        P.exit_dir     = [];        % report-only by default; T3 sets Mike's
        P.exit_tol_deg = 1.0;
        P.exit_wt      = 1e4;       % exit-direction residual weight, nm/deg
                                    % (an EQUALITY constraint solved in the
                                    % least-squares sense from S3 on -- a
                                    % boolean wall would freeze any stage
                                    % that STARTS outside tolerance)
        P.clear_m      = [0.050 0.035];
