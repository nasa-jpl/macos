    RINGS   = 1;          % 1 -> 7 segments, 2 -> 19
    GRID    = 'Hex';      % segment tiling / ray grid (Hex|Pie|...)
    GAP     = 50;         % inter-segment gap, mm
    MEASCFG = 2;          % edge sensors: 1 inner edges, 2 all adjacencies
    NF      = 6;          % MET fiducials on the hub (Dave: >=3, likely 6)
    R_FID   = 590;        % hub fiducial ring radius, mm -- fiducials must
                          % MOUNT ON M2 near its edge (~25 mm inside the
                          % 615 mm ApVec rim; no structure beyond)
    R_EXTRA = 100;        % aft ("M3") launcher ring radius, mm -- must hug
                          % that body (elt 11 is a Return with no ApVec/
                          % lMon to auto-size from; aft bench is ~90 mm)
    EDGE_OFF = 5;         % launcher clearance outward of the segment edge, mm
    SIG_ROT   = 1e-6;     % prior (deploy) uncertainty: rad per rot DOF
    SIG_TRANS = 1e-6;     % ... metres per trans DOF
    SIG_EDGE  = 1e-9;     % edge-sensor noise, metres
    SIG_MET   = 1e-9;     % MET gauge noise, metres
    NMC     = 200;        % Monte-Carlo draws for the cross-check
    MODEL   = 512;
