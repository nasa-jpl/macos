    % ====================  USER DESIGN CHOICES  ==========================
    D      = 6.605;                          % aperture (m) -- j18mono parent
    R      = [15.879722 1.778913 3.016227];  % |radii| M1 / M2(convex) / M3
    TBET   = [7.169041556 7.965313479];      % spacings (m)
    LAM    = 1.0e-6;                         % score at 1 um
    DPAST  = 1.2;                            % M4 past the TMA focus (m); > f4!
    R4     = 1.5;                            % M4 radius (concave; f4 = 0.75 m)
    PATCH_RAD_ASEC = 15;                     % coronagraph patch RADIUS (arcsec)
    SECT_MARGIN    = 0.01;                   % clearance margin (fraction of D)
    % =====================================================================
