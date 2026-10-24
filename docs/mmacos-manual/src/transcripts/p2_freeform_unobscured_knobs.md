    % ====================  KNOBS  ========================================
    D          = 8.0;        % aperture diameter (m)   -- e5mono heritage
    LAM        = 0.5e-6;     % 500 nm center wavelength: VISIBLE telescope
    MODEL      = 256;        % engine model size
    GRIDN      = 41;         % circular ray-grid points (~1300 rays)
    FOV_ARCMIN = 1.0;        % design HALF-field (arcmin) -- the visible DL
                             % bar is 35 nm RMS; hold a tighter field than
                             % sz_tma's 1-um +-2' and let [7] show the cost
    DIFFRACTION_LIMIT = 0.071;

    % --- 0th-order all-sphere layout (|radii| m, spacings m, fold tilts deg)
    % M1 f/3.2 (slower than the centered j18's f/1.2 -- gentle AOI), convex
    % M2 slows the cone to a real f/21 intermediate image between M2 and M3
    % (the metrology-injection point), M3 reimages it clear of the beam.
    R    = [51.534, 8.871, 3.0];   % M1 concave, M2 CONVEX, M3 concave
    TBET = [22.0,   28.0];         % M1->M2, M2->M3 along the folded chief
    TILT = [-7.2,   8.46,  12.0];  % fold about x (M1/M2 = e5mono)

    % --- correction: e5mono's BornWolf Zernike set (piston dropped) ------
    MODES = [3 4 5 9 10 11 12 13 19 20 21 22 23 24 25];
    ZTYPE = 'BornWolf';
    ITERS = 200;
    % =====================================================================
