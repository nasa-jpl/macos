    % ====================  KNOBS  ========================================
    D          = 8.0;        % aperture diameter (m)  -- e5mono
    LAM        = 1e-6;       % wavelength (m)         -- 1 um
    MODEL      = 256;        % engine model size
    GRIDN      = 41;         % circular ray-grid points (-> ~1300 rays)
    FOV_ARCMIN = 2.0;        % design HALF-field (arcmin)
    DIFFRACTION_LIMIT = 0.07;% "well-corrected" RMS WFE (waves)

    % --- 0th-order all-sphere layout (|radii| m, spacings m, fold tilts deg)
    R    = [51.534, 8.871, 3.0];   % M1 concave, M2 convex, M3 concave
    TBET = [22.0,   28.0];         % M1->M2, M2->M3 along the folded chief ray
    TILT = [-7.2,   8.46,  12.0];  % fold about x (M1/M2 reproduce e5mono)

    % --- correction: e5mono's BornWolf Zernike set (piston dropped) ------
    MODES = [3 4 5 9 10 11 12 13 19 20 21 22 23 24 25];
    ZTYPE = 'BornWolf';
    ITERS = 200;
    % =====================================================================
