    % ====================  USER DESIGN CHOICES  ==========================
    o.EFL_MM     = 100;          % focal length (mm)
    o.FNUM       = 5;            % f-number
    o.HALF_FIELD = 20;           % semi-field (deg)
    o.CROWN      = 'SK16';       % outer elements (Taylor's classic pair: SK16 / F4)
    o.FLINT      = 'F4';         % central element
    o.PETZVAL    = 0.2;          % Petzval sum x EFL (0 = flat; 0.2-0.3 is the classic compromise)
    o.T_MM       = [5 2 5];      % centre thicknesses (mm)
    o.LAMBDA_NM  = 587.5618;     % design wavelength (d); color is set by the seed
    o.CAL_FIELDS = [0 0.7 1];    % CALIB fields as fractions of HALF_FIELD
    o.N_CYCLES   = 3;            % CALIB cycles (a radii round + a spacings round each)
    o.FLOOR_MM   = 1.0;          % smallest axial gap, glass or air (mm)
    % ---------------------------------------------------------------------
    o.NGRID  = 45;
    o.MODEL  = 128;
    o.OUTDIR = fullfile(exdir, 'out');
    % =====================================================================
