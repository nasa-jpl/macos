    % ====================  USER DESIGN CHOICES  ==========================
    o.EFL_MM   = 100;          % focal length at d (mm)
    o.FNUM     = 5;            % f-number (stop at the first surface)
    o.CROWN    = 'N-BK7';      % front element (positive)
    o.FLINT    = 'F2';         % rear element (negative); any catalog name
    o.T_MM     = [6 2.5];      % centre thicknesses crown, flint (mm)
    o.LAM_NM   = [486.1327 587.5618 656.2725];   % F, d, C (the achromat's pair is F and C)
    o.VENDOR_PART = 'PAC095AR.14';               % a lens_catalog part ('' = skip stage 5)
    % ---------------------------------------------------------------------
    o.NGRID    = 63;
    o.MODEL    = 128;
    o.OUTDIR   = fullfile(exdir, 'out');
    % =====================================================================
