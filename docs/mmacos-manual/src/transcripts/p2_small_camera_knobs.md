    % ====================  USER DESIGN CHOICES  ==========================
    o.DECK         = '';           % the deck ('' = small_camera.in beside this runner)
    o.D_MM         = 350;          % aperture (mm)
    o.EFL_MM       = 3500;         % camera focal length (mm)
    o.PRIMARY_FNUM = NaN;          % primary f-number (NaN = the deck's)
    o.PIXEL_UM     = 13;           % detector pitch (um)
    o.NPIX         = 2048;         % pixels across the field (a line array)
    o.BAND_NM      = [500 700 900];% band; the middle entry is the design wavelength
    o.GLASSES      = {};           % corrector + filter glasses ({} = the catalog matches of the deck's indices)
    o.CAL_FIELDS   = [0 0.5 1];    % CALIB fields as fractions of the half field
    o.N_ROUNDS     = 4;            % stage-4 rounds (CALIB, then the first-order re-solve; the best is kept)
    % ---------------------------------------------------------------------
    o.NGRID  = 33;
    o.MODEL  = 128;
    o.OUTDIR = fullfile(exdir, 'out');
    % =====================================================================
