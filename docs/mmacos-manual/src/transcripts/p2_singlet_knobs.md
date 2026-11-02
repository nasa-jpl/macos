    % ====================  USER DESIGN CHOICES  ==========================
    o.EFL_MM    = 100;        % focal length (mm)
    o.FNUM      = 5;          % f-number EFL / EPD (the stop at the first surface)
    o.GLASS     = 'N-BK7';    % any name in macos_f90/macos_glass_list.txt
    o.LAMBDA_NM = 587.5618;   % the design wavelength (He d line)
    o.THICK_MM  = 4;          % centre thickness (mm)
    o.OBJECT_MM = Inf;        % object distance in front of the lens (Inf = collimated)
    % ---------------------------------------------------------------------
    o.NGRID     = 63;         % rays across the entrance pupil
    o.MODEL     = 128;        % engine model size
    o.NSCAN     = 41;         % lenses in the bending scan (0 = skip the scan and its figure)
    o.OUTDIR    = fullfile(exdir, 'out');
    % =====================================================================
