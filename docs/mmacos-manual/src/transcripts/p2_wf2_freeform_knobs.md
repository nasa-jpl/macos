    % ====================  USER DESIGN CHOICES  ==========================
    APERTURE_M   = 1.0;        % aperture diameter D (m)
    PRIMARY_FNUM = 2.0;        % M1 f/#
    SYSTEM_FNUM  = 8.0;        % system f/#  (mag = SYSTEM/PRIMARY)
    BFD_M        = 0.25;       % back focal distance: focus this far behind M1 (m)
    FOV_ARCMIN   = 7.0;        % design HALF-field (arcmin): optimize + evaluate over +-this
    FREEFORM_ELT = 2;          % mirror(s) to make freeform (2 = secondary, away from pupil).
                               %   Try [1 2] for both; M1 is at the pupil (field-constant).
    FREEFORM_MODES = [5 6 11]; % ANSI Zernike DOF: 5 defocus, 6 astig, 11 spherical
    LAM          = 633e-9;     % wavelength (m)
    MODEL        = 256;        % diffraction grid model size
    DIFFRACTION_LIMIT = 0.071; % ~lambda/14 (waves RMS) -- the Marechal criterion
    % =====================================================================
