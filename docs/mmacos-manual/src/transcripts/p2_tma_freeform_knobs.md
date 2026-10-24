    % ====================  USER DESIGN CHOICES  ==========================
    APERTURE_M   = 1.0;        % aperture diameter D (m)
    % Coaxial Korsch TMA seed (validated proof_korsch f/8 fixture, scales with D):
    M1_RADIUS_M  = 8.0;  M1_SPACING_M = 3.0;   % M1 radius + M1->M2 spacing
    M2_RADIUS_M  = 2.0;  M2_SPACING_M = 4.5;   % M2 radius + M2->M3 spacing
    M3_RADIUS_M  = 4.0;                        % M3 radius (M3->focus derived)
    FOV_ARCMIN   = 3.0;        % design HALF-field (arcmin): optimize + evaluate over +-this
    FREEFORM_ELT = [2 3];      % mirrors to make freeform (M2, M3 -- both off the pupil)
    FREEFORM_MODES = [5 6 7 8 11];   % ANSI Zernike DOF: defocus, astig, coma, spherical
    LAM          = 633e-9;     % wavelength (m)
    MODEL        = 256;        % diffraction grid model size
    DIFFRACTION_LIMIT = 0.071; % ~lambda/14 (waves RMS) -- the Marechal criterion
    % =====================================================================
