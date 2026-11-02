    % ====================  USER DESIGN CHOICES  ==========================
    APERTURE_M   = 1.0;     % aperture diameter D (m)
    PRIMARY_FNUM = 2.0;     % M1 f/#
    SYSTEM_FNUM  = 10.0;    % system f/#
    BFD_M        = 0.30;    % back focal distance: focus this far behind M1 (m)
    FOV_ARCMIN   = 4.0;     % evaluate WFE out to this field half-angle (arcmin)
    BALANCE_FIELD = true;   % field-balance the conics (radius+conic) for a WIDER
                            % corrected field -- giving up a little on-axis WFE to
                            % gain performance across the field is usually the better
                            % deal for a real instrument than a perfect-on-axis
                            % aplanat.  This balanced design is the saved deliverable.
    BALANCE_FOV  = 4.0;     % balance over +-this; keep it SHORT of FOV_ARCMIN
    % ---------------------------------------------------------------------
    LAM   = 633e-9;         % wavelength (m)
    MODEL = 256;            % diffraction grid model size
    DIFFRACTION_LIMIT = 0.07;   % ~lambda/14 (waves)
    % =====================================================================
