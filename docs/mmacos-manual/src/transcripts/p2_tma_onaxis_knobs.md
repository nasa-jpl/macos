    % ====================  USER DESIGN CHOICES  ==========================
    APERTURE_M    = 1.0;    % aperture diameter D (m)
    PRIMARY_FNUM  = 1.5;    % M1 f/#  (f1 = PRIMARY_FNUM*D, R1 = 2*f1)
    SYSTEM_FNUM   = 8.0;    % system f/#  (EFL = SYSTEM_FNUM*D)
    SECONDARY_MAG = 8.0;    % Cassegrain feed magnification (> 1)
    INT_FOCUS_D   = -0.125; % intermediate-focus z, units of D (between M1 & M2)
    M3_BEHIND_D   = 0.6;    % tertiary this far behind the primary, units of D
    FOV_ARCMIN    = 2.0;    % half-field optimized/evaluated about the bias
    BIAS_SWEEP_ARCMIN = [1 2 3 4 6];   % off-axis biases to explore (arcmin)
    % ---------------------------------------------------------------------
    LAM   = 633e-9;         % wavelength (m)
    MODEL = 256;            % diffraction grid model size
    DIFFRACTION_LIMIT = 0.07;   % ~lambda/14 RMS (waves)
    % =====================================================================
