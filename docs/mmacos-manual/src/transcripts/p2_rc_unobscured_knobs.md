    % ====================  USER DESIGN CHOICES  ==========================
    APERTURE_M   = 1.0;     % aperture diameter D (m)
    PRIMARY_FNUM = 1.5;     % M1 f/#    (FASTER -> smaller secondary -> less off-axis)
    SYSTEM_FNUM  = 20.0;    % system f/# (SLOWER -> smaller secondary -> less off-axis)
    OFF_AXIS_M   = [];      % off-axis decenter (m).  [] = AUTO (minimum to clear);
                            % a number sets it explicitly (>= the minimum).
    % ---- secondary choices ----------------------------------------------
    BFD_M        = 0.30;    % back focal distance: focus this far BEHIND M1 (m)
    FOV_ARCMIN   = 0.5;     % design HALF-field (arcmin): optimize + evaluate over +-this
    LAM          = 633e-9;  % wavelength (m)
    MODEL        = 256;     % diffraction grid model size
    DIFFRACTION_LIMIT = 0.07;   % ~lambda/14 (waves)
    % =====================================================================
