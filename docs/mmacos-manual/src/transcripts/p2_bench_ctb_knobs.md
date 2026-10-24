    P.F_OAP = [2500 1524 1143 1350 675 635 635 762];  % DST2R parent focal lengths
    P.AOI   = 5;                                       % OAP/DM angle of incidence, deg
    P.R_OAP = 75;                                       % OAP clear-aperture radius, mm
    P.R_DM  = 22.5;                                     % DM aperture radius, mm (stop)
    P.FILL  = 0.95;                                     % source NA fills FILL of DM
    P.WAVLEN= 500e-6;                                   % 500 nm in mm
    P.r     = P.F_OAP / cosd(P.AOI)^2;                  % conjugate distances (mm)
    % free (collimated-space) leg lengths -- set generously, no overlap:
    P.L_DM1=600; P.L_DM2=500; P.L_O2=700; P.L_APOD=350; P.L_O4=500;
    P.L_LYOT=350; P.L_O6=450; P.L_BACK=300; P.L_O8=450;
    P.AP    = P.R_DM*P.FILL / P.r(1);                   % source numerical aperture
    % aluminium coating for the polarization study (index n - i*kappa):
    P.AL_N = 1.2;  P.AL_K = 7.0;  P.AL_T = 1e-4;        % ~100 nm Al
