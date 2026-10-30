    % ---- the DM and the scale ----------------------------------------
    NACT = 96;  PITCH = 1.0;  DM_MM = NACT*PITCH;      % 96 mm aperture
    s    = DM_MM/56;                                   % uniform scale vs v1
    LAM  = 6.328e-4;  QWP = 0.25;  THETAS = [0 45 90 135];
    MODEL = 1024;  NGRID = 385;  N_G = 384;  DX_G = 0.28;  % 3.57 px/actuator
    % NGRID 385 (Dave 2026-09-03, 'a detector capable of resolving the
    % imagery, up to 1 Mpix square'): the wavefront seed inherits the ray
    % grid's sampling, so the pupil image spans ~NGRID px at the detector.
    % 63 px gave detector Nyquist 31.5 cyc/pupil -- BLIND at the 96x96 DM's
    % 48 (the run of record for that finding is tg96_report_ng63.txt).
    % 385 px -> ~4 px/actuator, Nyquist ~192 cyc/pupil, inside the 1024^2
    % (1 Mpix) grid with FFT padding.
    % MODEL 1024 is LOAD-BEARING: mGridMat is NOT the model size -- 512 caps
    % grids at 256 (macos_param.txt) and the nGridMat= parse is unguarded
    % (BRIEF_gridmat_guard.md), so a 384 grid on model 512 corrupts the heap.
    POKE  = 50e-6;                                     % mm (50 nm commands)
