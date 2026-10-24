    P.inst = struct( ...
        'type', "offner", ...       % "offner" (concentric ring-field 1:1
        ...                         % relay -- no tilted powered surfaces,
        ...                         % zero Seidel over the ring; the biased
        ...                         % patch is a ring arc, and several
        ...                         % small-field instruments can pick off
        ...                         % arcs of the same ring) | "zigzag"
        ...                         % (the tilted-sphere bench relay)
        'offner_R',   2.0, ...      % Offner concave radius (m); convex = R/2
        'offner_h',   0.25, ...     % ring radius (m): object/image offset
        'fov_arcmin', 2.0, ...      % widened half-field target
        'field_center', "auto", ... % "auto" (Dave 2026-07-18): after the
        ...                         % pass-1 solve, point the nominal chief
        ...                         % at the y-CENTROID of the region where
        ...                         % raw WFE < field_center_thresh on the
        ...                         % +-map_fov map, and re-solve once (the
        ...                         % [4f] scan then verifies optimality).
        ...                         % "manual" keeps field_dy_arcmin as-is.
        ...                         % set_field_bias is a +y scalar, so
        ...                         % only the y-centroid is adoptable; the
        ...                         % x-centroid is reported.
        'field_center_thresh', 0.02, ...  % waves, the "good region" bar
        'map_fov_arcmin', 3, ...    % WFE-map half-width (Dave: +-3')
        'map_n', 13, ...            % map sampling (13x13 = 0.5' pitch)
        'field_dy_arcmin', -0.7, ...% STARTING science-field-center shift
        ...                         % (+y bias
        ...                         % units) applied to the s1 bias for the
        ...                         % instrument stage: the s2 WFE map's
        ...                         % sweet spot sat below (0,0) (Dave
        ...                         % 2026-07-18: center at (0,-0.7)' or
        ...                         % lower).  BOTH candidates were solved
        ...                         % end-to-end and kept (s2_variants/):
        ...                         %  -0.70' -> best FIELD EDGE: worst +-2'
        ...                         %    0.0231 -tilt, 2' ring Strehl 0.965
        ...                         %  -1.05' -> FLATTER interior (Strehl
        ...                         %    >=0.985 through 1.5') but the 2'
        ...                         %    edge pays (0.038 -tilt, 0.945)
        ...                         % Adopted -0.7 on the worst-field
        ...                         % criterion; switch here if the
        ...                         % instruments live inside the ring.
        ...                         % The runner prints a center-scan table
        ...                         % after the solve; clearance is
        ...                         % re-checked at the shifted bias.
        'dpast_m',    0.45, ...     % M4 corrector past the telescope focus
        'R_m',        [], ...      % [] = DERIVE: M4 weak (20 m); M5 from
        ...                         % the collimator condition f5 = dpast +
        ...                         % leg1 (exact, whatever the telescope
        ...                         % conjugates); M6 = M5 (unit mag)
        'legs_m',     [1.5 2.0], ...       % M4->M5, M5->M6 (M6->FP derives)
        'tilt_deg',   [5 -5 5], ...        % bench zigzag tilts (astig ~
        ...                         % tilt^2; 6 deg OVERLAPS the in/out
        ...                         % beams at these legs and diverges)
        'modes',      [], ...       % relay Zernike modes ([] -> P.modes)
        'nfield_svd', 5, ...        % NxN dense field grid for the SVD solve
        ...                         % stages + final scoring (CALIB caps at 12)
        'max_iters_ff', 150);       % joint freeform field-solve iterations
