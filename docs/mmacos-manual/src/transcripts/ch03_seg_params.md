    P.seg = struct( ...
        'variants', ["pie" "hex2"], ... % s3 builds ALL of these ->
        ...                         % e2e_pie.in / e2e_hex2.in
        'variant',  "pie", ...      % the one stages 4-6 consume (one-knob
        ...                         % switch; both artifacts always exist)
        'gap_mm',         25, ...   % inter-segment gap (BaseUnits are m in
        ...                         % this family -- the runner converts)
        'dofs',           6, ...    % per-segment DOFs
        'meas_config',    2, ...    % edge sensors: 1 inner edges, 2 all
        'emit_apertures', true, ... % declare physical PolyApVec boundaries
        'ap_pad',         0, ...    % 0 = physical edge (gap rays clip)
        'grid_npts',      128, ...  % segmented-source sampling (e5-corpus
        ...                         % density; the parent's 41 gives only
        ...                         % ~50 rays/segment at 2 rings -- Dave
        ...                         % 2026-07-18: sampling seemed coarse)
        'model_size',     512);     % segmented stages need the big model
