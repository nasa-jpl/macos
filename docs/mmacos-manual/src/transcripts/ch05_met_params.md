    P.met = struct( ...
        'nf',               6, ...  % fiducials on the hub (>=3, likely 6)
        'fid_rim_inset_mm', 25, ... % fiducial ring inset inside the hub rim
        'edge_off_mm',      5, ...  % launcher clearance outward of seg edge
        'min_sep_mm',       50, ... % min separation between ANY 2 launchers
        'aft_ring_r_mm',    232, ...% aft launcher ring radius (on the FOLD,
        ...                         % the on-axis bench body behind M1: its
        ...                         % beams reach the M2 fiducials THROUGH
        ...                         % the M1 hole -- a ring at M3's lateral
        ...                         % position has no hole line of sight).
        ...                         % Dave 2026-07-19: the stable aft
        ...                         % structure extends to the SM RADIUS
        ...                         % (M2 lMon = 0.232 m) -- maximizes the
        ...                         % ring lever arm.
        'hole_r_mm',        232, ...% EFFECTIVE central-hole radius for the
        ...                         % aft-beam clearance check: the SM
        ...                         % shadow radius (Dave: the hole will be
        ...                         % that large anyway -- the SM obscures
        ...                         % more than the current 0.122 m hole).
        ...                         % Enlarging the REAL Rx hole to match is
        ...                         % a P.tel change cascading s1->s5.
        ...                         % SCENARIO (Dave 2026-07-19): MET tracks
        ...                         % the post-WFSC DRIFT state, not deploy
        ...                         % errors -- HWO-class drift is << 1 nm
        ...                         % between updates, and the MET accuracy
        ...                         % roadmap runs to ~1 pm.  The HEADLINE
        ...                         % run_met numbers are dimensionless
        ...                         % (WEM + floor fraction) and invariant
        ...                         % to these scales; the sigmas below set
        ...                         % only the absolute forecast column.
        'sig_rot',   1e-10, ...     % prior drift sigma, rad per rot DOF
        'sig_trans', 1e-10, ...     % ... metres per trans DOF (0.1 nm class)
        'sig_edge',  1e-11, ...     % edge-sensor noise, metres (10 pm class)
        'sig_met',   1e-12, ...     % MET gauge noise, metres (1 pm target)
        'mc_draws',  200);          % Monte-Carlo acceptance draws
