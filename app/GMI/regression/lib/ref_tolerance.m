function [tol, L] = ref_tolerance(rx, opts)
% REF_TOLERANCE  Tolerance for a comparison against a stored reference.
%
%   [tol, L] = ref_tolerance(rx, opts)
%
% Inputs:
%   rx    -- 'e2e_pie' or 'e5hex1'
%   opts  -- opts.ulp = how many units of round-off of the path are allowed
%
% Returns:
%   tol   -- opts.ulp * eps * L: the largest |a - b| allowed, in the base
%            units of the prescription (m for e2e_pie, mm for e5hex1);
%            eps is the machine epsilon, 2.22e-16
%   L     -- the chief-ray path length of the prescription, same units
%
% Why the tolerance follows the path.  An OPD is the small difference of
% two optical paths, so it is good to the round-off of the PATH, not of
% its own size: eps * L is the round-off of a number the size of L.
% (eps * L, not MATLAB's eps(L): eps(L) jumps by two at each power of
% two, and e2e_pie's 31.96 m sits just below one.)
% The references were made on Linux x86-64.  The first Mac (arm64) run,
% 2026-11-18, differed from them by 8 and 12 of these units on e2e_pie
% and by 8 and 6 on e5hex1 -- the same thing in each prescription's own
% units, though as plain numbers they are 8.5e-14 and 1.2e-10.  A fixed 1e-12
% passed the first and failed the second only because one deck is in
% metres and the other in millimetres.  (Dave's rule of 2026-10-19: no
% pin on digits the two architectures do not share.)
%
% L is the engine's own "Total Length" for the chief ray traced to the
% last element (CLI: old <deck>; srt <nElt>; 0,0), read 2026-11-18.  It
% sets a scale, so its first three figures are all that matter; re-read
% it only if a deck's layout changes.

    switch rx
        case 'e2e_pie'
            L = 31.9612075129113;        % m,  e2e_pie_met.in, element 17
        case 'e5hex1'
            L = 6.91208464052923e4;      % mm, Rx_e5hex1.in,   element 13
        otherwise
            error('ref_tolerance:rx', 'no path length recorded for "%s"', rx);
    end
    tol = opts.ulp * eps * L;
end
