function [pass, max_diff, msg] = compare_within(a, b, tol, label)
% COMPARE_WITHIN  Hard absolute-tolerance equality check.
%
%   [pass, max_diff, msg] = compare_within(a, b, tol, label)
%
% Inputs:
%   a, b   -- arrays of equal shape (numeric, real)
%   tol    -- max allowed |a - b| (absolute)
%   label  -- short string for the failure message
%
% Returns:
%   pass     -- true if max(|a - b|) <= tol
%   max_diff -- max(|a - b|) (the actual measurement)
%   msg      -- one-line description suitable for printing
%
% The check itself is absolute.  The caller chooses tol: 0 for two calls
% of the same build (they must agree bit for bit), and
% lib/ref_tolerance.m for a comparison against a stored reference, which
% scales with the prescription's path length so that it holds on a
% second architecture.

    if ~isequal(size(a), size(b))
        pass = false;
        max_diff = NaN;
        msg = sprintf('[%s] shape mismatch: %s vs %s', label, ...
                      mat2str(size(a)), mat2str(size(b)));
        return
    end

    d = abs(a(:) - b(:));
    max_diff = max(d);
    if isempty(max_diff)
        max_diff = 0;
    end

    pass = max_diff <= tol;
    if pass
        msg = sprintf('[%s] PASS  max|a-b| = %.3e (tol %.3e)', ...
                      label, max_diff, tol);
    else
        msg = sprintf('[%s] FAIL  max|a-b| = %.3e > tol %.3e', ...
                      label, max_diff, tol);
    end
end
