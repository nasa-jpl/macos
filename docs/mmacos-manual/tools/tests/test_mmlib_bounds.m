% test_mmlib_bounds.m -- mm_zero / mm_floor fail CLOSED: a value that breaks its bound errors instead
% of printing '< b' / '> b' (2026-10-23: ch06a printed 'w < 0.01' for a measured 0.945, because the
% helpers printed the bound unchecked and the regen script checks a pin only on runs AFTER it is
% recorded).  Red on the pre-fix helpers (2/4), green after (4/4).
%   matlab -batch "run('tools/tests/test_mmlib_bounds.m')"     (from docs/mmacos-manual)
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'transcripts', 'mmlib'));
ok = 0;
try, mm_zero('t', 'w 0.945 x', 'w (\S+) x', 1e-2); catch e, ok = ok + strcmp(e.identifier, 'mm_zero:bound'); end
try, mm_floor('t', 'c 5 x', 'c (\S+) x', 1e6); catch e, ok = ok + strcmp(e.identifier, 'mm_floor:bound'); end
s = mm_zero('t', 'w 1e-9 x', 'w (\S+) x', 1e-2);  ok = ok + strcmp(s, 'w < 0.01 x');
s = mm_floor('t', 'c 5e7 x', 'c (\S+) x', 1e6);  ok = ok + strcmp(s, 'c > 1e+06 x');
fprintf('MMLIB_BOUNDS %d/4\n', ok);
assert(ok == 4, 'mm_zero / mm_floor bound gate: %d/4', ok);
