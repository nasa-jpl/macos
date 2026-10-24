% ch04_jacobians.m -- Chapter 4 (Sensitivities): the e2e flow's stage 4, run UNCHANGED in the thread
% directory (MMMAN_THREAD) on the e2e_pie.in ch03_segmentation.m left there.  Stage 4 harvests the
% linear model's three channels over the +-2' science field (center + 4 corners) through the general
% runner design/runners/run_sensitivities.m: dwdx (rigid body, 6 DOF per element), dwdz (segment-
% local Zernike figure modes), dwdgrid (per-segment grid influence functions).  Emits the stage's
% sensitivity report, the figures, and the numbers the prose quotes, parsed from the report's
% conditioning table (so they are the numbers the reader sees).  Model size 512: a MATLAB of its own.
%
% Pinned as BOUNDS, not values (dw_dx's finite-difference floor: ~1e-6 abs, up to 3e-5 on ~0
% entries, tests/tDwDx.m; arm64 flips the sign of such entries), one line each:
%   ch04_dwdx_svmin   < 5e-5  dwdx's smallest kept singular value sits at the floor
%   ch04_seg1_rz      < 5e-5  center-segment clocking moves (almost) nothing: a ~0 column
%   ch04_dwdx_cond    > 1e6   sv_max over a floor-level sv_min: only its order is physical
%   ch04_dwdx_segcond > 1e6   the same for the 42 segment columns (segment clocking is in them)
%   ch04_dwdx_deficit < 5     the rank cut sits between the floor and two exact nulls (the fold's Rx/Ty)
% Everything else (sizes, dwdz/dwdgrid rank and condition) is a signal: sig5 or printed.

root = getenv('MMACOS_ROOT');
work = getenv('MMMAN_THREAD');
if ~isfile(fullfile(work, 'e2e_pie.in'))
    error('ch04:thread', 'no stage-3 products in %s -- run ch03_segmentation first', work);
end
copyfile(fullfile(root, 'templates', '80_end_to_end', 'e2e', 's4_jacobians.m'), work);
addpath(fullfile(root, 'src'));  addpath(fullfile(root, 'design', 'src'));
addpath(fullfile(root, 'design', 'runners'));  addpath(fullfile(root, 'sensitivities'));

t0 = tic;
run(fullfile(work, 's4_jacobians.m'));
fprintf('stage 4 wall time %.0f s\n', toc(t0));
work = getenv('MMMAN_THREAD');

% ---- the report, with its numerical zeros printed as bounds -------------------------------------
% dw_dx is a finite difference: its noise floor is ~1e-6 absolute, up to 3e-5 on the ~0 entries
% (tests/tDwDx.m, the 2026-10-09 re-pin and CCMac's arm64 gates).  So dwdx's smallest kept singular
% value and the near-inert rz column of the center segment are zeros at that floor (bound 5e-5, the
% tDwDx AbsTol), and the dwdx condition numbers, sv_max over that floor, are only known to exceed
% 1e6; dwdx's rank is decided at the floor too and is pinned as a bound on the rank deficit
% (how many of the 90 columns the cut drops).  dwdz and dwdgrid are well conditioned: signals.
rep = fileread(fullfile(work, 's4_report.txt'));
rep = mm_floor('ch04_dwdx_cond', rep, 'dwdx\s+\d+x\d+\s+\d+\s+\S+\s+\S+\s+(\S+)', 1e6);  % before the
rep = mm_zero('ch04_dwdx_svmin', rep, 'dwdx\s+\d+x\d+\s+\d+\s+\S+\s+(\S+)', 5e-5);     % zero, whose '< b' is 2 tokens
rep = mm_floor('ch04_dwdx_segcond', rep, 'segment-only: \S+\s+cond\+ (\S+)', 1e6);
rep = mm_zero('ch04_seg1_rz', rep, '\n\s+1\s+\S+\s+\S+\s+(\S+)\s+\S+\s+\S+\s+\S+\s*\n', 5e-5);
t = regexp(rep, 'dwdx\s+\d+x(\d+)\s+(\d+)\s', 'tokens', 'once');
mm_pin('ch04_dwdx_deficit', str2double(t{1}) - str2double(t{2}), 'abs:5', '< 5');
rep = regexprep(rep, '(dwdx\s+\d+x\d+\s+)\d+(\s)', '$1 <90$2', 'once');
mm_transcript('ch04_report', strtrim(rep));

mm_figure('ch04_svspec',        fullfile(work, 's4_svspec.png'));
mm_figure('ch04_dwdx_channels', fullfile(work, 's4_dwdx_channels.png'));
mm_figure('ch04_dwdz_channels', fullfile(work, 's4_dwdz_channels.png'));
mm_figure('ch04_dwdg_channels', fullfile(work, 's4_dwdgrid_channels.png'));

% ---- the conditioning table, as printed ---------------------------------------------------------
t = regexp(rep, '\n\s+dwdx\s+(\d+)x(\d+)', 'tokens', 'once');
mm_pin('ch04_rows', str2double(t{1}), 'sig5', t{1});
mm_pin('ch04_dwdx_cols', str2double(t{2}), 'sig5', t{2});
for ch = ["dwdz" "dwdgrid"]
    t = regexp(rep, "\n\s+" + ch + "\s+(\d+)x(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)", 'tokens', 'once');
    assert(~isempty(t), 'ch04: no conditioning row for %s', ch);
    mm_pin("ch04_" + ch + "_cols", str2double(t{2}), 'sig5', t{2});
    mm_pin("ch04_" + ch + "_rank", str2double(t{3}), 'sig5', t{3});
    mm_pin("ch04_" + ch + "_cond", str2double(t{6}), 'printed', t{6});
end
