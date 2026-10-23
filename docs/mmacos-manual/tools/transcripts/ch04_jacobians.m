% ch04_jacobians.m -- Chapter 4 (Sensitivities): the e2e flow's stage 4, run UNCHANGED in the thread
% directory (MMMAN_THREAD) on the e2e_pie.in ch03_segmentation.m left there.  Stage 4 harvests the
% linear model's three channels over the +-2' science field (center + 4 corners) through the general
% runner design/runners/run_sensitivities.m: dwdx (rigid body, 6 DOF per element), dwdz (segment-
% local Zernike figure modes), dwdgrid (per-segment grid influence functions).  Emits the stage's
% sensitivity report, the figures, and the numbers the prose quotes, parsed from the report's
% conditioning table (so they are the numbers the reader sees).  Model size 512: a MATLAB of its own.

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

rep = fileread(fullfile(work, 's4_report.txt'));
mm_listing('ch04_report', fullfile(work, 's4_report.txt'));

mm_figure('ch04_svspec',        fullfile(work, 's4_svspec.png'));
mm_figure('ch04_dwdx_channels', fullfile(work, 's4_dwdx_channels.png'));
mm_figure('ch04_dwdz_channels', fullfile(work, 's4_dwdz_channels.png'));
mm_figure('ch04_dwdg_channels', fullfile(work, 's4_dwdgrid_channels.png'));

% ---- the conditioning table, as printed ---------------------------------------------------------
for ch = ["dwdx" "dwdz" "dwdgrid"]
    t = regexp(rep, "\n\s+" + ch + "\s+(\d+)x(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)", 'tokens', 'once');
    assert(~isempty(t), 'ch04: no conditioning row for %s', ch);
    mm_pin("ch04_" + ch + "_rows", str2double(t{1}), 'sig5', t{1});
    mm_pin("ch04_" + ch + "_cols", str2double(t{2}), 'sig5', t{2});
    mm_pin("ch04_" + ch + "_rank", str2double(t{3}), 'sig5', t{3});
    mm_pin("ch04_" + ch + "_cond", str2double(t{6}), 'printed', t{6});
end
t = regexp(rep, 'segment-only: \S+\s+cond\+ (\S+)', 'tokens', 'once');
mm_pin('ch04_dwdx_segcond', str2double(t{1}), 'printed', t{1});
