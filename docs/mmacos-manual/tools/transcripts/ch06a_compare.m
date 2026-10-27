% ch06a_compare.m -- Chapter 6 (Closing the loop), part 1: the e2e flow's stage 6 in the thread
% directory (MMMAN_THREAD), on the products of ch03-ch05.  Stage 6 (design/runners/run_compare.m)
% pokes every degree of freedom of the linear model in turn (54 rigid, 56 figure, 42 grid) and
% compares the engine's response with the linear prediction, for the wavefront w and the two
% measurements, the MET beam lengths l and the edge sensors e.  Run UNCHANGED (it writes only beside
% itself).  The chapter prints the report's header and its [3] summary, not the per-poke table.
%
% Pinned as BOUNDS (agreement ratios at the noise of the comparison), one line each:
%   ch06_x_l   < 1e-5   worst rigid-body l: METcalc vs dldx, FD-noise class
%   ch06_x_e   < 1e-6   worst rigid-body e: finite-rotation kinematics vs dedx
% SIGNALS (printed), the worst w of each channel: on the live chain they are ~1, not noise --
% run_compare's engine truth takes the engine's OPD default ('chief', since 2026-10-10) while s4
% harvested at 'mean' (PLAN_CONSOLIDATION 1.2 (kkk): Seg 1 Tz +100 nm, engine 'chief' 188.9 nm rms,
% 'mean' 61.98, linear 61.98).  (A first version bounded x/z w at 1e-2 and the helpers printed the
% bound unchecked; they fail closed since dbc2aea2.)

root = getenv('MMACOS_ROOT');
work = getenv('MMMAN_THREAD');
for f = ["e2e_pie.in" "e2e_pieHx.m" "s4_jacobians.mat" "e2e_pie_met.mat" "s4_grid.in"]
    if ~isfile(fullfile(work, f))
        error('ch06a:thread', 'no %s in %s -- run ch03-ch05 first', f, work);
    end
end
copyfile(fullfile(root, 'templates', '80_end_to_end', 'e2e', 's6_compare.m'), work);
copyfile(fullfile(root, 'templates', '80_end_to_end', 'e2e', 'e2e_params.m'), work);   % the stage's knobs as they are now (P.piston, 1.2 (www))
addpath(fullfile(root, 'src'));  addpath(fullfile(root, 'design', 'src'));
addpath(fullfile(root, 'design', 'runners'));

t0 = tic;
run(fullfile(work, 's6_compare.m'));
fprintf('stage 6 wall time %.0f s\n', toc(t0));
work = getenv('MMMAN_THREAD');

rep = fileread(fullfile(work, 's6_report.txt'));
k2 = strfind(rep, sprintf('\n[2] '));  k3 = strfind(rep, sprintf('\n[3] '));
assert(numel(k2) == 1 && numel(k3) == 1, 'ch06a: report sections [2]/[3] not found once');
txt = [rep(1:k2) rep(k3+1:end)];               % the header, then the summary from [3] on
t = regexp(txt, '\n\s+x\s+\(\s*(\d+) pokes\): w (\S+) \(.*?\), l (\S+), e (\S+)\n', 'tokens', 'once');
assert(~isempty(t), 'ch06a: no x summary row');
mm_pin('ch06_x_pokes', str2double(t{1}), 'sig5', t{1});
t = regexp(txt, '\n\s+x\s+\(\s*\d+ pokes\): w (\S+) \(', 'tokens', 'once');
mm_pin('ch06_x_w', str2double(t{1}), 'printed', t{1});
txt = mm_zero('ch06_x_l', txt, '\n\s+x\s+\(\s*\d+ pokes\): w .*?\), l (\S+),', 1e-5);
txt = mm_zero('ch06_x_e', txt, '\n\s+x\s+\(\s*\d+ pokes\): w .*?\), l .*?, e (\S+)\n', 1e-6);
t = regexp(txt, '\n\s+z\s+\(\s*\d+ pokes\): w (\S+) \(', 'tokens', 'once');
mm_pin('ch06_z_w', str2double(t{1}), 'printed', t{1});
t = regexp(txt, '\n\s+grid\s+\(\s*(\d+) pokes\): w (\S+) \((.*?)\)', 'tokens', 'once');
mm_pin('ch06_g_pokes', str2double(t{1}), 'sig5', t{1});
mm_pin('ch06_g_w', str2double(t{2}), 'printed', t{2});
mm_transcript('ch06_compare', strtrim(txt));

fr = dir(fullfile(work, 's6_frames', '*.png'));  fr = {fr.name};
px = fr(~cellfun(@isempty, regexp(fr, '_Seg_2_Tz_')));
pg = fr(~cellfun(@isempty, regexp(fr, '_Seg_2_G5_')));
assert(numel(px) == 1 && numel(pg) == 1, 'ch06a: the Seg 2 Tz / G5 frames not found once');
mm_figure('ch06_poke_x', fullfile(work, 's6_frames', px{1}));
mm_figure('ch06_poke_g', fullfile(work, 's6_frames', pg{1}));
