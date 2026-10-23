% ch06b_simulate.m -- Chapter 6 (Closing the loop), part 2: the e2e flow's stage 7 in the thread
% directory (MMMAN_THREAD).  Stage 7 (design/runners/run_simulator.m) plays a 500 s disturbance
% history through the engine, uncorrected and corrected, twice: (A) a metrology-bias drift the loop
% cannot see, (B) a focus/astigmatism figure drift the rigid-body loop cannot sense; image-based
% wavefront control at 20 s, the metrology loop holding pose, a wavefront-maintenance recontrol at
% 400 s.  Seeded (rng(7)); run UNCHANGED (it writes only beside itself).  The chapter prints each
% scenario's setup ([0]-[3b]) and its final block ([5]), not the 51-row time table, and pins the
% corrected wavefront just before the recontrol (390 s), at it (400 s) and at the end (500 s).
% All pinned numbers are signals (printed): engine wavefronts of nanometres.

root = getenv('MMACOS_ROOT');
work = getenv('MMMAN_THREAD');
for f = ["e2e_pie.in" "e2e_pieHx.m" "s4_jacobians.mat" "e2e_pie_met.mat"]
    if ~isfile(fullfile(work, f))
        error('ch06b:thread', 'no %s in %s -- run ch03-ch05 first', f, work);
    end
end
copyfile(fullfile(root, 'templates', '80_end_to_end', 'e2e', 's7_simulate.m'), work);
addpath(fullfile(root, 'src'));  addpath(fullfile(root, 'design', 'src'));
addpath(fullfile(root, 'design', 'runners'));

t0 = tic;
run(fullfile(work, 's7_simulate.m'));
fprintf('stage 7 wall time %.0f s\n', toc(t0));
work = getenv('MMMAN_THREAD');

for s = ["A" "B"]
    rep = fileread(fullfile(work, "s7" + s + "_report.txt"));
    k4 = strfind(rep, sprintf('\n[4] '));  k5 = strfind(rep, sprintf('\n[5] '));
    assert(numel(k4) == 1 && numel(k5) == 1, 'ch06b: s7%s report sections [4]/[5] not found once', s);
    mm_transcript("ch06_sim" + s, strtrim([rep(1:k4) rep(k5+1:end)]));
    for tt = [390 400 500]
        r = regexp(rep, sprintf('\\n\\s+%d\\s+\\S+\\s+(\\S+)\\s', tt), 'tokens', 'once');
        assert(~isempty(r), 'ch06b: no t=%d row in s7%s', tt, s);
        mm_pin(sprintf('ch06_%s_corr_%d', lower(s), tt), str2double(r{1}), 'printed', r{1});
    end
    mm_figure("ch06_sim" + s + "_t390", fullfile(work, "s7" + s + "_sim_frames", "t039.png"));
    mm_figure("ch06_sim" + s + "_t500", fullfile(work, "s7" + s + "_sim_frames", "t050.png"));
end
