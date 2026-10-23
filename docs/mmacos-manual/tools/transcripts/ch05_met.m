% ch05_met.m -- Chapter 5 (Laser metrology): the e2e flow's stage 5 in the thread directory
% (MMMAN_THREAD), on the e2e_pie.in + e2e_pieHx.m of ch03 and the s4_jacobians.mat of ch04.  Stage 5
% (design/runners/run_met.m) builds the as-built Stewart truss (6 launchers per segment, 6 fiducials
% on the M2 rim, an aft ring on the fold), checks dldx engine-FD against analytic, scores sensing
% (edge sensors, MET, both) as dimensionless WEMs, optimizes the launcher layout per segment shape
% class, and runs a Monte-Carlo acceptance (seeded: rng(1)).
%
% ONE CUT from the template (asserted to match once, and said in the chapter): the script's last
% block copies the optimized layout into the SHARED preset library (mmacos/design/runners/presets/),
% a write into MACOS_resources the manual must not make.  Everything else runs as committed.
%
% Pinned as BOUNDS (numerical zeros), one line each:
%   ch05_dldx_rel  < 1e-5   dldx engine-FD vs analytic, max relative difference (agreement = FD noise)
% Everything else is a signal (sig5 or printed).

root = getenv('MMACOS_ROOT');
work = getenv('MMMAN_THREAD');
for f = ["e2e_pie.in" "e2e_pieHx.m" "s4_jacobians.mat"]
    if ~isfile(fullfile(work, f))
        error('ch05:thread', 'no %s in %s -- run ch03_segmentation and ch04_jacobians first', f, work);
    end
end
addpath(fullfile(root, 'src'));  addpath(fullfile(root, 'design', 'src'));
addpath(fullfile(root, 'design', 'runners'));

src = fileread(fullfile(root, 'templates', '80_end_to_end', 'e2e', 's5_met.m'));
cut = regexp(src, '% promote the optimized layout to the shared preset library.*?\nend\n', 'match');
assert(numel(cut) == 1, 'ch05: the preset-promotion block did not match exactly once');
fid = fopen(fullfile(work, 's5_met.m'), 'w');  fwrite(fid, strrep(src, cut{1}, ''));  fclose(fid);

pl = splitlines(string(fileread(fullfile(work, 'e2e_params.m'))));
k0 = find(startsWith(pl, "P.met"), 1);                 % the block, through its closing ');'
k1 = k0 - 1 + find(endsWith(strtrim(regexprep(pl(k0:end), '%.*', '')), ");"), 1);
mm_listing('ch05_met_params', fullfile(work, 'e2e_params.m'), [k0 k1]);

t0 = tic;
run(fullfile(work, 's5_met.m'));
fprintf('stage 5 wall time %.0f s\n', toc(t0));
work = getenv('MMMAN_THREAD');

rep = fileread(fullfile(work, 's5_report.txt'));
rep = mm_zero('ch05_dldx_rel', rep, 'engine-FD vs analytic: (\S+) max rel', 1e-5);
mm_transcript('ch05_report', strtrim(rep));

for f = ["met_layout" "metopt_layout" "met_metric" "met_view_rx"]
    mm_figure("ch05_" + f, fullfile(work, "e2e_pie_" + f + ".png"));
end

% ---- the numbers the prose quotes, as printed ---------------------------------------------------
t = regexp(rep, 'as-built truss: (\d+) beams', 'tokens', 'once');
mm_pin('ch05_beams', str2double(t{1}), 'sig5', t{1});
t = regexp(rep, 'dedx (\d+)x(\d+)', 'tokens', 'once');
mm_pin('ch05_edge_rows', str2double(t{1}), 'sig5', t{1});
mm_pin('ch05_dofs', str2double(t{2}), 'sig5', t{2});
t = regexp(rep, 'min margin (\S+) \(', 'tokens', 'once');
mm_pin('ch05_hole_margin', str2double(t{1}), 'printed', t{1});
for c = ["as-built" "optimized"]
    t = regexp(rep, "\n\s+" + c + "\s+\|\s+(\S+)\s+\S+\s+(\S+)\s+\|\s+(\S+)\s+\S+\s+(\S+)", 'tokens', 'once');
    assert(~isempty(t), 'ch05: no WEM row for %s', c);
    k = strrep(c, '-', '');
    mm_pin("ch05_" + k + "_met",  str2double(t{1}), 'printed', t{1});
    mm_pin("ch05_" + k + "_emet", str2double(t{3}), 'printed', t{3});
end
t = regexp(rep, 'Monte-Carlo \((\d+) draws\): (\S+) nm rms\s+\(analytic (\S+)\)', 'tokens', 'once');
mm_pin('ch05_mc_draws', str2double(t{1}), 'sig5', t{1});
mm_pin('ch05_mc_nm', str2double(t{2}), 'printed', t{2});
mm_pin('ch05_mc_analytic_nm', str2double(t{3}), 'printed', t{3});
