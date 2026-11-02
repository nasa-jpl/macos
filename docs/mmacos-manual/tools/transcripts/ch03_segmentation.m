% ch03_segmentation.m -- Chapter 3 (Segmentation): the e2e flow's stage 3, run UNCHANGED in the
% thread directory (MMMAN_THREAD) on the stage-2 prescription ch02a_imager.m left there.  Stage 3
% segments M1 two ways (a 1-ring pie, 7 segments; a 2-ring hex, 19) through the general runner
% design/runners/run_segmentation.m.  Emits the segmentation knobs, the stage report and the pie
% runner's report, a fresh-session reload of e2e_pie.in (the variant stages 4-6 consume), the
% figures, and the pinned numbers.  Model size 512: a MATLAB of its own.

root = getenv('MMACOS_ROOT');
work = getenv('MMMAN_THREAD');
if ~isfile(fullfile(work, 's2_instrument.in'))
    error('ch03:thread', 'no stage-2 products in %s -- run ch02a_imager first', work);
end
copyfile(fullfile(root, 'templates', '80_end_to_end', 'e2e', 's3_segmentation.m'), work);
addpath(fullfile(root, 'src'));  addpath(fullfile(root, 'design', 'src'));
addpath(fullfile(root, 'design', 'runners'));

pl = splitlines(string(fileread(fullfile(work, 'e2e_params.m'))));
k0 = find(startsWith(pl, "P.seg"), 1);                 % the block, through its closing ');'
k1 = k0 - 1 + find(endsWith(strtrim(regexprep(pl(k0:end), '%.*', '')), ");"), 1);
mm_listing('ch03_seg_params', fullfile(work, 'e2e_params.m'), [k0 k1]);

t0 = tic;
run(fullfile(work, 's3_segmentation.m'));
fprintf('stage 3 wall time %.0f s\n', toc(t0));
work = getenv('MMMAN_THREAD');

mm_listing('ch03_report', fullfile(work, 's3_report.txt'));
mm_listing('ch03_pie_report', fullfile(work, 'e2e_pie_seg_report.txt'));

% ---- the reader's check: reload the pie deck alone ----------------------------------------------
rx = fullfile(work, 'e2e_pie.in');
macos.init(512);
n  = macos.load_rx(rx);
s  = macos.trace(n);
ri = macos.get_ray_info(s.nRays);
np = nnz(logical(ri.ok_pass) & logical(ri.ok_trace));
chk = { ...
 '>> macos.init(512);'
 '>> n = macos.load_rx(''e2e_pie.in'');  s = macos.trace(n);'
 sprintf('   %d elements, %d of %d rays pass', n, np, s.nRays)
 sprintf('   rmsWFE %.4g m', s.rmsWFE)};
mm_transcript('ch03_reload', strjoin(string(chk), newline));

mm_figure('ch03_footprints_pie',  fullfile(work, 's3_footprints_pie.png'));
mm_figure('ch03_views_pie',       fullfile(work, 's3_views_pie.png'));
mm_figure('ch03_footprints_hex2', fullfile(work, 's3_footprints_hex2.png'));

S = load(fullfile(work, 's3_segmentation.mat'), 'S');  S = S.S;
for v = ["pie" "hex2"]
    q = S.(v);
    mm_pin("ch03_" + v + "_nseg",  q.seg.nseg,  'sig5', sprintf('%d', q.seg.nseg));
    mm_pin("ch03_" + v + "_width", q.seg.width, 'sig5', sprintf('%.3g m', q.seg.width));
    mm_pin("ch03_" + v + "_pass",  q.pass_ap,   'sig5', sprintf('%d', q.pass_ap));
    mm_pin("ch03_" + v + "_src",   q.src,       'sig5', sprintf('%d', q.src));
    mm_pin("ch03_" + v + "_clip",  q.pass_bare - q.pass_ap, 'sig5', sprintf('%d', q.pass_bare - q.pass_ap));
end
mm_pin('ch03_reload_nelt', n,  'sig5', sprintf('%d', n));
mm_pin('ch03_reload_pass', np, 'sig5', sprintf('%d', np));
