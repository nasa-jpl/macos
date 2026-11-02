% ch01_telescope.m -- the generated artefacts of Chapter 1 (Telescope), run by
% tools/regen_mmacos_transcripts.py in a scratch directory (MMMAN_WORK).
%
% It runs the e2e flow's stage 1 UNCHANGED: copies templates/80_end_to_end/e2e/{s1_telescope.m,
% e2e_params.m} into the thread directory (MMMAN_THREAD; the later stages' drivers consume what it
% leaves there) and runs the copy (the template writes its products beside itself, so the committed
% template is never touched); then it reloads the emitted prescription in a fresh engine
% session, the way a reader checks an assistant's work, and asks the engine for the first-order
% properties.  Everything the chapter prints is emitted from this run: the knobs (e2e_params.m), the
% design report, the head of the emitted .in, the reload check, the figures, and the pinned numbers.

root = getenv('MMACOS_ROOT');
work = getenv('MMMAN_THREAD');
tpl  = fullfile(root, 'templates', '80_end_to_end', 'e2e');
copyfile(fullfile(tpl, 's1_telescope.m'), work);
copyfile(fullfile(tpl, 'e2e_params.m'), work);
addpath(fullfile(root, 'src'));  addpath(fullfile(root, 'design', 'src'));

% ---- the knobs: e2e_params.m up to the instrument block (the telescope's own parameters) --------
pl = splitlines(string(fileread(fullfile(work, 'e2e_params.m'))));
k_inst = find(startsWith(pl, "P.inst"), 1);
mm_listing('ch01_params', fullfile(work, 'e2e_params.m'), [1 k_inst-1]);

% ---- the run ------------------------------------------------------------------------------------
t0 = tic;
run(fullfile(work, 's1_telescope.m'));
fprintf('stage 1 wall time %.0f s\n', toc(t0));       % in the log only, never printed
work = getenv('MMMAN_THREAD');                        % the stage script shares this workspace

mm_listing('ch01_report', fullfile(work, 's1_report.txt'));

% ---- the emitted prescription: the header and the first element ---------------------------------
rx = fullfile(work, 's1_telescope.in');
rl = splitlines(string(fileread(rx)));
k2 = find(startsWith(strtrim(rl), "iElt="), 2);       % up to the second element
mm_transcript('ch01_rx_head', strjoin(rl(1:k2(end)-1), newline));

% ---- the reader's check: reload the .in alone and ask the engine --------------------------------
macos.init(256);
n  = macos.load_rx(rx);
s  = macos.trace(n);
ri = macos.get_ray_info(s.nRays);
np = nnz(logical(ri.ok_pass) & logical(ri.ok_trace));
fo = macos.first_order_properties();
chk = { ...
 '>> macos.init(256);'
 '>> n = macos.load_rx(''s1_telescope.in'');  s = macos.trace(n);'
 sprintf('   %d elements, %d of %d rays pass', n, np, s.nRays)
 '>> p = macos.first_order_properties()'
 sprintf('   efl_baseunits  %.5g     (m: the deck''s BaseUnits)', fo.efl_baseunits)
 sprintf('   fno            %.5g', fo.fno)
 sprintf('   dpup_m         %.5g', fo.dpup_m)
 sprintf('   obscuration    %.4g', fo.obscuration)
 sprintf('   lamD_arcsec    %.5g', fo.lamD_arcsec)};
mm_transcript('ch01_reload', strjoin(string(chk), newline));

% ---- figures (the template's own files, unmodified) ---------------------------------------------
mm_figure('ch01_views',     fullfile(work, 's1_views.png'));
mm_figure('ch01_wfe_field', fullfile(work, 's1_wfe_field.png'));
mm_figure('ch01_fpmap',     fullfile(work, 's1_fpmap.png'));

% ---- the numbers the prose quotes ---------------------------------------------------------------
S = load(fullfile(work, 's1_telescope.mat'), 'rpt', 'bias', 'r_hole');
mm_pin('ch01_efl_m',    fo.efl_baseunits, 'sig5', sprintf('%.2f m', fo.efl_baseunits));
mm_pin('ch01_fno',      fo.fno,           'sig5', sprintf('f/%.2f', fo.fno));
mm_pin('ch01_lamD_mas', fo.lamD_arcsec*1e3, 'sig5', sprintf('%.2f mas', fo.lamD_arcsec*1e3));
mm_pin('ch01_bias_arcmin', S.bias,        'sig5', sprintf('%g''', S.bias));
mm_pin('ch01_hole_r_m', S.r_hole,         'sig5', sprintf('%.3f m', S.r_hole));
mm_pin('ch01_rays_pass', np,              'sig5', sprintf('%d', np));
mm_pin('ch01_rays',     s.nRays,          'sig5', sprintf('%d', s.nRays));
