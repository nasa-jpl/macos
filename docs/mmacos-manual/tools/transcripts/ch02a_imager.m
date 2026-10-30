% ch02a_imager.m -- Chapter 2, the thread's instrument: the e2e flow's stage 2 (the Offner ring-field
% relay behind the stage-1 telescope), run UNCHANGED in the thread directory (MMMAN_THREAD) on the
% products ch01_telescope.m left there.  Emits the instrument knobs, the stage-2 design report, a
% fresh-session reload of the emitted .in with the engine's first-order properties, the figures, pins.

root = getenv('MMACOS_ROOT');
work = getenv('MMMAN_THREAD');
if ~isfile(fullfile(work, 's1_telescope.mat'))
    error('ch02a:thread', 'no stage-1 products in %s -- run ch01_telescope first', work);
end
copyfile(fullfile(root, 'templates', '80_end_to_end', 'e2e', 's2_instrument.m'), work);
addpath(fullfile(root, 'src'));  addpath(fullfile(root, 'design', 'src'));

pl = splitlines(string(fileread(fullfile(work, 'e2e_params.m'))));
k0 = find(startsWith(pl, "P.inst"), 1);                % the block, through its closing ');'
k1 = k0 - 1 + find(endsWith(strtrim(regexprep(pl(k0:end), '%.*', '')), ");"), 1);
mm_listing('ch02_inst_params', fullfile(work, 'e2e_params.m'), [k0 k1]);

t0 = tic;
run(fullfile(work, 's2_instrument.m'));
fprintf('stage 2 wall time %.0f s\n', toc(t0));
work = getenv('MMMAN_THREAD');

mm_listing('ch02_report', fullfile(work, 's2_report.txt'));

rx = fullfile(work, 's2_instrument.in');
macos.init(256);
n  = macos.load_rx(rx);
s  = macos.trace(n);
ri = macos.get_ray_info(s.nRays);
np = nnz(logical(ri.ok_pass) & logical(ri.ok_trace));
fo = macos.first_order_properties();
chk = { ...
 '>> macos.init(256);'
 '>> n = macos.load_rx(''s2_instrument.in'');  s = macos.trace(n);'
 sprintf('   %d elements, %d of %d rays pass', n, np, s.nRays)
 '>> p = macos.first_order_properties()'
 sprintf('   efl_baseunits  %.5g', fo.efl_baseunits)
 sprintf('   fno            %.5g', fo.fno)};
mm_transcript('ch02_reload', strjoin(string(chk), newline));

mm_figure('ch02_views',     fullfile(work, 's2_views.png'));
mm_figure('ch02_wfe_field', fullfile(work, 's2_wfe_field.png'));

mm_pin('ch02_nelt',  n,               'sig5', sprintf('%d', n));
mm_pin('ch02_efl_m', fo.efl_baseunits, 'sig5', sprintf('%.2f m', fo.efl_baseunits));
mm_pin('ch02_fno',   fo.fno,          'sig5', sprintf('f/%.2f', fo.fno));
mm_pin('ch02_rays_pass', np,          'sig5', sprintf('%d', np));
mm_pin('ch02_rays',  s.nRays,         'sig5', sprintf('%d', s.nRays));
