% ch02d_imagedemo.m -- Chapter 2's counter-example: ImageDemo.in (the MACOS Manual's Appendix A
% example, after the Mars Observer Camera), a SMALL DEDICATED IMAGER -- telescope and instrument one
% design.  Loaded as distributed (with its Zernike terms; PLAN_MMACOS_MANUAL: the camera template
% strips them first), traced, and asked for its first-order properties.

root = getenv('MMACOS_ROOT');
addpath(fullfile(root, 'src'));
here = fileparts(fileparts(fileparts(mfilename('fullpath'))));       % docs/mmacos-manual
rx = fullfile(here, '..', 'macos-manual', 'examples', 'ImageDemo.in');
macos.init(256);
n  = macos.load_rx(rx);
s  = macos.trace(n);
ri = macos.get_ray_info(s.nRays);
np = nnz(logical(ri.ok_pass) & logical(ri.ok_trace));
fo = macos.first_order_properties();
mm_transcript('ch02_imagedemo', strjoin(string({ ...
 '>> macos.init(256);  n = macos.load_rx(''ImageDemo.in'');  s = macos.trace(n);'
 sprintf('   %d elements, %d of %d rays pass', n, np, s.nRays)
 '>> p = macos.first_order_properties()'
 sprintf('   efl_baseunits  %.5g', fo.efl_baseunits)
 sprintf('   fno            %.5g', fo.fno)
 sprintf('   dpup_m         %.5g', fo.dpup_m)
 sprintf('   obscuration    %.4g', fo.obscuration)}), newline));
mm_pin('ch02_imagedemo_nelt', n,       'sig5', sprintf('%d', n));
mm_pin('ch02_imagedemo_fno',  fo.fno,  'sig5', sprintf('f/%.3g', fo.fno));
mm_pin('ch02_imagedemo_obsc', fo.obscuration, 'sig5', sprintf('%.2f', fo.obscuration));
