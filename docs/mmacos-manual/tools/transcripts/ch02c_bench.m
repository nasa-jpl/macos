% ch02c_bench.m -- Chapter 2, a bench on its own: templates/40_benches/bench_layout, the generic
% folded test-arm bench built optic by optic with macos.design.Bench and optimized in three staged
% solves.  The script is copied to the scratch directory and the copy is run (it writes its
% products beside itself), in a MATLAB of its own.

root = getenv('MMACOS_ROOT');
work = getenv('MMMAN_WORK');
copyfile(fullfile(root, 'templates', '40_benches', 'bench_layout', 'example_bench_layout.m'), work);
addpath(fullfile(root, 'src'));  addpath(fullfile(root, 'design', 'src'));

txt = evalc('run(fullfile(work, ''example_bench_layout.m''));');
work = getenv('MMMAN_WORK');
txt = mm_clean(txt);
txt = mm_zero('ch02_bench_chief_agree', txt, 'chief agreement: max (\S+) mm', 1e-9);
txt = mm_zero('ch02_bench_spread', txt, 'ray spread \S+ -> (\S+) mrad', 1e-6);
txt = mm_zero('ch02_bench_shift', txt, 'beam shift \S+ -> (\S+) um', 1e-5);
mm_transcript('ch02_bench_run', strtrim(txt));
mm_figure('ch02_bench_view', fullfile(work, 'bench_layout_view_std.png'));
mm_figure('ch02_bench_params', fullfile(work, 'bench_layout_params.png'));
rl = splitlines(string(fileread(fullfile(work, 'bench_layout_opt.in'))));
mm_pin('ch02_bench_nelt', str2double(regexp(strjoin(rl(startsWith(strtrim(rl), "nElt=")), ''), '\d+', 'match', 'once')), 'sig5');
