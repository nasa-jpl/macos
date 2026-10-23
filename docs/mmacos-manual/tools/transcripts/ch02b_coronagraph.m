% ch02b_coronagraph.m -- Chapter 2, a coronagraph: templates/30_instruments/coro_planet_demo, the
% CoroExample recipe (the CLI's CoroExample.jou "pgp 4" pages: entrance pupil, the star occulted at the
% focal-plane mask, the Lyot pupil, the star at the detector).  Run UNCHANGED with its figures sent to
% the scratch directory (opts.outdir), in a MATLAB of its own (model size 512).  The chapter uses the
% surfaces page only: the demo's planet-injection page does not show a planet PSF on 2026-10-23
% (OPEN_QUESTIONS.md), so its figure and numbers are not printed.

root = getenv('MMACOS_ROOT'); work = getenv('MMMAN_WORK');
addpath(fullfile(root, 'src'));
addpath(fullfile(root, 'templates', '30_instruments', 'coro_planet_demo'));
txt = mm_clean(evalc('out = coro_planet_demo(''outdir'', work);'));
disp(txt); disp(out)                                   % the whole run, in the log
k = regexp(txt, '[^\n]*wrote coro_surfaces\.png', 'match', 'once');
mm_transcript('ch02_coro_run', strjoin({'>> out = coro_planet_demo();', strtrim(k)}, newline));
mm_figure('ch02_coro_surfaces', fullfile(work, 'coro_surfaces.png'));
