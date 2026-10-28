% ch08c_pol.m -- Chapter 8, what polarization costs the dark zone: templates/30_instruments/bench_ctb's
% ctb_efc_physics on a scratch copy of the template (mm_template), on its default chain (scalar vortex,
% charge 4, Lyot 0.60, no apodizer, N = 512, 3-15 lambda/D, one wavelength).  Two runs:
%   1. 'pol' false -- EFC on the scalar chain; measures the vortex chain's Jacobian (tag 'mono').
%   2. 'pol' true  -- the same loop with the coated train's Jones pupil (MgF2 over Al on all ten
%      reflectors, macos.jones_pupil at the exit pupil) applied as four per-component pupil screens,
%      unpolarized light; the Jacobian of run 1 is reused (the screens differ from identity by <4e-3).
% The printed 'pol floor' is the part one DM setting cannot null: the spread of the polarization
% components about their common field.  The figure is this driver's own: both contrast histories.
% Deck: ctb_dm.in from the committed ctb_dcr.in and its reference spheres (PLAN 1.2 (zzz)).

[here, check] = mm_template('30_instruments', 'bench_ctb');
work = getenv('MMMAN_WORK');
t0 = tic;
t1 = mm_clean(evalc('o1 = ctb_efc_physics(''pol'', false, ''tag'', ''mono'', ''outdir'', work);'));
jp = o1.jac;
t2 = mm_clean(evalc('o2 = ctb_efc_physics(''pol'', true, ''tag'', ''pol'', ''jac'', jp, ''outdir'', work);'));
fprintf('ctb_efc_physics wall time %.0f s\n', toc(t0));
check();
disp(t1); disp(t2)

sel = @(t) strjoin(string(regexp(t, '\[phys\] (iter 0|final|computing|iter \d+: no)[^\n]*', 'match')), newline);
mm_transcript('ch08_pol_run', strjoin(["% scalar vortex, no coatings", sel(t1), "", ...
                                      "% the same, coated train, unpolarized light", sel(t2)], newline));

f = figure('Visible', 'off', 'Position', [100 100 640 420]);
semilogy(0:numel(o1.contrast)-1, o1.contrast, 'o-', 'LineWidth', 1.5);  hold on
semilogy(0:numel(o2.contrast)-1, o2.contrast, 's-', 'LineWidth', 1.5);
yline(o2.pol_floor, 'k--', 'polarization floor', 'LabelHorizontalAlignment', 'left');
ylim([1e-16 1e-7]);  grid on;  xlabel('EFC iteration');  ylabel('dark-zone mean contrast, 3-15 \lambda/D');
legend('scalar vortex, no coatings', 'coated train, unpolarized', 'Location', 'northeast');
title('CTB vortex (charge 4, Lyot 0.60): EFC with and without polarization');
print(f, fullfile(work, 'ctb_pol_efc.png'), '-dpng', '-r150');  close(f);
mm_figure('ch08_pol', fullfile(work, 'ctb_pol_efc.png'));

e = @(v) regexprep(sprintf('%.1e', v), 'e-0?', 'e-');
mm_pin('ch08_vx_before', o1.c_before, 'sig5', e(o1.c_before));
mm_pin('ch08_vx_after',  o1.c_after,  'sig5', e(o1.c_after));
mm_pin('ch08_pol_before', o2.c_before, 'sig5', e(o2.c_before));
mm_pin('ch08_pol_after',  o2.c_after,  'sig5', e(o2.c_after));
mm_pin('ch08_pol_floor',  o2.pol_floor, 'sig5', e(o2.pol_floor));
