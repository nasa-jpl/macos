% ch08b_efc.m -- Chapter 8, closing the loop on the dark zone: templates/30_instruments/bench_ctb's
% ctb_efc on a scratch copy of the template (mm_template).  The copy carries no Jacobian (the 37 MB
% .mat is gitignored), so ctb_efc measures it first through ctb_dm_jacobian: every active actuator of
% both DMs poked through the full masked chain (ctb_dm.in, built from the committed ctb_dcr.in and its
% reference spheres -- PLAN 1.2 (zzz)).  Then the EFC loop, iterations capped at 25 (the loop stops
% itself when no regularization improves the measured contrast).  Mask configuration: the shipped
% hard occulter (r_fpm 2.70 lambda/D, Lyot 0.50), N = 512.

[here, check] = mm_template('30_instruments', 'bench_ctb');
work = getenv('MMMAN_WORK');
t0 = tic;
txt = mm_clean(evalc('out = ctb_efc(''niter'', 25, ''outdir'', work);'));
fprintf('ctb_efc wall time %.0f s\n', toc(t0));
check();
disp(txt)

ln = splitlines(string(txt));
keep = startsWith(ln, "[jac] dark zone") | startsWith(ln, "[jac] done") | startsWith(ln, "[efc]");
keep = keep & ~startsWith(ln, "[efc] figure") & ~contains(ln, "absent -- measuring");
% the Jacobian's wall clock moves on every run (13.9, 13.3, 12.6 min so far): masked
ln = regexprep(ln, "in [0-9.]+ min \([0-9.]+ s/poke\)", "in <t> min (<t> s/poke)");
mm_transcript('ch08_efc_run', strjoin(ln(keep), newline));
mm_figure('ch08_efc', fullfile(work, 'ctb_efc.png'));

nit = numel(out.contrast) - 1;
mm_pin('ch08_efc_before', out.c_before, 'sig5', strrep(sprintf('%.2e', out.c_before), 'e-0', 'e-'));
mm_pin('ch08_efc_after',  out.c_after,  'sig5', strrep(sprintf('%.2e', out.c_after),  'e-0', 'e-'));
mm_pin('ch08_efc_gain',   out.c_before/out.c_after, 'sig5', sprintf('%.0f', out.c_before/out.c_after));
mm_pin('ch08_efc_iter',   nit, 'sig5', sprintf('%d', nit));
mm_pin('ch08_efc_s1',     out.stroke_rms_nm(1), 'sig5', sprintf('%.1f nm', out.stroke_rms_nm(1)));
mm_pin('ch08_efc_s2',     out.stroke_rms_nm(2), 'sig5', sprintf('%.1f nm', out.stroke_rms_nm(2)));
nact = regexp(txt, 'G: \d+ px x (\d+) acts', 'tokens', 'once');
mm_pin('ch08_efc_nact',   str2double(nact{1}), 'sig5', nact{1});
