% ch08e_telescope.m -- Chapter 8, the coronagraph on the telescope: the e2e thread's stages 9, 10 and
% 11 (templates/80_end_to_end/e2e: s9_diffraction.m the coupled telescope + CTB diffraction model,
% s10_sensing.m pairwise probing on it, s11_cycles.m recompute-and-restart EFC cycles on it).
% Unlike the other drivers this one runs no stage: the three stages cost about eighteen hours
% (thirteen control matrices of 3520 chain runs each at model 512), so the chapter reads their committed
% products (s9_run.mat, s10_run.mat, s11_run.mat and the stages' own figures, PLAN_CONSOLIDATION 3a,
% CCL 2026-10-30) from the library's e2e directory and errors if they are not there -- the page is
% then not regenerated (the regenerator writes nothing).  It does run ctb_e2e_figures (bench_ctb),
% two minutes, for the instrument's views and the station walk (2026-11-02).  To refresh it, run the three stages in order
% (s9 -> s10 -> s11, from the e2e README) and regenerate this driver.  The numbers quoted in the prose
% are pinned from the .mat files; the listings are built from them (no dates, no paths).
% The library is read, never written.

root = getenv('MMACOS_ROOT');
e2e  = fullfile(root, 'templates', '80_end_to_end', 'e2e');
for f = {'s9_run.mat', 's10_run.mat', 's11_run.mat', 's9_fields.png', 's10_probe.png', 's11_cycles.png'}
    if ~isfile(fullfile(e2e, f{1}))
        error('ch08e:products', 'no %s in %s -- run the e2e stages 9-11 first (their README)', f{1}, e2e);
    end
end
addpath(e2e);
P  = e2e_params();  lam_nm = P.lambda_m * 1e9;
S9 = load(fullfile(e2e, 's9_run.mat'));
S10 = load(fullfile(e2e, 's10_run.mat'));
S11 = load(fullfile(e2e, 's11_run.mat'));
say = @(varargin) sprintf(varargin{:});
e2 = @(v) regexprep(sprintf('%.1e', v), 'e-0?', 'e-');
e3 = @(v) regexprep(sprintf('%.2e', v), 'e-0?', 'e-');
rad = @(nm) 2*pi*2*nm/lam_nm;                     % surface stroke nm -> wavefront phase, rad

% ---- stage 9: the coupled model -------------------------------------------------------------------
c = S9.S.coupled;  b = S9.S.bench;
L = strings(0);
L(end+1) = say('[s9] coupled model: static dark zone %.3e with the bench''s plain Lyot disk (bench alone %.3e)', c.contrast, b.contrast);
L(end+1) = say('[s9] the secondary''s shadow is %.2f %% of the pupil area; the same deck without it: %.3e', 100*c.shadow_frac, c.noobs);
L(end+1) = say('[s9] central Lyot blocker, pupil radii %s -> %s', mat2str(c.lyot_inner), strjoin(compose('%.3e', c.lyot_scan), ' '));
L(end+1) = say('[s9] the pick %.2f: static %.3e; flattening DM1 (%.4f -> %.4f waves rms on the pupil, %.2f nm rms on DM1): %.3e', ...
    c.lyot_pick, c.contrast_lyot, c.flat.opd_before_waves, c.flat.opd_after_waves, c.flat.stroke_rms_nm, c.flat.contrast);
L(end+1) = say('[s9] model limit: pupil axis ratio %.4f, PSF core axis ratio %.4f; static %.3e .. %.3e over the occulter bracket', ...
    c.pupil_ratio, c.core_ratio, c.bracket(1), c.bracket(2));
L(end+1) = say('[s9] EFC, control matrix at the flat mirrors: %.3e -> %.3e in %d iterations, strokes %.2f / %.2f nm rms', ...
    S9.efc_contrast(1), S9.efc_contrast(end), numel(S9.efc_contrast)-1, S9.strokes_nm(1), S9.strokes_nm(2));
L(end+1) = say('[s9] re-measured about that state: %.3e -> %.3e in %d more, strokes %.2f / %.2f nm rms', ...
    S9.efc_contrast_r1(1), S9.efc_contrast_r1(end), numel(S9.efc_contrast_r1)-1, S9.strokes_nm_r1(1), S9.strokes_nm_r1(2));

% ---- stage 10: pairwise probing on it -------------------------------------------------------------
d = S10.default;  s = S10.static;  g = S10.dug;  k = S10.control;
L(end+1) = say('[s10] the default amplitude rule: %.1f nm peak (%.2f rad), estimate %.3e against the field''s %.3e, floor %.2e, relative error %.3f, ratios %s', ...
    d.probe_nm, rad(d.probe_nm), d.c_est, d.c_true, d.floor, d.rel_err, strjoin(compose('%.3f', d.ratio), ' '));
for a = S10.amp
    L(end+1) = say('[s10] amplitude %.1fx = %6.2f nm: noiseless floor %.2e, relative error %.2e; at 1e10 photons floor %.2e, error %.2e', ...
        a.rel, a.nm, a.floor, a.err, a.floor_ph, a.err_ph);
end
L(end+1) = say('[s10] the pick %.1fx = %.1f nm (%.2f rad): relative error %.3f, floor %.2e', s.probe_rel, s.probe_nm, rad(s.probe_nm), s.rel_err, s.floor);
for p = S10.photons
    L(end+1) = say('[s10] photons %-6s: floor %.2e, relative error %.2e', ternary_(p.ph == 0, 'none', sprintf('%.0e', p.ph)), p.floor, p.rel);
end
L(end+1) = say('[s10] the dug state %.3e with the matrix measured about the flat state: relative error %.3f; at that matrix''s own point (%.3e): %.4f', ...
    g.c_true, g.rel_err, k.c_true, k.rel_err);
L(end+1) = say('[s10] EFC on the estimate, matrix fixed: %s', strjoin(compose('%.3e', S10.efc), ' '));
L(end+1) = say('[s10]   perfect sensing (stage 9):   %s', strjoin(compose('%.3e', S9.efc_contrast), ' '));
L(end+1) = say('[s10] re-measured on the estimate:  %s', strjoin(compose('%.3e', S10.efc_r1), ' '));
L(end+1) = say('[s10]   perfect sensing (stage 9):   %s', strjoin(compose('%.3e', S9.efc_contrast_r1), ' '));

% ---- stage 11: recompute-and-restart cycles -------------------------------------------------------
R = S11.R;  cfg = S11.cycles_cfg;
L(end+1) = say('[s11] cycles to %d, stop below x%.2f per cycle, goal %.1e; probe %.1fx the default', cfg.max, cfg.stop_gain, cfg.goal, S11.pick_rel);
L(end+1) = say('[s11] perfect sensing: %d cycles, %d iterations, %.3e -> %.3e, stop ''%s''', numel(R.perf.cycles), R.perf.total_iters, ...
    R.perf.cycles(1).c_start, R.perf.c_final, R.perf.stop_reason);
L(end+1) = say('[s11] probed:          %d cycles, %d iterations, %.3e -> %.3e, stop ''%s''', numel(R.sens.cycles), R.sens.total_iters, ...
    R.sens.cycles(1).c_start, R.sens.c_final, R.sens.stop_reason);
mm_transcript('ch08_tel_run', strjoin(L, newline));

C = strings(0);
C(end+1) = 'chain    cycle  matrix        start      end       iterations  gain   strokes nm rms   estimate error';
for cy = R.perf.cycles
    C(end+1) = say('perfect  %d      %-9s   %.3e  %.3e  %2d        x%5.2f  %5.2f / %5.2f', cy.k, ternary_(cy.jac_reused, 'reused', 'measured'), ...
        cy.c_start, cy.c_end, cy.iters, cy.gain, cy.stroke_rms_nm(1), cy.stroke_rms_nm(2));
end
for cy = R.sens.cycles
    C(end+1) = say('probed   %d      %-9s   %.3e  %.3e  %2d        x%5.2f  %5.2f / %5.2f    %.3f', cy.k, ternary_(cy.jac_reused, 'reused', 'measured'), ...
        cy.c_start, cy.c_end, cy.iters, cy.gain, cy.stroke_rms_nm(1), cy.stroke_rms_nm(2), cy.est_rel0);
end
mm_transcript('ch08_tel_cycles', strjoin(C, newline));

% ---- figures: the stages' own, and the two ctb_e2e_figures draws from their products --------------------
% (the folded instrument's views from the stage-8 deck, the field at eight planes of the stage-9 model;
%  about two minutes: one view_rx per panel and one masked chain pass at model 512)
addpath(fullfile(root, 'templates', '30_instruments', 'bench_ctb'));
figdir = tempname;  mkdir(figdir);
FX = ctb_e2e_figures('outdir', figdir);
mm_figure('ch08_tel_views',  FX.fig_views);
mm_figure('ch08_tel_walk',   FX.fig_walk);
mm_figure('ch08_tel_fields', fullfile(e2e, 's9_fields.png'));
mm_figure('ch08_tel_probe',  fullfile(e2e, 's10_probe.png'));
mm_figure('ch08_tel_cycles', fullfile(e2e, 's11_cycles.png'));

% ---- pins ---------------------------------------------------------------------------------------------
jt = [R.perf.cycles.jac_time_s, R.sens.cycles.jac_time_s];  jt = jt(jt > 60) / 60;
mm_pin('ch08_tel_static',    c.contrast,      'sig5', e2(c.contrast));
mm_pin('ch08_tel_bench',     b.contrast,      'sig5', e2(b.contrast));
mm_pin('ch08_tel_noobs',     c.noobs,         'sig5', e2(c.noobs));
mm_pin('ch08_tel_shadow',    100*c.shadow_frac, 'sig5', sprintf('%.1f %%', 100*c.shadow_frac));
mm_pin('ch08_tel_blocker',   c.lyot_pick,     'sig5', sprintf('%.2f', c.lyot_pick));
mm_pin('ch08_tel_lyot',      c.contrast_lyot, 'sig5', e2(c.contrast_lyot));
mm_pin('ch08_tel_thr',       100*(1 - (c.lyot_pick/0.5)^2), 'sig5', sprintf('%.0f %%', 100*(1 - (c.lyot_pick/0.5)^2)));
mm_pin('ch08_tel_wfe',       c.flat.opd_before_waves, 'sig5', sprintf('%.4f', c.flat.opd_before_waves));
mm_pin('ch08_tel_anam',      100*(c.pupil_ratio - 1), 'sig5', sprintf('%.1f %%', 100*(c.pupil_ratio - 1)));
mm_pin('ch08_tel_efc0',      S9.efc_contrast(end), 'sig5', e2(S9.efc_contrast(end)));
mm_pin('ch08_tel_efc0_it',   numel(S9.efc_contrast)-1, 'sig5', sprintf('%d', numel(S9.efc_contrast)-1));
mm_pin('ch08_tel_efc1',      S9.efc_contrast_r1(end), 'sig5', e2(S9.efc_contrast_r1(end)));
mm_pin('ch08_tel_efc1_it',   numel(S9.efc_contrast_r1)-1, 'sig5', sprintf('%d', numel(S9.efc_contrast_r1)-1));
mm_pin('ch08_tel_gain9',     c.contrast_lyot / S9.efc_contrast_r1(end), 'sig5', sprintf('%.0f', c.contrast_lyot / S9.efc_contrast_r1(end)));
mm_pin('ch08_tel_def_nm',    d.probe_nm,      'sig5', sprintf('%.0f nm', d.probe_nm));
mm_pin('ch08_tel_def_rad',   rad(d.probe_nm), 'sig5', sprintf('%.2f rad', rad(d.probe_nm)));
mm_pin('ch08_tel_def_rel',   d.rel_err,       'sig5', sprintf('%.0f %%', 100*d.rel_err));
mm_pin('ch08_tel_pick_nm',   s.probe_nm,      'sig5', sprintf('%.0f nm', s.probe_nm));
mm_pin('ch08_tel_pick_rel',  s.rel_err,       'sig5', sprintf('%.1f %%', 100*s.rel_err));
mm_pin('ch08_tel_lin_rel',   S10.amp(1).err,  'sig5', sprintf('%.1f %%', 100*S10.amp(1).err));
mm_pin('ch08_tel_ratio',     max(abs(s.ratio - 1)), 'sig5', sprintf('%.2f', max(abs(s.ratio - 1))));
mm_pin('ch08_tel_ph10_rel',  S10.photons(2).rel, 'sig5', sprintf('%.0f %%', 100*S10.photons(2).rel));
mm_pin('ch08_tel_ph8_floor', S10.photons(3).floor, 'sig5', e2(S10.photons(3).floor));
mm_pin('ch08_tel_dug_rel',   g.rel_err,       'sig5', sprintf('%.0f %%', 100*g.rel_err));
mm_pin('ch08_tel_ctl_rel',   k.rel_err,       'sig5', sprintf('%.1f %%', 100*k.rel_err));
mm_pin('ch08_tel_sefc0',     S10.efc(end),    'sig5', e2(S10.efc(end)));
mm_pin('ch08_tel_sefc1',     S10.efc_r1(end), 'sig5', e2(S10.efc_r1(end)));
mm_pin('ch08_tel_goal',      cfg.goal,        'sig5', e2(cfg.goal));
mm_pin('ch08_tel_jac_min',   mean(jt),        'sig5', sprintf('%.0f', mean(jt)));
mm_pin('ch08_tel_perf_nc',   numel(R.perf.cycles), 'sig5', sprintf('%d', numel(R.perf.cycles)));
mm_pin('ch08_tel_perf_it',   R.perf.total_iters, 'sig5', sprintf('%d', R.perf.total_iters));
mm_pin('ch08_tel_perf_end',  R.perf.c_final,  'sig5', e2(R.perf.c_final));
mm_pin('ch08_tel_perf_last', R.perf.cycles(end).gain, 'sig5', sprintf('%.0f', R.perf.cycles(end).gain));
mm_pin('ch08_tel_sens_nc',   numel(R.sens.cycles), 'sig5', sprintf('%d', numel(R.sens.cycles)));
mm_pin('ch08_tel_sens_it',   R.sens.total_iters, 'sig5', sprintf('%d', R.sens.total_iters));
mm_pin('ch08_tel_sens_end',  R.sens.c_final,  'sig5', e2(R.sens.c_final));
mm_pin('ch08_tel_sens_c3',   R.sens.cycles(end-1).c_end, 'sig5', e2(R.sens.cycles(end-1).c_end));
mm_pin('ch08_tel_est_last',  R.sens.cycles(end).est_rel0, 'sig5', sprintf('%.2f %%', 100*R.sens.cycles(end).est_rel0));
mm_pin('ch08_tel_est_max3',  max([R.sens.cycles(1:end-1).est_rel0]), 'sig5', sprintf('%.1f %%', 100*max([R.sens.cycles(1:end-1).est_rel0])));
mm_pin('ch08_tel_stroke',    max([R.perf.cycles.stroke_rms_nm, R.sens.cycles.stroke_rms_nm]), 'sig5', ...
    sprintf('%.0f nm', max([R.perf.cycles.stroke_rms_nm, R.sens.cycles.stroke_rms_nm])));

function out = ternary_(c, a, b)
    if c, out = a; else, out = b; end
end
