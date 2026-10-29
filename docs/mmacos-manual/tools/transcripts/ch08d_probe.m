% ch08d_probe.m -- Chapter 8, sensing the dark zone with pairwise probes: templates/30_instruments/
% bench_ctb's ctb_probe_estimate and ctb_efc('sensing', ...) on a scratch copy of the template
% (mm_template), on the shipped hard occulter (r_fpm 2.70 lambda/D, Lyot 0.50), N = 512, 3-15 lambda/D,
% the committed decks (ctb_dm.in from ctb_dcr.in, PLAN 1.2 (zzz)).
%   1. the Jacobian, measured live (ctb_dm_jacobian, central differences at 2 nm, not saved);
%   2. one noiseless estimate at the flat DMs, scored against the engine's true field;
%   3. a probe-amplitude sweep and a photon-count sweep (estimation floor against each);
%   4. the probe DM's columns re-measured the way ctb_dm_jacobian did before 2026-10-28 (one-sided
%      2 nm pokes), and the amplitude sweep again on them -- what the poke method costs;
%   5. EFC, 10 iterations: perfect sensing, then closed on the noiseless estimate, then on the
%      estimate with photon noise.
% The figures are this driver's own.  The true field is read only for the scores.

[here, check] = mm_template('30_instruments', 'bench_ctb');
work = getenv('MMMAN_WORK');
t0 = tic;
J  = ctb_dm_jacobian('save', false, 'verbose', false);
ch = ctb_chain('rx', fullfile(here, 'ctb_dm.in'), 'model_size', J.N);
L = strings(0);
say = @(varargin) sprintf(varargin{:});

% ---- 2. one noiseless estimate -----------------------------------------------------------------
e0 = ctb_probe_estimate(ch, J, {});
L(end+1) = say('[probe] flat DMs: true contrast %.3e, estimated %.3e', e0.c_true, e0.c_est);
L(end+1) = say('[probe] %d pairs, default amplitude %.2f nm peak: estimation floor %.2e, relative field error %.2e', ...
    size(e0.dE, 2), e0.probe_nm, e0.floor, e0.rel_err);
L(end+1) = say('[probe] diagnostic ratios %s; worst pixel conditioning %.1f', ...
    strjoin(compose('%.3f', e0.ratio), ' '), max(e0.cond2));
q = prctile(e0.err_px, [50 90 99 100]);
L(end+1) = say('[probe] per-pixel error |dE|^2/peak: median %.2e, 90%% %.2e, 99%% %.2e, max %.2e', q);

% ---- 3. sweeps ---------------------------------------------------------------------------------
amp = [0.03 0.1 0.3 1 3];
sw1 = struct('nm', {}, 'floor', {}, 'rel', {}, 'ratio', {});
for f = amp
    e = ctb_probe_estimate(ch, J, {}, 'probe_nm', f * e0.probe_nm);
    sw1(end+1) = struct('nm', e.probe_nm, 'floor', e.floor, 'rel', e.rel_err, 'ratio', e.ratio); %#ok<SAGROW>
    L(end+1) = say('[amp] %6.2f nm: floor %.2e, relative error %.2e, ratios %s', e.probe_nm, e.floor, ...
        e.rel_err, strjoin(compose('%.3f', e.ratio), ' '));
end
pho = [0 1e10 1e8 1e6];
sw2 = struct('ph', {}, 'floor', {}, 'rel', {});
for ph = pho
    e = ctb_probe_estimate(ch, J, {}, 'photons', ph);
    sw2(end+1) = struct('ph', ph, 'floor', e.floor, 'rel', e.rel_err); %#ok<SAGROW>
    L(end+1) = say('[photons] %-8s at the bare peak: floor %.2e, relative error %.2e', ...
        ternary_(ph == 0, 'none', sprintf('%.0e', ph)), e.floor, e.rel_err);
end

% ---- 4. the probe DM's columns by one-sided 2 nm pokes (the method before 2026-10-28) --------------
Js = J;  G = double(J.G);
udm = unique(J.col_dm, 'stable');  dmof = arrayfun(@(c) find(udm == c, 1), J.col_dm);
d1 = J.dm(1);
dm1 = ctb_dm('ielt', d1.ielt, 'ng', d1.ng, 'gdx_mm', d1.gdx_mm, 'nact', d1.nact, ...
    'pitch_mm', d1.pitch_mm, 'beam_d_mm', d1.beam_d_mm, 'coupling', d1.coupling);
d2 = J.dm(2);
dm2 = ctb_dm('ielt', d2.ielt, 'ng', d2.ng, 'gdx_mm', d2.gdx_mm, 'nact', d2.nact, ...
    'pitch_mm', d2.pitch_mm, 'beam_d_mm', d2.beam_d_mm, 'coupling', d2.coupling);
dm1.clear();  dm2.clear();
cols = find(dmof == 1);  h = J.h_mm;  bias = zeros(numel(cols), 1);
E0f = double(ch.run());  e0f = E0f(J.dz_idx);
for i = 1:numel(cols)
    a = zeros(d1.nact^2, 1);  a(J.col_act(cols(i))) = h;
    dm1.apply(a);  Ep = double(ch.run());
    g = (Ep(J.dz_idx) - e0f) / h;
    bias(i) = norm(g - G(:, cols(i))) / norm(G(:, cols(i)));
    G(:, cols(i)) = g;
end
dm1.clear();
Js.G = single(G);
L(end+1) = say('[jac] one-sided 2 nm columns vs the central differences, probe DM: median %.2e, max %.2e', ...
    median(bias), max(bias));
sw3 = struct('nm', {}, 'floor', {}, 'rel', {}, 'ratio', {});
for f = amp(1:4)
    e = ctb_probe_estimate(ch, Js, {}, 'probe_nm', f * e0.probe_nm);
    sw3(end+1) = struct('nm', e.probe_nm, 'floor', e.floor, 'rel', e.rel_err, 'ratio', e.ratio); %#ok<SAGROW>
    L(end+1) = say('[amp, one-sided] %6.2f nm: floor %.2e, relative error %.2e, ratios %s', e.probe_nm, ...
        e.floor, e.rel_err, strjoin(compose('%.3f', e.ratio), ' '));
end

% ---- 5. EFC --------------------------------------------------------------------------------------
PH = 1e8;
oT = ctb_efc('jac', J, 'niter', 10, 'save', false);
oE = ctb_efc('jac', J, 'niter', 10, 'save', false, 'sensing', struct());
oP = ctb_efc('jac', J, 'niter', 10, 'save', false, 'sensing', struct('photons', PH));
L(end+1) = say('[efc] perfect sensing:            %.3e -> %.3e in %d iterations', oT.c_before, oT.c_after, numel(oT.contrast)-1);
L(end+1) = say('[efc] on the noiseless estimate:  %.3e -> %.3e in %d iterations', oE.c_before, oE.c_after, numel(oE.contrast)-1);
L(end+1) = say('[efc] on the estimate, %.0e photons: %.3e -> %.3e in %d iterations', PH, oP.c_before, oP.c_after, numel(oP.contrast)-1);
fprintf('ch08d wall time %.0f s\n', toc(t0));
check();
disp(strjoin(L, newline))
mm_transcript('ch08_probe_run', strjoin(L, newline));

% ---- figures ---------------------------------------------------------------------------------------
f = figure('Visible', 'off', 'Position', [100 100 1100 400]);
subplot(1, 2, 1);
loglog([sw1.nm], [sw1.floor], 'o-', 'LineWidth', 1.5);  hold on
loglog([sw3.nm], [sw3.floor], 's-', 'LineWidth', 1.5);
yline(e0.c_true, 'k--', 'dark-zone contrast', 'LabelHorizontalAlignment', 'left');
grid on;  xlabel('probe amplitude, nm peak');  ylabel('estimation floor, mean |\DeltaE|^2 / peak');
legend('Jacobian as shipped (central differences)', 'probe DM by one-sided 2 nm pokes (before 2026-10-28)', 'Location', 'southeast');
title('Estimation floor against probe amplitude');
subplot(1, 2, 2);
semilogx([sw1.nm], max(abs(reshape([sw1.ratio], [], numel(sw1)) - 1), [], 1), 'o-', 'LineWidth', 1.5);  hold on
semilogx([sw3.nm], max(abs(reshape([sw3.ratio], [], numel(sw3)) - 1), [], 1), 's-', 'LineWidth', 1.5);
grid on;  xlabel('probe amplitude, nm peak');  ylabel('worst |ratio - 1|');
title('Probe-intensity diagnostic');
print(f, fullfile(work, 'ctb_probe_amp.png'), '-dpng', '-r150');  close(f);
mm_figure('ch08_probe_amp', fullfile(work, 'ctb_probe_amp.png'));

f = figure('Visible', 'off', 'Position', [100 100 640 420]);
phx = [sw2.ph];  phx(phx == 0) = 1e12;
loglog(phx, [sw2.floor], 'o-', 'LineWidth', 1.5);  hold on
yline(e0.c_true, 'k--', 'dark-zone contrast', 'LabelHorizontalAlignment', 'left');
grid on;  set(gca, 'XDir', 'reverse');
xticks(sort(phx));  xticklabels(compose('%.0e', sort(phx)));  xticklabels(regexprep(xticklabels, '1e\+12', 'none'));
xlabel('photons per image at the bare peak');  ylabel('estimation floor');
title(sprintf('Estimation floor against photon count (probe %.2f nm)', e0.probe_nm));
print(f, fullfile(work, 'ctb_probe_photons.png'), '-dpng', '-r150');  close(f);
mm_figure('ch08_probe_photons', fullfile(work, 'ctb_probe_photons.png'));

f = figure('Visible', 'off', 'Position', [100 100 700 440]);
n = @(o) 0:numel(o.contrast)-1;
semilogy(n(oT), oT.contrast, 'k-o', 'LineWidth', 1.5);  hold on
semilogy(n(oE), oE.contrast, 'b-s', 'LineWidth', 1.5);
semilogy(n(oE), oE.contrast_est, 'b:s', 'LineWidth', 1.2);
semilogy(n(oP), oP.contrast, 'r-^', 'LineWidth', 1.5);
semilogy(n(oP), oP.contrast_est, 'r:^', 'LineWidth', 1.2);
grid on;  xlabel('EFC iteration');  ylabel('dark-zone mean contrast, 3-15 \lambda/D');
legend('perfect sensing', 'probe estimate: true contrast', 'probe estimate: estimated', ...
    sprintf('%.0e photons: true contrast', PH), sprintf('%.0e photons: estimated', PH), 'Location', 'northeast');
title('CTB EFC closed on pairwise-probe estimates');
print(f, fullfile(work, 'ctb_probe_efc.png'), '-dpng', '-r150');  close(f);
mm_figure('ch08_probe_efc', fullfile(work, 'ctb_probe_efc.png'));

f = figure('Visible', 'off', 'Position', [100 100 1200 380]);
subplot(1, 3, 1);
histogram(log10(e0.err_px), 60);  grid on
xlabel('log_{10} |\DeltaE|^2 / peak, per pixel');  ylabel('pixels');  title('Per-pixel estimation error (no cut)');
img = @(v) accum_(v, e0.dz_idx, J.N, J.center_px, J.lamD_px, J.outer_lamD);
subplot(1, 3, 2);
[A, ax] = img(e0.amp_err);  imagesc(ax, ax, A.', 'AlphaData', ~isnan(A.'));  axis image xy;  colorbar;  clim([-0.2 0.2]);
title('amplitude error, |E_{est}|/|E_{true}| - 1');  xlabel('\lambda/D');  ylabel('\lambda/D');
subplot(1, 3, 3);
[A, ax] = img(e0.phase_err);  imagesc(ax, ax, A.', 'AlphaData', ~isnan(A.'));  axis image xy;  colorbar;  clim([-0.2 0.2]);
title('phase error, rad');  xlabel('\lambda/D');  ylabel('\lambda/D');
print(f, fullfile(work, 'ctb_probe_err.png'), '-dpng', '-r150');  close(f);
mm_figure('ch08_probe_err', fullfile(work, 'ctb_probe_err.png'));

% ---- pins ---------------------------------------------------------------------------------------------
e2 = @(v) regexprep(sprintf('%.1e', v), 'e-0?', 'e-');
mm_pin('ch08_probe_nm',     e0.probe_nm, 'sig5', sprintf('%.1f nm', e0.probe_nm));
mm_pin('ch08_probe_floor',  e0.floor,    'sig5', e2(e0.floor));
mm_pin('ch08_probe_rel',    e0.rel_err,  'sig5', sprintf('%.1f%%', 100*e0.rel_err));
mm_pin('ch08_probe_ratio',  max(abs(e0.ratio - 1)), 'sig5', sprintf('%.1f%%', 100*max(abs(e0.ratio - 1))));
mm_pin('ch08_probe_lin_rel', sw1(1).rel, 'sig5', e2(sw1(1).rel));
mm_pin('ch08_probe_os_rel',  sw3(1).rel, 'sig5', sprintf('%.1f%%', 100*sw3(1).rel));
mm_pin('ch08_probe_bias',    max(bias),  'sig5', sprintf('%.1f%%', 100*max(bias)));
mm_pin('ch08_probe_ph8',     sw2(3).floor, 'sig5', e2(sw2(3).floor));
mm_pin('ch08_probe_ph6',     sw2(4).floor, 'sig5', e2(sw2(4).floor));
mm_pin('ch08_probe_efcT',    oT.c_after, 'sig5', e2(oT.c_after));
mm_pin('ch08_probe_efcE',    oE.c_after, 'sig5', e2(oE.c_after));
mm_pin('ch08_probe_efcP',    oP.c_after, 'sig5', e2(oP.c_after));

function out = ternary_(c, a, b)
    if c, out = a; else, out = b; end
end

function [A, ax] = accum_(v, idx, N, c, lamD, outer)
    F = nan(N);  F(idx) = v;
    w = ceil((outer + 1) * lamD);  ix = c-w : c+w;
    A = F(ix, ix);  ax = (ix - c) / lamD;
end
