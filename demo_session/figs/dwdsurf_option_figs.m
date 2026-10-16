% dwdsurf_option_figs.m -- deck figures: ONE segment's Kr / Kc column on e5hex1
% under every driver option: opd_ref mean|chief, remove_ptt off|on, orient
% raw|xy, sign opl|wavefront.  Writes figs/dwdsurf_ref.png, figs/dwdsurf_opts.png
run('/home/dcr/dev/MACOS_resources/mmacos/mmacos_setup.m');
od = '/home/dcr/dev/macos/demo_session/figs';  if ~exist(od, 'dir'), mkdir(od); end
rx = '/home/dcr/dev/MACOS_resources/pymacos/tests/Rx/e5hex1.in';
m = macos.Session(128);  NG = 63;  PK = 3;
base = {'elts', PK, 'params', {'Kr','Kc'}, 'ngridpts', NG};
R = struct();
R.mean_xy   = macos.dw_dsurf(m, rx, base{:}, 'orient','xy', 'remove_ptt',false, 'opd_ref','mean');
R.chief_xy  = macos.dw_dsurf(m, rx, base{:}, 'orient','xy', 'remove_ptt',false, 'opd_ref','chief');
R.mean_ptt  = macos.dw_dsurf(m, rx, base{:}, 'orient','xy', 'remove_ptt',true,  'opd_ref','mean');
R.chief_ptt = macos.dw_dsurf(m, rx, base{:}, 'orient','xy', 'remove_ptt',true,  'opd_ref','chief');
R.mean_raw  = macos.dw_dsurf(m, rx, base{:}, 'orient','raw','remove_ptt',false, 'opd_ref','mean');
R.chief_wf  = macos.dw_dsurf(m, rx, base{:}, 'orient','xy', 'remove_ptt',false, 'opd_ref','chief', 'sign','wavefront');
save(fullfile(od, 'dwdsurf_option_maps.mat'), 'R');
ink = [11 11 11]/255;  ink2 = [82 81 78]/255;  surf_c = [1 1 1];
% diverging map: blue (#2a78d6) - white - orange (#eb6834), the validated pair
cb = [42 120 214]/255;  co = [235 104 52]/255;  nn = 128;
cmap = [ [linspace(cb(1),1,nn).' linspace(cb(2),1,nn).' linspace(cb(3),1,nn).']; ...
         [linspace(1,co(1),nn).' linspace(1,co(2),nn).' linspace(1,co(3),nn).'] ];
unf = @(o, c) local_unf(o, c);
% the poked segment's rays = the support of the chief-referenced, un-PTT'd column
SUPP = isfinite(R.chief_xy.dwds(:,1)) & abs(R.chief_xy.dwds(:,1)) > 1e-9*max(abs(R.chief_xy.dwds(:,1)));
OKM  = isfinite(R.chief_xy.dwds(:,1));
stat = @(v, o) local_stat(v, o, SUPP, OKM);
% ---- figure 1: the reference (Kr and Kc, mean vs chief) ----------------
f = figure('Color', surf_c, 'Position', [50 50 1400 700], 'Visible', 'off');
tl = tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
tl.Title.String = 'e5hex1, segment 2 (element 3) poked alone: the OPD reference decides what the other six segments read';
tl.Title.FontSize = 13;  tl.Title.Color = ink;
P = {R.mean_xy, 1, 'Kr column, opd\_ref = mean (engine default)'; R.chief_xy, 1, 'Kr column, opd\_ref = chief'; ...
     R.mean_xy, 2, 'Kc column, opd\_ref = mean'; R.chief_xy, 2, 'Kc column, opd\_ref = chief'};
for k = 1:4
    ax = nexttile;  [w, s] = deal(unf(P{k,1}, P{k,2}), stat(P{k,1}.dwds(:,P{k,2}), P{k,1}));
    imagesc(ax, w, 'AlphaData', ~isnan(w));  axis(ax, 'image');  axis(ax, 'xy');
    % colour scale = +-3 x the mean-reference constant of this column, so the
    % flat offset on the six unpoked segments is visible; the poked segment
    % saturates (its rms is 8-9x the constant)
    cst = abs(stat(R.mean_xy.dwds(:,P{k,2}), R.mean_xy).oth_mean);
    colormap(ax, cmap);  clim(ax, [-3*cst 3*cst]);  cbh = colorbar(ax);
    cbh.Label.String = 'dW/dp';
    title(ax, sprintf('%s\nunpoked segments: %.3g (spread %.1e); poked segment rms %.3g', P{k,3}, s.oth_mean, s.oth_std, s.pk_rms), ...
        'Color', ink, 'FontWeight', 'normal', 'FontSize', 10);
    set(ax, 'XColor', ink2, 'YColor', ink2, 'XTick', [], 'YTick', []);
end
exportgraphics(f, fullfile(od, 'dwdsurf_ref.png'), 'Resolution', 110, 'BackgroundColor', surf_c);
% ---- figure 2: the other options on the Kr column ------------------------
f = figure('Color', surf_c, 'Position', [50 50 1400 700], 'Visible', 'off');
tl = tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
tl.Title.String = 'The same Kr column under the remaining options: PTT removal, orientation, sign';
tl.Title.FontSize = 13;  tl.Title.Color = ink;
P = {R.mean_ptt, 1, 'remove\_ptt = true, opd\_ref = mean'; R.chief_ptt, 1, 'remove\_ptt = true, opd\_ref = chief'; ...
     R.mean_raw, 1, 'orient = raw (engine array; index 1 = global X), opd\_ref = mean'; R.chief_wf, 1, 'sign = wavefront (negated), opd\_ref = chief'};
STATSRC = {R.mean_ptt, R.chief_ptt, R.mean_xy, R.chief_wf};   % raw panel: stats from its xy twin
for k = 1:4
    ax = nexttile;  [w, s] = deal(unf(P{k,1}, P{k,2}), stat(STATSRC{k}.dwds(:,P{k,2}), STATSRC{k}));
    imagesc(ax, w, 'AlphaData', ~isnan(w));  axis(ax, 'image');  axis(ax, 'xy');
    cst = abs(stat(R.mean_xy.dwds(:,1), R.mean_xy).oth_mean);
    colormap(ax, cmap);  clim(ax, [-3*cst 3*cst]);  cbh = colorbar(ax);
    cbh.Label.String = 'dW/dp';
    title(ax, sprintf('%s\nunpoked segments: rms %.3g (mean %.3g); poked segment rms %.3g', P{k,3}, s.oth_rms, s.oth_mean, s.pk_rms), ...
        'Color', ink, 'FontWeight', 'normal', 'FontSize', 10);
    set(ax, 'XColor', ink2, 'YColor', ink2, 'XTick', [], 'YTick', []);
end
exportgraphics(f, fullfile(od, 'dwdsurf_opts.png'), 'Resolution', 110, 'BackgroundColor', surf_c);
fprintf('FIGS DONE\n');

function w = local_unf(o, c)
ix = o.indx;  w = nan(ix.size);  w(sub2ind(ix.size, ix.i, ix.j)) = o.dwds(:, c);
end
function s = local_stat(v, o, SUPP, OKM) %#ok<INUSL>
% statistics over the xy-layout column: the poked segment's rays = SUPP (the
% support of the chief-referenced un-PTT'd column), the rest = the six
% unpoked segments.  Callers pass the xy TWIN's column for the raw panel.
ok = OKM & isfinite(v);  supp = SUPP & ok;  oth = ok & ~supp;
s.oth_mean = mean(v(oth));  s.oth_std = std(v(oth));  s.oth_rms = rms(v(oth));  s.pk_rms = rms(v(supp));
end
