% dwdsurf_jwst_figs.m -- deck figures on Luis's deck: jwst_ote_designc.in through
% run_sensitivities (dwdsurf channel), ONE real segment (elt 5 = Seg2) poked in
% Kr / Kc, centre field, under every option.  Writes demo_session/figs/dwdsurf_jwst_ref.png
% and dwdsurf_jwst_opts.png + a stats text file.
run('/home/dcr/dev/MACOS_resources/mmacos/mmacos_setup.m');
od = '/home/dcr/dev/macos/demo_session/figs';
sc = '/tmp/claude-1000/-home-dcr-dev-macos-macos-f90/8399cf58-01e2-4984-a298-34d9f7d62691/scratchpad/s7/jwst_runs';
if ~exist(sc, 'dir'), mkdir(sc); end
RX = '/home/dcr/dev/MACOS_resources/mmacos/templates/50_sensitivities/zoom_5x5/jwst_ote_designc.in';
MODEL = 512;  NG = 63;  STOP = 25;  FOV = 2.90888e-4;  PK = 5;
base = {'fov_rad', FOV, 'channels', "dwdsurf", 'stop_elt', STOP, 'ngridpts', NG, ...
        'model_size', MODEL, 'surf_params', {'Kr','Kc'}, 'elts', PK, 'out_dir', sc, 'verbose', false};
runs = {'mean_xy',   {'orient','xy',  'opd_ref','mean',  'surf_remove_ptt',false}; ...
        'chief_xy',  {'orient','xy',  'opd_ref','chief', 'surf_remove_ptt',false}; ...
        'mean_ptt',  {'orient','xy',  'opd_ref','mean',  'surf_remove_ptt',true}; ...
        'chief_ptt', {'orient','xy',  'opd_ref','chief', 'surf_remove_ptt',true}; ...
        'mean_raw',  {'orient','raw', 'opd_ref','mean',  'surf_remove_ptt',false}; ...
        'chief_wf',  {'orient','xy',  'opd_ref','chief', 'surf_remove_ptt',false, 'sign','wavefront'}};
if exist(fullfile(sc, 'jwst_option_runs.mat'), 'file')
    load(fullfile(sc, 'jwst_option_runs.mat'), 'R');  fprintf('loaded saved runs\n');
else
    R = struct();
    for r = 1:size(runs,1)
        t0 = tic;
        art = run_sensitivities(RX, base{:}, runs{r,2}{:}, 'name', ['jw_' runs{r,1}]);
        R.(runs{r,1}) = art.os;
        fprintf('run %s done in %.1f min\n', runs{r,1}, toc(t0)/60);
    end
    save(fullfile(sc, 'jwst_option_runs.mat'), 'R', '-v7.3');
end
KR = find(strcmp(R.chief_xy.channel_names, sprintf('Elt %d Kr', PK)), 1);
KC = find(strcmp(R.chief_xy.channel_names, sprintf('Elt %d Kc', PK)), 1);
assert(~isempty(KR) && ~isempty(KC), 'Seg2 columns not found');
fprintf('Seg2 columns: Kr %d, Kc %d of %d channels\n', KR, KC, numel(R.chief_xy.channel_names));
os = R.chief_xy;
fprintf('indxall fields: %s\n', strjoin(fieldnames(os.indxall), ','));
ft = os.field_table;  kC = find(ft(:,1) == 0 & ft(:,2) == 0, 1);
fprintf('fields %d, centre field index %d; channels %s\n', size(ft,1), kC, strjoin(os.channel_names, ','));
% rows of the centre field in the stacked Jacobian
[rowsC, ix, ixC] = centre_rows(os, kC);  ix = ixC;
fprintf('centre-field rows %d of %d; map %s\n', nnz(rowsC), size(os.dwdxall,1), mat2str(ix.size));
% source frame + segment position (orientation anchor)
macos.init(MODEL);  macos.load_rx(RX);
c = macos.get_src_csys();  rp = macos.get_elt_rpt(PK);
tk = regexp(fileread(RX), 'ChfRayDir=\s*(\S+)\s+(\S+)\s+(\S+)', 'tokens', 'once');
cz = cellfun(@(t) str2double(strrep(t, 'D', 'E')), tk).';
fprintf('xGrid %s yGrid %s ChfRayDir %s -> triple product %+.3f; Seg2 (elt %d) RptElt %s mm\n', ...
    mat2str(c.xDir.',4), mat2str(c.yDir.',4), mat2str(cz(:).',4), dot(c.xDir(:), cross(c.yDir(:), cz(:))), PK, mat2str(rp.',5));
% support of the poke = nonzero px of the chief-ref un-PTT'd Kr column (centre field)
vC = os.dwdxall(rowsC, KR);  okm = isfinite(vC);
SUPP = okm & abs(vC) > 1e-9*max(abs(vC(okm)));
[ii, jj] = deal(ixC.i(:), ixC.j(:));
fprintf('poke support: %d px of %d; centroid (i=%.1f, j=%.1f) of a %dx%d map (centre %.1f,%.1f) -> in the xy map the response sits at row offset %+.1f, col offset %+.1f\n', ...
    nnz(SUPP), nnz(okm), mean(ii(SUPP)), mean(jj(SUPP)), ix.size(1), ix.size(2), (ix.size(1)+1)/2, (ix.size(2)+1)/2, ...
    mean(ii(SUPP)) - (ix.size(1)+1)/2, mean(jj(SUPP)) - (ix.size(2)+1)/2);
% raw-orientation support centroid for the caption
osr = R.mean_raw;  [rowsR, ~, ixr] = centre_rows(osr, kC);
% raw column: identify the support by the same pixel SET as xy (transpose)
fprintf('raw map %s; (raw stores index 1 along global X)\n', mat2str(ixr.size));
ink = [11 11 11]/255;  ink2 = [82 81 78]/255;  surf_c = [1 1 1];
cmap = jet(256);          % MATLAB 'jet', to match Luis's plots (Dave 2026-09-09)
stat = @(v) local_stat(v, SUPP, okm);
LITM = unf_field(R.mean_xy, kC, KR, true(size(R.mean_xy.per_field_w_nom_2d{1})));  LIT = isfinite(LITM) & LITM ~= 0;
fprintf('lit pupil mask: %d px (centre-field rows %d)\n', nnz(LIT), nnz(rowsC));
dmc = R.mean_xy.dwdxall(rowsC, KR) - R.chief_xy.dwdxall(rowsC, KR);   % should be one constant
okd = isfinite(dmc);
fprintf('mean - chief (Kr column, all px): mean %.6e std %.3e -> %.2e relative; on the poked support std %.3e, on the others %.3e\n', ...
    mean(dmc(okd)), std(dmc(okd)), std(dmc(okd))/abs(mean(dmc(okd))), std(dmc(okd & SUPP)), std(dmc(okd & ~SUPP)));
fid = fopen(fullfile(od, 'dwdsurf_jwst_stats.txt'), 'w');
fprintf(fid, 'jwst_ote_designc, run_sensitivities dwdsurf, elt %d (Seg2) Kr/Kc, centre field, model %d ng %d stop %d fov %.5g\n', PK, MODEL, NG, STOP, FOV);
fprintf(fid, 'mean - chief (Kr, all px): mean %.6e std %.3e (%.2e relative)\n', mean(dmc(okd)), std(dmc(okd)), std(dmc(okd))/abs(mean(dmc(okd))));
fprintf(fid, 'poke support %d of %d px (%.4f); Seg2 RptElt %s mm; frame triple product %+.3f; support centroid row %+.1f col %+.1f from centre\n', nnz(SUPP), nnz(okm), nnz(SUPP)/nnz(okm), mat2str(rp.',5), dot(c.xDir(:), cross(c.yDir(:), cz(:))), mean(ii(SUPP)) - (ix.size(1)+1)/2, mean(jj(SUPP)) - (ix.size(2)+1)/2);
% ---- figure 1: the reference ---------------------------------------------
f = figure('Color', surf_c, 'Position', [50 50 1400 720], 'Visible', 'off');
tl = tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
tl.Title.String = 'jwst\_ote\_designc, Seg2 (element 5) poked alone through run\_sensitivities: the OPD reference decides what the other 17 segments read';
tl.Title.FontSize = 13;  tl.Title.Color = ink;
P = {R.mean_xy, KR, 'Kr column, opd\_ref = mean (engine default)'; R.chief_xy, KR, 'Kr column, opd\_ref = chief'; ...
     R.mean_xy, KC, 'Kc column, opd\_ref = mean'; R.chief_xy, KC, 'Kc column, opd\_ref = chief'};
for k = 1:4
    o = P{k,1};  cc = P{k,2};  [rr, ~] = centre_rows(o, kC);
    v = o.dwdxall(rr, cc);  w = unf_field(o, kC, cc, LIT);
    s = stat(v);  cst = abs(stat(R.mean_xy.dwdxall(rr, cc)).oth_mean);
    ax = nexttile;  imagesc(ax, w, 'AlphaData', ~isnan(w));  axis(ax, 'image');  axis(ax, 'xy');
    colormap(ax, cmap);  cbh = colorbar(ax);  cbh.Label.String = 'dW/dp';
    title(ax, sprintf('%s\nunpoked segments: %.3g (spread %.1e); poked segment rms %.3g', P{k,3}, s.oth_mean, s.oth_std, s.pk_rms), ...
        'Color', ink, 'FontWeight', 'normal', 'FontSize', 10);
    set(ax, 'XColor', ink2, 'YColor', ink2, 'XTick', [], 'YTick', []);
    fprintf(fid, 'REF %-40s unpoked mean %.6e std %.3e rms %.6e | poked rms %.6e | leak/poked %.4f\n', strrep(P{k,3},'\\',''), s.oth_mean, s.oth_std, s.oth_rms, s.pk_rms, abs(s.oth_mean)/s.pk_rms);
end
exportgraphics(f, fullfile(od, 'dwdsurf_jwst_ref.png'), 'Resolution', 110, 'BackgroundColor', surf_c);
% ---- figure 2: the other options (Kr) --------------------------------------
f = figure('Color', surf_c, 'Position', [50 50 1400 720], 'Visible', 'off');
tl = tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
tl.Title.String = 'The same Kr column under the remaining options: PTT removal, orientation, sign';
tl.Title.FontSize = 13;  tl.Title.Color = ink;
P = {R.mean_ptt, 'surf\_remove\_ptt = true, opd\_ref = mean', R.mean_ptt; R.chief_ptt, 'surf\_remove\_ptt = true, opd\_ref = chief', R.chief_ptt; ...
     R.mean_raw, 'orient = raw (engine array: index 1 along global X), opd\_ref = mean', R.mean_xy; R.chief_wf, 'sign = wavefront (negated), opd\_ref = chief', R.chief_wf};
cst = abs(stat(R.mean_xy.dwdxall(rowsC, KR)).oth_mean);
for k = 1:4
    o = P{k,1};  [rr, ~] = centre_rows(o, kC);
    v = o.dwdxall(rr, KR);  w = unf_field(o, kC, KR, LIT);
    [rs, ~] = centre_rows(P{k,3}, kC);  s = stat(P{k,3}.dwdxall(rs, KR));   % raw panel: stats from the xy twin
    ax = nexttile;  imagesc(ax, w, 'AlphaData', ~isnan(w));  axis(ax, 'image');  axis(ax, 'xy');
    colormap(ax, cmap);  cbh = colorbar(ax);  cbh.Label.String = 'dW/dp';
    title(ax, sprintf('%s\nunpoked segments: rms %.3g (mean %.3g); poked segment rms %.3g', P{k,2}, s.oth_rms, s.oth_mean, s.pk_rms), ...
        'Color', ink, 'FontWeight', 'normal', 'FontSize', 10);
    set(ax, 'XColor', ink2, 'YColor', ink2, 'XTick', [], 'YTick', []);
    fprintf(fid, 'OPT %-60s unpoked mean %.6e rms %.6e | poked rms %.6e\n', strrep(P{k,2},'\\',''), s.oth_mean, s.oth_rms, s.pk_rms);
end
exportgraphics(f, fullfile(od, 'dwdsurf_jwst_opts.png'), 'Resolution', 110, 'BackgroundColor', surf_c);
fclose(fid);
fprintf('JWST FIGS DONE\n');

function w = unf_field(os, kC, col, lit)
% the stacked Jacobian is scattered onto the TILED FIELD CANVAS (dw_multi_core
% [stack]); macos.v2m on indxall rebuilds the canvas, the centre field's tile
% sits at field_table(kC, 3:4) (tile_row, tile_col), one nominal-map size each.
% LIT = the pupil mask (rays present); exact zeros INSIDE it are data (the
% chief-referenced unpoked segments), so they are kept, not blanked.
canv = macos.v2m(os.dwdxall(:, col), os.indxall);
n = size(os.per_field_w_nom_2d{1}, 1);
tr = os.field_table(kC, 3);  tc = os.field_table(kC, 4);
w = canv((tr-1)*n + (1:n), (tc-1)*n + (1:n));
w(~lit) = NaN;
end
function [rows, ix, ixk] = centre_rows(os, kC)
% rows of the centre field in the stacked Jacobian = the rows whose canvas
% pixel falls inside the centre tile
ix = os.indxall;
n = size(os.per_field_w_nom_2d{1}, 1);
tr = os.field_table(kC, 3);  tc = os.field_table(kC, 4);
rows = ix.i(:) > (tr-1)*n & ix.i(:) <= tr*n & ix.j(:) > (tc-1)*n & ix.j(:) <= tc*n;
ixk = struct('i', ix.i(rows) - (tr-1)*n, 'j', ix.j(rows) - (tc-1)*n, 'size', [n n]);
end
function s = local_stat(v, SUPP, okm)
ok = okm & isfinite(v);  supp = SUPP & ok;  oth = ok & ~supp;
s.oth_mean = mean(v(oth));  s.oth_std = std(v(oth));  s.oth_rms = rms(v(oth));  s.pk_rms = rms(v(supp));
end
