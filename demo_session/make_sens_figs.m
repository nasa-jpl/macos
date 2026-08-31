% MAKE_SENS_FIGS  Readable sensitivity-channel figures for deck slides
% (Dave 2026-08-31: the all-channel stack is unreadable; show 1-2 segments,
% centre-field + 5-fov, for each of dw/dx, dw/dz, dw/dgrid).
%
% Data: the committed r3 harvest (e2e6m_r2/r3_sens.mat) from stage
% r3_sensitivities -> run_sensitivities -> macos.dw_dx_multi /
% dw_dz_zernike_multi / dw_dgrid_multi.  Canvas via macos.v2m is a 3x3
% quincunx mosaic (centre + 4 corner fields), 255 px tiles; centre-only =
% the middle tile.  Display conventions from plot_dw_channels: NaN outside
% the pupil, piston removed, per-panel autoscale.
%
% Run: MACOS_HOME=~/dev/macos/macos_f90 matlab -batch "run('make_sens_figs.m')"
% Writes figs/fig_dwdx_read.png, fig_dwdz_read.png, fig_dwdgrid_read.png.

run('/home/dcr/dev/MACOS_res_dev/mmacos/mmacos_setup.m');
S = load(['/home/dcr/dev/MACOS_res_dev/mmacos/templates/80_end_to_end/' ...
          'e2e6m_r2/r3_sens.mat']);
here = fileparts(mfilename('fullpath'));

fam(1) = struct('key',"ox", 'chan',{["Elt 1 Rx","Elt 8 Rx"]}, ...
                'png',"fig_dwdx_read.png",  'ttl',"dw/dx — rigid body");
fam(2) = struct('key',"oz", 'chan',{["Elt 1 MonZern4","Elt 8 MonZern4"]}, ...
                'png',"fig_dwdz_read.png",  'ttl',"dw/dz — segment figure (MonZernike)");
fam(3) = struct('key',"og", 'chan',{["Elt 1 Grid4","Elt 8 Grid4"]}, ...
                'png',"fig_dwdgrid_read.png", 'ttl',"dw/dgrid — segment influence grid");

tl0 = 256; tl1 = 510;   % centre tile of the 765x765 quincunx
for f = 1:3
    o  = S.(fam(f).key);
    cn = string(o.channel_names);
    fig = figure('Position',[40 40 1560 440], 'Color','w', 'Visible','off');
    t = tiledlayout(fig, 1, 4, 'Padding','compact', 'TileSpacing','compact');
    for j = 1:2
        k = find(cn == fam(f).chan(j), 1);
        assert(~isempty(k), 'channel %s not found', fam(f).chan(j));
        M = macos.v2m(o.dwdxall(:,k), o.indxall);
        M(M == 0) = NaN;
        M = M - mean(M(:), 'omitnan');
        % quincunx sanity: edge-middle tiles empty, centre tile occupied
        assert(all(isnan(M(tl0:tl1, 1:255)), 'all') && ...
               any(~isnan(M(tl0:tl1, tl0:tl1)), 'all'), 'not a quincunx');
        C = M(tl0:tl1, tl0:tl1);
        cl = prctile(M(~isnan(M)), [1 99]);   % edge spikes must not own the scale
        if cl(1) >= cl(2), cl = [min(M(:),[],'omitnan') max(M(:),[],'omitnan')]; end
        ax = nexttile(t);
        h = imagesc(ax, C);  set(h, 'AlphaData', ~isnan(C));
        axis(ax, 'image', 'off');  clim(ax, cl);
        title(ax, sprintf('%s — centre field', fam(f).chan(j)), 'FontSize', 12);
        ax = nexttile(t);
        h = imagesc(ax, M);  set(h, 'AlphaData', ~isnan(M));
        axis(ax, 'image', 'off');  clim(ax, cl);  hold(ax, 'on');
        for g = [255.5 510.5]
            plot(ax, [g g], [0.5 765.5], '-', 'Color', [0.85 0.85 0.85]);
            plot(ax, [0.5 765.5], [g g], '-', 'Color', [0.85 0.85 0.85]);
        end
        title(ax, 'all 5 fields (centre + 4 corners)', 'FontSize', 12);
    end
    title(t, sprintf('%s — two of the %d channels, OPD at the coronagraph exit pupil', ...
          fam(f).ttl, numel(cn)), 'FontSize', 13, 'FontWeight', 'bold');
    exportgraphics(fig, fullfile(here, 'figs', char(fam(f).png)), 'Resolution', 150);
    fprintf('wrote %s\n', fam(f).png);
end
