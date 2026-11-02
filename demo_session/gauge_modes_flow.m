function gauge_modes_flow(out)
%GAUGE_MODES_FLOW  The DM surface gauge's operating modes as a flow diagram.
%   gauge_modes_flow('figs/gauge_modes_flow.png') draws the four modes
%   (ground flat, image-based phase retrieval, capture, closed-loop hold),
%   the recalibration events, and each reading's capture limit on the
%   arrow into the capture mode (BRIEF_gauge_deck.md sections 7.1 and
%   11.1).  1800 px wide, 17 pt type (the deck's figure rule).  Drawn in
%   MATLAB because graphviz is not on the box the deck is built on.
if nargin < 1, out = fullfile(fileparts(mfilename('fullpath')), 'figs', 'gauge_modes_flow.png'); end
W = 1800;  H = 1000;
f = figure('Visible', 'off', 'Units', 'pixels', 'Position', [0 0 W H], 'Color', 'w');
ax = axes(f, 'Position', [0 0 1 1], 'XLim', [0 W], 'YLim', [0 H], 'YDir', 'normal', 'Visible', 'off');
hold(ax, 'on');
fs = 16;  fsh = 19;  fss = 15;
blue = [0.20 0.35 0.65];  gray = [0.40 0.40 0.40];  green = [0.15 0.50 0.25];  red = [0.70 0.15 0.15];
arrow_ = @(x1, y1, x2, y2, c) quiver(ax, x1, y1, x2 - x1, y2 - y1, 0, 'Color', c, 'LineWidth', 2.5, 'MaxHeadSize', 0.5);
% ---- the four mode boxes, left to right -------------------------------
bw = 400;  bh = 270;  y0 = 430;  xs = [30 470 910 1350];
titles = {'Mode 0   Ground flat', 'Mode 1   Phase retrieval', 'Mode 2   Capture', 'Mode 3   Closed-loop hold'};
froms  = {'on the ground: ~0 WFE', 'from 50-100 nm of surface', 'from where the retrieval stops', 'from ~30 nm of surface'};
bodies = {sprintf('the ground-calibrated voltage map\napplied; on orbit the residual is\nlaunch, gravity release and thermal:\n100-200 nm WFE (50-100 nm of surface)'), ...
          sprintf('the WFS&C loop''s own focal-plane\nphase retrieval, no gauge; wraps at\nhalf a wave of high-spatial-frequency\nWFE (one wave before the gauges do)'), ...
          sprintf('an externally referenced reading,\nunwrapped, the matrix re-measured on\nthe surface as it moves; or the sensor\nwith a second color'), ...
          sprintf('the sensor at picometers: stepped or\nvector Zernike, or the pinhole;\ngain 0.5; the matrix measured on the\nheld surface')};
for k = 1:4
    rectangle(ax, 'Position', [xs(k) y0 bw bh], 'Curvature', 0.08, 'EdgeColor', blue, 'LineWidth', 2.5, 'FaceColor', [0.95 0.97 1.0]);
    text(ax, xs(k) + bw/2, y0 + bh - 28, titles{k}, 'FontSize', fsh, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Color', blue);
    text(ax, xs(k) + bw/2, y0 + bh - 62, froms{k}, 'FontSize', fss, 'FontAngle', 'italic', 'HorizontalAlignment', 'center', 'Color', gray);
    text(ax, xs(k) + 12, y0 + bh/2 - 40, bodies{k}, 'FontSize', fs - 0.5, 'VerticalAlignment', 'middle');
end
for k = 1:3
    arrow_(xs(k) + bw + 3, y0 + bh/2, xs(k+1) - 4, y0 + bh/2, gray);
end
% ---- what captures into mode 2: the readings' limits (section 7.1) ------
cx = xs(3) + bw/2;
arrow_(cx, y0 - 66, cx, y0 - 5, gray);
text(ax, cx, y0 - 92, 'Which reading captures, from how far: the largest starting surface brought to 3 pm', ...
    'FontSize', fs, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
caps = {'Interferometer, lens rig:  150-300 nm of surface, unwrapping alone', ...
        'P/SRI:  100 nm, unwrapping + recalibration on the way down', ...
        'Pinhole with a shutter frame:  100 nm, unwrapping + recalibration', ...
        'Stepped / vector Zernike:  30 / 60 nm  (self-referenced: no reach beyond)', ...
        'Interferometer on the OAP rig:  does not capture'};
cols = {green, green, green, red, red};
for k = 1:numel(caps)
    yk = y0 - 140 - 46*(k-1);
    text(ax, cx - 470, yk, ['\bullet   ' caps{k}], 'FontSize', fs, 'Color', cols{k});
end
% ---- recalibration events (above the modes) ------------------------------
yr = y0 + bh + 80;  rh = 90;
rectangle(ax, 'Position', [xs(1) yr bw rh], 'Curvature', 0.5, 'EdgeColor', green, 'LineWidth', 2, 'FaceColor', [0.93 0.98 0.94]);
text(ax, xs(1) + bw/2, yr + rh/2, sprintf('The flat re-taken\n(the ground reference)'), 'FontSize', fss, 'HorizontalAlignment', 'center', 'Color', green);
arrow_(xs(1) + bw/2, yr - 2, xs(1) + bw/2, y0 + bh + 5, green);
rx = xs(3) + 60;  rw = xs(4) + bw - 60 - rx;
rectangle(ax, 'Position', [rx yr rw rh], 'Curvature', 0.5, 'EdgeColor', green, 'LineWidth', 2, 'FaceColor', [0.93 0.98 0.94]);
text(ax, rx + rw/2, yr + rh/2, sprintf('Recalibration: the response matrix re-measured through the sensor\non the current surface (photon cost 5-40x of a measurement)'), ...
    'FontSize', fss, 'HorizontalAlignment', 'center', 'Color', green);
arrow_(xs(3) + bw/2, yr - 2, xs(3) + bw/2, y0 + bh + 5, green);
text(ax, xs(3) + bw/2 - 10, y0 + bh + 40, 'every few cycles of the descent', 'FontSize', fss - 1, 'Color', green, 'HorizontalAlignment', 'right');
arrow_(xs(4) + bw/2, yr - 2, xs(4) + bw/2, y0 + bh + 5, green);
text(ax, xs(4) + bw/2 - 10, y0 + bh + 40, 'when the calibration ages (gain off by 10%)', 'FontSize', fss - 1, 'Color', green, 'HorizontalAlignment', 'right');
% ---- title --------------------------------------------------------------
text(ax, W/2, H - 45, 'DM surface gauge: the operating modes, and what carries the surface between them', ...
    'FontSize', 22, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
set(f, 'PaperPositionMode', 'auto');
print(f, out, '-dpng', '-r96');
close(f);
fprintf('wrote %s\n', out);
end
