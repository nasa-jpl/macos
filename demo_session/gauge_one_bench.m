function gauge_one_bench(out)
%GAUGE_ONE_BENCH  One bench, all modes: the switching schematic.
%   Draws the shared front end, the interferometer's reference arm (shutter),
%   the mask seat (translating substrate), the vector tail (plate in/out,
%   cube, two cameras) and the flip-in P/SRI arm, with each switch marked
%   and the modes it selects (BRIEF_gauge_deck.md section 11.3).  A
%   schematic, not a ray trace: the engine-drawn layouts of each rig are on
%   their own slides.  1800 px, 16-19 pt.
if nargin < 1, out = fullfile(fileparts(mfilename('fullpath')), 'figs', 'gauge_one_bench.png'); end
W = 1800;  H = 1000;
f = figure('Visible', 'off', 'Units', 'pixels', 'Position', [0 0 W H], 'Color', 'w');
ax = axes(f, 'Position', [0 0 1 1], 'XLim', [0 W], 'YLim', [0 H], 'YDir', 'normal', 'Visible', 'off');
hold(ax, 'on');
fs = 16;  fss = 14;
blue = [0.20 0.35 0.65];  gray = [0.35 0.35 0.35];  green = [0.15 0.50 0.25];  beam = [0.85 0.30 0.20];
box_ = @(x, y, w, h, c, fc) rectangle(ax, 'Position', [x y w h], 'Curvature', 0.15, 'EdgeColor', c, 'LineWidth', 2.2, 'FaceColor', fc);
lab_ = @(x, y, s, c, varargin) text(ax, x, y, s, 'FontSize', fs, 'HorizontalAlignment', 'center', 'Color', c, varargin{:});
seg_ = @(x1, y1, x2, y2, c, w) plot(ax, [x1 x2], [y1 y2], '-', 'Color', c, 'LineWidth', w);
arr_ = @(x1, y1, x2, y2, c) quiver(ax, x1, y1, x2-x1, y2-y1, 0, 'Color', c, 'LineWidth', 2.2, 'MaxHeadSize', 0.6);
% ---- the main line: source -> L1 -> BS -> DM (y = 600) ----------------------
y = 600;
box_(40, y-40, 130, 80, blue, [0.95 0.97 1]);   lab_(105, y, sprintf('laser\n632.8 nm'), blue);
seg_(170, y, 300, y, beam, 3);
box_(300, y-40, 90, 80, blue, [0.95 0.97 1]);   lab_(345, y, 'L1', blue);
seg_(390, y, 560, y, beam, 3);
box_(560, y-45, 60, 90, blue, [0.95 0.97 1]);   text(ax, 545, y-70, sprintf('plate splitter\n7 deg'), 'FontSize', fs, 'Color', blue, 'HorizontalAlignment', 'right');
seg_(620, y, 900, y, beam, 3);
lab_(760, y+24, '700 mm leg', gray, 'FontSize', fss);
box_(900, y-70, 70, 140, blue, [0.95 0.97 1]);  lab_(935, y, sprintf('DM\n96x96'), blue);
% ---- the reference arm, up from the splitter -------------------------------
seg_(590, y+45, 590, 760, beam, 3);
box_(560, 760, 60, 60, green, [0.93 0.98 0.94]);  lab_(590, 790, 'S', green, 'FontWeight', 'bold');
lab_(640, 790, sprintf('shutter: open = interferometer\nclosed = every sensor mode'), green, 'FontSize', fss, 'HorizontalAlignment', 'left');
seg_(590, 820, 590, 860, beam, 3);
box_(500, 860, 180, 44, blue, [0.95 0.97 1]);  lab_(590, 882, 'reference flat on PZT', blue, 'FontSize', fss);
% ---- the detector tail, down from the splitter ------------------------------
seg_(590, y-45, 590, 470, beam, 3);
box_(560, 400, 60, 70, blue, [0.95 0.97 1]);  lab_(590, 435, 'L2', blue);
seg_(590, 400, 590, 330, beam, 3);
% the mask seat (translating)
box_(440, 270, 300, 60, green, [0.93 0.98 0.94]);
lab_(590, 300, sprintf('mask seat: one substrate,\ntranslated'), green, 'FontWeight', 'bold', 'FontSize', fss);
text(ax, 425, 300, sprintf('dimple | metasurface |\npinhole | clear window'), 'FontSize', fss, 'Color', green, 'HorizontalAlignment', 'right');
text(ax, 755, 300, sprintf('= Zernike | vector Zernike |\npinhole | interferometer, retrieval'), 'FontSize', fss, 'Color', green, 'HorizontalAlignment', 'left');
seg_(590, 270, 590, 200, beam, 3);
box_(560, 130, 60, 70, blue, [0.95 0.97 1]);  lab_(590, 165, sprintf('field\nlens'), blue, 'FontSize', fss);
seg_(590, 130, 590, 80, beam, 3);
% the vector tail: plate (in/out), cube, two cameras
box_(560, 40, 60, 40, green, [0.93 0.98 0.94]);  lab_(590, 60, 'QWP', green, 'FontWeight', 'bold', 'FontSize', fss);
lab_(450, 60, sprintf('plate in = vector split\nout = one camera'), green, 'FontSize', fss);
seg_(620, 60, 760, 60, beam, 3);
box_(760, 20, 80, 80, blue, [0.95 0.97 1]);  lab_(800, 60, sprintf('PBS\ncube'), blue, 'FontSize', fss);
seg_(840, 60, 960, 60, beam, 3);  box_(960, 35, 110, 50, blue, [0.95 0.97 1]);  lab_(1015, 60, 'camera A', blue, 'FontSize', fss);
seg_(800, 100, 800, 160, beam, 3);  box_(745, 160, 110, 50, blue, [0.95 0.97 1]);  lab_(800, 185, 'camera B', blue, 'FontSize', fss);
lab_(1015, 100, 'pupil image', gray, 'FontSize', fss);
% ---- the flip-in pickoff and the P/SRI arm (to the right of the tail) -------
box_(1120, 400, 60, 70, green, [0.93 0.98 0.94]);  lab_(1150, 435, 'flip', green, 'FontWeight', 'bold', 'FontSize', fss);
lab_(1150, 490, sprintf('flip-in pickoff plate:\nin = P/SRI, out = all others'), green, 'FontSize', fss);
seg_(620, 435, 1120, 435, beam, 2);   % the pickoff takes light from the tail before the mask seat
plot(ax, 620, 435, 'o', 'MarkerSize', 7, 'MarkerFaceColor', beam, 'MarkerEdgeColor', beam);
xs = [1240 1360 1480 1630];  ws = [90 90 120 110];  names = {'Lr1', sprintf('pinhole\n+ Lr2'), sprintf('waveguide +\nphase shifter'), sprintf('BS3 +\ncamera C')};
for k = 1:4
    box_(xs(k), 400, ws(k), 70, blue, [0.95 0.97 1]);  lab_(xs(k)+ws(k)/2, 435, names{k}, blue, 'FontSize', fss);
    if k < 4, seg_(xs(k)+ws(k), 435, xs(k+1), 435, beam, 3); end
end
seg_(1180, 435, 1240, 435, beam, 3);
lab_(1490, 360, sprintf('the P/SRI reference arm on its own breadboard\n(folds and compensator not drawn)'), gray, 'FontSize', fss);
% ---- the interferometer forms note ------------------------------------------
box_(1120, 620, 640, 130, green, [0.93 0.98 0.94]);
text(ax, 1140, 685, sprintf(['interferometer forms on this bench:\n' ...
    '  PZT four-step: the shutter open, nothing else\n' ...
    '  polarization snapshot: in-arm quarter-wave plates in, analyzer in\n' ...
    '  hybrid: both; the v2 cemented cube does not switch in (backup)']), 'FontSize', fss, 'Color', green, 'VerticalAlignment', 'middle');
% ---- legend + title -----------------------------------------------------------
box_(1120, 800, 30, 30, green, [0.93 0.98 0.94]);  text(ax, 1165, 815, 'a switch: inserted, removed, translated or shuttered -- nothing realigned', 'FontSize', fss, 'Color', green);
box_(1120, 850, 30, 30, blue, [0.95 0.97 1]);  text(ax, 1165, 865, 'a part that stays built', 'FontSize', fss, 'Color', blue);
text(ax, W/2, H - 30, 'One bench, every mode: what switches', 'FontSize', 22, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
set(f, 'PaperPositionMode', 'auto');
print(f, out, '-dpng', '-r96');
close(f);
fprintf('wrote %s\n', out);
end
