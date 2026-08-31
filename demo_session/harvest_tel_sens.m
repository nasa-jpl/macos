% HARVEST_TEL_SENS  Small telescope-only sensitivity harvest for the deck's
% control-basis slides (Dave 2026-08-31: the 5-fov panels on the coronagraph
% train show its field stop vignetting the corner fields -- expected, the
% coronagraph is narrow-field; show the telescope's own full-field channels
% instead, and label the coronagraph channel at its centre field).
%
% Harvests dw/dx + dw/dz + dw/dgrid for segments E1 and E8 on the IMAGING
% LEG deck s3_imager_full.in -- the full-field instrument (same telescope
% trunk as the r2 train; s2_segmented's hand-built EP-return sphere kills
% corner fields even under per-field FEX, so the imaging leg, which ends
% at a plain focal plane, carries the 5-fov panels) -- via the same run_sensitivities
% supervisor and P.sn settings as stage r3.  Saves sens_tel.mat here
% (NOT committed -- ~50 MB class; the committed figures are the artifact).
% Also attributes the coronagraph-train corner-field vignetting to its
% element, by per-plane ray survival at a corner field.

run('/home/dcr/dev/MACOS_res_dev/mmacos/mmacos_setup.m');
MM = '/home/dcr/dev/MACOS_res_dev/mmacos';
R2 = fullfile(MM, 'templates', '80_end_to_end', 'e2e6m_r2');
E1 = fullfile(MM, 'templates', '80_end_to_end', 'e2e6m');
addpath(R2);
addpath(fullfile(MM, 'design', 'runners'));
addpath(fullfile(MM, 'design', 'src'));
addpath(fullfile(MM, 'sensitivities'));
here = fileparts(mfilename('fullpath'));

P = e2e6m_r2_params(struct());
macos.init(P.sn.model);

art = run_sensitivities(string(fullfile(E1, 's3_imager_full.in')), ...
        'fov_rad', deg2rad(P.tel.fov_arcmin/60), ...
        'elts', [1; 8], ...
        'zmodes_fig',  P.sn.zmodes_fig, ...
        'zmodes_grid', P.sn.zmodes_grid, ...
        'ng', P.sn.ng, ...
        'model_size', P.sn.model, ...
        'out_dir', string(tempdir), 'name', "tel", ...
        'verbose', false);
ox = art.ox;  oz = art.oz;  og = art.og;
save(fullfile(here, 'sens_tel.mat'), 'ox', 'oz', 'og');
fprintf('saved sens_tel.mat: ox %d, oz %d, og %d channels\n', ...
        size(ox.dwdxall,2), size(oz.dwdxall,2), size(og.dwdxall,2));

% ---- attribute the coronagraph-train corner-field vignetting ------------
deck = fullfile(R2, 'r1_seg_prop.in');
txt = fileread(deck);
g3 = @(key) subsref(sscanf(char(subsref(regexp(txt, [key '=\s*([^\n]*)'], ...
     'tokens','once'), struct('type','{}','subs',{{1}}))), '%f'), ...
     struct('type','()','subs',{{1:3}}));
cd0 = g3('ChfRayDir');  cp0 = g3('ChfRayPos');  ap = g3('ApStop');
stand = dot(ap - cp0, cd0);
f = deg2rad(P.tel.fov_arcmin/60) * [1 1];        % a corner of the box
bx = asin(cd0(1)) + f(1);  by = asin(cd0(2)) + f(2);
cdir = [sin(bx); sin(by); sqrt(max(0, 1 - sin(bx)^2 - sin(by)^2))];
cpos = ap - stand*cdir;
v3 = @(x) sprintf('%+.15e %+.15e %+.15e', x(1), x(2), x(3));
s = regexprep(txt, '(ChfRayDir=\s*)[^\n]*', ['$1' v3(cdir)]);
s = regexprep(s,   '(ChfRayPos=\s*)[^\n]*', ['$1' v3(cpos)]);
s = regexprep(s,   '(yGrid=\s*)[^\n]*',     ['$1' v3([0; cos(by); -sin(by)])]);
tmp = [tempname '.in'];
fid = fopen(tmp, 'w');  fprintf(fid, '%s', s);  fclose(fid);
macos.load_rx(tmp);
macos.ray_hist('on');  t = macos.trace();  hi = macos.ray_hist(t.nRays);
macos.ray_hist('off');  delete(tmp);
nm = regexp(txt, 'EltName=\s*(\S+)', 'tokens');
nm = cellfun(@(c) c{1}, nm, 'UniformOutput', false);
nE = numel(nm);  off = size(hi.P, 3) - nE;
fr = squeeze(sum(hi.ok(:, off+1:end), 1)) / size(hi.ok, 1);
fprintf('corner-field ray survival along the coronagraph train:\n');
prev = fr(1);
for k = 1:nE
    if fr(k) < prev - 0.02
        fprintf('  %-14s %5.1f%%  (drop %.1f%%)\n', nm{k}, 100*fr(k), ...
                100*(prev - fr(k)));
    end
    prev = fr(k);
end
fprintf('  final survival %.1f%%\n', 100*fr(end));
