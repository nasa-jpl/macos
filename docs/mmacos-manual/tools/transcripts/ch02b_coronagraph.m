% ch02b_coronagraph.m -- Chapter 2, the coronagraph channel: the CTB (templates/30_instruments/
% bench_ctb) built as a Bench behind the thread's stage-2 telescope (MMMAN_THREAD/s2_instrument.in)
% and spliced onto it with append_rx.  PLAN 1.2 (xxx), case C, ruled 2026-10-27: the telescope keeps
% its arcminute field; the coronagraph takes an arcsecond sub-field about the stage-2 field center.
% example_ctb.m's PARAMETER BLOCK and its rot_axis fold helper are TAKEN from the template file at run
% time (the template is unchanged); one number differs from it, OAP1's focal length (f/# x the CTB's
% DM pupil diameter), and one spacing, OAP1 -> DM1 (the telescope stop's image, measured here).
% The figures draw the FOLDED instrument: e2e stage 8's recorded pick (e2e_params P.coro.pick,
% PLAN_CONSOLIDATION 3a, CCL 2026-10-29) replayed on this thread's telescope by the library's
% ctb_fold_build / ctb_fold_measure, and the unfolded case C by the same code for the comparison.
% Products land in the scratch directory (MMMAN_WORK); the thread is read, never written.

root = getenv('MMACOS_ROOT');  work = getenv('MMMAN_WORK');  thread = getenv('MMMAN_THREAD');
tel  = fullfile(thread, 's2_instrument.in');
if ~isfile(tel), error('ch02b:thread', 'no %s -- run ch02a_imager first', tel); end
ctb  = fullfile(root, 'templates', '30_instruments', 'bench_ctb');
addpath(fullfile(root, 'src'));  addpath(fullfile(root, 'design', 'src'));
cd(work);

% ---- example_ctb's parameter block and fold helper, taken from the template --------------------
ex = fullfile(ctb, 'example_ctb.m');
el = splitlines(string(fileread(ex)));
k0 = find(contains(el, "PARAMETER BLOCK"), 1) + 2;            % past the '% ====' rule
k1 = k0 - 1 + find(strtrim(el(k0:end)) == "", 1) - 1;          % through the line before the blank
mm_listing('ch02_ctb_params', ex, [k0 k1]);
eval(char(strjoin(el(k0:k1), newline)));                       % defines P (lengths mm)
f0 = find(startsWith(el, "function d = rot_axis"), 1);
f1 = f0 - 1 + find(el(f0:end) == "end", 1);
writelines(el(f0:f1), fullfile(work, 'rot_axis.m'));           % the template's helper, on the path (cwd)
nd = @(d, k) rot_axis(d, [0;0;1], (1 - 2*mod(k,2))*(180 - 2*P.AOI));   % example_ctb's nextdir
dmt = fileread(fullfile(ctb, 'ctb_dm.m'));                     % the DM lattice defaults
nact  = str2double(regexp(dmt, 'opts\.nact\s+\(1,1\)[^=]*=\s*(\d+)', 'tokens', 'once'));
beamd = str2double(regexp(dmt, 'opts\.beam_d_mm\s+\(1,1\)[^=]*=\s*([\d.]+)', 'tokens', 'once'));
pitch = beamd/nact*1e-3;                                       % m
mm = 1e-3;  Dpup = P.R_DM*P.FILL*mm;  rpup = Dpup/2;           % example_ctb: R_DM*FILL lands as the DIAMETER
EDGE_LD = 8;                                                   % the CTB's dark-zone edge, lambda/D (deck_ctb.md)

% ---- the telescope, as the coronagraph sees it -------------------------------------------------
macos.init(256);
nT = macos.load_rx(tel);
r1 = macos.get_ray_info(macos.trace(1).nRays);
sT = macos.trace(nT);  rT = macos.get_ray_info(sT.nRays);
ok = r1.ok_pass & r1.ok_trace & rT.ok_pass & rT.ok_trace;
pF = rT.pos(:,1);  dF = rT.dir(:,1)/norm(rT.dir(:,1));
th = acos(min(1, dF.' * (rT.dir(:,ok)./vecnorm(rT.dir(:,ok)))));
fno = 1/(2*sin(max(th)));                                      % the marginal ray's f/#
u = cross(dF, [0;0;1]);  u = u/norm(u);  v = cross(dF, u);
a = [u v].' * (rT.dir(:,ok)./vecnorm(rT.dir(:,ok)));
Pp = r1.pos(1:2, ok) - r1.pos(1:2, 1);
A = a / [Pp; ones(1, size(Pp,2))];
efl = sort(1 ./ svd(A(:,1:2)));                                % pupil-side EFL per principal axis
txt = fileread(tel);
Dtel = str2double(regexp(txt, 'Aperture=\s*(\S+)', 'tokens', 'once'));
dB = str2double(regexp(txt, 'ChfRayDir=\s*(\S+)\s+(\S+)\s+(\S+)', 'tokens', 'once')).';
pB = str2double(regexp(txt, 'ChfRayPos=\s*(\S+)\s+(\S+)\s+(\S+)', 'tokens', 'once')).';
dB = dB/norm(dB);  L0 = norm(pB);
assert(norm(pB + L0*dB) < 1e-9*L0, 'ch02b: ChfRayPos is not ApStop(0) - L*ChfRayDir');
bias = atan2(dB(2), dB(3)) * 180/pi*60;                        % arcmin, the stage-2 field center
ld = P.WAVLEN*mm / Dtel * 180/pi*3600;                         % lambda/D, arcsec
rep = fileread(fullfile(thread, 's2_report.txt'));
xpd = str2double(regexp(rep, '(\S+) m from the image', 'tokens', 'once'));

% ---- the CTB as built, then OAP1 matched, then DM1 moved to the stop image ---------------------
fOAP1 = fno*Dpup;
[fa, ia] = bench_(P, P.F_OAP(1)*mm, P.L_DM1*mm, pF, dF, fno, nd, tel, work, 'as_built');
macos.init(256);  macos.load_rx(fa);
nA = macos.trace(ia.DM1).nRays;  rA = macos.get_ray_info(nA);  okA = nnz(rA.ok_pass & rA.ok_trace);
[fb, ib] = bench_(P, fOAP1, P.L_DM1*mm, pF, dF, fno, nd, tel, work, 'matched');
s_img = stop_image_(fileread(fb), dB, L0, ib.OAP1, 3, work);
LDM1 = round(s_img, 4);
[fc, ic] = bench_(P, fOAP1, LDM1, pF, dF, fno, nd, tel, work, 'C');

macos.init(256);  nF = macos.load_rx(fc);
row = strings(0);  names = {'OAP1','DM1','DM2','Lyot','FPA'};  ie = [ic.OAP1 ic.DM1 ic.DM2 ic.Lyot nF];
for j = 1:numel(ie)
    s = macos.trace(ie(j));  ri = macos.get_ray_info(s.nRays);
    okk = ri.ok_pass & ri.ok_trace;  dc = ri.dir(:,1)/norm(ri.dir(:,1));
    w = ri.pos(:,okk) - ri.pos(:,1);  w = w - dc*(dc.'*w);
    row(end+1) = sprintf('   %-5s %4d of %4d rays, beam radius %6.2f mm', names{j}, nnz(okk), s.nRays, max(vecnorm(w))*1e3); %#ok<AGROW>
    if j == 2
        [~, ~, V] = svd(w.', 'econ');  q = V(:,1:2).' * w;
        ell = sort(max(q, [], 2) - min(q, [], 2), 'descend') * 1e3;   % DM1 footprint, mm
        np1 = nnz(okk);
    end
    if j == numel(ie), npF = nnz(okk); nR = s.nRays; end
end
macos.trace(nF);  W = macos.opd();  M = macos.opd_mask();
wfe = std(W(M)) / (P.WAVLEN*mm);

% ---- the pupil walk across the coronagraph's field ---------------------------------------------
RS = [EDGE_LD*ld 0.5 1 3 10 30 60];                            % ring radii, arcsec
wk = walk_(fileread(fc), dB, L0, RS, [ic.DM1 ic.DM2], work);

% ---- the fold: e2e stage 8's recorded pick, replayed on this thread's telescope -----------------
% The same chain by the library's builder (ctb_fold_build: example_ctb's block via
% ctb_template_params), unfolded first as a check against the hand build above, then folded by the
% stage's pick.  The measure is the stage's (ctb_fold_measure: rays, wavefront, AOI, the walk on the
% ring, deck_clearance over the coronagraph's own field, the hardware union about the launch axis).
addpath(fullfile(root, 'templates', '80_end_to_end', 'e2e'));   % e2e_params: P.coro (the pick), the hole
PE = e2e_params();  Q = PE.coro;  rp = Q.pick;
holes = struct('M1', max(PE.hole_min_r_m, 0.232));           % the primary's through-hole: a crossing inside it is the design
Pc = ctb_template_params();
Kc = struct('dm1_m', LDM1, 'fno', fno, 'sign0', -1, 'model', 256);   % case C in the chapter's zigzag (first turn -)
Rc = ctb_fold_build(tel, fullfile(work, 'libC'), Pc, Kc);
Mc = ctb_fold_measure(Rc, 'holes', holes, 'leg_from', 3, 'ring_as', Q.ring_as, 'clear_as', Q.clear_as, 'quiet', true);
assert(Mc.rays(end) == npF && abs(Mc.wfe_waves - wfe) < 1e-6, 'ch02b: the library''s case C is not the hand build''s');
Kp = struct('fold0_m', rp.fold_m, 'fold0_tilt_deg', rp.tilt_deg, 'plane_deg', Q.plane_deg, 'sign0', rp.sign, ...
            'flip_at', rp.flip_at, 'relay_scale', rp.scale, 'dm1_m', LDM1, 'fno', fno, 'model', 256);
Rf = ctb_fold_build(tel, fullfile(work, 'fold'), Pc, Kp);
Mf = ctb_fold_measure(Rf, 'holes', holes, 'leg_from', 3, 'ring_as', Q.ring_as, 'clear_as', Q.clear_as, 'quiet', true);
sh = shroud_deck(tel, struct('shroud_D_m', 0, 'model', 256), 'model', 256);   % the telescope's own hardware union
Rgate = sh.D/2 + Q.shroud_grow_m;
margin = Rgate - Mf.shroud.ctb_max_r;
zc = [Mf.clear.bodies([Mf.clear.bodies.elt] > Rf.n_base).z];
wdev = Mf.walk.dm2_m/Mc.walk.dm2_m - 1;
e = @(nm) Rf.ix.(nm) - Rf.n_base;                               % element numbers in the bench alone

% the figures: the folded train (view_std) and the folded bench in world views
macos.init(256);  macos.load_rx(Rf.full_in);
fvw = macos.view_std('title', 'The stage-2 telescope with the CTB folded behind it', 'visible', false);
print(fvw, fullfile(work, 'ctb_coupled_view.png'), '-dpng', '-r150');  close(fvw);
fvb = figure('Visible', 'off', 'Position', [100 100 1600 640], 'Color', 'w');
ax1 = subplot(1, 2, 1);  macos.load_rx(Rf.bench_in);
macos.view_rx('ax', ax1, 'view', [0 0], 'title', 'the folded CTB alone: elevation, z up (the primary''s vertex at z = 0)');
xlabel(ax1, 'x  [m]');  zlabel(ax1, 'z  [m]');
ax2 = subplot(1, 2, 2);  macos.load_rx(Rf.full_in);
macos.view_rx('ax', ax2, 'view', [0 90], 'title', sprintf('plan, looking down z: inside the %.1f m envelope (dashed)', 2*Rgate));
hold(ax2, 'on');  tq = linspace(0, 2*pi, 361);
plot3(ax2, Rgate*cos(tq), Rgate*sin(tq), zeros(size(tq)), 'k--', 'LineWidth', 1.2);
xlim(ax2, 1.1*Rgate*[-1 1]);  ylim(ax2, 1.1*Rgate*[-1 1]);  xlabel(ax2, 'x  [m]');  ylabel(ax2, 'y  [m]');
print(fvb, fullfile(work, 'ctb_bench_view.png'), '-dpng', '-r150');  close(fvb);

% ---- the listing ---------------------------------------------------------------------------
t = strings(0);
t(end+1) = "the telescope (s2_instrument.in), as the coronagraph sees it";
t(end+1) = sprintf('   %d elements, %d of %d rays; field center %+.4f'' (the chief ray''s offset)', nT, nnz(ok), sT.nRays, bias);
t(end+1) = sprintf('   marginal ray f/%.2f; along the two pupil axes f/%.2f x f/%.2f (EFL %.2f / %.2f m)', fno, efl/Dtel, efl);
t(end+1) = sprintf('   lambda/D %.2f mas; the CTB''s %d lambda/D dark-zone edge = %.3f"', ld*1e3, EDGE_LD, EDGE_LD*ld);
t(end+1) = "";
t(end+1) = sprintf('the CTB as built (OAP1 f = %.3f m): %d of %d rays reach DM1', P.F_OAP(1)*mm, okA, nA);
t(end+1) = sprintf('OAP1 matched, f = f/%.2f x %.3f mm = %.4f m', fno, Dpup/mm, fOAP1);
t(end+1) = sprintf('   the telescope stop images %.4f m after OAP1 (first order %.4f m); DM1 there, DM2 %.3f m past it', ...
                   s_img, fOAP1*(fOAP1 + xpd)/xpd, P.L_DM2*mm);
t = [t, row];
t(end+1) = sprintf('   DM1 beam %.2f x %.2f mm on the %d-actuator grid (%.3f mm across, pitch %.3f mm; mirror clear aperture %.1f mm)', ell, nact, beamd, pitch/mm, 2*P.R_DM);
t(end+1) = sprintf('   FPA rms wavefront %.4f waves at %.0f nm', wfe, P.WAVLEN*1e6);
t(end+1) = "";
t(end+1) = "pupil walk, 8 field points per ring about the field center (max footprint-centroid walk)";
t(end+1) = "   ring      DM1 mm  /r_pup   DM2 mm  /r_pup  /pitch   rays DM2      rays FPA";
for k = 1:numel(RS)
    t(end+1) = sprintf('   %6.3f"  %7.3f  %6.3f  %7.3f  %6.3f  %6.2f   %4d..%4d   %4d..%4d', RS(k), ...
        wk.walk(k,1)*1e3, wk.walk(k,1)/rpup, wk.walk(k,2)*1e3, wk.walk(k,2)/rpup, wk.walk(k,2)/pitch, ...
        wk.nmin(k,2), wk.nmax(k,2), wk.fmin(k), wk.fmax(k)); %#ok<AGROW>
end
mm_transcript('ch02_ctb_run', strjoin(t, newline));

t = strings(0);
t(end+1) = "the fold: e2e stage 8's pick (e2e_params P.coro.pick) on this telescope, by ctb_fold_build / ctb_fold_measure";
t(end+1) = sprintf('   F0, a flat %.2f m before the focus, deviation %g deg (AOI %.1f) out of the bench plane toward +z;', rp.fold_m, rp.tilt_deg, rp.tilt_deg/2);
t(end+1) = sprintf('   the accordion''s plane %g deg from the bench plane, first turn %+d, the pattern inverted after turn %d; relay scale %.2f', Q.plane_deg, rp.sign, rp.flip_at, rp.scale);
t(end+1) = sprintf('   unfolded (case C, the library''s build): CTB bodies to r = %.3f m; DM2 walk at %g" %.3f mm', Mc.shroud.ctb_max_r, Q.ring_as, Mc.walk.dm2_m*1e3);
for j = 1:numel(Mf.stations)
    t(end+1) = sprintf('   %-5s %4d of %4d rays, beam radius %6.2f mm', Mf.stations{j}, Mf.rays(j), Mf.nrays, Mf.beam_r(j)*1e3); %#ok<AGROW>
end
t(end+1) = sprintf('   FPA rms wavefront %.4f waves (unfolded %.4f); lambda/D at the FPM %.1f um (unfolded %.1f)', Mf.wfe_waves, Mc.wfe_waves, Rf.lamD_fpm_um, Rc.lamD_fpm_um);
a = "   AOI (deg):";  for k = 1:numel(Mf.aoi), a = a + sprintf(' %s %.1f', Mf.aoi(k).name, Mf.aoi(k).aoi_deg); end;  t(end+1) = a;
t(end+1) = sprintf('   walk at %g": DM2 %.3f mm (unfolded %.3f, %+.2f %%); rays DM2 %d..%d, FPA %d..%d', Q.ring_as, Mf.walk.dm2_m*1e3, Mc.walk.dm2_m*1e3, 100*wdev, Mf.walk.rays_dm2, Mf.walk.rays_fpa);
t(end+1) = sprintf('   envelope: the telescope''s own hardware union %.3f m = %.2f x D, gate radius %.3f m; CTB bodies to r = %.3f m, margin %+.3f m', sh.D, sh.D/Dtel, Rgate, Mf.shroud.ctb_max_r, margin);
t(end+1) = sprintf('   train %.3f m (the telescope alone %.3f); CTB bodies from z = %.3f to %.3f m behind the primary''s vertex (the bench at %.2f)', Mf.shroud.len, sh.len, Mf.zmin_ctb, max(zc), Q.z_floor_m);
t(end+1) = sprintf('   clearance (mm, signed): CTB vs telescope %+.1f (leg from %s vs %s); within the CTB %+.1f (%s vs %s), bodies and legs over the %g" field;', ...
                   Mf.clear.min*1e3, Mf.clear.worst.leg_name, Mf.clear.worst.body_name, Mf.clear.min_addon*1e3, Mf.clear.worst_addon.leg_name, Mf.clear.worst_addon.body_name, Q.clear_as);
t(end+1) = sprintf('   the telescope''s own worst pair %+.1f (%s vs %s), as built', Mf.clear.min_base*1e3, Mf.clear.worst_base.leg_name, Mf.clear.worst_base.body_name);
t(end+1) = sprintf('   in the bench alone: E%d F0, E%d OAP1, E%d DM1, E%d DM2, E%d FPM, E%d Lyot, E%d FPA', e('F0'), e('OAP1'), e('DM1'), e('DM2'), e('FPM'), e('Lyot'), e('FPA'));
mm_transcript('ch02_ctb_fold', strjoin(t, newline));
mm_figure('ch02_ctb_view', fullfile(work, 'ctb_coupled_view.png'));
mm_figure('ch02_ctb_bench', fullfile(work, 'ctb_bench_view.png'));

mm_pin('ch02_ctb_fno',     fno,        'sig5', sprintf('f/%.2f', fno));
mm_pin('ch02_ctb_fno_u',   efl(1)/Dtel, 'sig5', sprintf('f/%.2f', efl(1)/Dtel));
mm_pin('ch02_ctb_fno_v',   efl(2)/Dtel, 'sig5', sprintf('f/%.2f', efl(2)/Dtel));
mm_pin('ch02_ctb_bias',    bias,       'sig5', sprintf('%.2f′', bias));
mm_pin('ch02_ctb_asbuilt', okA,        'sig5', sprintf('%d', okA));
mm_pin('ch02_ctb_foap1',   fOAP1,      'sig5', sprintf('%.3f m', fOAP1));
mm_pin('ch02_ctb_simg',    s_img,      'sig5', sprintf('%.1f mm', s_img*1e3));
mm_pin('ch02_ctb_pass',    npF,        'sig5', sprintf('%d', npF));
mm_pin('ch02_ctb_rays',    nR,         'sig5', sprintf('%d', nR));
mm_pin('ch02_ctb_ell_a',   ell(1),     'sig5', sprintf('%.1f', ell(1)));
mm_pin('ch02_ctb_ell_b',   ell(2),     'sig5', sprintf('%.1f', ell(2)));
mm_pin('ch02_ctb_wfe',     wfe,        'sig5', sprintf('%.4f waves', wfe));
mm_pin('ch02_ctb_edge',    EDGE_LD*ld, 'sig5', sprintf('%.2f″', EDGE_LD*ld));
mm_pin('ch02_ctb_w2_edge', wk.walk(1,2)/pitch, 'sig5', sprintf('%.2f', wk.walk(1,2)/pitch));
mm_pin('ch02_ctb_w1_edge', wk.walk(1,1)*1e3,   'sig5', sprintf('%.3f mm', wk.walk(1,1)*1e3));
mm_pin('ch02_ctb_w2_10',   wk.walk(5,2)/rpup,  'sig5', sprintf('%.0f%%', 100*wk.walk(5,2)/rpup));
mm_pin('ch02_ctb_mag',     Dtel/Dpup,  'sig5', sprintf('%.0f', Dtel/Dpup));
mm_pin('ch02_ctb_np1',     np1,        'sig5', sprintf('%d', np1));
mm_pin('ch02_fold_reach',  Mc.shroud.ctb_max_r, 'sig5', sprintf('%.1f m', Mc.shroud.ctb_max_r));
mm_pin('ch02_fold_gate',   2*Rgate,    'sig5', sprintf('%.1f m', 2*Rgate));
mm_pin('ch02_fold_f0_m',   rp.fold_m,  'sig5', sprintf('%.2f m', rp.fold_m));
mm_pin('ch02_fold_f0_aoi', rp.tilt_deg/2, 'sig5', sprintf('%.0f°', rp.tilt_deg/2));
mm_pin('ch02_fold_rays',   Mf.rays(end), 'sig5', sprintf('%d', Mf.rays(end)));
mm_pin('ch02_fold_wfe',    Mf.wfe_waves, 'sig5', sprintf('%.4f waves', Mf.wfe_waves));
mm_pin('ch02_fold_wfe_dev', Mf.wfe_waves - Mc.wfe_waves, sprintf('abs:%g', Q.wfe_tol_waves), sprintf('< %g waves', Q.wfe_tol_waves));
mm_pin('ch02_fold_walk10', Mf.walk.dm2_m*1e3, 'sig5', sprintf('%.2f mm', Mf.walk.dm2_m*1e3));
mm_pin('ch02_fold_walk_dev', wdev,     'abs:0.01', '< 1 %');
mm_pin('ch02_fold_aoi_max', max([Mf.aoi.aoi_deg]), 'sig5', sprintf('%.0f°', max([Mf.aoi.aoi_deg])));
mm_pin('ch02_fold_zmax',   max(zc),    'sig5', sprintf('%.2f m', max(zc)));
mm_pin('ch02_fold_margin', margin,     'sig5', sprintf('%.3f m', margin));
mm_pin('ch02_fold_clear',  Mf.clear.min*1e3, 'sig5', sprintf('%.0f mm', Mf.clear.min*1e3));
mm_pin('ch02_fold_own',    Mf.clear.min_addon*1e3, 'sig5', sprintf('%.0f mm', Mf.clear.min_addon*1e3));
mm_pin('ch02_fold_fm',     Mf.clear.min_base*1e3, 'sig5', sprintf('%.0f mm', Mf.clear.min_base*1e3));
mm_pin('ch02_fold_len',    Mf.shroud.len, 'sig5', sprintf('%.2f m', Mf.shroud.len));
mm_pin('ch02_fold_tel_len', sh.len,    'sig5', sprintf('%.2f m', sh.len));

% =================================================================================================
function [full_in, ix] = bench_(P, f1, LDM1, pF, dF, fno, nd, tel, work, tag)
% example_ctb's build_ctb chain (P in mm -> m), seeded 0.30 m before the telescope focus, OAP1 f = F1,
% OAP1 -> DM1 = LDM1; spliced onto the telescope (its last three elements, FP_return/ExitPupil/FP, dropped)
    mm = 1e-3;  Fo = P.F_OAP*mm;  Fo(1) = f1;  r = Fo / cosd(P.AOI)^2;  back = 0.30;
    Ra = P.R_OAP*mm;  Rd = P.R_DM*mm;
    b = macos.design.Bench(['ctb_' tag], 'baseunits','m', 'pos', pF - back*dF, 'dir', dF, ...
            'wavelen', P.WAVLEN*mm, 'aperture', 2*atan(1/(2*fno)), 'ngridpts', 63, 'zsource', back);
    k = 0;  d0 = dF;
    k=k+1; O1 = b.add_oap(back + r(1), nd(d0,k), 'mode','collimate','f',Fo(1),'name','OAP1','aprad',Ra); d0=b.dir;
    k=k+1; iDM1 = b.add_mirror(LDM1,'out',nd(d0,k),'name','DM1','aprad',Rd); d0=b.dir;
    k=k+1; iDM2 = b.add_mirror(P.L_DM2*mm,'out',nd(d0,k),'name','DM2','aprad',Rd); d0=b.dir;
    k=k+1; b.add_oap(P.L_O2*mm, nd(d0,k),'mode','focus','f',Fo(2),'name','OAP2','aprad',Ra); d0=b.dir;
    b.add_reference(r(2), 'Focus23');
    k=k+1; b.add_oap(r(3), nd(d0,k),'mode','collimate','f',Fo(3),'name','OAP3','aprad',Ra); d0=b.dir;
    b.add_reference(P.L_APOD*mm, 'Apodizer');
    k=k+1; b.add_oap(P.L_O4*mm, nd(d0,k),'mode','focus','f',Fo(4),'name','OAP4','aprad',Ra); d0=b.dir;
    b.add_reference(r(4), 'FPM');
    k=k+1; b.add_oap(r(5), nd(d0,k),'mode','collimate','f',Fo(5),'name','OAP5','aprad',Ra); d0=b.dir;
    iLyot = b.add_reference(P.L_LYOT*mm, 'Lyot');
    k=k+1; b.add_oap(P.L_O6*mm, nd(d0,k),'mode','focus','f',Fo(6),'name','OAP6','aprad',Ra); d0=b.dir;
    b.add_reference(r(6), 'FieldStop');
    k=k+1; b.add_oap(r(7), nd(d0,k),'mode','collimate','f',Fo(7),'name','OAP7','aprad',Ra); d0=b.dir;
    b.add_reference(P.L_BACK*mm, 'Backend');
    k=k+1; b.add_oap(P.L_O8*mm, nd(d0,k),'mode','focus','f',Fo(8),'name','OAP8','aprad',Ra);
    b.add_detector(r(8), 'FPA');
    bench_in = fullfile(work, ['ctb_' tag '.in']);  b.emit(bench_in);
    full_in  = fullfile(work, ['e2e_ctb_' tag '.in']);
    info = macos.design.append_rx(tel, bench_in, full_in, 'drop_base_tail', 3);
    ix = struct('OAP1', info.n_base + O1.i, 'DM1', info.n_base + iDM1, 'DM2', info.n_base + iDM2, ...
                'Lyot', info.n_base + iLyot);
end

function f = field_deck_(txt, d0, L0, dxy, work, k)
% the coupled deck with its chief ray turned by DXY (rad) about the field center
    d = d0 + [dxy(1); dxy(2); 0];  d = d/norm(d);
    t2 = regexprep(txt, '(?m)^(\s*ChfRayDir=)[^\n]*', sprintf('$1  %.16E  %.16E  %.16E', d));
    t2 = regexprep(t2, '(?m)^(\s*ChfRayPos=)[^\n]*', sprintf('$1  %.16E  %.16E  %.16E', -L0*d));
    f = fullfile(work, sprintf('field_%02d.in', k));
    fid = fopen(f, 'w');  fwrite(fid, t2);  fclose(fid);
end

function s = stop_image_(txt, d0, L0, iO1, ring_as, work)
% where the field chief rays leaving OAP1 cross the field-center chief: the telescope stop's image
    F = [0 0; macos.design.field_ring(ring_as/60, 'units', 'arcmin', 'n', 8, 'inner', 0)];
    p = zeros(3, size(F,1));  q = p;
    for k = 1:size(F,1)
        macos.init(256);  macos.load_rx(field_deck_(txt, d0, L0, F(k,:), work, k));
        ri = macos.get_ray_info(macos.trace(iO1).nRays);
        p(:,k) = ri.pos(:,1);  q(:,k) = ri.dir(:,1)/norm(ri.dir(:,1));
    end
    sx = zeros(1, size(F,1) - 1);
    for k = 2:size(F,1)
        w0 = p(:,1) - p(:,k);  b = q(:,1).'*q(:,k);
        sx(k-1) = (b*(q(:,k).'*w0) - q(:,1).'*w0) / (1 - b^2);
    end
    s = median(sx);
end

function wk = walk_(txt, d0, L0, RS, eD, work)
% footprint-centroid walk on the elements ED across 8-point rings of radius RS (arcsec)
    c0 = zeros(3, numel(eD));
    wk.walk = zeros(numel(RS), numel(eD));  wk.nmin = wk.walk;  wk.nmax = wk.walk;
    wk.fmin = zeros(numel(RS), 1);  wk.fmax = wk.fmin;
    for kr = 0:numel(RS)
        if kr == 0, F = [0 0]; else, F = macos.design.field_ring(RS(kr)/60, 'units', 'arcmin', 'n', 8, 'inner', 0); end
        nk = zeros(size(F,1), numel(eD));  nf = zeros(size(F,1), 1);  wl = nk;
        for k = 1:size(F,1)
            macos.init(256);  n = macos.load_rx(field_deck_(txt, d0, L0, F(k,:), work, k));
            for j = 1:numel(eD)
                rj = macos.get_ray_info(macos.trace(eD(j)).nRays);  okj = rj.ok_pass & rj.ok_trace;
                nk(k,j) = nnz(okj);
                if nk(k,j) == 0, wl(k,j) = NaN; continue; end
                c = mean(rj.pos(:,okj), 2);
                if kr == 0, c0(:,j) = c; end
                wl(k,j) = norm(c - c0(:,j));
            end
            rn = macos.get_ray_info(macos.trace(n).nRays);  nf(k) = nnz(rn.ok_pass & rn.ok_trace);
        end
        if kr == 0, continue; end
        wk.walk(kr,:) = max(wl, [], 1);  wk.nmin(kr,:) = min(nk, [], 1);  wk.nmax(kr,:) = max(nk, [], 1);
        wk.fmin(kr) = min(nf);  wk.fmax(kr) = max(nf);
    end
end
