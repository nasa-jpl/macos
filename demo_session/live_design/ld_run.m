% ld_run.m -- LIVE DESIGN BEAT (rehearsal + demo runner)
% ===================================================================
%  An audience-spec two-mirror telescope, designed on the spot.
%  The audience gives THREE numbers; the presenter supplies the two
%  the audience didn't know they owed (primary speed, back focus) --
%  "the engineer's two knobs" -- and the design layer does the rest:
%  closed-form seed -> emitted Rx -> SMACOS-validated -> native
%  multi-field conic solve -> layout drawn.
%
%  Run:  matlab -batch "run('ld_run.m')"   from this directory.
% ===================================================================

% ---- the audience's three numbers --------------------------------
D_M           = 0.33;    % aperture diameter, metres
SYS_FNUM      = 16;      % system f/#
FULLFIELD_DEG = 0.5;     % full field of view, degrees

% ---- the presenter's two knobs (+ fixed choices) -----------------
FAMILY   = 'RC';         % aplanat: null spherical AND coma
PRI_FNUM = 2.0;          % primary speed (packaging: f1 = 0.66 m)
BFD_MM   = 300;          % back focal distance behind the primary vertex
MODEL    = 256; WVL = 633e-9;

t_all = tic;
addpath('~/dev/MACOS_res_dev/mmacos/src');
half_arcmin = FULLFIELD_DEG*60/2;              % 15' at 0.5 deg full field
dl_rms_m    = WVL/14;                          % Marechal DL criterion

fprintf('==================================================================\n');
fprintf(' LIVE DESIGN: %.2f m f/%g %s, %.2f deg full field\n', ...
        D_M, SYS_FNUM, FAMILY, FULLFIELD_DEG);
fprintf('==================================================================\n');

% -- Stage 1: intent -> closed-form design -> validated Rx ---------
t1 = tic;
t = macos.design.Telescope('family',FAMILY, ...
        'aperture_diameter_m',D_M, 'primary_fnum',PRI_FNUM, ...
        'system_fnum',SYS_FNUM, 'BFD_mm',BFD_MM, ...
        'model_size',MODEL, 'wavelength_m',WVL);
rx = t.build();
t.describe();
fprintf('  [stage 1: seed+build+validate  %.1f s]\n', toc(t1));

% -- Stage 2: multi-field conic solve (structure held fixed) -------
t2 = tic;
res = t.optimize('fields_arcmin', [half_arcmin/2, half_arcmin]);
fprintf('\n  field    before        after      (RMS WFE, nm; DL = %.1f nm)\n', dl_rms_m*1e9);
labels = {sprintf('%5.1f''',0), sprintf('%5.1f''',half_arcmin/2), sprintf('%5.1f''',half_arcmin)};
for k = 1:numel(res.wfe_before)
    verdict = 'DL'; if res.wfe_after(k) > dl_rms_m, verdict = 'NOT DL'; end
    fprintf('  %s   %9.2f   %9.2f    %s\n', labels{k}, ...
            res.wfe_before(k)*1e9, res.wfe_after(k)*1e9, verdict);
end
fprintf('  conics: K1 %.6f  K2 %.6f (seed K1 %.6f K2 %.6f)\n', ...
        res.conics(1), res.conics(2), t.spec.derived.K1, t.spec.derived.K2);
fprintf('  [stage 2: multi-field solve  %.1f s]\n', toc(t2));

% -- Stage 3: draw the design they named ---------------------------
t3 = tic;
t.view_layout('YZ', 'save','ld_layout.png', 'visible',false);
fprintf('  [stage 3: layout figure -> ld_layout.png  %.1f s]\n', toc(t3));

% -- Stage 4: the wall, answered -- add a third mirror -------------
%  Two conics null spherical + coma; NOTHING is left for field
%  astigmatism (~field^2).  A Korsch TMA's third conic supplies
%  exactly that DOF.  Same Cassegrain pair the audience's numbers
%  derived, plus a 1:1 tertiary relay (preserves f/16), Seidel-seeded.
t4 = tic;
d = t.spec.derived;
fprintf('\n Stage 4 -- the third mirror (Korsch TMA on the same front end)\n');
tm = macos.design.Telescope('family','TMA', ...
        'aperture_diameter_m',D_M, 'model_size',MODEL, 'wavelength_m',WVL);
tm.add_mirror('M1','radius_m',abs(d.R1), 'spacing_after_m',abs(d.sep));
tm.add_mirror('M2','radius_m',abs(d.R2), 'spacing_after_m',abs(d.sep)+BFD_MM/1e3+0.4);
tm.add_mirror('M3','radius_m',0.4, 'spacing_after','derive');   % 1:1 relay of the Cass focus
tm.build();
res3 = tm.optimize('fields_arcmin', [half_arcmin/2, half_arcmin]);
fprintf('\n  field    RC (2 conics)   TMA (3 conics)   (RMS WFE, nm; DL = %.1f nm)\n', dl_rms_m*1e9);
for k = 1:numel(res3.wfe_after)
    v2 = 'NOT DL'; if res.wfe_after(k)  <= dl_rms_m, v2 = 'DL'; end
    v3 = 'NOT DL'; if res3.wfe_after(k) <= dl_rms_m, v3 = 'DL'; end
    fprintf('  %s    %9.2f %-7s %9.2f %-7s\n', labels{k}, ...
            res.wfe_after(k)*1e9, v2, res3.wfe_after(k)*1e9, v3);
end
fprintf('  TMA conics: K1 %.4f  K2 %.4f  K3 %.4f\n', res3.conics(1:3));
fprintf('  [stage 4: TMA build+solve  %.1f s]\n', toc(t4));

t5 = tic;
tm.view_layout('YZ', 'save','ld_layout_tma.png', 'visible',false);
fprintf('  [stage 5: TMA layout -> ld_layout_tma.png  %.1f s]\n', toc(t5));

% -- Stage 6: keep both --------------------------------------------
t.save('ld_design.in');      t.save_spec('ld_design.mat');
tm.save('ld_design_tma.in'); tm.save_spec('ld_design_tma.mat');
fprintf('\n  saved: ld_design{,_tma}.in + .mat\n');
fprintf('  TOTAL %.1f s\n', toc(t_all));
