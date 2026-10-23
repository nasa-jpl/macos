% ch06c_opdref.m -- Chapter 6's check of its own comparison: one poke, segment 1 (elt 1) piston
% +100 nm, on the thread's e2e_pie_metopt.in (the deck run_compare uses), the center-field OPD change
% at the compare's wavefront element under each OPD reference, beside the Section 4 linear column
% (harvested at 'mean').  The reader's check of PLAN_CONSOLIDATION 1.2 (kkk).  Model size 512.
% Note the macos.modify() after the perturb: without it the next trace reuses the cached geometry
% and the change reads zero (macos.channels.RigidBodyChannel does the same).
% Signals (printed), except ch06_ref_mean_mean < 1e-9 nm: the mean-referenced map's own mean, a zero.

root = getenv('MMACOS_ROOT');
work = getenv('MMMAN_THREAD');
addpath(fullfile(root, 'src'));
rx = fullfile(work, 'e2e_pie_metopt.in');
macos.init(512);
L = {'>> macos.init(512);'};
for ref = ["chief" "mean"]
    n = macos.load_rx(rx);  macos.opd_ref(char(ref));
    macos.trace(n-1);  W0 = macos.opd();  m = W0 ~= 0;
    macos.perturb(1, 'translation', [0;0;1e-7]);  macos.modify();
    macos.trace(n-1);  W1 = macos.opd();
    d = (W1(m) - W0(m)) * 1e9;                              % nm (BaseUnits m)
    L{end+1} = sprintf('>> n = macos.load_rx(''e2e_pie_metopt.in'');  macos.opd_ref(''%s'');', ref); %#ok<SAGROW>
    L{end+1} = '>> macos.trace(n-1);  W0 = macos.opd();  macos.perturb(1, ''translation'', [0;0;1e-7]);  macos.modify();';
    L{end+1} = '>> macos.trace(n-1);  d = macos.opd() - W0;';
    L{end+1} = sprintf('   %-5s  rms %.4g nm   mean %.4g nm   rms about the mean %.4g nm', ...
        ref, sqrt(mean(d.^2)), mean(d), std(d, 1));
    if ref == "mean"                                    % the mean-referenced map's mean is a numerical zero
        L{end} = mm_zero('ch06_ref_mean_mean', L{end}, 'mean (\S+) nm   rms about', 1e-9);
    end
    mm_pin("ch06_ref_" + ref + "_rms", sqrt(mean(d.^2)), 'printed', sprintf('%.4g', sqrt(mean(d.^2))));
    if ref == "chief"
        mm_pin('ch06_ref_chief_mean', mean(d), 'printed', sprintf('%.4g', mean(d)));
    end
end
J = load(fullfile(work, 's4_jacobians.mat'), 'ox');  ox = J.ox;
i = find(strcmp(ox.channel_names, 'Elt 1 Tz'), 1);  icf = find(strcmp(ox.field_names, 'C'), 1);
c = ox.per_field_dwdx{icf}(:, i) * 1e-7 * 1e9;  c = c(isfinite(c));
L{end+1} = '>> the Section 4 column (Elt 1 Tz, center field, harvested at ''mean'') x 100 nm:';
L{end+1} = sprintf('   linear rms %.4g nm', sqrt(mean(c.^2)));
mm_pin('ch06_ref_linear_rms', sqrt(mean(c.^2)), 'printed', sprintf('%.4g', sqrt(mean(c.^2))));
mm_transcript('ch06_opdref', strjoin(string(L), newline));
