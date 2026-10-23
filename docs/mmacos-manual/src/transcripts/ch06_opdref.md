    >> macos.init(512);
    >> n = macos.load_rx('e2e_pie_metopt.in');  macos.opd_ref('chief');
    >> macos.trace(n-1);  W0 = macos.opd();  macos.perturb(1, 'translation', [0;0;1e-7]);  macos.modify();
    >> macos.trace(n-1);  d = macos.opd() - W0;
       chief  rms 188.9 nm   mean 178.4 nm   rms about the mean 61.98 nm
    >> n = macos.load_rx('e2e_pie_metopt.in');  macos.opd_ref('mean');
    >> macos.trace(n-1);  W0 = macos.opd();  macos.perturb(1, 'translation', [0;0;1e-7]);  macos.modify();
    >> macos.trace(n-1);  d = macos.opd() - W0;
       mean   rms 61.98 nm   mean < 1e-09 nm   rms about the mean 61.98 nm
    >> the Section 4 column (Elt 1 Tz, center field, harvested at 'mean') x 100 nm:
       linear rms 61.98 nm
