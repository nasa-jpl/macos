% ch08a_contrast.m -- Chapter 8, the dark zone before control: templates/30_instruments/bench_ctb's
% ctb_mask_compare, run unchanged on a scratch copy of the template (mm_template).  Six literature
% focal-plane-mask families on the same committed CTB deck (ctb_dcr.in, its reference spheres as
% committed -- PLAN 1.2 (zzz)), the same 1024 grid, the same 3-15 lambda/D annulus and the same
% normalization (the bare on-axis peak).  Pins: each family's dark-zone mean and throughput proxy.

[here, check] = mm_template('30_instruments', 'bench_ctb');
work = getenv('MMMAN_WORK');
t0 = tic;
txt = mm_clean(evalc('out = ctb_mask_compare(''outdir'', work);'));
fprintf('ctb_mask_compare wall time %.0f s\n', toc(t0));
check();
disp(txt)

k0 = regexp(txt, '\[cmp\] N=', 'once');
k1 = regexp(txt, 'throughput = off-axis proxy\)', 'end', 'once');
mm_transcript('ch08_masks_run', strtrim(txt(k0:k1)));
mm_figure('ch08_masks', fullfile(work, 'ctb_mask_compare.png'));

key = {'hard', 'blc', 'aplc', 'vortex', 'roddier', 'dualzone'};
lab = {'hard occulter', 'band-limited', 'APLC', 'vortex \(matched\)', 'Roddier', 'dual-zone'};
for k = 1:numel(key)
    t = regexp(txt, [lab{k} '[^|]*\|\s*(\S+)\s*\|\s*\S+\s*\|\s*\S+\s*\|\s*(\S+)%'], 'tokens', 'once');
    assert(~isempty(t), 'ch08a: no table row for %s', lab{k});
    v = str2double(t{1});  th = str2double(t{2});
    mm_pin(['ch08_dz_' key{k}], v, 'sig5', strrep(sprintf('%.1e', v), 'e-0', 'e-'));
    mm_pin(['ch08_th_' key{k}], th, 'sig5', sprintf('%.0f%%', th));
end
