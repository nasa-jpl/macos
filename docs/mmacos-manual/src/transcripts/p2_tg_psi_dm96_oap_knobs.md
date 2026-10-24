    % ---- the bench (macos.design.twyman_green options; s = 96/56 applied
    %      in tg96_run so the whole rig scales uniformly off the 56 mm v1) --
    P.bench.polarizing = true;
    P.bench.BS_AOI     = 22.5;         % Dave 2026-09-15: pinned (the Stage-A solve's 7 deg cleared only the
                                       % end bodies; the node parts need >= 22.5: dmg_bench_clearance)
    P.bench.D_RECOMB   = 150;          % physical mm (NOT scaled by s): the recomb plane and the output
    P.bench.D_RC_L2    = 55;           % optics 150 mm behind the splitter, L2 at 150 + 55 = 205 as before
    P.bench.F1 = 500;   P.bench.F2 = 250;      % *s in the runner
