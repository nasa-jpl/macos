    P.tag = 'cophi';
    P.outdir = '';                                   % '' = <this dir>/runs/<tag>
    P.stages = {'bench'};                            % bench | battery | noise | loop | capture | format | figs (+ tg96_run's figs / clearance)
    P.bench.ref_arm = 'lo';                          % the flat's seat -> the traced LO arm (twyman_green 'ref_arm','lo'; ruling 1)
    P.cophi.record_tail = struct('lens', 'redo96_lens_tail.mat', 'oap', 'redo96_oap_tail.mat');   % the deck's IFO rows ran on these (RUNBOOK_mac_cycle3b)
    P.cophi.oap_record = true;                       % 'bench.optics','oap' => the MIRROR RIG OF RECORD (POL_IN source, AOIs 20/25, sides +1/-1, D_RC_L2 125)
    P.bench.LO_DX   = [0;0;0];                       % fiber-tip offset, global mm (the LO drift knob: lateral = tilt d/F1, axial = defocus)
