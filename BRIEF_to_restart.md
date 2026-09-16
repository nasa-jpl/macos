# BRIEF for TO after the restart (2026-09-16)

CCL for TO (Dave relays).  Branch `dev-candidate`, both repos, shared tree on
this 32 GB box.  Both repos were PUSHED at 08:45 (macos 3f50c3d, resources
d28fdeb); everything of yours and mine is on origin.  The box is idle: no
MATLAB, no sequencer, no watcher.  You are a fresh session; nothing has to be
recovered from memory, only read.

## Read, in this order (10 minutes)

1. `BRIEF_to_gauge_close.md` section 0 (the rules) and section 0b (the CTB
   retraction: the DM model is RIGHT, your probe settled it; the real finding
   is the point-source Aperture convention; Dave's ruling: docs now -- DONE by
   CCL, resources d28fdeb -- and regeneration at the intended 42.75 mm beam
   LATER, after the gauge queue below is drained).
2. The status table at the top of `tg_psi_dm96_oap/REPORT_reflective.md`
   (your own; items 0, 1, 2(b), 5, 6 done; 2(a), 3, 4, 7 open) and its
   "What is running" table (everything through the CTB probe done; the chain
   stopped at the one manual step).
3. `runs/*.log` tails for exit codes only if something looks off.  Do NOT
   re-read sections 4.1-4.4 of the report (retracted attributions).

## The queue, in order, each item whole and committed before the next

| # | step | closes | wall time |
|---|---|---|---|
| 1 | apply the staged patch: `cd runs && python3 apply_tg96_pending.py` (the wrap stage, the camera line, `MASK_SUB`, the stations width; dry-run clean against copies) | -- | 2 min |
| 2 | `runs/item2bseq.sh` (detached): the `wrap` stage on both rigs (`wrapoap`, `wraplens`: dA vs dD, n_cross) + both rigs' station figures at 1800 px | item 2(a): the 120 / 240 nm mechanism; item 6's lens-rig figure | ~45 min |
| 3 | `runs/item4bseq.sh` (detached): thk22's retune WITH the advisory gate + its gate run | item 4's interferometer rows with the 10 mm plates (the 20.09 nm null is the tuner's; the ROWS are what goes on a slide) | ~1 h |
| 4 | item 3: the tuner's verify measure over the lattice (`dmg_act_fit`, not `interp2` at one pixel); the same two legs as the test (objwin3 must be refused, lens_tail accepted); README line | the detector-leg tuner | ~2 h of you, 20 min of runs |
| 5 | item 7 steps 1-3 as written in `BRIEF_to_gauge_close.md` section 7 and your `bench_ctb/REPORT_field_servo.md` (step 1's prescription is written: dichroic at Apodizer_Pst, 300 mm lens F/9.4, 2 lam/D dimple 11.8 um; step 2's expectation on the traced beam: 47% conversion at the actuator Nyquist, 50% at 16.5 cycles -- the CTB's DMs separate at the top of the band) | the coronagraph field servo | 2-3 days |
| 6 | Dave's CTB item 2: regenerate at the intended 42.75 mm beam (`Aperture = 2*NA` in example_ctb, regenerate the decks through the generator, then `ctb_study` re-derives jac/efc/relin/physics/bandwidth/vvc) -- ONLY after 1-5, and Dave says when | the CTB at the sheet's fill | ~1 day of one MATLAB + a day of you |

Steps 2 and 3 are sequence scripts already in the tree; chain them with a
waiter on the wrapper's `] exit` line (your `closefinal.sh` pattern), never on
a sequencer's echo, and never edit a script that is running.

## What CCL did with your overnight results (so you do not redo it)

- Gauge deck: 47 slides (macos 6937893): your servo, descent, wrap, vector
  verdict, pinhole, polarization-at-the-built-angle, camera and thickness
  numbers are on the slides; a new "Interferometer, station by station"
  slide uses `oapifol2_stations.png` as it is; a new "The reflective front
  end, measured" slide carries the record / redesigned / lens table.  When
  steps 2-3 land, send CCL the two lines (the rows with the 10 mm plates;
  the wrap mechanism's verdict) and the 1800 px station figures.
- CTB docs: README (47% fill, the Aperture convention, the chain as
  diameters), `example_ctb.m` label, deck_ctb slides 1 and 9, a dated
  section in CTB_PROP_STATUS.md.  No generator or deck changed.
- Memory and briefs corrected to the traced beam; the 64x64 / 1 mm DM ruling
  is withdrawn with its premise.
- `dmg_loop` warns once above 1e-2 mm of start_rms (tDmgLoop 16/16); your
  sequencer rule is recorded.

## Rules that bit this week, restated

ONE model-1024 MATLAB on this box, every wrapper waits; kill by PID; no
`--amend`, no push (Dave pushes); commit per item; report as you go with the
status table on top; a document can be read wrongly -- trace it (the probe
was the right instrument, and it is why the retraction is clean).
