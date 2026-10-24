# BRIEF (TO, from a cleared state): why the redo bench's descent stopped converging

CCL for TO, 2026-09-18.  One question, one run, a verdict rule that separates two
named causes.  Nothing here touches the deck or the report -- CCL owns those this
session.  Do not push.

## Read-up (in this order, then stop reading and run)

1. `~/dev/macos/CLAUDE.md` (root rules) -- you are on the Linux box, engine tree
   `~/dev/macos`, resources `~/dev/MACOS_resources`, both on `dev-candidate`.
2. This brief.  You do NOT need `CURRENT_SLICE.md` for this job; the state you need
   is in section 2 below.
3. `~/dev/macos/RUNBOOK_mac_cycle3c.md` section "Lane M" -- the invocation this run
   is a variant of.  Read it for shape only; the exact command is in section 3 here.

Background you may want AFTER the run is launched, not before:
`MACOS_resources/mmacos/templates/40_benches/tg_psi_dm96_oap/REPORT_bench_realism.md`
section 8 (package A: what changed on this bench and why).

## 1. The question

The Twyman-Green bench was rebuilt in package A (collimator re-solved, substrates in,
the beam opened).  On the rebuilt bench the DM servo's **descent stopped converging**,
on BOTH rigs, with settings identical to the record's runs.  Two causes are on the
table and the evidence on disk cannot separate them:

- **(a) the operating point.**  The rebuilt bench's flat-DM null is 59 nm (lens) /
  74 nm (mirror) against the record's ~10 nm.  That static pattern sits on top of the
  descent's starting surface and may push the opening reading past the raw four-step
  fold (+-158 nm of surface).  If so, the descent is being asked to start outside its
  capture range and the bench's *gain* is fine.
- **(b) the control set.**  The 48 mm aperture on the DM now stops the beam: the
  calibration reports **3680 lit actuators** against the record's **7548**, with the
  per-actuator window 27 -> 37 px.  The Tikhonov weight `battery.matrix_lam` is
  unchanged at 1e-3 (relative to the median column energy of J).  If the inversion is
  now over-damped for the new lit set, the loop gain is genuinely lower and the fix is
  the regularization, not the start.

These predict opposite things at a SMALL start, which is what makes one run decisive.

## 2. The evidence (so you can recognize the failure when you see it)

Mirror rig, `loop.gain` 0.50, 1e13 photons/cycle, `recal never` -- identical settings
in both columns.  `rho` is the fitted per-cycle contraction (1 - gG), so rho 0.50 means
an effective loop gain G ~ 1.0 and rho 0.85 means G ~ 0.30.

| start | run | k(3pm) | r(K) pm | rho |
|---|---|---|---|---|
| 100 nm | record `runs/oapdesc2` | 17 | 2.34 | 0.502 |
| 200 nm | record `runs/oapdesc2` | 19 | 2.35 | 0.515 |
| 100 nm | redo `runs/redo96_oapdesc` | **never** | **105.6** | **0.851** |
| 200 nm | redo `runs/redo96_oapdesc` | **never** | **487.2** | 0.851 |

The loop OPENS at the same place (r(1) 69.9 vs 69.8 nm) and then crawls.  The same
bench's single-poke rows are healthy (gain 1.0019, corr 0.9999 at a 30 nm base;
`runs/redo96_oapuw`), so the READING is not the problem -- which is why (b) has to be
about the inversion rather than the measurement.

## 3. The run

One job, three starts in ONE invocation so the comparison is internal: 20 and 30 nm
(comfortably inside any plausible capture) plus **100 nm as the in-run control**, which
must reproduce the failure above.  A difference between two separate runs would be
worth nothing.

`loop.start_rms` is in **mm** (every gauge rms knob on this bench is), so 20 nm = 2e-5.
Each start gets its own calibration measured on its own starting surface.
`loop.unwrap` defaults to `'auto'`, which is ON whenever `start_rms` is set -- leave it.

```bash
cd ~/dev/MACOS_resources/mmacos/templates/40_benches/tg_psi_dm96_oap
export MACOS_HOME=$HOME/dev/macos/macos_f90
cp redo_oap_tail.mat descsmall_tail.mat          # the run MUST carry the redo bench's tail
OAP="'bench.optics','oap','bench.POL_IN','source','bench.D_RC_L2',125,'oap.OAP1_AOI',20,'oap.OAP2_AOI',25,'oap.OAP1_SIDE',1,'oap.OAP2_SIDE',-1"
nohup ./tg96_batch.sh descsmall "$OAP,'loop.start_rms',[2e-5 3e-5 1e-4],'loop.steps',[],'loop.drifts',{},'loop.floor',false,'loop.nph',1e13,'loop.recal_every',0,'stages',{'bench','loop','figs'}" </dev/null > runs/descsmall.nohup 2>&1 &
echo $! > runs/descsmall.pid
```

Do NOT set `TG96_NOWAIT` -- you want the wrapper's wait-for-another-MATLAB loop armed
on this box.  Expect ~1 hour (the Mac did two starts in 39 min at model 1024 / 385 px).

**Before you launch:** `ps -C MATLAB -o pid,etime,rss --no-headers` must come back
empty.  This box has 30 GB and one model-1024 job is ~11 GB; two have taken it down.

**The first thing to check when it finishes** -- the run is void without it:

```bash
grep 'Tail:' runs/descsmall/descsmall_report.txt
```

It must name `descsmall_tail.mat` and print `null 73.856 nm`.  A `Tail:` line naming
`oap_tail.mat` or "geometrically-scaled seed" means the copy above was skipped and the
run is on the RECORD's bench, not the redo bench -- discard it and start again.

Then:

```bash
grep -A6 'start     N/cyc' runs/descsmall/descsmall_report.txt
```

## 4. The verdict rule (decide this from the table, do not improvise)

- **rho ~ 0.5 and k(3pm) in the mid-teens at 20 and 30 nm, while the 100 nm row still
  fails** -> cause **(a)**.  The bench's gain is fine; the descent's capture range has
  shrunk because the bench's own 74 nm null eats it.  Say so and STOP -- the follow-on
  is a design question for Dave (start the descent inside capture, or recalibrate on
  the surface), not another run.
- **rho ~ 0.85 at every start, including 20 nm** -> cause **(b)**.  The start is
  irrelevant, so the inversion is over-damped for the new 3680-actuator lit set.  Then
  and only then, run leg 2 below.
- **Anything else** (e.g. rho recovers only at 20 nm, or the 100 nm control now
  converges) -> do not force it into either box.  Report the table as measured and say
  which prediction it breaks.

## 5. Leg 2 -- ONLY if the verdict is (b)

Same run, same tag plus `lam`, with the Tikhonov weight dropped two decades.  The
existing reg sweep on this bench (`runs/redo96_oapuw`, stage E) shows the dense-random
row at 0.9744 for lambda 1e-3 and 0.9869 for 1e-5, with the dark columns closing from
0.9519 to 0.9872 -- so 1e-5 is the end of the sweep that already exists, not a guess.

```bash
cp redo_oap_tail.mat descsmalllam_tail.mat
nohup ./tg96_batch.sh descsmalllam "$OAP,'battery.matrix_lam',1e-5,'loop.start_rms',[2e-5 3e-5 1e-4],'loop.steps',[],'loop.drifts',{},'loop.floor',false,'loop.nph',1e13,'loop.recal_every',0,'stages',{'bench','loop','figs'}" </dev/null > runs/descsmalllam.nohup 2>&1 &
echo $! > runs/descsmalllam.pid
```

If rho returns to ~0.5 with the weight alone, the regularization is the whole story and
the number to hand back is the lambda that does it.

## 6. Rules in force

- **Do not push.**  Commit locally if you like; Dave reviews before anything leaves.
- **Do not touch `~/dev/macos/demo_session/`** or `REPORT_bench_realism.md` -- CCL is
  editing the deck and the report in parallel and we will clobber each other.  Your
  output is the run directory plus a report back; CCL folds it.
- **One MATLAB at a time on this box.**  Kill by PID (`kill $(cat runs/descsmall.pid)`),
  never `pkill -f MATLAB` -- a name pattern self-matches and has killed another lane's
  run on this box before.
- The resources tree is shared with other lanes.  No `git commit --amend`, no rebase of
  anything you did not write, and expect files you did not touch to change under you.
- Bit-identical reruns are the norm here; if a repeat of the same tag gives different
  numbers, that is itself the finding.

## 7. Report back

Three things, in this order:

1. The `Tail:` line (the gate).
2. The ladder table, verbatim, all three starts.
3. One sentence: **(a)**, **(b)**, or "neither, and here is what it breaks" -- plus, if
   leg 2 ran, the lambda that restores rho.

Put the run under `runs/descsmall/` (and `runs/descsmalllam/`) where it lands by
default.  Nothing else needs writing.
