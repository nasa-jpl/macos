# RUNBOOK -- Mac cycle 3c: the rest of package C on the redo bench (2026-09-17)

CCL for Dave, to execute on the Mac by hand.  Cycle 3b gave the rows on the redo bench
(lens grid row 0.993, null 59 nm; mirrors 0.981 / 26.5 nm; the sensors at 385
unchanged).  This cycle runs everything else the deck quotes for the interferometer,
on both rigs, and the sensors' photons and servo on the mirror rig -- each the record's
own invocation on the redo bench.  Two lanes, one terminal each, two MATLABs at a time
on the 64 GB box: lane M is the mirror rig, lane L the lens rig; each lane runs its
jobs in sequence.  About 12 hours per lane; start both, come back to the report lines.

## Setup (once, in each terminal)

```bash
cd ~/dev/MACOS_resources && git pull
export MACOS_HOME=$HOME/dev/macos/macos_f90
cd ~/dev/MACOS_resources/mmacos/templates/40_benches/tg_psi_dm96_oap
grep -c 'SRC_AT_FOCUS = true' tg96_params.m       # must print 1
ls redo_lens_tail.mat redo_oap_tail.mat           # both must exist
# the runner finds a tail by its TAG and falls back to the RECORD's tuned tail, so every tag of this
# cycle gets its copy of the redo tail:
for t in redo96_oaploop redo96_oapdesc redo96_oapuw redo96_oapwrap redo96_oapstn; do cp redo_oap_tail.mat ${t}_tail.mat; done
for t in redo96_lensloop redo96_lensdesc redo96_lensuw redo96_lenswrap redo96_lensstn; do cp redo_lens_tail.mat ${t}_tail.mat; done
OAP="'bench.optics','oap','bench.POL_IN','source','bench.D_RC_L2',125,'oap.OAP1_AOI',20,'oap.OAP2_AOI',25,'oap.OAP1_SIDE',1,'oap.OAP2_SIDE',-1"
```
Every `Tail:` line in every report must name that run's own `redo96_*_tail.mat`; a run
whose `Tail:` line says `oap_tail.mat`, `lens_tail.mat` or "geometrically-scaled seed"
is on the wrong tail and is discarded.

## Lane M -- the mirror rig (terminal 1; ~12 h in sequence)

```bash
# M1 servo (the record: oapifol2, ~5.5 h): 3 pm from N photons per cycle, noise-only and under the 2 pm walk; the thermal floor
TG96_NOWAIT=1 ./tg96_batch.sh redo96_oaploop "$OAP,'stages',{'bench','loop','figs'}" ; \
# M2 descent (oapdesc2, ~2 h): from 100 and 200 nm of surface, unwrap on, matrix measured at the start
TG96_NOWAIT=1 ./tg96_batch.sh redo96_oapdesc "$OAP,'loop.start_rms',[1e-4 2e-4],'loop.steps',[],'loop.drifts',{},'loop.floor',false,'loop.nph',1e13,'loop.recal_every',0,'stages',{'bench','loop','figs'}" ; \
# M3 unwrapped capture ladder (oapuw2, ~1 h): the single-site row on a growing base with unwrapping
TG96_NOWAIT=1 ./tg96_batch.sh redo96_oapuw "$OAP,'battery.unwrap',true,'battery.rows',{'base/single'},'stages',{'bench','battery'}" ; \
# M4 the wrap stage (wrapoap, ~1 h): where the raw four-step folds on this bench (the record: 60-120 nm on both rigs)
TG96_NOWAIT=1 ./tg96_batch.sh redo96_oapwrap "$OAP,'stages',{'bench','wrap'}" ; \
# M5 the station figure (stnoap, ~1 h): the signals along the train, flat and on the 30 nm surface, 1800 px
TG96_NOWAIT=1 ./tg96_batch.sh redo96_oapstn "$OAP,'stages',{'bench','figs'}"
```
Then the sensors on the mirror rig (the cycle-2 jobs on the redo bench; the seat 0.632
explicit, as cycle 3b):
```bash
cd ~/dev/MACOS_resources/mmacos/templates/40_benches/zwfs_dm96
OAPZ="'bench.optics','oap','bench.OAP1_AOI',20,'bench.OAP2_AOI',25,'bench.OAP1_SIDE',1,'bench.OAP2_SIDE',-1,'bench.MASK_TRIM',0.632,'mask.v_arm','engine','battery.calib_surface','base'"
# M6 photons for 1 pm (oapnoise22, ~40 min)
ZWFS_NOWAIT=1 ./zwfs_batch.sh redo96_oapnoise "$OAPZ,'MODEL',1024,'NGRID',193,'readings',{'S','V','P'},'noise.readings',{'S','V','P'},'stages',{'bench','noise'},'noise.nstates',10.^(8:2:14),'noise.nreal',6,'battery.base_rms',30e-6" ; \
# M7 the sensors' servo + descent (oapcap22, ~3 h)
ZWFS_NOWAIT=1 ./zwfs_batch.sh redo96_oapcap "$OAPZ,'MODEL',1024,'NGRID',193,'readings',{'S','V','P'},'loop.readings',{'S','V','P'},'stages',{'bench','loop','figs'},'loop.start_rms',[60 100]*1e-6,'loop.nph',[1e13 1e15],'loop.drifts',{},'loop.floor',false,'loop.steps',[],'loop.K',60,'loop.recal_list',[0 10],'loop.unwrap',true"
```

## Lane L -- the lens rig (terminal 2; ~11 h in sequence)

```bash
cd ~/dev/MACOS_resources/mmacos/templates/40_benches/tg_psi_dm96_oap
TG96_NOWAIT=1 ./tg96_batch.sh redo96_lensloop "'stages',{'bench','loop','figs'}" ; \
TG96_NOWAIT=1 ./tg96_batch.sh redo96_lensdesc "'loop.start_rms',[1e-4 2e-4],'loop.steps',[],'loop.drifts',{},'loop.floor',false,'loop.nph',1e13,'loop.recal_every',0,'stages',{'bench','loop','figs'}" ; \
TG96_NOWAIT=1 ./tg96_batch.sh redo96_lensuw "'battery.unwrap',true,'battery.rows',{'base/single'},'stages',{'bench','battery'}" ; \
TG96_NOWAIT=1 ./tg96_batch.sh redo96_lenswrap "'stages',{'bench','wrap'}" ; \
TG96_NOWAIT=1 ./tg96_batch.sh redo96_lensstn "'stages',{'bench','figs'}"
```
(The lens rig's record for these: loop_lens, descent_lens, lensuw2, wraplens, stnlens --
the last with its open 62 nm misregistration, which the redo bench should close: the
pupil-image bowl was the size it predicted.)

## Report lines, per job (send as they land)

```bash
cd ~/dev/MACOS_resources/mmacos/templates/40_benches/tg_psi_dm96_oap
grep -h 'Tail:' runs/redo96_*/redo96_*_report.txt
grep -h 'photons per cycle\|thermal\|floor' runs/redo96_oaploop/redo96_oaploop_report.txt runs/redo96_lensloop/redo96_lensloop_report.txt | head -12
grep -h 'DESCENT\|k(3pm)\|reaches\|rho' runs/redo96_oapdesc/redo96_oapdesc_report.txt runs/redo96_lensdesc/redo96_lensdesc_report.txt | head -16
grep -h 'unwrap\|capture\|beyond-fold' runs/redo96_oapuw/redo96_oapuw_report.txt runs/redo96_lensuw/redo96_lensuw_report.txt | head -12
grep -h 'wrap\|fold' runs/redo96_oapwrap/redo96_oapwrap_report.txt runs/redo96_lenswrap/redo96_lenswrap_report.txt | head -10
grep -h 'station\|residual vs the engine' runs/redo96_oapstn/redo96_oapstn_report.txt runs/redo96_lensstn/redo96_lensstn_report.txt | head -6
cd ../zwfs_dm96
grep -h 'photons for 1 pm\|1 pm' runs/redo96_oapnoise/redo96_oapnoise_report.txt | head -8
grep -h 'photons per cycle to hold\|CONVERGES' -A4 runs/redo96_oapcap/redo96_oapcap_report.txt | head -16
```
What the deck expects (the record on the old bench, to be beaten or matched): servo
3 pm from 5.4e12 photons per cycle noise-only and 1.7e13 under the 2 pm walk (mirrors),
5.5e12 / 2.0e13 (lens), thermal floor ~10 pm; descent to 3 pm from 200 nm in ~19 cycles;
the raw four-step wraps between 60 and 120 nm on both rigs; the mirror rig's station
figure 626 pm on the 30 nm surface; the sensors 3 pm from 2-8e12 under the walk and
photons for 1 pm 6-9e13.

## Harvest (each lane, when done)

```bash
cd ~/dev/MACOS_resources
git add mmacos/templates/40_benches/tg_psi_dm96_oap/runs/redo96_* mmacos/templates/40_benches/zwfs_dm96/runs/redo96_* mmacos/templates/40_benches/tg_psi_dm96_oap/redo96_*_tail.mat
git commit -m "Mac cycle 3c: <tags> (exit codes)" && git pull --rebase && git push
```

## Failure signatures

- A `Tail:` line naming the record's tail: the copy loop was skipped for that tag.
- `exit 1` inside the loop stage after hours: send the last 30 lines of that job's
  `.log`; the other lane is independent.
- Memory: two model-1024 jobs are ~24 GB; do not start a third.
