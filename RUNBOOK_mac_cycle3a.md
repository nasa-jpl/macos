# RUNBOOK -- Mac cycle 3a: the substrates in, the full 96 mm beam, the seed tail -- both rigs (2026-09-17)

CCL for Dave, to execute on the Mac by hand.  Four jobs, two at a time (64 GB).
What it answers before TO touches the tuner: do the rows and the null hold on the
bench we will build -- the record's optics, the decided substrates (10 mm splitter and
compensator, 2 mm fused-silica plates under every polarizing element, the 2 mm mask
plate, 4 mm lens edges, protected aluminum on the parabolas), the DM as the stop (the
full 96 mm beam), and the SEED detector-leg geometry (the field lens 10.8 mm past the
focus, which images the DM flat) -- on both rigs, with the pupil stage on each emitted
deck.  The lens rig's collimator is still the record's (the redo fixes it in package A);
these are the baseline the redo's rows are gated against.  All of this is in the sheet
defaults now (resources commit of 2026-09-17, "substrates baked in"), so the jobs pass
only the rig and the stages.

## Setup (pull first: the sheets changed)

```bash
cd ~/dev/MACOS_resources && git pull
export MACOS_HOME=$HOME/dev/macos/macos_f90
cd ~/dev/MACOS_resources/mmacos/templates/40_benches/tg_psi_dm96_oap
OAP="'bench.optics','oap','bench.POL_IN','source','bench.SRC_AT_FOCUS',true,'bench.D_RC_L2',125,'oap.OAP1_AOI',20,'oap.OAP2_AOI',25,'oap.OAP1_SIDE',1,'oap.OAP2_SIDE',-1,'clear.BODY',struct('Baffle',50,'Detector',50,'TestOptic',90,'PZT',60)"
```
(`OAP` is TO's redesigned mirror rig, exactly as `oapifo2` ran it; the coating now comes
from the sheet: protected aluminum.)

## Job A -- interferometer, lens rig: rows on the 30 nm surface, seed tail (~2 h, 12 GB)

```bash
TG96_NOWAIT=1 nohup ./tg96_batch.sh sub96_lens "'bench.tail_from_mat',false,'stages',{'bench','deck','figs'}" > runs/sub96_lens.nohup 2>&1 &
```
Watch: `tail -3 runs/sub96_lens.log`; done when it prints `[tg96_batch] exit 0`.
Report head to send: `grep -A12 'deck: rows' runs/sub96_lens/sub96_lens_report.txt | head -20`
and the flat null line: `grep -i 'null' runs/sub96_lens/sub96_lens_report.txt | head -3`.

## Job B -- interferometer, mirror rig: the same (~2 h, 12 GB; start beside A)

```bash
TG96_NOWAIT=1 nohup ./tg96_batch.sh sub96_oap "$OAP,'bench.tail_from_mat',false,'stages',{'bench','deck','figs'}" > runs/sub96_oap.nohup 2>&1 &
```
Same watch and report lines with `sub96_oap`.

## Job C -- the pupil stage on both emitted decks (~15 min, 3 GB; after A and B)

```bash
TG96_NOWAIT=1 ./tg96_pupil_batch.sh lens "'tool','sim','deck','$PWD/runs/sub96_lens/sub96_lens_test.in','tag','pupilsim_sub96_lens'"
TG96_NOWAIT=1 ./tg96_pupil_batch.sh oap  "'tool','sim','deck','$PWD/runs/sub96_oap/sub96_oap_test.in','tag','pupilsim_sub96_oap'"
```
Report lines to send: `grep -h 'BEST PUPIL\|working surface 3\|pupil distortion' runs/pupilsim_sub96_*/pupilsim_sub96_*_report.txt`.
Expected: the seed tail's numbers of REPORT_bench_realism 7.1 (Nyquist gain >= 0.998
worst, distortion < 0.01 mm rms) on both rigs, now with the substrates in.

## Job D -- the sensors on the mirror rig, substrates in, 385 px (~3.5 h, 20 GB; after A or B frees memory)

```bash
cd ~/dev/MACOS_resources/mmacos/templates/40_benches/zwfs_dm96
OAPZ="'bench.optics','oap','bench.OAP1_AOI',20,'bench.OAP2_AOI',25,'bench.OAP1_SIDE',1,'bench.OAP2_SIDE',-1,'bench.SRC_AT_FOCUS',true,'bench.MASK_TRIM',0,'mask.v_arm','engine','battery.calib_surface','base'"
ZWFS_NOWAIT=1 nohup ./zwfs_batch.sh sub96_oapsens "$OAPZ,'MODEL',2048,'NGRID',385,'param_file','macos_param_2048.txt','readings',{'S','V','P'},'stages',{'bench','battery'},'battery.rows',{'base/single','base/grid','base/rand'}" > runs/sub96_oapsens.nohup 2>&1 &
```
(the cycle-2 `oapsens385` job with the sheet's new defaults: the mask plate is in the
converging beam now, the plates under the polarizing elements, the coating protected
aluminum.)  Report lines: `grep -h 'capture range\|rand30/single10\|rand30/grid\|camera: pupil' runs/sub96_oapsens/sub96_oapsens_report.txt`.

## Harvest

```bash
cd ~/dev/MACOS_resources
git add mmacos/templates/40_benches/tg_psi_dm96_oap/runs/sub96_lens mmacos/templates/40_benches/tg_psi_dm96_oap/runs/sub96_oap mmacos/templates/40_benches/tg_psi_dm96_oap/runs/pupilsim_sub96_lens mmacos/templates/40_benches/tg_psi_dm96_oap/runs/pupilsim_sub96_oap mmacos/templates/40_benches/zwfs_dm96/runs/sub96_oapsens
git commit -m "Mac cycle 3a: <tags> (exit codes)" && git pull --rebase && git push
```
Send the four report heads; CCL folds them (deck: the rows with the substrates on
both rigs, the seed tail's null, the sensors' rows at 385 with the mask plate) and
they become the baseline the redo's package C is gated against.

## What a failure looks like, and what to do

- `exit 127` at launch: the wrapper's `flock` guard is in (resources 174ef6b); if it
  recurs, send the first 5 lines of the `.log`.
- A job that dies in stage B with a clearance assertion: the 10 mm plates moved a node
  part inside the 25 mm margin; send the `Stage A` lines of the report and stop that
  rig's jobs (the other rig's are independent).
- Memory: `ps aux | grep MATLAB` -- A and B together are ~24 GB, D ~20 GB; never three
  model-1024-or-more jobs at once on the 64 GB box.
