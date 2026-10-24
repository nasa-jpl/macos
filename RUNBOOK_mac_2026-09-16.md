# RUNBOOK -- the Mac as an execution box (Dave runs, 2026-09-16)

CCL for Dave.  CCMac's session is gone; its box (64 GB) can run TWO
model-1024 MATLABs at once (measured 11.9 GB each).  Everything below is
already scripted in the tree and pushed (macos a61f090, resources d28fdeb);
nothing needs editing on the Mac.  The wrappers detect macOS (no systemd cap)
and run a plain `matlab -batch`.  Each job writes `runs/<tag>.log` with an
exit line and `runs/<tag>/<tag>_report.txt`.

## 0. Preflight (5 min, once)

```bash
cd ~/dev/MACOS_resources && git pull
cd ~/dev/macos && git pull            # docs only; the mex on the Mac stays valid
export MACOS_HOME=$HOME/dev/macos/macos_f90
cd ~/dev/MACOS_resources/mmacos
matlab -batch "run('mmacos_setup.m'); macos.init(128); disp('mex ok'); exit(0)"
```
`mex ok` = go.  If the mex fails to load, stop and tell CCL (the Mac build is
CCMac's; MAC_PORT.md has the gfortran arms).

## Job A -- the interferometer's rows with the 10 mm plates (~40 min)

TO's queued `item4bseq.sh`: the lens rig's detector-leg retune with 10 mm
splitter and compensator and 4 mm lens edges, through the ADVISORY gate,
then the gate run (the single-actuator row on the 30 nm surface).  Closes
brief item 4's "does the 20 nm null cost the gauge?".

```bash
cd ~/dev/MACOS_resources/mmacos/templates/40_benches/tg_psi_dm96_oap
nohup runs/item4bseq.sh > runs/item4bseq.nohup 2>&1 &
```
Done when `grep '] exit' runs/thk22.log` shows `exit 0` (the retune's log is
`runs/thk22_tail.log`).  Numbers to send back: the retune's null line in
`runs/thk22_tail.log` (`null ... nm`; the tuner's winner was 20.09 nm), and
in `runs/thk22/thk22_report.txt` the `base/single` row (gain / err / floor
-- the lens rig of record reads 0.9916 / 2.2 pm).

## Job B -- the sensors' servo on the redesigned mirror rig (~2 h)

The deck's servo table has the interferometer on both front ends but the
three sensors only on the lens rig.  This is the oapsens22 bench (OAP1 20 /
OAP2 25 deg, fed at its focus, bare Al, the engine's arm maps) with the loop
stage: photons per cycle to hold 3 pm for S, V and P, noise-only / walk /
ramp.  Run it BESIDE job A (the bypass is for the Mac only):

```bash
cd ~/dev/MACOS_resources/mmacos/templates/40_benches/zwfs_dm96
ZWFS_NOWAIT=1 nohup ./zwfs_batch.sh oaploop22 "'bench.optics','oap','bench.OAP1_AOI',20,'bench.OAP2_AOI',25,'bench.OAP1_SIDE',1,'bench.OAP2_SIDE',-1,'bench.SRC_AT_FOCUS',true,'bench.MASK_TRIM',0,'bench.coat_oap','bareAl','MODEL',1024,'NGRID',193,'readings',{'S','V','P'},'loop.readings',{'S','V','P'},'stages',{'bench','loop','figs'},'battery.calib_surface','base','mask.v_arm','engine'" > runs/oaploop22.nohup 2>&1 &
```
Done when `grep '] exit' runs/oaploop22.log` shows `exit 0`.  Numbers to
send back: in `runs/oaploop22/oaploop22_report.txt` the block "photons per
cycle to hold 3.0 pm rms" (three readings x noise-only / walk / thermal) --
the lens rig of record: V 1.5e12 / 5.3e12, P 2.3e12 / 7.0e12, S 2.6e12 /
7.5e12.  The figure `oaploop22_loop.png` goes on the deck as is.

## Watching, and what "running" looks like

```bash
tail -3 runs/<tag>.log            # the wrapper's log; the report grows in runs/<tag>/
ps -ef | grep -c "MATLAB -batch"  # 2 while both run
```
A job that dies leaves no `] exit` line; send CCL the last 30 lines of the
log.  Do not start a third model-1024 job.

## Hand-back

When a job's exit line is 0:
```bash
cd ~/dev/MACOS_resources
git add mmacos/templates/40_benches/tg_psi_dm96_oap/runs/thk22 mmacos/templates/40_benches/tg_psi_dm96_oap/runs/thk22_tail.log mmacos/templates/40_benches/tg_psi_dm96_oap/thk22_tail.mat   # job A
git add mmacos/templates/40_benches/zwfs_dm96/runs/oaploop22                                                                            # job B
git commit -m "Mac runs: <tag> (exit 0)" && git push
```
Then on this box `git pull` and CCL folds the numbers into the deck.  The
`.mat` under runs/ are large but this tree tracks them for the record; the
`.nohup` files are not tracked (leave them).

## Not on the Mac (yet)

The staged tg96 patch, `item2bseq` and item 3 change source files TO's queue
here depends on -- they stay on this box.  The CTB regeneration at the
intended 42.75 mm beam is a day of one MATLAB and is the Mac's next job once
CCL has written its sequence (generator fix, deck regeneration, the study);
Dave says when.
