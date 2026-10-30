# RUNBOOK -- Mac cycle 2 (Dave runs, 2026-09-16 afternoon)

Three scripted runs that fill the deck's remaining mirror-rig gaps.  All on
the redesigned reflective bench (the oapsens22 settings: OAP1 20 / OAP2 25
deg, fed at its focus, bare Al, the engine's arm maps).  The wrappers now
skip `flock` on macOS (resources 3c88b46, pushed) -- pull first.  Jobs A and
B run side by side (12 GB each); job C alone or beside one of them (20 GB).

## 0. Preflight (2 min)

```bash
cd ~/dev/MACOS_resources && git pull
export MACOS_HOME=$HOME/dev/macos/macos_f90
cd ~/dev/MACOS_resources/mmacos/templates/40_benches/zwfs_dm96
OAP="'bench.optics','oap','bench.OAP1_AOI',20,'bench.OAP2_AOI',25,'bench.OAP1_SIDE',1,'bench.OAP2_SIDE',-1,'bench.SRC_AT_FOCUS',true,'bench.MASK_TRIM',0,'bench.coat_oap','bareAl','mask.v_arm','engine','battery.calib_surface','base'"
```
(`OAP` is a shell variable holding the bench settings; each job below pastes
it in.  Keep the same terminal, or re-run these three lines in a new one.)

## Job A -- photons for 1 pm on the mirror rig (~40 min, 12 GB)

The deck's "Photons for 1 pm" slide has the sensors on the lens rig only.
Noise stage: the readings priced over 1e8 to 1e14 photons per measurement on
the 30 nm surface.

```bash
ZWFS_NOWAIT=1 nohup ./zwfs_batch.sh oapnoise22 "$OAP,'MODEL',1024,'NGRID',193,'readings',{'S','V','P'},'noise.readings',{'S','V','P'},'stages',{'bench','noise'},'noise.nstates',10.^(8:2:14),'noise.nreal',6,'battery.base_rms',30e-6" > runs/oapnoise22.nohup 2>&1 &
```
Send back: `runs/oapnoise22/oapnoise22_report.txt`, the "photons for 1 pm"
lines per reading (lens rig of record: pinhole 9.8e13, vector 6.1e13,
stepped 8.8e13 on the surface).

## Job B -- capturing the initial figure with the sensors on the mirror rig (~2.5 h, 12 GB)

The capture slide has the sensors' descents on the lens rig only.  Starts of
60 and 100 nm rms of surface, unwrapping on, the matrix re-measured every 10
cycles, 1e13 and 1e15 photons per cycle.

```bash
ZWFS_NOWAIT=1 nohup ./zwfs_batch.sh oapcap22 "$OAP,'MODEL',1024,'NGRID',193,'readings',{'S','V','P'},'loop.readings',{'S','V','P'},'stages',{'bench','loop','figs'},'loop.start_rms',[60 100]*1e-6,'loop.nph',[1e13 1e15],'loop.drifts',{},'loop.floor',false,'loop.steps',[],'loop.K',60,'loop.recal_list',[0 10],'loop.unwrap',true" > runs/oapcap22.nohup 2>&1 &
```
(`loop.start_rms` is in mm: `[60 100]*1e-6` is 60 and 100 nm.  A bare 60
would ask for a 60 mm surface and the loop now warns.)
Send back: the "DESCENT LADDER" table in `runs/oapcap22/oapcap22_report.txt`
(lens rig of record: vector and pinhole capture from 60 nm; the pinhole with
its shutter frame from 100 with unwrap + recal).

## Job C -- the mirror rig at the compliant sampling (~3 h, 20 GB; start after A finishes)

The deck's footnote "most rows at 193 pixels; the 385 checks read the same"
is measured on the lens rig only.  Model 2048, 385 rays, the battery rows.

```bash
ZWFS_NOWAIT=1 nohup ./zwfs_batch.sh oapsens385 "$OAP,'MODEL',2048,'NGRID',385,'param_file','macos_param_2048.txt','readings',{'S','V','P'},'stages',{'bench','battery'},'battery.rows',{'base/single','base/grid','base/rand'}" > runs/oapsens385.nohup 2>&1 &
```
A model-2048 MATLAB may exit with code 137 AFTER writing everything (a known
exit-time segfault); the report and .mat are complete if the report's last
line is "run complete".  Send back: the three rows (single / grid / dense) for
S, V, P against oapsens22's 193-ray rows.

## Watching

```bash
tail -2 runs/<tag>.log
ps -ef | grep -c "MATLAB -batch"     # 2 while A and B run
```

## Hand-back (after each exit line)

```bash
cd ~/dev/MACOS_resources
git add mmacos/templates/40_benches/zwfs_dm96/runs/oapnoise22 mmacos/templates/40_benches/zwfs_dm96/runs/oapcap22 mmacos/templates/40_benches/zwfs_dm96/runs/oapsens385
git commit -m "Mac cycle 2: <tags> (exit codes)" && git pull --rebase && git push
```
Add only the directories that exist yet; a second commit later is fine.
