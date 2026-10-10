# RUNBOOK -- Mac cycle 3b: package C's first cut on the redo bench (2026-09-17)

CCL for Dave, to execute on the Mac by hand.  Package A passed both rigs (TO, 13:30:
the bench collimated for real, the tails tuned on the reading, re-emitted under the
tags `redo_lens` and `redo_oap`; REPORT_bench_realism section 8).  This cycle runs the
first of package C's rows on that bench: the interferometer's DECK stage on both rigs
(rows on the 30 nm surface with the matrix on it, the null, the capture ladder's aging
rung) and the sensors at 385 on the mirror rig.  Cycle 3a's `sub96_*` runs (the seed
tail, the record's collimator) are the baseline these are read against.

**Precondition: TO's package-A commits must be on origin** (they were local to the Linux
box at 13:30: `d1a3b57` and twelve before it).  The Mac's pull must bring
`tg96_params.m` with `SRC_AT_FOCUS true`, `L1_Kr 249.246`, `MASK_TRIM 1.2318`, and the
tail files `redo_lens_tail.mat`, `redo_oap_tail.mat` in `tg_psi_dm96_oap/`.  Check
before launching:

```bash
cd ~/dev/MACOS_resources && git pull
cd mmacos/templates/40_benches/tg_psi_dm96_oap
grep -n 'SRC_AT_FOCUS = true\|L1_Kr = 249' tg96_params.m | head -2      # both lines must print
ls redo_lens_tail.mat redo_oap_tail.mat                                  # both must exist; if not, STOP and say so
```

## Setup

```bash
export MACOS_HOME=$HOME/dev/macos/macos_f90
cd ~/dev/MACOS_resources/mmacos/templates/40_benches/tg_psi_dm96_oap
# the runner finds a tail by its TAG (<tag>_tail.mat) and falls back to <optics>_tail.mat = the RECORD's tuned tail,
# so the redo tails are copied under this cycle's tags:
cp redo_lens_tail.mat redo96_lens_tail.mat;  cp redo_oap_tail.mat redo96_oap_tail.mat
OAP="'bench.optics','oap','bench.POL_IN','source','bench.D_RC_L2',125,'oap.OAP1_AOI',20,'oap.OAP2_AOI',25,'oap.OAP1_SIDE',1,'oap.OAP2_SIDE',-1"
```
(`OAP` is exactly TO's `redo_oap` emission set; the collimation, seat and coating now come
from the sheet.)

## Job A -- interferometer, lens rig, redo bench: rows, null, aging rung (~2 h, 12 GB)

```bash
TG96_NOWAIT=1 nohup ./tg96_batch.sh redo96_lens "'stages',{'bench','deck','figs'}" > runs/redo96_lens.nohup 2>&1 &
```
Watch `tail -3 runs/redo96_lens.log`; done at `[tg96_batch] exit 0`.  Report lines:
```bash
grep -n 'Tail:\|flat-DM null\|single10\|grid@1nm\|capture range to\|camera: pupil' runs/redo96_lens/redo96_lens_report.txt
```
The `Tail:` line must name `redo96_lens_tail.mat` (not "geometrically-scaled seed", not
`lens_tail.mat`).  Gates against cycle 3a (`sub96_lens`: null 92 nm, single 1.003, grid
0.754): the null must fall well under a quarter wave of surface, the grid row at 1 nm
>= 0.98, the single row 0.99-1.01.

## Job B -- interferometer, mirror rig, redo bench (~2 h, 12 GB; beside A)

```bash
TG96_NOWAIT=1 nohup ./tg96_batch.sh redo96_oap "$OAP,'stages',{'bench','deck','figs'}" > runs/redo96_oap.nohup 2>&1 &
```
Same lines with `redo96_oap`; the `Tail:` line must name `redo96_oap_tail.mat`.  Cycle
3a's `sub96_oap`: null 26 nm, rows 0.986 / 0.981.

## Job C -- the sensors on the mirror rig at 385, redo bench (~3.5 h, 20 GB; after A or B)

```bash
cd ~/dev/MACOS_resources/mmacos/templates/40_benches/zwfs_dm96
OAPZ="'bench.optics','oap','bench.OAP1_AOI',20,'bench.OAP2_AOI',25,'bench.OAP1_SIDE',1,'bench.OAP2_SIDE',-1,'bench.MASK_TRIM',0.632,'mask.v_arm','engine','battery.calib_surface','base'"
ZWFS_NOWAIT=1 nohup ./zwfs_batch.sh redo96_oapsens "$OAPZ,'MODEL',2048,'NGRID',385,'param_file','macos_param_2048.txt','readings',{'S','V','P'},'stages',{'bench','battery'},'battery.rows',{'base/single','base/grid','base/rand'}" > runs/redo96_oapsens.nohup 2>&1 &
```
No `SRC_AT_FOCUS` on the line (the sheet carries the collimation).  **`MASK_TRIM` 0.632
explicitly:** the sheet's per-run seat scan (`'scan'`) converged on the Mac's first
launch to -5.93 mm with a mask-plane peak of 2e-6 (the lens record's neighborhood,
not a focus) and the run died; +0.632 mm is the mirror rig's seat with the plate in,
the value cycle 3a passed at.  The scan's start and range on the mirror rig are TO's
to fix (brief section 8).  Report lines:
```bash
grep -n 'focal spot\|camera: pupil\|single10\|grid@1nm\|capture range to\|G5' runs/redo96_oapsens/redo96_oapsens_report.txt
```
Cycle 3a (`sub96_oapsens`): rows 0.993 / 0.995, capture S 35 / V 51 / P 42 nm.

## Harvest

```bash
cd ~/dev/MACOS_resources
git add mmacos/templates/40_benches/tg_psi_dm96_oap/runs/redo96_lens mmacos/templates/40_benches/tg_psi_dm96_oap/runs/redo96_oap mmacos/templates/40_benches/zwfs_dm96/runs/redo96_oapsens mmacos/templates/40_benches/tg_psi_dm96_oap/redo96_lens_tail.mat mmacos/templates/40_benches/tg_psi_dm96_oap/redo96_oap_tail.mat
git commit -m "Mac cycle 3b: redo96_lens redo96_oap redo96_oapsens (exit codes)" && git pull --rebase && git push
```

## Next cycle (3c), once these pass: the rest of package C on the Mac

The capture ladders and their photon price (`cap`/`noise` stages), the servo and the
descent (`loop` stage: `loop.start_rms [60 100]*1e-6`, unwrap on), the station figures,
D1/D4, the plates and the polarization at the built angle -- both rigs, the mirror rig
first.  CCL writes it from cycle 3b's reports.

## Failure signatures

- `exit 1` within a minute: the sheet or a tail file did not arrive (see the
  precondition); send `tail -20 runs/<tag>.log`.
- `mask plane not focused` in job C: the seat scan did not run; send the bench lines
  of the report.
- A `Tail:` line naming `lens_tail.mat` / `oap_tail.mat`: the copy step was skipped; the
  run is on the record's tuned tail and must be discarded.
