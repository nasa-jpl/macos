# BRIEF for CCMac: JPL-private verification, round 2 -- rulings on IRIS + OPTIIX, three measurements, then the merge

From CCL for Dave, 2026-09-09.  Follows `BRIEF_ccmac_jpl_private_verification.md`
(round 1; amended today at `3cc4549` with the post-brief engine commit --
read its section 1 addendum first).  Cold-start: everything needed is here.

## 0. Where things stand

Both JPL decks are inside the expectation tables.  The merge gate -- the
supervisor axis, PRE vs POST on the same NEW engine -- passed on both:
every `reset_xp=false` Jacobian byte-identical, every other difference
error-parity or the designed emptyOPD softening.  The engine axis is
classified: IRIS FEX moved -0.35 mm (half the T/S split, expected), OPTIIX
FEX bit-identical on axis.  Dave's call, executed 2026-09-09: `sens-core`
MERGED into resources `dev-candidate` at `57a6ec0` and the branch deleted
(local + origin).  Run everything below on `dev-candidate` >= `57a6ec0`
(resources) with the macos `dev-candidate` engine >= `cdf8636`; the
PRE/POST pair no longer exists.  Three items stay open as measurements,
none of them a sens-core regression; they are below, cheapest first.

## 1. OPTIIX Test 2c -- RULED: expected (the finite-difference floor)

Your candidate: `reset_xp=true` Jacobians differ OLD<->NEW at 1.25e-7
relative (6.8e-5 absolute), uniform across fields including on-axis C
whose FEX is bit-identical; nominal w0 differs by 9.3e-13.

Ruling: expected.  Not the four-probe arithmetic (a bit-identical FEX
OUTPUT cannot propagate through the geometry it writes) but engine
STATE around the FEX call: `81d3308` (idempotent re-traces) landed after
round 1 was written and is in your NEW engine, not in OLD `82ced2b`.
The numbers are the floor's: with the supervisors' default step
delta = 1e-8, 6.8e-5 absolute = 1.4e-12 mm between the +/- traces, a few
ulp of the OPTIIX path; w0's 9.3e-13 is the same class.  The NEW Jacobian
is the clean one.  Three discriminators if you want it nailed (one MATLAB
each): `git rev-parse HEAD` of NEW contains `81d3308`; twice-trace
idempotency at field C (`trace`, `modify`, `trace`; max|w1-w2| = 0 on
NEW, a few ulp on OLD); the OLD-NEW difference falls as 1/delta at
delta 1e-7 and 1e-6.  Optional -- the ruling stands without them.

The stop-less OPTIIX run adds nothing (no stop -> `reset_xp=true` is
noStop error-parity by design; the rest is already identical).  Skip it.

## 2. IRIS `reset_xp=true` empty on axis under the stop -- NOT benign by construction; two numbers decide

The guard fires on `nnz(w_nom_2d) == 0`, and the no-ray sentinel is
also 0, so it cannot tell "0 rays survived at elt 55" from "every OPD
sample exactly 0".  Your reading needs the second; a real design does
not give 0.0 residual, so the first is the working hypothesis: the
FEX-placed sphere loses the beam on this deck.  On-axis is outside round
1's row 4 (that row is off-axis vignetting from too-wide fields).  PRE
errored identically, so it is not a sens-core regression -- but the
default `reset_xp=true` yields nothing on IRIS with a stop, and that
should be understood.  One MATLAB, NEW engine, POST tree, same model
size you used:

```matlab
nE = macos.load_rx(rx);  macos.stop(24);
s0 = macos.trace(nE-1);  r0 = macos.get_ray_status(s0.nRays);  w0 = macos.opd(nE-1);
i0 = macos.get_elt_info(nE-1);  k0 = macos.get_elt_kr(nE-1);   % vertex/psi/kr before
macos.fex(1);                                                   % exactly what reset_xp does
s1 = macos.trace(nE-1);  r1 = macos.get_ray_status(s1.nRays);  w1 = macos.opd(nE-1);
i1 = macos.get_elt_info(nE-1);  k1 = macos.get_elt_kr(nE-1);
fprintf('rays %d->%d | pass %d->%d | lost %d->%d | nnz(w) %d->%d | max|w| %.3g->%.3g | kr %.6g->%.6g | vpt moved %.3g\n', ...
  s0.nRays, s1.nRays, nnz(r0.status==0), nnz(r1.status==0), nnz(r0.status>=2), nnz(r1.status>=2), ...
  nnz(w0), nnz(w1), max(abs(w0(:))), max(abs(w1(:))), k0, k1, norm(i1.vpt(:)-i0.vpt(:)));
```
Read: rays lost after `fex` -> STOP-AND-REPORT (engine: the placed
sphere misses the beam on this deck; send those numbers).  Rays kept and
`nnz(w1)==0` with `max|w1|==0` -> your mechanism holds and the guard has
a blind spot (sens-core follow-up: key the empty test on the ray-pass
mask, not nnz) -- benign for the merge.

## 3. IRIS `save_rx` -> reload SIGSEGV -- not reproducible here; two cheap asks before the OLD-engine check

Here, 7 decks x 2 engines (gfortran, ifx) x {no stop, element stop},
load -> save -> reload -> OPD: OPD rms and ray loss identical pre/post in
all 28 runs, saves byte-identical across compilers.  Decks: five
FreeForm/grid (`tst_FF_fg/mg/g`, `FFSegDemoData`, `SegDemo3data`) and
two synthetic non-sequential ones built from `Rx_Cass_NS` (a grid on a
Conic NSReflector; a FreeForm+Zernike+grid NSReflector).  The stop is
not the operative ingredient in anything buildable here; the IRIS
factor is missing.  Numbers only, no deck content:

a. The same round trip WITHOUT a stop on iris_dp_ZGD.  If it also
   crashes, the stop is a red herring and the SAVE writer is the suspect
   (last changed July, `662e86e`) -- pre-existing, does not gate the merge.
b. A diff of original vs saved restricted to the six FreeForm elements'
   grid keys -- run in the deck directory:
   `diff <(grep -E 'iElt=|nGridMat|GridFile|GridSrfdx|lFF=|pData|xData|yData|zData|FFZern' <orig>.in) <(grep -E '...same...' <saved>.in)`
   First suspect: a `GridSrfdx=` or frame line the writer emits as ZERO
   where the original never set it -- a zero dx is a division by zero in
   FreeFormSrf's grid index, which is the SIGSEGV signature you saw.
c. The model size used, and whether the grids are 256 (mGridMat caps at
   256 regardless of model size).
d. THEN the OLD-engine check, only if (a) also crashes.

## 4. After the merge

The byte-identity gate is proven on the two most different decks in the
JPL set, so the remaining IRIS decks need only the intake card, the
Test 3 pupil-read check, and a note on whether `reset_xp=true` produces
rows (item 2's numbers on any deck where it does not).  No PRE/POST pair
once sens-core is on dev-candidate; the OLD/NEW engine axis is optional.

## 5. Two engine items found on this side, for the record (not yours)

- `Get_Values` (iosub.inc:3432) reads 36 bytes past its 220-char
  buffer: any deck with `ArrWaveLen=` / `ArrIndRef=` parses
  nondeterministically (the July fixture `tst_save_keys.in` loaded 3 of
  5 times here).  FIXED at macos `cdf8636` (gate 20/20 loads on both
  compilers) -- rebuild your NEW engine from that tip or later.  If an
  IRIS deck ever failed to load with "Bad real number" on the older
  engine, that was this, not you.
- `macos.trace(26)` then `trace(27)` on the jwst zoom deck returns a
  wrong first OPD (PLAN section 0); harvests never see it.
