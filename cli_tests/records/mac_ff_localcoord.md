# (ff) the arm64 COORD frame — source-level bisection in LocalCoord (CCMac, 2026-10-20)

**Finding: (ff) is a SOURCE defect in `LocalCoord` (tracesub.F), not a gfortran codegen fault and
not the suspected argument aliasing. The y-axis rotate branch omits the `DSUB(yLocal,RayPos,crLocalPos,3)`
that the x-axis rotate branch has, so `yLocal` is stale (a ~1e-16 noise vector) when `DUNITVEC` normalizes
it; the optimizer merely decides whether that noise normalizes to `(0,-1,0)` (correct, by luck) or to a
garbage direction (the 76.5° frame). A one-line fix makes -O2 arm64 correct; it regresses nothing.**

Platform: macOS 26.6 (Darwin 25), Apple Silicon arm64, GNU Fortran (Homebrew GCC 16.1.0) 16.1.0.
Clone: `~/mac_gate_3ead7473` (engine `3ead7473`). Deck: `CassWithExitPupil.in`, `coord 6` (COORD_ACOOR.jou).
All scratch in `/tmp/mac_ff/`.

## 1. The measurement — the first divergent line

`LocalCoord` builds each stop-frame axis by a differential chief ray: translate the chief by `delMag=1e-6`
along xGrid/yGrid, trace, `DSUB(axis,RayPos,crLocalPos)`; if that differential is degenerate (`DMAG<1e-3`)
fall back to a `thMag=1e-6` *rotation*. Instrumented a copy of tracesub.F to print, in `(A,3ES24.16)` to
unit 77, every intermediate; built the SAME instrumented copy two ways and ran `coord 6`:

| print (stage) | full **-O0** (correct) | **-O2** / tracesub-only-O0 (rotated) |
|---|---|---|
| `C_xLocal_t` (x translate, DSUB) | `-5.33e-17  0  0` | `-5.33e-17  0  2.22e-16` |
| **`C_v1_x` = DUNITVEC(xLocal)** | **`-1  0  0`** | **`-0.23345  0  0.97237`** ← first divergence |
| `E_v2_x` (x rotate) | `-1  0  0` | `-1  0  1.98e-11` |
| `F_xLocal_final` | `-1  0  0` | `-1  0  1.98e-11` (x ends correct both) |
| `G_yLocal_t` (y translate, DSUB) | `0 -5.33e-17  0` | `0 -5.33e-17  2.22e-16` |
| `G_v1_y` = DUNITVEC(yLocal) | `0 -1  0` | `0 -0.23345  0.97237` |
| **`I_v2_y` = DUNITVEC(yLocal) in the y ROTATE branch** | **`0 -1  0`** | **`0 -0.23345  0.97237`** ← the one that sticks |
| `J_yLocal_final` | `0 -1  0` | `0 -0.23345  0.97237` |
| `K_zLocal` | `0 0 1` | `4.62e-12  0.97237  0.23345` (76.5° about x) |

The DSUB `axis = RayPos - crLocalPos` of the *translate* differential is numerical zero: at -O0 the z-components
cancel exactly (`...0`), at -O2 they leave a `2.22e-16` residual. `DUNITVEC` (mathsub.F:267) normalizes any
vector with `DMAG>0`, so that `|v|~2.3e-16` noise becomes a unit garbage direction. Because `DMAG1<1e-3` the
code correctly enters the rotate branch — but it has already stored garbage in `v1`.

## 2. Why x survives and y does not — the asymmetry

- **x rotate branch** (tracesub.F ~3071–3072): `CALL DSUB(xLocal,RayPos,crLocalPos,3)` **then**
  `CALL DUNITVEC(v2,xLocal,DMAG2,3)`. `v2` is recomputed from the freshly-traced rotated ray → the real
  x-axis (`-1,0,…`), `DMAG2~2.2e-5 > DMAG1`, so `v2` wins. x ends correct on every build.
- **y rotate branch** (tracesub.F ~3124): `CALL DUNITVEC(v2,yLocal,DMAG2,3)` with **NO preceding
  `DSUB(yLocal,RayPos,crLocalPos,3)`**. `yLocal` is still the *translate*-stage near-zero vector, so `v2`
  is `DUNITVEC` of the same ~1e-16 noise. At -O0 that noise is `(0,-5.33e-17,0)` → `(0,-1,0)` (correct by
  luck); at -O2 it is `(0,-5.33e-17,2.22e-16)` → `(0,-0.23345,0.97237)` (garbage). `DMAG2>DMAG1`, garbage wins.

That single missing recompute is the whole defect. It is latent on x86_64 / -O0 only because the stale
noise there has a zero z-component and normalizes back to the axis; arm64 -O1+ leaves the 2.22e-16 z-residual
and the normalization tilts 76.5°.

## 3. Experiments CC asked for

- **(a) whole tree -O1(here -O2) except tracesub.F at -O0:** STILL ROTATED, bit-identical to full -O2
  (frame and every print). So the opt level of LocalCoord *itself* is not the lever — proof it is not an
  -O-sensitive transformation inside LocalCoord but the arithmetic of a callee (the DSUB/DUNITVEC in
  mathsub.F, i.e. whether `RayPos-crLocalPos` cancels to 0 or 2.22e-16) feeding a stale read. Consistent with
  a source bug the optimizer only *exposes*.
- **(b) perturb a LOCAL instead of the module array (`tDir`), both rotate sites:** STILL ROTATED, bit-identical.
  **CC's aliasing hypothesis is ruled out** — `ChfRayDir` passed as an actual argument while use-associated is
  not the cause; the frame does not depend on it.
- **-O0 honesty leg:** the instrumented copy at full -O0 prints `(0,0,1)` — the instrumentation did not move
  the fault.

## 4. The fix (candidate patch, verified on arm64 -O2)

Add the missing recompute the x-branch already has, in the y rotate branch before `DUNITVEC(v2,yLocal,...)`:

```
--- a/macos_f90/tracesub.F
+++ b/macos_f90/tracesub.F
@@ -3121,6 +3121,7 @@ (y-axis rotate branch, after the undefined-ray guard)
             ChfRayPos(:) = oldChfRayPos(:)
             return
 	  END IF
+	  CALL DSUB(yLocal,RayPos,crLocalPos,3)   ! recompute from the rotated-ray RayPos (the x-branch has this; its absence left yLocal stale)
 	  CALL DUNITVEC(v2,yLocal,DMAG2,3)
 	  IF (DMAG2.GT.DMAG1) THEN
 	    CALL DEQUATE(yLocal,v2,3)
```

Built release `-O2` on arm64 with this one line:
- `coord 6` → xLocal `(-1,0,1.98e-11)`, yLocal `(3.9e-22,-1,1.98e-11)`, zLocal `(1.98e-11,1.98e-11,1)` —
  the correct frame to arm64 round-off (the 1.98e-11 is the same residual the correct xLocal carries, not the rotation).
- **COORD_ACOOR.jou PASSES** on the fixed arm64 binary (its EXPECT is the Linux frame).
- **cmd 112/112** (0 genuine fails; 2 must-fail red by design) — COORD cleared, nothing else moved.
- **load 164/164, 0 movers** vs `load_gfortran_3ead7473.csv`.

## 5. For TO

This is TO's to land in the engine (I did not change the shared tree; the clone's tracesub.F is restored
to pristine). Two considerations beyond the one-liner:
1. The one-line `DSUB(yLocal,...)` is the direct, minimal, symmetric fix and is sufficient (verified above).
   It makes the y rotate branch mirror the x rotate branch exactly.
2. A deeper hardening, if wanted: the real trigger is `DUNITVEC` normalizing a `|v|~1e-16` noise vector as
   if it were a direction. The `DMAG1<1e-3` test is meant to reject exactly that degenerate differential, but
   the garbage `v1`/`v2` is produced *before* the test is consulted. A magnitude guard in `DUNITVEC`
   (leave YHAT=0, or return a flagged DMAG, when `DMAG` is below a floor) would make the degeneracy explicit
   rather than relying on the `DMAG2>DMAG1` comparison to discard it. Not needed for (ff); noted because the
   x-branch is correct only because its real `v2` (2.2e-5) beats the noise `v1` (2.3e-16) — a thin margin.

Linux (CC) can mirror steps 1–3 on x86_64 at -O1 to show the prints agree there; the first divergent line to
look for is `C_v1_x` / the DSUB z-residual (`0` on x86_64, `2.22e-16` on arm64).

## Files
`/tmp/mac_ff/`: `tracesub_orig.F` (pristine), `tracesub_instr.F` (instrumented), `tracesub_fix.F` (the patch),
`tracesub_b.F` (experiment b); `ff_trace_{O0tr,debug,O2}.txt` (the print streams); `ff_fix.patch`;
`coord_*_*.log` (the frames); `load_fix.log`, `cmd_fix_full.log` (the regression legs).
