# dev_optimization_surfsub — archived patches

Five commits by Norbert Sigrist (TraceMoreRays, `sigrist1000@gmail.com`)
that lived on the `dev_optimization_surfsub` branch of the macos repo
between 2026-01-23 and 2026-01-29.  The branch was never merged and
was deleted from origin after the contents were captured here.  These
patches are the durable archive — apply with `git am < <patch>.patch`
from the repo root, or `git apply` if you want a working-tree-only
preview.

Branch diverged from opt-dev's ancestor (`9abe5ce`, 2026-01-23) just
before all of the FreeForm SrfType=14, per-ray-status, and
STOP/PAUSE→`SetRayFail` work landed on opt-dev.  Direct `git am`
applies will conflict heavily against current opt-dev — these are
**reference patches**, not drop-in cherry-picks.

## Patch contents

| Patch | Lines | What it does | Conflict risk on opt-dev |
|---|---|---|---|
| 0001 (b1b167f) | 1257 | MonSrf rewrite: PAUSE removal + `ieee_arithmetic` / `ieee_exceptions` INF/NaN handling + redundant-computation cuts + Newton/Secant fast paths before Brent's | High — MonSrf has the FreeForm `SrfType_FreeForm = 14` dispatch on opt-dev |
| 0002 (ebcdc22) | 34 | Comment-only: conic-radicand-check documentation in MonSrf + ConSrf | Low |
| 0003 (49c94cd) | 410 | ConSrf updates: IEEE flag clear/check, near-tangent branch, cancellation-aware quadratic-root selection, `INTENT(IN/OUT)` modernization, `Conic_Sag` helper returning `ieee_value(quiet_nan)` on out-of-aperture | Medium — ConSrf evolved on opt-dev but not in the radicand region |
| 0004 (a9e687a) | 527 | Extract `SolveConicIntersection` as shared helper between ConSrf and MonSrf, inline-friendly | High — would tangle with FreeForm dispatch |
| 0005 (3573d7e) | 3466 | GridSrf large rewrite | Very high — opt-dev's jGridSrf mapping + FreeForm grid-component work touched everything |

## Worth-porting patterns (bounded scope, low conflict)

The following patterns from 0003 (49c94cd) are the most concrete
accuracy improvements and the cleanest to lift out by hand:

1. **IEEE-flag clear + post-check at routine boundaries** in ConSrf and
   MonSrf — catches silent NaN propagation that would otherwise sail
   through downstream `dot_product` / norm calls without diagnostic.
   ~10-20 lines per routine.

   ```fortran
   use, intrinsic :: ieee_arithmetic
   use, intrinsic :: ieee_exceptions
   ...
   ! Top of routine: clear flags
   call ieee_set_flag(ieee_invalid, .false.)
   call ieee_set_flag(ieee_overflow, .false.)
   call ieee_set_flag(ieee_divide_by_zero, .false.)
   ...
   ! Before LROK = .true.:
   call ieee_get_flag(ieee_invalid, ieee_invalid_flag)
   if (ieee_invalid_flag) go to 98   ! treat as miss
   ```

2. **Discriminant pre-check before `SQRT(k2)`** — returns
   `ieee_value(sag, ieee_quiet_nan)` instead of letting the runtime
   fault.  Different failure mode from what opt-dev's `SetRayFail`
   already handles (bracket/max-iter); the radicand-negative case is
   currently silent.

   ```fortran
   k2 = b*b - 4d0*a*c
   if (k2 < 0d0) go to 98   ! no real roots — surface doesn't exist here
   ```

3. **Cancellation-aware quadratic-root selection.** When
   `|b| ≈ √k2`, the standard quadratic formula loses precision.
   Switch to the conjugate form `2c / (b + sign(b)·√k2)` to preserve
   digits.  Textbook IEEE-arithmetic fix.

   ```fortran
   sqrtk2 = sqrt(k2)
   abs_b_minus_sqrtk2 = abs(abs(b) - sqrtk2)
   if (abs_b_minus_sqrtk2 < TOL_CANCEL * abs(b)) then
     ! Cancellation: use conjugate form
     L = 2d0 * c / (-b - sign(sqrtk2, b))
   else
     ! Standard quadratic formula
     ksqrt = -0.5d0 * (b + sign(sqrtk2, b))
     ...
   end if
   ```

4. **Near-tangent fallback.** When `k2 < TOL_TANGENT·b²` (≈14 digits
   between the two roots), use the double-root formula `L = -b/(2a)`
   instead of the standard two-root selection.  Avoids amplifying
   noise into the spread between roots that are numerically identical.

   ```fortran
   if (k2 < TOL_TANGENT * b * b) then
     L = -b / (2d0 * a)   ! double root — near-tangent
   else
     ...                   ! standard two-root case
   end if
   ```

`TOL_TANGENT`, `TOL_CANCEL`, `TOL_ZERO`, `TOL_NORMAL`, `TOL_GEOM` are
all defined in patch 0003's prelude on `Constants` (constants.f90)
— pull those definitions over too.

## Patterns **NOT** worth porting (superseded / too risky)

- **PAUSE removal in MonSrf** — opt-dev independently converted all
  PAUSE / STOP statements in surfsub to graceful `SetRayFail` returns
  via commit 51bbf1f (per-ray status tracking).  No re-do needed.
- **`SolveConicIntersection` extraction** — would conflict with
  opt-dev's `FreeFormSrf` / `MonomialEval` / `MonomialD2` module
  structure.  The duplication between ConSrf and MonSrf is small;
  not worth the integration cost.
- **GridSrf rewrite** — opt-dev's `jGridSrf` mapping and the FreeForm
  grid-component work overlap with most of the touched code.  Way
  too high a conflict surface for a refactor.

## How to apply a pattern

Working from a patch file in this directory:

```bash
# Inspect the patch
less docs/Archive/dev_optimization_surfsub/0003-Updated-ConSrf-...patch

# Try a clean apply (will fail due to conflicts)
git am --3way docs/Archive/dev_optimization_surfsub/0003-Updated-ConSrf-...patch

# Recommended: open the patch, locate the pattern you want, and edit
# the current ConSrf / MonSrf / etc. by hand.  Then add a smoke test
# that wasn't there before (e.g. a Rx that exercises the cancellation
# path).
```
