# SLSQP — Kraft's Sequential Least Squares Programming

Replaces NPSOL as MACOS's constrained-optimization back end.  Branch
`sls-dev` (branched off opt-dev 2026-06-04).

## Phase 0 — NPSOL usage survey

NPSOL is invoked in only one place: `design_cons_optim.F:382` via
`np_optim_dvr`, called from `macos_cmd_loop.inc:3662` (start_optim
path) and `macos_cmd_loop.inc:3965` (restore path).  `setbeam.inc`
also uses `npoptn` + `npsol` for the set-beam initialization.

### Problem shape

From `design_cons_optim.F:95`:

    Integer, parameter :: nclin=0, ncnln=0, nrowA=1, nrowj=1, nrowR=MAXVAR

* `nclin = 0`   — **no linear constraints**
* `ncnln = 0`   — **no nonlinear constraints**
* Bounds `bl`, `bu` are the only constraint

In other words, MACOS's "constrained" optimization is **bound-constrained
only**.  This is the simplest possible SLSQP regime: a smooth objective
+ box bounds.

### Objective callback (`npsol_funcobj`, design_cons_optim.F:422)

Signature: `(mode, n, a, obf, grd, nstat)`
* `a(n)`   in  — variable vector
* `obf`    out — scalar objective (squared-residual sum for WFE / Zern / Beam targets)
* `grd(n)` out — gradient (finite difference)

Body calls `smacos_compute_perturbed_np(...)` to evaluate the optical
system at the candidate `a`, then assembles `obf` and `grd` from the
target type (WFE / WFE_ZMODE / BEAM).

### NPSOL features MACOS does NOT use

* nonlinear constraints (`ncnln=0`)
* linear constraints  (`nclin=0`)
* warm starts (`COLD` mode is hardwired)
* multiplier returns (computed by NPSOL but not propagated back to
  caller in `np_optim_dvr`)
* `npoptn` runtime options other than verify-level, print-level,
  linesearch-tol (all set to suppress noise)

### NPSOL → SLSQP mapping

| NPSOL                              | SLSQP (Kraft, reverse-comm) |
|------------------------------------|-----------------------------|
| `npsol(n,...,bl,bu,funcobj,...)`   | `slsqp(0,0,1,n,x,xl,xu,...)` |
| `funcobj(mode,n,a,obf,grd,nstat)`  | driver loop computes `f,g` directly |
| `Verify level = -1` etc            | (drop — SLSQP doesn't need these) |
| nclin / ncnln / clamda             | (drop — bound-only) |

## Phase 1 — Vendor Kraft SLSQP

Source: Kraft's original Fortran translated by Schittkowski (DLR), as
distributed with SciPy under `scipy/optimize/slsqp/`.  Public domain.

Files to vendor under `macos_f90/slsqp/`:
* `slsqp_optmz.f`   (~1900 lines, main + lsq + ldp + nnls helpers)
* `slsqp.h.in` or equivalent constants

License notice: `LICENSE.txt` (public domain).

## Phase 2 — `slsqp_optim_dvr`

New file `macos_f90/design_slsqp_optim.F` mirroring
`design_cons_optim.F`.  Same external signature as `np_optim_dvr` so
the call sites in `macos_cmd_loop.inc` don't change semantically.

Internals: reverse-comm loop:
```fortran
mode = 0
do
  call slsqp(m=0, meq=0, la=1, n=len_a, x=aparams, xl=bl, xu=bu, &
             f=obf, c=cval, g=grd, a=cjac, acc=dopt_tol,           &
             iter=nitrs_dopt, mode=mode, w=w, l_w=l_w, jw=jw, l_jw=l_jw)
  select case (mode)
  case (1)   ! request f, g
    call smacos_compute_perturbed_sls(... aparams, obf, grd, ...)
  case (-1)  ! request c, a (constraint values + Jacobian; unused here)
    continue
  case (0)   ! converged
    exit
  case default  ! error
    rtn_flag = mode
    exit
  end select
end do
```

## Phase 3 — Dispatch swap

* `macos_cmd_loop.inc:3661/3964`: rename `np_optim_dvr` →
  `slsqp_optim_dvr`; drop `#ifdef USE_NPSOL` guards (SLSQP always
  available now)
* `setbeam.inc:103-105`: parallel changes (npoptn → SLSQP options)
* CMakeLists: drop `USE_NPSOL` option, drop `npsol_src_dir`, drop
  links to `libnpsol.a / liblapacklib.a / libblaslib.a`
* Build scripts: drop `npsol` argument from
  `makems.sh / makeall.sh / makegfortran.sh / makegmi.sh`
* `MACOS_resources` companion: drop `libnpsol.a / liblapacklib.a /
  libblaslib.a` auto-detect from `GMI/Makefile`

## Phase 4 — Regression

A/B vs NPSOL on every constrained Rx:
* `opt_example.in` (BeamOnly + bound-constrained rigid body)
* Any prescription in `ZGD_test_files/` that exercises constrained
  optimization
* Pymacos tests/test_calib_constrained.py (new, exercises bounds via
  programmatic setters)

Pass: final objective within `1e-6` of NPSOL result; iteration count
comparable order of magnitude.

## Phase 5 — Bindings (deferred)

`m.calib_set_bounds(...)` etc. — natural now that Phase 1d's
introspection plan can wrap SLSQP's bound-setting directly.
