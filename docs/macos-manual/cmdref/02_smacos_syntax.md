## SMACOS syntax conventions

SMACOS is the engine as a callable Fortran subroutine.  One entry
point executes any Part I command:

```fortran
CALL SMACOS(command, CARG, DARG, IARG, LARG, RARG,
&           OPDMat, RaySpot, RMSWFE, PixArray)
```

- `command` — `CHARACTER(LEN=132)`, the command string (same names
  and minimum-match rules as the CLI).
- `CARG(9)` — `CHARACTER(LEN=32)` arguments (file names, element
  names, option keywords).
- `DARG(9)` — `REAL*8` arguments; `IARG(9)` — `INTEGER`;
  `LARG` — `LOGICAL`; `RARG(9)` — `REAL*4`.
- `OPDMat(mpts,mpts)`, `RaySpot(mRay,2)`, `RMSWFE`,
  `PixArray(mPix,mPix)` — result buffers the engine fills for output
  commands (OPD, SPOT, PIXel, ...).

Arguments land in the same order the CLI would have prompted for
them: where the CLI asks "Enter prescription file name", SMACOS
expects it in `CARG(1)`; numeric prompts map to `DARG`/`IARG` slots in
prompt order.  The manual's Section 9 documents the mapping for the
commonly-scripted commands, and the CLI dialog itself (run the command
interactively once) is the reliable way to see what a command wants.

Before the first call the runtime must be initialized with the model
size — the programmatic equivalent of the CLI's startup prompt:

```fortran
call macos_init_all(256)
macos_init = .true.
```

A complete minimal driver (load, trace, OPD RMS) is listed in manual
Section 9.5; `smacos_dvr` in the build tree is a ready-made
interactive driver linked against `libsmacos.a`, useful for testing
call sequences before embedding them.

Two practical notes:

- SMACOS shares no state with any other copy of the engine in the
  process; sizes fixed at `macos_init_all` time bound all later work.
- Commands that would page output or draw graphics in the CLI run
  silently or write files instead; result data comes back through the
  buffer arguments.

The higher-level bindings (mmacos, pymacos) wrap this same entry
point plus the finer-grained routines of `macos_api_mod.F90`; if you
are writing new Fortran, consider calling `macos_api_mod` directly —
it is the supported programmatic surface and each Part II entry names
the routine involved.
