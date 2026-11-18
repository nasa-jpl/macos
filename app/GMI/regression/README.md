# GMI regression suite

Six lightweight tests that catch the most common GMI-breaking change
classes.  Each test loads a known Rx, calls the mex via the standard
`call_GMI.m` entry point, and compares the result against a committed
`.mat` reference.

> **GMI build choice (2026-05-25):** GMI builds with either ifx or
> gfortran.  Both pass all six regression tests with bit-identical
> numeric results, both exit MATLAB cleanly, and `run_regression.sh`
> trusts the exit code directly.
>
> - **gfortran (default):** `source ./makegmi.sh` from the repository
>   root.  Mex at `app/GMI/GMI.mexa64`.  Requires
>   `build_release_gfortran/` from `makems.sh release gfortran` first.
> - **ifx:** `source ./makegmi.sh ifx`.  Requires
>   `build_release/` from `makems.sh release` first.  app/GMI/Makefile
>   links the single-threaded Intel runtime (`-reentrancy=none`,
>   pulls in `libifcore.so.5` instead of `libifcoremt.so.5`); the
>   default multi-threaded variant kept worker threads alive across
>   `clear mex` unloads and SIGSEGV'd at MATLAB process exit by
>   waking those threads to a now-unmapped callback.  See
>   app/GMI/Makefile for the full back-story.

## Run

```bash
./run_regression.sh           # all tests, exits non-zero on any fail
./run_regression.sh --bootstrap   # (re)generate the .mat references
                                  # after an intentional change
```

Needs MATLAB (any `R20xx` under `/usr/local/MATLAB/`) and a built
`GMI.mexa64` in the parent directory.

## What's covered

| Test | Rx | Channel | Catches |
|---|---|---|---|
| 01 smoke (e2e_pie) | `e2e_pie_met` | nominal | mex won't load, missing symbol, signature regression |
| 02 nominal repro (e2e_pie) | e2e_pie | nominal x2 | `SetToNominalSettings` state-drift |
| 03 Zernike response (e2e_pie) | e2e_pie | Z4 on Elt 8 (Zernike-typed) | Zernike-channel apply path, `iNode=4..nzern+3` indexing, ELSE-branch reset |
| 04 smoke (e5hex1) | `Rx_e5hex1.in` (FreeForm chain) | nominal | basic FreeForm loading |
| 05 nominal repro (e5hex1) | e5hex1 | nominal x2 | FreeForm save/restore (`pFF/xFF/yFF/zFF`, `pData/xData/...`) |
| 06 Zernike response (e5hex1) | e5hex1 | Z4 on Elt 8 (Zernike-typed) | Zernike apply on a mixed FreeForm+Zernike chain |
| 07 tolerance refuses | none | no engine call | a tolerance loosened until it passes everything |

Tolerance against the stored references: `|a - b| <= 32 * eps * L`, where
`L` is the prescription's chief-ray path length and `eps` the machine
epsilon, so `eps * L` is the round-off of the path (`lib/ref_tolerance.m`):
2.3e-13 m for e2e_pie (L = 31.96 m) and 4.9e-10 mm for e5hex1
(L = 69,121 mm).  An OPD is
the difference of two paths, so it is good to the round-off of the path.
The references were made on Linux x86-64; the first Mac (arm64) run,
2026-11-18, differed from them by 8 and 12 of these units on e2e_pie and by
8 and 6 on e5hex1.  The fixed `1e-12` used until then passed the first and failed
the second only because one deck is in metres and the other in
millimetres.  The two calls of one build in tests 02 and 05 are still
compared at tolerance 0.

**References re-made 2026-11-17** (gfortran mex; the ifx mex reproduces
them bit for bit).  The previous set dated from before 2026-10-10, when
the engine's default OPD reference changed from the aperture mean to the
chief ray; GMI takes the engine default, so every map had moved by one
constant, its "Average OPD" (e5hex1: -8.2589e-06 mm; e2e_pie: -2.184e-10 m),
and tests 02, 03, 05 and 06 failed.  With the old reference forced back
(`UseChfRay4OPD= N` in the deck) the e2e_pie maps matched the old files to
3e-14; the e5hex1 maps matched to 8.7e-11 mm, six units in the last place
of the 66,573 mm path, which is round-off from engine changes since and
was above the fixed `1e-12` tolerance of the time.  The current references are chief-referenced.

## Files

```
regression/
├── README.md
├── run_regression.sh            entry point
├── regression_main.m            top-level runner
├── bootstrap_reference.m        regenerate references
├── tests/
│   ├── test01_smoke_e2e_pie.m
│   ├── test02_nominal_repro_e2e_pie.m
│   ├── test03_zern_response_e2e_pie.m
│   ├── test04_smoke_e5hex1.m
│   ├── test05_nominal_repro_e5hex1.m
│   ├── test06_zern_response_e5hex1.m
│   └── test07_tolerance_refuses.m
├── lib/
│   ├── compare_within.m         absolute-tolerance equality check
│   ├── ref_tolerance.m          the tolerance against a reference, from the path length
│   ├── init_e2e_pie.m           param struct for the e2e_pie Rx
│   └── init_e5hex1.m            param struct for the e5hex1 Rx
├── reference/                   committed `.mat` ground truth
│   ├── nominal_e2e_pie.mat
│   ├── zern_response_e2e_pie.mat
│   ├── nominal_e5hex1.mat
│   └── zern_response_e5hex1.mat
└── Rx/                          prescription files seen by SMACOS OLD
    ├── Rx_e5hex1.in
    └── e2e_pie_met.in          (copied from parent on first run)
```

## When a test fails

The summary table tells you which test + the max |diff| vs reference.

- **All tests fail with mex error**: rebuild GMI.  `source ./makegmi.sh` (gfortran) or `source ./makegmi.sh ifx`, from the repository root.
- **Just `nominal_repro_*` fails (round-trip non-zero)**: state-drift bug; check `ObtainNominalSettings` / `SetToNominalSettings` in `GMI.F` for fields not being snapshotted (e.g., the FreeForm `pFF/xFF/yFF/zFF` class of fix).
- **All tests fail with `vs reference`**: an intentional change has shifted the numbers.  Inspect the diff, regenerate via `./run_regression.sh --bootstrap`, commit the new `reference/*.mat`.
- **`zern_response_*` fails with "perturbation produced zero response"**: GMI is forcing `SrfType=8` on the perturbed element but the trace dispatch in `propsub.F` is silently no-op'ing.  Two known root causes (both fixed on release-candidate, watch for regression):
  - `ZernTypeL(iElt)=0` at the apply site — GMI must set a non-zero default when forcing SrfType=8 (see `GMI.F` ~line 1066).
  - Rx uses legacy "Malacara" naming that `ParseZernType` no longer recognizes — see the `Mala`-prefix alias in `elt_mod.F` ParseZernType.
- **One specific `*_response_*` fails with `vs reference`**: numeric drift in the perturbation-apply path.

## Adding a test

Drop a `testNN_*.m` in `tests/`, add it to the `tests` cell array in
`regression_main.m`, regenerate references via the bootstrap script,
commit both the new test and the new `reference/*.mat`.
