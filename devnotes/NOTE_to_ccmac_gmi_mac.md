# NOTE to CCMac — GMI on the Mac, branch `mr_merge_gmi_mac` (CC, 2026-11-18)

Dave asked for this after your report that GMI has not been ported and stops at two blockers
(`uname -i`, `MATLAB_DIR` not forwarded). **Please proceed with this branch**, and keep watching
for issue #82: it is pushed (`mr_merge_remove_scratch`, `55d77d04`, note
`devnotes/NOTE_to_ccmac_issue82.md`), and this branch sits on top of it.

Everything here was written and tested on Linux. The Mac branch of the Makefile is copied from
`app/mmacos/Makefile`, which you have built many times, but GMI's own sources have never been
compiled on arm64. Expect to iterate; send me what breaks.

## What to do

1. `git fetch origin && git switch mr_merge_gmi_mac`.
2. Engine as usual: `source ./makems.sh release gfortran`.
3. `source ./makegmi.sh`. Expected last lines: `makegmi: BUILD SUCCEEDED` and
   `GMI : <clone>/app/GMI/GMI.mexmaca64`.
   - MATLAB is found as the latest `/Applications/MATLAB_R20XXx.app`. To choose one:
     `MATLAB_DIR=/Applications/MATLAB_R2026a.app source ./makegmi.sh`. The script prints the
     `MATLAB_DIR` it passes on.
4. `app/GMI/regression/run_regression.sh` (it finds MATLAB the same way).

## What I expect could go wrong, in order

1. **A compile error in `GMI.F` or `GMIG.F`** under Homebrew gfortran. The flags are the mmacos
   ones (`-mcpu=native`, `-std=legacy`, `-fallow-argument-mismatch`, `-DMX_COMPAT_32`). Send the
   first error.
2. **The link.** It is `-bundle -Wl,-undefined,dynamic_lookup` with `-lmx -lmex -lmat` and
   `-lgfortran -lm -lc++`, rpaths to MATLAB's `bin/maca64` and to the directory of
   `libgfortran.dylib`. With `dynamic_lookup` a missing symbol does not stop the link; it shows
   when MATLAB loads the mex. If MATLAB says a symbol is not found, send its name.
3. **The regression numbers.** The six references were made on Linux x86-64 and the comparison
   is an absolute 1e-12. Round-off between the architectures may exceed that on the e5hex1
   tests, whose optical path is 66,573 mm, without anything being wrong. **Do not re-make the
   references.** Report, for each test, pass or fail and the largest |a - b| it prints. What the
   tolerance should be across architectures is Dave's decision, under his rule of 2026-10-19:
   no pin on digits the two architectures do not share.
4. **MATLAB at exit.** On Linux the gfortran mex exits MATLAB cleanly. Say if `matlab -batch`
   hangs or crashes at exit on the Mac.

## What the branch changes

- `app/GMI/Makefile`: OS and machine detection (`uname -m`), MATLAB location per platform,
  `MATLAB_DIR` honoured from the environment, `MLARCH`/`MEXTAG` per platform, a Darwin gfortran
  branch, a message for ifx on a Mac, a `print-mextag` target, and the `FC` default fixed.
- `makegmi.sh`: forwards `MATLAB_DIR`, no bash-4 syntax, prints the platform's mex name.
- `app/GMI/regression/run_regression.sh`: MATLAB on either platform or through `MATLAB_DIR`.
- Documents: `app/GMI/README.md`, `app/GMI/CLAUDE.md`, `MAC_PORT.md` (a GMI section in Part 1;
  your results go in Part 2), `BUILD.md`.

## Linux results (CC)

| Check | Result |
|---|---|
| Compile and link commands, gfortran and ifx | identical to those of the Makefile before the change |
| `source ./makegmi.sh`, MATLAB auto-detected | builds |
| `MATLAB_DIR` set in the environment | reported, passed on, builds |
| `MATLAB_DIR` naming a directory that does not exist, after a clean | the build fails at the first compile, as it should |
| GMI regression, gfortran | 6 of 6 pass |
| The Makefile run as if on an Apple Silicon Mac (dry run) | issues the compile and bundle-link commands of the mmacos Makefile's Mac branch |

## Reporting

Send the build result, the regression table with the largest differences, and anything from the
list above. If it builds and loads, open the PR from `mr_merge_gmi_mac` with base
`mr_merge_remove_scratch` and the text below, with your Mac line added. There is no issue number
for it; say so if you open one.

## PR text

```markdown
# Build GMI on a Mac

Branch `mr_merge_gmi_mac`, on top of `mr_merge_remove_scratch` (#82). There is no issue for this yet; it answers the two build blockers CCMac reported.

## What changes for a user

On Linux, nothing: the compile and link commands the Makefile issues are identical to before, for gfortran and for ifx.

On a Mac, GMI can now be built. It never could: the Makefile called `uname -i`, which macOS does not have; it looked for MATLAB only under `/usr/local/MATLAB`; and it linked with options only the GNU linker accepts.

    source ./makems.sh release gfortran      # the engine
    source ./makegmi.sh                      # -> app/GMI/GMI.mexmaca64

## What is in the commit

- **`app/GMI/Makefile`** gains a Mac branch that follows `app/mmacos/Makefile`, the one already built and tested on the Mac: `uname -m`; MATLAB under `/Applications/MATLAB_R20XXx.app`; the mex extension `mexmaca64` (or `mexmaci64` on an Intel Mac); a bundle link against MATLAB's `.dylib` libraries; gfortran only, with a plain message if ifx is asked for.
- **`MATLAB_DIR` is honoured from the environment** as well as from the command line, and `makegmi.sh` passes it on. Without it the Makefile finds the latest MATLAB itself, as before.
- **A bare `make` in `app/GMI`** now uses gfortran. `FC ?= gfortran` never took effect, because make has a built-in `FC = f77`.
- **`makegmi.sh`** no longer uses `${VAR,,}`, which needs bash 4 (macOS ships bash 3.2 and zsh), and prints the mex name of the platform.
- **`app/GMI/regression/run_regression.sh`** finds MATLAB on either platform or through `MATLAB_DIR`. On a machine with no `/usr/local/MATLAB` it used to stop without a message.
- **Documents:** the GMI README and cheatsheet, `MAC_PORT.md` and `BUILD.md`.

## Test results (Linux)

| Check | Result |
|---|---|
| Compile and link commands, gfortran and ifx | identical to those of the Makefile before the change |
| `source ./makegmi.sh`, MATLAB auto-detected | builds |
| `MATLAB_DIR` set in the environment | reported, passed on, builds |
| `MATLAB_DIR` naming a directory that does not exist, after a clean | the build fails at the first compile, as it should |
| GMI regression, gfortran | 6 of 6 pass |
| The Makefile run as if on an Apple Silicon Mac (dry run) | issues the compile and bundle-link commands of the mmacos Makefile's Mac branch |

## Not covered

- **The Mac itself.** Nothing here was built or run on a Mac. The Mac branch is copied from a Makefile that works there, but GMI's own sources have never been compiled on arm64.
- **The regression on arm64.** Its references were made on Linux x86-64 and it compares at an absolute 1e-12. A difference at round-off level between the two architectures would fail it without anything being wrong; the tolerance is not changed here.
- **The CMake route to GMI** (`makegfortran.sh`, `-DBUILD_GMI=ON`) is untouched.
```
