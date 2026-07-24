# Porting MACOS to macOS (Apple Silicon) — prep notes

Forward-looking guide for building the engine, the design layer, and the
bindings on an Apple Silicon Mac.  Nothing here is exercised on Linux;
the build files stay Linux/Windows-only until these arms are applied and
validated on real hardware (marked **[validate on first Mac build]**).

## The one big fact: gfortran end-to-end
- **gfortran is the compiler** on macOS.  `ifx` has no Apple Silicon
  (arm64) support, so it is not an option here — and gfortran is already
  the **MathWorks-supported Fortran compiler on Linux** (confirmed in the
  ifx bug-report exchange) and the project's gfortran build is clean
  (`makegfortran.sh`, `macos_api_mod` ported gfortran-clean).  ifx remains
  the opt-in *alternative* on Linux/Windows only.
- **The mmacos mex is a DIRECT gfortran link**, not MATLAB's `mex` driver
  (`MACOS_resources/mmacos/Makefile:122` links against `libmx`/`libmex`/
  `libmat` itself).  So MathWorks's per-platform "supported Fortran
  compiler" matrix is **moot** — there is no `mex -setup FORTRAN` gate to
  satisfy.  The Mac mex is a *link-recipe* port (extension + lib path +
  Darwin link flags), not a compiler-support problem.

## Toolchain setup
- **Homebrew**: `brew install gcc cmake cairo readline`.  Homebrew's `gcc`
  package provides `gfortran-NN` (versioned, e.g. `gfortran-14`) — there is
  usually no bare `gfortran` symlink, so pass it explicitly to CMake.
- **C compiler = Apple clang** (default) is fine *with* the implicit-decl
  downgrade (below); the gfortran+clang mix is the standard Mac
  scientific-computing setup and links cleanly.  (Homebrew `gcc-NN` for C
  also works if you prefer real GCC.)
- **MATLAB**: use the **native arm64** build (R2023b+), not the Rosetta
  x86_64 one — the mex arch must match (`maca64`).
- **Filesystem — clone onto a CASE-SENSITIVE APFS volume.**  macOS default
  APFS is case-INSENSITIVE; this tree relies on case-distinct names
  (CLAUDE.md: never `surfsub.f` beside `surfsub.F`) and they **collide on
  checkout** otherwise.  `diskutil apfs addVolume disk1 "Case-sensitive
  APFS" devwork` and clone there.
- **Debugger = `lldb`**, not gdb.  gdb on Apple Silicon is painful
  (weak arm64-macOS support + code-signing); lldb is Apple-native and
  works with gfortran binaries.  The "gdb to pinpoint engine bugs"
  workflow (memory) becomes lldb — same idea, different driver.
- **Editors**: VS Code (CMake Tools + Modern Fortran/`fortls`, all native
  on Mac) as the IDE + Claude Code (native; `.claude/` config + memory
  port as-is).  BBEdit optional for fast grep/quick edits.

## Engine + library (CMakeLists.txt)

### (a) Fortran compiler default
`CMakeLists.txt:8-13` searches for `ifx`; on Mac it won't be found, so
CMake's default Fortran detection picks gfortran — but Homebrew installs
`gfortran-14`, not bare `gfortran`.  Configure with it explicitly:
```
cmake -B build_release_gfortran -G Ninja -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_Fortran_COMPILER=gfortran-14 \
      -DBUILD_SMACOS=ON -DBUILD_MACOS=OFF   # library-only for the design layer
```

### (b) GNU compiler-flag arm — add the APPLE clang downgrade
In the `elseif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")` block
(`CMakeLists.txt:80-87`), append the clang flag.  `-DLNXOS` is KEPT on
macOS: the LNXOS sites are POSIX and Mac shares them.  **[validate on
first Mac build]** the three sites: `sunsub.F:117`, `utilsub_c.c:24,189`
(the C ones already co-check `__linux__`; clang auto-defines `__APPLE__`,
so add `|| defined(__APPLE__)` there only if a Mac path is needed).
```cmake
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  # Linux + macOS gfortran/gcc.  macOS/Apple-Silicon: gfortran is the only
  # path (ifx has no arm64) and the MathWorks-supported Linux compiler.
  # -DLNXOS kept on macOS: the LNXOS sites are POSIX (Mac shares them).
  set(BASE_FFLAGS "-cpp -ffixed-line-length-132 -fd-lines-as-comments -fallow-argument-mismatch -std=legacy -DPGPLOT -DLNXOS -fPIC")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -fbounds-check -fcheck=all")
  set(CMAKE_Fortran_FLAGS         "${BASE_FFLAGS}")
  set(CMAKE_C_FLAGS "-DLNXOS -fPIC")
  if(APPLE)
    # Apple clang ERRORS on the legacy C's implicit declarations -> downgrade.
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-implicit-function-declaration -Wno-implicit-int")
  endif()
endif()
```

### (c) RPATH — drop the GNU-ld-only flag on macOS
The existing `elseif(NOT WIN32 AND CMAKE_Fortran_COMPILER_ID MATCHES
"GNU")` RPATH arms (macos / smacos_dvr / GMI targets, ~`CMakeLists.txt:412,
479`) set `BUILD_RPATH/INSTALL_RPATH = ${CMAKE_Fortran_IMPLICIT_LINK_
DIRECTORIES}` (correct on Mac — that's the Homebrew gcc lib dir) but then
add `-Wl,--disable-new-dtags`, which **Mac's ld64 rejects**.  Guard it:
```cmake
  elseif(NOT WIN32 AND CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set_target_properties(<tgt> PROPERTIES
        BUILD_RPATH   "${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES}"
        INSTALL_RPATH "${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES}")
    if(NOT APPLE)
      target_link_options(<tgt> PRIVATE -Wl,--disable-new-dtags)  # GNU ld only
    endif()
```
On macOS CMake emits `@rpath` install names automatically from INSTALL_RPATH.

### (d) Giza / interactive graphics — DEFER
Giza = Cairo + X11; X11 isn't native on macOS (needs XQuartz).  But the
**library + design layer don't need it** — `libsmacos.a` links the
`pgplotdummy.F` stubs, headless.  Build `BUILD_SMACOS=ON BUILD_MACOS=OFF`
first; port the Giza/XQuartz path only when the interactive `macos` CLI
is wanted on Mac.

## mmacos mex (MACOS_resources/mmacos/Makefile)

### (a) Arch / MEX extension (replace `Makefile:45-51`)
`uname -i` is invalid on macOS — use `uname -m`:
```makefile
# Architecture / MEX extension
OS      := $(shell uname -s)
MACHINE := $(shell uname -m)
ifeq ($(OS),Darwin)
  ifeq ($(MACHINE),arm64)
    ARCH := maca64
    MEXTAG := mexmaca64
  else
    ARCH := maci64
    MEXTAG := mexmaci64
  endif
else
  ARCH := glnxa64
  MEXTAG := mexa64
endif
```
(Then use `$(ARCH)` for the MATLAB `bin/`/`extern/lib/` subdir instead of
the hardcoded `glnxa64` in LDFLAGS/POST_LIBS.)

### (b) Darwin gfortran branch (add after the `else ifeq ($(FC),gfortran)`
arm, `Makefile:98-112`).  Mac uses a **bundle** with dynamic_lookup, not
`-shared`/version-script/`--no-undefined`; MATLAB libs are `.dylib` (plain
`-lmx`, not GNU-ld `-l:libmx.so`); `-mcpu=native` not `-march=native`;
`-lc++` not `-lstdc++`:
```makefile
else ifeq ($(OS),Darwin)            # macOS / Apple Silicon -- gfortran only
  GFLIB := $(shell dirname $$(gfortran -print-file-name=libgfortran.dylib))
  FFLAGS := -fPIC -cpp -ffixed-line-length-132 \
            -fno-omit-frame-pointer -O2 -mcpu=native \
            -fallow-argument-mismatch -std=legacy \
            -DMX_COMPAT_32 \
            -I$(MATLAB_DIR)/extern/include \
            -I$(SMACOS_OBJS) \
            -J$(SRCPATH)
  LDFLAGS := -fPIC -O -bundle -Wl,-undefined,dynamic_lookup \
             -Wl,-rpath,$(MATLAB_DIR)/bin/$(ARCH) \
             -Wl,-rpath,$(GFLIB)
  POST_LIBS := -L$(MATLAB_DIR)/bin/$(ARCH) -lmx -lmex -lmat \
               -L$(GFLIB) -lgfortran -lm -lc++
```
Notes / **[validate on first Mac build]**:
- The `-reentrancy=none` ifx thread-at-exit fix is irrelevant here (no
  Intel runtime) — one fewer worry.
- The `clear mex` / batch-exit hang is platform-independent; keep the
  `exit(0)` rule in batch scripts.
- If `mexFunction` isn't found at load, add an exported-symbols list
  (`-Wl,-exported_symbols_list,<file with _mexFunction>`); usually
  `-undefined dynamic_lookup` on a `-bundle` is enough.

## pymacos on Mac — the smoother binding
f2py + Homebrew gfortran is well-trodden on macOS.  Build `libsmacos.a`
(headless, above), then build pymacos with `FC=gfortran-14`,
`-DMACOS_BUILD_DIR=<...>/build_release_gfortran`.  Because the design
layer is state-as-data (PLAN_DESIGN_LAYER §3) and the Python port is cheap
by construction, **leading with pymacos on Mac sidesteps the mex link
work entirely** — decide up front whether Mac == mmacos or pymacos-first.

## Validation checklist (when a Mac is in hand)
1. `libsmacos.a` builds headless (gfortran-14 + clang, `BUILD_SMACOS=ON
   BUILD_MACOS=OFF`); `nm` shows `__macos_api_mod_MOD_init`.
2. Review the 3 LNXOS sites for Mac correctness (sunsub.F, utilsub_c.c).
3. pymacos f2py build + `pytest` (the no-MATLAB path) green.
4. mmacos mex links (`.mexmaca64`) and `run_mmacos_tests.sh fast` green
   (seat-check already opt-in via MM_SEAT_CHECK).
5. lldb a deliberate fault to confirm the debug workflow.
6. (later) Giza/XQuartz for the interactive `macos` CLI.

## VALIDATED ON HARDWARE — 2026-07-22 (Apple M-series, macOS 26.5, MATLAB R2024a)

**Everything below actually built + ran on Dave's Mac.** The guide above was
mostly right; the deltas are recorded here.

### Toolchain actually installed
- `brew install gcc cmake ninja` → **GCC 16** (so `gfortran-16`, and Homebrew
  *does* provide a bare `gfortran` symlink → gcc 16). cmake 4.4.0, ninja 1.13.2.
- MATLAB **R2024a is arm64** (`bin/maca64/*.dylib`, `extern/lib/maca64/fexport.map`).
  Uses **JPL network license** (`cae-lm-mlm*.jpl.nasa.gov:7282`) — first `matlab`
  launch is slow (ServiceHost init); not a hang.
- pymacos needs **Python ≥3.13** (`find_package(Python 3.13 REQUIRED)`): installed
  `brew install python@3.13` → venv `~/.venv/pymacos` + `numpy 2.5.1`.
- `~/dev` is case-INSENSITIVE; no `.f`/`.F` collisions today, built fine as-is.

### Build files edited (all `if(APPLE)` / `Darwin`-guarded — Linux/Win unchanged)
- **`macos/CMakeLists.txt`**: GNU flag block — `if(APPLE)` adds
  `-Wno-implicit-function-declaration -Wno-implicit-int` to `CMAKE_C_FLAGS`.
  smacos_dvr GNU RPATH arm — `-Wl,--disable-new-dtags` wrapped in `if(NOT APPLE)`.
  (Compiler id is now **`IntelLLVM`**, not `Intel`; GNU block ~L128; giza
  `add_subdirectory` is unconditional but **XQuartz is present** so it configures
  fine, and the `smacos` LIBRARY links only pgplotdummy — no giza edit needed.)
- **`mmacos/CMakeLists.txt`**: `elseif(APPLE)` platform arm (maca64 / `.mexmaca64`
  / plain `-lmx -lmex -lmat`); `add_library(mmacos MODULE …)` on APPLE (MODULE =
  `-bundle` automatically); Apple link opts = `-Wl,-undefined,dynamic_lookup` +
  MATLAB rpath (drop version-script/no-undefined/dtags); `-mcpu=native` on arm64.
- **`mmacos/Makefile`** (the path `run_mmacos_tests.sh` actually drives): Darwin
  MATLAB location (`/Applications/MATLAB_R*.app`), `uname -m` arch → maca64/mexmaca64,
  `MACOS_BUILD_DIR` default → `build_release_gfortran`, a Darwin gfortran arm
  (bundle+dynamic_lookup, `-lc++`, `-mcpu=native`, libgfortran on rpath), and
  `FC` default via `ifeq ($(origin FC),default)` (GNU make's built-in `FC=f77`
  otherwise defeats `?=`). Added a `print-mextag` helper target.
- **`mmacos/run_mmacos_tests.sh`**: Darwin MATLAB detection; ext-aware mex path
  via `make print-mextag`; and **`export MACOS_HOME=…/macos_f90`** (see gotcha).
- **`pymacos/src/cmake/CMakeLists.txt`**: `elseif(APPLE)`→MACOS64; gfortran
  compiler arm; a gfortran `FFlags` arm (`-ffree-form -O2 -fPIC -cpp
  -fallow-argument-mismatch -std=legacy -mcpu=native` — the Intel arm's
  `-fpp/-qmkl/-names/-assume` are invalid); and the preprocess custom-command
  uses `-E -cpp -P -o <out>` on Mac (bare `-P` makes gfortran *compile*, needing
  kinds.mod — Intel's `-P` preprocesses-only).
- **`pymacos/src/cmake/source/pymacos_f2py.f90`**: two mixed-length char-array
  constructors (`(/'m','cm',…/)`) → `[character(len=4):: …]` (gfortran rejects
  mixed lengths; ifx tolerated).

### THE load-bearing gotcha (cost the most time): `macos_param.txt` + MEX = SIGSEGV
`init` (→`macos_init_all`) loads `macos_param.txt`, searched in: cwd, exe dir,
`$MACOS_HOME`, then a compile-time `build_loc` (literal `'BUILD_LOC'` — never
resolves). `matlab -batch` runs from the mmacos root, which has no
`macos_param.txt`. **A MISSING param file makes the engine abort fatally, and
inside the MEX host that abort is a hard SIGSEGV** (crash dump on the MCR
thread, empty fault stack — looks like a codegen bug, is not). Fix: point
`MACOS_HOME` at `macos_f90/`. A standalone Fortran driver shows a clean
"could not be found, bye!" for the same condition — only the MEX turns it into
a segfault. First engine call (`init`) is where it dies, so it masquerades as
"the whole mex is broken."

### Results
- **libsmacos.a**: clean build, gfortran-16 + Apple clang, 155 `macos_api_mod`
  symbols, 42 `.mod`. `build_release_gfortran/` (release) + `build_debug_gfortran/`
  (`-fcheck=all`, used to diagnose the param-file issue).
- **mmacos**: `src/mmacos.mexmaca64` (Mach-O arm64 bundle, exports `_mexfunction_`
  = exactly what maca64 `fexport.map` wants). `./run_mmacos_tests.sh fast` =
  **226 pass, 0 real fail, 3 skipped** (tRunCompare assumes `SegMirMaker` built
  under `build_release_ifx/` — impossible on arm64; not a port defect).
- **pymacos**: `pymacosf90.cpython-313-darwin.so` builds; import + init + load +
  trace + opd on `Rx_Cass_FarField.in` → RMS OPD 1.29e-12 m, 47084 rays. Full
  `pytest` not run this pass.
- Debugger: `lldb` (not gdb) is the Mac tool; a debug libsmacos + tiny standalone
  driver is the fastest way to isolate engine faults *outside* MATLAB.

### macos INTERACTIVE CLI — also built (2026-07-22, one fix)
The full `macos` exe (needs giza→cairo→pixman graphics) builds + runs on arm64.
Only ONE fix was needed beyond the library: the vendored **pixman** was hardcoded
x86-64 (`pixman-x86.c` carries x86 inline asm; `%eax`/`=a` constraints clang
rejects on arm64). Fix in `macos_f90/giza/src/subprojects/pixman/CMakeLists.txt`:
arch-gate via `CMAKE_SYSTEM_PROCESSOR` — on non-x86, don't `#define USE_SSE2/
USE_SSSE3` in the generated config.h and drop `pixman-sse2.c`/`pixman-ssse3.c`
from the sources. `pixman-x86.c` is fully `#if`-guarded on those macros, so it
compiles to an empty TU; pixman falls back to its generic C paths (correct, just
no SIMD accel — ARM NEON left off, fine for occasional CLI plotting). cairo/giza/
pgplot then compiled clean. **XQuartz must be installed** (it is) for the X11
window; giza's `find_package(X11 REQUIRED)` resolves against `/opt/X11`.
Verified: `bin/macos` (Mach-O arm64) loads Rx_Cass_FarField, traces, OPD RMS
1.286e-12 m (== pymacos), reaches the giza `/xw` device prompt.
Build: add `-DBUILD_MACOS=ON` to the library configure below, `--target macos`.
(Bundled `libreadline.a` is absent → arrow-key history disabled, non-fatal.)

### Build commands (copy-paste, verified)
```sh
# library (release)
cd ~/dev/macos && cmake -B build_release_gfortran -G Ninja -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_Fortran_COMPILER=gfortran-16 \
  -DBUILD_SMACOS=ON -DBUILD_MACOS=OFF -DBUILD_GMI=OFF -DBUILD_MMACOS=OFF -DBUILD_SMACOS_DVR=OFF
cmake --build build_release_gfortran --target smacos

# mmacos mex (Makefile path — outputs src/mmacos.mexmaca64)
cd ~/dev/MACOS_resources/mmacos && make FC=gfortran        # MACOS_BUILD_DIR defaults to build_release_gfortran
MACOS_HOME=~/dev/macos/macos_f90 ./run_mmacos_tests.sh fast # (runner exports MACOS_HOME itself)

# pymacos (into ~/.venv/pymacos with numpy>=2.2)
cd ~/dev/MACOS_resources/pymacos/src/cmake
PATH="$HOME/.venv/pymacos/bin:$PATH" cmake -B build_mac -S . \
  -DCMAKE_Fortran_COMPILER=gfortran-16 \
  -DMACOS_BUILD_DIR=$HOME/dev/macos/build_release_gfortran \
  -DPython_EXECUTABLE=$HOME/.venv/pymacos/bin/python
PATH="$HOME/.venv/pymacos/bin:$PATH" cmake --build build_mac
```
