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
