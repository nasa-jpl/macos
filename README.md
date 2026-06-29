# MACOS — Modeling and Analysis for Controlled Optical Systems

Jet Propulsion Laboratory

MACOS is an optical ray-tracing application. SMACOS is its static library
counterpart. Graphics are provided by **Giza**, a Cairo-based drop-in
replacement for PGPLOT.

> **PGPLOT removed on release-candidate.** Giza is now the only PGPLOT-API
> provider. The `[giza|pgplot]` argument previously accepted by build scripts
> is gone.
>
> **SLSQP is the default constrained-optimization back end** (Kraft 1988,
> BSD-licensed, vendored under `macos_f90/slsqp/`). It builds and runs out
> of the box — no extra source tree, no license. The legacy NPSOL path is
> still available behind `-DUSE_NPSOL=ON` for sites that have the Stanford
> SOL source; default builds do not require it.

---

## Table of Contents
- [Directory Layout](#directory-layout)
- [Clone the Repositories](#clone-the-repositories)
- [Linux](#linux)
- [Windows](#windows)
- [Build Options Reference](#build-options-reference)
- [Build Outputs Reference](#build-outputs-reference)
- [Building pymacos (Python Interface)](#building-pymacos-python-interface)
- [Building mmacos (MATLAB MEX Interface)](#building-mmacos-matlab-mex-interface)
- [Troubleshooting](#troubleshooting)

---

## Directory Layout

Both repositories must be siblings under the same parent directory.
Replace `<username>` with your login name throughout.

```
MACOS_resources/                   (top-level resources repo)
    GMI/                           (GMI MATLAB mex interface)
    mmacos/                        (mmacos MATLAB mex interface -- CMakeLists.txt here)
    pymacos/
        src/
            macos/                 (main MACOS/SMACOS source tree -- CMakeLists.txt here)
                macos_f90/         (Fortran source files)
                build/             (CMake build output -- generated)
```

---

## Clone the Repositories

```bash
mkdir -p ~/dev && cd ~/dev

git clone git@github.com:nasa-jpl/macos.git
git -C macos checkout joint-dev

git clone git@github.com:nasa-jpl/MACOS_resources.git
git -C MACOS_resources checkout dr-dev
```

---

## Linux

### Prerequisites

| Tool | Purpose |
|------|---------|
| Intel oneAPI HPC Toolkit | Fortran compiler (`ifx`) and C compiler (`icx`/`gcc`) |
| gcc / g++ | C sources and Giza graphics library |
| CMake ≥ 3.21 | Build system |
| X11 development headers | Interactive graphics window |
| MATLAB *(optional)* | Required only for `GMI.mexa64` or `mmacos.mexa64` |

Install X11 headers on Ubuntu/Debian:
```bash
sudo apt install libx11-dev
```

For a gfortran build use `makegfortran.sh` instead — Intel oneAPI is not required.

### Setup

1. Install Intel oneAPI HPC Toolkit from
   [intel.com/oneapi](https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit.html)

2. Source the Intel environment (build scripts do this automatically):
   ```bash
   source /opt/intel/oneapi/setvars.sh
   ```

3. Build the bundled readline library (enables arrow-key history in macos):
   ```bash
   cd macos_f90/readline-8.2
   ./configure && make
   cd ../..
   ```

### Building with Scripts

All scripts run from `~/dev/macos/` and accept `[debug|release]` options (default: release).

#### Build all targets at once

```bash
cd ~/dev/macos
source ./makeall.sh              # ifx, release (most common)
source ./makeall.sh debug        # ifx, debug
source ./makeall.sh npsol        # ifx, release, NPSOL also (opt-in; SLSQP always on)
source ./makegfortran.sh         # gfortran, release
```

#### Build individual targets

| Script | Targets built |
|--------|--------------|
| `source ./makems.sh` | `macos` + `libsmacos.a` |
| `source ./makesd.sh` | `smacos_dvr` (also builds macos/smacos if needed) |
| `source ./makegmi.sh` | `GMI.mexa64` (requires macos + smacos built first) |

All three accept the same `[debug|release]` options. `makesd.sh` and
`makegmi.sh` must use matching options to the preceding `makems.sh`
(or `makeall.sh`) call so they target the same build directory.

### Building with CMake Directly

#### Build macos + smacos (default targets)
Builds the interactive `macos` executable and the `libsmacos.a` static library.
GMI and mmacos are **not** included unless explicitly enabled with `-DBUILD_GMI=ON` / `-DBUILD_MMACOS=ON`.
```bash
rm -rf build_release
cmake -B build_release -G "Ninja" -DCMAKE_BUILD_TYPE=Release
cmake --build build_release -j$(nproc)
```

#### Build macos + smacos + mmacos (MATLAB MEX interface)
Builds `macos`, `libsmacos.a`, and `mmacos.mexa64` in one step.
The build system auto-detects MATLAB under `/usr/local/MATLAB/R*`.
`mmacos.mexa64` is output to `build_release/bin/`.
```bash
rm -rf build_release
cmake -B build_release -G "Ninja"   -DCMAKE_BUILD_TYPE=Release   -DBUILD_MACOS=ON   -DBUILD_SMACOS=ON   -DBUILD_MMACOS=ON
cmake --build build_release -j$(nproc)
```

#### Build macos + smacos + GMI (MATLAB GMI mex interface)
Builds `macos`, `libsmacos.a`, and `GMI.mexa64` in one step.
GMI is the legacy MATLAB interface; mmacos is the newer replacement.
`GMI.mexa64` is output to the GMI source directory.
```bash
rm -rf build_release
cmake -B build_release -G "Ninja"   -DCMAKE_BUILD_TYPE=Release   -DBUILD_MACOS=ON   -DBUILD_SMACOS=ON   -DBUILD_GMI=ON   -DMATLAB_ROOT=/usr/local/MATLAB/R2025b
cmake --build build_release -j$(nproc)
```

#### Build macos + smacos + GMI + mmacos (all MATLAB interfaces)
Builds everything including both MATLAB interfaces in one step.
```bash
rm -rf build_release
cmake -B build_release -G "Ninja"   -DCMAKE_BUILD_TYPE=Release   -DBUILD_MACOS=ON   -DBUILD_SMACOS=ON   -DBUILD_GMI=ON   -DBUILD_MMACOS=ON   -DMATLAB_ROOT=/usr/local/MATLAB/R2025b
cmake --build build_release -j$(nproc)
```

#### Build only smacos library
```bash
rm -rf build_release
cmake -B build_release -G "Ninja" -DCMAKE_BUILD_TYPE=Release -DBUILD_MACOS=OFF -DBUILD_SMACOS=ON
cmake --build build_release -j$(nproc)
```

#### Build only macos executable
```bash
rm -rf build_release
cmake -B build_release -G "Ninja" -DCMAKE_BUILD_TYPE=Release -DBUILD_MACOS=ON -DBUILD_SMACOS=OFF
cmake --build build_release -j$(nproc)
```

#### Build macos + smacos + smacos_dvr
```bash
rm -rf build_release
cmake -B build_release -G "Ninja" -DCMAKE_BUILD_TYPE=Release -DBUILD_SMACOS_DVR=ON
cmake --build build_release -j$(nproc)
```

#### Build smacos lean (minimal, no FITS, no graphics)
```bash
rm -rf build_lean
cmake -B build_lean -DBUILD_SMACOS_LEAN=ON
cmake --build build_lean -j$(nproc)
```

#### Build with NPSOL constrained optimization (opt-in)
SLSQP is the default constrained-optim back end and is always built in;
you only need this if you also want NPSOL available for A/B comparison or
to honor existing journals that hardwired the NPSOL code path.
Requires a separately-licensed Stanford SOL NPSOL source tree at
`macos_f90/npsol/`.
```bash
rm -rf build_release_npsol
cmake -B build_release_npsol -G "Ninja" -DCMAKE_BUILD_TYPE=Release -DUSE_NPSOL=ON
cmake --build build_release_npsol -j$(nproc)
```

#### Build with GMI MATLAB mex interface
```bash
rm -rf build_release
cmake -B build_release -G "Ninja" -DCMAKE_BUILD_TYPE=Release -DBUILD_GMI=ON \
      -DMATLAB_ROOT=/usr/local/MATLAB/R2025b
cmake --build build_release -j$(nproc)
```

#### Debug build
```bash
rm -rf build_debug
cmake -B build_debug -G "Ninja" -DCMAKE_BUILD_TYPE=Debug
cmake --build build_debug -j$(nproc)
```

### Rebuilding After Source Changes

Re-run the same script — CMake recompiles only changed files.
To force a clean rebuild delete the build directory first:

```bash
rm -rf ~/dev/macos/build_release
source ~/dev/macos/makeall.sh
```

### Shell Aliases

Add to `~/.bashrc` (or `~/.bash_aliases`), then `source ~/.bashrc`:

```bash
alias macos='~/dev/macos/build_release/bin/macos'
alias smacos_dvr='~/dev/macos/build_release/bin/smacos_dvr'
```

For a debug build:
```bash
alias macos='~/dev/macos/build_debug/bin/macos'
alias smacos_dvr='~/dev/macos/build_debug/bin/smacos_dvr'
```

For MATLAB, add the GMI directory to the MATLAB path in `startup.m` or via
the MATLAB **Set Path** dialog:

```matlab
addpath('~/dev/MACOS_resources/GMI')
```

### Linux Build Output Locations

| Target | Location |
|--------|---------|
| `macos` executable | `build_release/bin/macos` |
| `smacos_dvr` executable | `build_release/bin/smacos_dvr` |
| `libsmacos.a` static library | `build_release/lib/libsmacos.a` |
| `GMI.mexa64` MATLAB mex | `~/dev/MACOS_resources/GMI/GMI.mexa64` |
| `mmacos.mexa64` MATLAB mex | `build_release/bin/mmacos.mexa64` |

### Linux Build Directory Names

| Command | Build directory |
|---------|----------------|
| `makeall.sh` | `build_release/` |
| `makeall.sh debug` | `build_debug/` |
| `makeall.sh npsol` | `build_release_npsol/` |
| `makeall.sh debug npsol` | `build_debug_npsol/` |
| `makegfortran.sh` | `build_release_gfortran/` |
| `makegfortran.sh debug` | `build_debug_gfortran/` |

---

## Windows

### Prerequisites

| Tool | Confirmed working | Confirmed NOT working | Purpose |
|------|-------------------|------------------------|---------|
| Intel oneAPI HPC Toolkit | **2025.3.3 (icx) / 2025.3.2 (ifx)** | **2026.0** (clang-22 `mmintrin.h` regression — see Troubleshooting) | Fortran (`ifx`) and C (`icx`) compilers |
| Visual Studio | **2022 Community** (stable) | **2026 Insiders** (preview channel, paired with Intel 2026.0 above) | Linker and Windows SDK |
| CMake | **4.1.1** | **4.2.x, 4.3.x** (fail to set `CMAKE_Fortran_PREPROCESS_SOURCE` for `ifx` during ABI detection) | Build system |
| MATLAB | R2024b or later | — | Required only for `mmacos.mexw64` / `GMI.mexw64` |

> All other version combinations (other Intel 2025.x point releases, VS 2022
> non-Community editions, CMake versions outside the two ranges above) are
> **untested** — not confirmed working, not confirmed broken.

> **Verified working versions (exact, tested combination — only this
> combination has actually been confirmed; other version numbers below
> are untested, not "known good"):**
> - `icx` **2025.3.3** / `ifx` **2025.3.2** (Intel oneAPI HPC Toolkit,
>   installed under `compiler\2025.3`, with `latest` symlinked to it).
>   **No other Intel oneAPI version has been verified** — not 2024.2,
>   not other 2025.x point releases, not 2026.x. Use this exact toolkit
>   release unless you've independently confirmed a different one works.
> - **Visual Studio 2022 Community** (stable release — NOT a preview/Insiders
>   channel), with the **"Desktop development with C++"** workload installed
>   (provides `cl.exe`, MSVC toolset, Windows SDK)
> - **CMake 4.1.1**
>
> **Confirmed broken:** Intel oneAPI **2026.0** (bundles clang 22) paired
> with Visual Studio **2026 Insiders** reproducibly fails building
> Giza/Cairo on Windows — see
> [Troubleshooting: mmintrin.h errors](#troubleshooting) below. This is a
> regression in Intel's bundled clang-22 headers in clang-cl/MSVC-compatible
> mode (`mmintrin.h` `__builtin_shufflevector` errors), not something fixable
> from this project's CMakeLists.txt. Earlier drafts of this README
> recommended 2026.0; that guidance was never actually verified against a
> clean install and should be disregarded.
>
> Separately: CMake 4.2.x and 4.3.x have a known incompatibility with newer
> ifx releases regardless of version (missing `CMAKE_Fortran_PREPROCESS_SOURCE`
> during ABI detection). CMake 4.1.1 avoids this entirely.

### Setup

1. Install [Visual Studio 2022 Community](https://visualstudio.microsoft.com/vs/older-downloads/)
   (the stable release — avoid Insiders/preview channels), selecting
   **"Desktop development with C++"** under the Workloads tab. This installs
   `cl.exe`, the MSVC toolset, and the Windows SDK, all of which `ifx`/`icx`
   depend on even when you never invoke `cl.exe` directly. **Verify after
   install** by opening a fresh terminal and running `where cl` — if it
   doesn't return a path, the C++ workload didn't install correctly and
   everything downstream (including Fortran compiler identification in
   CMake) will fail in confusing ways.

2. Install [Intel oneAPI HPC Toolkit](https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit.html)
   — specifically the **2025.3** release (icx 2025.3.3 / ifx 2025.3.2 is the
   exact build verified above). Avoid 2026.0, which has a confirmed
   Giza/Cairo build failure (see Troubleshooting). Other versions are
   untested — if you must use a different one, expect to troubleshoot.

3. Install [CMake 4.1.1](https://cmake.org/files/v4.1/) — select
   **"Add CMake to system PATH"** during install. Avoid CMake 4.2.x/4.3.x.

4. Open the **Intel oneAPI Command Prompt** from the Start menu
   (this sets up `ifx` and `icx` on your PATH automatically — it also runs
   Visual Studio's `vcvars64.bat` internally, which is why step 1 must be
   done first and fully working).

5. **Verify the full toolchain before building anything**, in a fresh
   Intel oneAPI Command Prompt:
   ```cmd
   icx --version
   ifx --version
   cmake --version
   where cl
   ```
   All four must succeed and show the versions above. If `where cl` fails
   here even though it worked in step 1, your Intel oneAPI install may not
   be detecting Visual Studio correctly — see
   [Troubleshooting: Visual Studio was not found](#troubleshooting).

### Building with the Batch Script

Run from the **Intel oneAPI Command Prompt**:

```cmd
makems.bat              REM release (default)
makems.bat debug        REM debug
makems.bat npsol        REM release + NPSOL
makems.bat debug npsol  REM debug + NPSOL
```

### Building with CMake Directly

#### Build macos + smacos (default targets)
Builds `macos.exe` and `smacos.lib`. mmacos is **not** included unless
explicitly enabled with `-DBUILD_MMACOS=ON`.
```cmd
rmdir /s /q build
cmake -B build -G "Ninja" -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

#### Build macos + smacos + mmacos (MATLAB MEX interface)
Builds `macos.exe`, `smacos.lib`, and `mmacos.mexw64` in one step.
The build system auto-detects MATLAB under `C:\Program Files\MATLAB\R*`.
`mmacos.mexw64` is output to `build\bin\`.
Run from the **Intel oneAPI Command Prompt**:
```cmd
rmdir /s /q build
cmake -B build -G Ninja ^
  -DCMAKE_BUILD_TYPE=Release ^
  -DBUILD_MACOS=ON ^
  -DBUILD_SMACOS=ON ^
  -DBUILD_MMACOS=ON
cmake --build build
```

> **Note:** GMI is also supported on Windows — see the GMI build example below.

#### Build macos + smacos + GMI (MATLAB GMI mex interface)
Builds `macos.exe`, `smacos.lib`, and `GMI.mexw64` in one step.
Run from the **Intel oneAPI Command Prompt**:
```cmd
rmdir /s /q build
cmake -B build -G Ninja ^
  -DCMAKE_BUILD_TYPE=Release ^
  -DBUILD_MACOS=ON ^
  -DBUILD_SMACOS=ON ^
  -DBUILD_GMI=ON
cmake --build build
```

#### Build macos + smacos + GMI + mmacos (all MATLAB interfaces)
Builds everything including both MATLAB interfaces in one step.
Run from the **Intel oneAPI Command Prompt**:
```cmd
rmdir /s /q build
cmake -B build -G Ninja ^
  -DCMAKE_BUILD_TYPE=Release ^
  -DBUILD_MACOS=ON ^
  -DBUILD_SMACOS=ON ^
  -DBUILD_GMI=ON ^
  -DBUILD_MMACOS=ON
cmake --build build
```

#### Build only smacos library
```cmd
rmdir /s /q build
cmake -B build -G "Ninja" -DCMAKE_BUILD_TYPE=Release -DBUILD_MACOS=OFF -DBUILD_SMACOS=ON
cmake --build build
```

#### Build only macos executable
```cmd
rmdir /s /q build
cmake -B build -G "Ninja" -DCMAKE_BUILD_TYPE=Release -DBUILD_MACOS=ON -DBUILD_SMACOS=OFF
cmake --build build
```

#### Build smacos lean (minimal, no FITS, no graphics)
```cmd
rmdir /s /q build_lean
cmake -B build_lean -DBUILD_SMACOS_LEAN=ON
cmake --build build_lean
```

#### Build with NPSOL constrained optimization (opt-in)
SLSQP is the default and is always built; only enable NPSOL if you have
the Stanford SOL source available at `macos_f90/npsol/`.
```cmd
rmdir /s /q build_release_npsol
cmake -B build_release_npsol -G "Ninja" -DCMAKE_BUILD_TYPE=Release -DUSE_NPSOL=ON
cmake --build build_release_npsol
```

#### Debug build
```cmd
rmdir /s /q build_debug
cmake -B build_debug -G "Ninja" -DCMAKE_BUILD_TYPE=Debug
cmake --build build_debug
```

### Rebuilding After Source Changes

```cmd
rmdir /s /q build
cmake -B build -G "Ninja" -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

### Windows Build Output Locations

| Target | Location |
|--------|---------|
| `macos.exe` executable | `build\bin\macos.exe` |
| `smacos.lib` static library | `build\lib\smacos.lib` |
| `GMI.mexw64` MATLAB GMI mex | `build\bin\GMI.mexw64` |
| `mmacos.mexw64` MATLAB mex | `build\bin\mmacos.mexw64` |

> **Note:** smacos_dvr is not supported on Windows. GMI and mmacos are both supported.

---

## Build Options Reference

Pass these to cmake with `-D<OPTION>=ON/OFF`:

| Option | Default | Description |
|--------|---------|-------------|
| `BUILD_MACOS` | ON | Build the interactive `macos` executable |
| `BUILD_SMACOS` | ON | Build the `libsmacos.a` / `smacos.lib` static library |
| `BUILD_SMACOS_LEAN` | OFF | Minimal smacos (no FITS, no graphics). Mutually exclusive with `BUILD_SMACOS` |
| `BUILD_SMACOS_DVR` | OFF | Build the `smacos_dvr` test driver (Linux only) |
| `BUILD_GMI` | OFF | Build the MATLAB GMI mex interface (Linux + Windows) |
| `BUILD_MMACOS` | OFF | Build the MATLAB mmacos mex interface (Linux + Windows) |
| `MMACOS_SRC_DIR` | auto-detected | Path to mmacos source directory containing `mmacos_mex.F` |
| `MMACOS_FC` | (inherit) | Override Fortran compiler for mmacos only (e.g. `gfortran` while main build uses `ifx`) |
| `USE_NPSOL` | OFF | Also build the NPSOL constrained-optim back end (opt-in). SLSQP is always built and is the default. Requires `macos_f90/npsol/` source tree |
| `MATLAB_ROOT` | `/usr/local/MATLAB/R2025b` | MATLAB installation path (for GMI build) |

---

## Build Outputs Reference

| File | Platform | Description |
|------|----------|-------------|
| `build_release/bin/macos` | Linux | Interactive MACOS executable |
| `build/bin/macos.exe` | Windows | Interactive MACOS executable |
| `build_release/lib/libsmacos.a` | Linux | SMACOS static library |
| `build/lib/smacos.lib` | Windows | SMACOS static library |
| `build_release/bin/smacos_dvr` | Linux | SMACOS test driver |
| `build_release/bin/GMI.mexa64` | Linux | MATLAB GMI mex interface |
| `build\bin\GMI.mexw64` | Windows | MATLAB GMI mex interface |
| `build_release/bin/mmacos.mexa64` | Linux | MATLAB mmacos mex interface |
| `build/bin/mmacos.mexw64` | Windows | MATLAB mmacos mex interface |

---

## Building pymacos (Python Interface)

pymacos is the Python interface to SMACOS. It must be built **after** SMACOS.

### Linux

#### Prerequisites
- Python 3.13 or later with NumPy 2.2 or later
- Intel oneAPI (ifx, icx)
- SMACOS already built in `build_release/`

#### Build
```bash
cd ~/dev/MACOS_resources/pymacos/src/cmake
rm -rf build && mkdir build && cd build
cmake -DCMAKE_C_COMPILER=icx -DCMAKE_CXX_COMPILER=icpx -DCMAKE_Fortran_COMPILER=ifx -S ..
make
```

Output: `pymacos/src/pymacos/pymacosf90.cpXXX-linux.so`

Run tests from `pymacos/src`:
```bash
cd ~/dev/MACOS_resources/pymacos/src
python -m pytest ../tests/ -v
```

### Windows

#### Prerequisites
- Python 3.13 or later with NumPy 2.2 or later
- Intel oneAPI Command Prompt (ifx)
- SMACOS already built in `macos\build\`

#### Build
```cmd
cd pymacos\src\cmake
rmdir /s /q build
mkdir build
cd build
cmake -G "NMake Makefiles" -S ..
nmake
```

Output: `pymacos\src\pymacos\pymacosf90.cp314-win_amd64.pyd`

Run tests from `pymacos\src`:
```cmd
cd pymacos\src
python -m pytest ..\tests\ -v
```

> **Note:** Tests must be run from `pymacos\src` so Python can find the `pymacos` package.

---

## Building mmacos (MATLAB MEX Interface)

mmacos is the MATLAB MEX interface to SMACOS. It is the MATLAB sibling of
pymacos — both wrap the same `macos_api_mod` module compiled into
`libsmacos.a` / `smacos.lib`. This means a single SMACOS library backs both
Python and MATLAB simultaneously.

mmacos is built from its own `CMakeLists.txt` located in
`MACOS_resources/mmacos/`. It can be built two ways:

- **Integrated**: called automatically from the main `CMakeLists.txt` via
  `BUILD_MMACOS=ON` — SMACOS and mmacos are built in one step
- **Standalone**: invoked separately from the `mmacos/` directory, pointing
  at a pre-built SMACOS tree via `MACOS_BUILD_DIR`

### Linux

#### Prerequisites
- MATLAB R2024b or later installed under `/usr/local/MATLAB/`
- Intel oneAPI (ifx) or gfortran
- SMACOS already built (or build together with `BUILD_MMACOS=ON`)

#### Build mmacos together with macos + smacos (recommended)

This is the simplest approach — one cmake invocation builds everything:

```bash
cd ~/dev/macos
rm -rf build_release
cmake -B build_release -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_MACOS=ON \
  -DBUILD_SMACOS=ON \
  -DBUILD_MMACOS=ON
cmake --build build_release -j$(nproc)
```

The build system automatically locates MATLAB under `/usr/local/MATLAB/R*`
and picks the latest installed version. To specify a particular version:

```bash
cmake -B build_release -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_MACOS=ON \
  -DBUILD_SMACOS=ON \
  -DBUILD_MMACOS=ON \
  -DMATLAB_DIR=/usr/local/MATLAB/R2025b
cmake --build build_release -j$(nproc)
```

#### Build mmacos standalone (after SMACOS is already built)

Use this when you want to rebuild mmacos without rebuilding SMACOS.
`-DMMACOS_SRC_DIR` is only needed if mmacos is not at the standard
relative path (`../../../mmacos` from `pymacos/src/macos/`):

```bash
cd ~/dev/MACOS_resources/mmacos
rm -rf build
cmake -B build -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DMACOS_BUILD_DIR=~/dev/macos/build_release
cmake --build build -j$(nproc)
```

#### Use gfortran for mmacos (while the main build uses ifx)

The Makefile default is gfortran for mmacos (it exits MATLAB more cleanly
than ifx on some Linux configurations). To use gfortran for mmacos only
while keeping ifx for the rest:

```bash
cmake -B build_release -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_MACOS=ON \
  -DBUILD_SMACOS=ON \
  -DBUILD_MMACOS=ON \
  -DMMACOS_FC=gfortran
cmake --build build_release -j$(nproc)
```

#### Output

`build_release/bin/mmacos.mexa64`

### Windows

#### Prerequisites
- MATLAB R2024b or later installed under `C:\Program Files\MATLAB\`
- Intel oneAPI HPC Toolkit 2026.0 (ifx is the only supported compiler on Windows)
- Visual Studio 2022 or 2026 with "Desktop development with C++"
- CMake 4.1.1
- SMACOS already built (or build together with `BUILD_MMACOS=ON`)
- Open an **Intel oneAPI Command Prompt** (sets up ifx on PATH)

#### Build mmacos together with macos + smacos (recommended)

Run from the **Intel oneAPI Command Prompt**:

```cmd
cd C:\path\to\macos\pymacos\src\macos
rmdir /s /q build
cmake -B build -G Ninja ^
  -DCMAKE_BUILD_TYPE=Release ^
  -DBUILD_MACOS=ON ^
  -DBUILD_SMACOS=ON ^
  -DBUILD_MMACOS=ON
cmake --build build
```

Replace `C:\path\to\` with your actual paths. The build system automatically
locates MATLAB under `C:\Program Files\MATLAB\R*` and picks the latest version.

To specify a particular MATLAB version:

```cmd
cmake -B build -G Ninja ^
  -DCMAKE_BUILD_TYPE=Release ^
  -DBUILD_MACOS=ON ^
  -DBUILD_SMACOS=ON ^
  -DBUILD_MMACOS=ON
  -DMATLAB_DIR="C:\Program Files\MATLAB\R2025b"
cmake --build build
```

#### Build mmacos standalone (after SMACOS is already built)

```cmd
cd C:\path\to\MACOS_resources\mmacos
rmdir /s /q build
cmake -B build -G Ninja ^
  -DCMAKE_BUILD_TYPE=Release ^
  -DMACOS_BUILD_DIR=C:\path\to\macos\pymacos\src\macos\build
cmake --build build
```

#### Output

`build\bin\mmacos.mexw64`

### Using mmacos in MATLAB

After building, add the output directory to the MATLAB path and set the
working directory to a folder containing `macos_param.txt`:

```matlab
% Add mmacos to MATLAB path
addpath('C:\path\to\build\bin')   % Windows
addpath('~/dev/macos/build_release/bin')  % Linux

% Set working directory to folder containing macos_param.txt
% macos_param.txt defines MACOS model dimensions and is required
% by mmacos('init', ...) before any ray tracing can be performed
cd('C:\path\to\MACOS_resources\mmacos')   % Windows
cd('~/dev/MACOS_resources/mmacos')        % Linux

% Initialize MACOS with a model size (number of elements)
mmacos('init', 128)

% Load a prescription file
mmacos('load_rx', 'path/to/prescription.in')

% Run a trace
mmacos('trace_rays')

% Get OPD map
opd = mmacos('opd', 11)    % OPD at element 11

% Get complex field
cf = mmacos('complex_field', 9)   % Complex wavefront at element 9
```

> **Important:** `macos_param.txt` must be present in the MATLAB current
> working directory before calling `mmacos('init', ...)`. This file defines
> MACOS model parameters (array sizes, element limits, etc.). Without it,
> MATLAB will exit without an error message.

### Running the mmacos smoke test

```matlab
% Quick smoke test using a Cassegrain prescription
cd('C:\path\to\MACOS_resources\mmacos')
addpath('C:\path\to\build\bin')
test_mmacos('path/to/Rx_Cass_FarField.in')
exit(0)   % explicit exit needed to prevent MATLAB hang in -batch mode
```

Or via CTest (after building):

```cmd
ctest --test-dir build -R smoke
ctest --test-dir build -R unittest
```

### mmacos compiler notes

| Platform | Default compiler | Alternative |
|----------|-----------------|-------------|
| Linux | gfortran (matches Makefile default; cleaner MATLAB exit) | ifx via `-DMMACOS_FC=ifx` |
| Windows | ifx (only supported option via oneAPI) | — |

On Linux, gfortran is preferred for mmacos because it exits MATLAB more
cleanly than ifx when the MEX is unloaded. Both produce numerically
identical results.

---

## Troubleshooting

**`CMAKE_Fortran_PREPROCESS_SOURCE` error on Windows**
Two distinct causes produce this same error message — check both:

1. **CMake version**: CMake 4.2.x/4.3.x fail to auto-populate this variable
   during Fortran ABI detection for recent `ifx` releases. Use **CMake 4.1.1**.
2. **Missing `cl.exe`** (more common cause): `ifx`'s ABI/identification probe
   silently fails if Visual Studio's C++ workload isn't installed or isn't
   being found by `setvars.bat`. Symptoms: the configure log shows
   `-- The Fortran compiler identification is unknown` *before* the
   `CMAKE_Fortran_PREPROCESS_SOURCE` error. Fix: run `where cl` in your
   terminal — if it returns nothing, install the **"Desktop development
   with C++"** workload via the Visual Studio Installer (Modify → Workloads
   → Desktop and Mobile), then open a fresh Intel oneAPI Command Prompt.

**`WARNING: Visual Studio was not found in a standard install location`**
(printed by Intel's `setvars.bat`) — Intel's setup script only recognizes
VS install paths matching `2017`/`2019`/`2022` by default. If your VS
install lives elsewhere (e.g. a preview/Insiders build under a numeric
folder like `...\Microsoft Visual Studio\18\...` instead of a year),
set the install-dir variable yourself before calling `setvars.bat`:
```cmd
set "VS2022INSTALLDIR=C:\Program Files\Microsoft Visual Studio\<your-folder>\<Edition>"
"C:\Program Files (x86)\Intel\oneAPI\setvars.bat" --force
```
To avoid retyping this every session, set it permanently with `setx`:
```cmd
setx VS2022INSTALLDIR "C:\Program Files\Microsoft Visual Studio\<your-folder>\<Edition>"
```
Then open a **new** terminal for it to take effect.

**`mmintrin.h` / `__builtin_shufflevector` errors when building Cairo/Giza on Windows**
This is a **confirmed regression in Intel oneAPI 2026.0's bundled clang-22
backend** when operating in clang-cl/MSVC-compatible mode — not something
fixable via compiler flags (`/arch:SSE2` etc. do not help; the flag reaches
the compiler but the headers are internally broken regardless). The
project's verified-working toolchain is **Intel oneAPI 2025.3.x**, which
does not bundle clang 22 and does not exhibit this bug.
**Fix: install Intel oneAPI 2025.3 instead of 2026.0** (side-by-side
install is supported; point `setvars.bat`/your terminal at `compiler\2025.3`
or whichever version's `latest` symlink resolves there). Re-verify with
`icx --version` — it should report `2025.3.x`, not `2026.0`.

**`pgplot_mod` not found**
`pgplotdummy.F` must be present in the `macos_f90/` directory. It provides
stub graphics module definitions needed by the smacos library.

**`ifx: command not found`**
Source the Intel oneAPI environment first:
- Linux: `source /opt/intel/oneapi/setvars.sh`
- Windows: Open the Intel oneAPI Command Prompt from the Start menu

**`libreadline.a` not found warning (Linux)**
Arrow-key history in the macos prompt will be disabled. To fix:
```bash
cd macos_f90/readline-8.2 && ./configure && make
```

**MATLAB root not found (GMI or mmacos build)**
Set the MATLAB installation path explicitly:
```bash
# Linux
cmake -B build_release -DMATLAB_DIR=/usr/local/MATLAB/R2025b \
      -DBUILD_MMACOS=ON -DCMAKE_BUILD_TYPE=Release
```
```cmd
REM Windows
cmake -B build -DMATLAB_DIR="C:\Program Files\MATLAB\R2025b" ^
      -DBUILD_MMACOS=ON -DCMAKE_BUILD_TYPE=Release
```

**`LNK2019` unresolved external symbols on Windows**
Ensure you are using the latest `CMakeLists.txt` which includes
`windows_stubs.c` generated at configure time to provide Windows-compatible
stubs for Linux-only functions.

**mmacos: MATLAB exits without error after `mmacos('init', ...)`**
`macos_param.txt` is missing from the MATLAB current working directory.
This file is required by `param_mod_init` during SMACOS initialization.
Copy it from `MACOS_resources/mmacos/macos_param.txt` to your working
directory, or `cd` to the mmacos directory before calling `init`:
```matlab
cd('path/to/MACOS_resources/mmacos')
mmacos('init', 128)
```

**mmacos: "Gateway function is missing" on Windows**
The MEX file was built with an older `mmacos/CMakeLists.txt`. Rebuild with
the latest version which correctly exports `mexFunction` via the linker.

**mmacos: "unknown command" for valid commands**
The string dispatch in `mexFunction` is not receiving the command correctly.
This is a known issue with some compiler/MATLAB version combinations. Ensure
you are using the latest `mmacos/CMakeLists.txt` which uses `MX_COMPAT_32`
and the correct MATLAB API symbol aliases.

**`error #7002: Error in opening the compiled module file` / `error #6580: Name in only-list does not exist`**
A source file does `use some_mod, only: SomeRoutine`, but `some_mod`'s
`.F`/`.F90` file isn't listed in `COMMON_SOURCES` (or `SMACOS_ONLY_SOURCES`/
`MACOS_ONLY_SOURCES`) in the main `CMakeLists.txt`, so its `.mod` file was
never built. This typically happens when a newer branch adds a new Fortran
module that the `CMakeLists.txt` hasn't been updated to include yet. Find
the missing file and add it to `COMMON_SOURCES`:
```cmd
dir /s /b macos_f90\\*MODULENAME*
```
(replace `MODULENAME` with the actual module name). Then add
`${smacos_src_dir}/MODULENAME.F` to the `COMMON_SOURCES` list, anywhere
before the targets that consume it (CMake's Fortran dependency scanner
generally orders compilation correctly once a file is present in the
sources list — exact position rarely matters). If you're tracking a branch
that adds new modules over time, periodically diff the contents of the
`macos_f90` folder's `.F` and `.F90` files against what's listed in
`COMMON_SOURCES` to catch additions before they surface as build errors
one at a time.

**`mmacos_mex.F not found` during configure**
Set `-DMMACOS_SRC_DIR` to the full path of the `mmacos/` directory:
```cmd
cmake -B build -DMMACOS_SRC_DIR=C:\path\to\MACOS_resources\mmacos
```
Newer `mmacos` checkouts moved Fortran sources into `mmacos/src/` (mirroring
pymacos's `src/pymacos` layout) instead of keeping them directly in `mmacos/`.
Both this project's main `CMakeLists.txt` and `mmacos/CMakeLists.txt` check
both locations automatically (`mmacos/src/mmacos_mex.F` first, falling back
to `mmacos/mmacos_mex.F`), so this should resolve itself on a current
checkout. If you still hit this error, confirm `mmacos_mex.F` actually
exists at one of those two locations in your checkout.

**MATLAB closes silently with no crash dump when calling `mmacos`/`macos.Session`, even though the MEX loads without a "module not found" error**
This is a different failure mode than the DLL-not-found case above — no
error dialog, no crash dump in `%LOCALAPPDATA%\Temp\1\matlab_crash_dump*`,
MATLAB just exits. Two distinct causes produce this identical symptom;
check both:

1. **Mixed Intel Fortran runtime libraries** (Windows only). Run:
   ```cmd
   dumpbin /dependents build\bin\mmacos.mexw64 | findstr /i ifcore
   ```
   If this shows **both** `libifcoremd.dll` and `libifcorert.dll`, the MEX
   is linking two incompatible Intel runtime variants — one multi-threaded
   DLL runtime, one single-threaded static runtime — and loading both in
   the same process causes a silent exit during runtime initialization,
   before any of your code runs. This is caused by a stray
   `/reentrancy:none` compile flag on the **Windows** branch of
   `mmacos/CMakeLists.txt`'s Fortran compile options. `/reentrancy:none`
   is correct and required on **Linux** (see the `mmacos` Makefile's own
   `LDFLAGS`, which uses `-reentrancy=none` deliberately to avoid a
   different SIGSEGV-at-exit issue) but must **not** appear in the Windows
   `target_compile_options` block, since `smacos.lib`/`fitslib.lib` (built
   by the main `CMakeLists.txt`) link the standard multi-threaded DLL
   runtime. Compare against `GMI.mexw64`'s dependencies (`dumpbin
   /dependents build\bin\GMI.mexw64`) — a working GMI build shows only
   `libifcoremd.dll`; use that as your reference for what "clean" looks
   like.
2. **`macos_param.txt` missing** — see the entry above. Confirm with
   `exist('path/to/macos_param.txt', 'file')` in MATLAB before assuming
   the runtime-mixing cause applies; this file's absence produces the
   exact same symptom (a `stop` statement inside `macos_io_failure`/
   `param_mod_init` that exits MATLAB with no dialog and no dump).

Diagnostic tip: a genuine access-violation crash always writes a dump file
to `%LOCALAPPDATA%\Temp\1\matlab_crash_dump.<pid>-1`. If MATLAB closes and
no *new* dump appears (check the timestamp), the cause is one of the two
silent-exit paths above, not a memory-access crash — don't spend time on
crash-dump frame analysis for this symptom.

**Newer `mmacos` checkout: `macos.Session` / package functions not found**
`macos.Session(...)` is a MATLAB package-qualified call resolved via a
`+macos/Session.m` file. For MATLAB to find it, the **parent** directory
of the `+macos` folder must be on the MATLAB path — typically
`mmacos/src/`, not `mmacos/` itself and not the CMake `build/bin/` output
folder. You need **both** on the path: `build/bin` (for the compiled
`mmacos.mexw64`) and `mmacos/src` (for the `+macos` package `.m` files):
```matlab
addpath('path/to/pymacos/src/macos/build/bin')
addpath('path/to/MACOS_resources/mmacos/src')
rehash
m = macos.Session(256);
```

**`which <function> -all` lists a stale MEX file ahead of a freshly-built one**
MATLAB's `which -all` ordering depends on path order and current directory,
not file modification time — an old build artifact left in an `examples/`
or test folder can shadow a fresh rebuild if its containing folder happens
to resolve first. Use `addpath(..., '-begin')` to force your build output
to the front of the search path, or simply delete/rename stale `.mexw64`
copies you find via `which <name> -all` so they can't be picked up by
accident in a later session.
