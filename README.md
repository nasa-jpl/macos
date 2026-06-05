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
- [Troubleshooting](#troubleshooting)

---

## Directory Layout

Both repositories must be siblings under the same parent directory.
Replace `<username>` with your login name throughout.

```
/home/<username>/dev/
├── macos/            ← main MACOS/SMACOS source tree
└── MACOS_resources/  ← GMI mex interface and SegMirMaker
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
| MATLAB *(optional)* | Required only for `GMI.mexa64` |

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

#### Build everything (macos + smacos)
```bash
rm -rf build_release
cmake -B build_release -G "Ninja" -DCMAKE_BUILD_TYPE=Release
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

| Tool | Version | Purpose |
|------|---------|---------|
| Intel oneAPI HPC Toolkit | 2026.0 or later | Fortran (`ifx`) and C (`icx`) compilers |
| Visual Studio | 2022 or 2026 Community | Linker and Windows SDK |
| CMake | 4.1.1 confirmed working | Build system |

> **Note:** CMake 4.2.x and 4.3.x have a known incompatibility with ifx 2026.x.
> CMake 4.1.1 + ifx 2026.0 is the confirmed working combination.

### Setup

1. Install [Visual Studio 2022 Community](https://visualstudio.microsoft.com/vs/older-downloads/)
   or Visual Studio 2026, selecting **"Desktop development with C++"**

2. Install [Intel oneAPI HPC Toolkit 2026.0](https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit.html)

3. Install [CMake](https://cmake.org/download/) — select
   **"Add CMake to system PATH"** during install

4. Open the **Intel oneAPI Command Prompt** from the Start menu
   (this sets up `ifx` and `icx` on your PATH automatically)

### Building with the Batch Script

Run from the **Intel oneAPI Command Prompt**:

```cmd
makems.bat              REM release (default)
makems.bat debug        REM debug
makems.bat npsol        REM release + NPSOL
makems.bat debug npsol  REM debug + NPSOL
```

### Building with CMake Directly

#### Build everything (macos + smacos)
```cmd
rmdir /s /q build
cmake -B build -G "Ninja" -DCMAKE_BUILD_TYPE=Release
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

> **Note:** GMI and smacos_dvr are not supported on Windows.

---

## Build Options Reference

Pass these to cmake with `-D<OPTION>=ON/OFF`:

| Option | Default | Description |
|--------|---------|-------------|
| `BUILD_MACOS` | ON | Build the interactive `macos` executable |
| `BUILD_SMACOS` | ON | Build the `libsmacos.a` static library |
| `BUILD_SMACOS_LEAN` | OFF | Minimal smacos (no FITS, no graphics). Mutually exclusive with `BUILD_SMACOS` |
| `BUILD_SMACOS_DVR` | OFF | Build the `smacos_dvr` test driver (Linux only) |
| `BUILD_GMI` | OFF | Build the MATLAB GMI mex interface (Linux only) |
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
| `~/dev/MACOS_resources/GMI/GMI.mexa64` | Linux | MATLAB mex interface |

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

## Troubleshooting

**`CMAKE_Fortran_PREPROCESS_SOURCE` error on Windows**
Caused by CMake 4.2.x + ifx 2025.3.x incompatibility. Update Intel oneAPI
to 2026.0 or later. CMake 4.1.1 + ifx 2026.0 is confirmed working.

**`mmintrin.h` errors when building Cairo/Giza on Windows**
Caused by icx version incompatibility with the Windows SDK. Update Intel
oneAPI to 2026.0 or later and ensure Visual Studio 2022 or 2026 is installed.

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

**MATLAB root not found (GMI build)**
Set the MATLAB installation path:
```bash
cmake -B build_release -DMATLAB_ROOT=/usr/local/MATLAB/R2025b \
      -DBUILD_GMI=ON -DCMAKE_BUILD_TYPE=Release
```

**`LNK2019` unresolved external symbols on Windows**
Ensure you are using the latest `CMakeLists.txt` which includes
`windows_stubs.c` generated at configure time to provide Windows-compatible
stubs for Linux-only functions.
