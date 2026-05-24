# MACOS — Modeling and Analysis for Controlled Optical Systems

Jet Propulsion Laboratory

---

## Prerequisites

- **Intel oneAPI** (`ifx` compiler) — install from intel.com/oneapi.
  The build scripts source `/opt/intel/oneapi/setvars.sh` automatically.
- **gcc / g++** — for C sources and the Giza graphics library.
- **CMake ≥ 3.20** — verify with `cmake --version`.
- **X11 development headers** — `sudo apt install libx11-dev` on Ubuntu/Debian.
- **MATLAB** *(optional, for `GMI.mexa64`)* — any recent release under `/usr/local/MATLAB/`.

For a gfortran build use `makegfortran.sh` instead; Intel oneAPI is not required.

---

## 1. Directory layout

Both repositories must be siblings under the same parent directory.
Replace `<username>` with your login name throughout.

```
/home/<username>/dev/
├── macos/            ← main MACOS/SMACOS source tree
└── MACOS_resources/  ← GMI mex interface and SegMirMaker
```

---

## 2. Clone the repositories

```bash
mkdir -p ~/dev && cd ~/dev

git clone git@github.com:nasa-jpl/macos.git
git -C macos checkout joint-dev

git clone git@github.com:nasa-jpl/MACOS_resources.git
git -C MACOS_resources checkout dr-dev
```

---

## 3. Build

All scripts run from `~/dev/macos/` and accept options
`[debug|release] [gfortran]` in any order
(defaults: release, ifx).

> **PGPLOT removed on release-candidate.**  Giza is now the only
> PGPLOT-API provider.  The `[giza|pgplot]` argument the build
> scripts previously accepted is gone.

### Build all four targets at once

```bash
cd ~/dev/macos
source ./makeall.sh              # ifx, release (most common)
source ./makeall.sh debug        # ifx, debug
source ./makegfortran.sh         # gfortran, release
```

### Build individual targets

| Script | Targets built |
|---|---|
| `source ./makems.sh` | `macos` + `libsmacos.a` |
| `source ./makesd.sh` | `smacos_dvr` (also builds macos/smacos if needed) |
| `source ./makegmi.sh` | `GMI.mexa64` (requires macos+smacos built first) |

All three accept the same `[debug|release] [gfortran]` options.
`makesd.sh` and `makegmi.sh` must use matching options to the preceding
`makems.sh` (or `makeall.sh`) call so they target the same build directory.

---

## 4. Build outputs

After a default `source ./makeall.sh` (ifx, release):

| Target | Location |
|---|---|
| `macos` executable | `~/dev/macos/build_release/bin/macos` |
| `smacos_dvr` executable | `~/dev/macos/build_release/bin/smacos_dvr` |
| `libsmacos.a` static library | `~/dev/macos/build_release/lib/libsmacos.a` |
| `GMI.mexa64` MATLAB mex | `~/dev/MACOS_resources/GMI/GMI.mexa64` |

Build directory names for other combinations:

| Script + options | Build directory |
|---|---|
| `makeall.sh` | `build_release/` |
| `makeall.sh debug` | `build_debug/` |
| `makegfortran.sh` | `build_release_gfortran/` |
| `makegfortran.sh debug` | `build_debug_gfortran/` |

---

## 5. Shell aliases

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

---

## 6. Rebuilding after source changes

Re-run the same script; CMake recompiles only changed files.
To force a clean rebuild delete the build directory first:

```bash
rm -rf ~/dev/macos/build_release
source ~/dev/macos/makeall.sh
```
