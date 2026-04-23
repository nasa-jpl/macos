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

Run one script from the `macos/` directory — it builds all four targets.

### Intel ifx (recommended)

```bash
cd ~/dev/macos
source ./makejoint.sh
```

### gfortran alternative

```bash
cd ~/dev/macos
source ./makegfortran.sh
```

Both scripts accept optional `debug` and `pgplot` arguments (order-independent):

| Invocation | Compiler | Graphics | Optimisation |
|---|---|---|---|
| `source ./makejoint.sh` | ifx | Giza | -O2 (default) |
| `source ./makejoint.sh debug` | ifx | Giza | -O0 -check all |
| `source ./makejoint.sh pgplot` | ifx | PGPLOT | -O2 |
| `source ./makegfortran.sh` | gfortran | Giza | -O2 |
| `source ./makegfortran.sh debug` | gfortran | Giza | -O0 |

---

## 4. Build outputs

After a default `source ./makejoint.sh` (ifx, Giza, release):

| Target | Location |
|---|---|
| `macos` executable | `~/dev/macos/build_release_giza/bin/macos` |
| `smacos_dvr` executable | `~/dev/macos/build_release_giza/bin/smacos_dvr` |
| `libsmacos.a` static library | `~/dev/macos/build_release_giza/lib/libsmacos.a` |
| `GMI.mexa64` MATLAB mex | `~/dev/MACOS_resources/GMI/GMI.mexa64` |

Build directory names for other combinations:

| Script + options | Build directory |
|---|---|
| `makejoint.sh` | `build_release_giza/` |
| `makejoint.sh debug` | `build_debug_giza/` |
| `makejoint.sh pgplot` | `build_release_pgplot/` |
| `makegfortran.sh` | `build_release_giza_gfortran/` |
| `makegfortran.sh debug` | `build_debug_giza_gfortran/` |

---

## 5. Shell aliases

Add to `~/.bashrc` (or `~/.bash_aliases`), then `source ~/.bashrc`:

```bash
alias macos='~/dev/macos/build_release_giza/bin/macos'
alias smacos_dvr='~/dev/macos/build_release_giza/bin/smacos_dvr'
```

For a debug build:

```bash
alias macos='~/dev/macos/build_debug_giza/bin/macos'
alias smacos_dvr='~/dev/macos/build_debug_giza/bin/smacos_dvr'
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
rm -rf ~/dev/macos/build_release_giza
source ~/dev/macos/makejoint.sh
```
