# HOW TO COMPILE — MACOS / SMACOS (Linux)

A short, copy-paste walkthrough that takes you from **nothing** to a
**running `macos`** in three steps. Everything compiles from one command;
you do **not** need to learn CMake.

> Busy? The whole thing is:
> ```bash
> mkdir -p ~/dev && cd ~/dev                                  # 1. make a parent dir
> git clone git@github.com:nasa-jpl/macos.git                 #    pull both repos
> git clone git@github.com:nasa-jpl/MACOS_resources.git       #    (side by side)
> git -C macos checkout opt-dev && git -C MACOS_resources checkout opt-dev
> cd macos && source ./makeall.sh                             # 2. build everything
> ./build_release/bin/macos                                   # 3. run it
> ```
> The rest of this file just explains each line. For the full option matrix
> (Windows, CMake-direct, lean/NPSOL/pymacos builds) see [README.md](README.md).

---

## A fresh clone is pure source — the first build makes *everything*

There are **no pre-compiled libraries** in the repository. The first
`source ./makeall.sh` builds the **whole dependency stack** from source, in
order, before it builds MACOS itself:

| Built on first run | What it is | How |
|--------------------|-----------|-----|
| **readline** | arrow-key history at the `MACOS>` prompt | autoconf `./configure && make` — the script runs it for you, **once** |
| **zlib + libpng + pixman + Cairo** | the 2-D graphics back end (all **bundled** under `macos_f90/giza/src/subprojects/`) | compiled by CMake |
| **Giza** | the PGPLOT-API graphics layer (drop-in PGPLOT replacement) | compiled by CMake |
| **SLSQP** | the constrained-optimization solver (vendored, BSD) | compiled by CMake |
| **libsmacos.a → macos / smacos_dvr / GMI** | the engine, library, driver, MATLAB mex | compiled by CMake (GMI via its Makefile) |

The **only** graphics piece that comes from your system is **X11** — Cairo,
libpng, zlib and pixman are all in the repo. That is why the prerequisite
list below is so short.

> **So the first build is slow** (several minutes — it's compiling Cairo and
> readline too). Every build *after* that is incremental and fast: CMake
> recompiles only the files you changed.

---

## What you end up with

`source ./makeall.sh` produces **all four** targets in one shot:

| Target | What it is | Where it lands |
|--------|-----------|----------------|
| `macos` | the interactive ray-tracer | `build_release/bin/macos` |
| `libsmacos.a` | the static library (backs pymacos / mmacos / GMI) | `build_release/lib/libsmacos.a` |
| `smacos_dvr` | the batch/scripting test driver | `build_release/bin/smacos_dvr` |
| `GMI.mexa64` | the MATLAB interface (skipped if no MATLAB) | `~/dev/MACOS_resources/GMI/GMI.mexa64` |

---

## Step 0 — Install the prerequisites (one time)

Pick **one** Fortran compiler:

- **Intel oneAPI HPC Toolkit** (`ifx`) — the default, fastest binary.
  Free download: [intel.com/oneapi](https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit.html)
- **gfortran** — no Intel install needed. Use `makegfortran.sh` instead of
  `makeall.sh` everywhere below. This is also the supported compiler for the
  MATLAB (GMI / mmacos) interface.

Plus the common tools (Ubuntu/Debian shown):

```bash
sudo apt install build-essential cmake gfortran git libx11-dev
```

- `cmake` ≥ 3.21, `make`, `gcc`/`g++` — the C sources and the bundled graphics
  stack need them. **The C compiler must be `gcc`** (the build scripts set this
  for you).
- `libx11-dev` — the interactive graphics window. This is the **only**
  graphics dev package you need; Cairo/libpng/zlib/pixman are bundled and built
  from source, so do **not** apt-install them.
- MATLAB is **optional** — needed only for `GMI.mexa64` (auto-detected under
  `/usr/local/MATLAB/`; the build just skips GMI if it's missing).

You do **not** need to build readline, set `LD_LIBRARY_PATH`, or source the
Intel environment by hand — the build scripts do all of that.

---

## Step 1 — Pull the code (two repos, side by side)

MACOS is **two repositories** that must sit **side by side under one parent
directory**, both checked out to the **same branch name**:

```
~/dev/
├── macos/            ← engine + build scripts (this repo)
└── MACOS_resources/  ← GMI / mmacos / pymacos / SegMirMaker
```

### First time (fresh clone)

```bash
mkdir -p ~/dev && cd ~/dev

git clone git@github.com:nasa-jpl/macos.git
git clone git@github.com:nasa-jpl/MACOS_resources.git

# Check BOTH out to the SAME branch.  opt-dev = current release target.
git -C macos            checkout opt-dev
git -C MACOS_resources  checkout opt-dev
```

### Already cloned (just get the latest)

```bash
git -C ~/dev/macos            pull
git -C ~/dev/MACOS_resources  pull
```

> **The two repos must be on the same branch** (`opt-dev` with `opt-dev`,
> `sls-dev` with `sls-dev`). The bindings compile against the engine's Fortran
> module files; a mismatched pair points at the wrong modules and the API
> looks "missing" (e.g. unresolved `macos_api_mod_mp_*` symbols at link, or
> `init` absent from the module layer). If a build fails right after a `git
> pull`, check both branches match first:
> ```bash
> git -C ~/dev/macos branch --show-current
> git -C ~/dev/MACOS_resources branch --show-current
> ```

---

## Step 2 — Build everything (one command)

From inside `~/dev/macos`:

```bash
cd ~/dev/macos
source ./makeall.sh
```

That one command:
1. sources the Intel oneAPI environment for you (silently),
2. bootstraps the bundled readline library the **first** time only,
3. configures + compiles the bundled graphics stack, SLSQP, and then macos,
   smacos, smacos_dvr, and GMI.

The first run takes several minutes (it's building Cairo + readline too); it
finishes with a `BUILD SUCCEEDED` block listing each output path.

> **Use `source`, not `./makeall.sh`.** The leading `source ` (or `.`) runs the
> script in *your current shell* so its environment setup sticks. Running it as
> `./makeall.sh` spawns a throwaway subshell and the Intel/library setup is
> lost. This is the single most common "it won't build / won't run" mistake.

### No Intel compiler? Use gfortran (one word changes)

```bash
cd ~/dev/macos
source ./makegfortran.sh
```

Same four targets, output in `build_release_gfortran/` instead of
`build_release/`. Nothing else differs.

---

## Step 3 — Run it

The binary is self-contained (readline and the graphics stack are baked in;
Intel/gfortran runtime libraries are found automatically), so just run it:

```bash
~/dev/macos/build_release/bin/macos             # ifx build
# or, for the gfortran build:
~/dev/macos/build_release_gfortran/bin/macos
```

Make it a one-word command by adding an alias:

```bash
echo "alias macos='~/dev/macos/build_release/bin/macos'" >> ~/.bashrc
source ~/.bashrc
# now just type:  macos
```

For the MATLAB interface, add GMI to your MATLAB path (in `startup.m` or via
**Set Path**):

```matlab
addpath('~/dev/MACOS_resources/GMI')
```

---

## Which script do I run?

`makeall.sh` is almost always what you want. The others are shortcuts for
when you only need part of the tree (faster rebuilds):

| Script | Builds | Use it when… |
|--------|--------|--------------|
| `source ./makeall.sh` | macos + smacos + smacos_dvr + GMI | **default** — first build, or you want everything |
| `source ./makems.sh` | macos + smacos | you only touched engine source (fastest) |
| `source ./makesd.sh` | smacos_dvr | you only need the batch driver |
| `source ./makegmi.sh` | GMI.mexa64 | you only need to relink the MATLAB mex |
| `source ./makegfortran.sh` | all four, via gfortran | you don't have Intel oneAPI (or are debugging an ifx-only issue) |

Every script takes the same options, in any order:

| Option | Effect | Build directory |
|--------|--------|-----------------|
| *(none)* | `ifx`, optimized | `build_release/` |
| `debug` | `-O0 -check all` | `build_debug/` |
| `gfortran` | use gfortran (`makems`/`makesd` only) | `…_gfortran/` |
| `npsol` | also build the opt-in NPSOL solver (SLSQP is always on) | `…_npsol/` |

Examples: `source ./makems.sh debug` · `source ./makems.sh release gfortran`

> `makesd.sh` and `makegmi.sh` read the library built by the **preceding**
> `makeall`/`makems` call, so give them the **same** options (compiler +
> debug/release) so they look in the same build directory.

---

## Rebuilding after you edit source

Just re-run the **same** script — CMake recompiles only the files that
changed (pure `.m` MATLAB-layer edits need no rebuild at all):

```bash
source ./makems.sh          # engine edit → fast incremental rebuild
```

Force a clean rebuild by deleting the build directory first (this also forces
the bundled libraries to rebuild, so it's slow again):

```bash
rm -rf ~/dev/macos/build_release
source ~/dev/macos/makeall.sh
```

---

## If it doesn't work — the five usual causes

1. **`ifx: command not found`** — Intel oneAPI isn't installed or wasn't
   sourced. Either install it, or just build with gfortran:
   `source ./makegfortran.sh`.
2. **Ran it as `./makeall.sh`** — the environment didn't stick. Re-run with
   `source ./makeall.sh`.
3. **`macos_api_mod` symbols missing / unresolved `macos_api_mod_mp_*` at
   link** — the two repos are on **different branches**, or `MACOS_resources`
   isn't a sibling of `macos`. Put them side by side under `~/dev/` and check
   `git branch --show-current` matches in both.
4. **GMI / smacos_dvr says "module directory not found"** — build the engine
   first (`makeall.sh` or `makems.sh`), with the **same** compiler, so the
   `mod_smacos/` modules exist in the build directory it looks in.
5. **X11 / graphics build error** — install the dev headers:
   `sudo apt install libx11-dev`. (Cairo/png/zlib are bundled — don't go
   hunting for those packages.)

For more (Windows, CMake-direct, lean/NPSOL builds, pymacos, and the full
troubleshooting list) see [README.md](README.md).
