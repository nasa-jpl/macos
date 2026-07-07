# MACOS Command Reference

*One optical engine, four ways to drive it.*

This document is the command-level reference for the MACOS optical
modeling system and its programming interfaces.  It complements the
MACOS Manual (concepts, modeling guidance, worked examples) with a
catalog: every command and function, one at a time, with its syntax.

## The four surfaces

**MACOS** is a Fortran optical analysis engine: exact ray tracing,
differential ray tracing (linear optical models), and diffraction
propagation, built for *system-level* modeling — optics on dynamic
structures, segmented and deformable mirrors, controlled optics —
rather than lens design.  There are four ways to use it:

1. **The `macos` CLI** — an interactive command interpreter.  You type
   commands (about 140 of them; Part I of this reference) at the
   `MACOS>` prompt; the engine prompts for whatever arguments it
   needs.  Sessions can be recorded to journal (`.jou`) files and
   replayed with `EXEcute`.  This is the exploratory surface and the
   place to learn what the engine can do.

2. **SMACOS** — the same engine as a Fortran subroutine
   (`libsmacos.a`).  Your program issues the *same command strings*
   programmatically and passes arguments through a fixed buffer list.
   This is the integration surface for custom Fortran drivers and is
   what every higher-level binding is built on.  Same command catalog
   (Part I); calling convention in the SMACOS chapter.

3. **mmacos** — a MATLAB package (`+macos`) wrapping the engine
   through a single MEX gateway.  Functions like `macos.load_rx`,
   `macos.trace`, `macos.opd`, `macos.perturb`, `macos.dw_dx` give
   MATLAB-native access (matrices in, matrices out) plus higher-level
   drivers the CLI does not have: multi-channel Jacobians, the CALIB
   optimizer interface, segment/Zernike grid-basis generators, and the
   parametric telescope **design layer** (`macos.design.*`).
   Parts II and III.

4. **pymacos** — a Python package wrapping the same engine via f2py,
   with numpy arrays and (mostly) the same function vocabulary as
   mmacos.  Part II documents both side by side; where only one
   language has a function, the entry says so.

All four share one prescription format (`.in` files), one state model,
and one set of units conventions, because underneath they are the same
compiled engine.  The bindings (3, 4) go through a common wrapper
layer, `macos_f90/macos_api_mod.F90` — each Part II entry names the
underlying routine, which is the ground truth for argument meanings.

## Where things live

| Piece | Location | Build |
|---|---|---|
| engine + CLI (`macos`) | `nasa-jpl/macos` repo, `macos_f90/` | `source ./makeall.sh` (cmake; gfortran or Intel ifx) |
| SMACOS library | same repo | built by `makeall.sh`/`makems.sh` (`libsmacos.a`) |
| mmacos | `MACOS_resources/mmacos/` | see its README (MEX build against libsmacos) |
| pymacos | `MACOS_resources/pymacos/` | see its README (f2py build) |

## Prerequisite tags used throughout

- **\[Rx\]** — needs a loaded prescription (`OLD` or `NEW` first).
- **\[BLD\]** — needs a built linear model (`BUild` first).
- **\[DIFF\]** — needs a propagated wavefront (`PROpagate` first).

## Status of this document

Part I and Part II skeletons are **generated from the engine and the
binding sources** (`tools/gen_cmdref.py`), so the catalog cannot drift
from the code; entries marked *TODO* are awaiting expanded prose.
Part III (design layer, channels) documents packages that are under
active development and may change faster than the rest.
