# MACOS Command Reference

*One optical engine, four ways to drive it.*

This document is the command-level reference for the MACOS optical
modeling system and its programming interfaces: a catalog of every
command and function, one at a time, with its syntax.

## Companion document: the MACOS Manual

This reference is a lookup document; the **MACOS Manual** is the
learning document.  The Manual covers the concepts this catalog only
names — geometry and coordinate conventions (VptElt, PsiElt, the
global frame), prescription structure, ray-grid and source setup,
reference surfaces and diffraction modeling, linear optical models,
and SMACOS integration — and works complete examples end to end
(Appendix A pairs each example prescription with a journal file that
drives a full analysis).

The Manual lives beside this document in the repo
(`docs/macos-manual/`, built as `macosMan.pdf` — see the README there)
and is referenced throughout this catalog as "manual §N".  If you are
new to MACOS: read Manual Sections 1–4 first, run one Appendix A
example, and use this reference from then on.

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

## Where the examples are

Worked examples are not enumerated in this reference; they live with
each surface:

- **CLI / engine** — `docs/macos-manual/examples/`: the manual's
  Appendix A prescriptions, each `<name>.in` paired with a
  `<name>.jou` journal that drives a complete analysis
  (`EXEcute <name>` after loading).
- **mmacos** — `MACOS_resources/mmacos/examples/`; further runnable
  demos sit in `mmacos/design/` (telescope design layer) and
  `mmacos/sensitivities/` (Jacobian / dw_* drivers).
- **pymacos** — `MACOS_resources/pymacos/tests/`: the test suite
  doubles as usage examples.

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
