## APPENDIX D: Working with AI Agents

This appendix is written for an AI agent that is asked to run, extend or document MACOS, and for the person who sets one to work. It says where the authoritative facts are, how to drive the program so that its answers can be trusted, which conventions bite, and the rules this project holds its own agents to. Nothing here is specific to one agent product.

### Where the facts live

-   **The Command Reference** (`docs/macos-manual/cmdref/`, rendered as macosCmdRef) is the authority on every command: its minimum abbreviation, its syntax, its prompts and their defaults, its binding in MATLAB and Python, and the journal in the test suite that exercises it. When this manual's prose and a transcript disagree, the transcript made on the current build wins, and the Command Reference is updated from it.
-   **This manual** explains the models: the prescription language of Section 4 (Table 13 is the design-optimization keyword reference), the ray-trace, diffraction, imaging and linear-model commands of Sections 5 to 8, SMACOS in Section 9, worked examples in Section 10, polarization in Section 11. The example prescriptions of Appendix A are in `docs/macos-manual/examples/` and load as printed.
-   **The notes beside the code.** `CLAUDE.md` at the top of each repository, and in `macos_f90/`, `macos_f90/giza/`, `macos_f90/slsqp/` and `MACOS_resources/pymacos/`, record the engine's conventions and the defects that were closed, each with the test that pins it; `macos_f90/notes/` holds the engine record by topic. They are maintained with the code and are the first thing to read before changing it.
-   **The test suite is the specification.** `cli_tests/` runs the interactive program on every prescription in the repositories (the load gate) and one journal per command (the command gate), on both compilers; `MACOS_resources/mmacos/tests/` and `pymacos/tests/` are the binding suites. A behaviour that no test pins is not yet a promise.

### Driving the program

-   Run the interactive program on a pseudo-terminal with `TERM=dumb`; under a real terminal type the line editor turns on bracketed paste and prefixes every line with an escape sequence. Send one line per prompt. A prompt ends in `MACOS>`, in `]:` (a default in brackets) or in `?`; the next line answers it, and a line that happens to be a command name is read as the answer, not as a command.
-   The first prompt asks the model size (the largest ray grid the session will hold: 128, 256, 512, 1024, 2048). A prescription whose `nGridpts` exceeds it is clamped with a message; a diffraction grid is four times the ray grid, and the coronagraph example of Section 7 wants 512.
-   A journal (`cli_tests/cmd/*.jou`) is the reproducible form of a run: the lines sent, with `# EXPECT:` regular expressions the transcript must match and `# EXPECT-NOT:` ones it must not, written from a transcript on the current build, never from prose. `cli_tests/README.md` is the recipe, including the shell-escape and plot-device quirks.
-   Every number an agent reports comes from a run on the build it names, with the prescription and the commands that produced it. The engine prints its OPD reference (chief ray or aperture mean), its ray-loss counts by cause and element, and its propagation summaries; read them.
-   SAVE writes the prescription as loaded and modified; a saved prescription reloads to the same trace, and a second SAVE is byte-identical. Use the round trip to check that a change to a prescription means what was intended.

### Conventions that bite

-   `KrElt=` carries the sign of the vertex radius with a mirror's concave side facing the incoming beam negative; `KcElt=` is the conic constant; `psiElt=` points against the incoming beam on a reflector and along it on a focal plane or reference surface. The output-coordinate block (`nOutCord=`, `Tout=`) ends a prescription; when it is absent the program proposes the beam frame at the last element and asks.
-   A ray trace never changes the prescription; a trace reads it and writes trace state. The design-optimization keywords of Table 13 are prescription state, written by SAVE and announced at load.
-   Session options (POLARIZATION, VECTOR, CENTROID, LNEG, CIR, the OPD reference) survive loads; the bindings report them at every load.
-   The OPD map's zero is not "no ray": a valid ray at the chief ray's path reads exactly 0. Validity is the mask the engine keeps (`opd_mask` in the bindings), not a test on the value.
-   A numerical zero (a collimated system's RMS near 1e-12, an on-axis coordinate near 1e-11) has no shared digits across machines or architectures. Pin a real-signal value to five significant figures; bound a numerical zero in magnitude; never pin its digits.
-   Masks for diffraction are `Element= Obscuring`; a Reference surface with obscurations clips rays only. Reference and Return surfaces are passive: they are read, never traced through.

### Working rules

1.  Every fix ships with a test, and the test has a leg that fails on the build before the fix. A test that cannot fail measures nothing.
2.  Measure before you guard. Before adding a check, find the value's other writers; a refusal added at parse time once broke two hundred tests because the clamp already existed downstream.
3.  Guards warn and count; they do not refuse behaviour a long-standing tool tolerated. A new refusal needs the maintainer's decision.
4.  A change that can move a trace runs the whole corpus on the old and the new build before it is committed, and every prescription that moved is named and explained.
5.  Read the frame before calling a mismatch: angles and lengths depend on the coordinate frame they are quoted in.
6.  Look at the raw per-ray data before blaming a solver: a histogram of the per-ray OPD finds a wrong conic sheet that an outlier cut hides.
7.  An "exact" replica of a kernel that passes a test can hide a broken kernel; diff the engine's own intermediate quantity instead.
8.  Write the resolution of an oddity where the next reader will find it, at the time it is resolved.
9.  Report outcomes faithfully: if a test fails, say so with its output; if a step was skipped, say that; a number comes from a run, not from memory.
