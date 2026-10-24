### `design_layer_api`: Design-layer API examples (`+macos/+design`)

Runnable worked examples for the MATLAB design layer -- the four
scripts that introduce the two front-ends.  The bench / coronagraph /
interferometer rows below moved to their own thread directories in the
2026-08 reorganization (see `../../00_INDEX.md`);
their descriptions are kept here until the Week-2 rewalk absorbs this
file into that index.  Two front-ends,
one analysis core (PLAN_DESIGN_LAYER §1.0): **import** an existing
prescription (`macos.design.System`, Sprint 2A-i) or **build** one from
design intent (`macos.design.Telescope`, Sprint 2A-ii).  Both feed the
same `vary → evaluate → optimize` surface and the bitwise-verified
Phase 7 sensitivity machinery.  Each script ends in `exit(0)` (batch-mode
rule, `../../../CLAUDE.md`).

The builder examples derive geometry from closed-form optics
(`../../../../optical_design/` reference + fixtures); the import examples run
against the committed `e5hex1.in` in `../../50_sensitivities/e5hex1/` and run
unchanged on a CodeV/Zemax-converted prescription (the expected dominant
entry point, §1.0).

Run one:

```bash
matlab -batch "run('$(pwd)/example_telescope_design.m')"
```

These are distinct from `../../50_sensitivities/`
(lower-level `macos.dw_dx` / `macos.dw_dz_zernike` driver demos); the
design-layer examples wrap that machinery behind the `+macos/+design`
API.  The telescope-design *progression* (rc_* -> tma_* -> freeform) is
the rest of `../`.

**Run it.**

```matlab
cd('<path-to>/mmacos/templates/10_telescopes/design_layer_api');
run(fullfile(pwd,'example_telescope_design.m'));
```

The README's design-template example, `example_telescope_design.m`, the builder front end from design intent to a saved deliverable. Its three siblings are the import front end (`example_sensitivities_from_rx.m`, `example_align_from_rx.m`) and a focused alignment loop.

**What it prints:**

<!-- transcript: p2_design_layer_api_run -->

**Its figure:** none: the example prints its stages and saves the deliverable (`rc_6m_f20.in`, `.mat`).

**The number it reproduces.** From design intent (6 m, f/20) the builder derives a Ritchey-Chretien with K1 = %pin:p2_design_layer_api_k1% and K2 = %pin:p2_design_layer_api_k2%, diffraction-limited on axis (its rms wavefront error is at the floor of the arithmetic, printed as a bound). A 2 mm despace of the secondary costs %pin:p2_design_layer_api_pert% m rms, and the optimizer returns it to nominal, again to the floor.

**Adapt it.** *"Using templates/10_telescopes/design_layer_api/example_telescope_design.m, build a 4 m f/15 Ritchey-Chretien instead, despace the secondary by 1 mm and let the optimizer recover it, and report the conics and the recovered wavefront error."*

