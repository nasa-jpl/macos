# BRIEF — the bindings program (TO): plan §1.4a (pymacos covers the api; one generator, three targets; one default per option) then §1.4b (the cmdref BINDING column; the missing CLI commands through the api)

CC for TO, 2026-10-16, on Dave's word (he chose this over the budget tool's follow-ons).  You
start from a cleared session; this file is the whole hand-off.  Read, in order:
`macos/CLAUDE.md` (root rules; the build summary points at `BUILD.md`), `macos/PLAN_CONSOLIDATION.md`
§1.4a and §1.4b (the program as planned, with Dave's rulings inline — read them before you touch
anything), `macos_f90/CLAUDE.md` (a 54-line index) and then `macos_f90/notes/api_bindings.md`,
`MACOS_resources/mmacos/CLAUDE.md` (the mex side) and `MACOS_resources/pymacos/CLAUDE.md` (the
f2py side; the venv, the build, the test layout), the generator
`mmacos/src/gen_mex_wrappers.py` (it reads `macos_api_mod.F90` as the single source of truth —
memory `reference_binding_codegen_source`), and `pymacos/src/cmake/source/pymacos_f2py.f90` +
`pymacos/src/pymacos/macos.py` (what a hand-written binding looks like today).  Both repos
`dev-candidate`, pushed through macos `8f1be77` / resources `920eb58`; state SHAs when you commit.
Engine changes are NOT expected; if an api routine needs a change to be bindable, say so first.

## A. §1.4a — pymacos covers the api (one generator, three targets)

A.1 **Re-measure the gap.**  List every public routine of `macos_api_mod.F90` from the source and
    every routine pymacos binds (f2py layer + veneer); the 2026-10-10 count was 136 of 146 with
    `deadband_notes_get`, `fwd_root_flips_get`, `draw_rays_cmd/_get/3d_get`, `elt_info_get`,
    `elt_obs_get`, `elt_z`, `elt_srf_zrn_*`, `jmat_elt_get`, `polelt_get/_set`, `ray_hist_set`,
    `ray_pos_hist_get`, `ray_status_get`, `src_seg_get`, `met_calc/_geom_get/_get`, `translate`,
    `xps_cmd` missing, plus whatever landed since (`modes_get`, `calib_stats_get`,
    `grid_idx_ovf_get`, …).  Commit the table as `cli_tests/records/api_binding_gap_<sha>.md`.
A.2 **One generator, three targets.**  Extend `gen_mex_wrappers.py` to emit, from the same parse
    of the api source, (1) the mex helpers it already emits, (2) the f2py layer, (3) the Python
    veneer stubs (docstring from the Fortran comment, argument marshalling by the same rules the
    mex uses).  Generated code is regenerated, never hand-edited; hand-written veneers that add
    behaviour stay in `macos.py` and call the generated layer.  Gate: the generator reproduces
    today's mex helpers byte-for-byte (the refactor-safety leg) before it emits anything new.
A.3 **The missing wrappers through it**, in the order of use: `ray_status_get` (the design layer's
    ray-loss guard), `modes_get`, `calib_stats_get`, the `elt_*` getters, `xps_cmd`, `met_*`,
    `draw_rays_*`, `polelt_*` / `jmat_elt_get`, the rest.  Each with a pytest that exercises it on
    a real deck (`meaningful tests`: no toy inputs that skip the code path).
A.4 **The coverage gate**: a pytest that lists the api routines FROM THE SOURCE and asserts each
    has a binding; must-fail = remove one binding and watch it go red.  Same gate for mmacos
    (`tApiCoverage`) if it does not exist.
A.5 **A side-tree pymacos build recipe** so a pre-fix `.so` red leg is routine (the gap recorded on
    `test_opd_ref`): write it in `pymacos/CLAUDE.md`, use it once in this brief.
A.6 **One default per option (Dave 2026-10-10).**  Every session/user option the api exposes —
    OPD reference, FEX/SXP axis, `spcOption` / spot centre, `ffcut`, the obscuration option, the
    five session modes, `remove_ptt` in the harvests, grid orientation, model size, `opd_mask` —
    has ONE default set in the engine at allocation; the three front ends inherit it.
    Deliverables: a cmdref table (option, engine default, how each front end reads / sets it, the
    override) generated from the api source by the same generator; `tDefaults` (mmacos) +
    `test_defaults.py` (pymacos) reading the same deck and asserting the same answers for each
    default and its override; the rule in `macos_f90/notes/api_bindings.md` (and one line in the
    index): a veneer never carries a default an api getter can supply.  Expect to find drift —
    the history is pymacos `vpt_center=True` vs the CLI's chief-centred spot, the api's chief FEX
    axis vs the CLI's centroid, the drivers' `'mean'` vs the manual's chief; each is a measured
    row in the commit, each fixed by removing the veneer's restatement, not by adding one.

## B. §1.4b — the cmdref BINDING column, then the missing commands through the api

B.1 The cmdref catalog (`docs/macos-manual/cmdref/`, generated from `macos_help.inc`) gains a
    BINDING column filled from the api source by the generator: `SMACOS command` / `direct
    setter` (name it) / `CLI-only` / `CLI-only, legacy` / `missing`.  **Dave's ruling:** the
    macos-internal LINEAR-MODEL family (LIN, LNEG, LOPD, LSPOT, LRESET, LPERTURB, LPREAD, LPI,
    LQW, LL, LLO, LOS, LS) and the differential commands (DOPD, DOPL, DVAR, DFOV, DADD, MVAR,
    AVAR) are `CLI-only, legacy` with the reason in the column (superseded by the `dw_dx` route;
    the analytic forms cover only some surfaces) — not binding targets.  BUILD / DMBUILD /
    PARTIALS the same.
B.2 The genuinely missing group, through the api, in order of use: **single-ray tracing first**
    (RAY / SRAY / RRAY / PRAY / SRT / SRTRACE / SEGRAYTRACE / SINT — the gauge benches want it),
    the Zernike / analysis set second (ZCOEF / ZRM / SZCO / MAP / ELTS; `zern_fit` is on
    `opdValid` since 2026-10-12), the image-plane post-processing third (MTF, BLUR, NOISE, …:
    veneer functions over the field the bindings already return, one implementation per
    binding; MULTISPEC / CADD through `compose`).  Each new api routine: added in
    `macos_api_mod.F90` following the buffer-wrapper template, generated into all three targets,
    gated in mmacos AND pymacos on a real deck.  New SMACOS-dispatched commands go in
    `macos_cmd_loop.inc` + `smacosutil.F`'s LoadStack, never `macos_ops.F`
    (`notes/pupil_fex_stop.md`, the SXP entry).

## C. How we work

One item in flight; its gate first, red before green (the side-tree `.so` / a side worktree for
the mex); fast suite + the pymacos suite after each landing (pymacos runs from `tests/`);
corpus A/B only if an engine file moves (then from the clean worktree `~/dev/macos_cli_base`);
commits local with SHA + branch to CC; CC gates; Dave says push.  Never relink the shared mex
while another MATLAB runs (`ps -C MATLAB`); one model-1024 MATLAB on the box.  Numbers only from
output you have read.  The generator's output is reviewed by diff, not by eye.

## D. Report

One message to CC per landed item: SHAs (both repos), the gate numbers (coverage count, the
byte-identical regeneration, the defaults table's drift rows), suite counts.  Anything needing
Dave (an api signature change, a default whose value is a judgment call) goes to him through
CC with the measurement attached.
