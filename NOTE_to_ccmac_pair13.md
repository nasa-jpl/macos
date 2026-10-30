# NOTE to CCMac — your next pair: macos c6c83d55 + resources 446bd45 (CC, 2026-10-24)

On origin (macos tip 7ddf1d32 carries c6c83d55's records above it; resources unchanged since your pair 12).  Engine
CHANGED by ONE line -- the full 7-check clean-clone gate, with one place to look.  Dave is away until ~10-26; CC
holds his standing push word for gated work.

**What changed (CC, PLAN_CONSOLIDATION 1.2 (ooo), found by macos-fd on the tg_psi_dm templates):** POLARIZATION is
a session mode and has survived a load since 2026-10-11 (PLAN 4 item 2), but `src_mod_init_vars` zeroed the source
Jones vector `Ex0`/`Ey0` on every load, so polarization-then-load traced a ZERO field -- INTENSITY summed to exactly
0, `macos.intensity` all zeros, the templates' detector centroid 0/0.  `reinitialise_variables` (macos_init.F) now
keeps `Ex0`/`Ey0` across the load when `ifPol` is in force; a deck's own `Ex0Ey0=` line still overrides.  Journal
`POL_across_load.jou` (new): the switch before and after the load read the same intensity sum on the manual's
Cassegrain (6.4607758D+10 on Linux, a signal -- pinned to five figures by the rule); red pre-fix.

**Linux (both compilers, clean worktree):** load 160/160, ab 08b92085 -> c6c83d55 0 movers, cmd 120/120 + 2
must-fail.  pymacos side tree: tg_test.in 0.185955 in both orders.  On the relinked mmacos mex the two tg_psi_dm
templates pass every gate (v1: 7.4786 deg arm rotation, closure corr 0.998434, 0.304 nm; v2: 0.999573, 0.183 nm).

**Expected on arm64:** cmd 120/120 + 2 must-fail (119 + POL_across_load); load 165/165, 0 movers vs your gate-12
base; mmacos fast 711/0; pymacos 7023/0/1 skip.  **The one place to look:** POL_across_load's intensity sum is a
five-figure pin on a diffraction propagation (65536 grid points) -- if arm64 differs past the fifth figure, send the
two sums and it becomes a bound; and after your gate, if you have a MATLAB window, run the two tg_psi_dm
examples (templates/90_polarization) once: they were refusal pages until this fix, and their records above are
the first the Mac has had a chance at.

**Also on origin, no engine effect:** the mmacos manual's Part II design catalog (docs/mmacos-manual, 109 pages,
35 template pages generated from live runs) and the findings (nnn)-(rrr) in PLAN 1.2 for Dave.
