# Mac release gate (4) — gfortran, Apple Silicon

Pair: macos `36699d3` + MACOS_resources `629a2d2`. The macos tip advanced past
pair 3's `725e0389`: ImageDemo repaired, the SRS near-field sign fix (gg), the
ORS admission + RSTRACE fixes (ii), the manual rebuilt with a `transcripts`
gate. **Resources `629a2d2` is unchanged from pair 3**, so the binding-side
(mmacos/pymacos) residual is identical to pair 3 by construction. Host: macOS
26.6, gfortran 16.1.0 + Apple clang 21, MATLAB R2024a arm64. Run 2026-10-19,
clean clone. Baselines: Linux `bba2e362` (load 160/160, cmd 103/103).

## Counts (check — Mac — Linux — differs?)
- 1 clean clone — both at SHA — — — no
- 2 build release+debug — OK, arm64, 189 syms — OK — ifx n/a
- 3a load (161) — 161/161 ok, rt 154/7, 1 ulp — 160/160 — 2 edge-ray only (ImageDemo now MATCHES Linux)
- 3b cmd — 96/103 pass, 2 must-fail red — 103/103 — 7 fail (same as pair 3)
- 4 mmacos fast — 650 pass / 7 fail, subprocess filtered — all pass — 7 (same as pair 3)
- 5 pymacos — 6939 pass / 2 fail / 2 skipped(darwin) — pass — 2 (same as pair 3)
- 6 /png + /xw — render, glyph `1x10-10` — pass — none
- 7 self-containment — engine CLI + pymacos import with ~/dev renamed — — — none

## What this pair verified
- **The gg/ii engine fixes work on arm64:** the new ORS/RSTRACE/SRS and BEAM journals (the cmd corpus grew 102→103) all PASS on Mac; load is clean (161/161) with only the 2 known edge-ray decks moving; ImageDemo now matches Linux (578/163).
- **COORD (ff) persists** — the engine delta did not touch `LocalCoord`, so `coord 6` on CassWithExitPupil.in still prints the 76.5°-rotated frame on the -O2 release build (xLocal (-1,0,1.98e-11), yLocal (1.9e-11,-0.2335,0.9724), zLocal (4.6e-12,0.9724,0.2335)); correct under -O0 debug. Still TO's -O2 LocalCoord item.
- **(ee) guards hold** — the 3 mmacos subprocess tests Filtered-by-assumption, the 2 pymacos capfd tests skipped(darwin), with reason.

## Residual — unchanged from pair 3 (resources 629a2d2 identical); round-3 for CC
Exact Mac values are in `mac_gate_gfortran_d3f2777.md`. In brief:
- cmd 7: COORD_ACOOR (-O2 bug), EXPORT (BSD `wc -l` leading space), CALIB_SPOT_mixed (tilt E-05, journal pattern vs binding), GPERTURB/WLENS_SWL/ZABERR/OBS_POLYGON (near-zero/precision pins; perl + fixture fixes all work).
- mmacos 7 / pymacos 2: mixed_SPOT rms 2.7255e-9 (>1e-10), tDwDx near-zero element 1.42e-6 vs -2.49e-6 (>AbsTol 1e-6), tPolElement a ray >AbsTol 1e-15, tOpdRef chief-path count 0, tCalibStop stop-reason/n_iter, spot_centre brittle allclose (6305 vs -5938) still present, trial-loses-chief empty on arm64.

## Not run / notes
- ifx: no Apple-Silicon arm — gfortran only.
- Raw logs /tmp/mac_gate4_*.log; cmd transcripts cli_tests/records/cmd_gfortran_36699d3_transcripts/.
- No engine *defect* except (ff) the -O2 LocalCoord miscompile; the rest is arm64 round-off, the (ee) unit-6 mechanism (api-fixable), one BSD `wc`, and round-3 re-pins/sub-asserts.
