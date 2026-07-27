<!-- GENERATED FILE -- do not edit.
     Source: 00_frontmatter.md.in
     Numbers: generated/numbers.json (MATLAB driver)
     Regenerate: make polval-regen  (docs/macos-manual)
-->
# MACOS Polarization — Validation Report

**Scope.** This is the reviewer-facing companion to the polarization test
code: one evidence section per validation gate, each stating the claim, the
figure, the measured number, the truth it is measured against, and the test
that pins it in CI. It covers the work described in `PLAN_POLARIZATION.md`
Phases 0–2b (exposure of the engine's existing polarization physics, the
Jones-pupil layer and the polarization-aberration maps) and Phase 3a
Tranche 1 (vector propagation across a multi-leg chain).

It does **not** cover Phase 2c (coronagraph contrast floor), Phase 2d (the
interferometer track), Phase 3 (polarizer / waveplate / vector vortex),
Phase 3a Tranche 2 (Jones-through-chain), or Phase 4 (spatially variable
coatings). Those phases append their own evidence sections as they land —
see *Coverage and gaps* at the end.

## Provenance

| | |
|---|---|
| Generated | 2026-07-26 16:58:55 |
| Engine (`macos`) | `b19e7a6` on `pol-core` |
| Bindings (`MACOS_resources`) | `3c9a42b` on `pol-core` |
| MATLAB | R2026a |
| Model size | 128 |
| Host | red-river |
| Figure / number consistency | consistent (7 figures, all older than numbers.json) |

## How this document is produced

Every figure and every number below is generated. The prose sources
(`polval/*.md.in`) contain no numeric literals at all — only substitution
placeholders of the form `@@`*NAME*`@@`. One command refreshes the lot:

```
cd docs/macos-manual && make polval-regen   # measure + render
make polval                                 # build docx / HTML
make polval-pdf                             # and PDF
```

`polval-regen` runs the MATLAB driver
(`MACOS_resources/mmacos/tools/pol_validation_report/`), which re-runs every
validation case, writes `media/*.png` and `generated/numbers.json`, and then
substitutes those numbers into the prose. Rendering **fails** if any token is
unresolved, so a half-updated report cannot be built, and the provenance row
above stamps the exact engine and binding commits the numbers came from. A
stale figure beside a fresh number is not something a reader has to catch by
eye.

### Numbers this driver cannot measure

Some gates run in another language binding, under another compiler, in
another repository, or against a pre-fix engine binary that no longer exists
in the tree. They are **not** omitted, and they are **not** quietly presented
as if they were regenerated: each is marked *(external, captured DATE)* and
its producing command is recorded in `external.json` beside the driver.
Everything else on this page was measured by the run stamped above.

### Reading the numbers

Two conventions matter for interpreting round-off-level results:

* **Mean and variation are reported separately, never conflated.** A uniform
  diattenuation or retardance across the pupil is a state change — after
  folds it also absorbs the system's geometric rotation — and is not an
  aberration. Only the *variation* drives a contrast floor or a phase-shifting
  interferometry systematic. Conflating the two into one RMS is the most
  likely way to produce a plausible-looking wrong number.
* **At round-off, the statistic can be louder than the physics.** Where a
  quantity is uniform to a few times machine epsilon, an ordinary `mean()`
  over ~10⁴ pupil points accumulates more summation error than the spread
  being measured. Where that happens it is called out explicitly and the
  floor itself is published alongside the result, so the reader can see which
  of the two a number is reporting. The transmission-uniformity result in the
  unitarity gate is the worked case.
