# Figure rescue — status & audit log

## 2026-10-19 — every engine-made figure regenerated from the current macos (Dave)

26 figures are now PRODUCED BY THE ENGINE on the /png device and assembled by
`tools/regen_engine_figs.py` (its SPEC table: figure -> page directory -> pages side by
side) from the journals in `tools/figs/` -- `FIG03 FIG04 FIG06 FIG10 FIG31 FIG33_34
FIG36_37 FIG39 FIG47_55 FIG56 FIG58 FIG60_64`, each the manual's own example macro
(`examples/*.jou`) or its printed dialog, run by the CLI suite's driver so a failed
EXPECT stops the figure: 3, 4, 6, 10, 31, 33, 34, 36, 37, 39, 47-49, 51-58, 60-64, written
as `src/media/fig_engine_NN.png`.  The 3.2-era `image*.png` and the two rescued crops they
replace (31, 56, 58 were engine-made already) stay in media/ unreferenced until the next
cleanup.  Kept as artwork (diagrams and schematics, not engine output): 1, 2, 5, 7-9,
11-30, 32, 35, 38, 40-45, 50, 59, 65, 66.  Figure 46 (beam profiles) is remade on
CassWithExitPupil.in, its 25-element beam-expander deck being absent from the tree.  Changes beside the pictures: figure 34 is the ROW plot its caption
always named (the 3.2 dialog said WIRE); figures 60-62 are made at the detector (element
6, where the 3.2 pictures were taken) and their dialogs re-emitted; the Section 5 OPD
example's RMS is today's 1.545e-2 (the 3.2 text's 1.609e-2).  Each journal copies its
pages into its own `$TMP/figNN` (the gate's scratch directory is shared by every journal
of a run; the first regeneration mixed pages).  The axis labels' multiplication sign
(finding 1.2 (t)) is fixed in giza the same day: the vendored cairo has no FreeType, its built-in
font is ASCII only, so the escape now draws `x` -- the labels read `-1x10^-10`, not a box.

34 composite line-art figures (FrameMaker-era, exploded to fragments by
the PDF->docx conversion) were rasterized from `docs/macosMan3.2.pdf`
(vector-clean) into `src/media/fig_rescued_NN.png` and inserted before
their captions; fragment debris near the captions was removed (audit
trail below).

## 2026-10-17 — three figures that had no image, hand-cropped

Figures 5 (exit pupil for Cassegrain.in, p. 18), 26 (conic mirror
parameters, p. 65) and 29 (obscuration types, p. 74) existed in the
Markdown only as label debris: `pdftoppm -r 150` of the page from
`docs/macosMan3.2.pdf`, cropped with PIL, saved as
`fig_rescued_05/26/29.png`, inserted before the caption.  The label
debris around figures 5, 6, 11, 12, 18, 26, 27, 29, 30, 65 and 66 was
removed by hand (36 tiny glyph images on 15 lines, the remaining
"equation images" of the README), with the running-header artifacts
('Background', 'Commands for Building Models', 'Element Surface Types',
'Ray-Trace Commands') and the `{.underline}` wrappers on four prose
lines.  The two-panel OPD figure 6 keeps its panel labels as an italic
line under the image.  GOOD is now 34.  Later the same day: figure 59 (complex plane, p. 132)
cropped from the PDF; figures 31 (DRAW views of ImageDemo.in), 56 (GBLUR
example) and 58 (Luneberg gain) have NO artwork in the 3.2 PDF either —
their captions sit over blank space, which is why every crop of them
caught text — so they are PRODUCED BY THE ENGINE: `tools/figs/FIG31.jou`,
`FIG56.jou`, `FIG58.jou` through `cli_tests/cli_cmd_gate.py` on the /png
device, assembled by `tools/regen_engine_figs.py` (31: the two views side
by side; 56: the two composites; 58: the column cut of the gain).  Nothing
is pending.

## Review status (2026-07-06)

- GOOD (31): figs 1, 2, 8, 9, 11-18, 21, 22, 24, 25, 27, 28, 30, 32,
  35, 38, 40-45, 50, 65, 66 — clean artwork; a few carry one faint
  line of stray body text at the top edge (cosmetic).
- NEEDS HAND-TIGHTENING (3):
  - fig 31 (DRAW views of ImageDemo.in): crop caught the draw DIALOG
    and prose instead of the two view plots.
  - fig 56 (GBLUR example): crop is transcript text; artwork not
    located in either source PDF page.
  - fig 59 (complex plane): usable diagram at bottom but crop spans
    unrelated content (incl. FIGURE 58 caption).
- To re-crop one figure: tools/rescue_figures.py --redo <pdf> <bbox> N
  (or hand-crop from a page render and overwrite the PNG; the md
  reference stays valid).

## Debris-removal audit

## FIGURE 1 (02_technical_overview.md)
- trimmed imgs: ![](media/image1.png){width="0.10551181102362205in" height="...
- dropped label: '2.5'

## FIGURE 2 (02_technical_overview.md)
- trimmed imgs: Ellipsoidal mirror Paraboloidal mirror...
- dropped label: 'Ellipsoidal mirror Paraboloidal mirror'

## FIGURE 8 (02_technical_overview.md)

## FIGURE 9 (02_technical_overview.md)
- dropped label: '**Physical Optics: Diffraction**'
- trimmed imgs: L...
- dropped label: 'L'

## FIGURE 11 (04_describing_optical_systems.md)
- dropped label: 'zSource=1d22'
- dropped label: 'Collimated beam'
- dropped label: 'zSource<0 zSource>0'
- dropped label: 'Diverging beam Converging beam'

## FIGURE 12 (04_describing_optical_systems.md)
- dropped label: 'Diverging beam Converging beam'
- dropped label: '(zSource = ∞)'
- dropped label: '(zSource \\< 0)'
- dropped label: '(zSource \\>0)'

## FIGURE 13 (04_describing_optical_systems.md)
- dropped label: 'GridType=Circular GridType=Rectangular'
- dropped label: 'GridType=Flower'
- dropped label: 'GridType=Hex GridType=Pie'

## FIGURE 14 (04_describing_optical_systems.md)
- dropped label: 'xGrid'
- dropped label: 'zGrid'
- dropped label: 'Chief ray'
- dropped label: 'xGrid'
- dropped label: 'nGridPts=7'

## FIGURE 15 (04_describing_optical_systems.md)
- trimmed imgs: Secondary Mirror...
- dropped label: 'Secondary Mirror'
- dropped label: 'Detector'
- dropped label: 'Reference Surface'
- dropped label: 'Primary Mirror'

## FIGURE 16 (04_describing_optical_systems.md)
- dropped label: '0.15'
- dropped label: '1.5'
- dropped label: 'Z'
- trimmed imgs: 3.0...
- dropped label: '3.0'

## FIGURE 17 (04_describing_optical_systems.md)
- trimmed imgs: This ray hits the elements out of sequence, becomes undefine...
- dropped label: 'lens_back'
- dropped label: 'lens_front'

## FIGURE 18 (04_describing_optical_systems.md)
- dropped label: 'sArrayWidth'
- trimmed imgs: Top View, hex array (LensArrayTyp=2)...

## FIGURE 21 (04_describing_optical_systems.md)
- dropped: ![](media/image19.png){width="0.13897637795275591in" height=

## FIGURE 22 (04_describing_optical_systems.md)
- dropped: ![](media/image24.png){width="0.17913385826771652in" height=

## FIGURE 25 (04_describing_optical_systems.md)
- dropped label: 'Ray B'
- dropped label: 'Ray A'
- dropped label: 'Face1'
- dropped label: 'Face2'

## FIGURE 27 (04_describing_optical_systems.md)
- dropped: ![](media/image2.png){width="0.11023622047244094in" height="
- dropped label: 'xMon'
- dropped label: 'mirror segment'

## FIGURE 28 (04_describing_optical_systems.md)
- dropped label: 'Enter Element 1 Data:'

## FIGURE 30 (05_ray_trace_analysis.md)
- dropped label: 'RayPos3'
- dropped label: 'z'
- dropped label: 'y'
- dropped label: 'MACOS\\>**draw**'

## FIGURE 31 (05_ray_trace_analysis.md)
- dropped label: 'MACOS\\>**draw**'

## FIGURE 32 (05_ray_trace_analysis.md)
- dropped: ![](media/image24.png){width="0.18346456692913385in" height=
- trimmed imgs: Misses element Must go backwards to hit next element....

## FIGURE 35 (05_ray_trace_analysis.md)
- dropped label: 'Element'
- dropped label: 'Projection plane'

## FIGURE 38 (05_ray_trace_analysis.md)
- dropped label: '*New source position, angle*'

## FIGURE 40 (06_diffraction_analysis.md)
- dropped label: 'Secondary'
- dropped label: 'Spiders'
- dropped label: 'Fold mirror'

## FIGURE 41 (06_diffraction_analysis.md)
- dropped label: 'zElt4'
- dropped label: 'zElt3'
- dropped label: 'zElt8 zElt11'
- dropped label: 'zElt10'
- dropped label: 'zElt13 D 14'

## FIGURE 42 (06_diffraction_analysis.md)
- dropped label: 'iElt 1 PropType'

## FIGURE 43 (06_diffraction_analysis.md)
- dropped label: 'start'
- dropped label: '[z1]{.underline}'
- dropped label: 'start'

## FIGURE 44 (06_diffraction_analysis.md)
- dropped label: 'iElt=5'
- dropped label: 'iElt=8'
- dropped label: 'Example A Example B'

## FIGURE 45 (06_diffraction_analysis.md)
- dropped label: 'z1 z2'
- trimmed imgs: Through-focus beam z1 \> 0, z2 \< 0...
- dropped: ![](media/image85.png){width="7.716535433070866e-2in" height

## FIGURE 50 (07_beam_propagation_imaging.md)
- dropped label: '**Beam Propagation Commands**'

## FIGURE 56 (07_beam_propagation_imaging.md)

## FIGURE 59 (07_beam_propagation_imaging.md)
- trimmed imgs: imaginary...
- dropped label: 'imaginary'
- trimmed imgs: (a + b⋅i)...
- dropped label: '(a + b⋅i)'
- dropped label: 'b'
- dropped label: 'a real'

## FIGURE 65 (08_differential_ray_tracing.md)
- dropped label: 'RayDir'
- dropped label: 'RayPos'
- dropped label: 'Ray-trace'
- dropped label: 'Ray state'

## FIGURE 66 (08_differential_ray_tracing.md)
- dropped label: 'dγ'
- dropped label: 'dri'
- dropped label: 'Perturbed surface'
- dropped label: 'ri'

## FIGURE 24 (04_describing_optical_systems.md)
- dropped label: 'SegYgrid'
- dropped: ![](media/image54.png){width="0.2826388888888889in" height="
- dropped label: 'SegXgrid'
- dropped label: 'GridType=Flower'
