# BRIEF: dwd* plots large enough to be interpretable at any segment count

From CCL for Dave, 2026-09-10.  Extends `BRIEF_dwd_plot_pagination.md`
(2026-08-28, queued, never executed).  Dave: "change the way dwd* data
is plotted to make the plots large enough to be interpretable, even
with large numbers of segments.  This will mean many more plots and
pages."  Cold start: read `MACOS_resources/mmacos/CLAUDE.md`, memory
`project_opd_conventions` (rounds 2-4), `feedback_deck_plots_unmodified`,
`feedback_luis_facing_figures_jet`, `feedback_look_at_users_own_figures_first`,
and the three plotters below.  Work on `dev-candidate`; commit locally
and report SHAs; push only when Dave says push (`feedback_push_only_on_review`).

## What exists (MACOS_resources/mmacos/sensitivities/)

| file | what it draws today | the problem |
|---|---|---|
| `plot_dw_channels.m` | ONE sheet, every channel a subplot; each subplot = the multi-field canvas (5 fields tiled) | 42 subplots on jwst dwdsurf (18 segs + hub + SM/TM, x Kr/Kc) are postage stamps; dwdx = 138 channels; unreadable |
| `plot_dw_per_element.m` | one page per element (or group); subplots = that element's channels; modes 'center' (centre field, per-field cell) and 'multi' (the tiled canvas via `indxall`); 1400x950 px, r140 | readable for 2-6 channels; a 5-field canvas in a 6-DOF grid is small; pages go to `<name>_pages/` |
| `plot_opd_canvas.m` | the OPDall canvas | fine at 5x5; unbounded for big 'grid','NxM' field sets |

Call site: `run_sensitivities.m` ~lines 566-580 (`plot_dw_channels` per
channel kind + `plot_dw_per_element` per mode in `opts.per_element`,
default `["center" "multi"]`).  Reconstruction rules that MUST be kept:
the stacked multi Jacobian lives on a TILED FIELD CANVAS, rebuilt with
`macos.v2m(J(:,c), out.indxall)`; a per-field map uses
`per_field_indx(out, k)` (orientation-aware; the 2026-09-10 fix -- never
rebuild an index with m2v on a map an 'orient' option may have transposed).
Gate: `tRunSensitivities/test_per_element_page_index_follows_orient_xy`.

## Rulings that bind the design (Dave, 2026-09-09/10)

- Deck plots are the tool's OWN output, unmodified: no re-render, no
  clim, no bespoke masks.  So the TOOL must produce interpretable pages.
- Colormap `jet`, autoscaled; exact zeros not drawn (white); piston
  removed on the per-element pages (existing convention).
- The mean-referenced / PTT-removed columns are the same correct data
  under a different convention -- no "leak" language anywhere.
- Guards warn, never error; existing filenames for small decks stay
  unchanged (small-deck READMEs and artifacts reference them).

## The task

1. **Size first, count second.**  Fix a MINIMUM panel size (e.g. each
   OPD map >= 3.5 in at 140 dpi on a 16:9 page) and derive the panels
   per page from it -- never the other way round.  A 5-field canvas
   panel needs more width than a single-field map: size by the canvas
   tile count.
2. **Paginate on element boundaries** (never split one element's DOF /
   mode / parameter block).  Section by `out.kind` (element / group /
   source blocks).  File names `<name>_<ch>_channels_p01.png` ...; keep
   `<name>_<ch>_channels.png` when one page suffices.  Titles carry the
   page index and the element range.
3. **Per-element pages at full size**: one page per element AND per
   field mode as now, but with the panel size from (1); when an
   element's channels exceed a page (dwdgrid mode counts, group blocks),
   continue onto `_p02`.  Offer a per-FIELD page mode ('field': one page
   per element per field, single-field maps at maximum size) for the
   segment-count regime.
4. **An index sheet** per harvest: one contact sheet listing every page
   (element id, kind, channels, file) so a 200-page output is navigable;
   the report's `figures:` line points at it.
5. Options in `run_sensitivities` and the two plotters: `'panel_in'`
   (minimum panel size), `'max_per_page'`, `'per_element'` gains
   `"field"`; defaults reproduce today's small-deck output byte-for-byte
   (e5hex1 templates) -- assert that in a test.
6. Gates in `tRunSensitivities`: page count and per-panel size on the
   jwst zoom fixture at 63 rays (cheap), element-boundary pagination,
   the orient-xy identity on every page mode, small-deck byte-identity.
7. Re-run `templates/50_sensitivities/zoom_5x5/run_dwdsurf_5zoom_5fov.m`
   and `run_dwdx_5zoom_5fov.m` and LOOK at the pages (Read the PNGs);
   the acceptance is visual: a 19-segment dwdx harvest whose every
   segment's 6-DOF block is readable without zooming.

## Scale to expect

jwst zoom: dwdx 138 channels x 5 fields x 5 zoom configurations; dwdsurf
42; dwdz 19 elts x n modes; dwdgrid 19 x ~6.  Hundreds of pages per
harvest at full size -- that is the point; keep them in `<name>_pages/`
and index them.  Runtime is plotting only (no traces): batch-safe, but
print at r140 not r300 and keep figures 'Visible','off'.
