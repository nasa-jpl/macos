### The functions no template calls

Counted by a word-boundary search of every template script for each function's name: `design/src`, 40 of 75; `sensitivities`, 12 of 24. They are documented here, not deleted (the Part I report's default); a runner or a test may still call them.

**`design/src`** (the design layer's helpers, 40 of 75):

| Function | What its help line says |
|---|---|
| `afocal_first_order` | Paraxial marginal AND chief trace of an N-mirror train. |
| `afocal_ladder_deck` | The 3-rung afocal reference ladder on a COMMITTED deck. |
| `afocal_plane_opl` | Exact OPL from the source to a flat reference. |
| `afocal_refs` | The flat wavefront reference for an AFOCAL exit beam, plus |
| `afocal_rungs` | The three reference-freedom rungs of the AFOCAL WFE metric. |
| `afocal_score_psf` | Score an AFOCAL deck with a focus, without touching it. |
| `afocal_wfe_deck` | The afocal WFE metric evaluated on a COMMITTED .in deck. |
| `calib_blind_fields` | Which CALIB fields did the solve not see?  Gated on RAYS, not on the size of the seed's WFE. |
| `chain_aim` | Launch POINT on a plane whose ray along d hits surface iStop at its vertex. |
| `chain_bundle` | The multi-field, multi-lambda COLLIMATED bundle through a chain. |
| `chain_footprints` | Per-surface beam footprint of a bundle, in the APERTURE frame. |
| `chain_trace` | Sequential exact 3-D trace through a surface chain. |
| `dyson_layout` | Closed-form seed + exact-trace verification of a concentric Dyson. |
| `dyson_scaling` | The block radius a concentric Dyson needs for a given field. |
| `e2e_geom` | Telescope + spectrometer as ONE chain, the grating its stop. |
| `focal_surface` | Measure a system's best-focus SURFACE and emit the deck |
| `optic_footprints` | Trace all fields x all configs; per-optic footprint geometry. |
| `spectrometer_layout_fig` | Deck-standard layout of a spectrometer_geom chain. |
| `spectrometer_maps_fig` | Performance maps of one scored deck, deck standard. |
| `spectrometer_score_chain` | The same spectrometer metrics from the exact chain. |
| `spectrometer_sens` | The tolerance ladder of a spectrometer deck: one perturbation at a time, scored by the engine. |
| `spectrometer_sens_defaults` | The addendum-49 tolerance list for a Dyson deck (spectrometer_rx emission). |
| `spectrometer_sens_dn` | The index change the lens_index row actually applies, over wavelength: the Sellmeier B |
| `spectrometer_slit_loss` | The slit-width diffraction loss, measured vs sinc^2. |
| `spectrometer_wave` | The propagation twin of spectrometer_score. |
| `spectrometer_wave_fig` | The propagation twin's evidence, deck standard. |
| `strict_sphere_opl` | Exact OPL from the source to a reference sphere. |
| `telescope_clear_quick` | A cheap clearance wall for the telescope solve. |
| `telescope_cloud` | The spectrometer's bodies as a point cloud for the telescope's clearance wall. |
| `telescope_geom` | The fore-optics chain that feeds the slit: an off-axis section of a coaxial three-mirror anastigmat. |
| `telescope_ladder` | The telescope's departure ladder on the exact chain, under the slit scorer. |
| `telescope_maps_fig` | The telescope's score at the slit, deck standard. |
| `telescope_score` | The telescope scored at the slit in the ENGINE. |
| `telescope_score_chain` | The telescope scored at the slit on the exact chain. |
| `telescope_seed` | First-order PNP three-mirror layout + Seidel conics, closed form. |
| `tma_screen` | First-order (paraxial, engine-free) clearance screen of a coaxial three-mirror train. |
| `tms_clear` | Beam-leg vs body clearances of a two-mirror TMS chain (tms_geom), exact rays, signed. |
| `tms_firstorder` | General first order of the two-mirror TMS: EFL held by eliminating R1, the stop telecentric. |
| `tms_geom` | The two-mirror modified Schwarzschild (TMS) that feeds the slit, as an exact chain. |
| `tms_paraxial` | Closed-form first order of the two-mirror modified Schwarzschild (flat, telecentric). |

**`sensitivities`** (the sensitivity drivers, 12 of 24):

| Function | What its help line says |
|---|---|
| `dw_block_keys` | Per-channel PAGE BLOCK: which element / group / source a |
| `dw_canvas_tiles` | [rows cols] of FIELD TILES in a harvest's canvas. |
| `dw_draw_map` | One dW map in one panel, in the pinned display convention. |
| `dw_page_axes` | Axes for one panel slot of a DW_PAGE_LAYOUT page. |
| `dw_page_fig` | An invisible figure whose PRINTED size is FIG_IN inches. |
| `dw_page_layout` | Size-FIRST page plan for the dw_d* channel plotters. |
| `dw_page_title` | Page header at a FIXED distance from the top of the page. |
| `dw_panel_label` | Short panel title from a channel name. |
| `per_field_indx` | m2v bookkeeping for ONE field of a dw_d*_multi harvest, in |
| `plot_dw_index` | Contact sheet: every channel as a thumbnail, page-labelled. |
| `run_dwdx_multi_zoom` | -- multi-CONFIGURATION x multi-field dw/dx (GENERIC). |
| `write_page_index` | The harvest's page index: one line per figure written. |

