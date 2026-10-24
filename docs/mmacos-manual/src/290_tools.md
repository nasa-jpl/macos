## SECTION 9: The Tool Templates

The other seventeen templates are tools, not designs: examples of an interface, a method, or a driver that
the design templates and the runners call. Two are the entry points of Section 1's thread and get a page each;
the table lists all seventeen, each described by its own README's title or, where it has none, by its main
script's header. A last table lists the library functions no template calls.

### The seventeen tools

| Template | Section | What it is | Described by |
|---|---|---|---|
| `design_layer_api` | 10_telescopes | Design-layer API examples (`+macos/+design`) | README |
| `pupil_id` | 10_telescopes | beyond FEX: the exit pupil as a SURFACE | README |
| `grid_surface` | 20_segmentation | GRID-DATA SURFACE API walkthrough on FFSegDemoAll.in | `grid_surface_example.m` header |
| `coro_experiments` | 30_instruments | Coronagraph division-of-labor experiments (Sprint 1, E1–E4) | README |
| `dm_gauge_lib` | 40_benches | the ONE copy of the DM-gauge scoring machinery | README |
| `e5hex1` | 50_sensitivities | graphical dw/d(grid-data) sensitivity demo + verification | `run_dwdgrid.m` header |
| `e5hex1_budget` | 50_sensitivities | the live-formula WFE calculator (Excel or LibreOffice; regenerate with `run_budget_e5hex1`). | README |
| `e5hex2_refzern` | 50_sensitivities | conforming Reference → GS Zernike basis → dw/d(grid) | README |
| `gen_segment_gridmat` | 50_sensitivities | per-segment Zernike GridMat generator | README |
| `run_dwdgrid_multi` | 50_sensitivities | multi-field dw/d(grid-data) sensitivity (example) | `run_dwdgrid_multi.m` header |
| `run_dwdgrid_multi_multisegbasis` | 50_sensitivities | dw/d(grid) with a PER-SEGMENT basis | `run_dwdgrid_multi_multisegbasis.m` header |
| `run_dwdgrid_multi_singlesegbasis` | 50_sensitivities | dw/d(grid) with ONE shared basis | `run_dwdgrid_multi_singlesegbasis.m` header |
| `run_dwdsurf_multi` | 50_sensitivities | multi-field dw/d(Kr,Kc) POWERED-SURFACE sensitivity (example) | `run_dwdsurf_multi.m` header |
| `run_dwdx_multi` | 50_sensitivities | multi-field dw/dx RIGID-BODY sensitivity (example) | `run_dwdx_multi.m` header |
| `run_dwdz_multi` | 50_sensitivities | multi-field dw/dz ZERNIKE-COEFFICIENT sensitivity (example) | `run_dwdz_multi.m` header |
| `zoom_5x5` | 50_sensitivities | multi-configuration × multi-field sensitivities | README |
| `view_rx_demo` | 60_visualization | the general prescription visualizer | README |

`40_benches/vsg_wip` is a planning folder (a VSG2 model plan and review notes), not runnable; it has no page.


