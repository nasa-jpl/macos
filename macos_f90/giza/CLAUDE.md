# Giza / PGPLOT-API graphics — subsystem cheatsheet

> **Post-compaction / resume:** if you are resuming any plot/display work (`macos_f90/giza/`, `pgplotsub.F`, GRAY/CONTOUR/SLICE/WIRE/SPOTDIAG, IMGMODE, wedge/annotation, DRAW data-only), re-read
> THIS file — nested CLAUDE.md files are NOT auto-injected after
> compaction; they reload only when CC next reads a file in this
> directory. Root rules, conventions, and the subsystem index live in
> `../../CLAUDE.md`.

*Sections below are lifted verbatim from the former monolithic root
CLAUDE.md. Move text, don't paraphrase — engine gotchas are exact.*

---

> Note: `draw_rays` (DRAW data-only) straddles `macos_api_mod` — kept
> here because it guards the Giza render path. RPATH/build sections stay
> in root.

## draw_rays — DRAW's ray bundle as DATA (data-only mode)
The design layer's MATLAB layout viewer needs DRAW's *real* ray bundle
(`DrawRayVec`: per-ray, per-surface positions) WITHOUT opening a Giza
device.  Mechanism (additive — `DrawDataOnly` defaults `.FALSE.`, so
normal interactive DRAW is byte-unchanged):
- `src_mod.F`: `DrawDataOnly` flag + capture arrays `DrawRayVec_save` /
  `DrawEltVec_save` / `nDrawElt_save` + inputs `DrawDataPlane` /
  `DrawDataStartElt` / `DrawDataEndElt`.
- DRAW handler (`macos_cmd_loop.inc`): when `DrawDataOnly`, take the
  plane (1=YZ/2=XZ/3=XY) + elt range from the module vars instead of
  `IACCEPT_S`/`CACCEPT` (the SMACOS LoadStack pushes **3 ints + 1 char**
  but the handler reads only 2 ints + 1 char — the plane CARG would be
  misread off the stack, hence module vars).  After the trace it copies
  the bundle into the save arrays, and the `GRAINI` + `DRAW`/`DRAWOUT`
  render is guarded `.AND.(.NOT.DrawDataOnly)`.
- `macos_api_mod.F90`: `draw_rays_cmd(OK,plane,iStart,iEnd,nDEmax,nRay)`
  sets the module vars + `DrawDataOnly` + `CALL SMACOS('DRAW',…)` +
  resets the flag; `draw_rays_get(OK,RayU,RayV,EltVec,nEltPerRay,nDE,nDR)`
  copies the saved bundle out (codegen-friendly rank-2).
DRAW natively supports plane (`pgplotDrawPlane`) + first/last-elt
slicing, so off-axis "XY after the SM, PM hidden" deconfliction views
come for free.  Landed 2026-06-19 (commit f3e98e5; see
PLAN_DESIGN_LAYER §8 Sprint-4 layout realizability).

## Display polarity + plot annotations
- `plot_mode_mod.F` (registered early in `COMMON_SOURCES`) holds session
  state:
  - `ifAstroMode` (default `.TRUE.` = astronomy / negative polarity,
    large→dark — matches legacy PGPlot users' expectation)
  - Bottom-of-plot annotation: `annotOn`, `annotLine1`, `annotLine2`
  - `wedgeLabel` (default `'pixel value'`) — text under the color wedge
  - `ClearAnnot()` resets all three (annotOn + both lines + wedgeLabel)
    after the draw routine emits them.
- `IMGMODE` command toggles polarity interactively (prompts `NEG|POS`,
  accepts `ASTRO|CONV` as synonyms). Must be dispatched BEFORE the
  `LCMP(command,'IMG',3)` branch in `macos_cmd_loop.inc` — placed at
  ~line 4362 for that reason (otherwise "IMGMODE" is absorbed by the
  3-char 'IMG' match and treated as a plot command requiring `ifLoad`).

### Raster inversion — the key Giza gotcha
Giza's `giza_render_gray` (backing `PGGRAY`) internally calls
`giza_set_colour_table_gray()` and forcibly resets its own ramp,
**ignoring any `PGSCR` setup the caller does** (see
`macos_f90/giza/src/giza-render.c:368`). You cannot use `PGGRAY` for
polarity control. The fix in `pgplotsub.F:GRAY`:
- Render the image via `PGIMAG` + a 2-point `PGCTAB` we build ourselves
  (both the color and gray paths).
- For gray: `CTAB_R/G/B = [1.0, 0.0]` (white→black) under astro mode,
  `[0.0, 1.0]` under POS.
- For color: walk `CInd` in reverse under astro mode (palette flipped
  end-to-end; hue order preserved).

### The second Giza gotcha: `PGCTAB` clobbers CI 0 and 1
`giza_set_colour_table` (called by `PGCTAB`) ends with
`_giza_set_range_from_colour_table(_giza_colour_index_min,
_giza_colour_index_max)` — that range starts at **0**. A 2-point ramp
written with `L=[0,1]` overwrites CI 0 and CI 1 with the ramp endpoints,
clobbering the bg=white / fg=black that axes, tick labels, wedge axis,
and `PGMTXT` rely on. Symptoms (all debugged this way): axis-label text
invisible on first plot, annotation visible but axes ticks gone on a
second plot, etc.
- Fix: **reserve CI 0=white and CI 1=black** via explicit `PGSCR` at
  start of GRAY, then `PGSCIR(16, 255)` pushes the image ramp into CI
  16..255 so `PGCTAB` can't touch CI 0/1. Always use `PGSCI(1)` (never
  `PGSCI(ICILO)`) for axes/box/labels/wedge/annotation.
- Re-assert `PGSCI(1)` before `PGWEDG` and before `PGMTXT` — the
  intervening `PGCTAB`/`PGIMAG` may shift the current CI.

### Non-raster plots — always black ink on white
`CONTOUR`, `SLICE`, `WIRE`, `SPOTDIAG`, `PLOTCOL` all set
`PGSCR(0,1,1,1); PGSCR(1,0,0,0); PGSCI(1)` at entry — overrides any
`IMGMODE` setup. Harmless under PGPlot (same as its default).

### Labels, titles, wedge
- `PGLABEL(' ', ' ', CTITLE)` in GRAY suppresses the placeholder
  `'X-Axis'` / `'Y-Axis'` arg strings (the 76 callers that pass those
  literals are untouched — GRAY just ignores them). Only the title
  survives above the plot.
- `PGMTXT('B', 3.5, 0.5, 0.5, line)` centers the annotation under
  the x-axis (coord=0.5 is middle, fjust=0.5 is center-justify).
- Tick numeric labels come from `PGENV(...,axis=0)` → internally
  `PGBOX('BCNST', ...)` — the `N` option prints numeric values at
  each major tick. Keep as-is; for OPD the tick values are projected
  dimensions on the reference surface (usually the OPD reference,
  e.g. Elt 9 in SegDemo3) which is sometimes meaningful.
- Titles:
  - OPD → `'OPD, Surface N'` (no `File=` suffix)
  - INT / PIX → `'<kind>, Surface N, File=<filnam>'` (via `SrfOut`
    and `PixOut` in `utilsub.F`; `File=` kept because the in-file
    stretch variant — `LOG10 Intensity` etc. — can distinguish runs)
- Wedge labels (`wedgeLabel` in `plot_mode_mod`):
  - OPD with BaseUnits → `'OPD (<cUnit>)'` (e.g., `OPD (m)`)
  - OPD without BaseUnits → `'OPD'`
  - INT → `'Intensity'`
  - PIX → `'Pixel value'`
  - Others → default `'pixel value'`

### Annotation lines (bottom of plot)
Command handler sets `annotLine1/2` + `annotOn=.TRUE.` before the draw
call. `GRAY` and `SPOTDIAG` emit via `PGMTXT` after `PGLABEL`, then
`ClearAnnot()`.
- OPD → `OPD=<rms> <cUnit> RMS, <pv> <cUnit> P-V` (from `RMSWFE`,
  `WFEPV`, `cUnit`; falls back to no-unit form if `.NOT.BaseUnits_FLG`)
- SPOT → `RMS spot radius=<r> <cUnit>` — computed from `RaySpot`
  about the centroid in the SPOT handler (not a pre-existing var).
- INT → `Peak pixel=<MaxInt>` (from `Ca2Int` output).
- PIX → `Peak pixel=<maxval(PixArray)>`.
- Phase/Amplitude/etc. intentionally left unannotated.

## Giza plot improvements (macos_f90/giza/)
- Supersampling: off-screen image surface is 2× (GIZA_XW_SUPERSAMPLE) in giza-driver-xw.c;
  `cairo_surface_set_device_scale(surface, 2, 2)` keeps MACOS logical coordinates at 1×
  while rendering at 2× physical pixels.  Downsample to 1× pixmap on flush, then XCopyArea
  or cairo-scale to the window (handles user resize).  Text/font antialiasing: GRAY hinting
  enabled in giza-drivers.c; image data uses CAIRO_FILTER_NEAREST in giza-render.c
  (sharp pixel boundaries for OPD/PIX/INT raster arrays).
- Idle window responsiveness: mhist.c sets rl_event_hook to `_macos_rl_event_hook` which
  (via weak `giza_process_events`) calls `_giza_process_events_xw` — drains ConfigureNotify,
  Expose, and WM_DELETE_WINDOW while readline blocks at the prompt.
- Close button (WM_DELETE_WINDOW): unmaps the window and sets XW[id].window_hidden=1;
  the next `_giza_flush_device_xw` re-maps it.  Does NOT call `giza_close_device()` —
  MACOS tracks its own `ifPlot` flag and skips `GRAINI`/`PGBEGIN` after first plot.
  Closing the device would strand MACOS with spam of "No device open" errors.
- `giza_set_colour_representation` (giza-colour-index.c) no longer requires an open
  device — `colourIndex[]` is static global, matching PGPLOT's pgscr semantics.
  `giza_set_colour_index(ci)` is still guarded (uses Cairo context).


## smacos_dvr graphics initialization
- smacos_dvr.F calls `SMACOS()` directly; unlike the interactive macos.F command loop
  it never opens a Giza device. Any plot-producing SMACOS command (OPD, SPOT, INT, ...)
  then emits "No device open" for every Giza call. Fix: after `macos_init_all(modelSize)`,
  set `nPgPanel = 1` and call `GRAINI` (opens device via `PGBEGIN(0,'?',1,1)`).
- `macos_init_all` does NOT include `macos_init.inc`, so `nPgPanel` is not defaulted to 1
  — it stays at 0 and `GRAINI`'s panel-layout if/elseif (1/2/3/4) silently skips
  `PGBEGIN`. Consumers outside the interactive path must set `nPgPanel` themselves.

## Giza GRAY CR/CG/CB array bounds
- macos_f90/pgplotsub.F GRAY subroutine caches default color representations for
  restoration when exiting color mode. It calls `PGQCIR(ICILO,ICIHI)` then loops
  `Do J=ICILO,ICIHI; CALL PGQCR(J,CR(J),CG(J),CB(J))`. Under Giza, `PGQCIR` returns
  ICILO=0 (GIZA_COLOUR_INDEX_MIN), whereas classic PGPLOT starts at 1 — so CR/CG/CB
  must be declared 0-based: `REAL, Save :: CR(0:511),CG(0:511),CB(0:511)`. Same total
  size (512), but bounds now cover index 0. Otherwise the first-entry gray image crashes
  on CR(0) out-of-bounds write (hard to debug inside giza since graphics libs are
  typically stripped).


## Plot routine viewport (pgplotsub.F PANEL_ENV)
All raster/contour/spot plot routines in pgplotsub.F (CONTOUR,
SPOTDIAG, LINSPOTDIAG, SLICE, GRAY, WIRE, PLOTCOL) call `PANEL_ENV`
instead of `PGENV` directly. The helper:

- In single-panel mode (`nPgPanel <= 1`) calls `PGENV` unchanged.
- In multi-panel mode (`pgp 2/3/4`) expands PGENV into the explicit
  sequence `PGPAGE + PGVSTD + PGQVP + PGSVP + PGSWIN/PGWNAD + PGBOX`
  with the viewport's horizontal extent shrunk to 78% of the standard
  panel viewport. The 22% remaining width inside the panel is where
  PGWEDG (the colorbar in GRAY) lives without crashing into the next
  panel's left tick labels.
- Subtle bug found during development: `PGENV` implicitly calls
  `PGPAGE` to advance to the next panel. Custom replacements MUST
  include the explicit `CALL PGPAGE`; without it, every plot lands
  in the same panel.
- Sanity-checks the input bounds (`SaneBounds` contained function)
  for NaN, reversed (x1 >= x2), and overflow (|x2-x1| > 1e30) — emits
  `** Error: inconsistent propagation parameters` then lets the plot
  proceed with the original bounds so the historical giza
  `giza_set_window: Invalid arguments` is still produced. Catches the
  pathology where a FarField reference surface has `zElt = INF` and
  the propagation math returns garbage extents.

DRAW (system-layout sketch) still uses `PGENV` directly — it's a
single-plot routine, no multi-panel concern.

## GRAINI graphics-device caching (pgplotsub.F)
`GRAINI` is called once per `ifPlot=.FALSE.` event (e.g. on the
first plot after `PGP` changes layout, or after `NEWPAGE`). The
first call passes `'?'` to `PGBEGIN` so the user is prompted for
the graphics device (`/xw`, `/png`, etc.). After PGBEGIN returns,
`PGQINF('DEV/TYPE', dev_str, dlen)` queries the user's choice and
caches it in a `SAVE`'d local. Subsequent GRAINI calls pass the
cached `dev_str` to PGBEGIN instead of `'?'`, so giza doesn't
re-prompt every time the layout changes. With giza the new PGBEGIN
opens an additional window (it doesn't close the previous one);
that's intentional and supports the multi-window-history workflow
in CoroExample.jou-style scripts.

