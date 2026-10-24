# mmacos manual, Part II: report for Dave (TO, 2026-10-24)

The form of the Part I report (REPORT_PART_I.md): the SHAs, what was generated, the page count, the
open questions in one list, the findings by letter. All commits are on `dev-candidate` in `~/dev/macos`;
CCL reviewed each batch and pushed it on your standing word. Nothing in MACOS_resources was changed:
every template ran unchanged on a copy, and the library's `git status` was compared before and after
every run, failed runs included. A run that moved it would have been refused; none did.

## Part II: the template catalog and the tools, 116 pages in all

`make all pdf` builds `docs/mmacos-manual/build/mmacosMan.pdf`: Part I's 45 pages, then Section 8 (the
catalog) and Section 9 (the tools). Figures 20 to 50.

| Section | What it holds | Commits |
|---|---|---|
| 8 | One page per DESIGN template, all 35, in the index's order: what it is (its README, or generated text), the knobs, what it prints, its figure, the number it reproduces, and a request you can paste to an assistant to adapt it | batches 91cff24d, 40b07151, 777b1c45, 16901f6b, 9188938b, 1b2ff44e |
| 9 | The 17 tool templates in one table; pages for `design_layer_api` and `pupil_id`; the 52 library functions no template calls (`design/src` 40 of 75, `sensitivities` 12 of 24), each with its help line; `vsg_wip` as a line | d595fefb |

Machinery (fe0db3fa onward): one page generator, `tools/gen_part2.py`, over one runner,
`tools/transcripts/mmlib/mm_page.m`. Each page has a spec in `tools/part2/specs/<name>.json`. The only
typed text is the caption, the number sentence and the prompt. The listings, figures and numbers come
from the run, through Part I's gates (pins to the printed digit or 5 significant figures; numerical
zeros as checked bounds).

## How the pages stand

- **Runs that reproduce their record** (most of them), e.g.:
  - `tg_psi_dm96`: its committed report, line for line;
  - `bench_ctb`: 0.001423 waves, against the README's 0.0014;
  - `tg_psi_dm`: the 7.4786° arm rotation and the 11.7 % high gain.
- **Runs that print today's result beside a record they no longer reach**, each citing its letter:
  - `offset_imager` (lll);
  - `e2e6m` (qqq);
  - `e5_seg` and `e5_pie`, live against stale READMEs;
  - `tma_onaxis`, whose README table predates the engine.
- **Refusals, printed as they run today:**
  - `e2e2` (ppp) and `rc_offaxis` (rrr);
  - `zwfs_dm96` and `pdi_dm96` (nnn) show the bare call refusing, then run the sheet's own value,
    `bench.coat_oap 'none'`, and say so.

  A refusal page fails its own build if the template stops refusing, so a repair forces the page over to
  a real run.
- **`coro_walkthrough`** prints what it does today and cites (ccc).
- **Generated README text:** the 15 designs without a README get it in `docs/mmacos-manual/part2/<name>.md`,
  generated from the entry script's header. Each page says so. The templates were not touched.
- **Reused from Part I:** `e2e` (Sections 1 to 6), and `coro_planet_demo` and `bench_layout` (Section 2).

## The defaults (overturn any in one line)

1. **A template that refuses today gets a page that shows the refusal**, verbatim, and cites its letter.
   Where the README itself gives a value that runs, the page runs that value and says so. No value is
   invented (CCL's ruling).
2. **README text for the 15 designs without one lives in the manual** (`docs/mmacos-manual/part2/`), not in
   the templates.
3. **Heavy campaign benches run their bench stages only** (`cophi_dm96`, `tg_psi_dm96_oap`, `zwfs_dm96`,
   `pdi_dm96`), not their hours-long batteries. Each page's run line says so. `tg96.m` (ten minutes)
   runs whole.
4. **Templates Part I already ran are not rerun** (`e2e`, `coro_planet_demo`, `bench_layout`). Their pages
   point to Part I's generated runs and pins.
5. **The 52 uncalled functions are documented, not deleted.** The table is the input to your pruning
   decision.

## Open questions, in one list

1. **The 15 generated README texts.** Promote them into the templates (a library change), or keep them in
   the manual?
2. **(nnn)** A bare `zwfs_run` refuses its own defaults (`coat_oap 'protectedAl'` on the lens rig), and
   `pdi_run` inherits it. Warn, like `tg96_run`, or change the default?
3. **(ppp)** e2e2's stage 1 still gates `seidel_seed` against the fixture that res 3e96ad8 un-pinned. Should
   the gate be re-targeted the way `tE2E2Axial` was?
4. **(rrr)** `rc_offaxis`'s final design cuts a foreign beam at M2 (−15.7 mm). Fix the design, or retire it
   behind `rc_unobscured`?
5. **(lll)** `offset_imager()`'s defaults are not `run_t3`'s inputs. Which call is the template's run of
   record?
6. **(qqq)** e2e6m's telescope stage departs from its committed report after S0, with its files unchanged
   since 2026-08-25. Should it be traced?
7. **(mmm)** tma_unobscured's closing report reads NaN: all rays are obscured after `add_pupil`.
8. **Figure defects, the same class as (hhh)–(jjj):**
   - `tg96_run`'s layout legend says the test arm is blue and draws it red;
   - `pupil_id`'s Zernike chart puts its axis multiplier under its title.

   The captions say so.
9. **(sss) One `exit(0)` in an example** (`design_layer_api/example_sensitivities_from_rx.m:49`), against
   your no-exit rule.

Carried from Part I and still open: the thread's budget spec, (kkk) (the OPD reference for the harvest),
(bbb) (regenerate e2e's committed outputs), (ccc)/(ddd) (the two coronagraph templates), and the
authorship block.

## The findings by letter (PLAN_CONSOLIDATION 1.2; measured, none fixed by the manual)

- **(lll)** `offset_imager()` with its defaults matches the rodgers3 record at S1 (37.8 nm), then departs (S5 12415.1 against 118.2 nm); the record is `run_t3`'s.
- **(mmm)** `tma_unobscured`'s design report is NaN: all 1185 rays are obscured after `add_pupil` (1135 at the FP, 50 at M3).
- **(nnn)** `zwfs_run`'s coat_oap guard refuses its defaults (zwfs_run.m:195); `pdi_run` inherits it.
- **(ooo)** With POLARIZATION on across `load_rx`, `macos.intensity` read all zeros. **FIXED by CCL (macos c6c83d55, pushed; gate passed).** Both tg_psi_dm pages now print real runs.
- **(ppp)** e2e2's s1 gates `seidel_seed` against the fixture 3e96ad8 un-pinned (max|dK| 3.03e-1).
- **(qqq)** e2e6m's telescope stage departs from its committed report after S0 (f/11.76 against f/25.39; the dense map fails).
- **(rrr)** `rc_offaxis` stops at its own clipping assert (M2 clearance −15.7 mm, 11 rays).
- **(sss)** `design_layer_api/example_sensitivities_from_rx.m:49` ends in `exit(0)`, and its README says all four examples do.
