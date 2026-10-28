# mmacos manual, Part I round 2: report for Dave (TO, 2026-10-28)

This report takes the Part I report's form (REPORT_PART_I.md): the SHAs, what was regenerated, the page count, the open
questions in one list, and the findings by letter. The brief was `BRIEF_to_mmacos_manual_round2.md`, your rulings
of 2026-10-26. All commits are on `dev-candidate`. CCL reviewed and pushed A to D on your word. The Part II reruns
commit is local.

## What landed

| Item | What | macos | resources |
|---|---|---|---|
| A.1 (kkk) | The sensitivity harvest and `run_compare`'s engine truth share the engine's chief-ray reference; each harvest stamps the reference it ran at | ef29af59 | 1f32955, 6efff8f |
| A.2 (fff) | The zero gap-ray clip is diagnosed: there is no library defect. ColSource's acceptance radius changed in 36a5f5aa (12520 → 12128 source rays) | ef29af59 | none |
| A.3/A.4 (hhh) (iii), ch1b | The MET chart shows every bar; the hub inset names the deck's units; the report labels its plate-scale EFL | ef29af59 | 692deab |
| (www), from the repin | Analyses remove wavefront piston by default, because the chief-referenced harvest puts a piston in every column | 0d36137d, a94b0b76 | 2b70b8d, 6556b85 |
| B (bbb) (ggg) | The e2e template's committed outputs (356 files) were removed with a verdict for every one; the MET preset was promoted; the fingerprints were re-pinned; chapters 1–7 were regenerated and re-pinned | 995f51be, a94b0b76, ca65dc71, 48ef7779 | 458f398 |
| C | Chapter 2's coronagraph is now case C: the CTB behind the e2e telescope, with OAP1 at f = 0.3819 m. The chapter shows the ladder, the pupil walk and the arcsecond field. 18 pins | 130027f9 | none |
| D | New Section 8, coronagraph analyses on the CTB: the mask families, EFC on a Jacobian measured live, and vortex EFC with and without the coated train's Jones screens. 24 pins, each equal to the August record | 41fd3356 | none |
| Reruns | Every Part II page was rerun, and each template's committed outputs got a verdict | d1a15f43 | none |

`make all pdf`: 125 pages. Sections were renumbered: 8 is the CTB chapter, 9 Part II, 10 the tools, 11 Part III. Figures 21–23 are new, and Part II's figures are 24–54.

## The Part II reruns

All 33 pages in the catalog and tools sections ran on scratch copies of their templates (`mm_page`; the library's
`git status` was checked after every run). The listings, pins and page figures all held.
- The 29 re-rendered page figures came back pixel-equal and were restored to HEAD.
- Three pages failed only on wall-clock seconds printed in their listings: `tg_psi_dm`, `tg_psi_dm_v2` and `cophi_dm96`. Their specs now mask those seconds, and the three pages were re-pinned. No pin moved.

The record is `cli_tests/records/run_artifacts_part2_reruns_2026-10-28.tsv`. It has one row per committed output file
of each of the 22 templates that commit outputs (the other 11 commit none), compared against the page's run. Nothing was removed.

| Verdict | Files | Meaning |
|---|---|---|
| exact | 25 | The run rewrote the file and it matched: bytes, masked timestamps and paths, pixels, or numbers to 1e-12 abs / 1e-9 rel |
| not-exercised | 2530 | The page's run did not rewrite the file, so this rerun does not test it |
| diagnosed | 43 | The file differs and the cause is named with a SHA and a number |
| held | 35 | The file differs and the cause is not established; the candidates are named |

The diagnosed ones are stale outputs left behind by later library or engine fixes:
- `bench_ifo` and `bench_ifo_dm`: refractor extinction (cafc53c);
- `e5_pie`: ColSource (36a5f5aa);
- `e5_seg`: edge sensors and the grid frames (b41503d, 20d7f03);
- `tg_psi_dm96`, and the `tg_psi_dm` and `tg_psi_dm_v2` decks: QWP placement (9b2181c). Only the QWP vertex moves;
- `tg_psi_dm96_oap`, lens run: beam radius (1c63f7d, ba86a28);
- `tma_longslit` and `e2e6m_r2`: their parameter structs gained fields;
- `bench_ifo_pol`: numerical zeros only.

`e2e2` reproduced exactly.

The held ones:
- **`tma_onaxis`, `tma_offaxis`, `sz_tma`.** Their committed decks differ from today's design-layer build. The scripts have not changed since, so the change is somewhere in `macos.design`, which I have not bisected.
- **`bench_ctb`, geometry layer.** The staged solve lands about 5e-6 relative away. The candidate is the LM stopping rule (84a91903). The committed PNGs are exactly 4/3 today's size; the candidate is a HiDPI host.
- **`e2e6m`, S1/S2.** This is the open round-1 finding (qqq): today's S1 solves to f/11.76, where the committed run reached f/25.39.
- **`tg_psi_dm`'s mat and two figures.** The orientation search now picks pixel axes `x=-col,y=-row` where the August run picked `x=+col,y=+row`. That is the same map turned 180° in pixel indices: the correlation is unchanged at 0.998434, the README's value. The flip predates 10-24, since the page committed then already printed it. The candidates are the engine's frame changes since August; not bisected.

None of these is removed: the record is the verdict, and fixing each one is a separate decision for you.

## Gate

No library change in the reruns commit. The changes are the regen tool's `MMMAN_KEEP` option (keep a passing run's
scratch tree for the audit), three spec masks, the three re-pinned listings, the record, and OPEN_QUESTIONS.

## Open questions, in one list

1. **The three held design-layer pages** (`tma_onaxis`, `tma_offaxis`, `sz_tma`). Should I bisect `macos.design` and regenerate their committed decks, or record the difference and leave them?
2. **`bench_ctb`'s geometry outputs** are 5e-6 relative from today's solve. Should I regenerate them with (zzz)'s re-pin, since that re-pin touches the same template?
3. **`e2e6m` (qqq)** is still open from round 1.
4. **`tg_psi_dm`'s 180° pixel-axis flip.** The science is unchanged. Should I find which change caused it?
5. Carried: the authorship block, `elt_csys`, the perturb frame (PLAN_CONSOLIDATION 1.4a).

Ruled on 2026-10-28 and queued in this order:
1. (zzz): re-pin the CTB hand decks to the generator.
2. Three library fixes: the offset_imager README, `add_pupil`'s fex default, and `run_budget` without a spec.
3. (aaaa): `ctb_efc` checks its Jacobian cache.
4. Pairwise probing: first on the CTB, then the fold, then the coupled model, then pairwise probing on the coupled system.

## The findings by letter (PLAN_CONSOLIDATION 1.2)

- **(uuu)** e2e s1's `rms_waves` masks validity with `W ~= 0`. **(vvv)** The same mask appears in 34 library files. The sweep is entered in 3a.
- **(www)** The chief-referenced harvest puts a piston in every column, and three consumers assumed piston-free columns. Fixed.
- **(xxx)** `bench_ctb` takes the e2e exit beam only if OAP1's focal length changes (f = f/# × Dpup = 0.3819 m). The exit beam is anamorphic: f/17.69 × f/18.47.
- **(yyy)** The committed MET preset differs from today's live search (the class-1 angles). It was promoted in B.
- **(zzz)** The generator's FEX reference-sphere radii have drifted from the committed CTB decks (tCtbProp fails). Ruled RE-PIN.
- **(aaaa)** `ctb_efc` reads its Jacobian cache by file name only. Ruled FIX. Section 8 states the workaround until the fix lands.

## One thing I got wrong, caught before commit

Writing Section 8, I first said that `ctb_efc` checks its cached Jacobian the way its siblings do. It does not. I
corrected the prose before commit, and the gap became (aaaa). Also, in a draft of this audit I counted the pdi and
zwfs files as exact when their runs had not rewritten them. The record now labels such files not-exercised.
