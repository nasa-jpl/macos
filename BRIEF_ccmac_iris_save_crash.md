# BRIEF — the IRIS `save_rx` → reload SIGSEGV on real ZrnGrData grids (CCMac)

CC for CCMac, 2026-10-10, on Dave's word.  You are on the government cloud with the
IRIS deck, so you can chase this where we cannot; the deliverable that comes back to
the public tree must be free of it.  Read `macos/PLAN_CONSOLIDATION.md` §0 first (the
ten rules — the definition of done), then `macos/REPORT_iris_save_crash.md` (your own
round-2/3 record of this defect) and the IRIS lines of `macos/PLAN.md` §0 (line 29).
Branch `dev-candidate` in both repos; engine tip **47cd323 or later** (TO's phase-1
items are all in: the Rx reader refuses short/non-numeric values cleanly instead of
crashing, line endings are normalized at open, every parser STOP exits through label
99 — rebuild before you start, so you are not chasing a ghost those fixed).

## The defect, as recorded

`iris_dp_ZGD` loads and traces (12760 rays).  `save_rx` it, reload the saved deck,
and the trace SIGSEGVs in `SFFSrf → FreeFormSrf → MonGridSrf → FindSrf → CTRACE` at
one of the six real `ZrnGrData` elements (iElt 17/19/21/35/37/39).  Stop-independent.
The phantom-grid half (38 elements with `nGridMat= 99` / `GridFile= None` gaining a
frame on SAVE) was FIXED 2026-09-09 (`cda178e`); this half was not.  Your attribution:
the original deck's `GridFile=` lines carried a TAB + comment that the GridFile-tab bug
turned into "file not found", so the original NEVER traced those grids; `save_rx`
cleans the name, the reload loads the grid for real, and the live-grid trace hits a
corrupted rewrite of the real-grid block.  Eliminated already: a pre-existing
live-grid trace bug (the original with the grid activated by hand traces fine), the
missing file, the added `lData`.  The public `SegDemo3data.in` (6 grids) round-trips
clean, so the save path is not universally broken.  Unverified candidates from the
report: the `GridSrfOrder= 3` bicubic edge stencil overrunning `GridMat` at the rim;
the rewritten `GridFile` name/path; `nGridMat` against the file's real dimensions; the
`pData/xData/yData/zData` frame written for a ZrnGrData (SrfType 13) element.

## What to do, in order

1. **Reproduce on the current engine, both compilers**, with a debug build (`source
   ./makems.sh debug` and `... debug gfortran`; `-check all` on ifx, bounds checks on
   gfortran).  gfortran's checker has found the real line every time this year — start
   there.  **Three outcomes, each informative:** (i) it still crashes — chase it as
   below; (ii) it is now a clean refusal naming a key (`** Rx load refused: <key>`)
   — the reader has found the corrupted rewrite for you, that key is the bisection's
   answer and the writer is the fault; (iii) it loads and traces — one of last
   night's fixes closed it as a side effect: bisect across TO's kept pre-fix binaries
   (`~/dev/macos/build_release_prefix` … `_prefix8`, one per item, both compilers; on
   the cloud rebuild them from the SHAs in PLAN_CONSOLIDATION §1.1) to NAME the fix,
   and still build the public reproducer (step 4) red on the binary before it.
   Use the pty-driven CLI, not MATLAB, for the first reproduction (a mex crash
   kills MATLAB and hides the stack): `macos/cli_tests/cli_load_gate.py <debug binary>
   <list with the saved deck> out.csv --workdir keep` gives you load → opd → save →
   reload → save in two processes, and `gdb` on the second (the GDB-first rule in the
   memory: `printf | macos` hangs on readline; drive it through a pty).
2. **Bisect the saved deck, not the code.**  Take the saved IRIS deck and the original;
   diff them element by element for the six ZrnGrData elements (keys, order, the
   grid frame, `GridSrfdx`, `nGridMat`, `GridSrfOrder`, `GridFile`).  Then make the
   SAVED deck load by reverting ONE difference at a time until it traces — the first
   revert that cures it names the key.  That is faster than reading the writer.
3. **Name the engine fault** with the debug build's own line and index: the
   out-of-range index, the stale pointer, the frame that never got set.  Only then
   change code.  Keep the fix minimal and in the engine (`iosub.inc` writer and/or the
   grid surface routines in `surfsub.F`); if the fix is a guard, it WARNS and counts,
   it does not refuse (rule 5) — a hard refusal needs Dave.
4. **Build the PUBLIC reproducer.**  A synthetic deck in `macos/ZGD_test_files/`
   (`tst_zrngr_roundtrip.in` + its grid file) with ZrnGrData elements of the SAME shape
   as IRIS's — same `GridSrfOrder`, the same `nGridMat` against the file's dimensions,
   the same frame idiom, a grid that is live from the first load — that crashes the
   PRE-fix engine the same way and round-trips after.  No IRIS numbers, no IRIS
   geometry, no IRIS names: a Cassegrain-class deck with the grids bolted on is fine.
   If the synthetic deck does NOT crash pre-fix, you have not captured the trigger
   yet — keep bisecting (step 2) until it does.  This is the gate's red leg and the
   only thing that can live in the repo.
5. **Gates**, all red on the pre-fix binary (keep it: `build_release_prefix*/` is the
   convention): the CLI round trip of the synthetic deck on both compilers (the
   `cli_tests` load gate records it: status ok, roundtrip `same`); `mmacos/tests/
   tZrnGrRoundTrip.m` (SUITE_FAST) = load → save_rx → reload → trace, with the
   pre-fix crash as a must-fail leg run in a `matlab -batch` subprocess asserting the
   exit code; a pymacos twin.  Then the IRIS deck itself, on your side only, as the
   confirmation that the public reproducer stands for it — report its numbers
   (rays surviving, RMS OPD before and after) in the report, never the deck.
6. **Corpus A/B** if the fix touches `surfsub.F`, `iosub.inc` or `msmacosio.inc`
   (rule 3): `cli_tests/run_cli_tests.sh corpus` from a clean worktree at the pre-fix
   and post-fix commits (`~/dev/macos_cli_base` pattern in `cli_tests/README.md`), every
   moved deck named in a committed record (TO's `cli_tests/records/ab_6a_moved_*.txt`
   is the form).  The fast suite and the pymacos subset after.
7. **Deliverables:** the engine commit(s), the fixture + three gates, the A/B record,
   `macos/REPORT_iris_save_crash.md` updated with the mechanism and the closure (keep
   the existing text; add the section), the PLAN §0 line 29 checkbox, and a one-line
   pointer in `macos_f90/CLAUDE.md`'s "GridSrf" cluster.  Commit locally with SHA +
   branch stated; Dave says when to push.

## Rules that bite here

- The IRIS deck, its grid files and its numbers never enter either repo, a commit
  message, a test name or a report table.  "An instrument deck with six ZrnGrData
  elements" is the most you write.
- A fix ships with a gate that is RED on the pre-fix binary (rule 1); a synthetic
  reproducer that does not crash pre-fix is not a gate.
- Guards warn once per run and count (rule 5); no new refusal without Dave.
- Nothing in `elt_mod_init_vars` that is a user option (rule 6); ≤ 72 columns (8).
- Never relink the shared mex while another MATLAB runs; on the cloud that is your
  own, but say so in the report.
- One item in flight; report to CC (this channel or a NOTE_to_ccl file) with the
  gate's red/green counts on both compilers and the A/B's moved list.
