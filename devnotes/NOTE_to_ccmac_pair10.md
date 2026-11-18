# NOTE to CCMac — your next pair: macos b7556b81 + resources d436b9c (CC, 2026-10-22)

Both on origin dev-candidate now.  Engine CHANGED (TO's cleanup round 1 and three of Dave's rulings): the full
7-check clean-clone gate.  Thirty-one macos commits above your gate-9 base e474a137; the Mac line was declared free
of engine defects at gate 9 -- this pair tests that it stays so.

**Engine changes since ad4161fa (all TO's, each with a journal red on its pre-fix binary; Linux gates PASS both
compilers at each step, corpus 1468/1468, the plan's 1.2 lines name them):**
- (k) DOPD after ZABERR / ZRM reports the real dOPD, not the -1.0 sentinel (f743b307).
- (m) LSPOT honours the NONE plot type (968c36a6); (i) RESET resets the plot type (82ee9d7f).
- (v) api `save_rx` no longer prints `** Unknown command = YES` (98d65f73; resources twins 5a39f47).
- (dd) an unknown `PropType=` is REFUSED at load by name (a239d5e6) -- ImageDemo.in's `Annulus` fixed.
- (f) `GridFile=` is looked up beside the deck (0f9983eb): **16 corpus movers on Linux, all grid decks that were
  silently loading without their grids**, each attributed (the pre-fix binary run in the deck's own folder gives
  the tip numbers) -- `cli_tests/records/ab_gridfile_moved_0f9983eb.txt`; expect the same 16 on arm64, same numbers.
- (zz) coating keys: 65 slots, refused by name above; `Coating=` matches its full name so `CoatIndx=` is reachable
  (b9cdb843) -- the old 20-slot store was written out of bounds at `nCoat= 21`.
- (aaa) the CLI tokenizer: 128 tokens per 256-character line, refused by name beyond (623898d4); the old 10-slot
  store was written out of bounds by any `!` shell line over ten words.  `LOH_tokens.jou` (new) -- on macOS check
  the refusal text and that the next line is read.
- The dead `sv_surfsub.F`, `wk_surfsub.F`, `tst_zern_wf.F` removed (4a140f9a); nothing built them.
- Harness: `pymacos/tests/proper_compare/conftest.py` reads PYMACOS_SRC (d436b9c) -- a whole-suite run used to
  import the SHARED .so; `test_side_tree_is_loaded.py` guards it.  Your clean-clone runs were never affected
  (one tree), but note it for any side-tree work.

**Release prep in the same push (no engine effect):** `release-exclude.txt` amended in both repos (the CLAUDE.md
set, `macos_f90/notes/`, MAC_PORT.md and the optical_design guides SHIP; `PLAN_*.md` added), `tools/release_strip.py`
(dry run by default), `RELEASE_EXCLUSIONS_REVIEW.md`, `RUNBOOK_release_moves_2026-10.md` -- the last two are for
Friday's meeting; the runbook's step 2 describes the re-clone everyone does after the history rewrite, you included.

**Expected on arm64:** cmd 117/117 + 2 must-fail (116 of gate 9 + LOH_tokens; no (ff)-class miss remains), load
160/160 with the 16 GridFile movers vs your gate-9 base (named, same numbers as Linux), mmacos fast 656/0 + the
filtered leg, pymacos 7007 + the new side-tree guard = 7008 / 0 / 2 skips (TO's Linux side tree read 7007 / 1 skip).
If the GridFile movers differ in number or value on arm64, that is the one thing to send.
