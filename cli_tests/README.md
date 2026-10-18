# cli_tests — the CLI regression suite (PLAN_CONSOLIDATION §2)

The interactive `macos` CLI shares its parsers and command loop with mmacos and
pymacos but had no regression suite of its own; the bindings have ~570 fast
tests and pymacos ~6600.  This directory is that suite.  It drives the CLI through
a pty (readline needs a tty; the model-size prompt comes first) and writes RECORDS
— the record is the result, the script prints the summary.

```
./run_cli_tests.sh load        # the core corpus (114 decks), both compilers  (~minutes)
./run_cli_tests.sh corpus      # every committed prescription (~1270)         (~an hour)
./run_cli_tests.sh extended    # + untracked sandbox decks (informational)
./run_cli_tests.sh cmd         # the command journals (cmd/*.jou), both compilers
./run_cli_tests.sh manual      # load + re-emit the manual's examples (docs)
./run_cli_tests.sh ab <pre> <post> [load|corpus] [allow-file]   # the A/B of two records
./make_corpus.sh               # regenerate the three deck lists
```

## Layers

| layer | what | status |
|---|---|---|
| `load` (§2.1) | per deck: LOAD → `opd nElt` → SAVE; new process: LOAD the save → SAVE; byte identity of the second save against the first | built 2026-10-09 |
| `cmd` (§2.2) | one pty-driven journal per command-reference entry (`cli_cmd_gate.py`, `cmd/*.jou`): `# EXPECT:` / `# EXPECT-NOT:` regexes asserted on the transcript, `# COPY: <deck> [<subdir>]` copies a deck into `$TMP` first (anything that writes beside the deck), `$ROOT $RES $TMP $DECK` expand; `must_fail/*.jou` must NOT pass | built 2026-10-09 |
| `ab` (§2.3) | `ab_records.py`: the load-gate record at a PRE SHA against the record at a POST SHA, deck by deck (status, nElt, nPass, lost, RMS, P-V, avg, roundtrip, as printed); every mover written to `records/ab_<kind>_<compiler>_<pre>_<post>.txt`; exit 1 on a mover not named in the allow file or on an ok → not-ok deck.  The builds are the manual step (one clean worktree, one SHA at a time); the records are what gets committed | built 2026-10-10 |

## Writing a journal (Phase B, 2026-10-16)

- Draft the journal with only `# EXPECT-NOT: forrtl`, run it through `cli_cmd_gate.py <bin> <dir> <out.csv>`
  from a scratch directory, read `<out>_transcripts/<name>.txt`, then write the `# EXPECT:` lines from the
  transcript -- never from the manual or the cmdref prose (both have described behaviour that the engine
  does not have: BWK, TPR, SYSPROP after OPD were all found this way).
- **List-directed output differs between compilers**: ifx prints `0.000000000000000E+000` and
  `1.000000000000000E-003`, gfortran `0.0000000000000000` and `1.0000000000000000E-003`; formatted
  (`D` edit) output is identical.  Match a list-directed zero as `0\.0+(E\+000)?`, a mantissa to the
  digits both print (`-4\.41659775810825\d*E-003`), never the whole token.
- A prompt answer that happens to be a command name (`beam`, `no`) is matched by the cmdref's suite
  index only if it passes the dispatcher's minimum match; `cd`, `!` and the like are genuine commands.
- The shell escape `!` runs the rest of the line through `system()` (256 characters since 6ccd86d; it
  was cut at 130 before); `$DECK` is a COPY, so `! sed -i ... $DECK` is the way to vary a deck.
- The driver runs the CLI with `TERM=dumb` (`cli_load_gate.Cli`): under a real terminal type
  readline 8 turns on bracketed paste and every output line after a prompt starts with
  `ESC[?2004l`, so the `^`-anchored EXPECTs (SHELL_LONG, SESSION_INFO, ...) miss -- 24 journals
  failed that way from an `xterm-256color` shell on 2026-10-18.
- A `!` line that starts with `cd` is the CLI's own CD (the rest of the line is dropped): give a
  tool absolute paths instead (`atmos_map_check.py --dir $TMP`).
- Make the gate non-vacuous: pair the asserted state with its contrast (SPCENTER's two spot means,
  BWK's refusal then its walk, the FOV count before and after).
- A TEXT / BINARY / FITS output that already exists asks `exists.  Overwrite? [YES]:` — the next
  journal line answers it (a `! awk` line there is eaten as the answer); put a `yes` before it.  An
  external program the CLI waits for must exit by itself: `vi` with the file name `+q!` does.
- A deck that needs a data file beside it (`ObjFile=`, `.srfN.txt`) gets the file written by a
  `! awk` line after `cd $TMP`, and the key added with `! sed -i '1i Key= file' $DECK`; the gate's
  one journal line is 256 characters.

## Files

- `cli_load_gate.py` — the driver.  One CLI process per leg so a crash at load
  cannot take the rest down; prompt-driven, never sleeps; answers the validator's
  "Pick a different file" re-prompt with `q` (→ `load_fail`) and the plot-device
  prompt with `/null`.  Statuses: `ok`, `load_fail` (validator or parser refusal),
  `crash_load` / `crash_opd` (the process died — the HOST-KILLER class when the
  engine is a mex), `timeout_*`, `missing`.  Round trip: `same` / `differ` /
  `save1_fail` / `save2_fail`.
- `judge_load.py` — the verdict: FAIL on a must-fail deck that loads, on a NEW
  failure against the previous record of the same kind, or on a round trip that
  was `same` and is not; pre-existing failures are named, counted, not fatal.
  Lists the decks whose OPD / nPass / lost moved against the previous record.
- `make_corpus.sh` — `corpus_core.txt` (ZGD_test_files, the manual's examples,
  mmacos and pymacos `tests/Rx`, GMI `test_ff`, segmirmaker), `corpus_committed.txt`
  (every tracked `.in` that declares `nElt=` or `Element=`, both repos),
  `corpus_extended.txt` (+ `$MACOS_SANDBOX`, default `~/dev/MACOS_sandbox`; gitignored —
  machine-local).  The two committed lists carry this box's absolute paths; re-run
  `make_corpus.sh` on another machine before the first record.
- `must_fail/` — the non-vacuity decks (§2.4), regenerated into `list.txt` each run.
  `missing_value.in`: `nGridpts=` with no value — the validator refuses it.
  `segdemo3_broken.in`: the ZGD broken-deck fixture (excluded from the corpus lists
  so it counts once, here).  History: `short_zerncoef.in` (a `ZernCoef=` block one
  value short) crashed both CLIs until TO's item 1 (`ReadCoefBlock`, macos 1d318ac)
  made it pad and warn; it moved to `ZGD_test_files/tst_short_zerncoef.in` as a
  positive case the same day — the pattern for every must-fail deck a fix closes.
- `records/` — `<kind>_<compiler>_<sha>[+dirty].csv` + `.log`, `mustfail_*.csv`,
  `<kind>_<compiler>_latest.txt` (the last summary).  Committed: the records of
  tagged runs; the `.log` files are gitignored.
- `ab_records.py` — the A/B of two records (above).  Non-vacuity: `ab d6e2fc5 47cd323
  load records/ab_6a_moved_47cd323.txt` passes (the 6a round-off movers, all named);
  without the allow file it fails naming them.  Its first run over the core found one
  round-off mover TO's 406-deck record had missed (`Rx_AsphShort_base.in`, gfortran only).
- `ab_corpus_seed.py` — the FwdRoot corpus A/B driver (2026-10-04) this grew from.

## Conventions

- Both compilers, always; gfortran's stricter runtime has found real bugs.  The
  binaries come from a CLEAN worktree at a committed SHA (`MACOS_CLI_TREE`, default
  `~/dev/macos_cli_base`; `git worktree add ../macos_cli_base <sha>`, then
  `makems.sh release` and `makems.sh release gfortran` there) — never from the
  shared working tree, whose build dirs another lane may be rebuilding mid-edit
  (that happened on day one).  `MACOS_CLI_TREE=.` uses this tree, tagged `+dirty`.
- Model size 256 (`--model`); `--deadline` 90 s per command (the manual's `SegDemo.in`
  is the known slow loader).  A `timeout_load` is not always the engine: the three
  `pupilsim_redo_lens_*` decks in the 10-09 corpus record were a 138-character path
  cut at the CLI's 130-character command line (fixed b90887e; `cmd/LONG_PATH.jou`).
- The suite runs the CLI only — never MATLAB — so it never blocks a mex relink.
- A deck that fails is a finding, not a test defect, until shown otherwise; the
  reproducer goes into `ZGD_test_files/` and the item into PLAN §0.
