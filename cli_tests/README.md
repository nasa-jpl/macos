# cli_tests — the CLI regression suite (PLAN_CONSOLIDATION §2)

The interactive `macos` CLI shares its parsers and command loop with mmacos and
pymacos but had no regression suite of its own; the bindings have ~570 fast
tests and pymacos ~6600.  This directory is that suite.  It drives the CLI through
a pty (readline needs a tty; the model-size prompt comes first) and writes RECORDS
— the record is the result, the script prints the summary.

```
./run_cli_tests.sh load        # the core corpus (102 decks), both compilers  (~minutes)
./run_cli_tests.sh corpus      # every committed prescription (~1270)         (~an hour)
./run_cli_tests.sh extended    # + untracked sandbox decks (informational)
./make_corpus.sh               # regenerate the three deck lists
```

## Layers

| layer | what | status |
|---|---|---|
| `load` (§2.1) | per deck: LOAD → `opd nElt` → SAVE; new process: LOAD the save → SAVE; byte identity of the second save against the first | built 2026-10-09 |
| `cmd` (§2.2) | one pty-driven journal per command-reference entry, `# EXPECT:` lines asserted, a must-fail twin where the command has a known failure mode | next |
| `ab` (§2.3) | the previous tagged binary in `build_ab/` against this one on the `corpus` record; any changed deck named in the commit | next |

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
- `must_fail/` — the non-vacuity decks (§2.4).  `short_zerncoef.in`: a `ZernCoef=`
  block one value short — a bare `READ(VALUE,*)` abort that KILLS the CLI (PLAN §0
  line 270; TO's item 1); when that fix lands this deck will load and must move to
  the positive corpus with its padded-and-warned expectation.  `missing_value.in`:
  `nGridpts=` with no value — the validator refuses it; this one is permanent.
- `records/` — `<kind>_<compiler>_<sha>[+dirty].csv` + `.log`, `mustfail_*.csv`,
  `<kind>_<compiler>_latest.txt` (the last summary).  Committed: the records of
  tagged runs; the `.log` files are gitignored.
- `ab_corpus_seed.py` — the FwdRoot corpus A/B driver (2026-10-04) this grew from.

## Conventions

- Both compilers, always; gfortran's stricter runtime has found real bugs.  The
  binaries come from a CLEAN worktree at a committed SHA (`MACOS_CLI_TREE`, default
  `~/dev/macos_cli_base`; `git worktree add ../macos_cli_base <sha>`, then
  `makems.sh release` and `makems.sh release gfortran` there) — never from the
  shared working tree, whose build dirs another lane may be rebuilding mid-edit
  (that happened on day one).  `MACOS_CLI_TREE=.` uses this tree, tagged `+dirty`.
- Model size 256 (`--model`); `--deadline` 90 s per command (the manual's `SegDemo.in`
  is the known slow loader).
- The suite runs the CLI only — never MATLAB — so it never blocks a mex relink.
- A deck that fails is a finding, not a test defect, until shown otherwise; the
  reproducer goes into `ZGD_test_files/` and the item into PLAN §0.
