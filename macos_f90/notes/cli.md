# MACOS engine notes — The interactive CLI: readline sub-prompts, prompt loops under SMACOS, LOG

> **Read this file before touching:** `macosio.F`, `mhist.c`, `macos_cmd_loop.inc` prompts.  These notes are the engine record
> (what moved, the measurements, the gates); they are NOT auto-injected — the index is
> `macos_f90/CLAUDE.md`, the root rules `../../CLAUDE.md`.  Text is moved verbatim from the
> former monolithic cheatsheet (2026-10-15 split); add new entries here, not to the index.

## CLI sub-prompts via mhist cache (macosio.F + mhist.c)
The old fragile pattern was: Fortran ACCEPT routine WRITE-s a prompt
(no trailing newline), then calls READ_LOH which goes through MHIST →
`readline(empty)` for the input. readline couldn't see the prompt and
mis-managed cursor state, producing the recurring "MACOS> the new
element data?" overwrite class of bugs.

Replaced with a small cache in mhist.c:
- `sub_prompt_buf[256]` (file-scope, OUTSIDE the READLINE_LIBRARY
  guard so non-readline builds still link)
- `set_sub_prompt_(s, slen)` / `clear_sub_prompt_()` — Fortran-callable
- `mhist_` for sub-prompts (mp[0]==' ') passes `sub_prompt_buf` to
  readline, so readline renders the prompt itself and owns cursor
  state authoritatively. Setter appends one trailing space to match
  the legacy `' ',A,' [',A,']: '` format the old WRITE produced.

CACCEPT, IACCEPT, DACCEPT, RACCEPT in macosio.F (NOT in smacosio.F —
that one has no readline) now:
1. Build the full prompt string from TEXT + `[default]:`.
2. If `read_len(pstack) > 0` (token will come from the existing
   buffer, no MHIST call), WRITE the prompt directly with format
   `(A,$)` — matches legacy behavior for multi-token lines like
   `stop elt 1 0,0`.
3. Else `CALL set_sub_prompt(...)` then `CALL READ_LOH(...)` and
   `CALL clear_sub_prompt` afterwards.

Use `LEN_TRIM(prompt_str)` for the slice length, NOT `ICLEN` —
`ICLEN` returns the length of the first *non-blank* prefix (length
of the first word), which gives 4 for ` Enter...` etc. That bug ate
several debug cycles; don't repeat.

When adding new ACCEPT-style prompt routines, follow the same
pattern. New top-level commands that need their own prompt should
either reuse the existing ACCEPT routines or push their prompt
through `set_sub_prompt` directly.

Non-readline builds (Windows always; any clone whose
`readline-8.2/libreadline.a` isn't built) do NOT render readline
prompts: PARSE_LOH's `#else` branch is a bare `READ(*)`, so the
cached sub-prompt must be printed explicitly — `print_sub_prompt_`
(mhist.c, outside the guard) called before that READ (e81b235).
Without it, `pert` <ENTER> blocks on invisible input ("pert freezes",
Luis 2026-08-06); full one-line commands still work, which disguises
it. cmake now auto-builds the bundled readline when the .a is
missing; `-DMACOS_FORCE_NO_READLINE=ON` is the test hook for the
fallback path.

## IACCEPT_S reprompt-loop sweep (macos_cmd_loop.inc, 2026-07-16)
Interactive commands that validated an element/ray id and jumped BACK
to their ACCEPT prompt on failure (`GO TO <prompt label>`) hang or
misbehave under SMACOS batch: IACCEPT_S reads the exhausted stack and
returns the DEFAULT, so a sticky default (`IACCEPT_S(iElt, iElt, ...)`
— the OPD case) spins forever, and a fixed default (1, iCurWFElt,
nElt−1) silently substitutes the wrong element.  ALL such validation
reprompts now abort with `'<CMD>: command aborted.'` + `GO TO 1`:
BUILD, ORS, SRS (both prompts), FEXIT, SXP, XPS, SPOT, MVAR, DVAR,
GPERTURB, LPERTURB, and AVAR — whose old SMACOS branch called `stop`,
killing the host process when the engine is a mex/library.
Deliberately KEPT as loops: CTRACE / SEGRAYTRACE "0=quit" multi-entry
prompts (default 0 self-terminates on an exhausted stack), the
`#ifdef MACOS_CMD` MRESET size prompt (interactive-only), and the
Zernike-range prompt (sets valid defaults before looping, so it
self-heals).  When adding a NEW command prompt, never re-enter the
ACCEPT on validation failure — abort to label 1.

## LOG command cleanup (macos_cmd_loop.inc)
- The `LOG` interactive command (log10-intensity wavefront display, not
  transcript logging -- which is `JOU`/`JOURNAL`) previously dumped an
  always-on debug print (`*****B4 SrfOut: mdttl=`) plus a 5 MB ASCII file
  (`IntLog.txt`) to the cwd on every invocation, gated by an `#if 1`.
  Both removed; diagnostic block gated `#if 0`.  `LOG` syntax and 3-char
  min-match unchanged -- existing `.jou` scripts using `log <iElt>` keep
  working.


## `ifRayPrint` is global state (B.2, 2026-10-17)

`ifRayPrint` (macos_mod) decides whether CTRACE prints each ray's segments.  CTRACE sets it
`.FALSE.` on return (tracesub.F), so a handler that traces more than once in one command must set
it before each CTRACE (SEGRAYTRACE set it once before its loop; a second segment traced silently).
A helper must restore the caller's value, never force one: FSR (utilsub.F; RRAY / PRAY / FSR) ended
with `ifRayPrint=.TRUE.`, and the next OPD printed every ray of the grid (tens of thousands of lines
on Rx_Cass_FarField) -- it now saves and restores.  TPR clears it itself, so a journal leg for this
puts the OPD before any TPR.  FindMarginalRay and getRefRays (collim.inc) still set `.TRUE.` on
return; EFL's caller resets it (measured: no print after EFL), CALIB's was not measured.

## The `!` shell line and a leading `cd` (found 2026-10-18)
`! cd $TMP && python3 tool.py ...` does NOT run the tool: the `!` handler
treats a line that starts with `cd` as the CLI's own CD command (chdir, then
"Current Folder:") and drops the rest of the line.  In journals give the
tool absolute paths (or a `--dir` option) and use one `mv` per `!` line.

## The startup model-size prompt (PLAN_CONSOLIDATION 1.2 (yy), 2026-10-20)
`set_macos_problem_size` (macos_mod.F, `#ifdef CMACOS`) asks before readline or any LOG is up:
a plain `read(*,'(A)')`, not the ACCEPT routines, so none of the sub-prompt machinery applies.
It was a list-directed READ of the answer with no `err=`: a word, `256abc` or a number too long
for an integer killed the CLI (forrtl 59 / 64; gfortran "Bad integer" / "Integer overflow");
`256.5` crashed gfortran and read as 256 under ifx; `-1` was the internal no-answer sentinel and
gave the default; `512 x` gave 512.  `problem_size` was not reset before the retry and the buffer
was not cleared, so a refused size then a blank line or end of input re-read the stale answer
forever (2.3 M lines in 5 s on a pipe).  Now: blank = the default, otherwise the trimmed answer
must be at most five digits and one of the sizes; anything else is refused by name, at most
`mTries` = 5 times; end of input and the fifth refusal `STOP 1` with a named line.  MRESET is
not affected (IACCEPT_S re-asks on a non-integer).  The api path never reaches this prompt:
`macos_init_all` (param_mod) still `stop`s on a bad size, but pymacos `init` and mmacos
`macos.init` refuse before the call (pymacos's list omits 8192).  Gate
`cli_tests/cmd/STARTUP_size.jou` (the gate's `# MODEL:` takes `a | b | c`, answers in turn);
red on the pre-fix binary (dies before MACOS>).

## DOPD, LSPOT, RESET (PLAN_CONSOLIDATION 1.2 (k) (m) (i), 2026-10-22)
- **DOPD of identical OPDs** read `RMS dOPD = -1.0`: `RMS2` (mathsub.F) returns -1 when NO cell
  of the map is non-zero -- its sentinel -- and two identical OPD runs give an all-zero
  difference.  The finding blamed ZABERR / ZRM; measured, `opd; opd; dopd` did it too.  DOPD now
  prints 0 and "(the two OPDs are identical)".  Not changed: `RMS2` takes its mean and RMS over the
  NON-ZERO cells only, so an exactly-zero valid cell (e.g. the chief-ray reference) is left out --
  the "zero is not no ray" rule; changing it would move every DOPD value.
- **LSPOT** honours the NONE plot type (`iPloTyp == 0`) as SPOT does; it called SpotDiag
  regardless (eight `%giza ... No device open` lines).  Under a file plot type (5/6/8/9) it still
  calls SpotDiag, not SPOTOUT as SPOT does -- not changed.
- **RESET** sets the plot type back to its start-up value, GRAY (`iPloTyp = 2`, ArchInit).
  Gates: DOPD_identical.jou, LINEAR_MODEL.jou (`EXPECT-NOT: No device open`), SESSION_INFO.jou.
