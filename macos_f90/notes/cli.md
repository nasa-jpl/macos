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

