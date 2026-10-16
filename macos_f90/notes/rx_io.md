# MACOS engine notes — Prescription I/O: parser chains, comments, keys, defaults, validator, trace-must-not-write-Rx

> **Read this file before touching:** `msmacosio.inc`, `iosub.inc`, `macosio.F`, `smacosio.F`, `validate_prescription.F90`, SAVE.  These notes are the engine record
> (what moved, the measurements, the gates); they are NOT auto-injected — the index is
> `macos_f90/CLAUDE.md`, the root rules `../../CLAUDE.md`.  Text is moved verbatim from the
> former monolithic cheatsheet (2026-10-15 split); add new entries here, not to the index.

## A trace must not rewrite the prescription: the IndRef stamps (2026-10-09/10)
CTRACE's Return (PLAN item 3) and Obscuring (1.1b 3b) branches wrote the RUNNING medium into the
element's own `IndRef`.  Nothing reads either element's `IndRef` as a medium (ReturnSrf / ObsSrf use
`na`; CPROPAGATE resets its medium only after a Refractor or LensArray), so the stamps only mutated the
Rx: an obscuration inside glass reported / SAVEd the glass index after any trace (ZGD
`tst_obs_in_glass.in`, gate `cli_tests/cmd/IndRef_stamp.jou`).  Both removed.  The REFRACTOR branch
never had a stamp: its `nb` (the element's IndRef) is never assigned; a variable-index element
delivers its per-ray index through `nc`/`ifnc` -> CTRACE's `NewIndRef`.  Measured: SAVE without a
trace vs SAVE after `opd nElt` on 403 corpus decks -- `IndRef` 0 changed.  The same measurement found
two OTHER trace-time mutations, CLOSED 2026-10-12: `xObs` recomputed (106 decks) -- the obscuration
frame is now built ONCE at load by `OrthoObsFrameLoad` (math_mod; end of ChkDf2: an xObs within the
dead band of the perpendicular is KEPT, so SAVE -> load -> SAVE stays a fixed point -- rebuilding it
unconditionally turned e2e_pie's round trip 'ulp') and the traces' `OrthoObsFrame` keeps it unless it
moves by more than round-off; and `Surface= Interpolated` SAVEd as `Conic` (4 legacy decks)
-- an Interpolated surface with no data traces as its base conic through `SrfTypeTr` (elt_mod, set
per trace, passed at the 22 dispatch call sites); `SrfType` is never written by a trace.  Re-run:
`cli_tests/records/trace_rewrites.md` (0 and 0); gate `cli_tests/cmd/TRACE_keeps_rx.jou`.  **Rule:
a trace reads the Rx and writes TRACE state; if it needs a modified copy of an Rx value, keep the
copy (a `*Tr` array, a local), never overwrite the Rx array.**

## Prescription I/O notes (msmacosio.inc / iosub.inc)
- The parser has two chains: label 50 (header/global) and label 61 (per-element).
  FreeForm keywords must be in BOTH chains. The element chain (61) now includes:
  lData, lFF, pFF, xFF, yFF, zFF, FFZernType, nFFZernCoef, FFZernModes, FFZernCoef,
  MonZernType, nMonZernCoef, MonZernModes, MonZernCoef.
- nFFZernCoef/FFZernModes/nMonZernCoef/MonZernModes reset per-element at loop start
  (before LOHIN2 call).
- Output (PrtSingleEltInfo in iosub.inc): pData/xData/yData/zData written for
  FreeForm only when nGridMat > 0.
- ChkDf2: FreeForm with nGridMat=0 silently satisfies pData/xData/yData/zData/lData
  (no spurious "Default used" warnings).
- **An unknown HEADER key is reported (2026-10-11, PLAN 3.2).**  The header chain's
  final `ELSE` (label-50 chain, just above label 60) was a bare `GO TO 50`: an unknown or
  misspelled header key vanished with NO message, while the element chain has always
  printed `MBFile6: Unidentified string`.  That is how `OptMxItrs=` (9 chars against a
  10-char `OptMaxItrs` match) silently lost the user's cap.  Now `HdrUnkNote` (elt_mod)
  prints `MBFile6: Unidentified header string = <key>` once per key per load (list reset
  at the top of msmacosio.inc) -- warns, never refuses; a "key" not starting with a letter
  is skipped (a UTF-8 BOM alone on line 1: six FFSegDemo decks).  Corpus, engine truth:
  see the PLAN 3.2 line.  Gate `cli_tests/cmd/Rx_header_unknown.jou`.  **Prefix trap that
  the warning cannot see:** a short `LCMP` match CAPTURES look-alikes -- `OptFOVWt` is
  matched on 5 chars, so `OptFooKey= 1` is read as FOV weights.  When adding a header key,
  check that no shorter prefix above it already swallows it (and vice versa).
  **Since 2026-10-12 every active header match is the FULL key name** (23 lengthened, the
  table is `cli_tests/records/hdr_prefix_audit.md`; no deck in the 1353-deck extended corpus nor
  any binding writer used a shorter spelling), so `OptFooKey=` is reported, not read as
  `OptFOVWt`.  New header keys: match the full name.  Element-chain matches were not audited.
- **`OptAsph=` above its element's `VarElt=`** prints `** OptAsph= at element N IGNORED:
  OptAsph= must follow VarElt=` (the parse keys on `isVarElt`, which only `VarElt=` sets;
  `OptZern=` sets it itself).  Gates `cli_tests/cmd/OptAsph_order.jou`,
  `Rx_header_unknown.jou` (+ `OptFooKey=`), pymacos `test_rx_parse_warnings.py` (capfd).

## Deferred PolyObsVec/Poly3DObsVec projection (msmacosio.inc + iosub.inc)
- For polygon obscuration, the 3D-vertex forms (`PolyObsVec=`/`Poly3DObsVec=`)
  used to require xObs to be parsed BEFORE the polygon line, because
  `SetCvxPolyObsVtx` projects the 3D vertices to element-local 2D at parse
  time and needs xObs/yObs to do so. Wrong order produced a hard error.
- Now the order doesn't matter. On a 3D-vertex polygon line:
  - If `xObs_FLG` is true, project immediately (legacy behaviour).
  - Else, stash raw 3D vertices into `PolyObsVtx3DPending(3, mPolySide, mObs,
    mElt)` and the count into `nPolyObsVtx3DPending(mObs, mElt)` (both in
    elt_mod), defer the projection.
- xObs parser, after setting xObs_FLG, walks the pending-list for the
  current element and runs `SetCvxPolyObsVtx + SetCvxPolyObsBound` on any
  stashed vertices.
- ChkDf2 (iosub.inc): after the existing psiElt-based xObs default kicks in
  for elements that need it, runs the same pending-resolution pass. If
  xObs is still unavailable for an element with pending vertices, ELEM_OK
  is set false with a targeted message naming iElt and iObs.
- Error message text now references both accepted keyword names —
  `PolyObsVec/Poly3DObsVec` — not the never-implemented `PolyObs3DVec`.
- Note: `PolyApVec=`/`Poly3DApVec=` have the same xObs-ordering wrinkle
  but were left as-is (still error on bad order). Same deferral pattern
  would extend to them if needed.

## MOD command: empty-value crash guard (macosio.F MOD_LOH)
- `mod ngridpts = 33` (spaces around `=`) used to crash with
  `forrtl: severe (24): end-of-file during read, unit -5`. CACCEPT in
  MOD_LOH only consumes the first whitespace-separated token, so VALUE
  ends up empty and the subsequent `READ(VALUE,*) nGridpts` (which has
  no ERR=/END= label) hits EOF on the empty string.
- Guard added right after GET_EQ: if `LEN_TRIM(VALUE) == 0` and the
  command isn't HELP, print a clear warning ("assignment must be
  \"<var>=<value>\" with no spaces around \"=\"") and drain remaining
  buffered tokens (`read_len(pstack)=0, read_cur(pstack)=0`) so the
  leftover `=` and `33` tokens don't cascade into more bogus commands.
- Smaller scope than fixing the parser to slurp multi-token assignments.
  Users get one targeted message and prompt continues. `mod ngridpts=33`
  (no spaces) still works as before.
- Also dropped the cosmetic `WRITE(*,*)` blank line in the QUIT branch
  of MOD_LOH — left a stray empty line between the MOD prompt and
  `MACOS>` on every exit.

## Rx comments: block, whole-line and in-line (2026-09-29, Scott's report)
Three comment forms, three different fates -- know which is which before
touching any of them:
- **Block comments `/* ... */` and `CommentBegin ... CommentEnd`** are a
  FIRST-TOKEN rule in `GET_EQ` (iosub.inc ~:2670): the marker must start
  its line; label 20 is GET_EQ's own read-next-line loop, so a block is
  consumed within one call and the top-of-call `LInCommentMode=.FALSE.`
  never sees it.  **The parser has always handled them.**  What was broken
  (both fixed 223a6ff): the Phase-1 VALIDATOR ran first and knew only
  `%`/`!`, so a `/*` line (no `=`) was read as a continuation row of the
  last multi-row key and the deck refused ("blank line inside multi-row
  block" -- eac5mono.in at 149 for a `/*` at 151); and SAVE LOST every
  block, because capture saw only `%` lines.  Now the validator carries
  the same first-token mode, and `RxCommentCaptureRaw` keeps the opening,
  interior (incl. `%` lines and blanks) and closing lines verbatim so the
  SAVEd deck re-enters the block on reload.
- **A one-line `/* x */` is NOT supported** (first-token rule: `*/` is
  never a first token) -- the parser would swallow the rest of the file.
  The validator now refuses it as "comment block never closed" rather than
  let that happen.  Same for a genuinely unterminated block.
- **Whole-line `%` comments** round-trip through SAVE (`RxCommentCapture`,
  PLAN sec.0 item 3), furniture-filtered; capture is capped at
  `mRxComment` (500) lines and overflow is now counted and reported ONCE on
  the console at SAVE, not silent.
- **In-line `Key= value  % note`** is READ (the `%` ends the value,
  iosub.inc ~:2607) and, by the recorded decision, NOT preserved by SAVE
  (the fragile case behind the GridFile-tab bug).  The validator now judges
  `Key= % note` as the EMPTY value the parser will see.
- `!` is a comment char to the VALIDATOR only -- the parser has no `!`
  branch and would tokenize a `!` line as an unknown keyword.  Pre-existing,
  left alone.
- The CLI, mmacos and pymacos share GET_EQ and the validator, so one fix
  covers all three.  Gates: `ZGD_test_files/tst_block_comment.in` (CLI, both
  compilers, SAVE->load->SAVE byte-identical, gfortran==ifx) and
  `mmacos/tests/tRxBlockComment` (the binding path, with a live-keyword
  negative control).  The only round-trip diff on eac5mono.in is a 1-ulp
  psiElt print wobble on a tilted unit vector, reproduced on the pre-fix
  binary -- pre-existing, not comments.

## `GlassElt=` was DEAD engine-wide until 2026-09-30 (catalog wiped per load)
The glass catalog (`GlassName`/`GlassTable`, ~200 Sellmeier rows from
`glass_builtin.f90`) is loaded ONCE at start-up (`rl_macos_glass.inc` in the
CLI's model-size reset, `smacos_glass.inc` at SMACOS first entry).  Every Rx
load runs `reinitialise_variables()` -> `elt_mod_init_vars()`, which blanked
BOTH catalog arrays before the parser's `GlassElt=` lookup (msmacosio.inc
~:2796), so every glass element kept its written `IndRef` and traced as AIR
-- CLI, mmacos and pymacos alike.  Nothing caught it: no Rx in the corpus
used a glass name (PLAN sec.0 item 4 said so, and "no user" was never read
as "untested").  Found by TO's dyson5 gates.  Fix (macos b000390): the two
catalog arrays are blanked once in `elt_mod_init` (allocation); per-element
`GlassElt`/`GlassCoef` stay Rx state and are still reset per load.  **Rule:
`elt_mod_init_vars` resets PRESCRIPTION state only; anything loaded once at
start-up (catalogs, tables) is blanked at allocation.**  Gate: CLI A/B on
`mmacos/tests/Rx/Rx_GlassPlate.in` (pre-fix: IndRef 1.0 and OPD 0 at every
wavelength; post-fix: Malitson coefficients SAVEd, OPD 6.861e-4 / 6.679e-4 /
6.491e-4 m at 0.5/1.0/2.0 um) + `tGlassDispersion` / `tGratingImmersed`
(mmacos, SUITE_FAST).  Related facts pinned by TO: the Grating branch takes
the EXIT index from the grating element as written (an immersed grating
carries its glass on its own element); `ray_info_get` returns the OUTGOING
direction at the element despite its "before surface" comment; CaF2 is in
the table since 0ca61c1.

## Short multi-value lines: `AsphCoef=` / `AnaCoef=` pad with zero (2026-10-01)
A `AsphCoef=` line with fewer values than `nAsphCoef` (default 4) was an
UNCAUGHT end-of-file in the list-directed internal READ -- a Fortran runtime
abort that kills the HOST when the engine is a mex / .so (TO, dyson5 beat 3;
same class as the `mod ngridpts = 33` and `stop obj 0 0 0` EOF crashes).
Now `ReadRealsPad` (elt_mod, shared by the CLI and SMACOS parsers) reads as
many values as the line carries, pads the rest with 0 and prints ONE line
naming the element and the count -- fewer coefficients than slots means the
higher orders are zero (guards warn, not error).  Applied at the four
AsphCoef read sites (first line, continuation groups, the incomplete group)
and to `AnaCoef=`.  The other multi-value keywords (MonCoef, ZernCoef, ...)
still use bare `READ(VALUE,*)` -- same hazard, not yet guarded; extend
`ReadRealsPad` there when one bites.

## Prescription validator (validate_prescription.F90)
- Phase-1 pre-validator: `validate_prescription_mod%ValidatePrescription
  (filename, ios, msg)` runs before MBFile6 opens the .in file. Pure character
  I/O: catches missing values after `Key=` (inline or continuation), bad keys,
  file-not-found, and open errors. Also catches blank lines inside multi-row
  continuation blocks (e.g. a stray blank in the middle of a `Tout=` matrix),
  which trip msmacosio.inc's fixed-format READ even when every individual
  line is well-formed. Returns ios=0 on success; ios/=0 with a
  bare `msg` (no filename prefix — caller adds context).
- Per-key empty-value exceptions (KeyEq helper inside the module): `EltName`,
  `BaseUnits`, `WaveUnits` may have empty values. Older prescriptions leave
  these blank and the parser tolerates it. Add more keys to the `.OR.` chain
  in `ValidatePrescription` if Norbert flags additional ones.
- Wired into MBFile6 in macosio.F (interactive: re-prompt on failure via
  `GO TO 43`) and smacosio.F (SMACOS: clean abort with LOAD_SUCCESS=.FALSE.,
  `GO TO 99`). Both paths run the validator only when iopt=1 (existing file
  being loaded).
- Top-level `VALIDATE` command added in macos_cmd_loop.inc — same subroutine,
  no state change, just prints `<file>: OK` or `<file>: <msg>`.
- Phase 2 (deferred): IOSTAT-armoring of msmacosio.inc parser internals;
  enum-value validation against the type tables.

## Test prescriptions (ZGD_test_files/)
- All FreeForm test prescriptions use nGridpts=11 (gives 89 rays for Circular grid).
  Reason: model 256 has bRay=500 (BUILD limit) and mDrawRay=101 (nDrawElt array size);
  nGridpts=11 gives 89 rays which is under both limits.
  Grid-using prescriptions (tst_FF_g/mg/fg/refl) also need model 256 for mGridMat=256.
- Obscured rays now increment nBadRays so WARN fires: tracesub.F and propsub.F both
  have `nBadRays=nBadRays+1` inside the `IF (.NOT.L1(iRay))` / `IF (.NOT.LRayPass(iRay))`
  blocks that call SetRayFail(RayStat_Obscured).

