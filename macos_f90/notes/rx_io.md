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
  prints the header message once per key per load (list reset at the top of msmacosio.inc)
  -- warns, never refuses; a "key" not starting with a letter is skipped (a UTF-8 BOM alone
  on line 1: six FFSegDemo decks).  Both messages are in Manual §4 "Editing .In-Files".
  Gate `cli_tests/cmd/Rx_header_unknown.jou`.  **Prefix trap that
  the warning cannot see:** a short `LCMP` match CAPTURES look-alikes -- `OptFOVWt` is
  matched on 5 chars, so `OptFooKey= 1` is read as FOV weights.  When adding a header key,
  check that no shorter prefix above it already swallows it (and vice versa).
  **Since 2026-10-12 every active header match is the FULL key name** (23 lengthened, the
  table is `cli_tests/records/hdr_prefix_audit.md`; no deck in the 1353-deck extended corpus nor
  any binding writer used a shorter spelling), so `OptFooKey=` is reported, not read as
  `OptFOVWt`.  New header keys: match the full name.  Element-chain matches were not audited.
- **`OptAsph=` above its element's `VarElt=`** is reported and ignored (Manual Table 13; the
  parse keys on `isVarElt`, which only `VarElt=` sets; `OptZern=` sets it itself).  Gates `cli_tests/cmd/OptAsph_order.jou`,
  `Rx_header_unknown.jou` (+ `OptFooKey=`), pymacos `test_rx_parse_warnings.py` (capfd).

## Deferred PolyObsVec/Poly3DObsVec projection (msmacosio.inc + iosub.inc)
> User-facing: Manual §4 "Obscuration types and the ObsVec= line" (the type table, the three polygon forms, the Neg rule; 2026-10-17).  Gate `cli_tests/cmd/OBS_POLYGON.jou`.

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
User rule (Manual §4 MODIFY): no spaces around `=`.  Mechanism: CACCEPT in MOD_LOH consumes
only the first whitespace-separated token, so `mod ngridpts = 33` left VALUE empty and the bare
`READ(VALUE,*) nGridpts` (no ERR=/END=) hit EOF -- `forrtl: severe (24)`.  The guard after
GET_EQ (`LEN_TRIM(VALUE) == 0`, not HELP) prints one message and drains the buffered tokens
(`read_len(pstack)=0, read_cur(pstack)=0`) so the leftover `=` and `33` do not cascade into
bogus commands.  Deliberately narrower than a multi-token parser.  (The cosmetic blank
`WRITE(*,*)` in MOD_LOH's QUIT branch was dropped at the same time.)

## Rx comments: block, whole-line and in-line (2026-09-29, Scott's report)
The three forms and their fates as the USER sees them are Manual §4 "Editing .In-Files"
(comments bullet) and the VALIDATE paragraph.  Engine record:
- Block comments (`/* ... */`, `CommentBegin ... CommentEnd`) are a FIRST-TOKEN rule in
  `GET_EQ` (iosub.inc ~:2670): label 20 is GET_EQ's own read-next-line loop, so a block is
  consumed within one call and the top-of-call `LInCommentMode=.FALSE.` never sees it.  The
  parser always handled them; what 223a6ff fixed was the Phase-1 VALIDATOR (knew only `%`/`!`,
  read a `/*` line as a continuation row -> "blank line inside multi-row block", eac5mono.in
  149/151) and SAVE (capture saw only `%` lines, so every block was LOST).  Now the validator
  carries the same first-token mode and `RxCommentCaptureRaw` keeps opening, interior (incl.
  `%` lines and blanks) and closing lines verbatim.
- A one-line `/* x */` is unsupported by construction (`*/` is never a first token; the parser
  would swallow the rest of the file) -- the validator refuses it as "comment block never
  closed", as it does a genuinely unterminated block.
- Whole-line `%`: `RxCommentCapture` (PLAN sec.0 item 3), furniture-filtered, capped at
  `mRxComment` (500) with the overflow counted and reported ONCE at SAVE.
- In-line `Key= value  % note`: the `%` ends the value (iosub.inc ~:2607); by the recorded
  decision NOT preserved by SAVE (the fragile case behind the GridFile-tab bug).  The validator
  judges `Key= % note` as the EMPTY value the parser will see.
- `!` is a comment char to the VALIDATOR only -- the parser has no `!` branch.  Pre-existing.
- CLI, mmacos and pymacos share GET_EQ and the validator.  Gates: `ZGD_test_files/
  tst_block_comment.in` (CLI, both compilers, SAVE->load->SAVE byte-identical, gfortran==ifx),
  `mmacos/tests/tRxBlockComment` (live-keyword negative control).  eac5mono.in's only round-trip
  diff is the pre-existing 1-ulp psiElt print wobble.

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
(Manual §4, aspheric surfaces; guards warn, not error).  Applied at the four
AsphCoef read sites (first line, continuation groups, the incomplete group)
and to `AnaCoef=`.  The other multi-value keywords (MonCoef, ZernCoef, ...)
still use bare `READ(VALUE,*)` -- same hazard, not yet guarded; extend
`ReadRealsPad` there when one bites.

## Prescription validator (validate_prescription.F90)
What it refuses and accepts, for the user: Manual §4 VALIDATE.  Wiring: Phase-1 pre-validator
`validate_prescription_mod%ValidatePrescription(filename, ios, msg)` runs before MBFile6 opens
the .in file -- pure character I/O; ios=0 on success, else a bare `msg` (the caller adds the
file name).  Per-key empty-value exceptions are the `KeyEq` `.OR.` chain (`EltName`,
`BaseUnits`, `WaveUnits`; add there if Norbert flags more).  Wired into MBFile6 in macosio.F
(interactive: re-prompt via `GO TO 43`) and smacosio.F (SMACOS: `LOAD_SUCCESS=.FALSE.`, `GO TO
99`), only when iopt=1 (an existing file).  The top-level `VALIDATE` command
(macos_cmd_loop.inc) is the same subroutine with no state change.  Phase 2 (deferred):
IOSTAT-armoring of the msmacosio.inc parser internals; enum-value validation against the type
tables.  Line endings: CR / CRLF decks load from a temp LF copy made at OPEN (TO 2026-10-09,
item 5).

## Test prescriptions (ZGD_test_files/)
- All FreeForm test prescriptions use nGridpts=11 (gives 89 rays for Circular grid).
  Reason: model 256 has bRay=500 (BUILD limit) and mDrawRay=101 (nDrawElt array size);
  nGridpts=11 gives 89 rays which is under both limits.
  Grid-using prescriptions (tst_FF_g/mg/fg/refl) also need model 256 for mGridMat=256.
- Obscured rays now increment nBadRays so WARN fires: tracesub.F and propsub.F both
  have `nBadRays=nBadRays+1` inside the `IF (.NOT.L1(iRay))` / `IF (.NOT.LRayPass(iRay))`
  blocks that call SetRayFail(RayStat_Obscured).


## Three loader rulings of 2026-10-18 (PLAN_CONSOLIDATION 1.2 (d), (u) and section 4 item 1)
- **A numeric `PropType=` resolves to its name at load; `0` IS Geometric** (Dave).  The element chain
  (msmacosio.inc) matched the names and the numerals 1..15 but had no `0` branch: a legacy
  `PropType= 0` (old_Rx/KeckDemo.in, mars.in) left the allocation value 0 with the given-flag set, so
  ChkDf2 never defaulted it, the trace saw "not 1" (not Geometric) and SAVE printed `PropTypeName(0)`
  (bytes past the table on ifx, blanks on gfortran).  Now `0` -> 1, an unrecognized value is NAMED
  (`** PropType= Bogus at element 8 not recognized: the default applies`) and left to ChkDf2, and the
  SAVE writer guards the index (Geometric for anything outside 1..15).  Note the chain's 1-char
  numeral matches: `99` reads as 9 (GeomUpdate) -- a longstanding prefix quirk, not changed.
  Fixture `ZGD_test_files/tst_proptype_numeric.in` (0..15 in element order); gate `PROPTYPE_numeric.jou`.
- **A deck with no output-coordinate block loads** (it failed with no message: the label-61 line read
  was `END=99`).  Now `END=69`: end of file (or an `EOF` line) inside the element chain closes the last
  element as `nOutCord=` would (`LOutCordAbsent`), the two `70` blocks (macosio.F / smacosio.F) skip
  the nOutCord/Tout reads (`GO TO 72`), and ChkDf3's Tout default is `OutCordCandidate` (iosub.inc):
  the BEAM frame at the last element -- z = unit(VptElt(nElt) - VptElt(nElt-1)) (the source
  ChfRayDir for one element or coincident vertices; NOT psiElt, whose sign runs against the beam on
  a reflector and along it on a FocalPlane -- Cassegrain's FP has psi = +z with the beam +z), x from
  xObs or the global axis least aligned with z, y = z x x.  `LOutCordCandidate` makes the CLI ask
  (`OutCordQuery`, macosio.F, under MACOS_CMD in macos_cmd_loop.inc: YES keeps, NO reads nOutCord and
  the rows with the proposal as defaults); SMACOS/bindings take it and print the same lines.  The
  old cyclic default (`-psiElt(3), -psiElt(1), -psiElt(2)`) was orthonormal only for an axis-aligned
  element and no corpus deck reached it (all 1464 carry `Tout=`).  The load gate answers the prompt
  (`cli_load_gate.load`).  Fixture `tst_no_outcord.in`; gate `OUTCORD_missing.jou`.
- **SAVE writes the CALIB setup** (Dave: "SAVE the CALIB settings and have them transfer, loudly").
  `PrtCalibHeader` (source block, called from PrtSourceInfo before the comment emit) and
  `PrtCalibElt(iElt)` (element block) in iosub.inc write every Table 13 key whose value differs from
  its load default; the reader's counting rules shape the output: `ChfRayDir/ChfRayPos` and `Wavelen`
  are field 1 / wavelength 1, so fields 2..n go out as `OptChfRayDir=/OptChfRayPos=` pairs and
  wavelengths 2..n as the existing `ArrWaveLen=`; `OptFOVWt=` after the pairs; `OptTgtElt=` only when
  not nElt-1 (LOAD resolves the -1 sentinel); `OptRayGrid=` only when opt_npts differs from
  nGridpts/2-1 (the given-flag is a local of the two reader hosts); `VarCons=`/`OptAsph=` after their
  `VarElt=`; beam keys at the element they name, `OptBeamWt=`/`OptBeamCentroid=` with the first beam
  key.  `CalibInForceNote` (macos_cmd_loop.inc, after ModesInForceNote, CLI and bindings) prints one
  line when a load reads a setup.  Gate `CALIB_save_roundtrip.jou` (tst_calib_mixed,
  opt_example_constrained: saved keys, the notice, re-save byte-identical); the load gate's round trip
  now covers every CALIB fixture's keys.  `SAVE_KEYWORD_AUDIT.md`'s Opt* family is closed.
