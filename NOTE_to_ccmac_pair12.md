# NOTE to CCMac — your next pair: macos 08b92085 + resources 446bd45 (CC, 2026-10-23)

Both on origin dev-candidate (macos tip 668542ec carries 08b92085's records and plan line above it).  Engine CHANGED
(TO's bindings round): the full 7-check clean-clone gate.  Dave is away until ~2026-10-26; CC holds his standing
push word for gated work.

**What changed (TO; PLAN_CONSOLIDATION 1.4a/1.4b):**
- **SINT through the api** (`sint_cmd`), with three fixes in `SetIntSrf`: a header longer than the file kept one
  unread row as a data point (P-V 1.91 -> 4.13 um; ifx and gfortran differed) -- the rows read are used now, with a
  warning; the missing-file prompt looped forever under SMACOS (an empty stack re-parsed 'SINT') -- the bindings no
  longer hang, the SMACOS build does not prompt; the CLI looked for `<deck>.in.srfN` and SMACOS for `<deck>.srfN`
  -- each tries the other spelling second.  Journal `SINT_short_file.jou` (new).
- **Veneer defaults fixed on the engine-value rule** (an omitted argument = the ENGINE's value): mmacos `srs` link
  default false -> true (macos.srs had been a no-op since the SRS fix), `set_src_size` no longer clears the
  obscuration, `beam('cos')` no longer uses radius 1 silently, `opd_psf` model 512 -> 256; pymacos `setEltSrfZern`'s
  default call no longer raises; `tSrsBugFlatZ` un-vacated.  Two left for Dave: the elt_csys update flag (pymacos
  True vs engine/mmacos False) and the perturb frame (veneers LOCAL vs CLI GLOBAL).
- api coverage confirmed closed: 164/164 public routines bound, enforced by `test_api_coverage` / `tApiCoverage`.
- **Suite membership**: tSint, tVeneerDefaults, tMessages, tSaveNewFileQuiet, tRxLoadRulings1018 join the mmacos
  fast suite -- the last three were in no suite before (you ran them by name at gate 11).  cmdref regenerated.

**Linux (TO, both compilers):** load 160/160; cmd all + 2 must-fail; corpus 1468/1468; corpus A/B 492998db ->
08b92085 0 movers (records e7a0d6a2); pymacos full 7023 / 1 skip; mmacos fast 711/0.

**Expected on arm64:** cmd 119/119 + 2 must-fail (118 + SINT_short_file); load 165/165, 0 movers vs your gate-11
base; mmacos fast 711/0 (the five new classes included); pymacos 7023 / 0 / 1 skip.  The one line to read closely:
SINT's header-row fix is the compiler-sensitive one on Linux (ifx vs gfortran differed BEFORE the fix); if arm64
reads SINT_short_file differently, send the P-V numbers.

**Not for you:** the mmacos manual Part I is in progress in TO's lane (no engine change); the IODC pick waits for
Dave.  Your IODC record babb2c01 is merged and on origin -- thank you.
