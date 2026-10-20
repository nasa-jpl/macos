# NOTE to CCMac — your next pair: macos 3ead7473 + resources 00697be (CC, 2026-10-20, rewritten; supersedes the 1c867405 pair)

Both on origin dev-candidate.  Engine CHANGED since 9a971d25 (gate 6): the full 7-check clean-clone gate.

**What changed, all ruled by Dave 2026-10-20 (PLAN_CONSOLIDATION.md 1.2 (tt)-(yy), section 5 for the Linux gates):**
- (tt)/(vv) macos 7cd22b63 (CC): `ZernType= NormAnnularNoll` refused by name at load, MOD and the api setters; an
  ExtFringe coefficient in modes 47-49 / 58-64 refused by mode.  6fda7711: the OUTPUT fit too -- `ZCOZernType=
  NormAnnular` refused at load and SZCO's 'Ann' answer refused (zern_fit's annular basis is the same ratio-blind
  converter).  Journal `ZERN_refused.jou` (new, 6 legs).  The three SAVE-keys fixtures now carry `ZernType= NormNoll`
  and `ZCOZernType= NormHex` (0 movers on Linux).
- (uu) macos 548b0da7 (macos-fd): NormHex OPENED -- `ZernType_MaxMode(NormHex)=45`, `hzdef.inc` validated 45/45
  against a textbook Gram-Schmidt over the hexagon (`cli_tests/tools/hzdef_check.py`); its 313 literals were REAL(4)
  (1.8e-6 off at mode 38), now double.  Journal `ZERN_normhex.jou` (new): SRT sag = Mahajan & Dai mode 10 to every
  digit -- if arm64 differs in the last digits, send the two numbers (a 1e-12-class pin may need the 5-figure rule).
- (xx) macos 8f43e952 (CC): a coefficient past a type's mode count refused by mode at every face (ChkDf2, MOD --
  whose index bound was a 'ToDo' no-op -- and the FF/Mon sparse mode lists bounded by their arrays at load).
  Journal `ZERN_mode_count.jou` (new).
- (yy) macos 49c87e13 (macos-fd): the CLI startup model-size prompt refuses an answer it does not recognize and asks
  again (at most 5 times; end of input stops with exit status 1) -- it crashed before.  Journal `STARTUP_size.jou`
  (new; `cli_cmd_gate.py` gained `# MODEL: a | b | c`).  On macOS check the shell leg's exit status reads 1 there too.
- (ww) resources b80a431 (macos-fd): pymacos's `setEltSrfZern*` limits read the generated constants; new
  `test_zern_veneer_limits.py` (4).  `test_zernike_orderings.py` (b554792) expects the refusals and checks NormHex.
- Manual 58d6759a / 6dc22c0e: Appendix E regenerated, Section 3.3 / 4, cmdref, mmacos help.

**Expected on arm64:** cmd 112/114 + 2 must-fail (COORD_ACOOR = (ff) the only miss; 110 + ZERN_refused, ZERN_normhex,
ZERN_mode_count, STARTUP_size = 114 with the must-fails counted as before), mmacos fast 656/0 + the filtered leg,
pymacos 6941 + 6 + 4 = 6951 / 0 / 2 skips.  Linux: CC's gates at 7cd22b63, 39275005, 3ead7473 (806bf5e4) and TO's at
548b0da7, 49c87e13 all PASS, both compilers, records in `cli_tests/records/`.

**(ff)** unchanged, TO's.  Nothing else open for you.
