# NOTE to CCMac — your next pair: macos 1c867405 + resources b554792 (CC, 2026-10-20)

Both on origin dev-candidate.  Engine CHANGED since 9a971d25, so the full 7-check clean-clone gate; the binary
will not be byte-identical to gate 6's.

**What changed (Dave's rulings on macos-fd's four Zernike findings, PLAN_CONSOLIDATION 1.2 (tt)-(ww)):**
- (tt)/(vv), macos 7cd22b63 (CC): `ZernType= NormAnnularNoll` refused by name at load, MOD and the api setters;
  an ExtFringe coefficient in modes 47-49 / 58-64 refused by mode.  Journal `ZERN_refused.jou` (new).  The three
  SAVE-keys fixtures now carry `NormNoll` (0 movers on Linux).
- (uu), macos 548b0da7 (macos-fd): NormHex OPENED -- `ZernType_MaxMode(NormHex)=45`, `hzdef.inc` validated 45/45
  against a textbook Gram-Schmidt over the hexagon (`cli_tests/tools/hzdef_check.py`), and its 313 literals were
  REAL(4) (1.8e-6 off at mode 38), now double.  Journal `ZERN_normhex.jou` (new): SRT sag = Mahajan & Dai mode 10
  to every digit -- if arm64 differs in the last digits there, send the two numbers (a 1e-12-class pin may need
  the 5-figure rule).
- (ww), resources b80a431 (macos-fd): pymacos's `setEltSrfZern*` limits now read the generated constants
  (MZERNTYPE / MZERNMODES); new `test_zern_veneer_limits.py` (4 tests).  `test_zernike_orderings.py` (b554792)
  now expects the two refusals and checks NormHex: 20 tests.
- Manual 58d6759a: Appendix E regenerated, Section 4 / cmdref / mmacos help say "refused at load".

**Expected on arm64:** cmd 109/111 + 2 must-fail (COORD_ACOOR = (ff) the only miss; 110 journals + ZERN_refused
+ ZERN_normhex = 111 with must-fails counted separately as before), mmacos fast 656/0 + the filtered leg, pymacos
6941 + 6 (zernike_orderings 14 -> 20) + 4 (veneer limits) = 6951 / 0 / 2 skips.  Linux gates: CC at 7cd22b63 and
macos-fd at 548b0da7, both PASS, both compilers, records in `cli_tests/records/`.

**Open for Dave, not for you:** (xx) a MOD of a coefficient beyond a type's mode count is accepted silently (three
faces); `ZCOZernType= NormAnnular` (output fitting) is not refused.  (ff) unchanged, TO's.
