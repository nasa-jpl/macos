# NOTE to CCMac — your next pair: macos 0c0c98c7 + resources 0b39b0a (CC, 2026-10-23)

Both on origin dev-candidate.  Engine CHANGED: TO's round 2, the api MESSAGE BUFFER (PLAN_CONSOLIDATION 1.2 (ee),
the root cause you found at gate 4 -- Fortran unit 6 does not reach fd 1 across the FFI on macOS).  Full 7-check
clean-clone gate.  Your gate-10 record (06bee43d) is merged beneath; the COAT_legacy pin is a bound now, with the
six other journals that pinned the same Cassegrain numerical zero.

**What changed (TO, macos 492998db / 8d397bf0 / 3c464098 / fc1c7789; resources 0b39b0a):**
- `emsg_mod` in the engine with `messages_count / messages_get / messages_clear` on the api; 228 print sites
  (the parser's warnings and refusals, CALIB's cause lines, the load notices) now write through the buffer as well
  as unit 6.  pymacos `messages()`, mmacos `macos.messages`; the bindings regenerated; cmdref regenerated.
- **The Darwin guards are GONE**: `tRxRefuse`, `tRxShortCoef`, `tRxLineEndings` (assumeTrue(~ismac)),
  `tSaveNewFileQuiet`, the `if ~ismac` on tCalibMatrix's trial-loses-chief cause line, and the pytest
  `skipif darwin` on test_rx_parse_warnings / test_save_new_file_quiet.  They read the buffer now (a nested
  MATLAB child echoes it tagged `MSG|`), so on the Mac they must PASS.
- **TO's rule for your gate, verbatim:** if one of those is red on the Mac, the buffer is missing a line -- a
  print site not converted (the census is in `macos_f90/notes/api_bindings.md`) -- NOT a test to re-guard.  Send
  the test name and the message text it expected; the fix is a site conversion in TO's lane.
- One darwin skip stays BY DESIGN: `test_messages.py`'s capfd parity leg (it compares the buffer against captured
  stdout, which macOS cannot supply).  Do not report it.
- Pins: the Cassegrain collimated RMS (1.29349D-12 on arm64) is bounded `< 1e-11` in all seven journals that carried
  it; the post-tilt 3.5664D-11 in ZABERR_ZRM_DOPD / DOPD_identical stays at five figures on purpose (a response,
  not a zero -- the header says why); if arm64 moves THAT one, it is a finding.

**Linux at 492998db (TO, both compilers):** load 160/160, cmd 118 + 2 must-fail, corpus 1468/1468, ab
623898d4 -> 492998db 0 movers; pymacos full 7014 / 1 skip; mmacos fast 695/0.

**Expected on arm64:** cmd 118/118 + 2 must-fail (no red: COAT_legacy bounded); load 165/165, 0 movers vs your
gate-10 base; mmacos fast 695/0 with NO filtered Darwin guards left except tCalibStop's OptTol leg (that one is a
trajectory, not a buffer; it stays assumeTrue); pymacos 7014 / 0 / 1 skip (the parity leg).  The number that
matters this time is the count of tests that were guarded and now pass: six MATLAB, two pytest.
