# NOTE to CCMac — your next pair: macos ad4161fa + resources 00697be (CC, 2026-10-21)

Both on origin dev-candidate once Dave's word lands this stack.  Engine CHANGED since 3ead7473 by ONE line -- yours:
the full 7-check clean-clone gate, and this time the result you care about most is COORD_ACOOR.

**What changed:** PLAN_CONSOLIDATION 1.2 (ff) FIXED (macos 0efc2dec): `LocalCoord`'s y rotate branch now recomputes
`yLocal` from the rotated ray's `RayPos` before `DUNITVEC`, exactly your candidate line from
`cli_tests/records/mac_ff_localcoord.md`, with a comment crediting the bisection.  Linux was wrong too, which no
gate had seen: on NON-coaxial decks the stale residual is real on x86_64 (opt_example.in's focal-plane frame 1.4
degrees off its normal, e5hex1.in's 12 degrees, on every build) -- the coaxial Cassegrain of COORD_ACOOR hid it
with an exactly-zero residual.  New journal `cli_tests/cmd/COORD_offaxis_frame.jou` pins both decks' zLocal at five
figures (= the element normal; e5hex1 with the beam-side sign).  Linux gate at 0efc2dec: load 160/160, ab 0 movers,
cmd 113/113 + 2 must-fail; pymacos side tree 6963/0.  Your 1.1b closure (82c80e1b) is merged beneath.

**Expected on arm64:** cmd 115/115 + 2 must-fail -- COORD_ACOOR passes for the first time on the Mac (your -O2
verification already showed it) and COORD_offaxis_frame should read the same five figures as Linux; if its last
digits differ there, send them, but the pins are already 5-figure.  load 164/164, 0 movers (none on Linux); mmacos
fast 656/0 + the filtered leg; pymacos 6961 + 0 new = 6961 / 0 / 2 skips.  With that, the Mac line has NO engine
miss left: (ff) was the last.  Update MAC_PORT.md's Mac-only deltas accordingly in your record.

Thank you for the bisection -- the method (instrument, diff -O0 vs -O1, first divergent intermediate) is now in
`notes/pupil_fex_stop.md` as the way to read a "compiler-only" fault.
