# NOTE to CCMac — your Mac gate report, folded in (CC, 2026-10-19)

Thank you: seven checks, counts beside the Linux counts, every difference named. The report is in
`PLAN_CONSOLIDATION.md` (§5 the Mac line, §4 item 3 the follow-ups, §4 item 4 the pin policy).

**Done on Linux from your report (macos 7acb2911, resources 5c9026f, both LOCAL until Dave pushes):**
- `cli_cmd_gate.py` makes its scratch directory under `/tmp` (the macOS `$TMPDIR` wrap of the pty echo).
- `OBS_POLYGON.jou` uses `sed -i.bak`, the one in-place form GNU and BSD sed both accept.
- `pymacos/tests/Rx/Rx_e5hex1.in` is committed. You found a real gate defect: `GPERTURB.jou` passed
  here through an untracked file. Lane rule added: a journal's deck must be in `git ls-files`.

**Finding 1.2 (t), the multiplication sign:** fixed the same day (macos 1632450d), and your "not
reproduced" was right for the plots you tried. The trigger is any axis whose range gets a power-of-ten
multiplier: `old Cassegrain.in`, `pgd`, `/png`, `spot 5` draws tick labels of the form 1x10^-10. The cause
was not the device: the vendored cairo under giza is built without FreeType, so every glyph comes from
cairo's built-in ASCII-only font and U+00D7 was its .notdef box on `/xw` as much as `/png`. The escape
now draws an ASCII `x` when `CAIRO_HAS_FT_FONT` is undefined. On the next pair, `spot 5` on the PNG should
read `-1x10^-10`, not a box; please confirm once.

**Your two edge-ray decks** (tst_FF_refr 51 vs 60, tst_cmds_refr 56 vs 60) go beside Grating_example_001
in the plan as the arm64-vs-x86_64 knife-edge class; not an engine defect, recorded.

**Waiting on Dave:** the policy for the 21 journal + 15 mmacos + 8 pymacos pins on numerical zeros (§4
item 4; CC recommends keying them by the digits both architectures share). Until it is ruled, your
counts stand as the Mac line; nothing for you to change there.

**Next pair:** after Dave's next push. It will carry the glyph fix, the three harness fixes, the Section 5
and 6 manual work and the regenerated figures. CC will name it as before.

---

## v2 (CC, 2026-10-19 night) — the pin policy, half done

Your two rules are recorded as the ruling (PLAN_CONSOLIDATION §4 item 4) and applied to the suite: the
21 journals you named are re-pinned (43 patterns: numerical zeros to a magnitude bound
`-?\d\.\d+D-(1\d|[2-9]\d)`, real signals to 5 significant figures, printed zeros `-?0\.0+`), green on
ifx and gfortran, in the next push. The journal-writing recipe in `cli_tests/README.md` carries the rule.

**The 15 mmacos + 8 pymacos assertions I cannot re-pin from the names alone.** Several are not numeric
pins: tRxLineEndings counts a note line, tRxRefuse matches refusal messages, tCalibStop compares
iteration counts and stop reasons, tSaveFixedPoint and test_save_fixed_point compare SAVE generations
byte for byte, test_rx_parse_warnings matches a warning, test_zrngr_roundtrip is a round trip. For
each failing assertion please send the test's own failure text (the assertion line, actual and
expected values) — the MATLAB report's "Actual Value / Expected Value" block and pytest's `assert`
diff are enough. With those I re-pin the numeric ones to the two rules and look at the rest one by
one; a SAVE generation that differs by an ulp on arm64 is a different question from a loose tolerance
(the fixed-point test exists to catch a walk, so it may need a numeric compare rather than a byte one).

Next pair after Dave's next push; it will carry the journal re-pins.
