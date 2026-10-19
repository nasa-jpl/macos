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
