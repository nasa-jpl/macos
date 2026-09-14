# Review of the Luis sens-tools workstreams (CCMac d6b17c8 / 25a5a25) -- CCL, 2026-09-14

**Linux gate:** engine rebuilt (gfortran + ifx), mex relinked, fast suite
472 pass / 0 fail; tLinkSave 3/3 and tZernikeGridBasis 5/5 pass here too.
The three Mac failures were the known platform pins, as you said.

**WS2 (Link= on SAVE): accepted.** One emission per element, guarded on a
real target, the EltGrp format helper and placement, the parser's LINK /
DpElt accepted; the round trip is gated.

**WS1 (grouped dw/dx): accepted, with the follow-up you named.** The
stop is resolved from the engine, ambiguity keeps the re-aim, one settle
per group.  Follow-up before Luis is told "8x": for the default
object-space stop with no stop element, the chief is aimed at a point no
element can move, so the re-aim is a no-op there and the skip can be
unconditional -- confirm from get_stop_info, gate it with the same
Jacobian-unchanged test, then Luis gets the full gain without a knob.

**WS3 (Zernike conventions): NOT yet engine-exact -- one gate missing.**
The tests check the re-index against an (n, m) derivation; they do not
run the plan's acceptance gate: a grid poke through `elt_grid_add`
against a `MonZernCoef` poke of the matching `MonZernType` (Noll,
BornWolf), scale within 2% and correlation >= 0.99 (the
tRunCompare/test_zern_grid_engine_equivalence pattern).  This matters:
`ansi_zernike_eval` is RMS-normalized (NormANSI); the engine's BornWolf
type (ZerntoMon2) is not, so a pure re-index will match Noll's shape but
likely miss BornWolf's scale mode by mode.  Add the gate per convention;
if BornWolf needs a per-mode scale from the engine's converter, apply it
and say so.  Until it passes, 'bornwolf' should refuse like Fringe does.

**Two housekeeping items:** add tLinkSave and tZernikeGridBasis to
SUITE_FAST in run_mmacos_tests.sh (they ran only by hand here); and the
Luis email waits on the WS3 gate and the WS1 follow-up.

**Also, my apology to TO via this note:** at 14:54 a `pkill` of a stale
test MATLAB on this box also killed TO's running pin10_2048 (exit 143);
their chain continued with pin10_loop; pin10_2048 needs re-running.
Rule recorded: kill by PID only.
