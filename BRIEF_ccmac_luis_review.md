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

---

## Round 2 (CCL, 2026-09-14 evening): on 4772c83 / c44d179

**Linux gate on your push:** tDwDxGroups 15/15, tZernikeGridBasis 5/5,
tLinkSave 3/3, both tRunCompare zern_grid gates pass; fast suite result
in the slice.  WS1 object-space skip accepted (the catch branch now means
"no element stop", which is the only way get_stop_info fails on a live
session).  'bornwolf' documented as NormBornWolf in the veneer header, so
my BornWolf-scale worry does not apply -- agreed.  Housekeeping accepted.

**WS3 gate: passes, but it cannot see a Noll / Born & Wolf swap.**  Noll
index 8 and Born & Wolf index 8 both land on ANSI slot 9 (your two perm
tables: `perm(9) = 8` in each), so at mode 8 the two conventions are the
same map and the same poke -- which is why your gate reports identical
numbers for both (1.0101 / 0.9938 here).  Replayed with the pairings
crossed (`scratchpad/ws3_teeth*.m`, model 512, ng 128, sampling 31):

| map \ poke      | NormANSI | NormNoll | NormBornWolf |
|-----------------|----------|----------|--------------|
| ansi, mode 8    | 1.039 / 0.969 | 0.049 / 0.053 | 0.049 / 0.053 |
| noll, mode 8    | 0.058 / 0.053 | 1.010 / 0.994 | 1.010 / 0.994 |
| bornwolf, mode 8| 0.058 / 0.053 | 1.010 / 0.994 | 1.010 / 0.994 |
| noll, mode 7    |  --      | 1.039 / 0.969 | -0.089 / -0.109 |
| bornwolf, mode 7|  --      | -0.135 / -0.128 | 1.035 / 0.964 |

(scale / correlation).  Mode 7 separates them (Noll 7 = ANSI 8 coma,
Born & Wolf 7 = ANSI 10 trefoil) and the code is RIGHT there too -- but
the correct pairing reads 1.039 / 0.969, outside your 2% / 0.99, so
switching the gate to mode 7 as written would fail the right answer.

**The threshold is grid discretization, not a convention error.**  On a
NormANSI deck, every mode 4..10 against its own poke:

| ng / sampling | scale range   | corr range     |
|---------------|---------------|----------------|
| 128 / 31 (the gate) | 1.010 .. 1.039 | 0.964 .. 0.995 |
| 128 / 63      | 1.013 .. 1.031 | 0.981 .. 0.993 |
| 256 / 31      | 0.993 .. 1.008 | 1.0000 (all)   |
| 256 / 63      | 1.003 .. 1.018 | 0.984 .. 0.997 |

At 256 points with the rays coarser than the grid pixel the bilinear
facets vanish and every mode agrees to 1% with correlation 1.0000.

**Fix (gate only, ~20 lines):** `ng 256`, sampling 31; modes {4, 7, 8}
per convention (4 and 7 are where Noll and Born & Wolf differ; 8 is
where they coincide -- keep it as the "same slot" check); thresholds 1%
/ 0.999; and the NEGATIVE control in the same test: the noll map against
the NormBornWolf poke (and the reverse) at mode 7 must give |corr| < 0.5,
so the gate proves it can see a swap.  Then the Luis email is unblocked.
Note for the email: Noll and Born & Wolf coincide at index 8 (and 1..3),
so a user who checked only those would see no difference between them.
