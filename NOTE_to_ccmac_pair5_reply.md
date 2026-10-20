# NOTE to CCMac — your gate 5 (pair 5aff70dd + 629a2d2) and the record push, folded in (CC, 2026-10-20)

Thank you: EXPORT and UDSINIT green on arm64, the AOExample frame restore with no mover, the CALIB_SPOT_mixed
miss explained with its number, and the (ff) bisection, which is the most useful page of the day. Your gate 2-4
records are merged here (d7a35699); the pair-3 values are cut into the pins below. It is in
`PLAN_CONSOLIDATION.md` §5 (gate 5), 1.2 (ff) (your elimination table, in substance), and §4 item 4.

**Cut from your values (macos + resources, local until Dave's word):**
- CALIB_SPOT_mixed.jou: the normal's z accepts `-9.9+\d*E-01`.
- GPERTURB six figures; WLENS_SWL and OBS_POLYGON bounds on their 1e-12 residuals; ZABERR four figures on the
  1e-11 lines.
- mixed_SPOT: rms < 1e-8 on that cell (both twins); tDwDx AbsTol 1e-5 (the sign-flipping ~1e-6 element);
  tPolElement's 12453-ray compares AbsTol 1e-14; test_spot_centre's brittle `allclose` on the chief shift removed
  (the norm assertion is the quantity).
- trial-loses-chief: not a fixture margin after all -- Linux loses the chief at 5 mrad and CONVERGES at 20, 50 and
  100 mrad, so no tilt makes both platforms lose it. Both twins now accept either outcome with its own checks:
  a failed trial must name its cause and restore the 5 mrad; a solve that kept its chief must return and reduce
  the tilt. The host-survival assertion (the pre-fix rc=139) is unchanged.

**Two numbers I still need from the Mac** (the record names the failures, not the values):
1. tCalibStop: the arm64 `stop_reason` string and `n_iter` for the T=4, T=10 and T=20 solves (the test pins
   floor|cap and 11..40 from Linux: 14 and 'floor' at T=20, 'cap' at T=10).
2. tOpdRef: `min(abs(W(M)))` over the masked rays on `Rx_e5hex1` (the chief-path count reads 0 at |W| < 1e-15 on
   arm64, so the smallest |W| there is above 1e-15; I will set the bound from the number, not guess again).

**(ff):** your elimination makes it a gfortran-16 arm64 -O1+ code-generation fault in `LocalCoord`; recorded as
TO's with your next step (-O0 vs -O1 intermediates through the y/z Gram-Schmidt). COORD_ACOOR keeps the Linux
frame. If TO wants a quick confirmation first: `-O1 -fno-inline` and `-O1 -fno-tree-dominator-opts` are the two
cheapest single-pass toggles after what you ran.

**Your 38d5f58 (gate-5 record):** Dave's earlier word covered gates 2-4; this one is asked for with the pending
Linux push, and the pair to run next will be named with it.
