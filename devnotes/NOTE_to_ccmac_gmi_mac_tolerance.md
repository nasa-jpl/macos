# NOTE to CCMac — GMI regression tolerance, branch `mr_merge_gmi_mac` (CC, 2026-11-18)

Thank you for the first Mac build of GMI. Dave has decided the tolerance question, and the change
is a second commit on the same branch. Please pull, rerun the regression, and open the PR.

## What to do

1. `git fetch origin && git switch mr_merge_gmi_mac && git pull`. No rebuild is needed: only the
   regression's MATLAB files and two documents changed.
2. `app/GMI/regression/run_regression.sh`.
3. Expected on the Mac: **7 passed, 0 failed (of 7)**, with these lines:

       tolerance vs reference: 32 x eps x (chief-ray path length)
         e2e_pie 2.271e-13 m, e5hex1 4.911e-10 mm (absolute |a-b|)

   and your earlier differences inside them: about 5.7e-14 and 8.5e-14 for tests 02 and 03, about
   1.2e-10 and 8.7e-11 for tests 05 and 06.
4. If it is green, open the PR from `mr_merge_gmi_mac` with base `mr_merge_remove_scratch`, using
   the PR text at the end of `devnotes/NOTE_to_ccmac_gmi_mac.md` (updated for this change), with
   your Mac line added. No issue number is needed for it; Dave has said so. The same holds for
   the manuals branch.
5. Record the Mac result in `MAC_PORT.md`, Part 2. Part 1 already has a short account.

## What changed

- **The rule.** The regression compared against its references at a fixed 1e-12. It now allows
  32 x eps x L, where L is the prescription's chief-ray path length and eps the machine epsilon
  (`app/GMI/regression/lib/ref_tolerance.m`). An OPD is the difference of two paths, so it is
  good to the round-off of the path.
- **Why your numbers settled it.** In units of eps x L your four differences are 8.0 and 12.0
  (e2e_pie, L = 31.96 m) and 7.6 and 5.7 (e5hex1, L = 69,121 mm). The two prescriptions disagree
  with Linux by the same amount; the fixed 1e-12 passed one and failed the other only because one
  deck is in metres and the other in millimetres.
- **What stays exact.** The two calls of one build, in tests 02 and 05, are still compared at
  tolerance 0.
- **A seventh test**, `test07_tolerance_refuses`, needs no engine: half the tolerance must pass
  and twice the tolerance must fail, for both prescriptions.
- **The path lengths** are the engine's own "Total Length" for the chief ray traced to the last
  element of each deck.

## If it is not green

Send the summary table again with the largest differences. Do not re-make the references and do
not change the 32. A difference above 32 units on the Mac would mean the margin (2.7 times on
e2e_pie, 4.2 times on e5hex1 against your first run) is too thin, and that is Dave's call again.

## Still to watch

Issue #82 (`mr_merge_remove_scratch`, note `devnotes/NOTE_to_ccmac_issue82.md`) and the manuals
branch (`mr_merge_manuals`, note `devnotes/NOTE_to_ccmac_manuals.md`) have no PR yet.
