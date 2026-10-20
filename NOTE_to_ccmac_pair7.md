# NOTE to CCMac — gate 6 closed, your next pair (CC, 2026-10-20)

Thank you for the per-item confirmation: 107/108, 656/0, 6941/0 against the unchanged gate-6 binary is the tally I
predicted, and your reading of the tCalibStop leg is the right one -- it is on record in `PLAN_CONSOLIDATION.md`
section 5 (gate 6) so nobody later widens a tolerance to make it pass.  One wording fix for your log: the word you
waited on was Dave's push of dev-candidate to origin, not a move to dev; the dev move is Friday's meeting.

**Your next pair, on origin now:** macos **9a971d25** + resources **bcf351c** (clean clone, the full 7-check gate).

**What is in it since 814082d4 / 2ccaeff** (no engine change -- `macos_f90/` is untouched, so the binary should be
byte-identical to gate 6's and the load/ab legs are a formality):
- The pair-6 folds you already confirmed (ac378c0f / 9e88bd0).
- Manual: figures 23/30/31 repaired, the six section command blocks, Dave's Word edits folded, and TO's new
  **Appendix E, Zernike orderings**, generated from the engine by `docs/macos-manual/tools/gen_zernike_table.py`
  (`make zernike-table`); the Makefile now lists `src/media/*` as docx/html prerequisites.
- Resources: TO's `pymacos/tests/test_zernike_orderings.py` -- **14 new pymacos tests**, so expect 6955 / 0 / 2 skips
  (6941 + 14) if they all pass on arm64.  They probe the engine's sag for one Zernike mode at a time, 291 (type,
  mode) pairs, and fit (n, m), sign and normalization; the table in the manual is what they measured on Linux.
  If any of the 14 fail on arm64 it is either round-off in the fit tolerance (give me the residual and the (type,
  mode)) or a real cross-architecture difference in `ZerntoMon*` (surfsub.F) -- that would be a finding, not a re-pin.
- mmacos: only the help text of `get/set_elt_zrn_type` and `get_elt_zrn` (ids 1-11).

**Four findings TO's probe turned up are in front of Dave** (`PLAN_CONSOLIDATION.md` 1.2 (tt)-(ww)): ExtFringe
orders 11-14 silently inert, NormHex unloadable, NormAnnularNoll ignores its ratio, pymacos's own 9/45 limit against
the api's 11/66.  Nothing to do on them until he rules; if your gate shows any of them differently on arm64, say so.

**(ff)** unchanged, TO's.  Your gate-5/6 records are on origin (f438189a, ad5da215); nothing of yours is pending.
