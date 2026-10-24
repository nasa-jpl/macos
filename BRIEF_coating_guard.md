# BRIEF: `Coating=` parser unguarded against mCoat — bound check + gate

_Small ENGINE fix, post-demo queue (Dave sequences; any agent).
Origin: TO's tg_psi_dm_v2 build (delivery log at the foot of
`BRIEF_tg_ifo_v2.md`, §ENGINE FINDING), 2026-08-29.  Found the hard
way: an 11-layer MacNeille stack._

## The defect

`msmacosio.inc:2753` reads `EltCoat(iElt)` and loops
`Do k=1,EltCoat(iElt)` writing `IndRefArr(k,iElt)` /
`ExtincArr(k,iElt)` / `EltCoatThk(k,iElt)`, then `IndRefArr(i+1,iElt)`
— with NO bound check — while `IndRefArr` is `(0:mCoat,mElt)` and
`EltCoatThk` is `(mCoat,mElt)` with `mCoat = 10` (`elt_mod.F:330`).
An 11-layer `Coating=` block LOADS WITHOUT COMPLAINT and writes past
the array ends; the only visible symptom is `coat_get` failing
afterwards (it does guard `nLayer > mCoat`).  Memory-corruption
class: adjacent elt_mod arrays are silently clobbered.

## The fix

Reject `EltCoat > mCoat` at parse time with a named message
(element index, requested layers, mCoat ceiling), the way the
validator handles other bad values — check BOTH the interactive and
SMACOS paths (msmacosio.inc is shared, but confirm `coat_set` too:
it may carry its own guard; verify rather than assume).  Consider the
same audit for the Model-B arrays (`nCoatElt`/`CoatThkElt`) while in
there — same pattern risk.

## Gates

- A 10-layer stack loads and round-trips (boundary, not off-by-one).
- An 11-layer stack is REFUSED with the named message on load — CLI
  and binding paths both.
- Non-vacuity: the refusing gate shown to pass (i.e. load silently)
  against the pre-fix engine.
- `pbs_macneille.m` carries a MATLAB-side assert meanwhile — leave it
  in place even after the fix (defense in depth), but its comment
  gains the brief reference.

## Notes

- ChkDf2/validator conventions: see `validate_prescription` and the
  §0 hygiene cluster in `macos_f90/CLAUDE.md`.
- Engine repo work → dev branch; cherry-pick to opt-dev per the
  branch model.  Gates in the mmacos suite need a fixture .in with a
  deliberately over-long coating.
