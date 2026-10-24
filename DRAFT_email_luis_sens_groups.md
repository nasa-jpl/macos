Draft reply to Luis — for Dave to send (2026-09-14)
====================================================

Luis,

Thanks — all three are real, and two of them are on us to fix. Quick answers,
then what we're doing about each.

1) Grouped dw/dx running many hours vs ~54 min for the individual elements
--------------------------------------------------------------------------
You're right that it's abnormal, and it isn't the engine. The engine's group
perturb (GPERTURB) is just frame math over the group members — it's as cheap as
perturbing one element. The cost is in the MATLAB layer's grouped-channel object:
for every finite-difference poke it re-aims the chief ray through the object-space
stop (a full ChiefRayAiming + source-ray-grid rebuild), and if the group contains
the focal plane it also re-traces and re-finds the exit pupil. That runs three
times per Jacobian column (the +delta, the -delta, and the restore) — so on a
segmented 512 model a single group column costs roughly 45x what an element column
costs. That's the whole 8x: 3 groups can dominate 21 elements even though it's far
fewer columns.

The per-element path deliberately does none of that, which is why the 21-element
run is fast.

Fix (done, tested, on dev-candidate): the grouped channel now skips that
per-poke re-aim whenever it can't change the answer, and leaves the session at
the correct nominal aim afterward. The Jacobian is numerically unchanged —
verified by an A/B test (gate on vs off) — so it's purely a speed fix. Two
cases, both now covered:
 - Object-space stop (the default): the chief ray is aimed from the fixed
   source through a fixed point in space that no optic can move, so the re-aim
   is a genuine no-op — it's skipped for every group. This is your case, and it
   removes the whole per-poke cost, no knob needed.
 - Element-defined stop (ApStop= / macos.stop(elt)): the re-aim is skipped for
   any group sitting entirely downstream of the stop element (a rigid move
   there can't change which ray hits the stop), and kept for a group at/upstream
   of it.
Expect the group columns to drop back to roughly element-column speed — the 7 h
run should come back toward ~1 h.

If you're on the current binary and want the speedup right now without the
rebuild, run_sensitivities(..., 'group_stop_mode','none') also skips the re-aim;
just note that 'none' is only *exact* for groups downstream of (or with) an
object-space stop — with an element stop upstream it would change the numbers,
so prefer the rebuilt default, which decides correctly per group.

2) Link= iElt for grouped perturbation
--------------------------------------
Good news: it still works exactly as you remember, everywhere — the CLI, GMI, and
the MATLAB tools all share the same reader and perturb path, and Link=/DpElt= is
live (not deprecated or #ifdef'd out). Put "Link= 11" on element 12 and perturbing
11 rigidly drags 12, as it always did. In the MATLAB layer you don't need anything
special: if the loaded prescription carries Link=, perturb_elt on the master drags
the linked element automatically.

One gotcha worth knowing: SAVE does not currently write Link= back out, so if you
load a linked .in, perturb, SAVE, and reload, the linkage is gone in the saved
file. The linkage is intact for as long as the model stays loaded — it's only the
writer that drops it. We're adding Link= emission to SAVE so linked prescriptions
survive a round trip. Until that lands, keep your Link= in the source .in rather
than relying on a re-saved copy.

(Two structural notes, unchanged from the original design: only one linked element
per master, and Link= must reference a positive element id.)

3) Zernike basis types for non-segment surfaces
------------------------------------------------
Your observation was correct: the MATLAB *basis-map generator* for a plain
(non-segment) surface produced ANSI only. That's now fixed — the generator
(macos.zernike_grid_basis) takes a convention argument and produces engine-exact
Noll and Born & Wolf maps for any surface; dw_dgrid takes a 'zconv' option and
records the convention on the saved map so it's self-documenting.

Two things to separate:
 - Coefficient level (works today, all engine types 1-10): declare MonZernType=
   Noll (or BornWolf, the Norm* variants) in the Rx, or set it programmatically,
   and the finite-difference Zernike channel operates in that basis. The engine
   realizes those correctly — Noll included. (Only ExtFringe has no converter
   engine-side.)
 - Sampled grid/influence-map level (the dw_dgrid use case): now covered for
   ANSI, Noll and Born & Wolf. These are matched to the engine so a grid poke
   equals a MonZernCoef poke of the matching MonZernType (gated against the
   engine, per convention). Fringe / NormHex / annular-Noll are not offered yet
   (their orderings/radial polynomials need more work) — asking for one errors
   rather than returning a wrong basis.

One heads-up if you compare conventions: Noll and Born & Wolf coincide at
indices 1-3 AND at index 8 (both map to the same ANSI mode there), so a check
that only exercised those modes would see NO difference between them. They first
diverge at index 4 and clearly at index 7 (Noll 7 is coma, Born & Wolf 7 is
trefoil) — use one of those to tell them apart.

All three are in and tested on dev-candidate (the group speedup and the Zernike
generator on the MATLAB side; the Link=-on-SAVE fix is an engine change, so it
needs the rebuilt binary). Pull dev-candidate and rebuild when you get a chance.

Keep the emails coming — this batch was useful.

Dave
