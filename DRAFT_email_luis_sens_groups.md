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

Fix (done, in review): the grouped channel now (a) stops doing the re-aim on the
restore step, where it's pure waste (a universal ~1/3 saving), and (b) skips the
re-aim entirely when the whole group sits downstream of the aperture stop (moving
it rigidly can't change the chief-ray aim). The session is left at the correct
nominal aim afterward, and the Jacobian is numerically unchanged — verified by an
A/B test (gate on vs off) — so it's purely a speed fix. Interior/downstream group
columns drop back to roughly element-column speed.

One important detail on getting the FULL (b) speedup: the code resolves the stop
element from the engine, and that only works when the deck defines an ELEMENT
stop (an "ApStop=" in the header, or macos.stop(elt)). Under a pure object-space
stop (the default group_stop_mode='obj'), the engine can't report a stop element,
so the gate stays safe-but-conservative and you get only the ~1/3 (part a). If
your FSM decks have a defined stop element, run the groups with
group_stop_mode='elt' and that stop element to get the full skip; if you know a
group is entirely downstream of the stop, group_stop_mode='none' is exact and
fastest.

Interim workaround you can use right now, no rebuild:
    run_sensitivities(..., 'group_stop_mode','none', 'group_fp_mode','none')
Important caveat: 'group_stop_mode','none' is only *correct* for groups that are
entirely downstream of the aperture stop. If a group contains or sits upstream of
the stop, turning it off will change those numbers — so use it only for the
downstream groups, and keep the default for any group at/before the stop.

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
Your observation is correct: the MATLAB *basis-map generator* for a plain
(non-segment) surface currently produces ANSI only. The recent self-contained Noll
work went into the per-segment basis builders, not the general non-segment
generator, so Noll/Fringe/Born-Wolf sampled maps aren't available there yet.

Two things to separate:
 - Coefficient level (works today, all engine types 1-10): declare MonZernType=
   Noll (or Fringe, BornWolf, the Norm* variants) in the Rx, or set it
   programmatically, and the finite-difference Zernike channel operates in that
   basis. The engine realizes those correctly — Noll included. (Only ExtFringe has
   no converter engine-side.)
 - Sampled grid/influence-map level (the dw_dgrid use case): this is the gap. We're
   extending the non-segment generator to emit engine-exact maps in the other
   conventions — matched to the engine so a grid poke equals a MonZernType
   coefficient poke — and we'll label each saved map with its convention so it's
   self-documenting. We'll validate it in both grid orientations (the coma/trefoil
   orientation you've run into before).

Priority order on our side: the group speedup first (that's the one biting you
now), then the non-segment Zernike bases, then the SAVE/Link= fix (engine change,
needs a rebuild).

Keep the emails coming — this batch was useful.

Dave
