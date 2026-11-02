# BRIEF for TO (Opus): capture is a wrap problem -- unwrap the differential, then the start-rms ladder

From CCL for Dave, 2026-09-13.  Follow-up to `BRIEF_to_gauge_deck.md`
after your report (PF through the two decks; the shared knobs 14/14;
pdi_dm96; the descent finding).  Dave: push is OK -- push your five
commits now (MACOS_resources 5d99c20, e873703, 215a452, 2e0ab46; macos
e8c51d0); CCMac is blocked on the shared `dmg_loop` knobs and unblocks
the moment they are on origin.  Say when pushed.

## The finding, accepted

From a 100 nm rms surface the differential to the set point is ~2 rad
rms of phase, so the wrapped phase difference every reading returns is
wrapped whatever its absolute range: S and P diverge, PF barely moves.
Capture is a wrap problem, not a gain problem -- and it is the same
problem for every approach, the interferometer's four-step included
(its phase wraps at the same +-pi).  The re-measured matrix still holds
the gain at 120-160 nm (runs/cap385_b*); the wrap is what stops the
descent.

## Deliverable 8: unwrap the differential

1. A two-dimensional phase unwrapper in `dm_gauge_lib` (`dmg_unwrap.m`):
   least-squares (DCT-based, Ghiglia & Romero 1994) on the lit mask, with
   the off-mask pixels excluded, returning the unwrapped phase and a
   residue count (the number of 2 pi inconsistencies it could not
   resolve).  ~40 lines; no toolbox dependency (the release gate: no
   external paths).
2. Apply it to the wrapped differential of every phase reading (S, V,
   P, PF: `stepdiff`, `diffV`, the PDI `diff`) BEFORE the estimator,
   switched by one knob `battery.unwrap` (default off, so every record
   reproduces) and turned on by the loop when `loop.start_rms` is set;
   the linear reading L has no wrap and is untouched.  The DM's surface
   is smooth at the pixel scale (4 px per actuator at 385 rays), so the
   limit moves from the wrap (+-158 nm of surface) to the pixel gradient
   (adjacent pixels must differ by less than pi: several hundred nm rms).
3. Gates in tDmgLoop: a synthetic wrapped ramp + a random smooth surface
   of 1.5 waves unwrap to the truth to 1e-12 with zero residues; a
   gradient above pi per pixel is reported as residues, not silently
   wrong; `battery.unwrap` off reproduces the V3 record bit-for-bit
   (v3dev G4 = 0.296 pm at model 512 / NGRID 65 / grid 256 / 0.42,
   `dm_use 2`, `reg.mode record`, readings V).

## Deliverable 9: the start-rms ladder, both ways

Your gseq3 ladder, run with and without the unwrapper: start 30 / 60 /
100 / 150 / 200 / 300 nm rms, matrix measured at the start, gain 0.5,
`recal_every` 10 and never, photons per cycle 1e13 and 1e15, K 60,
readings L, S, V, P, PF (385 rays for the ladder; 193 is acceptable if
the box is the limit -- say which).  Report per reading: the largest
start that converges, cycles to 10 nm and to 3 pm, the final residual,
the residue count at the first cycle, and the photons per cycle it
took.  That table is the deck's capture slide.  The recalibrated
matrix's photon cost (5-40x at 120-160 nm) stays on top of it; state
both numbers.

Then the within-scan drift runs (`loop.intra`) as briefed, and the
reference-arm walk; CCMac mirrors the unwrapper in `tg96_run`'s
differential (`fsdiff_`) for the interferometer's descent -- keep
`dmg_unwrap` general (a masked phase map in, an unwrapped map out).

## Report

Append to `pdi_dm96/REPORT_gauge_pdi.md`: the unwrapper's gates, the
ladder table both ways, departures flagged, run tags.  Push when
tDmgLoop and the fast suite are green.
