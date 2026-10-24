# BRIEF for Terminal Opus: the "center-channel speckle" in dwdsurf_5zoom_5fov (quick exploration)

From CCL for Dave, 2026-09-08.  Scope is deliberately small: dw/dsurf,
centre zoom, centre field, one deck.  Budget half a day; stop at the
first measured attribution and write it up.  CCL stays clear for the
Fang work; report to Dave.

## 0. What Luis saw (as relayed by Dave)

Running the zoom x field sensitivity harvest on
`~/dev/MACOS_sandbox/sens_noise/jwst_ote_designc.in` (his 9-zoom x 5-fov
`dwdsurf_9zoom_5fov_rx_sens.mat`; viewer scripts in the same directory),
a low-level SPECKLE pattern appears in the dw/dsurf columns, FIXED for a
given (field, zoom), in the "center channel" only.  Deck facts: it is
the committed zoom fixture `templates/50_sensitivities/zoom_5x5/
jwst_ote_designc.in` plus three header lines -- `UseChfRay4OPD= Y`
(chief-ray OPD reference instead of the aperture mean), `ApStop=`
(object-space stop; the committed runner uses `stop_elt 25`, the FSM),
`PgplotImage= Color`.  Runner settings: model 512, NGRIDPTS 63, FOV
2.90888e-4 rad, stop elt 25, dw_dsurf default delta 1e-6 (Kr in
BaseUnits = mm), central differences.

## 1. Read first (cold start)

`MACOS_resources/mmacos/CLAUDE.md`; `mmacos/doc/SENSITIVITY_TOOLS.md`;
`templates/50_sensitivities/zoom_5x5/README.md` -- in particular the
paragraph on **element 4 (CenterSegment)**: a VIRTUAL element that passes
only the chief-ray sliver (~2.5 % of the beam) and exists to carry the
chief ray and reference the PM.  Its dw/dx columns are a documented NULL
FLOOR (norms ~1e-4 vs ~200 for a real segment; "the dead elt-4 null
floor 5.30e-05") and `flag_zero_norm_channels` + `drop_channels` remove
it from the dw/dx rung number-free.  dw/dsurf is a Kr/Kc-per-element
rung and elt 4 is a Conic with finite Kr, so **dw/dsurf still carries an
'Elt 4 Kr' / 'Elt 4 Kc' channel** whose true response is ~zero.  Also
read the FEX and pupil sections of `macos/macos_f90/CLAUDE.md`
("FEX probe is FRAME-INDEPENDENT", "Element STOP preserves the source
frame's HANDEDNESS") -- the engine changed on 2026-09-08; build the
current `dev-candidate` engine (`44fc362`+) and relink the mex before
measuring anything, and say which engine you used.

## 2. Step 0 -- which "center"?  (do this before any hypothesis)

"Center channel" has two readings and they lead to different work:
(a) the CenterSegment ELEMENT's channel (elt 4 Kr / Kc) in every
(zoom, field) block, or (b) the centre (zoom, field) BLOCK for every
element.  Settle it with Luis's own viewer on a regenerated block --
do NOT rerun the 9x5 harvest; the single-DOF driver is enough:

```matlab
run('<mmacos>/mmacos_setup.m'); macos.init(512);
rx = '~/dev/MACOS_sandbox/sens_noise/jwst_ote_designc.in';
m = macos.Session(512);
o = macos.dw_dsurf(m, rx, 'elts', [4; 5; 23; 24], 'verbose', true);  % elt 4 = CenterSegment, 5 = a real segment, 23/24 = SM/TM
% o.dwds is Nw x Ns, o.channel_names names the columns, o.w_nom_2d + mask give the canvas
```
(centre zoom = no FSM rotation = the deck as loaded; centre field = the
deck's nominal chief.)  Reconstruct each column onto the canvas
(`macos.v2m` with the mask from `w_nom_2d ~= 0`) and look: if the
speckle is in the elt-4 columns and not in elt 5 / SM / TM, it is reading
(a); if it is in every column of this block, reading (b).  Record the
column norms: elt 4 vs elt 5 (expect ~1e-4 vs O(1e2)-class per mm of Kr).

## 3. If reading (a): a null-floor channel (expected outcome)

The elt-4 response is ~zero, so its column IS the finite-difference
noise floor: the OPD is a difference of ~1e4-mm path lengths, whose
round-off (~1e-11 mm) divided by 2*delta = 2e-6 mm gives ~5e-6 mm per
mm-of-Kr of speckle -- fixed per (field, zoom) because the ray geometry,
hence the round-off pattern, is fixed.  Two measurements close it:
- delta scan: `'delta'` in {1e-7, 1e-6, 1e-5, 1e-4}, `'method'`
  'central' and 'forward'.  Round-off noise scales as 1/delta and is
  independent of method; a geometric artefact does not scale.  Report
  the rms of the elt-4 column vs delta, and of the elt-5 column (its
  SIGNAL must be delta-independent -- that is the check the delta is
  in the linear regime).
- the zero-norm test: does the dw/dsurf rung apply the same
  `flag_zero_norm_channels` / `drop_channels` treatment as dw/dx?  If
  not, that is the fix to recommend (number-free, by response norm),
  mirroring the README's dw/dx rule; do not hard-code elt 4.
Write it up with the numbers and stop.

## 4. If reading (b): a block-wide artefact (the interesting case)

Then something in the centre block differs from the others.  Measure in
this order, one hypothesis per measurement, and stop at the first hit:
1. Repeatability.  Run the identical `dw_dsurf` twice in ONE process and
   once more in a NEW process.  Identical columns -> deterministic
   (round-off or geometry); different -> jitter between traces.  If
   different, the difference pattern is the clue: compare it to the
   nominal OPD gradient (a re-aim / regrid jitter shows as grad(W)
   times a per-call shift -- one fixed spatial pattern, varying only in
   scale across columns, which is exactly "fixed to field and zoom").
   Check whether the stop re-aim (`Computed StopPos` prints,
   ChiefRayAiming's convergence residual) runs between the +/-delta
   traces; if it does, its tolerance is the suspect.
2. delta scan as in section 3.  1/delta -> round-off; flat -> geometry.
3. Mask flicker.  After each +delta / -delta trace, `macos.get_ray_status`
   and compare the pass masks: rays whose status flips between the two
   traces produce speckle at those pixels.  The centre field puts the
   hex ray grid symmetric on the segment gaps, so edge rays are the
   natural suspects.  Overlay the speckle pixels on the segment-gap map
   (a piston poke of each segment gives its footprint).
4. OPD reference.  Rerun with `UseChfRay4OPD= Y` removed from a copy of
   the deck (aperture-mean reference) and with it kept.  A chief-ray
   reference turns any jitter of ray 1's own path into a per-trace
   PISTON (uniform, not speckle); if the speckle changes with the
   reference, say how.
5. Pupil placement.  The columns are differences against a reference
   sphere held fixed through the pokes; confirm FEX is not re-run
   between the +/-delta traces (the FEX prints `EP crossing from 4
   probes` would appear twice per channel).

## 5. What NOT to do

No engine changes, no new guards, no tolerance changes in tests.  If
the attribution points at engine behaviour (re-aim tolerance, grid
regeneration, status flicker), write the measurement and the proposed
fix; Dave and CCL decide.  Do not regenerate the committed zoom_5x5
artefacts.  Deck names and numbers only in the write-up (the deck is
Luis's copy; the committed fixture is public).

## 6. Deliverable

`macos/REPORT_sens_noise_center.md`: which "center" it is (section 2
evidence), the attribution with its measurement (rms vs delta table, or
the repeatability / mask-flicker / reference result), and one
recommendation.  Keep the MATLAB you used under
`MACOS_sandbox/sens_noise/` (not committed) and cite it from the report
with the deck's provenance line.  Record the resolved oddity in the
zoom_5x5 README's oddities section at resolution time (the standing
rule), and tell Luis in one paragraph.
