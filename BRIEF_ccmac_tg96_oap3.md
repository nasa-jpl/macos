# BRIEF for CCMac: tg96-oap review of f46585d -- accepted in substance; six items before the merge

From CCL for Dave, 2026-09-11.  Reviewed `origin/tg96-oap` at `f46585d`
(`REPORT_oap.md`, `tg96_run.m` stage MATRIX / PLACE / D4, `tg96_place.m`,
the `dmg_frame` extension).  Written for a cold start; pull both repos'
`origin/dev-candidate` first (Dave pushes macos dev-candidate with the
loop addendum; MACOS_resources dev-candidate carries `dmg_loop`).

## Verdict

The last mile is real.  The affine placement is the right route and the
matrix calibration is a faithful lift (mean-reference, the rank-one
piston term, Cholesky; I read `stage_matrix_`, `build_J_`,
`est_matrix_tg` line by line against the ZWFS original).  The shared-lib
change (`dmg_frame` gains an optional third output) is backward
compatible; `zwfs_run` still calls it with two.  f46585d touches neither
the bench nor `tBench`, so the 'lens' byte-identity verified on
2026-09-11 stands.

**Reproduced here (Linux, mexa64, model 1024), lens rig, matrix mode:**

| row | CCMac (Mac) | here, tuned tail (tag lens) | here, UNTUNED tail (tag vlens) |
|---|---|---|---|
| flat-DM null | 0.1345 nm | 0.1345 nm | 9.111 nm |
| Stage C single 150 nm | 0.9968 / 98.8 pm | 0.9968 / 98.8 pm | 0.9928 / 21.7 pm |
| flat / single 10 nm | 0.9916 / 2.2 pm | 0.9916 / 2.2 pm | 0.9939 / 2.5 pm |
| flat / random 10 nm | 0.9894 | 0.9893 / 169 pm | 0.9948 / 135 pm |
| modal transfer 0.7-68 cyc/pupil | 0.96-0.99, x-talk < 6 % | 0.96-0.99, x-talk < 6.3 % | 1.00-0.93, x-talk < 2 % |
| break ladder | holds to 120, breaks 240 | holds to 120 (1.03 / 19 pm), breaks 240 (1.93 / 394 pm) | holds to 60, breaks 120 (gain -9.0: a phase wrap) |

Every D2 number reproduces to the printed digit on the tuned bench.  The
matrix beats the kernel record on the lens (0.9654 / 21 pm) on both
machines and even on the untuned bench -- which is itself evidence for
the D4 claim that the null is common-mode in the differential.

## Items (numbers first, then what to do)

1. **The D2-D4 evidence is not in the commit.**  `REPORT_oap.md` cites
   `runs/dev{lens,oap}{,mx}/`; f46585d carries none of them, and the
   committed `runs/lens/lens_report.txt` is the KERNEL record run
   (0.9654 / 21 pm).  Dave's standing rule: an uncommitted artifact is an
   unverifiable claim.  Commit the run reports + `.mat` for the lens
   matrix run, the OAP matrix run, the D1 gate runs and D4 (the
   `zwfs_dm96/runs/<tag>/` pattern; `zwfs_flat.txt`-class bulk is
   gitignored there, do the same for `tail_flat.txt` copies).

2. **The tuned tail is keyed by the run TAG** (`stage_B_`:
   `<tag>_tail.mat`).  Any tag other than `lens` / `oap` silently runs
   the geometrically-scaled seed: my `vlens` run got the 9.111 nm
   pre-retune null instead of 0.134, and everything downstream on a
   different bench.  Key the tail by `P.bench.optics` (`lens_tail.mat`,
   `oap_tail.mat`), keep `<tag>_tail.mat` as an override, and print the
   null beside its expected value with a WARN when it is 10x off.

3. **The OAP dense-pattern loss (random 10 nm 0.7486 / 4848 pm; modal
   cross-talk ~0.42) is attributed to physical astigmatism, but two
   software mechanisms are not excluded, and each is one run:**
   (a) *Window truncation.*  An elongated response's tail outside its
   own +/-half-step window lands in its neighbors' windows in a dense
   pattern and is assigned to them; a single poke has no neighbors.
   That is exactly the single 0.99 / dense 0.75 signature.  Test: run
   the OAP matrix with `matrix_step` 16 (windows twice as wide, 256
   states) -- or, better and cheaper, assign every mask pixel to its
   NEAREST poked actuator of that multiplexed frame (Voronoi cells:
   nothing truncated, nothing double-counted; a 10-line change in
   `build_J_`).  If the dense gain moves toward 0.99, the fold cost was
   the window, not the optics.
   (b) *Regularization shrinking the dim columns.*  ~25 % of actuators
   image dark on the OAP; their column energies are small and
   `matrix_lam` = 1e-3 of the MEDIAN energy shrinks exactly those.  Test:
   `matrix_lam` 1e-4 and 1e-5 on the OAP; report the dense-random gain
   over the bright actuators and over the dark 25 % separately.
   Only what survives (a) and (b) is the same-plane fold's physical cost.
   The 98.9 % "containment" number does not settle (a): say whether it
   is a fraction of actuators or a fraction of response ENERGY -- only
   the latter bounds the truncation.

4. **"Worst at centre" for D1 wants a picture.**  A pupil imaged by an
   off-axis parabola should not be darkest at its centre.  Plot the D1
   error and the column norm over the actuator grid for both rigs; an
   optical cause is smooth and symmetric about the fold plane, a
   chief-pixel-reference leak is a dip at the centre.  (Mean-referencing
   after a chief-pixel subtraction removes only the constant; check the
   four-step phase is referenced BEFORE any wrap.)

5. **The non-vacuity line is stated for the wrong rig.**  Here the best
   axis-aligned parity+scale map ALSO passes the lens D1 gate (100 %
   within 2 px, median 0.14 px) -- as it must: the lens mapping IS
   axis-aligned.  The check is meaningful on the OAP only; report it
   there with its number (what fraction the parity map reaches vs the
   affine's 72.8 %).

6. **`calib_surface 'base'` on the ladder** (S10 doctrine: the matrix
   measured on the working surface) -- one lens run and one OAP run;
   the 240 nm break is the flat-calibrated matrix.  Note the break here
   at 120 nm on the untuned bench reads gain -9.0 with corr -0.68: a
   four-step phase WRAP, not a gradual loss -- worth a wrap guard (the
   differential map's range vs lambda/2) so the runner says "wrapped"
   instead of printing a negative gain.

## Then: deliverable 7, the closed-loop hold metric

`BRIEF_ccmac_tg96_oap2.md`, addendum (written 2026-09-11 after your
merge; you have not seen it until Dave pushes macos dev-candidate).
Merge `origin/dev-candidate` (MACOS_resources) into `tg96-oap` again to
get `dm_gauge_lib/dmg_loop.m` + `tests/tDmgLoop.m`, then the `tg96_run`
stage 'loop' exactly as `zwfs_run`'s `stage_loop_`: four handles + lit,
seed 77, walk 2 pm / thermal 5 pm, steps 1 and 10 nm, photons
{1e12..1e15} PER MEASUREMENT (Dave: say "measurement", never "state"),
g 0.5, K 60, `hold_spec` 3 pm.  The ZWFS record to compare against is
`zwfs_dm96/README.md` S11 (L and S: 3 pm at ~2e12 noise-only, ~7.4e12
under the walk; S no fixed error; I+ diverges).  Lens rig first.

## Constraints (unchanged)

No engine work; 'lens' byte-identical (tBench); commits on `tg96-oap`;
push when items 1-2 are in (they are hygiene); Dave orders the merge.
One MODEL-1024 MATLAB at a time.

## Appendix: the two reproduction reports (this box, 2026-09-11)

### tag lens (tuned tail)

```
Stage MATRIX -- measured response matrix dw/da (lens rig, step 8, sign same):
  flat-DM null: 0.1345 nm rms surface (134.5 pm) -- the same-plane-fold arm difference
  J: 64 states, 3228 columns (lit), window 41 px, reg lambda 1.00e-03 (of median col energy)
  Stage C single actuator @150nm: gain 0.9968, off-target floor 98.8 pm
Stage D -- modal transfer (matrix estimator):
  mode(p,q)  cyc/pup     gain  cross-talk
    1,1          0.7   0.9877      0.0049
    2,2          1.4   0.9908      0.0028
    4,4          2.8   0.9884      0.0065
    8,8          5.7   0.9887      0.0049
   16,16        11.3   0.9886      0.0090
   24,24        17.0   0.9888      0.0158
   32,32        22.6   0.9889      0.0284
   48,48        33.9   0.9895      0.0563
   64,64        45.3   0.9909      0.0048
   80,80        56.6   0.9848      0.0043
   96,96        67.9   0.9639      0.0106
   48,0         24.0   0.9918      0.0621
Stage E -- differential rows (matrix estimator; actuator space, pm):
  base           deviation          gain   resid pm     corr
  flat           single 10nm      0.9916        2.2   1.0000
  flat           random 10nm      0.9893      169.2   0.9999
  random 16nm    single 10nm      0.9969        2.7   0.9999
  random 16nm    random 10nm      0.9894      231.3   0.9998
Stage E break ladder -- single 10nm differential vs base rms (lens):
  base rms       gain   floor pm     corr
      30 nm   1.0013        4.5   0.9997
      60 nm   1.0103        9.3   0.9987
     120 nm   1.0258       19.3   0.9946
     240 nm   1.9271      394.0   0.6872
     480 nm   1.5849      419.1   0.5661

wrote lens_report.txt + lens.mat + figures in /home/dcr/dev/MACOS_resources_wt_oap/mmacos/templates/40_benches/tg_psi_dm96_oap/runs/lens
```

### tag vlens (untuned tail -- item 2)

```
Stage MATRIX -- measured response matrix dw/da (lens rig, step 8, sign same):
  flat-DM null: 9.1111 nm rms surface (9111.1 pm) -- the same-plane-fold arm difference
  J: 64 states, 3268 columns (lit), window 41 px, reg lambda 1.00e-03 (of median col energy)
  Stage C single actuator @150nm: gain 0.9928, off-target floor 21.7 pm
Stage D -- modal transfer (matrix estimator):
  mode(p,q)  cyc/pup     gain  cross-talk
    1,1          0.7   1.0029      0.0019
    2,2          1.4   1.0024      0.0007
    4,4          2.8   1.0017      0.0025
    8,8          5.7   1.0019      0.0020
   16,16        11.3   1.0021      0.0030
   24,24        17.0   1.0019      0.0051
   32,32        22.6   1.0017      0.0093
   48,48        33.9   1.0002      0.0168
   64,64        45.3   0.9944      0.0024
   80,80        56.6   0.9695      0.0070
   96,96        67.9   0.9256      0.0192
   48,0         24.0   1.0006      0.0184
Stage E -- differential rows (matrix estimator; actuator space, pm):
  base           deviation          gain   resid pm     corr
  flat           single 10nm      0.9939        2.5   0.9999
  flat           random 10nm      0.9948      134.7   0.9999
  random 16nm    single 10nm      0.9922        2.7   0.9999
  random 16nm    random 10nm      0.9948      141.3   0.9999
Stage E break ladder -- single 10nm differential vs base rms (lens):
  base rms       gain   floor pm     corr
      30 nm   0.9906        3.0   0.9999
      60 nm   0.9878        3.7   0.9998
     120 nm  -9.0161     2427.7  -0.6843
     240 nm   0.5771      737.5   0.1364
     480 nm  -2.9720      861.4  -0.7144

wrote vlens_report.txt + vlens.mat + figures in /home/dcr/dev/MACOS_resources_wt_oap/mmacos/templates/40_benches/tg_psi_dm96_oap/runs/vlens
```

## Addendum 2026-09-12: 662e76a received -- items 1-6 closed; the D5 reframe; one modelling run before it is a design number

Items 1-6 are closed as asked (evidence committed as pruned results;
tail keyed by optics with a null WARN; Voronoi byte-identical -> not
truncation; bright/dark split -> the dark 25 % are dark, not
regularized; the D1 picture; the non-vacuity stated on the OAP and the
lens line labelled vacuous; calib_surface 'base' + the wrap guard).
`origin/tg96-oap` merges CLEAN into dev-candidate as it stands today
(V1 included).  Good work.

**The D5 reframe is accepted as measured, not yet as physics.**  The
picture says the dark columns are a bowed vertical BAND through the
pupil centre (d1_picture.png, column norm ~0 in the band, ~5.5 outside),
i.e. a fringe-VISIBILITY null: in the band the test arm's polarization
after the two folds is orthogonal to the reference arm's at the
analyzer.  A perfect conductor has r_s = -1, r_p = +1 at every angle
(zero retardance between s and p in the ray-following basis), and BARE
aluminium at 5-9 deg is within a degree or two of that -- so if the
null is the conductor's, bare metal has it too, and a coating whose
retardance is only a few degrees at 9 deg filling it from 0.05 to 0.82
means the null is a knife-edge in the polarization train, not a fold
cost at all.  Two runs settle which, and both are one command each:

1. **Bare aluminium** on both OAPs: `coat_set` with the Al index at
   632.8 nm (n = 1.373, k = 7.62, Rakic 1998 -- the
   `pol_external_anchor` tool carries van Harten's table if you prefer
   theirs).  Same rows as D5.  If the band persists, the uncoated D3
   numbers ARE the bare-metal numbers and the coating is a design
   parameter (state its retardance at both OAP AOIs and its chromatic
   slope); if it fills in, the "ideal reflector" idiom was the
   idealization and D3 is retired.
2. **The Jones pupil of the test arm** for the three cases (ideal /
   bare Al / protected Al) with `macos.jones_pupil` + `pol_maps` (the
   Phase-2 binding-side tools; double-pole basis): the retardance map
   and the per-pixel fringe visibility.  Report the visibility at the
   band centre and at the pupil edge for each case -- numbers, beside
   the dense-random gain.  That is the mechanism statement the report
   needs instead of "polarization artifact".

If the band is real physics, the cheaper fix is in the train, not the
mirror: re-solve the test-arm QWP azimuth for the folded arm (the
`tg_psi` runner solves the waveplate azimuths for the lens rig; the OAP
rig has a different arm rotation) before reaching for a coating.  Say
which you did.

**D7 (the loop).**  Blocked on Dave pushing MACOS_resources
dev-candidate (`dmg_loop` + `tDmgLoop`; also the V1 vector reading, so
`stage_loop_` now has four readings to mirror).  When it lands: merge
dev-candidate into tg96-oap (clean today), the stage exactly as
`zwfs_run` `stage_loop_`, seed 77, photons per MEASUREMENT, lens rig
first, the S11 table as the comparison; add the polarized-pair ZWFS row
(V1: 3 pm from 1.5e12 noise-only / 5.3e12 walk) beside S.

**Merge order (my recommendation to Dave):** merge tg96-oap into
dev-candidate now -- it is additive (new template dir, `dmg_frame`'s
optional third output, `twyman_green` 'optics', tBench cases) and
conflict-free -- and do D7 on dev-candidate directly, so the loop code
is edited in one tree.
