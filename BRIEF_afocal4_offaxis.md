# BRIEF: afocal4 OFF-AXIS SEEDED — is the FAMILY the wall? (24 h box)

**For TO, cold start.  Read first:** `MACOS_res_dev/mmacos/challenges/
afocal4/` — `RESULTS.md` §§ CLEARING (C.0–C.9) + DESCENT (D.1–D.7),
`descent/README.md`, and the off-axis probe wrap (Task 0 below).
Additive-only: new dir `challenges/afocal4/offaxis/`; nothing existing
overwritten.  Commits local on `MACOS_res_dev` dev-candidate, **no push**.
Box: **24 h from handoff** (Dave 2026-09-01).

## The ruling this brief executes

The descent found 71 nm unreachable in the coaxial family at any mirror
count (48× off at N=7, wavefront-only, every DOF free).  The rigid-body
off-axis probe then measured — correctly and honestly — only that **the
coaxial point is a local optimum under rigid-body perturbation**: with
±15° / ±300 mm available, the solver used 0.92° / 3.1 mm at most, the
wavefront-only arm landed in a *worse* basin (5157 vs 3842 nm) and broke
M to 1.03 % (the pre-flagged cheat).  That statement is about the
SOLVER'S BASIN, not the family.  Dave's ruling: **seed off-axis, don't
perturb there.  24 hours.**

## Task 0 — wrap the probe into the record (first hour)

- Append the probe as its own short record entry (RESULTS.md § OFFAXIS
  O.0 or DESCENT D.8 — your call, one findable place): the table, the
  two readings, and the explicit statement that "off-axis doesn't help"
  is NOT established — only "rigid-body perturbation from the coaxial
  point returns to it."  The M-cheat guard and the basin-scatter detail
  are worth their sentences; they are why the probe is trustworthy.
- Commit the signed-off **§S4b.4 correction** — Dave's sign-off is
  recorded at the foot of `BRIEF_afocal4_wall.md` (2026-08-31); it is
  still listed as held.  Close it.

## Task 1 — build a genuinely off-axis seed (the heart of the slice)

Seed the family, not a perturbation of the parent:

- **The classic seed is the off-axis Mersenne**: confocal parabola pair
  (f1/f2 = 30) used as OFF-AXIS SECTIONS — exact 30×, exact collimation
  BY CONSTRUCTION, unobscured by geometry.  Then a third (and if needed
  fourth) mirror for field correction and the interface-pupil work.
  Verify the seed's M and collimation at 1e-6-class BEFORE any metric
  is taken (the descent's rule: a cold closure is a specification, not
  a design — one descent probe traced M = 40.45 against paraxial 30).
- **Two build routes exist; use whichever closes first, name which:**
  (a) off-axis PARENT sections — parent conics + decentered apertures,
  the e2e6m-OAP idiom; (b) the design layer's tilted-component path —
  `freeform_unobscured` / `sz_tma` (3+n sphere+Zernike front end) and
  Telescope's Bauer unobscuring fold (`resolve_nmirror_fold_`, triggers
  on any nonzero `tilt_deg`).  Known trap on route (b):
  `realize_apertures` had a GLOBAL-XY→LOCAL-ApVec frame bug (sz_tma
  arc) — verify apertures land on the traced footprints before scoring.
- Multiple seeds (≥2 off-axis geometries) — the probe just demonstrated
  basin scatter; one seed is one basin.

## Task 2 — solve and score under the SAME requirement set

- Reuse the committed metric machinery: `afocal4_score` on the emitted
  deck, interface plane posed on the traced chief (recover `iface` from
  **zElt**, never vertex geometry — C.9 trap), wander at the refit
  plane with the as-emitted number beside it.
- Phase the requirements the way the descent priced them: first the
  **wavefront-only floor** (the family question — beside the coaxial
  3842 nm and the 71 nm target), then **add the interface-pupil ladder**
  and report the off-axis family's own pupil price (the coaxial price
  is 2.7× at 343 mm; is the off-axis family's cheaper?).
- **M and collimation printed on every result** (the probe's guard,
  now standing); converged rounds before anything is quoted (central
  differences — the wall slice's 17 %-low-gradient lesson moved every
  delivered number 19–36 %).

## Task 3 — realizability gates on anything quoted

- `afocal4_union` (body-vs-beam over the FIELD BOX, hull never a disk)
  + `afocal4_pack` on every design you'd put in a table.  An off-axis
  train trades obscuration for tilt — the failure mode moves from
  "collimator in the beam" to AOI growth and package girth; report max
  chief AOI per mirror against the 15° standing rule, and say where it
  is broken rather than hiding the design.
- Layout render per quoted design (review gate before numbers).

## Deliverables (end of box, whatever state)

1. RESULTS.md § OFFAXIS: probe wrap + seed geometries + the wavefront
   floor(s) beside coaxial 3842 / target 71, + the pupil price if
   reached.  Every resolved oddity written at resolution time.
2. Gate `tAfocal4Offaxis` (non-vacuous: fails a coaxial deck on the
   off-axis-seed check, passes the delivered seed).
3. Delivery log at the foot of this brief.  Partial results are fine —
   every rung is a finished floor; say which rung the box closed on.

**Success is an answer either way**: "off-axis seeded reaches X nm under
the same set" — X ≪ 3842 says the family was the wall; X ≈ 3842 says
the requirement set is the wall and the spec conversation goes back to
Dave/Mike with both families priced.

---

## DELIVERY LOG — TO, 2026-09-01

**Status at time of writing: Tasks 0, 1, 3 delivered; Task 2 wavefront arm in
flight (4 solves), pupil arm queued behind it.** Everything below is committed
locally on `MACOS_res_dev` dev-candidate (`4c32475`, `f7fb85d`, `01821fd`).
**Not pushed.** Additive: new `challenges/afocal4/offaxis/`; the one change
outside it is an optional `D.decenter` on `descent_build` that defaults to 0
and takes the identical code path when 0.

### Task 0 — the probe wrapped, the correction closed ✅

`RESULTS.md` § D.8 records the rigid-body probe with its limit stated
explicitly ("the coaxial point is a LOCAL OPTIMUM UNDER RIGID-BODY
PERTURBATION… it does **not** establish that an off-axis family cannot meet
the requirement set"), the M-cheat guard and the basin-scatter detail. The
signed-off **§ S4b.4 correction** is committed verbatim and the stale "wording
is with Dave" pointer retired.

### Task 1 — genuinely off-axis seeds ✅ (route (a), named)

**Route (a), off-axis parent sections with a decentered pupil.** The insight
that made it cheap: *an off-axis section is the same paraxial system with a
decentered pupil* — same powers, same spacings, so the afocal/magnification/
pupil closures and the entire descent machinery apply unchanged, and going
off-axis is a deck edit rather than a new closure.

Route (b) was not used: `Telescope.resolve_nmirror_fold_` is **focal-only**
(it resolves through `seidel_seed`'s `t_focus`), so it would need afocal
support added first. Recorded rather than attempted.

**Seed exactness, verified before any metric:** a confocal parabola pair is
afocal and 30× *for a beam entering anywhere on it*, and the engine agrees —
**collimation 0.00 µrad at decenters of 0, 0.6, 0.8 and 1.0 m**, and **exactly
0.00 nm** wavefront at true on-axis. This is the one seed in the arc that
needs no convergence before it can be trusted.

Geometries delivered (≥ 2, per the brief): Cassegrain and Gregorian Mersennes
at three focal lengths × five decenters; the committed 4-mirror and the descent
N = 5 / N = 7 rungs decentered; and a Cassegrain+Gregorian cascade.

### Task 2 — solve and score

**Delivered:** the bare-Mersenne floor with the pupil requirement dropped
(§ O.6) — 284× target at best — and, more importantly, **what that number
actually is** (§ O.6b).

**The finding of the slice.** Rung 2 does not remove power, and **60–99.7 % of
every design's wavefront variance in this study is power** — for an afocal
system, output collimation varying across the field, i.e. **field curvature**.
Confirmed by law: rung 2 runs 3666 → 18749 nm across the box, factor 5.11
where θ² predicts 5.90, and exactly 0 on axis. **It reverses § O.6's own
reading**: the off-axis Mersenne looks worst at rung 2 (20168 nm) and has the
**best rung-3 number in the study** (2778 nm), from two mirrors with no free
parameter — beating the coaxial *seven*-mirror floor. Its 20 µm is about
having two mirrors, not about being off-axis.

**In flight:** wavefront floors vs decenter at N = 4/5/7 (h = 0 control +
h = 0.55, 1.0) and a Mersenne-seeded arm in a different basin. Round-1
starts reproduce the descent's recorded starts exactly (10407 / 10775 / 7894),
which is the control the brief asked for. **The interface-pupil ladder is
queued behind them** — the pupil price is defined relative to the wavefront
floor, so it is sequenced, not skipped.

### Task 3 — realizability ✅

`run_offaxis_gates` scores `descent_require`'s full set plus two columns the
coaxial tables have no counterpart for: **chief AOI per mirror** against the
15° rule (decentering moves the chief off *every* parent axis at once, where a
fold spends its angle at one station) and the **parent radius** each element is
cut from — the off-axis family's own bill. Layout render per deck via
`macos.view_rx`, which draws the fitted off-axis *sections*. Validated on the
committed design: AOI 0.60/1.99/2.11/**12.84°**, union −79.9 mm, behind_m1
+1323 mm — all matching the record.

### Gate ✅

`tAfocal4Offaxis`, 6 tests, registered in `SUITE_FREEFORM`, **6/6**.
Non-vacuous by construction: a **coaxial deck fails** the off-axis-seed check,
and the widened-aperture trap is asserted to still lose rays.

### Things that were wrong and were fixed

1. **The measuring pass must REMOVE apertures, not widen them.** Widening to
   hold a decentered 1 m pupil asks the engine to intersect a |Kr| = 0.083 m
   parabola ~75 radii from its vertex: **0 of 1185 rays** survive, against
   1185 of 1185 with `ApType=None`.
2. **My loss count was reading the wrong flag.** `ray_hist`'s `ok` is
   *geometric* validity — an obscured ray keeps a valid intersection — so a
   fully vignetted beam read as lossless. Now `ok_trace .AND. ok_pass`.
3. **A test asserted a remembered mechanism on a deck where it does not fire**
   (the bare pair, not the N = 5 seed), and a nested-escaping slip meant its
   widening regex matched nothing. Both fixed; the trap is now measured, not
   recalled.
4. **The "Petzval sum" routine does not compute a Petzval sum** — MACOS emits
   `KrElt = −|R|` for every mirror and carries convexity in the geometry.
   Caught by the corpus itself (cass_greg and cass_cass return an identical
   11.600 while their defocus differs 47 %), renamed, and no Petzval claim made
   from it. **Corrected before it reached the record.**
5. **3841.8 nm is the N = 4 floor with tilt withheld**, not the seven-mirror
   floor (3424.2 nm, 48×). Two citations corrected at the point of use.

### For Dave — the question this changes

The slice was framed as *"is the FAMILY the wall?"* The measurements say the
binding term is neither the family nor the étendue (H = 2.18e-3 m·rad, ~4363
resolvable points per dimension; LSST's is 59× larger with three mirrors) but
**field curvature**, and that at 30× the final compressor mirror is
necessarily the smallest and most strongly curved in the train — the best
design in the study has the smallest curvature-magnitude sum of any deck
measured. One hypothesis about curing it (mixing Cassegrain and Gregorian
stages) was built and **refuted**, with its confound named.

**The open design question, and it is a new one:** whether a strongly curved
mirror of opposite sign *in the compressed beam near the exit* — where such a
mirror is at least small — can take the power term down. That is a different
move from adding mirrors to the front end, which is what the descent measured
and what bought 11 %.
