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
