# Focal-surface work — the FEX/SXP curved-`iElt+1` radius fix, and the focal-surface finder

_TO, 2026-08-28.  Tasking: `macos/BRIEF_focal_surface.md` (Dave ruled the
fix YES 2026-08-28).  Evidence base: `macos/REPORT_wnom_cli_ab.md`.
Engine change: `macos_f90/tracesub.F` (+189/−3).  Drivers, logs and the
corpus scan: `~/dev/MACOS_sandbox/xp_tst/fs_fix/`.  All lengths mm._

---

## Verdict

The fix is in and every gate is green.  On the zoom fixture the corner
nominals drop **4.26/4.37/3.09/3.23e-5 → 7.33/7.96/11.47/9.92e-6 mm**,
matching the A/B report's V3 column to **≤1.4e-5 relative**, with the
centre field unchanged.  Every flat-`iElt+1` deck is **bit-identical** —
written sphere and OPD map, `0.000e+00` on all four measures.  FEX ≡ SXP
to `0.000e+00` everywhere.

Two things in the brief turned out to be **wrong, and the correction
matters**:

1. **keckFF's `Kr = −0.0200` does not exist.**  The engine reports
   `Kr = −1e22` (flat) for keckFF's element 44.  So do iris_dp
   (not −2516.07) and eac2_7seg (not −306).  Those values came from
   text-parsing the `.in`, which mis-indexes elements on several decks
   (§5).  There is no degenerate-Kr deck to protect against.
2. **The blast radius is 11 of 208, not 9 of 149, and every one of them
   is the JWST OTE focal sphere** `Kr = −3017.56` (§5).  The earlier list
   (iris_dp ×3, keckFF, e5hex2) was produced by the same mis-indexing.

One thing the brief did **not** anticipate, and the fix now guards
(§4): a `Surface= Flat` element may carry a small `KrElt`, so gating on
`|Kr|` alone — which is what the brief's spec literally says — mis-reads
a flat detector as curved.  Measured cost of the naive form: **22.25 mm**
of spurious EP radius on e5hex1.

---

## 1. What changed

`macos_f90/tracesub.F`.  New module procedure **`FEXConicLeg`**, called
by **both** `FEX` and `SXP` at the point where each computes its
plane-leg distance.  The plane leg is computed exactly as before and
handed in; the routine returns a refined distance and a flag, and the
caller only overwrites when the flag is set.

```
        IF (LFPSrf) THEN
          ...print element, its Kr, and the plane−surface gap...
          fexT = fexTs
        ELSE IF (LFPTyp) THEN     ...curved but not a conic-base type...
        ELSE IF (curved)  THEN    ...chief misses the surface...
        END IF
```

Downstream is untouched: the telecentric guard, the beam-footprint
autoswitch and the Rx-order flag all still see `zp` and behave as before.

The geometry.  With `prel = o − Vpt(iNxt)`, `itpsi = d·psi`,
`psitprel = psi·prel`, `mpr = |prel|²`, the conic
`|p|² + 2Kr(psi·p) + Kc(psi·p)² = 0` gives

```
  a = 1 + Kc·itpsi²
  b = 2(Kr·itpsi + Kc·itpsi·psitprel + d·prel)
  c = (2Kr + Kc·psitprel)·psitprel + mpr
```

— the same coefficients `ConSrf` builds, so the surface being
intersected is provably the one the trace uses.

**Five decisions worth stating, each with its reason:**

| Decision | Why |
|---|---|
| Solve locally, **not** via `ConSrf` | `ConSrf` picks its root by `\|L²−mpr\|` proximity with no flow-of-light sense — the trap that produced the CENTROID branch's far-root garbage.  Here the root **nearest the plane leg** is taken, which is well posed: the sag correction is `h²/2R` while the two roots are ~`2R` apart. |
| Gate on **`SrfType`**, not on `Kr` | A `Surface= Flat` element may declare a small `KrElt`.  §4 measures what a Kr-only gate costs. |
| Refine only the **conic-base family** (Conic/Aspheric/Monomial/Zernike/the four grid composites) | Toric, Anamorphic, Interpolated and UserDefined give `Kr` other meanings.  Those keep the plane leg and print why.  No corpus deck is in that class. |
| `Kr = 0` keeps the plane leg | A zero-radius conic is a degenerate **point**, not a surface — `ConSrf` kills every ray on such an element.  Without this guard the near-tangent root returns the foot of the perpendicular from the vertex (measured: 1.14e-4 mm off, on the jwst deck). |
| Base conic only; any grid/FreeForm **figure** on `iNxt` ignored | µm-class against a sub-mm sag.  Stated in the header comment, not silent. |

---

## 2. Primary gate — the zoom fixture

`templates/50_sensitivities/zoom_5x5/jwst_ote_designc`, model 128,
ng 63, stop 25, OPD at elt 27, fields = deck `ChfRayDir + [dx;dy;0]`
normalised at ±2.90888e-4 rad.  Leg B4 (the CLI order: fresh load, field,
stop re-aimed, FEX).  Fit conventions are `wstats.m`, unchanged from the
A/B report (`a4` = coefficient of `2ρ²−1`; `rms` divides by N).

| fld | pre rms | post rms | ratio | pre a4 | post a4 | pre EP rad | post EP rad |
|---|---|---|---|---|---|---|---|
| C  | 6.8464e-06 | 6.8464e-06 | 1.000 | −8.1878e-07 | −8.1878e-07 | −3017.549071 | −3017.549071 |
| UL | 4.2603e-05 | **7.3329e-06** | 5.810 | −8.7713e-05 | −1.2405e-05 | −3018.083186 | −3017.596198 |
| UR | 4.3669e-05 | **7.9558e-06** | 5.489 | −9.0005e-05 | −1.4822e-05 | −3018.093178 | −3017.606978 |
| LL | 3.0923e-05 | **1.1470e-05** | 2.696 | −6.1105e-05 | +1.3972e-05 | −3017.978623 | −3017.491020 |
| LR | 3.2269e-05 | **9.9154e-06** | 3.254 | −6.4476e-05 | +1.0718e-05 | −3017.989875 | −3017.501917 |

Against the A/B report's V3 target column
(7.3329 / 7.9558 / 11.470 / 9.9154 e-6) the post-fix rms agrees to
**5.9e-6 / 2.5e-6 / 1.4e-5 / 8.0e-7 relative**, and the radii to
**2.3e-6 / 2.2e-6 / 1.5e-8 / 3.0e-6 mm** (the report quoted 5 decimals,
which is the whole residual).  Radii land on the V2/V3 values, not V1.

**Centre field unchanged**: 6.8464458386e-06 → 6.8464458483e-06,
**1.4e-9 relative**, from a 1.8e-12 mm radius round-off.  Not bit-equal,
and it should not be: the surface leg is computed by a different
expression that happens to agree to 1e-12 when the sag is zero.

**Non-vacuity** is the pre-fix engine itself, rebuilt from
`git checkout macos_f90/tracesub.F` in this session — it reproduces the
A/B report's table-1 and V1 numbers digit for digit
(B1 C 6.8464e-06, UL 4.2616e-05, UR 4.3682e-05, LL 3.0832e-05,
LR 3.2212e-05).

The **residual sign flip** on `a4` at LL/LR (−6.1e-05 → +1.4e-05) is the
expected signature: the pre-fix number was dominated by the reference
error, the post-fix number is the deck's own field curvature, and the
two have opposite sign at those corners.  Focus removal barely moves the
post-fix corners (`rms_pttf` 3.70/3.31/5.06/4.93e-06 vs `rms` 7.33/7.96/
11.47/9.92e-06) — i.e. what is left is no longer focus-dominated.

---

## 3. Flat-`iElt+1` bit-identity

Pre vs post, same recipe, model 1024, ng 63, FEX at `nElt−1`, OPD there:

| deck | Kr(next) | \|Δrad\| | \|Δvpt\| | \|Δpsi\| | max\|ΔW\| | verdict |
|---|---|---|---|---|---|---|
| e5hex1 | −1e22 | 0.000e+00 | 0.000e+00 | 0.000e+00 | 0.000e+00 | **bit-identical** |
| SegDemo3conic | −1e22 | 0.000e+00 | 0.000e+00 | 0.000e+00 | 0.000e+00 | **bit-identical** |
| iris_dp_conic | −1e22 | 0.000e+00 | 0.000e+00 | 0.000e+00 | 0.000e+00 | **bit-identical** |
| keckFF | −1e22 | 0.000e+00 | 0.000e+00 | 0.000e+00 | 0.000e+00 | **bit-identical** |
| eac2_7seg | −1e22 | 0.000e+00 | 0.000e+00 | 0.000e+00 | 0.000e+00 | **bit-identical** |
| synthFlatKr | −1e22 | 0.000e+00 | 0.000e+00 | 0.000e+00 | 0.000e+00 | **bit-identical** |
| jwst_zoom | −3017.56 | 1.36e-12 | 0 | 0 | 3.64e-12 | curved — moves off axis |
| j18sc | −3017.56 | 4.55e-13 | 0 | 0 | 3.66e-12 | curved — moves off axis |
| j18mono | −3017.56 | 1.36e-12 | 0 | 0 | 3.66e-12 | curved — moves off axis |

The three curved decks are probed at their **nominal on-axis** field,
where the sag is zero — 1e-12 is the whole change, which is the
mechanism's own prediction and a second confirmation of it.  The zoom
deck's off-axis behaviour is §2.

**FEX ≡ SXP**: `|FEXrad − SXPrad| = 0.000e+00` on all nine decks.

`e5hex2` belongs to the flat class but is **not** gated: it carries no
`ApStop=` and `macos.stop` fails on it (`stop_info_set failed`, elements
1–4), so FEX is unreachable on that deck.  Noted, not chased.

---

## 4. The `SrfType` gate, and why `|Kr|` alone is not enough

`docs/macos-manual/examples/SegDemo.in` declares
`Surface= Flat` with `KrElt= 0d0`.  A Kr-only gate reads that as curved.

The Rx path cannot actually deliver it — **`msmacosio.inc:1338` and
`ChkDf2` (`iosub.inc:1105`) both force `KrElt = −1d22` whenever
`Surface= Flat`**, unconditionally.  Measured: a hand-edited e5hex1 with
`KrElt= −1.0E+02` on its Flat element 13 (`synth_flat_smallKr.in`) comes
back from `load_rx` with `Kr = −1e22`, and its FEX radius is identical to
the unmodified deck's.  So that fixture is vacuous as a control, and it
is reported here as the parser finding it is.

The API path **can** deliver it: `elt_kr` writes `KrElt` directly and
leaves `SrfType` alone — and that is precisely the path the Part-2
design layer will use.  Driving it from there:

| engine | e5hex1 FEX radius after `set_elt_kr(13, −100)` |
|---|---|
| shipped (SrfType gate) | **−2548.0018521245** — unchanged from baseline |
| Kr-only gating (built for this measurement) | **−2570.2492488727** |

**22.25 mm of spurious EP radius on a flat detector**, from the
implementation the brief's wording describes.  The counterfactual engine
was built by deleting the `SrfType_Flat` return and neutralising the
conic-base whitelist; it is reproducible from the same two edits.

Across the deck set, forcing `Kr = −100` and `Kr = 0` on the next element
leaves the FEX radius **exactly unchanged** on every flat deck, and falls
back to the plane leg on the curved ones.

---

## 5. Blast radius — and a correction to the earlier count

**The engine is the only reliable source here.**  Text-parsing `.in`
files mis-indexes elements: several decks declare an `nElt` that
disagrees with their `Element=` block count (`e5hex2` declares 24, the
engine reports 25; `eac2_7seg` declares 47 with 48 blocks), so the
element table shifts and the wrong element's `Kr` gets attributed.  That
is how keckFF, iris_dp and eac2_7seg landed on the affected list in the
A/B report and in `BRIEF_focal_surface.md`.

Re-done through the engine (`scan_engine.m` → `scan_engine.csv`: per
deck, `num_elt`, `get_elt_info(nElt−1).elt_id`, `get_elt_kr(nElt)`):

```
decks attempted                          : 406
  loaded ok                              : 396
  load failed / crashed the engine       : 10   (pre-existing, §7)
  with a FEX-writable EP at nElt-1       : 208
    iElt+1 CURVED -> MOVED BY THE FIX    :  11
    iElt+1 flat -> bit-identical         : 197
```

All 11 are the same JWST OTE focal sphere, `Kr = −3017.56`:

| nElt | deck |
|---|---|
| 28 | `mmacos/templates/50_sensitivities/zoom_5x5/jwst_ote_designc.in` |
| 28 | `mmacos/templates/50_sensitivities/zoom_5x5/dwdgrid_5zoom_5fov_jwst_ote_designc_grid.in` |
| 60 | `MACOS_sandbox/old_Rx/j18dcWithStop.in` |
| 7  | `MACOS_sandbox/old_Rx/j18mono.in` |
| 28 | `MACOS_sandbox/old_Rx/j18sc.in` |
| 28 | `MACOS_sandbox/xp_tst/jwst_ote_designc.in` + the 5 `ab_cli/cli_*.in` copies the A/B run generated |

Corrected engine facts for the three decks the brief flagged:

| deck | brief / A/B report said | engine says |
|---|---|---|
| keckFF | `Kr(next) = −0.0200` (a "20 µm artefact") | `Kr = −1e22`, **flat** — no fallback needed, nothing to flag |
| iris_dp ×3 | `Kr = −2516.07` | `Kr = −1e22`, **flat** |
| eac2_7seg | `Kr = −306` | `Kr = −1e22`, **flat** |
| e5hex2 | `Kr = −2548` | `Kr = −1e22`, **flat** |

Also worth knowing: **no corpus deck has a non-conic-base curved
`iElt+1`**, and **none pairs `Surface= Flat` with a small `Kr`** — both
guarded branches are dead in the corpus today and live only for the API
path and for future decks.

---

## 6. Regression

| Gate | Result |
|---|---|
| `makems.sh release` (ifx) | clean |
| `makegfortran.sh release` | clean |
| mmacos mex relink | clean |
| GMI regression | **6/6 pass**, `vs-ref = 0.000e+00` — bit-identical to the committed reference |
| `run_mmacos_tests.sh fast` | **419 pass, 2 fail** — both fails PRE-EXISTING, see §6.1 |
| `run_mmacos_tests.sh tFocalSurface` (new) | **11/11 pass** (Part 2) |
| zoom_5x5 committed template numbers | moved in the 5th digit; regenerated, §6.2 |

### 6.1 The two suite failures are pre-existing

`tRunCompare/test_poke_agreement_end_to_end` and
`tRunCompare/test_run_simulator_time_history` fail.  They are **not**
this fix's doing, and that is measured rather than argued: the pre-fix
engine (rebuilt from `git checkout macos_f90/tracesub.F`) fails the same
two tests with **bit-identical** actual values —

| assertion | pre-fix | post-fix |
|---|---|---|
| `max w_rel < 0.05` | 9.990072415433924e+02 | 9.990072415433924e+02 |
| `max abs(u) > 1e-7` | 1.523371601413485e-08 | 1.523371601413485e-08 |
| `rms_wfe_unc(2) > 500` | 6.286616661310122e+02 | 6.286616661310122e+02 |

The fixture (`e5mono` → `segment_rx` → `pie.in`) is a flat-`iElt+1` deck,
so the fix is a structural no-op on it, consistent with the identity.
Flagged for Dave; not chased.

### 6.2 The committed zoom_5x5 numbers moved — regenerated

The brief asked to confirm no quoted digit moves and to regen + explain
if one does.  **One did.**

**Attribution first, and it is measured.**  Re-running the committed
driver on the **pre-fix** engine reproduces the committed
`dwdx_..._sens_report.txt` **byte for byte** (zero diff).  So the
artifacts were not stale, and the shift is caused by this fix and
nothing else.

What moved, on the fixed engine:

| quantity | committed | regenerated | rel |
|---|---|---|---|
| Elt 5 `Rz` column RMS | 2.1582e-02 | 2.1583e-02 | 5e-5 |
| group / Elt 5 `Rz` ratio | 11515.6161 | 11515.3882 | 2e-5 |
| elt-4 (dead) null floor | 5.30e-05 | 5.29e-05 | 2e-3 |
| `cfg zUL` smallest sv | 3.577e-05 | 3.293e-05 | 8e-2 |
| segment-only `cond+` | 2.525e+16 | 6.247e+16 | 2.5x |

Everything the README quotes except the `Rz` row is unchanged.  The
mechanism is the one the brief predicted: the dw/dx columns are
DIFFERENCES against a reference re-found once per field and then held
across that field's pokes, so the 0.49 mm reference change very nearly
cancels and only ~1e-5 relative survives.  The last two rows look
alarming and are not: the smallest singular value of a rank-deficient
cond-1e10 matrix, and the condition number of a cond-1e16 one, are the
most perturbation-sensitive numbers in the report — an 8% move in
`sv_min` under a 1e-5 matrix perturbation is expected behaviour, not a
finding.  The elt-4 row is the deliberately-dead obscured element whose
columns are pure null floor.

**Done:** the five tracked dwdx artifacts were regenerated **in the
template directory** (a first regen in a scratch dir leaked sandbox
paths into the report header — discarded), the `_resume/` checkpoints
the run created were removed, and the README's `Rz` row was updated with
a dated note carrying the attribution.

**Not done, deliberately:** the `dwdz` / `dwdsurf` / `dwdgrid` artifacts
in the same directory.  Same mechanism, same size, but they are separate
multi-hour runs and the README quotes none of their numbers.  Flagged
here so it is a decision, not an omission.

---

## 7. Loose ends found while measuring (none are this fix's doing)

- **10 decks fail or crash the engine at LOAD.**  `load_rx` never calls
  FEX, so none of this is the fix's doing; recorded because the corpus
  scan had to route around it.  Hard SIGSEGV: `old_Rx/ape.in`,
  `old_Rx/mars.in`, `old_Rx/nngx.in`, `old_Rx/wfpc3.in`,
  `old_Rx/wideAngle.in`, `ZGD_test_files/tst_save_keys.in`.  Clean
  `load_rx failed`: `docs/macos-manual/examples/SegDemo.in`,
  `segmirmaker/test_in/xx.in`, `old_Rx/eac1_opt_met_vcs.in`,
  `old_Rx/HubbleTour.in`.  A *documented manual example* that no longer
  loads (`SegDemo.in`) is the one worth Dave's attention.
- **`e5hex2` cannot set a stop** (`stop_info_set failed` at elements
  1–4) and carries no `ApStop=`, so FEX/SXP are unreachable on it.
- **The FEX radius still doubles as the far-field PROPAGATION distance**,
  where a tangent *plane* target is genuinely right.  This slice moved
  the OPD reference.  If a physical-optics leg ever needs the plane
  distance separately, the two uses must be split.  Unchanged from the
  A/B report; flagged again because the fix makes them differ.
- **The legacy `zp_iEm1` leg is not a general substitute.**  It equals
  the surface leg on the jwst deck only because elt 26 is coincident
  with elt 28.

---

## 8. Re-running

| What | How |
|---|---|
| primary gate | `cd ~/dev/MACOS_sandbox/xp_tst/fs_fix && matlab -batch "gate_zoom('post')"` |
| multi-deck probe | `matlab -batch "gate_decks('post')"` |
| pre/post diff | `matlab -batch "cmp_gates"` (needs both `gate_*_pre.mat` and `_post.mat`) |
| pre-fix engine | `git checkout macos_f90/tracesub.F`, `makegfortran.sh release`, `rm mmacos/src/mmacos.mexa64 && make FC=gfortran`, then the `('pre')` legs |
| corpus scan | `./scan_all.sh` → `scan_engine.csv` (resumes past the crashers) |

**Build gotcha that cost a cycle here:** the mmacos `Makefile` does not
treat `libsmacos.a` as a dependency, so an engine rebuild alone leaves
the mex linked against the OLD library and the gates silently measure
the previous engine.  `rm src/mmacos.mexa64` before `make` after EVERY
engine change (memory `reference_build_test_workflow` says so; I skipped
it once and got a wrong non-vacuity number out of it).

---
---

# Part 2 — the focal-surface finder

_The fix references `iElt+1`'s DECLARED surface, so the deck must declare
the right one.  This measures where the system actually focuses._

## P2.0 Verdict

`macos.design.focal_surface` (`mmacos/design/src/focal_surface.m`) works
end to end and the loop closes: measure the best-focus surface → emit the
replacement geometry → on the emitted deck the FIXED engine's own FEX
radius lands on the **focus-nulling** radius, so the corner nominals fall
to the A/B report's V4 column **to ≤2e-4 relative**.  That column is the
one the report had to construct by hand; the deck now produces it.

The other deliverable answer is a negative one, and it is the honest
result: **the deck author was right.**  The fitted radius is
**3027.31 ± 14.19 mm** against the deck's declared **3017.56** — 0.7σ.
At this system's own ±1′ field, a ±1′ cloud simply cannot separate them.

## P2.1 Surface

`focal_surface(rx, ...)` — a plain function in `design/src`, the
`pupil_find` convention (that folder is on the path, not a package).

| option | meaning |
|---|---|
| `'fov_rad'` + `'grid'` `'NxM'`\|`'cross5'`, or `'fields'` K×2 | the field set; default `3x3` |
| `'configs'` | a `macos.design.configs_from_table` schedule, or `[]` |
| `'fit'` | **`'sphere'` \| `'plane'` — BY CHOICE.  The other model's residual is PRINTED, never substituted.** |
| `'stop_elt'` `'xp_elt'` `'ngridpts'` `'model_size'` `'init'` | as `run_sensitivities` |
| `'elts'` | which elements receive the fit; default = auto-detect |
| `'write'` | path for the emitted `.in` (`[]` = measure only); refuses to overwrite the input |
| `'dR'` `'verify'` `'verbose'` | slope calibration offsets, verify pass, printing |

**One implementation note that is not cosmetic:** the supervisors' config
machinery (`config_axis`) is `private` to `+macos` and therefore
unreachable from `design/src`.  Rather than duplicate its
snapshot/apply/undo, `focal_surface` **reloads the Rx for every
(config,field) point** — the reload IS the restore, so there is no
restore to get wrong, and no field can inherit the previous field's
written exit pupil.

## P2.2 The cloud point, and why it is fix-independent

Per point: `macos.fex`, then — vertex and axis HELD — calibrate
`d(a4)/d(rad)` over `'dR'` and solve for the radius that nulls `a4`.
The cloud point is that sphere's **centre**, `vpt − R_null·psi`: the
field's best-focus image point, the **medial** focus.

Two properties make this the right quantity, and both are measured:

* **It does not depend on where FEX started.**  Between the pre-fix and
  post-fix engines the FEX radius at a corner moves by **0.487 mm**,
  while `R_null` moves by at most **1.1e-4 mm** — a rejection factor of
  ~4400.  (`fs_pre.mat` / `fs_post.mat`.)
* **It does not depend on the stop order.**  The A/B report measured that
  the CLI order writes a radius up to 0.76 mm different from the runner
  order while holding the sphere CENTRE to 8e-5 mm.  The centre is what
  this returns.

The verify pass is not decoration: `a4` at the solved radius comes back
at ~1.0e-07 against a starting `a4` of up to 1.5e-05, so the linear solve
really does null the term it claims to.

## P2.3 The jwst measurement

Model 256, ng 63, stop 25, 3×3 grid at ±2.90888e-4 rad, sphere fit:

```
  centre = [0.475155 1.147528 -1057.300905]   radius = 3027.308461
  1-sigma: centre [2.68e-01 1.43e+00 1.41e+01]  radius 1.419e+01   cond 1.399e+09
  deck declared Kr = -3017.560611  ->  fitted |R| = 3027.308461  (delta -9.748)
  residual rms 1.696e-03   max 2.399e-03   (9 points)
  the OTHER model (plane), for comparison only: rms 1.619e-01  max 3.237e-01
  cloud extent [76.8152 76.3171 9.1089]
```

Read it in this order:

* **The surface is genuinely curved.**  The plane fit is **95×** worse
  (1.62e-1 vs 1.70e-3 mm).  `'sphere'` was still a caller CHOICE — that
  is the rule — but it is the defensible one here.
* **The fitted radius is NOT significantly different from the deck's.**
  9.75 mm delta against a 14.19 mm 1σ.  Reporting the delta without the
  σ would be the easy mistake and it would be wrong.
* **Why σ is that large, and why more numerical care would not help.**  A
  ±1′ cloud is a ~110 mm patch of a ~3000 mm sphere, so the sag across it
  is only ~0.49 mm.  The fit residual (1.70e-3 mm) is real structure —
  the astigmatic tangential/sagittal split and higher-order field
  curvature, which no single sphere can absorb — so the radius inherits
  an amplification of about `R/sag ≈ 6000`.  Tightening the `a4` solve
  cannot move this; only a wider or denser field can, and ±1′ **is** this
  system's field.
* **The same amplification explains a number that looks inconsistent and
  is not:** the pre-fix and post-fix fitted radii differ by 0.578 mm
  although no cloud point moved more than 1.1e-4 mm.  6000 × 1e-4 ≈ 0.6.
* **`cloud extent [76.8 76.3 9.1]`** — the 9.1 mm of "z" is mostly the
  focal surface's 5.7° tilt against the global frame (77 mm × sin 5.7° =
  7.6 mm), not sag.  Do not read it as depth.
* **Configuration motion: not measured here.**  The 3×3 run carried no
  zoom schedule (`config_motion = 0` is the trivial one-config answer).
  The plumbing is in and gated; running the 5-zoom schedule is a
  ~5× longer scan and is left as the obvious next call.

## P2.4 Emission

Auto-detection finds **elements 26 and 28** — the focal Return and the
detector Reference, which share the declared sphere — and correctly
excludes element 27, the exit pupil, whose radius differs.  The emitted
file is a NEW `.in`; only `KrElt/KcElt/psiElt/VptElt/RptElt` inside those
two blocks change, byte-for-byte everywhere else:

```
-   KrElt=  -3017.560611d0
+   KrElt=  -3.0273084612E+03
-  psiElt=   1.790323563D-02  9.857893159D-02 -9.949681746D-01
+  psiElt=   1.8144115136E-02 1.0078031896E-01 -9.9474324245E-01
```

The normal's **hemisphere is COPIED** from the element's existing `psi`
(dot product 0.9998), never re-derived — the `pupil_find` psi-hemisphere
defect is the cautionary tale.  The `Kr` SIGN convention is copied too.

## P2.5 The verification loop — the point of the exercise

Same recipe as the Part-1 gate, on the emitted deck, fixed engine:

| fld | plane-leg (pre-fix) | V3 declared surface | **emitted** | V4 a4-null | emit/V3 | emit/V4 |
|---|---|---|---|---|---|---|
| C  | 6.8464e-06 | 6.8464e-06 | 6.8355e-06 | 6.8360e-06 | 0.9984 | **0.9999** |
| UL | 4.2603e-05 | 7.3329e-06 | **4.2196e-06** | 4.2186e-06 | 0.5754 | **1.0002** |
| UR | 4.3669e-05 | 7.9558e-06 | **3.4374e-06** | 3.4378e-06 | 0.4321 | **0.9999** |
| LL | 3.0923e-05 | 1.1470e-05 | **9.2542e-06** | 9.2545e-06 | 0.8068 | **1.0000** |
| LR | 3.2269e-05 | 9.9154e-06 | **8.4488e-06** | 8.4469e-06 | 0.8521 | **1.0002** |

The chain closes to **≤2e-4 relative** at every field.  Corners fall a
further 0.43–0.85× below the declared-surface result and the centre is
unchanged.  The tell that the focus term is actually gone is that
removing focus now buys almost nothing: at UL `rms_ptt` = 3.6890e-06 and
`rms_pttf` = 3.6868e-06, a 0.06 % difference (pre-fix the same field went
4.2560e-05 → 3.8049e-06, an 11× drop).  What remains is real aberration
plus the astigmatic split a single sphere cannot absorb.

Total journey at UR: **4.3669e-05 → 3.4374e-06**, a factor **12.7**, of
which 5.5× is the engine fix and 2.3× is declaring the measured surface.

## P2.6 Gates

`mmacos/tests/tFocalSurface.m` — **11/11 pass**.  One scan in
`TestClassSetup`; the methods cover: cloud completeness + the verify
residual; a tight regression pin on all nine `R_null`; the A/B V4
cross-check **with its tolerance explained** (see below); fix-
independence; sphere-vs-plane; the fitted radius inside its own σ; both
focal elements found and elt 27 excluded; psi-hemisphere and Kr-sign
copying; the emitted deck loading with its FEX radius following the fit;
the emitted deck nulling the corner focus; and the overwrite refusal.

**One test I had to fix before believing it.**  The brief's gate (1) asks
that the cloud reproduce the report's V4 radii "~1e-3 relative".  At
3017 mm that tolerance is **3 mm** — and the runner-vs-CLI stop order
offsets the radius by a systematic **±0.7 mm**, so the assertion passes
while being nearly vacuous, and my first draft additionally compared
against a mixed row of report and measured values.  It is now split: a
**tight** pin (1e-7) on the measured radii, and a separate coarse test
that asserts the offset is present, bounded by its documented 0.76 mm
class, and carries the systematic **sign pattern** (−y corners long, +y
corners short) rather than being scatter.

Not built: the synthetic flat-field control the brief lists as gate (3).
No corpus deck has a genuinely flat focal surface at `nElt` **and** a
usable stop, and a hand-built one would test the fitter against data I
also authored.  The plane branch is exercised instead by the `other`
model on every run (P2.3) and by `resolve_fields_`/emission unit paths.
Stated as a gap rather than papered over.

