## SECTION 10: Worked Examples

Four complete workflows, each a runner that exists in the companion
repository `MACOS_resources` and prints its own table.  They are run from
MATLAB with the mmacos bindings (Section 9 describes the engine they call);
every one reads a single parameter file, writes `<tag>_<stage>.txt` beside
itself, and keeps that text as its record, so the test of this section is to
run the command and compare.  The numbers quoted here were produced on
2026-10-17 from the committed sources.

```matlab
run('<path-to>/mmacos/mmacos_setup.m');     % once per session: the +macos package and the mex
```

### 10.1 A sensitivity workflow and an error budget (`e5hex1_budget`)

`templates/50_sensitivities/e5hex1_budget/run_budget_e5hex1.m` harvests the
rigid-body sensitivity Jacobian of the public seven-segment telescope
`e5hex1.in` — `macos.dw_dx_multi`, five field points, six degrees of freedom
per optic, model size 128 — and turns it into an optical error budget after
Redding 2014: a properly formed RSS budget is a covariance analysis, and a
Monte-Carlo draw gives the same answer.  The error groups, their
distributions and the control set are a table (`budget_spec_e5hex1.m`): the
seven primary segments and the secondary are controlled in six degrees of
freedom each, 48 actuators, by the pseudo-inverse controller; a wavefront
sensing row and one actuator row per control group are added.

```matlab
cd <path-to>/mmacos/templates/50_sensitivities/e5hex1_budget
run_budget_e5hex1
```

The run takes a few minutes and leaves the roll-up table
(`e5hex1_budget.csv`), the text report, a per-group bar chart and two
workbooks, one of them (`e5hex1_wfe_budget.xlsx`) a live-formula calculator
in the HabEx form whose white cells are the inputs.

From the report (`e5hex1_budget_report.txt`), with every error rigid-body on
this deck:

| line | WFE, nm rms |
|---|---|
| uncontrolled (pre-WFC) total | 317.4 |
| post-WFC total, covariance | 3.156 |
| post-WFC total, Monte-Carlo, 1000 draws with the controller, sensing and actuator draws in the loop | 3.179 |
| wavefront-sensing row (3 nm white noise, of which only the controllable part leaks: 48 of 51 225 rays) | 0.092 |

The actuator rows dominate the corrected budget (primary piston 2.60 nm,
primary tilt 1.28 nm, secondary tilt 1.25 nm).  The covariance and the
Monte-Carlo totals agree to 0.7 %, the point of the method.  The 2026-10-17
run reproduced the record of the template's README and of the regression
test `tRunBudget/test_wfc_stage_cov_equals_mc` to every printed digit.

### 10.2 A telescope design: the long-slit TMA front end (`tma_longslit`)

`templates/10_telescopes/tma_longslit/tma_longslit_run.m` designs a
three-mirror telescope that feeds a straight slit telecentrically with a
bounded cone, the SBG VSWIR zig-zag form: a weak concave M1, a small convex
M2 that is the stop, M3 carrying the power, the slit 300 mm behind M3.  The
default instrument is the dyson5 3k front end: focal length 330 mm, aperture
183 mm, F/1.8, a ±4.7° strip onto a 54 mm slit.  Two stages run in seconds:

```matlab
addpath('<path-to>/mmacos/templates/10_telescopes/tma_longslit');
OUT = tma_longslit_run({'first_order','section'});
```

`first_order` is the paraxial zig-zag; `section` emits the three mirrors as
off-axis conic sections (the first-order result forces this: a tilted sphere
cannot supply the same power in and across the fold plane), traces the strip
through the engine on a 41-point grid and measures it.  Its table, from
`tls_section.txt`:

| quantity | measured | spec |
|---|---|---|
| plate scale, local / to the edge of the strip | 329.99 / 329.87 mm | 330.0 mm |
| image of the ±4.7° strip | 54.24 mm | 54.0 mm slit |
| F/# at the slit, along / across the strip, every field | 1.92–1.93 / 1.95–1.96 | 1.8 cone bound |
| chief-ray angle to the slit normal at the strip edge (telecentricity) | 0.093° | — |
| every ray passes | 1.000 | — |
| worst clearance, footprint plus a 5 mm mount (M3 → slit leg against M2) | +10.2 mm | > 0 |

The rms spot at this stage is millimetres (2230 µm at the strip edge): the
conic sections only set the geometry.  The `figure` stage, an engine-traced
Levenberg–Marquardt ladder of 15–80 minutes a rung, brings the figure in; the
`e2e` stage joins the best rung to the spectrometer of Section 10.3.  The
2026-10-17 run reproduced the 2026-10-07 record to the printed digit; the
only moving column was the chief-ray miss of the M2 pole, at 1e-7 µm.

### 10.3 A spectrometer design: the Dyson challenge (`dyson5`)

`challenges/dyson5/dyson5_run.m` designs an EMIT-class push-broom imaging
spectrometer from a public specification: F/1.8 in air-equivalent image
space, a 3000 × 500 pixel array at 18 µm (a 54 mm slit, 9 mm of spectrum),
380–2500 nm, smile and keystone under 0.1 pixel.  Stage `s0` is engine-free
and runs in under a minute; it sets the classical concentric Dyson at the
spec and finds the scaling law the rest of the design obeys.

```matlab
addpath('<path-to>/mmacos/challenges/dyson5');
OUT = dyson5_run();                      % stage s0 at the spec (then s1, the deck)
```

From `dyson5_s0_scaling.txt` (slit centre offset 6 mm, the slit corner at
27.66 mm from the axis):

- At the Dyson condition the blur grows as the fourth power of the field
  height (fitted exponent 4.17, a fifth-order residual) and falls as the cube
  of the block radius (fitted −3.18): blur ~ h⁴ / r³.
- At the slit corner the 0.25-pixel blur budget (4.5 µm rms) is met only at
  a block radius of 208 mm or more: grating radius 671 mm, a 257 mm grating.
  The verdict is to report the law and not add elements yet; the later
  stages (the meniscus, the fold prism, the telescope join) are the beats
  recorded in the challenge's README.

The stored record of this stage was taken on 2026-09-30 with the slit centre
offset then in force, 8 mm (corner at 28.16 mm, 213 mm required); the law and
its exponents are the same.

### 10.4 A gauge bench: COPHI on the TG96 interferometer (`cophi_dm96`)

`templates/40_benches/cophi_dm96/cophi_run.m` builds the heterodyne camera
interferometer as a deformable-mirror surface gauge on the bench of record of
the TG96 campaign: the reference flat's seat is taken by a traced, fiber-fed
local-oscillator arm, so there is no analyzer and no quarter-wave plate and
the polarization gain is absent by construction.  The bench stage solves the
clearances, builds and emits both arms as prescriptions (`<tag>_test.in`,
`<tag>_ref.in`) and registers them at the camera: pixel pitch, the chief-ray
pixel, the magnifications from the DM and from the seat, the local
oscillator's footprint over the test arm's pupil image and its own wavefront
there, and the null phase on the flat DM.

```matlab
cd <path-to>/mmacos/templates/40_benches/cophi_dm96
OUT = cophi_run('bench.optics','oap','tag','bench_oap','stages',{'bench'});   % the mirror rig of record
```

The run takes about a minute at model size 1024 (a 385-point ray grid) and
writes `runs/bench_oap/bench_oap_cophi_report.txt`.  Its registration block:

| quantity | measured | pass criterion |
|---|---|---|
| chief-ray pixel, LO minus test arm | 0.000 px | < 1 px |
| LO footprint centroid offset from the test pupil image | 1.625 px, inside the LO's 33.9 px overfill | < overfill |
| LO coverage of the test pupil (77 877 px), flatness | 1.0000, 0.21 % rms | > 0.99 |
| null on the flat DM, interference phase as surface after piston | 0.782 nm rms (tilt 0.647, defocus 0.169, rest 0.405) | < 1 nm |
| fringe contrast over the test pupil | 0.9817 median | — |
| magnification, DM to camera / seat to camera | 10.5576 / 10.5600 | — |

The 1.6 px is the test pupil (the image of the DM stop) sitting off the chief
ray, a property of the test arm, and the 0.65 nm of tilt in the null is the
LO collimator's beam direction against the parabola's; both are fixed maps the
frozen null subtracts, as on the record's interferometer (whose mirror-rig
null is 12.9 nm).  The 2026-10-17 run reproduced the README's record row of
the mirror rig to every printed digit.
