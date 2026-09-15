# A DM surface gauge inside a coronagraph -- the options

CCL for Dave, 2026-09-15.  The ask: the CTB example and the e2e6m space
version show the packaging constraints; one option is a normal-incidence
gauge face-on to the DMs, another comes in at an angle for clearance using
the ZWFS and/or PDI.  Look at the CTB layout and develop a list of options.

Every number below comes from the two committed decks through
`demo_session/gauge_in_coro_clearance.py` (pure Python; no MATLAB run, TO
owns the box today).  The tool walks the chief ray through each deck so the
off-axis parabolas are placed at their POLES (the vertices in the decks are
the parent parabolas' -- TO's reflective item 1 lesson; the same mistake
here would have put OAP1 758 mm from DM1 instead of 600), then scans a gauge
beam about each DM's normal.  Nothing was re-solved; the coronagraphs are as
committed.

## 1. The two packages

| | CTB bench (`bench_ctb/ctb_dcr.in`) | e2e6m space relay (`e2e6m_r2/r1_seg_d040_full.in`) |
|---|---|---|
| DM clear radius / beam radius | 22.5 / 21.4 mm | 30 / 23.75 mm (47.5 mm beam) |
| actuators, pitch | 32 across, 0.67 mm | 32 across, 1.48 mm |
| science incidence on each DM | 4.9 / 5.0 deg | 6.0 / 6.0 deg |
| OAP1 -> DM1 / DM1 -> DM2 / DM2 -> OAP2 | 600 / 500 / 700 mm | 250 / 400 / 200 mm |
| OAP body radius (assumed here) | 75 mm (README) | 60 mm (2x the beam; not in the deck) |
| out-of-plane room | the whole bench is in one plane (bodies span +/-20 mm in z): above and below the table is EMPTY | a 3-D package; the telescope's M3 -> OAP1 beam and the M2 -> M3 beam pass within 0.5 m of the DM pair |
| DM pupil | DM1 is the pupil; DM2 is 500 mm downstream of it | same, 400 mm |

Rules used: mounts +8 mm on every body, clearance >= 25 mm (the gauge
bench's rule), gauge beam = the DM's clear radius, the gauge's IN and OUT
legs at +theta and -theta about the DM normal in the plane of azimuth phi.
Two numbers per case: the gauge BEAM against every body over 700 mm, and the
gauge's first OPTIC (a body of beam radius + mount) at 300 or 500 mm from the
DM against every body and every science beam.  Light crossing light is not a
collision; hardware in a beam is.

## 2. The scan (best azimuth at each angle; negative = interference)

**CTB, DM1** (DM2 is mirror-symmetric):

| theta | beam clearance | binds | optic at 300 mm | optic at 500 mm |
|---|---|---|---|---|
| 0 (face-on) | **-53** | OAP1 body | -25 (in the OAP1 -> DM1 beam) | -8 (DM2) |
| 6 | -24 | OAP1 | -10 | +16 |
| 8 | -7 | OAP1 | -2 | +30 |
| 10 | +11 | OAP1 | +7 | +41 |
| 12 | +30 | OAP1 | +16 | +54 |
| **15** | **+58** | OAP1 | **+31** | +76 |
| 20 | +106 | OAP1 | +55 | +115 |

The best azimuth is 90 deg at every angle: OUT OF THE BENCH PLANE.  Face-on
is blocked twice -- by OAP1's 150 mm body, which sits 5 deg off the DM normal
at 600 mm (lateral 52 mm against 104 needed), and by DM2 at 5 deg the other
way (44 mm against 52).  A gauge that comes in from above the table at 15 deg
and leaves below it clears everything by 31 mm with its first optic 300 mm
from the DM; at 10 deg it needs the optic at 500 mm.  The lever if 15 deg is
too steep is the DM fold angle itself: the record's 5 deg puts both OAP1 and
DM2 right on the normal.

**e2e6m space relay, DM1** (DM2 mirror-symmetric, binding on OAP2):

| theta | beam clearance | binds | optic at 300 mm | optic at 500 mm |
|---|---|---|---|---|
| 0 (face-on) | **-72** | OAP1 body | -40 (OAP1) | +26 |
| 10 | -48 | OAP1 | -23 | +59 |
| 15 | -29 | OAP1 | -7 | +90 |
| 20 | -9 | OAP1 | +13 | +122 |
| 25 | +6 | OAP2 | +33 | +145 (M2 -> M3 beam) |
| 30 | +19 | OAP2 | +54 | +102 (M2 -> M3 beam) |

**No angle to 30 deg gives an external gauge beam 25 mm of clearance in the
space relay.**  OAP1 is 250 mm from DM1 at 6 deg off its normal, so its body
(assumed 60 mm) sits inside any beam that leaves the DM within ~25 deg of
the normal, and past 25 deg OAP2 takes over from the other side.  The
assumed OAP radius is the lever: at 45 mm the 25-deg beam clears by ~21 mm.
Either way an external gauge on the flight relay is a 25-35 deg periscope
folded through the telescope's own beams, and the relay would be designed
around it, not fitted with it.

## 3. The options

**A. Face-on external gauge (normal incidence).**  The cleanest reading: the
gauge deck's bench as is, the DM in the retro seat.  Does not fit either
package as committed: the DM fold puts OAP1 (and the other DM) on the
normal.  Fits only if the coronagraph's DM fold opens to ~10 deg AND the
DM -> OAP leg is long (CTB-class, 600 mm+); never in the space relay.  A
flip-in fold in the science leg during calibration is the same thing with
the science beam blocked -- fine on the ground, not in flight.

**B. External gauge out of the science plane, oblique (ZWFS / PDI).**  The
CTB's free direction.  A collimated gauge beam from above the table at
theta, out below it, the gauge's focuser + mask + camera on a riser
(300 mm at 15 deg = 78 mm above the table).  One gauge per DM, each its own
pupil, each measuring its DM alone (job 1, 2 and 3 per mirror, exactly the
deck's readings).  Costs, all small and all modelable in the existing
runner with the DM tilted: OPD sensitivity 2 cos(theta) (0.966 at 15 deg);
the footprint is an ellipse 1/cos(theta) long (3.5%) -- the pupil image
elongates and the matrix measured on the surface absorbs it; the DM's own
surface is seen foreshortened along one axis, so the actuator grid is not
square on the camera (again the matrix).  The self-referenced sensors fit
this best: one beam in, one out, no reference arm to route; the
interferometer's reference flat would sit on the riser too.  Non-common
path with the science beam entirely, which is the point: it measures the
mirror, not the coronagraph.  Stray light: the gauge wavelength must be out
of the science band (a long-pass at the science camera, or a 1.5 um gauge);
scatter from the DM at 2 theta into the science cone is the DM's BRDF times
a very small solid angle -- a number to compute, not a worry.  Does NOT fit
the space relay (section 2).

**C. Gauge injected into the science path, read at the first focus after
the DMs (the coronagraph's own light path, near-normal incidence).**  The
science beam already hits each DM at 5-6 deg.  Inject the gauge at the
source (CTB: a second fiber on a dichroic at the source; flight: the
pseudo-star / calibration source at the telescope focus) and read it at
Focus23 (CTB) or after OAP2 (space): a ZWFS or pinhole mask at that focus,
the pupil reimaged by OAP3 onto a gauge camera behind a dichroic pickoff at
the apodizer conjugate.  Nothing new enters the DM neighborhood; the pickoff
is one plate in the collimated beam after OAP3, or the FPM's own reflected
light if the FPM is reflective (the standard LOWFS seat).  What it measures:
the SUM of the two DMs at DM1's conjugate, with DM2 seen through 400-500 mm
of Fresnel propagation.  This is the reading the coronagraph's wavefront
control wants anyway; per-DM surfaces come out only as far as the two mirrors
are separable (next paragraph).  Adds to the science train: the dichroic (a
plate in the collimated beam: 0.03 wave W040-class for 2 mm at F/4, a
polarization term at the vector sensor), the gauge's non-common path
downstream of the pickoff (the reference-arm walk lesson: a piston, benign
for a DM servo, not for absolute field work).

**Separability of DM1 and DM2 in one conjugate, the number that decides
C.**  A phase sinusoid of period L on DM2 converts to amplitude at DM1's
plane by sin(pi lambda z / L^2).  CTB (0.67 mm pitch, 500 mm, 550 nm):
47% at the actuator Nyquist period, 93% at one cycle per actuator, 50% at
33 cycles across the beam.  So on the CTB the two mirrors separate at the
actuator scale (a complex-amplitude reading -- the vector pair with its clear
frame, the pinhole, the interferometer -- sees DM1 as phase and DM2's
actuator print-through mostly as amplitude) and do NOT separate below ~30
cycles per aperture, where both are phase.  Space relay (1.48 mm pitch,
400 mm, 500 nm): 7% at Nyquist; 50% only at 43 cycles across the beam,
beyond the actuator Nyquist of 16.  **In the flight relay a single-conjugate
reading cannot tell DM2 from DM1 at any controllable frequency.**  Per-DM
knowledge there needs either a second conjugate (reimage DM2 onto a second
camera behind the same pickoff -- one more lens, and the two complex fields
over-determine the two phase screens) or differential attribution (poke one
DM, the matrix column says which -- the gauge deck's multiplexed matrix
already carries both DMs' actuators as columns; the deck's "complex
amplitude for two-mirror control" slide is this).

**D. The coronagraph's own sensors as the gauge (LOWFS / dark-hole
estimation).**  The FPM-rejected-light ZWFS (LOWFS) and the science-camera
pairwise probing already measure the field the DMs make.  They are the
capture-and-hold instrument for the SUM at the science wavelength; they are
not a surface gauge (no absolute per-mirror figure, low order only for the
LOWFS, probe-limited for the dark hole).  Listed because option C at its
minimum IS the LOWFS with a picometer-class reading in its seat, and the
deck's numbers (photons per measurement, capture range) transfer to it
directly.

**E. Gauge through a dichroic AT the DM leg (normal incidence without a
free normal).**  A dichroic plate in the DM1 -> DM2 leg at 45 deg would let
a gauge beam reach DM1 face-on and return through the same plate.  Fits the
CTB's 500 mm leg (a 60 mm plate 200 mm from DM1; the scan's DM2/OAP1
bindings vanish because the gauge comes from the side).  Costs the science
beam a tilted plate at a pupil-adjacent station -- astigmatism and a
polarization term in the collimated beam, ghosts into the dark hole -- the
one place a coronagraph least wants glass.  Ground calibration only, if at
all; ruled out for flight.

**F. Two gauges time-multiplexed on one out-of-plane periscope.**  Option B
with one gauge head serving both DMs by a fold that switches legs (or two
heads sharing a laser and a camera through a fiber switch).  Halves the
hardware, serializes the two readings; each DM's servo cycle is then shared.
Worth it only if the riser space is the constraint.

**G. Open the coronagraph's DM fold for the gauge (a relay design knob).**
The record's 5-6 deg DM incidence is what puts OAP1 and the other DM on the
normal.  At 10-12 deg incidence on the CTB the face-on beam clears DM2 and
OAP1 without a riser; the price is the relay's own polarization and the
off-axis parabolas' fold angles growing with it (the reflective front end's
"no fold coma once fed at its focus" applies to the relay's OAPs too, so the
price is polarization, not blur).  For the space relay the DM -> OAP legs
(200-250 mm) are the constraint, not the angle; lengthening them is a
volume trade against the 8 m shroud.

## 4. What I would model next (in order)

1. **B on the CTB:** the gauge deck's bench with the DM at 15 deg incidence
   (out of plane), S / V / P through the existing runner (`bench.DM_TILT`
   is a one-line addition: the retro seat tilted, the return leg folded);
   the rows on the 30 nm surface, the elliptical-footprint effect through
   the matrix, the photons.  One model-1024 run per reading.
2. **C's separability, measured not estimated:** the two-DM complex-amplitude
   reading on the CTB deck (both DMs as grid surfaces, the vector pair with
   its clear frame at Focus23 reimaged by OAP3), the multiplexed matrix over
   both DMs' actuators, the cross-talk between DM1 and DM2 columns vs spatial
   frequency.  Confirms or refutes the 33-cycle crossover.
3. **The space relay's OAP body sizes from the design** (not assumed) and a
   re-scan; then the second-conjugate variant of C (reimage DM2), since B
   does not fit there.
4. The stray-light number for B (DM BRDF at 2 theta into the science
   acceptance) and the dichroic's contrast cost for C -- both are inputs
   the coronagraph side has to accept.

## 5. Tool and assumptions

`demo_session/gauge_in_coro_clearance.py` (run with pymacos's venv python
for numpy).  Assumptions to replace with design values: OAP body radii (75
CTB from its README; 60 space, not in the deck); mask bodies 25 mm; M2/M3
bodies 400 mm at their vertices (far from the DMs; only their beams matter
and those are placed from the chief walk); the gauge beam equals the DM's
clear radius; the gauge's first optic is a sphere of that radius + mount.
The scan does not yet place the gauge's own tail (focuser, mask, camera)
beyond the first optic, nor the interferometer's reference arm; both go on
the riser and clear by construction on the CTB, and are moot in the space
relay until the OAP sizes are real.
