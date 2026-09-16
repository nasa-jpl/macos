# NOTE for CCL — close-out steps 1-4 landed (2026-09-16)

TO for CCL, via Dave. Branch `dev-candidate`, resources repo. Commits
`7aac1d8` (the staged patch) · `a935464` (step 2) · `2ffcd37` + `6b5141c`
(item 3) · `d43914a` + `03ebe27` (item 6) · `db7811d` (item 4). **Local, not
pushed** — Dave pushes.

## The two lines you asked for

**1. The wrap mechanism's verdict (item 2(a)).** The ladder break is the
**base reading WRAPPING**, and **both rigs do it at the same rung**. The
measured base rms saturates at 91 nm, and that value is analytic rather than
fitted: a four-step reading is the surface modulo λ/2 = 316.4 nm, so a base the
sensor cannot follow is uniform across that range with rms 316.4/√12 =
**91.34 nm**. Measured 90.2 / 92.1 / 91.9 (OAP) and 88.7 / 91.0 / 90.9 (lens)
at 120 / 240 / 480 nm; 30 and 60 nm track the command on both. **So "the
reflective rig has a smaller capture range" comes OFF the deck** — the capture
range is λ/2 of surface, set by the reading and not by the optics. The ladder's
own arithmetic is innocent (dA and dD agree pixel for pixel except at 2, 1 and
5 pixels of a ~10⁵ mask).

**2. The rows with the 10 mm plates (item 4).** `thk22_tail` reproduced its
winner — seed 75.2422 → **20.0938 nm**, 3.7× — and the advisory gate kept it,
so these rows describe a 20 nm-null bench (the earlier ones described a 75 nm
one). Stage C single actuator **0.9885**, off-target floor **40.0 pm**; modal
transfer ~1.00 from 0.7 to 45.3 cyc/pup and 0.9626 at 67.9, cross-talk ≤0.027;
differential rows 1.0023–1.0085, corr ≥0.9999; break ladder holds to 60 nm at a
**9.3 pm floor** and breaks at 120 — the same rung both rigs break at, so the
glass has not moved the capture range.

## The figures

`runs/stnoap/stnoap_stations.png` and `runs/stnlens/stnlens_stations.png`,
both **1800 × 560 px**. The OAP one is bit-identical in content to the
pre-patch `oapifol2` figure (0.00 pm flat, 626.43 pm on the 30 nm surface).

## Three things to know before they reach a slide

**The lens station figure's headline number is an OPEN defect.** It reads
**62 067 pm** where the OAP reads 626. Zero at flat on both rigs, √2 × the map
with structure — a lateral misregistration between the recovered map and the
engine field. It is the first lens station figure ever made, so not a
regression. The picture is sound (its first six panels are the tool's own
output); the seventh panel's number is not. **Do not put 626 and 62 067 side by
side** — that would assert a rig-to-rig quality difference that is not
established. `REPORT_bench_realism.md` §6.

**NEW, and it is slide-worthy on its own: the 10 mm plates read 13 % LOW.**
Measured base against commanded, 26.6/30 and 52.1/60, and in saturation 79.6
against the analytic 91.34 — one common factor moving all three. The record
2.6 mm bench tracks the command. The matrix absorbs a common scale, so this
costs **SNR and floor, not gain**, which is why every gain row looks healthy
and only the reading shows it. The old `max|h|` meter saturates at 1.00 from
120 nm up and could not have found it. §3.2.

**The tail gate now ENFORCES, and its criterion is RELATIVE to the seed**
(Dave's call, same day). It refuses a winner that reads worse than the
geometric seed it would fall back to, both measured through the same
estimator: `objwin3` ratio **0.1977 REFUSED**, `lens_tail` **0.9467
ACCEPTED**, `thk22_tail` **0.9930 ACCEPTED**. The two-leg test passes for the
first time. Two absolute thresholds had failed before it — the point-sample
measure inverted the verdict, and the lattice measure reads 8–19 % below the
battery (not the regularizer; the sweep is flat to ~1 %). A ratio of two
identically-estimated quantities cancels the bias instead of calibrating it.
Margin 0.047 at the tightest leg with exactly zero run-to-run spread. `4.9`.

## Queue state

Steps 1–4 of `BRIEF_to_restart` are done, item 3 included (the gate closed after Dave ruled on the criterion). **Step 5 (item 7, the coronagraph
field servo, 2–3 days) is not started**; step 6 (the CTB regeneration at the
intended beam) remains gated on Dave's word and on 1–5.
