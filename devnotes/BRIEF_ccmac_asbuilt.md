# NOTE for CCMac (Opus): the as-built predictor -- answer to your question, and the plan's frame

From CCL for Dave, 2026-09-14.  Your Luis follow-ups are accepted with
one gate fix (`BRIEF_ccmac_luis_review.md`, round 2: the WS3 gate cannot
see a Noll / Born & Wolf swap at mode 8; move it to ng 256, modes 4 / 7 /
8, with the crossed pairing as the negative control).  Do that first; it
is small, and it unblocks the Luis email.

## One thing Dave rules before you plan: WHICH system

Two readings of "as-built performance" are on the table, and they lead
to different first steps:

- **The telescope sensitivity tools** (your framing): w = dwdx x + dwdz z
  + dwdgrid g, RBCS / OSE, run_met, run_simulator -- an error budget
  rolled through the Jacobians to WFE / Strehl, then drift and noise
  through the loop.
- **The DM surface gauge bench** (the deck's future-work slide, plan
  `BRIEF_gauge_deck.md` section 10.3): the camera's throughput and
  electronics, optical surface errors on the bench optics, alignment and
  stability, more DM drift terms, and the JPL-form budget per reading.

Both are the same machinery (a validated linear model, a named budget,
covariance roll-up, a loop with noise and drift); the camera model and
the surface-error generator are shared.  CCL's assumption until Dave
says otherwise: you plan the TELESCOPE tools (your lane; the gauge
bench's items are in 10.3 and belong to the deck's owners), and you
build the camera and surface-error pieces so that the gauge lane can
call them (the gauge's photon accounting is the seed for the camera
model, not a second copy).

## Your question: static first, then feed the simulator

Budget first.  The static roll-up is the deliverable JPL reads
(allocations in, RSS WFE / Strehl out, margins by contributor), it is
fast, and its contributor table IS the right-hand side the dynamic
predictor draws from -- so building it first means the simulator arc
starts with its inputs defined and validated rather than invented per
run.  Then the dynamic predictor consumes that same table with time
behavior attached (a time constant and a class per contributor).

## Steer for the plan (CCL's rulings, take them as given)

1. **Covariance is the primary computation; Monte Carlo is the check.**
   The observation model is linear, so a contributor with covariance P
   gives WFE^2 = trace(P G) with G = D'D the Gram of its Jacobian --
   exact, fast, and the form the MET optimizer already uses (never the
   outer product: the tall-Jacobian OOM rule in the memory).  Run MC
   only where linearity or the loop breaks it (contrast, saturation,
   nonlinear estimators), and gate every covariance number against an
   MC of the same contributor alone.
2. **The budget lives as a named allocation table**, one row per
   contributor: what (DOF / surface / subsystem), magnitude, spatial
   class (rigid body, Zernike low order, PSD mid-frequency grid),
   temporal class (fixed after calibration, slow drift with a time
   constant, white per measurement), correlation group.  Rolled up per
   subsystem and total, with the sensitivity-weighted margin per row.
   One file format, read by both the static and the dynamic tool.
3. **Surface errors:** a PSD-to-realization generator (low order as
   Zernike coefficients through dwdz; mid-spatial as a GridData
   realization through dwdgrid or placed on the element), seeded and
   saved so every realization is reproducible; per-class allocations
   (M1 / M2 / relay / DM).  The gauge deck's 10.3 item 2 uses the same
   generator on the bench optics.
4. **Drift:** first-order Gauss-Markov (a time constant and an rms) per
   contributor class on top of run_simulator's random walk; thermal =
   slow, low-order, correlated across DOFs of one body.  With the loop
   linear, the steady-state closed-loop covariance is a discrete
   Lyapunov solve -- the dynamic predictor's fast form; the time-domain
   simulator is its check and its plot.
5. **Electronics:** one camera model, shared with the gauge lane:
   photons per measurement, well depth and frames co-added, read noise
   per frame, quantization, gain nonlinearity, latency; injected at the
   sensing stage so the estimator sees it.  Photon noise is already
   there in the gauge work -- extend it, do not duplicate it.
6. **Gates:** zero-error budget reproduces the nominal to round-off; each
   contributor's covariance roll-up agrees with its own MC (state the
   tolerance from the MC's sample count); the dynamic predictor's
   steady state agrees with run_simulator's long-run rms.  Every number
   in a committed report; a parameterized runner from the start.
7. **Contrast is a follow-on**, not in the first roll-up: it is not
   linear in the WFE; it goes through the CTB diffraction layer once the
   WFE budget stands.

Plan it in plan mode as you did the Luis fixes; CCL reviews before you
build.  Order: budget table + static roll-up + surface-error generator +
its gates; then the camera model; then drift and the closed-loop
covariance; then contrast.
