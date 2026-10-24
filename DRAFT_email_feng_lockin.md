To: Feng (cc Fang)
Subject: Re: COPHI analysis — the lock-in camera

Feng — agreed on both counts, and your suspicion is close to how we modeled it.

What the lock-in pixel does.  The heliCam line (Heliotis; the CSEM demodulation-pixel
lineage you bought in 2005) steers each pixel's photocurrent between two storage sites by
the demodulation reference — a square-wave, two-tap integration synchronous with the beat,
over thousands of cycles per output sample, with the DC rejected in the pixel before the
well.  That is the four-step done inside the pixel: I = (tap A − tap B) at one reference
phase, Q at the quadrature phase (the C4's "dual-phase" demodulator does both).  So the
firmware's I/Q IS the bucket algorithm, with two things different from a frame camera: the
buckets integrate over many cycles instead of one quarter-cycle each, and nothing but the
difference is stored, so the well bounds the fringe, not the DC.

What that means for our numbers.  Our lock-in rows carry exactly this: the two-tap
square-wave demodulation costs ×1.11 in phase noise (×1.23 in time) against an ideal
sinusoidal multiplier, and that factor is on every lock-in row of the analysis; the
per-pixel step error and trigger jitter terms of the frame-camera four-step apply inside the
pixel the same way (we gate them in the model).  The speed advantage we report (3.4 s per
picometer against 42 for the best frame camera, 3 pm held at 0.125 s cycles) comes from the
cycle count per output sample and the DC rejection, not from a better algorithm.

On not duplicating firmware and software: we would run the camera's own demodulator with
the beat on its reference input (the C4 has one) and take I/Q out at the frame rate; the
software then does only atan2(Q, I), the registration and the actuator fit.  Computing I/Q
in software from raw frames would need the frames at the beat rate, which no camera gives
at kHz — so the on-chip path is the only one that uses the device for what it is.

The questions I would put to Heliotis before we go further, which are the ones our model
cannot answer: (1) the tap structure — two taps with alternating reference phase, or four;
(2) whether the full well bounds each internal sample (their 1 MSa/s) or the integrated
output (5.6 kSa/s) — this decides whether the camera is photon-limited at the well or not,
a factor of ten in our numbers; (3) the per-pixel phase noise versus flux at a 100 ms
integration; (4) how the external reference is used — direct clock or an internal PLL, and
the phase delay's stability.  Happy to talk it through among ourselves first.

Dave
