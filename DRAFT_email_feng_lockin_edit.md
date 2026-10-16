To: Feng (cc Fang)
Subject: Re: COPHI analysis — the lock-in camera

Hi Feng!

That camera sounds like fun -- did it work OK?

I'm attaching updated slides, with more detail on the lock-in detector, showing more clearly how the demodulation reference is used in the pixels to extract the I and Q signals, which are ADCed on-chip, so we can (as you say) use their firmware. Do they use the 4-frame approach? Their pixels do a two-tap demodulation at 0/180° and 90/270°, the four-step done in the charge domain, so far as the datasheet shows.

The performance difference comes from integrating many cycles in the pixel instead of one quarter-cycle per frame. The larger effect is photon throughput: the lock-in's well empties every 1 ms sample, so a pixel can use ~5 × 10⁸ e⁻/s, where a 200 ke, 150 fps frame camera is capped at 3 × 10⁷ — about 17× the photons per second, which accounts for most of the 42 → 3.4 s per picometer lock-in improvement. The second effect is read noise: a frame camera pays r² at every bucket, four per beat, while the lock-in pays r² once per sample of N cycles — 1/N in variance per photon — which is why an assumed 100 e⁻ moves the lock-in only from 3.4 to 4.0 s per picometer. Photon noise per photon is the same in both.

Confirmation of these performance calcs depends on the answers to a couple questions we should ask the vendor:

1. Do the storage nodes have to hold the full DC of a sample before A − B is formed, or is the common part drained continuously so the nodes hold only the difference? In the first case the well caps the photons per sample; in the second it caps only the fringe term and the camera is ~10× faster at the same flux.
2. Is the 1 MSa/s internal sampling an analog accumulation read once per 5.6 kSa/s output sample, or a read per internal sample? Only the first keeps the read-noise dilution.
3. What is the phase noise per pixel versus flux at ~100 ms integration?
4. How is the external reference used -- direct clock or an internal PLL -- and what is the stability of the phase delay?

The revised BDs hopefully show how the 5 kHz beat locked to the AOM synthesizers goes to the camera's reference input, driving the pixels' demodulation, and our software does only atan2(Q, I), the registration and the actuator fit. Computing I/Q in software would need four frames per 5 kHz beat, 1/quadrant — so frame rate of 20 kfps at a 5 kHz beat — to match read noise.

There is also more detail on the resolution required to get best performance -- see slide 19. Bottom line, 3.3 pixels/actuator spacing gives the best performance, curves are flat beyond that.

Dave


