To: Feng (cc Fang)
Subject: Re: COPHI analysis — the lock-in camera

Hi Feng!

That camera sounds like fun -- did it work OK?

I'm attaching updated slides, with more detail on the lock-in detector, showing more clearly how the demodulation reference is used in the pixels to extract the I and Q signals, which are ADCed on-chip, so we can (as you say) use their firmware. They do use the 4-frame approach.

The performance differences come from integrating over many cycles instead of one quarter-cycle each, reducing read noise. Photon noise is the same in either case: every photoelectron contributes its shot noise whether the sum is formed in charge or in numbers, so σ_φ = π/(2γ√N_e) per pixel per sample holds for both (the π/2 is the square-wave multiplier's ×1.11 against a sine lock-in). What the charge-domain sum changes is the read noise: a frame camera pays r² at every bucket, four per beat, while the lock-in pays r² once per sample of N cycles -- 1/N in variance per photon. That is why an assumed 100 e⁻ of read noise moves the lock-in only from 3.4 to 4.0 s per picometer, where the Photonfocus's 220 e⁻, paid every frame, costs it 22 % at half well.

Two things set the ceiling on that gain, and they are the questions I would put to Heliotis -- their datasheet cannot answer them from outside:

1. Do the storage nodes have to hold the full DC of a sample before A − B is formed, or is the common part drained continuously so the nodes hold only the difference? In the first case the well caps the photons per sample (our "well-capped" numbers, ~3.4 s per pm); in the second it caps only the fringe term and the camera is ~10× faster at the same flux. Draining the DC removes its charge, not its shot noise, so the gain is more photons per sample, not less noise per photon.
2. Is the 1 MSa/s internal sampling an analog accumulation read once per 5.6 kSa/s output sample, or a read per internal sample? Only the first keeps the read-noise dilution.
3. Phase noise per pixel versus flux at ~100 ms integration, measured.
4. How the external reference is used -- direct clock or an internal PLL -- and the stability of the phase delay.

On your point about not duplicating the demodulation: agreed. The beat (5 kHz, locked to the AOM synthesizers through the Rb clock) goes to the camera's reference input, the pixels demodulate, and our software does only atan2(Q, I), the registration and the actuator fit. Computing I/Q in software would need frames at the beat rate, which no camera gives at kHz.

Dave


