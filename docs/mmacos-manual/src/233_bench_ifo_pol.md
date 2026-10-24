### `bench_ifo_pol`: Bench IFO polarization: ray-level Jones scoring (slice 1 of 3)

Polarization-honest Twyman-Green: a real Al coating on the beam-splitter face, ifPol on, scored with the GATED Phase-2 layer (jones_pupil + pol_maps). The coated BS face is hit as an EXTERNAL reflection in the test arm (air->Al, element 'BSrefl') and as an INTERNAL reflection in the reference arm (glass->Al, two glass transits, element 'BScrefr') -- genuinely different Jones -- so the arm-differential diattenuation / retardance at recombination is a MEASURED result, and its pupil variation is the polarization-induced PSI phase-error contribution.

*This template has no README; the text above is generated from its entry script's header (`docs/mmacos-manual/part2/bench_ifo_pol.md`, `tools/gen_readme_part2.py`).*

**Run it.**

```matlab
cd('<path-to>/mmacos/templates/90_polarization/bench_ifo_pol');
run(fullfile(pwd,'example_bench_ifo_pol.m'));
```

**What it prints:**

<!-- transcript: p2_bench_ifo_pol_run -->

**Its figure:** none: slice 1 prints its gates; the figures beside it belong to slices 2 and 3 (`example_bench_ifo_pol_slice2.m`, `_slice3.m`).

**The number it reproduces.** The coated beamsplitter's reflection matches the Fresnel equations to the floor of the arithmetic at 45 degrees (gate 1), and the coating is inert with polarization off (gate 2). Between the arms at recombination the plate leaves a diattenuation of %pin:p2_bench_ifo_pol_dia% and a retardance of %pin:p2_bench_ifo_pol_ret% rad, uniform across the pupil, so it costs the co-polarized fringe almost nothing.

**Adapt it.** *"Using templates/90_polarization/bench_ifo_pol, run slice 2 to trade the beamsplitter's angle of incidence against clearance, and tell me the angle that keeps the arm-differential retardance under 0.05 rad."*

