### `bench_layout`: Laying out and optimizing an optical bench

A generic folded-bench worked example: build a complete test-bench train SEQUENTIALLY with the macos.design.Bench add-optic utilities, then optimize it in three physically-staged steps. The bench is the classic surface-gauge / interferometer test-arm topology: point source -> baffle (pupil stop) -> L1 (collimating singlet) -> beam-splitter REFLECT -> DM (test mirror, retro) -> beam-splitter TRANSMIT (through the tilted plate, real walk-off) -> Fold1 -> Fold2 -> L2 (imaging singlet) -> focal-mask Reference (internal focus) -> Detector at the DM-PUPIL IMAGE.

**Run it.**

```matlab
cd('<path-to>/mmacos/templates/40_benches/bench_layout');
run(fullfile(pwd,'example_bench_layout.m'));
```

**What it prints:** Section 2 ran it unchanged: the builder's chief ray against the engine's, the staged solve, and the optimized bench.

**Its figure:** the inputs are FIGURE 6 and the optimized bench FIGURE 7, in Section 2.

**The number it reproduces.** The optimized bench is %pin:ch02_bench_nelt% elements, and the builder's chief ray and the engine's agree to %pin:ch02_bench_chief_agree% (Section 2).

**Adapt it.** *"Using templates/40_benches/bench_layout, rebuild the test-arm bench for a 48 mm DM with the same fold angles, and tell me the optimized element count and whether the builder's and the engine's chief rays still agree."*

