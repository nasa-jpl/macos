### `coro_planet_demo`: Illustrated coronagraph + off-axis planet injection

`coro_planet_demo()` runs the canonical CoroExample recipe through the mmacos framework and writes two 4-panel figures (matching the CoroExample.jou "pgp 4" plot pages) into the example folder: coro_surfaces.png the coronagraph surfaces (jou: int 1/6/11/16) [entrance pupil | star occulted at FPM | pupil at Lyot | star at det] coro_planet.png the planet injection (jou: ffp + add) [suppressed star | planet alone | star+planet | difference]

**Run it.**

```matlab
cd('<path-to>/mmacos/templates/30_instruments/coro_planet_demo');
out = coro_planet_demo();
```

**What it prints:** Part I does not run it: Section 2's coronagraph is the CTB behind the thread's telescope. Run with its defaults (model 512, CoroExample.in), it writes its surfaces page and its planet page and prints only the file names.

**Its figure:** none printed. Its surfaces page is sound (the Lyot pupil's bright rim is the mask working); its planet page is finding (ddd), below.

**The number it reproduces.** The template prints no number: its product is two pictures. The surfaces page is sound (the Lyot pupil's bright rim is the mask working). The planet page is not: asked for a planet at sky-angle cosines [+2e-5 +2e-5], it composes a shifted, clipped speckle field rather than a planet image, finding (ddd), and Part I prints only the surfaces page.

**Adapt it.** *"Using templates/30_instruments/coro_planet_demo, find why the planet page shows no planet image by printing where the engine puts the planet's chief ray at the detector for the requested sky angle, before changing anything."*

