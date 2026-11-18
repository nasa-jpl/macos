**Dave's sent version, 2026-10-09** (CC's draft: `DRAFT_email_jim_round3.md`; the deck sent: `demo_session/deck_dyson_3k.pptx`, FINAL).

Hi Jim and Joe!

Thanks for the answers and the paper -- they settled the trade toward two 3k modules and pointed us at the telescope form, and the "3k module is not there yet" of the last note is no longer true. Both modules now meet Joe's specification end to end, with SRF at the two-pixel slit's own floor. Details in the attached deck; the summary follows.

We built the telescope as a design template rather than a one-off. It follows the SBG zig-zag three-mirror form from the paper's Fig. 4b, scaled to Joe's 330 mm, with the first order solved in closed form (stop at M2 makes the slit telecentric exactly) and the figure solved against rows for telecentricity, a bounded cone (no marginal ray faster than the spectrometer accepts, Jim's "no F/1.2 rays"), plate scale, and the chief-minus-centroid offset at the slit. The mirrors are off-axis conic sections with a pole-frame polynomial of degree 3–6; even aspheres alone did not reach pixel resolution on either module.

The 3k module, scored (engine ray traces on 18 µm pixels; no coatings, no tolerances; the paper's bounds converted to micrometers):

| criterion | Joe | paper | spectrometer alone (240 mm CaF2 Dyson) | telescope alone | together, end to end |
|---|---|---|---|---|---|
| smile | < 0.1 px | < 1.8 µm | 0.005 px | – | 0.022 px slit-filled (0.4 µm); 0.016 px point-source |
| keystone | < 0.1 px | < 1.8 µm | 0.006 px | – | 0.006 px (0.1 µm) |
| CRF FWHM | < 1.5 px | < 50.4 µm | 1.21 px | 1.02 px along the slit | 1.17 px (21 µm) |
| SRF FWHM | 1.5–2.0 px | < 64.8 µm | 2.024 px | 1.02 px across (ARF 18 µm) | 2.025 px (36.5 µm) |
| energy in one pixel | > 0.75 | – | 0.82 | 1.00 | 0.85 |
| marginal rays at the slit | – | F/# anamorphicity bounded | – | F/1.80–1.89, none below 1.7 | grating admits 98.9 % |
| clearance, worst pair | > 0 | – | +0.54 mm | +12 mm | +0.6 mm |

The 1.5k module was designed with same template: 1.02 px at every field of its strip in five rungs; end to end smile 0.009, keystone 0.009, CRF 1.03, SRF 2.024, 0.97 of the energy in a pixel, 98.7 % through the grating, clearance +0.6 mm. One deviation on the way: the first rung at the pixel put the M3-to-slit beam 0.3 mm inside M2's body, and a clearance wall row in the solve bought +2.4 mm with the image held. It improves on the three-asphere 1.5k of the last note in every line but one (keystone at one roll, 0.012 vs 0.009, both a tenth of the bound). The price is size: this form's M2 and M3 are full aperture, so the four-module block is 494 L / 21.4 kg of optics against the previous version's three-asphere at 300 L / 14.6 kg.

With a 2-pixel slit a perfect spectrometer gives SRF = 2.010–2.023 px over the band (slit, pixel, diffraction); ours sits 0.001 px above that, alone and behind the telescope -- 36.5 µm against your 64.8. We read Joe's "1.5–2.0" as a slit-width choice, not a miss.

So for two telescopes, spectrometers and detectors (3072 × 512, a standard format) against four with a 1.5k × 0.5k detector that is NRE; 609 L and 39 kg of optics for the pair against 494 L and 21 kg for the four; integration, test and calibration twice instead of four times -- your answer, and the paper's.

Our Dyson of record is the paper's Option A (aspheric CaF2 lens, spherical grating). Every large sensitivity is a focus term the detector's one-axis focus removes (grating radius the strongest: dR/R 1e-4 costs CRF +1.9 px uncompensated, 0.006 after a 76 µm refocus). Keystone is lateral color and sets the alignment tolerances with no compensator: for the paper's as-built 1.03 µm, lens and grating decenter along the slit 26 µm, grating tilt 35 µrad, clocking 88 µrad. Thermal: the aluminum air space (−0.22 px/K) and the grating radius (+0.14 px/K) are the open focus terms, which is your "CRF limits the temperature requirements." Option B (the conic grating) takes energy in a pixel from 0.82 to 0.99 and CRF from 1.21 to 1.03 but raises the keystone sensitivities 6–45 % -- it buys nominal margin, not alignment insensitivity, in our alignment-only ladder.

Throughput by band, including the air gap at the order-sorting filter: uncoated 0.82 in every band; a single quarter-wave MgF2 centered at 2.2 µm lifts the long band to 0.87 and costs nothing elsewhere; a two-layer tuned to 1 µm is better mid-band (0.92) and worse at both ends; a 1.8 µm blaze puts 0.89 of first-order efficiency in the long band (scalar estimate).

Questions:

1. The entrance beam. Joe's 183 mm at f 330 mm is F/1.80, so the telescope's cone along the slit is F/1.89 and cannot be "a little faster than the spectrometer" without a larger beam -- the paper's 192 mm. Is the aperture Joe's number or the paper's?
2. The slit and SRF. Is the 1.5–2.0 px range read against a 2-pixel slit (where this design sits at the floor), or does the project want a narrower slit at the throughput's expense?
3. The spectrometer's form. What drove Option B in the build -- fabrication, thermal, or the margin? Our alignment ladder does not reproduce a tolerance advantage for it.

The template (parameters, runner, records, test) is in the public MACOS_resources tree under templates/10_telescopes/tma_longslit; any long-slit front end can run through it. And as before, this round surfaced macos and mmacos items now fixed with tests, which is the main purpose here--

        -Dave (and CC)
