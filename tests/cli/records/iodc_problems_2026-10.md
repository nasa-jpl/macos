# IODC lens design problems — facts for the mmacos-manual challenge (CCMac, 2026-10-23)

Task: `NOTE_to_ccmac_iodc.md` — pull the IODC design-problem papers, extract the facts, recommend one to
enter first and one to hold as the stress test. For `PLAN_MMACOS_MANUAL.md` "Challenges from the literature".

## Sourcing — what was reachable from here, and what was not
**The SPIE proceedings full text is NOT reachable from this machine**, as the NOTE anticipated:
- `doi.org` → crossref → `chooser.crossref.org` landing loops (abstract metadata only, no body).
- `spiedigitallibrary.org` returns an **Incapsula bot-challenge** page (HTTP 200, ~1 KB, no content) to curl.
- `opg.optica.org` returns a bot-protection script; `researchgate.net` 403s. All aggregators block automation.
- The WebFetch tool's summarizer backend is unavailable in this gov environment (it errors on an unroutable
  model), and WebSearch is unsupported on this model — so no tool-side extraction either.
- **Akira Yabe's page IS reachable** — but only over **plain HTTP** (`http://www.yabeopt.com/iodc.htm`); the
  HTTPS host presents a wrong Sakura-hosting certificate (NOTE's warning) and 403s. Its text is quoted below.

So the statements below are **Yabe's wording** (a top-ranked entrant every year, authoritative-adjacent) plus
the NOTE's open-abstract facts. **Every line a design would be scored against — the exact merit maximum, the
WFE evaluation (fields/reference), the vignetting rule, the 2014 blank dimensions, the submission format —
must be confirmed against the SPIE proceedings via the JPL library before a challenge is built.** Points that
Yabe and the NOTE's abstracts already disagree on are flagged in place; take the proceedings as the record.

## Per problem (Yabe verbatim; NOTE abstract facts; gaps flagged)

### 2014 Kona — the Cinderella lens — Proc. SPIE 9293 (ITH... ; Juergens)
- **Statement (Yabe, verbatim):** "design a 100 mm focal length lens in which all the components of the lens
  can be manufactured from **ten Schott N-BK7 lens blanks 100 mm in diameter x 30 mm thick**. The lens is used
  **monochromatically at 587.56 nm**." (Yabe ranked third among 45 participants.)
- **NOTE abstract adds:** maximize EPD × semi-field with **RMS WFE ≤ 0.070 wave over the field**; 45 entries,
  13 countries; element counts 10–52; a winning entry quoted EPD 81.3 mm, semi-field 43.5°.
- **Glass:** FIXED (N-BK7 only), so CALIB's curvature/spacing/conic/asphere DOF cover the solution space; no
  glass choice. Monochromatic. 587.56 nm = the Fraunhofer d-line (helium).
- **GAPS for the proceedings:** the exact blank-accounting rule (how "manufacturable from ten blanks" is
  checked — total glass volume? per-element thickness ≤ 30 mm and diameter ≤ 100 mm?), whether aspheres/conics
  are allowed or spherical-only, the vignetting rule, how the WFE is evaluated (which fields, what reference),
  the merit's theoretical maximum, and the submission format (CODE V/Zemax → whether `rx_converter` is tested).

### 2010 Jackson Hole — the green lens — Proc. SPIE 7652 (Juergens)
- **Statement (Yabe, verbatim):** "design a 100 mm focal length lens in which **every optical surface has the
  same radius of curvature**, positive or negative, or is plano. The lens is used **monochromatically at
  532 nm** and is made of **only Schott N-BK7 glass**." (Yabe ranked fourth among 37 participants.)
- **NOTE abstract adds:** maximize semi-field × entrance-pupil diameter; 37 entries, 8 countries.
- **Glass:** FIXED (N-BK7). Monochromatic (532 nm). The hard constraint is the **shared-radius rule** — a single
  |R| value (sign free, or plano) across all surfaces: an unusual coupled constraint, not a native CALIB DOF
  (it needs one shared variable driving every curvature — a MATLAB driver around CALIB, not an engine feature).
- **GAPS:** the merit maximum, WFE/vignetting rules, submission format.

### 2021 virtual — the down under lens — Proc. SPIE 12078, paper 1207820, DOI 10.1117/12.2605550
- **Statement (Yabe, verbatim):** "design a **100 mm focal length lens where if the lens is flipped
  end-to-end but the radii are not reversed in sign, the lens has to perform the same as in its original
  form**. The lens is used **monochromatically**. The goal ... is to maximize the product of the entrance
  pupil diameter and the semi-field of view **while holding the RMS wavefront error to ≤ 0.070 wave within
  the field of view**."
- **Results (Yabe, verbatim):** "**43 entries from 13 different countries. Three different commercial lens
  design programs were used, along with four in-house programs.** The winning entry from **Damien Gawron** has
  an **entrance pupil diameter of 198.4083 mm and a semi-field of view of 90° for a merit function product of
  17,857 out of a maximum possible 18,000.**"
- **DISCREPANCY with the NOTE:** the NOTE (from the abstract) named the winner's numbers but not the name;
  Yabe gives the winner as **Damien Gawron** (merit 17,857 / 18,000, EPD 198.4083 mm, 90° semi-field). The
  NOTE's "198.4 mm / 90° / 17,857" matches Yabe's; the name is new. Take the proceedings as the record.
- **Constraint class:** the flip-without-sign-reversal symmetry is a bilateral constraint on the prescription
  (surface k must mirror surface N−k with equal signed radius) — expressible as a parameter-tying driver, not
  an engine feature. Monochromatic, free element count. 90° semi-field means extreme-angle ray survival (TIR,
  vignetting) — a good engine stress but hard to converge.
- **GAPS:** the EFL tolerance, the WFE reference/fields, the wavelength value, the 18,000 maximum's definition.

### 2023 Quebec — the Yin-Yang lens — Proc. SPIE 12798, paper 1279837, DOI 10.1117/12.3012062
- **Statement (Yabe, verbatim):** "design a lens which is **diffraction-limited over the field of view at
  550 nm but at the same image plane location has large on-axis RMS wavefront error at 540 nm**. The entrance
  pupil diameter and the semi-field of view are not specified, but are free to be selected by the designer so
  as to **maximize the product of them and the 540 nm on-axis RMS wavefront error**. **Only two different
  glasses are allowed to be used (Schott N-BK7 and Schott N-SF6).**"
- **Results (Yabe, verbatim):** "**24 entries were received from seven different countries. The winner is
  Arnaud Davenel of Safran Electronics & Defense in France.**" (Yabe himself ranked fifth of 24.)
- **DISCREPANCY with the NOTE:** the NOTE's criteria call 2023 "free glass"; it is **NOT free** — exactly two
  glasses (N-BK7, N-SF6) are permitted. That weakens 2023's value as the "glass-table path" stress test (the
  table is exercised, but only two entries). Its real stretch is the **two-wavelength opposed merit** (good at
  550, deliberately bad at 540, same image location) and the merit being a PRODUCT of three quantities
  (EPD × semi-field × 540 nm on-axis RMS WFE) — which CALIB's least-squares WFE solve cannot express directly;
  it needs an outer loop (a MATLAB driver), as the NOTE says.
- **Surface rule (NOTE abstract):** all refractive, rotationally symmetric, **spherical or plano only, no
  cemented surfaces, no intermediate images.** Confirm against the paper.
- **GAPS:** the diffraction-limited criterion (Maréchal? Strehl?), the field definition, the merit maximum,
  the submission format.

### 2017 Denver — the Centennial "100 lens" — Proc. SPIE 10590 — NOT A CANDIDATE
- **Statement (Yabe, verbatim):** "design a lens which includes a '100 lens' to commemorate the OSA's 100th
  anniversary. Traditional lens specifications ... are not specified; the goal ... is to **maximize the used
  diameters of the ball lenses** in the 100 lens while maintaining diffraction-limited optical performance."
- **Winner (Yabe):** Takeshi Akiyama, **54 lens elements, 19,443 mm long, 165 mm diameter, EFL 0.02 mm**, merit
  398.04 / 400. A commemorative puzzle (ball lenses, meter-scale), not a design our tools should chase. Exclude.

### Older (context only — Yabe verbatim, not candidates)
- **2006 lens shuffler** (SPIE 6342): split a 100 mm N-BK7 cylinder into two lenses that reshuffle to reform the
  cylinder; equal EFL/EPD/FOV, both diffraction-limited. (Yabe ranked first of 32.)
- **2002 diffractive simulator** (SPIE 4832, pp. 473–485): a glass design with the same focal-length-vs-wavelength
  dependence as a diffraction grating, 400–750 nm — chromatic overcorrection with abnormal glasses. Yabe posts
  **CODE V/Zemax files** (`diff.seq`, `diff.len`, `diff.zmx`) for his improved retry. (Ranked second of 42.)
- **1998 solid glass lens** (SPIE 3482, pp. 2–8): aberration control at three wavelengths, glass weight ≤ 1 kg —
  secondary-spectrum glass selection. Yabe posts `solid.seq / solid.len / solid.zmx`. (Ranked fifth of 41.)
- **1990 Monochromatic quartet** (ILDC, SPIE 1354, pp. 548–554): the global-optimization benchmark (two known
  best solutions). **1985 reversible lens** (SPIE 554): aberration control at magnifications −1/2 and −2.

## Prescriptions / submission format
- **None of the candidate years' winning prescriptions were obtainable from here.** Yabe posts CODE V `.seq` /
  Zemax `.zmx` / `.len` files only for his **1998 and 2002** retries — not for 2010/2014/2021/2023. Whether the
  proceedings print the top entries (several years do) must be checked at the JPL library; **do not copy any
  printed prescription into the repo** — cite the page and reproduce independently.
- That Yabe's archive is CODE V/Zemax confirms the entry format is a lens-design-program file → **`rx_converter`
  is on the comparison path** for any year we score against a published winner.

## Recommendation (confirms the NOTE's prior, with two corrections)

**Enter first: 2014 Cinderella.** Fixed glass (N-BK7 only), monochromatic (587.56 nm), 100 mm EFL, a large
published field (45 entries, winner ~EPD 81 mm / 43.5° semi-field to score against). The solution space is
exactly CALIB's: curvatures, spacings, conics, aspheres — no glass choice, one wavelength. The only non-native
piece is the **"ten blanks" manufacturability accounting**, which is a check applied to the result (a MATLAB
scorer over element thicknesses/diameters), not a solve constraint — so it fits the "state the metric before
the number" README rule cleanly. **Blocker to resolve first:** the exact blank-accounting rule and the WFE
evaluation, from the proceedings.

**Hold as the stress test: 2023 Yin-Yang** — but note the correction: it is **two-glass, not free-glass.** Its
stretch is real but different from the NOTE's framing: (1) **two wavelengths with opposed goals** (diffraction-
limited at 550, large on-axis WFE at 540, same image location) — macos traces per-wavelength fine, but the
*merit* is a product of EPD × semi-field × 540 nm-WFE that CALIB's least-squares cannot express; it needs an
**outer MATLAB loop** over EPD and field around the WFE solve (not an engine change); (2) spherical/plano only,
no cements, no intermediate images — constraints a driver enforces; (3) the glass-table path is exercised but
only across two glasses. It remains the hardest to express and the best test of the driver layer.

**Alternates:**
- **2010 green lens** — the simplest to converge (fixed N-BK7, 532 nm, 100 mm EFL) BUT the **shared-radius**
  rule (every surface the same |R|, or plano) is an unusual coupled constraint needing a one-variable-drives-all
  driver; a cleaner *engine* test than 2014 but a narrower optical result. A good second success story.
- **2021 down under** — the flip-symmetry constraint + 90° semi-field is the widest-field / extreme-angle
  engine stress (TIR, vignetting, ray survival), with a clear published target (Gawron, 17,857/18,000). Harder
  to converge than 2014; a strong alternate stress test if a geometric-symmetry driver is wanted over 2023's
  two-wavelength merit.

**Net:** 2014 to succeed, 2023 to stress — the NOTE's prior holds, with (a) the 2023 "free glass" corrected to
"two glasses (N-BK7, N-SF6)" and (b) the 2021 winner named (Damien Gawron). Dave picks; the chosen problem
becomes `mmacos/challenges/iodc<year>/` only on his word, with its metric confirmed from the proceedings first.

## For Dave / the JPL-library pull (the lines a scorer needs, per year chosen)
EFL tolerance; wavelength(s) exact; glass rule; surface types allowed (asphere/conic vs spherical/plano);
symmetry/coupling rule stated verbatim; field & pupil definitions; vignetting rule; the WFE limit AND how it is
evaluated (which fields, which reference); the merit function AND its theoretical maximum; packaging constraints
(2014's blank dimensions above); submission format; whether the top prescriptions are printed (page numbers) and
the software entrants used. Sources consulted: `http://www.yabeopt.com/iodc.htm` (plain HTTP); the DOIs
10.1117/12.2605550 (2021) and 10.1117/12.3012062 (2023); SPIE vols 7652/9293/10590/12078/12798 — bodies at the
JPL library.
