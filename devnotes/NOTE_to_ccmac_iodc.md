# NOTE to CCMac — the IODC lens design problems: pull, review, select (CC for Dave, 2026-10-23)

Dave wants a design challenge from the literature for the mmacos manual (`PLAN_MMACOS_MANUAL.md`, "Challenges
from the literature"), scored against published entries.  The IODC has posed one lens design problem per
conference since 1990; the statements and results are in the SPIE proceedings, which JPL's library reaches and we
do not from here (the open pages 403 the fetch).  Your task: pull the papers below, extract the facts into one
record, and recommend one problem to enter first and one to hold as the stress test.  Record:
`cli_tests/records/iodc_problems_2026-10.md` (plain facts; quote constraints verbatim; no prescriptions that the
proceedings restrict -- say where they are instead).  No engine or template work in this item.

## The candidates (what the open abstracts say; verify each line against the paper)
| Year | Paper | Problem, as far as the abstracts go |
|---|---|---|
| 2010 Jackson Hole | Juergens, "The 2010 IODC lens design problem: the green lens", Proc. SPIE 7652 | 100 mm EFL; EVERY surface shares one radius of curvature; maximize semi-field x entrance pupil diameter; 37 entries, 8 countries |
| 2014 Kohala | Juergens, "The 2014 IODC lens design problem: the Cinderella lens", Proc. SPIE 9293 | 100 mm EFL built from TEN N-BK7 glass blanks (fixed glass); maximize EPD x semi-field with RMS WFE <= 0.070 wave over the field; 45 entries, 13 countries; element counts 10-52; a winning entry quoted at EPD 81.3 mm, semi-field 43.5 deg (confirm it is this problem's) |
| 2017 Denver | Juergens, "The 2017 IODC lens design problem: the centennial lens", Proc. SPIE 10590 | the "100 lens": maximize the used diameters of ball lenses, diffraction-limited; winner 54 elements, 19.4 m long, EFL 0.02 mm -- a puzzle, not a candidate |
| 2021 virtual | Juergens, "The 2021 IODC lens design problem: the down under lens", Proc. SPIE 12078, 10.1117/12.2605550 | 100 mm EFL, monochromatic, must perform identically flipped end to end WITHOUT reversing the radii; maximize EPD x semi-field at RMS WFE <= 0.070 wave; 43 entries, 13 countries (the abstract); winner quoted 198.4 mm pupil, 90 deg semi-field, product 17,857 of 18,000 possible |
| 2023 Quebec | Juergens, "The 2023 IODC lens design problem: the Yin-Yang lens", Proc. SPIE 12798 paper 1279837, 10.1117/12.3012062 | all refractive, rotationally symmetric, spherical or plano only, no cemented surfaces, no intermediate images; diffraction-limited over the field at 550 nm AND bad on-axis image quality at 540 nm at the same image location; maximize EPD x semi-field; 24 entries, 7 countries; winner Arnaud Davenel (Safran) |
Older: 1998 Kona (a solid glass lens, no air spaces, Proc. SPIE 3482); 2002 Tucson (100 mm f/8, spherical
refractive only, mimicking a diffractive element's chromatic behaviour).  Akira Yabe's page
(www.yabeopt.com/iodc.htm) carries his own entries and possibly prescriptions -- fetch it from a browser; its
certificate is wrong for the host name.  Note that the web summaries of 2021 and 2023 contradict each other on
entry counts and winners (Yabe's page vs the SPIE abstracts); take the proceedings as the record.

## Extract, per problem
1. The statement verbatim: EFL, wavelength(s), glass rule, surface types allowed, symmetry, field and pupil
   definitions, vignetting rule, the wavefront limit and how it is evaluated (which fields, which reference), the
   merit and its maximum, packaging constraints (the 2014 blank dimensions!), the submission format (a CODE V /
   Zemax file? -- that decides whether `rx_converter` gets a test).
2. The results: entry count, countries, the winner's numbers, and whether the winning prescriptions are printed
   in the paper (several years print the top entries).  If printed, say so and give the page; do not copy them.
3. The software the entrants used, if stated (it says what MACOS is being compared with).

## The criteria for the recommendation (Dave's two aims)
- **To SUCCEED**: a problem whose solution space our tools cover -- fixed glass (CALIB varies curvatures,
  spacings, conics and aspheres, NOT glass choice), monochromatic, spherical, a moderate element count.  By the
  abstracts that is 2014 (ten N-BK7 blanks: the glass is given) and perhaps 2010.
- **To TEST MORE of macos**: free glass (the glass-table path the consolidation repaired), two wavelengths (2023),
  wide fields with vignetting and TIR at 40-90 degrees semi-field (2014's 43 deg, 2021's 90 deg), many surfaces
  (40+ DOF for CALIB's LM), the reversibility constraint (2021), and the converter for comparing against the
  published winners.
CC's prior: enter 2014 first (fixed glass, monochromatic, spherical, 45 published entries to score against) and
hold 2023 as the stress test (two wavelengths, free glass, a merit CALIB cannot express directly -- it needs an
outer loop over EPD and field around the least-squares WFE solve, which is a MATLAB driver, not an engine change).
Say whether the papers change that.

## Then
Dave picks; the chosen problem becomes `mmacos/challenges/iodc<year>/` with the README's scoring-convention rule
(state the metric before the number), a runner that scans EPD x field around CALIB, and the published winner as
the target.  Not before Dave's word.
