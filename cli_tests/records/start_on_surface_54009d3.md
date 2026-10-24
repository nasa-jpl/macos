# A ray that STARTS on a surface: is it called a miss? (PLAN_CONSOLIDATION 1.2, measured 2026-10-12)

**Measurement only -- no surfsub change (CC's ruling: none without Dave).**  Engine at macos 54009d3, an
INSTRUMENTED side build (not committed), ifx and gfortran, model 256, each deck: LOAD + `opd nElt`.

## Method
- Every FORWARD-branch `L < 0` rejection (`GO TO 98`) -- 24 sites: FlatSrf, ConSrf, AsphSrf, AsphGridSrf,
  FreeFormSrf, UDSrf (surfsub.F); ReturnSrf, FocalPln, RefSrf, PolElt, ObsSrf (elemsub.F); IntSrf (didesub.F) --
  records the root nearest zero (`|L|`, or `min(|kp|,|km|)` after `FwdRoot`) and the element scale
  `max(1, |pin|_inf, |pv|_inf)`.
- At CTRACE's Miss stamp (sequential and NS sites) a miss is IN THE BAND when that root is within
  **16 ulp of the scale** (`DeadBandUlp`), and ON A COINCIDENT PLANE when, in addition, the element's vertex
  lies on its predecessor's plane (`|psiElt(prev).(VptElt(e)-VptElt(prev))| <= 16 ulp`).
- Per trace: misses, band misses, coincident band misses, first / last element of a band miss.

## Non-vacuity (the instrument sees the known case)
`Grating_example_001.in` BEFORE its fixture fix (resources 82c58bb^, ExitPupil vertex on element 2's plane):

| compiler | nPass | misses | in band | coincident | element |
|---|---:|---:|---:|---:|---:|
| ifx | 7303 | 5150 | 43 | 43 | 3 |
| gfortran | 7295 | 5158 | 51 | 51 | 3 |

The band misses are EXACTLY the rays each compiler lost against the fixed deck (nPass 7346 on both:
7346 - 7303 = 43, 7346 - 7295 = 51) -- the round-off coin toss the 1.1a record named.  The fixed deck: 0.

## The corpus
| corpus | decks | compiler | decks with misses | decks with band misses |
|---|---:|---|---:|---:|
| committed (`cli_tests/corpus_committed.txt`, both repos) | 1279 | ifx | 20 | 0 |
| committed | 1279 | gfortran | 20 | 0 |
| the 406-deck A/B list (`decks_all.txt`) | 406 | both | 18 (ifx) | 1 -- the OLD Grating_example_001 copy in the stale `~/dev/MACOS_res_dev` worktree (43 / 51, the table above) |

Decks that lose rays to misses (committed corpus; misses ifx / gfortran; none in the band):

| deck | nElt | misses ifx | misses gfortran |
|---|---:|---:|---:|
| ~/dev/MACOS_resources/mmacos/challenges/afocal4/descent/afocal4_descent_N7c.in | 8 | 442 | 442 |
| ~/dev/MACOS_resources/mmacos/challenges/afocal4/offaxis/afocal4_OAM5_N5_h0.in | 6 | 448 | 448 |
| ~/dev/MACOS_resources/mmacos/challenges/afocal4/offaxis/decks/oa_greg_N5.in | 6 | 522 | 522 |
| ~/dev/MACOS_resources/mmacos/challenges/afocal4/offaxis/decks/om_greg_f1.25_h1.5.in | 3 | 923 | 923 |
| ~/dev/MACOS_resources/mmacos/challenges/afocal4/offaxis/decks/om_greg_f1.25_h1.in | 3 | 198 | 198 |
| ~/dev/MACOS_resources/mmacos/templates/10_telescopes/offset_imager/t5_redemption/t5r_s3.in | 5 | 535 | 535 |
| ~/dev/MACOS_resources/mmacos/templates/10_telescopes/offset_imager/t5_redemption/t5r_s4.in | 5 | 518 | 518 |
| ~/dev/MACOS_resources/mmacos/templates/10_telescopes/offset_imager/t5_redemption/t5r_s5.in | 5 | 504 | 504 |
| ~/dev/MACOS_resources/mmacos/templates/40_benches/tg_psi_dm96_oap/runs/s2s_glassmove/s2s_glassmove_chk_flat.in | 46 | 51073 | 51073 |
| ~/dev/MACOS_resources/mmacos/templates/40_benches/tg_psi_dm96_oap/runs/s2s_glassmove/s2s_glassmove_chk_sin.in | 46 | 51073 | 51073 |
| ~/dev/MACOS_resources/mmacos/templates/40_benches/tg_psi_dm96_oap/runs/s2s_glassmove/s2s_glassmove_deck.in | 46 | 51073 | 51073 |
| ~/dev/MACOS_resources/mmacos/tests/Rx/Rx_Cass_NS.in | 6 | 1 | 1 |
| ~/dev/MACOS_resources/mmacos/tests/Rx/Rx_GratingImmersed.in | 4 | 61 | 61 |
| ~/dev/MACOS_resources/mmacos/tests/Rx/Rx_Luneberg.in | 17 | 8 | 8 |
| ~/dev/MACOS_resources/pymacos/tests/Rx/Grating_example_001.in | 4 | 5107 | 5107 |
| ~/dev/MACOS_resources/pymacos/tests/Rx/Rx_GratingConcave_air.in | 4 | 61 | 61 |
| ~/dev/macos/ZGD_test_files/tst_zrngr_overflow.in | 6 | 1 | 1 |
| ~/dev/macos/ZGD_test_files/tst_zrngr_roundtrip.in | 6 | 1 | 1 |
| ~/dev/macos/docs/macos-manual/examples/Luneberg.in | 17 | 8 | 8 |
| ~/dev/macos/docs/macos-manual/examples/emitted/Luneberg.in | 17 | 8 | 8 |

## Reading, for Dave
**Zero.**  Outside the one fixture already corrected (82c58bb), no deck in the committed corpus loses a single
ray to a start-on-surface root, on either compiler -- including every committed deck in the Return/Return/FocalPlane idiom (119 such decks in the
514-deck census `macos_f90/CLAUDE.md` cites; not re-counted here).
That is the case for LEAVING the forward branches as they are (a ray at `L` within round-off of 0 on a
coincident element stays a miss), with the fixture rule as the guidance: an element whose vertex sits on its
predecessor's plane is a deck error (the standard FEX sphere centred on the focus is the fix, as 82c58bb did).
The alternative -- accept `|L| <= 16 ulp` as a hit at `L = 0` -- would change no committed deck and would
remove the compiler sensitivity for a future deck that repeats the Grating fixture's geometry; it is a
surface-routine change at 24 sites, so it is Dave's call.  TO's recommendation: leave it; if wanted, a WARN
counter of band misses (rule 5: warn and count) costs nothing and names the deck that repeats the pattern.
