# Header-chain key matches: the prefix audit (PLAN_CONSOLIDATION 1.1b, 2026-10-12)

`LCMP(VAR_NAM,'Key',n)` compares the first n characters, case-insensitively, and the FIRST match in the
chain wins.  99 active header matches (msmacosio.inc, label-50 chain; commented-out and `#if 0` lines
excluded).  **Sibling captures** (an earlier prefix that a LATER key name starts with): one, `JZLouChfRayDir`
(9) over `JZLouChfRayPos` -- both disabled code, left alone.  **Look-alike exposure**: 23 active matches used
fewer characters than the key name, so any key sharing the prefix was read as that key (the brief's
`OptFooKey=` -> FOV weights).  All 23 now match the FULL name.  Deck usage below = decks in the 1353-deck
extended corpus (`cli_tests/corpus_extended.txt`, sandbox included) declaring the key; NO deck and no
mmacos/pymacos writer uses a shorter or variant spelling of any of them, so no deck moves.

| line | key | was n | now n | decks using it |
|---:|---|---:|---:|---:|
| 133 | `Save3DApVec` | 9 | 11 | 0 |
| 213 | `PrtEltRots` | 6 | 10 | 0 |
| 221 | `OptTarget` | 7 | 9 | 7 |
| 248 | `OptRayGrid` | 6 | 10 | 1 |
| 263 | `OptZernModes` | 8 | 12 | 0 |
| 392 | `OptNomSens` | 9 | 10 | 0 |
| 415 | `OptSvdSvCut` | 9 | 11 | 0 |
| 421 | `OptFOVWt` | 5 | 8 | 1 |
| 483 | `ZCOZernType` | 10 | 11 | 3 |
| 503 | `UseChfRay4OPD` | 10 | 13 | 2 |
| 548 | `BaseUnits` | 7 | 9 | 1341 |
| 570 | `WaveUnits` | 7 | 9 | 1343 |
| 672 | `pData` | 4 | 5 | 375 |
| 681 | `xData` | 4 | 5 | 375 |
| 688 | `yData` | 4 | 5 | 375 |
| 695 | `zData` | 4 | 5 | 375 |
| 702 | `lData` | 4 | 5 | 3 |
| 950 | `drawXGrid` | 6 | 9 | 0 |
| 957 | `drawYGrid` | 6 | 9 | 0 |
| 1018 | `nRings` | 5 | 6 | 3 |
| 1025 | `nPetals` | 5 | 7 | 3 |
| 1038 | `radCtr` | 4 | 6 | 3 |
| 1053 | `FlowerSegSrc` | 11 | 12 | 0 |

Residual (unchanged): a token that EXTENDS a full key name (`OptFOVWtX=`) still matches it -- `LCMP`
does not compare lengths; the element chain's matches were not audited here (not in this item).
