## APPENDIX E: Zernike Orderings

A Zernike surface (Section 4, “Zernike Surfaces”) adds to its conicoid the excess sag Σ c~k~ Z~k~(R, A), where the coefficients c~k~ are the `ZernCoef=` values in BaseUnits, R is the radial distance in the xMon–yMon plane divided by lMon, and A is the angle from xMon toward yMon. The sag is added along zMon. `ZernType=` sets which Zernike term Z~k~ each index k names, and whether the term is normalized. The same names apply to `MonZernType=` and `FFZernType=` on FreeForm surfaces.

Every term is a radial polynomial in R times cos(mA) for m > 0 or sin(|m|A) for m < 0, and is identified by its radial order n and azimuthal order m. Table E.1 gives the index of each term in every ordering, and the polynomial itself.

-   **ANSI** is the default whenever a Zernike-carrying surface declares no `ZernType=`. It is the ANSI Z80.28 (OSA) ordering, numbered from 1: the term (n, m) has index (n(n + 2) + m)/2 + 1, so the index runs through each radial order with m increasing from −n to n. `Malacara` is accepted as its legacy name (the ordering the 3.2 manual credited to Malacara's *Optical Shop Testing*).
-   **BornWolf** runs through each radial order with m decreasing from n to −n. The column that the 3.2 manual's Table 8 labeled “Code V (Noll)” is this ordering, not Noll's.
-   **Noll** is Noll's 1976 ordering. Within a radial order, |m| increases, and of each cos/sin pair the even index is the cosine.
-   **Fringe** is the 37-term Fringe (University of Arizona) set, ordered by (n + |m|)/2 and then by decreasing |m|. As in CODE V, its 37th term is the order-12 spherical term, not the next term of the ordering. **ExtFringe** (Extended Fringe) continues the same ordering past term 36.
-   **NormANSI**, **NormBornWolf**, **NormNoll** and **NormFringe** use the same index map as their unnormalized types, with every term multiplied by its RMS normalization factor √(2(n + 1)) for m ≠ 0, or √(n + 1) for m = 0. Each normalized term has an RMS of 1 over the unit circle, so a coefficient is that term's RMS contribution.
-   **NormAnnularNoll** takes the annulus ratio as a trailing value (`ZernType= NormAnnularNoll 0.3`). **NormHex** names hexagon-normalized terms. What each traces as today is in the notes under Table E.1.

ANSI, BornWolf, Noll, ExtFringe and their normalized forms accept up to 66 coefficients. For every type but ExtFringe, 66 is every term through radial order 10 (for ExtFringe, see the notes). Fringe and NormFringe accept 37.

The table is measured, not transcribed. `tools/gen_zernike_table.py` (`make zernike-table`) loads a flat Zernike mirror through *pymacos* for each type and sets one coefficient to a unit value at a time. It then traces the ray grid and fits the sag where the rays meet the mirror against every Zernike term through radial order 14. A mode prints as a row when its sag is exactly one term with coefficient 1 (or the RMS factor, for a normalized type). The gate `pymacos/tests/test_zernike_orderings.py` evaluates each printed polynomial against the engine's sag for every column. It also checks the indices against the textbook rules for ANSI, Noll and Fringe, and against the 3.2 table, and confirms that a row with sin and cos exchanged fails.

<!-- transcript: zernike_orderings -->
