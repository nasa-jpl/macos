<!-- GENERATED FILE — do not edit by hand.
     Regenerate with: python tools/gen_appendix_c.py
     Source: macos_f90/macos_help.inc  -->

## APPENDIX C: Complete Command Reference

The following is the complete command list as printed by the HELP command, regenerated from macos_f90/macos_help.inc (single source of truth).  Casing convention: \<UPPER prefix\>\<lower tail\> — the uppercase part is the minimum-match abbreviation tested by the dispatcher.  Tags: \[Rx\] needs a loaded prescription; \[BLD\] needs a built linear model; \[DIFF\] needs a propagated wavefront.  For a full per-command catalog — dialogs, behavior notes, and the mmacos/pymacos programming interfaces — see the companion **MACOS Command Reference** (macosCmdRef.pdf; source in docs/macos-manual/cmdref/).

MACOS command summary

Casing: UPPER part is the minimum-match abbreviation.

Tags: [Rx]=needs OLD/NEW, [BLD]=needs BUILD,

```
       [DIFF]=needs PROpagate
```

#### SESSION & FILES

```
   Quit, EXIT, BYE - exit MACOS
   HELP        - this command list
   REset       - reset options & defaults to startup
   STatus      - what is loaded / built / propagated
   SUMmarize   - pretty-print the optical system
   ELTS        - list elements (type, surface)
   EXEcute     - run a .jou journal file
   JOUrnal     - generate a .jou journal file
   ECHo <text> - print <text> to stdout (annotate
                 console / journal output)
   MREset      - reset MACOS model size (128..8192)
   ! <cmd>     - run a shell command
   PWD, CD     - working directory shell helpers
   RX, LS, LL  - list .in files / current directory
   VI, EMAcs   - open file in editor
```

#### PRESCRIPTION I/O

```
   NEW         - start a new optical system
   OLD <file>  - load a .in prescription
   FID <id>    - load by file id (from RX listing)
   VALidate <file> - syntax-check a .in prescription
                     without loading it
   SAVe        - write current system to .in
   EXPort      - write select results to file
   SHOw <elt>  - print all data for one element [Rx]
   MODify <elt>- interactively edit one element  [Rx]
```

#### SOURCE & WAVELENGTH

```
   CHIefray    - use chief ray (not centroid) as
                 the FEX reference point
   WLENS       - list Rx optimization wavelengths
   SWL         - select wavelength by index
   MULtispec   - multi-color image stacking
   NFIlt, RFIlt, SFIlt - load / read / save
                 MULtispec filter data
   ATMosphere  - apply atmospheric phase screen
   SETC        - set design-optimization params
   SAOpt       - source-adjustment mode
```

#### RAY TRACING

```
   RAY <i>     - trace ray i; print state at each
                 surface                       [Rx]
   SEGRAYTrace <i> - trace a ray through segment i s
                 center (uses RptElt)          [Rx]
   PRAy, RRAy  - find pupil / marginal rays
   TPR         - trace the PRAy/RRAy rays
   MAP         - print ray-to-segment map
   COOrd <elt> - beam coordinates at element    [Rx]
   ACOor       - all-element coords (xyzLocal.txt)
   SYSprop     - first-order system properties  [Rx]
   EFL         - system effective focal length /
                 f-ratio / plate scale          [Rx]
   BWK         - chief-ray walk vs nominal      [Rx]
   FEXit       - find exit pupil
   SXP         - set exit pupil at an element [Rx]
   XPS         - exit-pupil surface fit over the
                 full ray grid              [Rx]
   FSR         - find source ray nearest an
                 aperture position         [Rx]
   CENter      - center beam in system stop
   STOp        - define system stop element
   CENTRoid    - use centroid as ref point for FEX
   FFP <x,y>   - find field point that places chief
                 ray at element coordinate
   PFP <x,y>   - find field point that places chief
                 ray at pixel coordinate
   FDP         - find detector plane orientation
   SPCenter    - spot-diagram centering option
```

#### BEHAVIOR TOGGLES

```
   LNEg/NOLNeg - allow/disallow negative pathlengths
   POLarized   - polarization tracing on
   NOPolarization - polarization off
   OBS         - spot-diagram obscuration option
   REGrid/NOREGrid - resample wavefront at element
                     during propagation
   ORS         - optimize reference surface
   SRS         - slave one ref surface to another
   NONe        - turn off plotting (no PGPLOT calls)
```

#### SURFACE DATA

```
   SINt        - set up interpolated-surface data
   UDSinit     - set up user-defined-surface data
```

#### PERTURBATION

```
   PERturb <elt>   - apply 6-DOF perturbation     [Rx]
   GPERturb <grp>  - perturb an element group     [Rx]
   PREad <file>    - apply perturbations from file
```

#### LINEAR MODEL

```
   BUild       - build linear model from chief-ray
                 ray partials                     [Rx]
   DMBuild     - build with deformable mirror     [Rx]
   PARtials <i>- print partial derivatives, ray i [BLD]
   LPErturb    - linear-model perturb              [BLD]
   LPRead      - linear perturb from file          [BLD]
   LREset      - zero linear perturbations         [BLD]
   LSPot, LOPd - linear-model spot / OPD           [BLD]
   LPIxilate   - linear-model pixelated image      [BLD]
```

#### DIFFRACTION

```
   PROpagate   - propagate diffraction wavefront [Rx]
   BEAm        - set beam type (Gauss, plane, ...)
   VECtor, SCAlar - vector vs scalar diffraction
                    (with POLarized)
   COMpose     - start a multi-object/color image
   ADD, DAdd   - add to / display composed image
   CADD        - add current intensity to composed
   NOIse, SEEd - image noise + RNG seed
   STRetch     - LIN / LOG / SQRT display stretch
```

#### OUTPUTS (RAY-TRACE & DIFFRACTION)

```
   OPD         - wavefront error map             [Rx]
   SPOt        - spot diagram                    [Rx]
   ZABerr      - Zernike OPD aberration          [Rx]
   ZCOef       - Zernike coefficient extraction  [Rx]
   LOS         - line-of-sight output
   METcalc     - metrology beam output
   INTensity   - wavefront intensity            [DIFF]
   PIXel       - pixelated intensity output     [DIFF]
   AMPlitude   - amplitude output               [DIFF]
   PHAse, PH   - phase output                   [DIFF]
   REAl        - real & imaginary wavefront     [DIFF]
   LOG         - log10 intensity                [DIFF]
   MTF, CMTF   - modulation transfer function   [DIFF]
   IMG         - load saved image data
```

#### PLOT STYLE

```
   GRAy        - gray-scale surface plot
   WIRe        - wireframe plot
   SLIce       - slice surface plot
   COLumn      - column line plot
   CONtour     - contour surface plot
   IMGmode     - toggle image polarity
                 (NEG=astro, POS=conventional)
   CIR, GIR    - color / gray image rendering
   PGP         - PGPLOT panel layout (1,2,3,4,9)
```

#### FILE OUTPUT

```
   TEXt        - ASCII print output
   BINary      - binary image file
   FITs, WFIts - FITS / wide-FITS image file
   MAT         - MATLAB .mat file
   GETMatlabmatrix - fetch named matrix into MATLAB
                     (smacos only)
```

#### WINDOW & POST-PROCESSING

```
   PLOcate     - set plot pixel-array location
   NOPLOC      - clear plot location
   PWIn        - plot window definition
   SZCo        - sizing / coordinate options
   GBS         - Gaussian beam-size post-processing
   BLUr, GBLur - blur / Gaussian blur
   GAIn        - apply gain map to intensity
   ODRaw, PGD  - outline drawing / PGPLOT diagnostics
   ROW         - extract / plot row from image
```

#### SYSTEM OPTIMIZATION

```
   AVAR <elt>  - add a variable element
   MVAR <elt>  - modify a variable element
   DVAR <elt>  - delete a variable element
   VARS        - list all variable elements
   AFOV <fov>  - add a field of view
   DFOV <fov>  - delete a field of view
   FOVS        - list all fields of view
   CALib       - run system optimization
   SFOV        - select field of view
```

#### MISC / DEBUG

```
   SRAy, SRT, SRTrace - older single-ray variants
                        (use RAY)
   DOPD, DOPL, ZRM     - debug / diagnostic prints
   JWST_v3d, Vis3d     - 3D rendering demo
```
