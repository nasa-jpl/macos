## SECTION 3: User Interface

We begin this section by discussing the mechanics of interacting with the MACOS application. Basic file-oriented functions are discussed, as are various graphics options. Memory requirements are discussed. The commands covered in this section are:

    Help EXPort JOUrnal EXEcute

### Starting and Stopping MACOS

The particular command used to start the MACOS application program depends—of course—on which machine you are using. If your machine is a…

-   Macintosh: open the folder containing the MACOS application and double-click on the MACOS 1.0 icon

-   Unix: type the command MACOS with appropriate path names to point to the MACOS executable file

On all machines, the response is the same. A window appears with the MACOS command dialog:

    Modeling and Analysis for Controlled Optical Systems MACOS Version 2.8, May 31, 1999
    This run is limited to 70 surfaces, 69 segments, 4226
    rays ( 441 for sensitivity calcs). Up to 3 diffraction planes are stored. Up to 10 interpolated surfaces of up to 255 data points are supported.
    Maximum pixel array size is 128 by 128. Diffraction grid size is 128 by 128.
    Glass table /home/comp/macosData/macos.glass read, with 202 glasses included. For a list of commands, type H

MACOS\>

MACOS\>**q**

The MACOS\> prompt indicates readiness of the MACOS application to take commands. To exit MACOS, type Quit or BYE to end a session. The response you get depends on whether a PGPLOT graphics window is open (Unix machines running XWindows only). If so, MACOS will ask you to hit return before it goes away. On Sun machines a warning that “non-standard” floating point results may have been produced will be printed. This can be ignored.

Type \<RETURN\> for next page:

### Help

Type HELP (or H) at the MACOS\> prompt to see the list of all commands. The output is paginated: hit return to advance to the next screen. Commands are grouped into 15 categories:

-   Session & files: QUit, END, HELP, REset, STatus, SUMmarize, ELTS, EXEcute, JOUrnal, MREset, !, PWD, CD, RX/LS/LL, VI/EMAcs

-   Prescription I/O: NEW, OLD, FID, VALidate, SAVe, EXPort, SHOw, MODify

-   Source & wavelength: CHIefray, WLENS, SWL, MULtispec, NFIlt/RFIlt/SFIlt, ATMosphere, SETC, SAOpt

-   Ray tracing: RAY, SEGraytrace, PRAy/RRAy/TPR, MAP, COOrd, ACOor, EFL, BWK, FEXit, FSR, CENter, STOp, CENTRoid, FFP, PFP, FDP, SPCenter

-   Behavior toggles: LNEg/NOLNeg, POLarized/NOPol, OBS, REGrid/NOREGrid, ORS, SRS, NONe

-   Surface data: SINt, UDSinit

-   Perturbation: PERturb, GPERturb, PREad

-   Linear model: BUild, DMBuild, PARtials, LPErturb, LPRead, LREset, LSPot, LOPd, LPIxilate

-   Diffraction: PROpagate, BEAm, VECtor, SCAlar, COMpose, ADD, DAdd, CADD, NOIse, SEEd, STRetch

-   Outputs: OPD, SPOt, ZABerr, ZCOef, LOS, METcalc, INTensity, PIXel, AMPlitude, PHAse/PH, REAl, LOG, MTF/CMTF, IMG

-   Plot style: GRAy, WIRe, SLIce, COLumn, CONtour, IMGmode, CIR, GIR, PGP

-   File output: TEXt, BINary, FITs/WFIts, MAT, GETMatlabmatrix

-   Window & post-processing: PLOcate, NOPLOC, PWIn, SZCo, GBS, BLUr/GBLur, GAIn, ODRaw, PGD, ROW

-   System optimization: AVAR, MVAR, DVAR, VARS, AFOV, DFOV, FOVS, CALib, SFOV

-   Misc / debug: SRAy, SRT, SRTrace, DOPD, DOPL, ZRM, JWST_v3d / Vis3d

Casing convention. Commands are shown as \<UPPER prefix\>\<lower tail\>, e.g. RAYtrace, SUMmarize, PERturb. The uppercase part is the minimum-match abbreviation tested by the dispatcher - typing just RAY, SUM, or PER (case-insensitive) is sufficient.

Prerequisite tags. Lines in the HELP output may end in:

\[Rx\] - needs a loaded prescription (run OLD or NEW first)

\[BLD\] - needs a built linear model (run BUild or DMBuild)

\[DIFF\] - needs a propagated wavefront (run PROpagate first)

Commands without a tag work in any state. A complete reference of all commands appears in Appendix C.

### MACOS Model Size Specification

The first thing the user will be asked to do after starting MACOS is to specify a model size. The model size determines MACOS computational accuray, such as the number of rays traced and image grid resolution, and obviously also has an impact on the overall computer simulation time. MACOS currently supports the following model sizes: 128, 256, 512, 1024 and 2048 mainly because the fast Fourier transorm routine is optimized for power-of-2 arrays.

The MACOS model size can be changed at any time during a MACOS session. Use the *mreset* command to set a new model size, as shown in the following example to set a model size to 256:

MACOS\> **mreset 256**

    or

MACOS\> **mreset**

    Enter new size (128, 256, 512, 1024, 2048): [128]: 256
    Modeling and Analysis for Controlled Optical Systems MACOS Version 3.2, February, 2003
    This run is limited to 70 surfaces, 69 segments, 16385
    rays ( 500 for sensitivity calcs). Up to 3 diffraction planes are stored. Up to 10 interpolated surfaces of up to 255 data points are supported.
    Maximum pixel array size is 256 by 256. Diffraction grid size is 256 by 256.

In previous versions of MACOS, separate model parameter files are used for different model sizes. For example, the model parameter file *param256.cmn* is used for the model size 256. In MACOS 3.2b, a single model parameter file called *macos_param.txt* is used for all model sizes. *macos_param.txt* contains parameter sets for all supported model sizes. To run MACOS, either *macos_param.txt* or a link to it must exist in the MACOS execution directory.

### Entering Commands

MACOS\>**old**

MACOS commands can be abreviated using the first letters of the command. The neccessary letters are denoted in capital letters. However, they do not need to be entered in capitals as MACOS is case insensitive.

MACOS prompts the user for the next input. For example, as shown below, to enter an old file, the user types OLD. MACOS then prompts for the filename.

    Enter file name: lens_array
    Input file lens_array.in being loaded. Default used for SegXgrid
    Optical Train Summary:
    Tracing rays past 4 elements.
    Aperture grid type = Circular with 512 grid points.

Elt 1: LensArray : Conic with Kc= 0.0000D+00, Kr=-8.0000D-01, nECo-ords = -6

Elt 2: Refractor : Conic with Kc= 0.0000D+00, Kr=-1.0000D+18, nECo-ords = -6

    Elt 3: FocalPlane : Flat with Kc= 0.0000D+00, Kr=-1.0000D+18, nECo-ords = -6

**Entering Commands**

    Elt 4: FocalPlane : Flat with Kc= 0.0000D+00, Kr=-1.0000D+18, nECo-ords = -6
    Number of output coordinates is 5
    Too many grid points. Resetting npts to 128 MACOS>

MACOS also allows commands to be stacked on a single input line. The commands must be separated by a space.

MACOS\>**old lens_array show 2**

    Enter file name: lens_array
    Input file lens_array.in being loaded. Default used for SegXgrid
    Optical Train Summary:
    Tracing rays past 4 elements.
    Aperture grid type = Circular with 512 grid points.

Elt 1: LensArray : Conic with Kc= 0.0000D+00, Kr=-8.0000D-01, nECo-ords = -6

Elt 2: Refractor : Conic with Kc= 0.0000D+00, Kr=-1.0000D+18, nECo-ords = -6

    Elt 3: FocalPlane : Flat with Kc= 0.0000D+00, Kr=-1.0000D+18, nECo-ords = -6
    Elt 4: FocalPlane : Flat with Kc= 0.0000D+00, Kr=-1.0000D+18, nECo-ords = -6
    Number of output coordinates is 5
    Too many grid points. Resetting npts to 128 MACOS> show
    Enter number of element (0=aperture): [0]: 2

+----------+--------------------+-------------------------------------+
| iElt=    | 2                  |                                     |
| EltName= |                    |                                     |
| Element= | rs Refractor       |                                     |
|          |                    |                                     |
| Surface= | Conic              |                                     |
+==========+====================+=====================================+
| fElt=    | 1.000000000D+18    |                                     |
+----------+--------------------+-------------------------------------+
| eElt=    | 0.000000000D+00    |                                     |
+----------+--------------------+-------------------------------------+
| KrElt=   | -1.000000000D+18   |                                     |
+----------+--------------------+-------------------------------------+
| KcElt=   | 0.000000000D+00    |                                     |
+----------+--------------------+-------------------------------------+
| psiElt=  | 0.000000000D+00    | 0.000000000D+00 -1.000000000D+00    |
+----------+--------------------+-------------------------------------+
| VptElt=  | 0.000000000D+00    | 0.000000000D+00 2.000000000D-01     |
+----------+--------------------+-------------------------------------+
| RptElt=  | 0.000000000D+00    | 0.000000000D+00 1.000000000D+00     |
+----------+--------------------+-------------------------------------+
| IndRef=  | 1.000000000D+00    |                                     |
+----------+--------------------+-------------------------------------+
| Extinc=  | 0.000000000D+00    |                                     |
+----------+--------------------+-------------------------------------+
| nObs=    | 0                  |                                     |
+----------+--------------------+-------------------------------------+
| ApType=  | None               |                                     |
+----------+--------------------+-------------------------------------+
| zElt=    | 0.000000000D+00    |                                     |
+----------+--------------------+-------------------------------------+
| P        | Geometric          |                                     |
| ropType= |                    |                                     |
+----------+--------------------+-------------------------------------+
| nECoord= | -6                 |                                     |
+----------+--------------------+-------------------------------------+

MACOS\>

MACOS\>**show;;**

MACOS provides defaults for nearly all requested parameters or values. These are usually computed from contextual information. The value of the default for a particular parameter is shown in the prompt line enclosed in square brackets. To indicate that the default is to be accepted, the user can simply hit return when the prompt is displayed. If typing ahead, a semicolon can also be used to indicate that a default should be accepted. In the example below, the default element number 0 is used to select the aperture.

    Enter number of element (0=aperture): [0]:
    ChfRayDir= 0.000000000D+00 0.000000000D+00 1.000000000D+00 ChfRayPos= 0.000000000D+00 0.000000000D+00 -1.000000000D+00 zSource= 1.000000000D+22
    IndRef= 1.000000000D+00 Extinc= 0.000000000D+00 Wavelen= 1.000000000D-04 Flux= 1.000000000D+00
    GridType= Circular Aperture= 5.000000000D-01 Obscratn= 0.000000000D+00
    nGridpts= 128
    xGrid= 1.000000000D+00 0.000000000D+00 0.000000000D+00 yGrid= 0.000000000D+00 1.000000000D+00 0.000000000D+00
    nElt= 4

MACOS\>

#### Interrupting MACOS (Ctrl-C)

Pressing Ctrl-C at the MACOS\> prompt - or during a long trace, build, or propagation - prints 'MACOS: interrupted' and exits with status 0. Internally MACOS installs a SIGINT handler at startup that uses \_exit(0) so it is async-signal-safe and fires regardless of what the runtime is doing. This replaces the Intel-runtime 'forrtl: severe (69)' traceback that earlier MACOS versions produced on Ctrl-C.

#### Empty input

Pressing return at the MACOS\> prompt with no command is a no-op - useful for clearing the input buffer or scrolling history without re-running the last command.

#### Unknown commands

Unknown commands print a single-line warning rather than the verbose multi-line traceback that v3.2 produced.

### Graphics Display

MACOS provides in-line graphical displays on Unix machines using the PGPLOT package from CalTech (see installation notes for information on how to obtain PGPLOT). PGPLOT is started up at the first invocation of a MACOS graphics function, such as INTensity. PGPLOT prompts the user to define the desired form for the graphical output:

MACOS\>**stretch**

    Enter image stretch type (LINEAR, LOG10, SQRT): [LINEAR]: log10

MACOS\>**intensity**

    Enter number of element where data is to be generated: [9]: 9
    Tracing 3210 rays and propagating 16384 grid points... Scalar FF Prop between Elt 8 and Elt 9 to WF 1:
    z1=-9.1481D+01 dx1= 1.0581D-01 min= 1.0581D-01 max= 1.0581D-01 dev= 2.6446D-06 lin=100.00%
    z2= 1.0000D+22 dx2=-1.0996D-03
    Compute time was 1.4453 sec
    Graphics device/type (? to see list, default /NULL): /xwindow
    Starting /opt/local/bin/pgxwin_server.

The normal selection for Unix machines is “/XWindow.” PGPLOT then sets up a graphics window in which plots are displayed. Once selected, the graphical output option cannot be changed without quitting and restarting MACOS.

### Off-Line Graphics Display

MACOS graphics can also be displayed and analyzed using other applications. This is the only display option for Macintosh. It is also extremely useful for users on all platforms who wish to study their results in ways not provided by MACOS.

One such application is Matlab®, which provides an easy-to-use, powerful environment for manipulating and displaying data in many forms. Matlab is a product of The Math Works, Inc., 24 Prime Park, Natick, MA, USA (508)653-1415. We have included Matlab “.m-file” functions that provide some of the standard MACOS display types. These scripts import a specified datafile, generate a plot according to specified plot and stretch types, and return the data values. They functions also provide templates for modification or translation to other languages. By keeping the display application (e.g. Matlab) running in one window, and MACOS running in another, nearly seamless generation and display of results can easily be accomplished.

Generating the data for off-line display is as simple as selecting a file output “plot type,” using the BINARY, FITS, TEXT or MATLAB commands (Section 7). Subsequent use of graphics-generating commands such as INTensity then create a file and store the requested data in that file. If the file already exists, the user is asked if the file is to be

**PostScript and GIF Graphics Output**

overwritten; if not, a new name is requested. Other applications can then access and display the data. Repeating the previous example:

MACOS\>**binary**

    Plot type set to BINARY MACOS>intensity
    Enter number of element where data is to be generated: [9]: 9
    Tracing 3210 rays and propagating 16384 grid points... Scalar FF Prop between Elt 8 and Elt 9 to WF 1:
    z1=-9.1481D+01 dx1= 1.0581D-01 min= 1.0581D-01 max= 1.0581D-01 dev= 2.6446D-06 lin=100.00%
    z2= 1.0000D+22 dx2=-1.0996D-03
    Compute time was 1.0508 sec Wavefront Propagation Data Summary:
    Wavelength= 1.6280000D-04; Source Flux= 1.0000000D+00 u-v plane diam= 1.3438007D+01 du= 1.0581108D-01
    x-y plane diam=-1.3965216D-01 dx=-1.0996233D-03
    Peak intensity= 1.1584734D-01; Peak occurs at i= 65, j= 65 Sum of intensity= 7.6908687D-01
    Writing Intensity, Elt=9
    DIRECT ACCESS File=tmo.int9 Data file is “tmo.int9”
    Binary array is of dimension 128 by 128. Magnitude of central pixel is 0.3161060810D-01

MACOS\>

The output file name (in this case, “tmo.int9”) is nominally generated from the prescription file name, with the datatype (“int” for intensity) and element number appended after a period. Fits, text and Matlab format files have “.fit,” “.txt,” or “.mat” appended as well, respectively.

The Matlab command required to read the file, generate the plot (at log10 stretch), and return the data for further analysis is:

    >>intensity=macosIMG(‘gray’,’tmo.int9’,’log10’);

Here intensity is the returned intensity data, in this case a 128 by 128 matrix, ‘gray’ is the plot type (gray scale), ’tmo.int9’ is the file name, and ’log10’ is the stretch option. The result is a plot functionally identical to what MACOS would produce with the same settings.

Other applications and environments that you may find useful for analyzing and displaying MACOS data include include SAOImage (available free from the Smithsonian Astrophysical Observatory for Unix only), Spyglass Transform, MatrixX, IDL, plotting packages such as DeltaGraph, and spreadsheets such as Excel.

### PostScript and GIF Graphics Output

MACOS provides PostScript and GIF graphical file output directly from PGPLOT (Unix only). When PGPLOT is started up, at the first invocation of a MACOS graphics function, it prompts the user to define the desired form for the graphical output. Typing

? prints a list of PGPLOT output options:

    Graphics device/type (? to see list, default /NULL): ?
    PGPLOT v5.0.0 Copyright 1994 California Institute of Technology Legal PGPLOT device types are:
    /GIF (Graphics Interchange Format file, landscape orientation)
    /VGIF (Graphics Interchange Format file, portrait orientation)
    /NULL (Null device, no output)
    /PS (PostScript file, landscape orientation)
    /VPS (PostScript file, portrait orientation)
    /CPS (Colour PostScript file, landscape orientation)
    /VCPS (Colour PostScript file, portrait orientation)
    /TEK4010 (Tektronix 4010 terminal)
    /RETRO (Retrographics VT640 Tek emulator)
    /XTERM (XTERM Tek terminal emulator)
    /XWINDOW (X window window@node:display.screen/xw)
    /XSERVE (A /XWINDOW window that persists for re-use)
    /XDISP (pgdisp or figdisp server)
    Graphics device/type (? to see list, default /NULL):

The GIF options generate separate GIF files for each graphic, named pgplot.gif, pgplot.gif_2, etc. in sequence.

The PS options generate a single file, pgplot.ps, containing all of the plots generated in a single session. It can be printed (on Unix machines) with the command:

    % lpr pgplot.ps

The file can be renamed if you wish to save it.

Screen capture utilities and applications, such as SAOimage or XV, can also be very useful in recording MACOS graphical outputs.

### MACOS Files

#### .In-File

An optical system prescription created in MACOS is saved in a “filename.in” file. This file stores the prescription in a text format which can be edited offline and reused later. The MACOS prescription files are generated by the NEW command (Section 4).

#### .Out-File

MACOS no longer generates a “filename.out” file for recording each result generated in a session.

#### Glass Tables

Glass data is provided in a text file, named macos.glass, which is located in a special MACOS data directory (see next section). The form and content of this file is discussed in Section 4.4.3.

Alternatively, glass data can be provided as part of the optical prescription “.in-file.”

#### MACOS Data Directory

Beginning with Version 2.8, certain MACOS data files, such as glass tables, are being placed in a special directory. On Unix systems, the path to this directory can be specified by the user by defining an environment variable, MACOS_DATA. This is done in response to the Unix shell prompt, as in:

    > setenv MACOS_DATA /home/macos/data/

**MACOS Files**

The above command establishes the directory /home/macos/data as the location of the data files for subsequent MACOS sessions. The user may specify any directory as MACOS_DATA.

If MACOS_DATA is not specified, 3 other directories will be searched for data files. In order these are:

1.  Directory macosData in the user’s home directory

2.  Directory /home/comp/macosData

3.  The current directory

If the required data file is not found in any of these directories, a warning message is issued. In this event, the user is required to supply all relevant data, such as glass coefficients, through the prescription file.

#### Data Export

The EXPort command writes certain *variables* to a file, in several different formats (Section 7 discusses writing *graphics* data to files). The variables are useful when linking linear models to thermal or structural codes (see Section 8).

The following variables can be exported:

-   CMAtrix: optical sensitivity matrix in global or local coordinates

-   CmDM: deformable mirror sensitivity matrix

-   RayPos: current position 3-vectors for each ray

-   RayDir: current ray direction 3-vectors for each ray

-   CumRayL: current optical pathlengths for each ray Files may be saved in the following formats:

-   BINary: data is saved to a Fortran binary file, filename.bin

-   TEXt: data is saved to an ASCII text file, filename.txt

-   MFILE: data is saved to a Matlab file, filename.m

The following is a sample EXPort session.

MACOS\>**build ;**

    Enter terminal element number: [1]:
    Tracing 12669 rays but not BUILDing linear model. BUILD limited to 500 rays. MOD npts and recompute. MACOS>export
    Current element is 1
    Enter EXPORT mode (BINary, TEXt, MFIle, NAStran, H or Quit): [BINARY]: tex
    Exporting data to lens_array.txt
    Export data for CHFRay only, ALLRays, or NOChfray? [CHFRay]: chfr
    EXPORT (CMAtrix, RayIndex, RAYPos, RAYDir, RAYL, H or Q): [QUIT]: raypos
    % RayPos 3-vectors for 1 rays at Element 1
    EXPORT (CMAtrix, RayIndex, RAYPos, RAYDir, RAYL, H or Q): [QUIT]: q MACOS>

### Using Macros

MACOS has a scripting or macro capability. This is provided by keystroke recording and play-back commands, JOUrnal and EXEcute. Macros can call other macros, up to 10 deep. Comments (anything following a “%” on a line) and blank lines are ignored. Looping, conditionals and other “language” capabilities are not currently supported.

#### Journal Files

The JOUrnal command is used to start and stop a journal of all the commands the user inputs to the terminal. Typing the JOUrnal command initiates keystroke recording, which continues until toggled off by typing the JOUrnal command a second time.

JOUrnal allows the user to keep track of the commands entered in a MACOS session. The resulting .jou-file can later be used as a macro command (see Section 3.9.2). The

.jou-file can also be edited and modified offline.

The JOUrnal command has one argument, the filename. The filename entered by the user will automatically have the file extension .jou appended, to denote that the file is a JOUrnal file. If the file already exists, then the contents of the file are deleted before the JOUrnalizing begins. It is imperative that the user select a unique filename, if is it desired that old JOUrnal files be kept.

#### Executing Macros

The EXEcute command allows the user to designate a .jou-file from which further MACOS commands are read. Each command (including other EXEcute commands) is executed in sequence, until an END OF FILE symbol is read from the .jou-file.

MACOS then returns to interactive mode and the user can once again enter commands from the keyboard.

The EXEcute command has one argument, the filename of the batch file. The filename must have the file extension .jou to denote that the file is an journal file.

The user can edit a file produced with the JOUrnal comand to produce a macro. For example, if the .jou-file was created with an OPD at the 10th element, it is trivial to modify it to calculate the OPD at the 14th element. The user only needs to change the 10 in the .jou-file to a 14 using a text editor.

Examples are provided in the macosDemo directory distributed with the MACOS code.

### Memory Requirements

The MACOS application currently uses a fixed memory allocation set at compile time. Different run modules are supplied for different sized problems, depending on the machine and available memory. Table 4 summarizes requirements for specific platforms.

**TABLE 4.** Memory Requirements

+---------------------+-----------------+-----------------------------+
| **executable**      | **Sun/Solaris** | **PowerMacintosh**          |
+=====================+=================+=============================+
| macos128            | 9.75 MB         |                             |
+---------------------+-----------------+-----------------------------+
| macos256            | 18.3 MB         |                             |
+---------------------+-----------------+-----------------------------+
| macos512            | 38.4 MB         |                             |
+---------------------+-----------------+-----------------------------+
| macos1024           | 115 MB          |                             |
+---------------------+-----------------+-----------------------------+
