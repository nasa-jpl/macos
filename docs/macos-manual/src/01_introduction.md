## SECTION 1: Introduction

This section introduces the features of MACOS and gives an overview of the manual. New features introduced with this version of the code are summarized.

### Modeling and Analysis for Controlled Optical Systems

The Modeling and Analysis for Controlled Optical Systems software provides tools for modeling and analyzing optical systems. It consists of a stand-alone application program called MACOS and a Fortran subroutine package called SMACOS. SMACOS provides the full functionality of the MACOS application in a form that can be directly integrated with other codes.

MACOS is an optical analysis and model generation program, focused on generating and exercising fast, accurate, and detailed models combining ray-trace, differential ray-trace, and diffraction techniques. It is not intended as a replacement for a traditional optical design code, though it does include optimization features that can be used for design.

MACOS provides unique modeling capability for *system-level* design and analysis tasks. SMACOS functions can be directly integrated into simulation programs, to provide optical outputs driven by structural and/or thermal dynamics, feedback controls, atmospheric effects, or other system perturbations – the GMI Matlab MEX function provides an example of how to do this using the SMACOS API. More commonly, MACOS is used to generate Linear Optical Models (LOMs), either directly in MACOS, or driving SMACOS from Matlab or Python. MACOS LOMs are easily integrated with standard structures, controls, and user written analysis tools to create end-to-end system models.

This release provides *mmacos* and *pymacos*, Matlab and Python code libraries that call SMACOS to perform many functions. The *mmacos* functions include designing optical systems, adding structural interfaces, segmenting mirrors, adding metrology systems, generating LOMs, and performing analysis, while *pymacos* is focused on regression testing. The *mmacos* library includes design templates for telescopes, spectrometers, lab benches, coronagraphs, interferometers and other optical devices.

This release also provides a utility for conversion of prescriptions developed in the CodeV optical design code, into MACOS prescriptions. Once converted, MACOS is used to add modeling features, including mechanical and control system interfaces, optical figure aberrations, segmented mirrors, actuated mirrors, deformable mirrors, and reference and return surfaces for diffraction calculations. Then LOMs can be generated using MACOS functions, *mmacos* sensitivity generators, or SMACOS calls in codes like GMI. By exporting the LOMs into external analysis environments like Matlab, or by linking SMACOS-based code drivers, MACOS optical modeling functions can be driven by other calculations, from thermal-structural or controls codes, to compute controlled optical system response. The integrated models are then exercised in end-to-end dynamic simulations, Monte Carlo analyses, covariance analyses, or optimization studies to determine integrated system-level performance and guide system engineering and design.

MACOS (and SMACOS) are built on computationally efficient general ray-trace capabilities. These support standard analyses of optical performance such as ray-traces, spot-diagrams, and wavefront (OPD) plots. MACOS ray traces can be used to provide the phases for single- or multi-plane diffraction models: single-plane for simple imaging Point-Spread Function (PSF) generation; multiplane for more complex instruments such as coronagraphs. MACOS also provides polarization modeling, through polarization ray-trace and vector diffraction options. For diffraction calculations, MACOS uses Fresnel propagators driven by exact phases computed by its ray-trace engine. MACOS supports simulation of multiple source and wavelength images. Detector pixilation can be set to model detector resolution limits.

### Overview of this Manual

This manual provides two services to the user:

-   A technical orientation outlining MACOS ray-trace, differential ray-trace, and diffraction modeling capabilities.

-   Detailed description of MACOS commands including reference material and numerous examples.

Section 2 of this manual briefly summarizes the technical areas covered by MACOS: geometric optics and physical optics. The purpose is to establish basic definitions and to introduce certain MACOS features. More detailed technical discussion is provided at the beginning of each subsequent section and in the references. Section 3 discusses the basic procedures of starting MACOS, getting help, entering commands, and using files. Sections 4 through 9 proceed in the sequence that a typical MACOS session might:

-   Entering the optical prescription from a drawing of the system (Section 4)

-   Ray-tracing the system (Section 5)

-   Setting up the model for diffraction analysis (Section 6)

-   Simulating images (Section 7)

-   Building and exporting linear models (Section 8)

-   Generating and using SMACOS models (Section 9)

Several example problems are worked in the manual to illustrate specific features of the code. Appendix A contains the prescriptions and journal files for the examples. Appendix B contains selected technical papers.

This manual is the *learning* document; its companion, the **MACOS Command Reference** (`macosCmdRef.pdf`, source in `docs/macos-manual/cmdref/`), is the *lookup* document: every interactive command, and every function of the MATLAB (*mmacos*) and Python (*pymacos*) programming interfaces, cataloged one at a time with its syntax. Appendix C of this manual lists the commands in summary form; the Command Reference is the full per-command catalog and the place to start if you are driving MACOS from MATLAB or Python.

Three companion tools are worth knowing from the outset. **SegMirMaker** (`MACOS_resources/segmirmaker/`) converts a monolithic primary mirror into a many-segment one, generating the segmented-mirror prescription and segment-influence files from a parent-mirror description (see §4.4.11). **GMI** (`MACOS_resources/GMI/`) is a MATLAB MEX front end built directly on the SMACOS subroutine interface — both a complete forward simulator and the reference example of a real SMACOS application (see Section 9). The rx_converter, `cv2macos.seq`, can be used in CodeV to generate a MACOS prescription – a “.in-file” – from a CodeV prescription.

### MACOS and AI Agents

MACOS is increasingly driven by AI coding agents as well as by people, and the parts of the system an agent uses are the documented ones: the interactive program, which an agent runs on a pseudo-terminal and feeds command by command exactly as a user types them (the command-line test suite in `cli_tests/` is built that way, one journal of commands and expected output per command); the MATLAB and Python bindings (`mmacos`, `pymacos`), whose functions are the programmatic surface; and the Command Reference, which states every command's syntax, prompts and defaults and is the authority when prose and transcript disagree. Appendix D collects what an agent needs to drive MACOS reliably: where the facts live, how to run and read the program, the conventions that trip newcomers, and the working rules the project holds its own agents to. The repositories also carry agent-facing notes beside the code (`CLAUDE.md` at the top of each repository and in each subsystem directory), which are maintained with the code and are part of the release.

### New Features of Version 4.00

MACOS 4.00 (April 2026) is a feature release built on the Fortran 90 / SMACOS infrastructure introduced in v3.2. It adds a new general-purpose composite surface (FreeForm), adds multiple different forms of Zernike polynomials describing for optical surfaces, modernizes the graphics layer, expands the SMACOS / GMI MATLAB interface, and adds new utility commands plus a separate segmented-mirror prescription generator.

A summary of changes since MACOS 3.2:

-   FreeForm composite surface (SrfType = 14). Single surface type combining a conic base with a Cartesian monomial sag (Mon), a free-form Zernike sag (FF), and an interpolated grid sag, all on independent coordinate frames. See §4.5.11.

-   Per-element Zernike inputs. Sparse "modes + coefs" form MonZernCoef / MonZernModes and FFZernCoef / FFZernModes is accepted directly in .in files; conversion to monomial form happens internally at trace time. See §4.5.4.

-   Per-ray status tracking. Every ray records why it failed (Miss / Obscured / Bracket / MaxIter) and at which element. The end-of-trace WARN summary breaks down failures by category and throttles bracket / iter messages after the first 20.

-   Graphics modernization. Built-in giza graphics backend replacing legacy PGPLOT (still selectable at build time). New IMGMODE command toggles between astronomical (negative, large-\>dark) and conventional (positive, large-\>light) display polarity. Colorbar wedge labels reflect the active STRETCH (LIN / LOG10 / SQRT). See §7.x.

-   SEGRAYTRACE command. Trace a single ray through the geometric center of a chosen segment element; print ray state at every downstream surface. Useful for verifying segmented-mirror prescriptions. See §5.2.

-   SegMirMaker external tool. Standalone utility, in MACOS_resources/segmirmaker/, that generates segmented-mirror .presc and Hx.m files. Replaces SMPGe (1992) and adds support for FreeForm parents. See §4.4.11.

-   GMI MATLAB interface upgrades. New 14th MEX argument pmonzern for FreeForm MonZernCoef perturbations (parallel to pzern). New param.monzernSrf element list. Nominal-state save/restore extended to cover FreeForm coordinate frames so per-element rigid-body perturbations no longer leak across GMI calls. See §9.

-   Graceful Ctrl-C. SIGINT installs a clean-exit handler at start-up; pressing Ctrl-C at the MACOS\> prompt prints 'MACOS: interrupted' and exits cleanly, replacing the Intel-runtime 'forrtl: severe (69)' traceback.

-   Design optimization (CALIB). The native multi-field optimizer, never described in the 3.2 manual, is documented: the CALib command and its companions (SETC, AVAR, MVAR, DVAR, VARS, AFOV, DFOV, SFOV, FOVS, WLENS) run a Levenberg-Marquardt least-squares solve over the rigid-body, radius, conic, aspheric and Zernike degrees of freedom of the variable elements, scored at up to 12 fields and 6 wavelengths against a wavefront, spot, Zernike-mode or beam target, with bounds through the SLSQP driver.  Section 4 (Table 13) is the keyword reference; the Command Reference carries the dialogs; the setup is written by SAVE and announced at every load.

-   Metrology (MET). The metrology model -- per-element metrology points and target surfaces (`nMetPos=`, `tMetElt=`, `EdgeSensors=`) and the METCALC command that computes every metrology beam length and the edge-sensor readings -- is in the Command Reference (METcalc, which lists the prescription keys); the MATLAB design package's `add_met`, `edge_sensors` and `dmet_dx` build and differentiate metrology layouts on a prescription.

-   Design-optimization keyword reference.  The `Opt*` / `VarElt=` prescription keywords that set up CALib are tabulated with their defaults in Section 4 (Table 13), together with the solver's stopping rule.

-   HELP command rebuilt. All \~140 dispatched commands are now listed (vs. \~75 in v3.2), grouped into 15 categories with one-line effect descriptions and prerequisite tags (\[Rx\] / \[BLD\] / \[DIFF\]). See §3.2 and Appendix C.
