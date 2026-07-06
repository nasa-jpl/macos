## SECTION 1: Introduction

This section introduces the features of MACOS and gives an overview of the manual. New features introduced with this version of the code are summarized.

### Modeling and Analysis for Controlled Optical Systems

The Modeling and Analysis for Controlled Optical Systems software provides tools for modeling and analyzing optical systems. It consists of a stand-alone application pro-gram called MACOS and a Fortran subroutine package called SMACOS. SMACOS provides the full functionality of the MACOS application in a form that can be directly integrated with other codes.

MACOS is an optical analysis and model generation program, as opposed to an optical design program. It is focused on generating and exercising fast, accurate, and detailed models combining ray-trace, differential ray-trace, and diffraction techniques. It does not include extensive internal optimization functions or other features to support the optical designer.

MACOS does provide unique modeling capability for *system-level* design and analysis tasks, however. Capabilities include modeling optics on dynamic structures, deformable optics, and controlled optics. These models are easily integrated with standard structures, controls, and user written analysis tools to create an end-to-end system model.

MACOS is commonly used in conjunction with a design code such as CodeV. The design code produces an initial optical prescription, which is converted into a MACOS prescription. Further features are added in MACOS, including mechanical and control system interfaces, optical figure aberrations, segmented mirrors, actuated mirrors, deformable mirrors, and reference and return surfaces for diffraction calculations. Using MACOS linear optical models (or SMACOS subroutines), models can be integrated with structures, dynamics, and control codes. These models are exercised in end-to-end dynamic simulations, Monte Carlo analyses, covariance analyses, or optimization studies to determine integrated system-level performance.

MACOS (and SMACOS) provide computationally efficient general ray-trace capabili-ties. These support standard analyses of optical performance such as ray-traces, spot-diagrams, and wavefront (OPD) plots. They also support MACOS diffraction calcula-tions.

MACOS has a differential ray-trace capability which generates linear matrix models of optical systems. These models are accurate for systems that undergo only small pertur-bations and are essential for tasks such as covariance analysis and control design. Mod-els can be exported to linear analysis programs such as Matlab where they can be exercised alone or in combination with structural, thermal, and control models.

For diffraction calculations, MACOS uses Fresnel propagators driven by exact phases computed by its ray-trace engine. Simulated images (point-spread functions) can be generated using single- or multi-plane diffraction. MACOS supports simulation of multiple source and wavelength images. Detector pixilation can be set to model detector resolution limits.

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

-   HELP command rebuilt. All \~140 dispatched commands are now listed (vs. \~75 in v3.2), grouped into 15 categories with one-line effect descriptions and prerequisite tags (\[Rx\] / \[BLD\] / \[DIFF\]). See §3.2 and Appendix C.

### New Features of Version 3.2 (historical)

-   Interpolated data surfaces using general X-Y-Z data are once again supported

-   Interpolated data surfaces using square X-Y grids have been added

-   Influence function-based surfaces are provided as a "user-defined" surface

-   Transmissive and reflective straight-ruled gratings were added

-   "Flower" segmentation was added

-   A gaussian blur function was added for simulating image jitter and other noise effects

-   SMACOS call line arguments were changed

-   OPMmat and RaySpot variables were changed to double precision

### Changes since MACOS 2.8 (historical)

MACOS 3.2 Beta Version (3.2b) is a **Fortran 90 version** of MACOS, developed using MACOS 2.86 as the baseline. From a user's point of view, the most important improvement in MACOS 3.2b is its ability to support runtime specification of the MACOS model size in both interactive mode and in subroutine MACOS (SMACOS). The usage of this feature is described in Section 3.3 for the interactive mode and in Section 9.3 for SMACOS. Taking advantage of the Fortran 90 compiler's strong type-checking, a number of type mismatches between calling and called subroutines in MACOS have been fixed. MACOS 3.2b has been tested on several optical system simulations for space telescope wavefront sensing and control.

MACOS 3.2b retains all improvements that MACOS 2.8 made over earlier versions of MACOS. For users familiar with MACOS before version 2.8, version 2.8 brings some improvements and some changes. New features include additional surface types, additional element types, a new segmentation option, and glass tables. Full backward compatibility with existing MACOS prescriptions and journal files has been preserved.

MACOS 2.8 made some changes in its interface to SMACOS functions. These include a changed SMACOS argument list; also some important arguments and internal arrays, such as OPMmat, are now REAL\*8 rather than REAL\*4 variables. In addition to these changes, a user program calling SMACOS in MACOS 3.2b needs to include relevant Fortran 90 modules instead of Fortran 77 common block files as a way of communicating with SMACOS. Programs calling SMACOS should be examined to ensure compatibility with Version 3.2.

A summary of changes since MACOS 2.8:

-   Toroid surface type was added
