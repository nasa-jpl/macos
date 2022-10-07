# macos
MACOS source distribution
Jet Propulsion Laboratory
09/02/2022
Updated 09/30/2022

Modeling and Analysis for Controlled Optical Systems

The MACOS source code files and Linux Makefiles can be used to

1) Compile a MACOS executabale program on a Linux workstation
2) Compile a SMACOS library smacos_lib.a to be linked to a Matlab
   mex function module. The SMACOS mex module provides an interface
   to MACOS operations in a Matlab session.

On a Linux machine, "make macos" command generates a MACOS executable;
"make smacos" command generates a SMACOS library. 

To clean out existing complied modules of MACOS and SMACOS, say for
re-compilations, do 
"make clean-macos" for MACOS
"make clean-smacos" for SMACOS

The Makefiles use Intel Fortran and C compilers on a Linux platform
for MACOS compilations. Also included are sub-folders of source files for
1) a pgplot library for graphics plotting
2) a fits library for fits format input and output

At JPL, the compilations can be done on 
mustang2@jpl.nasa.gov
