#!/bin/bash
## exit on the first error.
set -e

# source intel iccvars file which will be used to compile all libraries
# source /opt/intel-$intel_version/bin/compilervars.sh intel64
source /opt/intel/bin/compilervars.sh intel64
cd macos_f90

#--------------------------------------------------------------------------
# Compile npsol library (can be compiled with gfortran or ifort)
#--------------------------------------------------------------------------
# cd to npsol/blas and compile
cd npsol/blas
make -f Makefile_Intel clean; make -f Makefile_Intel
cd ../

# cd to lappack and compile
cd lapack; make -f Makefile_Intel clean; make -f Makefile_Intel
cd ../

# now compile npsol
make -f Makefile_Intel clean; make -f Makefile_Intel

# cd back to /macos_f90 folder
cd ../

#--------------------------------------------------------------------------
# Compile readline-8.2 (uses gnu gcc to compile)
#--------------------------------------------------------------------------
# cd to readline-8.2, configure, and compile
cd readline-8.2
./configure
make

# cd back to /macos_f90
cd ../

#--------------------------------------------------------------------------
# Compile pgplot (Uses Intel Compiler)
#
# Note: we need to change libpgplot.so to sv_libpgplot.so
#--------------------------------------------------------------------------
#cd to pgplot and compile
cd pgplot
make
#macos compilation expects sv_libpgplot.so
mv libpgplot.so sv_libpgplot.so 
cd ../

#--------------------------------------------------------------------------
# Compile fitsio 
#--------------------------------------------------------------------------
cd fits_build
bash cmp
bash amp
cd ../

#--------------------------------------------------------------------------
# Compile macos and smacos
#--------------------------------------------------------------------------
make clean-macos
make macos

#make clean-smacos
browse confirm wa
make smacos



#--------------------------------------------------------------------------
# Requirements
#--------------------------------------------------------------------------
# Luis used following versions of compilers for compiling macos on "marchen":
#bash-4.2$ ifort --version
#ifort (IFORT) 13.1.3 20130607
#Copyright (C) 1985-2013 Intel Corporation.  All rights reserved.

#bash-4.2$ gcc --version
#gcc (GCC) 4.8.5 20150623 (Red Hat 4.8.5-39)
#Copyright (C) 2015 Free Software Foundation, Inc.
#This is free software; see the source for copying conditions.  There is NO
#warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

#bash-4.2$ gfortran --version
#GNU Fortran (GCC) 4.8.5 20150623 (Red Hat 4.8.5-39)
#Copyright (C) 2015 Free Software Foundation, Inc.












