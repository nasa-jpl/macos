#!/bin/bash
## exit on the first error.
set -e

# source intel iccvars file which will be used to compile all libraries
# source /opt/intel-$intel_version/bin/compilervars.sh intel64
#source /opt/intel/bin/compilervars.sh intel64
# The following line assumes compilation occurs on s383 Network on mustang4 
# for a different computer make sure to specify the correct path for setvars.sh
source /opt/intel-14.0.0/oneapi/setvars.sh intel64 --force

# User must also set paths in order for compilation and macos to run as follows (
# user can uncomment lines below if user would not like to update the .aliases file):
# export LD_LIBRARY_PATH=/opt/intel/oneapi/mkl/latest/lib/intel64
# export LD_LIBRARY_PATH=/home/lmarchen/MACOS/work/git_macos/macos-dev_wrk/macos_f90/pgplot
# export PGPLOT_FONT=/home/lmarchen/MACOS/work/git_macos/macos-dev_wrk/macos_f90/pgplot/grfont.dat
# export PGPLOT_DIR=/home/lmarchen/MACOS/work/git_macos/macos-dev_wrk/macos_f90/pgplot
# export LD_LIBRARY_PATH=/home/lmarchen/MACOS/work/git_macos/macos-dev_wrk/macos_f90/readline-8.2
# export LD_LIBRARY_PATH=/opt/intel-14.0.0/oneapi/compiler/2023.1.0/linux/compiler/lib/intel64_lin

#--------------------------------------------------------------------------
# cd to macos source code directory: macos_f90
#--------------------------------------------------------------------------
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

make clean-smacos
make smacos

#--------------------------------------------------------------------------
# Requirements
#--------------------------------------------------------------------------
# Luis used following versions of compilers for compiling macos on "marchen":
#bash-4.2$ ifort --version
#ifort (IFORT) 2021.9.0 20230302
#Copyright (C) 1985-2023 Intel Corporation.  All rights reserved.
#or
#ifx (IFX) 2023.1.0 20230320
#Copyright (C) 1985-2023 Intel Corporation. All rights reserved.

#bash-4.2$ gcc --version
#gcc (Ubuntu 9.4.0-1ubuntu1~20.04.2) 9.4.0
#Copyright (C) 2019 Free Software Foundation, Inc.
#This is free software; see the source for copying conditions.  There is NO
#warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

#bash-4.2$ gfortran --version
#GNU Fortran (Ubuntu 9.4.0-1ubuntu1~20.04.2) 9.4.0
#Copyright (C) 2019 Free Software Foundation, Inc.
#This is free software; see the source for copying conditions.  There is NO
#warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.













