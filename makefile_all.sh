#!/bin/bash
## exit on the first error.
set -e

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
## Compilation Instructions
#
# Instructions to compile MACOS, SMACOS, and necessary libraries such as 
# npsol, pgplot, and readline. Assumption is user has cd to the MACOS folder.
# Here we set some global variable which can serve as reference inside
# the Makefile to avoid the need for user to change or hack the MACOS 
# Makefile.
#
# To Compile:
#
#   1. Setup Environment
#   2. Setup MACOS source path with pwd
#   3. Compile npsol
#   4. Compile pgplot
#   5. Compile readline-8.2
#   6. Compile macos
#   7. Compile smacos
#   8. Compile GMI

#-----------------------------------------------------------------------------
# 1. Set up environment
# 
#  Note: User must set the version of the Intel Compiler being used. User may
#  also need to tinker with Intel Compiler path here because the path to
#  compiler depends on the setup of the Intel Compiler by Admin. In this case
#  we are using oneAPI which requires oneAPI folder inside inte-$intel_version
#  to source oneAPI Intel Compiler variables.
#------------------------------------------------------------------------------
export intel_version=14.0.0
export intel64_lib=/opt/intel-$intel_version/oneapi/compiler/2023.1.0/linux/compiler/lib/intel64_lin
source /opt/intel-$intel_version/oneapi/setvars.sh intel64 --force

#------------------------------------------------------------------------------------------------------
# 2. Setup MACOS source path with pwd, and Matlab path
#
# The user should add the source line above, and the export lines below to 
# their .aliases or .bashrc file:
#
#  1. export PGPLOT_FONT=$macossrc_dir/grfont.dat
#  2. export PGPLOT_DIR=$macossrc_dir/pgplot
#  3. export LD_LIBRARY_PATH=$macossrc_dir/readline-8.2
#  4. export LD_LIBRARY_PATH=$intel64_lib
#------------------------------------------------------------------------------------------------------
cd macos_f90
# Our source files are in the current directory
export macossrc_dir=$(pwd)
export matlab_version=$(matlab -e | sed -n 's/MATLAB=//p')  #/usr/local/MATLAB/R2013a #$(which matlab)

export PGPLOT_FONT=$macossrc_dir/grfont.dat:"$PGPLOT_FONT"
export PGPLOT_DIR=$macossrc_dir/pgplot:"$PGPLOTDIR"
export LD_LIBRARY_PATH=$macossrc_dir/readline-8.2:"$LD_LIBRARY_PATH"
#export LD_LIBRARY_PATH=$intel64_lib:"$LD_LIBRARY_PATH"

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
# Compile GMI
#
# Note:
#
# 1. User must change matlab_version above, or leave it as is if matlab
#    has been setup as /usr/local/bin/matlab
#
#-------------------------------------------------------------------------
cd ../../GMI
make clean
#make


eval $(ifx -C -traceback -fstack-protector -c  -I$macossrc_dir -I$matlab_version/extern/include -I$matlab_version/simulink/include -nologo -fpic -fpp -132 -gen-interfaces -fp-model strict -fno-omit-frame-pointer -D__amd64 -module $macossrc_dir/SMACOS_OBJS/Linux-x86_64  -DGMI_SVN_REV="''" -DGMI_DATE="'2024-01-18'"  -DMX_COMPAT_32 -O2 -xHOST  "GMI.F")

eval $(ifx -C -traceback -fstack-protector -c  -I$macossrc_dir -I$matlab_version/extern/include -I$matlab_version/simulink/include -nologo -fpic -fpp -132 -gen-interfaces -fp-model strict -fno-omit-frame-pointer -D__amd64 -module $macossrc_dir/SMACOS_OBJS/Linux-x86_64  -DGMI_SVN_REV="''" -DGMI_DATE="'2024-01-18'"  -DMX_COMPAT_32 -O2 -xHOST  "GMIG.F")

eval $(ifx -C -traceback -fstack-protector -O -shared-intel -shared -Wl,--version-script,$matlab_version/extern/lib/glnxa64/fexport.map  -Wl,--no-undefined -o  "GMI.mexa64"  GMI.o GMIG.o   -Wl,-rpath-link,$matlab_version/bin/glnxa64 -L$matlab_version/bin/glnxa64  -l:libmx.so -l:libmex.so -lmat -L/opt/intel-14.0.0/oneapi/compiler/2023.1.0/linux/compiler/lib/intel64 -lirc -lm -lstdc++  $macossrc_dir/SMACOS_OBJS/Linux-x86_64/smacos_lib.a)




#--------------------------------------------------------------------------
# Requirements to Compile on a Linux System 
#--------------------------------------------------------------------------
# Luis used following versions of compilers for compiling macos on "marchen":
# 1. bash (run on bash, any version)
# 2. Fortran Intel Compiler oneAPI (2024 or earlier down to 2013) ifort or ifx (ifx early 2024 won't compile GMI)
#    User must install oneAPI HPC Toolbox to compile Fortran and Intel oneAPIBase Toolkit. 
# 3. GCC (4.2 or later, used to compile supporting libraries)
# 4. Matlab 2013a to Latest (2023b, GMI will compile with any of these version with either ifort or ifx 
#    with the exception of the ifx early 2024 Fortran Intel Compiler


