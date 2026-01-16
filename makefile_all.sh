#!/bin/bash
## exit on the first error.
set -e

# source intel iccvars file which will be used to compile all libraries
# source /opt/intel-$intel_version/bin/compilervars.sh intel64
source /opt/intel-14.0.0/oneapi/setvars.sh intel64 --force

#export LD_LIBRARY_PATH=/opt/intel/oneapi/mkl/latest/lib/intel64
export LD_LIBRARY_PATH=/opt/intel-14.0.0/oneapi/compiler/2023.1.0/linux/compiler/lib/intel64

export LD_LIBRARY_PATH=/home/lmarchen/MACOS/work/git_macos/macos-develop2/macos_f90/pgplot
export PGPLOT_FONT=/home/lmarchen/MACOS/work/git_macos/macos-develop2/macos_f90/pgplot/grfont.dat
export PGPLOT_DIR=/home/lmarchen/MACOS/work/git_macos/macos-develop2/macos_f90/pgplot
export LD_LIBRARY_PATH=/home/lmarchen/MACOS/work/git_macos/macos-develop2/macos_f90/readline-8.2
export LD_LIBRARY_PATH=/opt/intel-14.0.0/oneapi/compiler/2023.1.0/linux/compiler/lib/intel64


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
# Compile GMI
#
# Note:
#
# 1. User must change text /home/lmarchen/MACOS/work/git_macos/macos-develop2
#    to user's local path to macos source code
# 2. User must change /usr/local/MATLAB/R2023b if using different version of
#    Matlab to compile GMI or if location of Matlab is different
# 3. User must change /opt/intel-14.0.0/oneapi/compiler/2023.1.0/linux/compiler/lib/intel64_lin 
#    if using a different version of the compiler or if different from it.
#-------------------------------------------------------------------------
cd ../GMI



ifx -C -traceback -fstack-protector -c  -I/home/lmarchen/MACOS/work/git_macos/macos-develop2/macos_f90 -I/usr/local/MATLAB/R2023b/extern/include -I/usr/local/MATLAB/R2023b/simulink/include -nologo -fpic -fpp -132 -gen-interfaces -fp-model strict -fno-omit-frame-pointer -D__amd64 -module /home/lmarchen/MACOS/work/git_macos/macos-develop2/macos_f90/SMACOS_OBJS/Linux-x86_64  -DGMI_SVN_REV="''" -DGMI_DATE="'2024-01-18'"  -DMX_COMPAT_32 -O2 -xHOST  "GMI.F"

ifx -C -traceback -fstack-protector -c  -I/home/lmarchen/MACOS/work/git_macos/macos-develop2/macos_f90 -I/usr/local/MATLAB/R2023b/extern/include -I/usr/local/MATLAB/R2023b/simulink/include -nologo -fpic -fpp -132 -gen-interfaces -fp-model strict -fno-omit-frame-pointer -D__amd64 -module /home/lmarchen/MACOS/work/git_macos/macos-develop2/macos_f90/SMACOS_OBJS/Linux-x86_64  -DGMI_SVN_REV="''" -DGMI_DATE="'2024-01-18'"  -DMX_COMPAT_32 -O2 -xHOST  "GMIG.F"

ifx -C -traceback -fstack-protector -O -shared-intel -shared -Wl,--version-script,/usr/local/MATLAB/R2023b/extern/lib/glnxa64/fexport.map  -Wl,--no-undefined -o  "GMI.mexa64"  GMI.o GMIG.o   -Wl,-rpath-link,/usr/local/MATLAB/R2023b/bin/glnxa64 -L/usr/local/MATLAB/R2023b/bin/glnxa64  -l:libmx.so -l:libmex.so -lmat -L/opt/intel-14.0.0/oneapi/compiler/2023.1.0/linux/compiler/lib/intel64 -lirc -lm -lstdc++  /home/lmarchen/MACOS/work/git_macos/macos-develop2/macos_f90/SMACOS_OBJS/Linux-x86_64/smacos_lib.a




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


