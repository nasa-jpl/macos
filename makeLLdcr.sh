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
export intel64_lib="/opt/intel/oneapi/compiler/latest/lib"
source /opt/intel/oneapi/setvars.sh intel64 --force

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
cd ~/dev/macos/macos_f90
# Our source files are in the current directory
export macossrc_dir=$(pwd)
export intellib_dir=""
export matlab_version=/usr/local/MATLAB/R2025b

export PGPLOT_FONT=$macossrc_dir/grfont.dat     # :"$PGPLOT_FONT"
#export PGPLOT_DIR="/usr/local/pgplot":"$PGPLOT_DIR"
export PGPLOT_DIR=$macossrc_dir/pgplot/:"$PGPLOT_DIR"
export LD_LIBRARY_PATH=$macossrc_dir/readline-8.2:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH="/opt/intel/oneapi/compiler/latest/lib":$LD_LIBRARY_PATH
export LD_LIBRARY_PATH="/usr/local/pgplot":$LD_LIBRARY_PATH

#--------------------------------------------------------------------------
# Compile npsol library (can be compiled with gfortran or ifort)
#--------------------------------------------------------------------------
# cd to npsol/blas and compile
cd npsol/blas
make -f Makefile_Intel clean
make -f Makefile_Intel
cd ../

# cd to lappack and compile
cd lapack; 
make -f Makefile_Intel clean
make -f Makefile_Intel
cd ../

# now compile npsol
make -f Makefile_Intel clean
make -f Makefile_Intel

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
make clean
make

#macos compilation expects sv_libpgplot.so
cp libpgplot.so sv_libpgplot.so
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
# End
#--------------------------------------------------------------------------
cd ../
