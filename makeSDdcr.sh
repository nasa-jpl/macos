# Makefile.
# Usage: source ./makeSDdcr.sh [debug]
#   debug  — build with -O0 -check all for interactive debugging
# Requires prior make macos and make smacos (for .mod files and smacos_lib.a)

source /opt/intel/oneapi/setvars.sh intel64 --force

#------------------------------------------------------------------------------------------------------

if [ "$1" = "debug" ]; then
  export OPT="-O0 -check all"
  echo "*** Debug build: OPT=$OPT"
else
  unset OPT    # let Makefile default to -Ofast
fi

cd ~/dev/macos/macos_f90

export macossrc_dir=$(pwd)
export intellib_dir=""
export matlab_version=/usr/local/MATLAB/R2025b

export PGPLOT_FONT=$macossrc_dir/grfont.dat     # :"$PGPLOT_FONT"
#export PGPLOT_DIR="/usr/local/pgplot":"$PGPLOT_DIR"
export PGPLOT_DIR=$macossrc_dir/pgplot/:"$PGPLOT_DIR"
export LD_LIBRARY_PATH=$macossrc_dir/readline-8.2:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH="/opt/intel/oneapi/compiler/latest/lib":$LD_LIBRARY_PATH
export LD_LIBRARY_PATH="/usr/local/pgplot":$LD_LIBRARY_PATH

export intel64_lib="/opt/intel/oneapi/compiler/latest/lib"

#--------------------------------------------------------------------------
# Compile smacos-driver (builds in MACOS_OBJS; requires prior make macos)
#--------------------------------------------------------------------------

make smacos_dvr

#--------------------------------------------------------------------------
# End
#--------------------------------------------------------------------------
cd ~/dev/macos/