#!/usr/bin/env bash
# makegfortran.sh  — CMake build for macos + smacos + smacos_dvr + GMI using gfortran/gcc
#
# Usage:
#   source ./makegfortran.sh [giza|pgplot] [debug|release]
#
# Options (order-independent, defaults: giza, release):
#   giza    — build Giza instead of pre-built PGPLOT (default)
#   pgplot  — link pre-built libpgplot.a + X11
#   debug   — -O0 -fcheck=all
#   release — -O2 (default)
#
# Examples:
#   source ./makegfortran.sh                      # giza, release
#   source ./makegfortran.sh debug                # giza, debug
#   source ./makegfortran.sh pgplot               # pgplot, release
#   source ./makegfortran.sh pgplot debug         # pgplot, debug
#
# Output locations after a default (giza, release) build:
#   macos executable  : ./build_release_giza_gfortran/bin/macos
#   smacos library    : ./build_release_giza_gfortran/lib/libsmacos.a
#   smacos_dvr        : ./build_release_giza_gfortran/bin/smacos_dvr
#   GMI mex           : ./build_release_giza_gfortran/lib/GMI.mexa64
#
# For other combinations substitute the build directory name:
#   build_release_pgplot_gfortran, build_debug_giza_gfortran, build_debug_pgplot_gfortran

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# --- Parse arguments ---
PLOT_CHOICE="giza"
BUILD_TYPE="Release"

for arg in "$@"; do
  case "$arg" in
    giza)    PLOT_CHOICE="giza"    ;;
    pgplot)  PLOT_CHOICE="pgplot"  ;;
    debug)   BUILD_TYPE="Debug"    ;;
    release) BUILD_TYPE="Release"  ;;
    *)
      echo "makegfortran.sh: unknown option '$arg'"
      echo "Usage: source ./makegfortran.sh [giza|pgplot] [debug|release]"
      return 1 2>/dev/null || exit 1
      ;;
  esac
done

BUILD_DIR="${SCRIPT_DIR}/build_$(echo "${BUILD_TYPE}" | tr '[:upper:]' '[:lower:]')_${PLOT_CHOICE}_gfortran"
echo "makegfortran: BUILD_TYPE=${BUILD_TYPE}  PLOT=${PLOT_CHOICE}"
echo "makegfortran: BUILD_DIR=${BUILD_DIR}"

# --- LD_LIBRARY_PATH ---
READLINE_LIB="${SCRIPT_DIR}/macos_f90/readline-8.2/shlib"
if [ "$PLOT_CHOICE" = "pgplot" ]; then
  PGPLOT_DIR="${SCRIPT_DIR}/macos_f90/pgplot"
  export LD_LIBRARY_PATH="${READLINE_LIB}:${PGPLOT_DIR}:${LD_LIBRARY_PATH}"
else
  export LD_LIBRARY_PATH="${READLINE_LIB}:${LD_LIBRARY_PATH}"
fi

# --- USE_GIZA flag ---
if [ "$PLOT_CHOICE" = "giza" ]; then
  USE_GIZA_FLAG="ON"
else
  USE_GIZA_FLAG="OFF"
fi

# --- Configure ---
cmake -S "${SCRIPT_DIR}" -B "${BUILD_DIR}" \
  -DCMAKE_Fortran_COMPILER=gfortran \
  -DCMAKE_C_COMPILER=gcc \
  -DCMAKE_BUILD_TYPE="${BUILD_TYPE}" \
  -DUSE_GIZA="${USE_GIZA_FLAG}" \
  -DBUILD_MACOS=ON \
  -DBUILD_SMACOS=ON \
  -DBUILD_SMACOS_DVR=ON \
  -DBUILD_GMI=ON \
  || { echo "makegfortran: cmake configure FAILED"; return 1 2>/dev/null || exit 1; }

# --- Build ---
cmake --build "${BUILD_DIR}" -j"$(nproc)" \
  || { echo "makegfortran: cmake build FAILED"; return 1 2>/dev/null || exit 1; }

# --- Report ---
echo ""
echo "makegfortran: BUILD SUCCEEDED"
echo "  macos      : ${BUILD_DIR}/bin/macos"
echo "  smacos     : ${BUILD_DIR}/lib/libsmacos.a"
if [ -f "${BUILD_DIR}/bin/smacos_dvr" ]; then
  echo "  smacos_dvr : ${BUILD_DIR}/bin/smacos_dvr"
fi
if ls "${BUILD_DIR}"/lib/GMI.mexa64 > /dev/null 2>&1; then
  echo "  GMI        : ${BUILD_DIR}/lib/GMI.mexa64"
fi
