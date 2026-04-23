#!/usr/bin/env bash
# makejoint.sh  — CMake build for macos + smacos + smacos_dvr + GMI
#
# Usage:
#   source ./makejoint.sh [giza|pgplot] [debug|release]
#
# Options (order-independent, defaults: giza, release):
#   giza    — build Giza instead of pre-built PGPLOT (default)
#   pgplot  — link pre-built libpgplot.a + X11
#   debug   — -O0 -check all
#   release — -O2 (default)
#
# Examples:
#   source ./makejoint.sh                      # giza, release
#   source ./makejoint.sh debug                # giza, debug
#   source ./makejoint.sh pgplot               # pgplot, release
#   source ./makejoint.sh pgplot debug         # pgplot, debug
#
# Output locations after a default (giza, release) build:
#   macos executable  : ./build_release_giza/bin/macos
#   smacos library    : ./build_release_giza/lib/libsmacos.a
#   smacos_dvr        : ./build_release_giza/bin/smacos_dvr
#   GMI mex           : ./build_release_giza/lib/GMI.mexa64
#
# For other combinations substitute the build directory name:
#   build_release_pgplot, build_debug_giza, build_debug_pgplot

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
      echo "makejoint.sh: unknown option '$arg'"
      echo "Usage: source ./makejoint.sh [giza|pgplot] [debug|release]"
      return 1 2>/dev/null || exit 1
      ;;
  esac
done

BUILD_DIR="${SCRIPT_DIR}/build_$(echo "${BUILD_TYPE}" | tr '[:upper:]' '[:lower:]')_${PLOT_CHOICE}"
echo "makejoint: BUILD_TYPE=${BUILD_TYPE}  PLOT=${PLOT_CHOICE}"
echo "makejoint: BUILD_DIR=${BUILD_DIR}"

# --- Intel oneAPI environment ---
ONEAPI_SETVARS="/opt/intel/oneapi/setvars.sh"
if [ -f "$ONEAPI_SETVARS" ]; then
  source "$ONEAPI_SETVARS" --force > /dev/null 2>&1
else
  echo "makejoint: WARNING — Intel oneAPI setvars.sh not found; using system compilers"
fi

# --- LD_LIBRARY_PATH ---
READLINE_LIB="${SCRIPT_DIR}/macos_f90/readline-8.2/shlib"
INTEL_LIB="/opt/intel/oneapi/compiler/latest/lib"
if [ "$PLOT_CHOICE" = "pgplot" ]; then
  PGPLOT_DIR="${SCRIPT_DIR}/macos_f90/pgplot"
  export LD_LIBRARY_PATH="${READLINE_LIB}:${INTEL_LIB}:${PGPLOT_DIR}:${LD_LIBRARY_PATH}"
else
  export LD_LIBRARY_PATH="${READLINE_LIB}:${INTEL_LIB}:${LD_LIBRARY_PATH}"
fi

# --- USE_GIZA flag ---
if [ "$PLOT_CHOICE" = "giza" ]; then
  USE_GIZA_FLAG="ON"
else
  USE_GIZA_FLAG="OFF"
fi

# --- Configure ---
cmake -S "${SCRIPT_DIR}" -B "${BUILD_DIR}" \
  -DCMAKE_Fortran_COMPILER=ifx \
  -DCMAKE_C_COMPILER=gcc \
  -DCMAKE_BUILD_TYPE="${BUILD_TYPE}" \
  -DUSE_GIZA="${USE_GIZA_FLAG}" \
  -DBUILD_MACOS=ON \
  -DBUILD_SMACOS=ON \
  -DBUILD_SMACOS_DVR=ON \
  || { echo "makejoint: cmake configure FAILED"; return 1 2>/dev/null || exit 1; }

# --- Build macos / smacos / smacos_dvr ---
cmake --build "${BUILD_DIR}" -j"$(nproc)" \
  || { echo "makejoint: cmake build FAILED"; return 1 2>/dev/null || exit 1; }

# --- Build GMI via its own Makefile (handles MATLAB detection itself) ---
GMI_DIR="${SCRIPT_DIR}/../MACOS_resources/GMI"
GMI_MEX="${GMI_DIR}/GMI.mexa64"
if [ -f "${GMI_DIR}/Makefile" ]; then
  echo "makejoint: building GMI via ${GMI_DIR}/Makefile ..."
  make -C "${GMI_DIR}" macossrc_dir="${SCRIPT_DIR}/macos_f90" \
    || { echo "makejoint: WARNING — GMI build failed (continuing)"; }
else
  echo "makejoint: GMI Makefile not found at ${GMI_DIR}, skipping GMI build."
fi

# --- Report ---
echo ""
echo "makejoint: BUILD SUCCEEDED"
echo "  macos      : ${BUILD_DIR}/bin/macos"
echo "  smacos     : ${BUILD_DIR}/lib/libsmacos.a"
if [ -f "${BUILD_DIR}/bin/smacos_dvr" ]; then
  echo "  smacos_dvr : ${BUILD_DIR}/bin/smacos_dvr"
fi
if [ -f "${GMI_MEX}" ]; then
  echo "  GMI        : ${GMI_MEX}"
fi
