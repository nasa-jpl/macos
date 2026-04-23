#!/usr/bin/env bash
# makesd.sh  — CMake build for smacos_dvr (reuses existing macos/smacos build)
#
# Usage:
#   source ./makesd.sh [giza|pgplot] [debug|release] [gfortran]
#
# Options (order-independent, defaults: giza, release, ifx):
#   giza     — Giza graphics (default)
#   pgplot   — pre-built libpgplot.a + X11
#   debug    — -O0 -check all
#   release  — -O2 (default)
#   gfortran — use gfortran instead of ifx
#
# Requires macos and smacos to be already built in the same build directory
# (run makems.sh or makejoint.sh first).
#
# Output location (default build):
#   smacos_dvr : ./build_release_giza/bin/smacos_dvr

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# --- Parse arguments ---
PLOT_CHOICE="giza"
BUILD_TYPE="Release"
FC="ifx"

for arg in "$@"; do
  case "$arg" in
    giza)     PLOT_CHOICE="giza"    ;;
    pgplot)   PLOT_CHOICE="pgplot"  ;;
    debug)    BUILD_TYPE="Debug"    ;;
    release)  BUILD_TYPE="Release"  ;;
    ifx)      FC="ifx"              ;;
    gfortran) FC="gfortran"         ;;
    *)
      echo "makesd.sh: unknown option '$arg'"
      echo "Usage: source ./makesd.sh [giza|pgplot] [debug|release] [gfortran]"
      return 1 2>/dev/null || exit 1
      ;;
  esac
done

BUILD_TAG="$(echo "${BUILD_TYPE}" | tr '[:upper:]' '[:lower:]')_${PLOT_CHOICE}"
[ "${FC}" = "gfortran" ] && BUILD_TAG="${BUILD_TAG}_gfortran"
BUILD_DIR="${SCRIPT_DIR}/build_${BUILD_TAG}"

echo "makesd: FC=${FC}  BUILD_TYPE=${BUILD_TYPE}  PLOT=${PLOT_CHOICE}"
echo "makesd: BUILD_DIR=${BUILD_DIR}"

# --- Intel oneAPI environment (ifx builds) ---
ONEAPI_SETVARS="/opt/intel/oneapi/setvars.sh"
if [ "${FC}" = "ifx" ] && [ -f "${ONEAPI_SETVARS}" ]; then
  source "${ONEAPI_SETVARS}" --force > /dev/null 2>&1
fi

# --- LD_LIBRARY_PATH ---
READLINE_LIB="${SCRIPT_DIR}/macos_f90/readline-8.2/shlib"
INTEL_LIB="/opt/intel/oneapi/compiler/latest/lib"
if [ "${PLOT_CHOICE}" = "pgplot" ]; then
  export LD_LIBRARY_PATH="${READLINE_LIB}:${INTEL_LIB}:${SCRIPT_DIR}/macos_f90/pgplot:${LD_LIBRARY_PATH}"
else
  export LD_LIBRARY_PATH="${READLINE_LIB}:${INTEL_LIB}:${LD_LIBRARY_PATH}"
fi

# --- Configure (adds smacos_dvr to existing build; macos+smacos are no-ops if current) ---
cmake -S "${SCRIPT_DIR}" -B "${BUILD_DIR}" \
  -DCMAKE_Fortran_COMPILER="${FC}" \
  -DCMAKE_C_COMPILER=gcc \
  -DCMAKE_BUILD_TYPE="${BUILD_TYPE}" \
  -DUSE_GIZA="$([ "${PLOT_CHOICE}" = "giza" ] && echo ON || echo OFF)" \
  -DBUILD_MACOS=ON \
  -DBUILD_SMACOS=ON \
  -DBUILD_SMACOS_DVR=ON \
  || { echo "makesd: cmake configure FAILED"; return 1 2>/dev/null || exit 1; }

# --- Build ---
cmake --build "${BUILD_DIR}" -j"$(nproc)" \
  || { echo "makesd: cmake build FAILED"; return 1 2>/dev/null || exit 1; }

echo ""
echo "makesd: BUILD SUCCEEDED"
echo "  smacos_dvr : ${BUILD_DIR}/bin/smacos_dvr"
