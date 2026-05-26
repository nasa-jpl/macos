#!/usr/bin/env bash
# makegfortran.sh  — CMake build for macos + smacos + smacos_dvr + GMI using gfortran/gcc
#
# Usage:
#   source ./makegfortran.sh [debug|release]
#
# Options (order-independent, default: release):
#   debug   — -O0 -fcheck=all
#   release — -O2 (default)
#
# Examples:
#   source ./makegfortran.sh                      # release
#   source ./makegfortran.sh debug                # debug
#
# Output locations after a default release build:
#   macos executable  : ./build_release_gfortran/bin/macos
#   smacos library    : ./build_release_gfortran/lib/libsmacos.a
#   smacos_dvr        : ./build_release_gfortran/bin/smacos_dvr
#   GMI mex           : ./build_release_gfortran/lib/GMI.mexa64
#
# Note: PGPLOT support was removed on release-candidate; Giza is now
# the only PGPLOT-API provider.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# --- Parse arguments ---
BUILD_TYPE="Release"
USE_NPSOL="OFF"

for arg in "$@"; do
  case "$arg" in
    debug)   BUILD_TYPE="Debug"    ;;
    release) BUILD_TYPE="Release"  ;;
    npsol)   USE_NPSOL="ON"        ;;
    *)
      echo "makegfortran.sh: unknown option '$arg'"
      echo "Usage: source ./makegfortran.sh [debug|release] [npsol]"
      return 1 2>/dev/null || exit 1
      ;;
  esac
done

BUILD_TAG="$(echo "${BUILD_TYPE}" | tr '[:upper:]' '[:lower:]')_gfortran"
[ "${USE_NPSOL}" = "ON" ] && BUILD_TAG="${BUILD_TAG}_npsol"
BUILD_DIR="${SCRIPT_DIR}/build_${BUILD_TAG}"
echo "makegfortran: BUILD_TYPE=${BUILD_TYPE}  USE_NPSOL=${USE_NPSOL}"
echo "makegfortran: BUILD_DIR=${BUILD_DIR}"

# --- LD_LIBRARY_PATH ---
READLINE_LIB="${SCRIPT_DIR}/macos_f90/readline-8.2/shlib"
export LD_LIBRARY_PATH="${READLINE_LIB}:${LD_LIBRARY_PATH}"

# --- Configure ---
cmake -S "${SCRIPT_DIR}" -B "${BUILD_DIR}" \
  -DCMAKE_Fortran_COMPILER=gfortran \
  -DCMAKE_C_COMPILER=gcc \
  -DCMAKE_BUILD_TYPE="${BUILD_TYPE}" \
  -DBUILD_MACOS=ON \
  -DBUILD_SMACOS=ON \
  -DBUILD_SMACOS_DVR=ON \
  -DBUILD_GMI=ON \
  -DUSE_NPSOL="${USE_NPSOL}" \
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
