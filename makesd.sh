#!/usr/bin/env bash
# makesd.sh  — CMake build for smacos_dvr (reuses existing macos/smacos build)
#
# Usage:
#   source ./makesd.sh [debug|release] [gfortran]
#
# Options (order-independent, defaults: release, ifx):
#   debug    — -O0 -check all
#   release  — -O2 (default)
#   gfortran — use gfortran instead of ifx
#
# Requires macos and smacos to be already built in the same build directory
# (run makems.sh or makeall.sh first).
#
# Output location (default build):
#   smacos_dvr : ./build_release/bin/smacos_dvr
#
# Note: PGPLOT support was removed on release-candidate; Giza is the
# only PGPLOT-API provider now.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# --- Parse arguments ---
BUILD_TYPE="Release"
FC="ifx"

for arg in "$@"; do
  case "$arg" in
    debug)    BUILD_TYPE="Debug"    ;;
    release)  BUILD_TYPE="Release"  ;;
    ifx)      FC="ifx"              ;;
    gfortran) FC="gfortran"         ;;
    *)
      echo "makesd.sh: unknown option '$arg'"
      echo "Usage: source ./makesd.sh [debug|release] [gfortran]"
      return 1 2>/dev/null || exit 1
      ;;
  esac
done

BUILD_TAG="$(echo "${BUILD_TYPE}" | tr '[:upper:]' '[:lower:]')"
[ "${FC}" = "gfortran" ] && BUILD_TAG="${BUILD_TAG}_gfortran"
BUILD_DIR="${SCRIPT_DIR}/build_${BUILD_TAG}"

echo "makesd: FC=${FC}  BUILD_TYPE=${BUILD_TYPE}"
echo "makesd: BUILD_DIR=${BUILD_DIR}"

# --- Intel oneAPI environment (ifx builds) ---
ONEAPI_SETVARS="/opt/intel/oneapi/setvars.sh"
if [ "${FC}" = "ifx" ] && [ -f "${ONEAPI_SETVARS}" ]; then
  source "${ONEAPI_SETVARS}" --force > /dev/null 2>&1
fi

# --- LD_LIBRARY_PATH ---
READLINE_LIB="${SCRIPT_DIR}/macos_f90/readline-8.2/shlib"
INTEL_LIB="/opt/intel/oneapi/compiler/latest/lib"
export LD_LIBRARY_PATH="${READLINE_LIB}:${INTEL_LIB}:${LD_LIBRARY_PATH}"

# --- Generator: prefer Ninja if installed, else Unix Makefiles.  On a
# --- re-configure, reuse the existing cache's generator (avoids a mismatch). ---
if [ -f "${BUILD_DIR}/CMakeCache.txt" ]; then
  GEN_ARG=()
elif command -v ninja >/dev/null 2>&1 || command -v ninja-build >/dev/null 2>&1; then
  GEN_ARG=(-G Ninja)
else
  GEN_ARG=(-G "Unix Makefiles")
fi

# --- Configure (adds smacos_dvr to existing build; macos+smacos are no-ops if current) ---
cmake -S "${SCRIPT_DIR}" -B "${BUILD_DIR}" "${GEN_ARG[@]}" \
  -DCMAKE_Fortran_COMPILER="${FC}" \
  -DCMAKE_C_COMPILER=gcc \
  -DCMAKE_BUILD_TYPE="${BUILD_TYPE}" \
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
