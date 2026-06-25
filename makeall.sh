#!/usr/bin/env bash
# makeall.sh  — CMake build for macos + smacos + smacos_dvr + GMI
#
# Usage:
#   source ./makeall.sh [debug|release]
#
# Options (default: release):
#   debug   — -O0 -check all
#   release — -O2 (default)
#
# Examples:
#   source ./makeall.sh                      # release
#   source ./makeall.sh debug                # debug
#
# Output locations after a default (release) build:
#   macos executable  : ./build_release/bin/macos
#   smacos library    : ./build_release/lib/libsmacos.a
#   smacos_dvr        : ./build_release/bin/smacos_dvr
#   GMI mex           : ./build_release/lib/GMI.mexa64
#
# Note: PGPLOT support was removed on release-candidate.  Giza is the
# only PGPLOT-API provider now; the [giza|pgplot] arg is gone.

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
      echo "makeall.sh: unknown option '$arg'"
      echo "Usage: source ./makeall.sh [debug|release] [npsol]"
      return 1 2>/dev/null || exit 1
      ;;
  esac
done

BUILD_TAG="$(echo "${BUILD_TYPE}" | tr '[:upper:]' '[:lower:]')"
[ "${USE_NPSOL}" = "ON" ] && BUILD_TAG="${BUILD_TAG}_npsol"
BUILD_DIR="${SCRIPT_DIR}/build_${BUILD_TAG}"
echo "makeall: BUILD_TYPE=${BUILD_TYPE}  USE_NPSOL=${USE_NPSOL}"
echo "makeall: BUILD_DIR=${BUILD_DIR}"

# --- Intel oneAPI environment ---
ONEAPI_SETVARS="/opt/intel/oneapi/setvars.sh"
if [ -f "$ONEAPI_SETVARS" ]; then
  source "$ONEAPI_SETVARS" --force > /dev/null 2>&1
else
  echo "makeall: WARNING — Intel oneAPI setvars.sh not found; using system compilers"
fi

# --- Readline bootstrap (one-time, on first clone) ---
READLINE_DIR="${SCRIPT_DIR}/macos_f90/readline-8.2"
if [ ! -f "${READLINE_DIR}/libreadline.a" ]; then
  echo "makeall: bootstrapping bundled readline (one-time setup) ..."
  ( cd "${READLINE_DIR}" && ./configure -q && make -j"$(nproc)" ) \
    > /tmp/makeall_readline_bootstrap.log 2>&1 \
    || { echo "makeall: readline bootstrap FAILED — see /tmp/makeall_readline_bootstrap.log"
         return 1 2>/dev/null || exit 1; }
fi

# --- LD_LIBRARY_PATH (readline + Intel runtime; pgplot path removed) ---
READLINE_LIB="${READLINE_DIR}/shlib"
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

# --- Configure ---
cmake -S "${SCRIPT_DIR}" -B "${BUILD_DIR}" "${GEN_ARG[@]}" \
  -DCMAKE_Fortran_COMPILER=ifx \
  -DCMAKE_C_COMPILER=gcc \
  -DCMAKE_BUILD_TYPE="${BUILD_TYPE}" \
  -DBUILD_MACOS=ON \
  -DBUILD_SMACOS=ON \
  -DBUILD_SMACOS_DVR=ON \
  -DBUILD_GMI=ON \
  -DUSE_NPSOL="${USE_NPSOL}" \
  || { echo "makeall: cmake configure FAILED"; return 1 2>/dev/null || exit 1; }

# --- Build macos / smacos / smacos_dvr ---
cmake --build "${BUILD_DIR}" -j"$(nproc)" \
  || { echo "makeall: cmake build FAILED"; return 1 2>/dev/null || exit 1; }

# --- Build GMI via its own Makefile (handles MATLAB detection itself) ---
GMI_DIR="$(realpath "${SCRIPT_DIR}/../MACOS_resources/GMI" 2>/dev/null || (cd "${SCRIPT_DIR}/../MACOS_resources/GMI" && pwd))"
GMI_MEX="${GMI_DIR}/GMI.mexa64"
if [ -f "${GMI_DIR}/Makefile" ]; then
  echo "makeall: building GMI via ${GMI_DIR}/Makefile (FC=ifx) ..."
  # Force FC=ifx — makeall.sh is the ifx all-four path.  The GMI
  # Makefile defaults to gfortran; use makegfortran.sh for the
  # all-four gfortran build.
  make -C "${GMI_DIR}" FC=ifx \
    macossrc_dir="${SCRIPT_DIR}/macos_f90" \
    MACOS_BUILD_DIR="${BUILD_DIR}" \
    || { echo "makeall: WARNING — GMI build failed (continuing)"; }
else
  echo "makeall: GMI Makefile not found at ${GMI_DIR}, skipping GMI build."
fi

# --- Report ---
echo ""
echo "makeall: BUILD SUCCEEDED"
echo "  macos      : ${BUILD_DIR}/bin/macos"
echo "  smacos     : ${BUILD_DIR}/lib/libsmacos.a"
if [ -f "${BUILD_DIR}/bin/smacos_dvr" ]; then
  echo "  smacos_dvr : ${BUILD_DIR}/bin/smacos_dvr"
fi
if [ -f "${GMI_MEX}" ]; then
  echo "  GMI        : ${GMI_MEX}"
fi
