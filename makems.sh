#!/usr/bin/env bash
# makems.sh  — CMake build for macos + smacos (no smacos_dvr, no GMI)
#
# Usage:
#   source ./makems.sh [debug|release] [gfortran] [npsol]
#
# Options (order-independent, defaults: release, ifx, no-npsol):
#   debug    — -O0 -check all
#   release  — -O2 (default)
#   gfortran — use gfortran instead of ifx
#   npsol    — build with -DUSE_NPSOL=ON (constrained-optim path);
#              requires macos_f90/npsol/ tree (only on opt-dev branch).
#              Build dir gets _npsol suffix.
#
# Output locations (default build):
#   macos executable  : ./build_release/bin/macos
#   smacos library    : ./build_release/lib/libsmacos.a
#
# Note: PGPLOT support was removed on release-candidate.  Giza is the
# only PGPLOT-API provider; the [giza|pgplot] arg is gone.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# --- Parse arguments ---
BUILD_TYPE="Release"
FC="ifx"
USE_NPSOL="OFF"

for arg in "$@"; do
  case "$arg" in
    debug)    BUILD_TYPE="Debug"    ;;
    release)  BUILD_TYPE="Release"  ;;
    ifx)      FC="ifx"              ;;
    gfortran) FC="gfortran"         ;;
    npsol)    USE_NPSOL="ON"        ;;
    *)
      echo "makems.sh: unknown option '$arg'"
      echo "Usage: source ./makems.sh [debug|release] [gfortran] [npsol]"
      return 1 2>/dev/null || exit 1
      ;;
  esac
done

BUILD_TAG="$(echo "${BUILD_TYPE}" | tr '[:upper:]' '[:lower:]')"
[ "${FC}" = "gfortran" ]    && BUILD_TAG="${BUILD_TAG}_gfortran"
[ "${USE_NPSOL}" = "ON" ]   && BUILD_TAG="${BUILD_TAG}_npsol"
BUILD_DIR="${SCRIPT_DIR}/build_${BUILD_TAG}"

echo "makems: FC=${FC}  BUILD_TYPE=${BUILD_TYPE}  USE_NPSOL=${USE_NPSOL}"
echo "makems: BUILD_DIR=${BUILD_DIR}"

# --- Intel oneAPI environment (ifx builds) ---
ONEAPI_SETVARS="/opt/intel/oneapi/setvars.sh"
if [ "${FC}" = "ifx" ] && [ -f "${ONEAPI_SETVARS}" ]; then
  source "${ONEAPI_SETVARS}" --force > /dev/null 2>&1
fi

# --- Readline bootstrap (one-time, on first clone) ---
READLINE_DIR="${SCRIPT_DIR}/macos_f90/readline-8.2"
if [ ! -f "${READLINE_DIR}/libreadline.a" ]; then
  echo "makems: bootstrapping bundled readline (one-time setup) ..."
  ( cd "${READLINE_DIR}" && ./configure -q && make -j"$(nproc)" ) \
    > /tmp/makems_readline_bootstrap.log 2>&1 \
    || { echo "makems: readline bootstrap FAILED — see /tmp/makems_readline_bootstrap.log"
         return 1 2>/dev/null || exit 1; }
fi

# --- LD_LIBRARY_PATH ---
READLINE_LIB="${READLINE_DIR}/shlib"
INTEL_LIB="/opt/intel/oneapi/compiler/latest/lib"
export LD_LIBRARY_PATH="${READLINE_LIB}:${INTEL_LIB}:${LD_LIBRARY_PATH}"

# --- Configure ---
cmake -S "${SCRIPT_DIR}" -B "${BUILD_DIR}" \
  -DCMAKE_Fortran_COMPILER="${FC}" \
  -DCMAKE_C_COMPILER=gcc \
  -DCMAKE_BUILD_TYPE="${BUILD_TYPE}" \
  -DBUILD_MACOS=ON \
  -DBUILD_SMACOS=ON \
  -DUSE_NPSOL="${USE_NPSOL}" \
  || { echo "makems: cmake configure FAILED"; return 1 2>/dev/null || exit 1; }

# --- Build ---
cmake --build "${BUILD_DIR}" -j"$(nproc)" \
  || { echo "makems: cmake build FAILED"; return 1 2>/dev/null || exit 1; }

echo ""
echo "makems: BUILD SUCCEEDED"
echo "  macos  : ${BUILD_DIR}/bin/macos"
echo "  smacos : ${BUILD_DIR}/lib/libsmacos.a"
