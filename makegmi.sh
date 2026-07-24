#!/usr/bin/env bash
# makegmi.sh  — Build GMI.mexa64 via MACOS_resources/GMI/Makefile
#
# Usage:
#   source ./makegmi.sh [debug|release] [ifx|gfortran]
#
# Options (order-independent, defaults: release, gfortran):
#   debug    — use debug build module files
#   release  — use release build module files (default)
#   gfortran — compile GMI mex with gfortran (default).  Produces a
#              mex that exits MATLAB cleanly.
#   ifx      — compile GMI mex with ifx.  Mex SIGSEGVs at MATLAB
#              process-exit teardown (cosmetic; results correct).
#
# Either compiler choice picks up libsmacos.a from the matching
# build_release[_gfortran] tree, so build macos+smacos with the
# same compiler first (makems.sh for ifx, makegfortran.sh for
# gfortran).  MATLAB is auto-detected under /usr/local/MATLAB/.
#
# Output location:
#   GMI.mexa64 : ~/dev/MACOS_resources/GMI/GMI.mexa64
#
# Note: PGPLOT support was removed on release-candidate; Giza is now
# the only PGPLOT-API provider.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# --- Parse arguments ---
BUILD_TYPE="Release"
FC="gfortran"
USE_NPSOL="OFF"

for arg in "$@"; do
  case "$arg" in
    debug)    BUILD_TYPE="Debug"    ;;
    release)  BUILD_TYPE="Release"  ;;
    ifx)      FC="ifx"              ;;
    gfortran) FC="gfortran"         ;;
    npsol)    USE_NPSOL="ON"        ;;
    *)
      echo "makegmi.sh: unknown option '$arg'"
      echo "Usage: source ./makegmi.sh [debug|release] [ifx|gfortran] [npsol]"
      return 1 2>/dev/null || exit 1
      ;;
  esac
done

BUILD_TAG="$(echo "${BUILD_TYPE}" | tr '[:upper:]' '[:lower:]')"
[ "${FC}" = "gfortran" ]    && BUILD_TAG="${BUILD_TAG}_gfortran"
[ "${USE_NPSOL}" = "ON" ]   && BUILD_TAG="${BUILD_TAG}_npsol"
BUILD_DIR="${SCRIPT_DIR}/build_${BUILD_TAG}"

echo "makegmi: BUILD_TYPE=${BUILD_TYPE}  FC=${FC}  USE_NPSOL=${USE_NPSOL}"
echo "makegmi: MACOS_BUILD_DIR=${BUILD_DIR}"

# --- Intel oneAPI environment (only needed when FC=ifx) ---
if [ "${FC}" = "ifx" ]; then
  ONEAPI_SETVARS="/opt/intel/oneapi/setvars.sh"
  if [ -f "${ONEAPI_SETVARS}" ]; then
    source "${ONEAPI_SETVARS}" --force > /dev/null 2>&1
  fi
fi

# --- Locate MACOS_resources/GMI ---
GMI_DIR="$(realpath "${SCRIPT_DIR}/../MACOS_resources/GMI" 2>/dev/null || \
           (cd "${SCRIPT_DIR}/../MACOS_resources/GMI" && pwd))"

if [ ! -f "${GMI_DIR}/Makefile" ]; then
  echo "makegmi: ERROR — GMI Makefile not found at ${GMI_DIR}"
  echo "  Clone MACOS_resources alongside macos/: git clone git@github.com:nasa-jpl/MACOS_resources.git"
  return 1 2>/dev/null || exit 1
fi

# --- Check that the cmake build dir exists ---
if [ ! -d "${BUILD_DIR}/mod_smacos" ]; then
  echo "makegmi: ERROR — module directory not found: ${BUILD_DIR}/mod_smacos"
  if [ "${FC}" = "gfortran" ]; then
    echo "  Build macos+smacos first: source ./makegfortran.sh ${BUILD_TYPE,,}"
  else
    echo "  Build macos+smacos first: source ./makems.sh ${BUILD_TYPE,,}"
  fi
  return 1 2>/dev/null || exit 1
fi

# --- Build ---
echo "makegmi: building GMI via ${GMI_DIR}/Makefile ..."
make -C "${GMI_DIR}" \
  FC="${FC}" \
  macossrc_dir="${SCRIPT_DIR}/macos_f90" \
  MACOS_BUILD_DIR="${BUILD_DIR}" \
  || { echo "makegmi: GMI build FAILED"; return 1 2>/dev/null || exit 1; }

echo ""
echo "makegmi: BUILD SUCCEEDED"
echo "  GMI : ${GMI_DIR}/GMI.mexa64"
