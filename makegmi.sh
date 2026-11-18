#!/usr/bin/env bash
# makegmi.sh  — Build the GMI mex (GMI.mexa64 on Linux, GMI.mexmaca64 on an Apple Silicon Mac)
#               via app/GMI/Makefile
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
# gfortran).  On macOS only gfortran is supported (ifx has no macOS build).
#
# MATLAB is auto-detected: the latest release under /usr/local/MATLAB/ on
# Linux, /Applications/MATLAB_R20XXx.app on macOS.  To choose another one,
# set MATLAB_DIR in the environment before sourcing this script:
#   MATLAB_DIR=/Applications/MATLAB_R2025b.app source ./makegmi.sh
#
# Output location:
#   ./app/GMI/GMI.<mex extension of this platform>
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

# lower case without ${VAR,,}: that form needs bash 4, and macOS ships bash 3.2 and zsh
BUILD_TYPE_LC="$(echo "${BUILD_TYPE}" | tr '[:upper:]' '[:lower:]')"
BUILD_TAG="${BUILD_TYPE_LC}"
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

# --- Locate app/GMI/ ---
GMI_DIR="${SCRIPT_DIR}/app/GMI"

if [ ! -f "${GMI_DIR}/Makefile" ]; then
  echo "makegmi: ERROR — GMI Makefile not found at ${GMI_DIR}"
  return 1 2>/dev/null || exit 1
fi

# --- Check that the cmake build dir exists ---
if [ ! -d "${BUILD_DIR}/mod_smacos" ]; then
  echo "makegmi: ERROR — module directory not found: ${BUILD_DIR}/mod_smacos"
  if [ "${FC}" = "gfortran" ]; then
    echo "  Build macos+smacos first: source ./makems.sh ${BUILD_TYPE_LC} gfortran"
  else
    echo "  Build macos+smacos first: source ./makems.sh ${BUILD_TYPE_LC}"
  fi
  return 1 2>/dev/null || exit 1
fi

# --- Build ---
echo "makegmi: building GMI via ${GMI_DIR}/Makefile ..."
# MATLAB_DIR from the environment is passed on when it is set; otherwise the
# Makefile finds MATLAB itself.
GMI_MAKE_ARGS=(FC="${FC}" macossrc_dir="${SCRIPT_DIR}/src" MACOS_BUILD_DIR="${BUILD_DIR}")
if [ -n "${MATLAB_DIR:-}" ]; then
  echo "makegmi: MATLAB_DIR=${MATLAB_DIR}"
  GMI_MAKE_ARGS+=(MATLAB_DIR="${MATLAB_DIR}")
fi
make -C "${GMI_DIR}" "${GMI_MAKE_ARGS[@]}" \
  || { echo "makegmi: GMI build FAILED"; return 1 2>/dev/null || exit 1; }

GMI_MEXTAG="$(make -s -C "${GMI_DIR}" "${GMI_MAKE_ARGS[@]}" print-mextag 2>/dev/null)"
echo ""
echo "makegmi: BUILD SUCCEEDED"
echo "  GMI : ${GMI_DIR}/GMI.${GMI_MEXTAG:-mex}"
