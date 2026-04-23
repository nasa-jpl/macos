#!/usr/bin/env bash
# makegmi.sh  — Build GMI.mexa64 via MACOS_resources/GMI/Makefile
#
# Usage:
#   source ./makegmi.sh [giza|pgplot] [debug|release] [gfortran]
#
# Options (order-independent, defaults: giza, release, ifx):
#   giza     — use build_release_giza module files (default)
#   pgplot   — use build_release_pgplot module files
#   debug    — use debug build module files
#   release  — use release build module files (default)
#   gfortran — use build_*_gfortran module files
#              (Note: GMI itself is always compiled with ifx)
#
# Requires macos+smacos to be built first in the matching build directory.
# MATLAB is auto-detected under /usr/local/MATLAB/.
#
# Output location:
#   GMI.mexa64 : ~/dev/MACOS_resources/GMI/GMI.mexa64

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
      echo "makegmi.sh: unknown option '$arg'"
      echo "Usage: source ./makegmi.sh [giza|pgplot] [debug|release] [gfortran]"
      return 1 2>/dev/null || exit 1
      ;;
  esac
done

BUILD_TAG="$(echo "${BUILD_TYPE}" | tr '[:upper:]' '[:lower:]')_${PLOT_CHOICE}"
[ "${FC}" = "gfortran" ] && BUILD_TAG="${BUILD_TAG}_gfortran"
BUILD_DIR="${SCRIPT_DIR}/build_${BUILD_TAG}"

echo "makegmi: PLOT=${PLOT_CHOICE}  BUILD_TYPE=${BUILD_TYPE}  FC=${FC}"
echo "makegmi: MACOS_BUILD_DIR=${BUILD_DIR}"

# --- Intel oneAPI environment (GMI Makefile uses ifx) ---
ONEAPI_SETVARS="/opt/intel/oneapi/setvars.sh"
if [ -f "${ONEAPI_SETVARS}" ]; then
  source "${ONEAPI_SETVARS}" --force > /dev/null 2>&1
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
  echo "  Build macos+smacos first: source ./makems.sh${FC:+ ${FC}} ${PLOT_CHOICE} ${BUILD_TYPE,,}"
  return 1 2>/dev/null || exit 1
fi

# --- Build ---
echo "makegmi: building GMI via ${GMI_DIR}/Makefile ..."
make -C "${GMI_DIR}" \
  macossrc_dir="${SCRIPT_DIR}/macos_f90" \
  MACOS_BUILD_DIR="${BUILD_DIR}" \
  || { echo "makegmi: GMI build FAILED"; return 1 2>/dev/null || exit 1; }

echo ""
echo "makegmi: BUILD SUCCEEDED"
echo "  GMI : ${GMI_DIR}/GMI.mexa64"
