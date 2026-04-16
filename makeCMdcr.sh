# CMake build for macos + smacos + smacos_dvr
#
# Usage: source ./makeCMdcr.sh [debug] [gfortran]
#   debug    — build with -O0 and runtime checks
#   gfortran — use gfortran instead of ifx (default)
#   Both options can be combined in any order.

# ---- Parse options ----
BUILD_TYPE="Release"
FC_CHOICE="ifx"
for arg in "$@"; do
    case "$arg" in
        debug)    BUILD_TYPE="Debug" ;;
        gfortran) FC_CHOICE="gfortran" ;;
        *)        echo "Unknown option: $arg"; echo "Usage: source ./makeCMdcr.sh [debug] [gfortran]"; return 1 2>/dev/null || exit 1 ;;
    esac
done

# ---- Intel environment (needed even for gfortran link against Intel libs) ----
source /opt/intel/oneapi/setvars.sh intel64 --force 2>/dev/null

# ---- Paths ----
cd ~/dev/macos/macos_f90
export macossrc_dir=$(pwd)
export PGPLOT_FONT=$macossrc_dir/grfont.dat
export PGPLOT_DIR=$macossrc_dir/pgplot/:"$PGPLOT_DIR"
export LD_LIBRARY_PATH=$macossrc_dir/readline-8.2:/opt/intel/oneapi/compiler/latest/lib:/usr/local/pgplot:$LD_LIBRARY_PATH
export intel64_lib="/opt/intel/oneapi/compiler/latest/lib"

# ---- Configure + Build ----
BUILD_DIR="build_$(echo $BUILD_TYPE | tr A-Z a-z)_${FC_CHOICE}"

echo "*** CMake build: $BUILD_TYPE with $FC_CHOICE → $BUILD_DIR"

cmake -S . -B "$BUILD_DIR" \
    -DCMAKE_Fortran_COMPILER="$FC_CHOICE" \
    -DCMAKE_C_COMPILER=gcc \
    -DCMAKE_BUILD_TYPE="$BUILD_TYPE" \
    -DBUILD_SMACOS_DVR=ON

cmake --build "$BUILD_DIR" -j$(nproc)

if [ $? -eq 0 ]; then
    echo "*** Build succeeded: $BUILD_DIR/macos, $BUILD_DIR/smacos_dvr, $BUILD_DIR/libsmacos_lib.a"
    # Symlink for convenience
    ln -sf "$BUILD_DIR/macos" macos
    ln -sf "$BUILD_DIR/smacos_dvr" smacos_dvr
else
    echo "*** Build FAILED"
fi

cd ../
