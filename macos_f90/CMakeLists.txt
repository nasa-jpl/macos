#########################################################################
#  CMake build for MACOS / SMACOS / SMACOS_DVR / GMI
#
#  Usage:
#    mkdir build && cd build
#    cmake .. -DCMAKE_Fortran_COMPILER=ifx -DCMAKE_C_COMPILER=gcc
#    make -j$(nproc)
#
#  For gfortran:
#    cmake .. -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_C_COMPILER=gcc
#
#  Options:
#    -DCMAKE_BUILD_TYPE=Debug        # -O0 -g -check all (ifx) or -fbounds-check (gfortran)
#    -DCMAKE_BUILD_TYPE=Release      # -Ofast (default)
#    -DBUILD_GMI=ON                  # build GMI mex target (requires Matlab)
#    -DMATLAB_ROOT=/usr/local/MATLAB/R2025b
#    -DBUILD_SMACOS_DVR=ON          # build smacos_dvr executable
#########################################################################
cmake_minimum_required(VERSION 3.20)
project(MACOS VERSION 3.34 LANGUAGES Fortran C)

# ---- Options ----
option(BUILD_GMI        "Build GMI Matlab mex interface"  OFF)
option(BUILD_SMACOS_DVR "Build smacos_dvr executable"     OFF)

# Auto-detect the latest MATLAB under /usr/local/MATLAB (R20XXa/b).
# Override: -DMATLAB_ROOT=/usr/local/MATLAB/R2023a
if(NOT DEFINED MATLAB_ROOT OR MATLAB_ROOT STREQUAL "")
    file(GLOB _matlab_dirs LIST_DIRECTORIES true "/usr/local/MATLAB/R[0-9][0-9][0-9][0-9][ab]")
    if(_matlab_dirs)
        list(SORT _matlab_dirs COMPARE NATURAL ORDER ASCENDING)
        list(GET  _matlab_dirs -1 MATLAB_ROOT)
    endif()
endif()
set(MATLAB_ROOT "${MATLAB_ROOT}" CACHE PATH "Matlab installation root" FORCE)

# ---- Platform detection ----
if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    if(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64")
        set(MACOS_OS "LNXOS64")
    else()
        set(MACOS_OS "LNXOS32")
    endif()
elseif(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    set(MACOS_OS "Darwin")
endif()
set(MACOS_FFT "NR_FFT")

# ---- Compiler flags ----
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    # ifx / ifort
    set(MACOS_FORTRAN_FLAGS "-xhost -fPIC -132 -traceback")
    set(MACOS_FORTRAN_FLAGS_RELEASE "-Ofast")
    set(MACOS_FORTRAN_FLAGS_DEBUG   "-O0 -g -check all,noarg_temp_created")
    set(MACOS_LINK_FLAGS "-O3 -qopenmp")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    # gfortran
    set(MACOS_FORTRAN_FLAGS "-fPIC -ffixed-line-length-132 -fd-lines-as-comments -fallow-argument-mismatch -std=legacy")
    set(MACOS_FORTRAN_FLAGS_RELEASE "-Ofast")
    set(MACOS_FORTRAN_FLAGS_DEBUG   "-O0 -g -fbounds-check -fcheck=all")
    set(MACOS_LINK_FLAGS "-O3 -fopenmp")
endif()

set(CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} ${MACOS_FORTRAN_FLAGS}")
set(CMAKE_Fortran_FLAGS_RELEASE "${MACOS_FORTRAN_FLAGS_RELEASE}")
set(CMAKE_Fortran_FLAGS_DEBUG   "${MACOS_FORTRAN_FLAGS_DEBUG}")
set(CMAKE_EXE_LINKER_FLAGS      "${CMAKE_EXE_LINKER_FLAGS} ${MACOS_LINK_FLAGS}")

# C compiler flags (legacy C files, match original gcc -fPIC)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")

# Default to Release
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# ---- External pre-built libraries ----
# pgplot: expected at ${CMAKE_CURRENT_SOURCE_DIR}/pgplot/libpgplot.a (or .so)
# fitsio: expected at ${CMAKE_CURRENT_SOURCE_DIR}/fits_build/libfits.a
# These are built by their own build systems before running CMake.
# If not found, the build will fail at link time with a clear message.

find_library(PGPLOT_LIB pgplot
    PATHS ${CMAKE_CURRENT_SOURCE_DIR}/pgplot
    NO_DEFAULT_PATH)
if(NOT PGPLOT_LIB)
    message(WARNING "libpgplot not found in pgplot/ — build pgplot first (cd pgplot && make)")
endif()

find_library(FITS_LIB fits
    PATHS ${CMAKE_CURRENT_SOURCE_DIR}/fits_build
    NO_DEFAULT_PATH)
if(NOT FITS_LIB)
    message(WARNING "libfits not found in fits_build/ — build fitsio first (cd fits_build && bash cmp && bash amp)")
endif()

find_library(X11_LIB X11 PATHS /usr/lib64 /usr/lib/x86_64-linux-gnu /usr/lib)

# Intel runtime libraries (only needed for ifx)
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(INTEL_LIB_DIR "/opt/intel/oneapi/compiler/latest/lib" CACHE PATH
        "Intel compiler runtime library directory")
    find_library(LIBIMF    imf    PATHS ${INTEL_LIB_DIR} NO_DEFAULT_PATH)
    find_library(LIBIFPORT ifport PATHS ${INTEL_LIB_DIR} NO_DEFAULT_PATH)
    find_library(LIBIFCORE ifcore PATHS ${INTEL_LIB_DIR} NO_DEFAULT_PATH)
    set(INTEL_LIBS ${LIBIMF} ${LIBIFPORT} ${LIBIFCORE})
else()
    set(INTEL_LIBS "")
    set(INTEL_LIB_DIR "")
endif()

# ---- Preprocessor defines ----
# Common to all targets
set(COMMON_DEFINES ${MACOS_OS} ${MACOS_FFT} MT0)
# BUILD_LOC is set per-target below

# ---- NPSOL sub-library ----
add_subdirectory(npsol)

# ---- Source file lists ----
# Module sources (order matters for Fortran module dependencies)
set(MOD_SRCS
    kinds.f90
    constants.f90
    dftsub.F
    lohpars_mod.F
    macos_IO.f90
    param_mod.F
    elt_mod.F
    nn_util.F
    mathsub.F
    macos_debug.F
    src_mod.F
    cfiles_mod.F
    usersub.F
    lsq.F
    loh_mod.F
    traceutil_mod.F
    nls.F
    macos_mod.F
    dopt_mod.F
    surfsub.F
    didesub.F
    elemsub.F
)

# C sources for macos
set(C_SRCS
    utilsub_c.c
    mhist.c
    splicubi2.c
)

# C sources for smacos (fk.c instead of mhist.c)
set(C2_SRCS
    utilsub_c.c
    fk.c
    splicubi2.c
)

# Fortran sources specific to macos executable
set(M_SRCS
    macos_vars_mod.F
    macosio.F
    macos_init.F
    funcsub.F
    pixsub.F
    sourcsub.F
    tracesub.F
    srtrace.F
    utilsub.F
    linsub.F
    propsub.F
    pgplotsub.F
    sunsub.F
    macos_ops.F
    zern_wf.F
    design_optim.F
    # design_cons_optim.F removed from build by rem_npsol PR (main):
    # the NPSOL-based constrained-optim path is dead code now; the
    # file stays in the tree but isn't compiled because its npsol
    # callees were deleted.  Levenberg-Marquardt via nls_optim_dvr
    # (inside design_optim_mod) is the replacement.
    nls.F
    stop_set.F
    psesub.F
    f90prt.F
)

# Fortran sources specific to smacos library
set(SM_SRCS
    smacos_vars_mod.F
    smacosutil.F
    smacosio.F
    macos_init.F
    funcsub.F
    pixsub.F
    sourcsub.F
    tracesub.F
    srtrace.F
    utilsub.F
    linsub.F
    propsub.F
    pgplotdummy.F
    sunsub.F
    macos_ops.F
    zern_wf.F
    design_optim.F
    # design_cons_optim.F removed -- see comment above.
    nls.F
    stop_set.F
    psesub.F
)

#########################################################################
#  Target 1: macos (interactive executable)
#########################################################################
add_executable(macos
    macos.F
    ${MOD_SRCS}
    ${M_SRCS}
    ${C_SRCS}
)

target_compile_definitions(macos PRIVATE
    ${COMMON_DEFINES}
    CMACOS
    BUILD_LOC='${CMAKE_CURRENT_SOURCE_DIR}'
)

# Fortran flags: -g for local compile, defines via target_compile_definitions
target_compile_options(macos PRIVATE
    $<$<COMPILE_LANGUAGE:Fortran>:-g>
)

target_link_libraries(macos PRIVATE
    npsol
    ${INTEL_LIBS}
    ${X11_LIB}
    ${PGPLOT_LIB}
    ${FITS_LIB}
)

# Place .mod files in a per-target directory to avoid collisions
set_target_properties(macos PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod_macos
)
target_include_directories(macos PRIVATE
    ${CMAKE_CURRENT_BINARY_DIR}/mod_macos
    ${CMAKE_CURRENT_SOURCE_DIR}
)

#########################################################################
#  Target 2: smacos_lib (static library)
#########################################################################
add_library(smacos_lib STATIC
    smacos.F
    ${MOD_SRCS}
    ${SM_SRCS}
    ${C2_SRCS}
)

target_compile_definitions(smacos_lib PRIVATE
    ${COMMON_DEFINES}
    CSMACOS
    BUILD_LOC='${CMAKE_CURRENT_SOURCE_DIR}'
)

target_compile_options(smacos_lib PRIVATE
    $<$<COMPILE_LANGUAGE:Fortran>:-g>
)

# smacos_lib also includes npsol objects in the final archive
# We handle this by linking npsol into the final consumers instead
set_target_properties(smacos_lib PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod_smacos
    OUTPUT_NAME smacos_lib
)
target_include_directories(smacos_lib PRIVATE
    ${CMAKE_CURRENT_BINARY_DIR}/mod_smacos
    ${CMAKE_CURRENT_SOURCE_DIR}
)

#########################################################################
#  Target 3: smacos_dvr (driver executable, optional)
#########################################################################
if(BUILD_SMACOS_DVR)
    # smacos_dvr is compiled with -DCMACOS and needs the macos-flavored
    # .mod files (e.g. macos_mod has ifPGColor only under CMACOS).
    # It links against smacos_lib plus macosio, pgplotsub, macos_vars_mod
    # from the macos build.  We use the macos target's module directory.
    # These sources are re-compiled with -DCMACOS (the macos_mod in
    # smacos_lib.a was compiled with -DCSMACOS and lacks CMACOS symbols).
    # The linker resolves these objects before the archive, matching the
    # original Makefile which listed them explicitly before smacos_lib.a.
    add_executable(smacos_dvr
        smacos_dvr.F
        macos_mod.F
        macos_vars_mod.F
        macosio.F
        pgplotsub.F
    )

    target_compile_definitions(smacos_dvr PRIVATE
        ${COMMON_DEFINES}
        CMACOS
        BUILD_LOC='${CMAKE_CURRENT_SOURCE_DIR}'
    )

    target_compile_options(smacos_dvr PRIVATE -g)

    # Use macos target's .mod files (CMACOS flavor), not smacos (CSMACOS)
    set_target_properties(smacos_dvr PROPERTIES
        Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod_macos
    )
    target_include_directories(smacos_dvr PRIVATE
        ${CMAKE_CURRENT_BINARY_DIR}/mod_macos
        ${CMAKE_CURRENT_BINARY_DIR}/mod_smacos
        ${CMAKE_CURRENT_SOURCE_DIR}
    )
    # Must build macos first (CMACOS .mod files) and smacos_lib (smacos_mod)
    add_dependencies(smacos_dvr macos smacos_lib)

    target_link_libraries(smacos_dvr PRIVATE
        smacos_lib
        npsol
        ${INTEL_LIBS}
        ${X11_LIB}
        ${PGPLOT_LIB}
    )
endif()

#########################################################################
#  Target 4: GMI (Matlab mex shared library, optional)
#########################################################################
if(BUILD_GMI)
    set(GMI_SRC_DIR "$ENV{HOME}/dev/MACOS_resources/GMI" CACHE PATH
        "Directory containing GMI.F and GMIG.F")

    if(NOT EXISTS "${GMI_SRC_DIR}/GMI.F")
        message(WARNING
            "BUILD_GMI=ON but GMI.F not found in ${GMI_SRC_DIR}. "
            "Set -DGMI_SRC_DIR=<path> to the MACOS_resources/GMI directory. "
            "Disabling GMI target.")
        set(BUILD_GMI OFF)
    elseif(NOT MATLAB_ROOT
           OR NOT EXISTS "${MATLAB_ROOT}/extern/include/mex.h"
           OR NOT EXISTS "${MATLAB_ROOT}/extern/include/fintrf.h")
        message(WARNING
            "BUILD_GMI=ON but MATLAB MEX/Fortran headers not found under "
            "${MATLAB_ROOT}/extern/include (need mex.h and fintrf.h). "
            "Set -DMATLAB_ROOT=/usr/local/MATLAB/R20XXx or install full MATLAB. "
            "Disabling GMI target.")
        set(BUILD_GMI OFF)
    endif()
endif()

if(BUILD_GMI)

    add_library(GMI SHARED
        ${GMI_SRC_DIR}/GMI.F
        ${GMI_SRC_DIR}/GMIG.F
    )

    set_target_properties(GMI PROPERTIES
        OUTPUT_NAME GMI
        SUFFIX ".mexa64"
        PREFIX ""
    )

    target_compile_definitions(GMI PRIVATE
        MX_COMPAT_32
    )

    target_compile_options(GMI PRIVATE
        -fPIC -fpp -132 -gen-interfaces -fp-model strict
        -fno-omit-frame-pointer -O2 -xHOST
        $<$<COMPILE_LANGUAGE:Fortran>:-check all -traceback -fstack-protector>
    )

    target_include_directories(GMI PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_CURRENT_BINARY_DIR}/mod_smacos
        ${MATLAB_ROOT}/extern/include
    )

    target_link_libraries(GMI PRIVATE
        smacos_lib
        npsol
        ${INTEL_LIBS}
        -L${MATLAB_ROOT}/bin/glnxa64
        -l:libmx.so -l:libmex.so -lmat
        -lm -lstdc++
    )

    target_link_options(GMI PRIVATE
        -Wl,--version-script,${MATLAB_ROOT}/extern/lib/glnxa64/fexport.map
        -Wl,--no-undefined
        -Wl,-rpath-link,${MATLAB_ROOT}/bin/glnxa64
    )
endif()

#########################################################################
#  Install (optional)
#########################################################################
install(TARGETS macos RUNTIME DESTINATION bin)
install(TARGETS smacos_lib ARCHIVE DESTINATION lib)
if(BUILD_SMACOS_DVR)
    install(TARGETS smacos_dvr RUNTIME DESTINATION bin)
endif()
