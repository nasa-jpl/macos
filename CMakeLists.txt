cmake_minimum_required(VERSION 3.21)

# ==============================================================================
# 1. Compiler Selection & Project Setup
# ==============================================================================
if(NOT DEFINED CMAKE_Fortran_COMPILER)
  find_program(IFX_EXECUTABLE ifx)
  set(CMAKE_Fortran_COMPILER "${IFX_EXECUTABLE}" CACHE FILEPATH "Fortran compiler")
endif()

project(SMACOS LANGUAGES Fortran C)

option(BUILD_MACOS       "Build the MACOS executable"          ON)
option(BUILD_SMACOS      "Build the SMACOS static library"     ON)
option(BUILD_SMACOS_LEAN "Build the SMACOS static library"     OFF)
option(BUILD_SMACOS_DVR  "Build the smacos_dvr test driver"    OFF)
option(BUILD_GMI         "Build the GMI Matlab mex interface"  OFF)
option(USE_GIZA          "Build and link Giza instead of pre-built PGPLOT" OFF)
set(MATLAB_ROOT "/usr/local/MATLAB/R2025b" CACHE PATH "Matlab installation root")

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)


if(BUILD_SMACOS_LEAN)
  # ---------------------------------------------
  # Safety net
  if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
    message(FATAL_ERROR
        "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there.\n"
    )
  endif()
  # ---------------------------------------------
  set(BUILD_MACOS       OFF)
  set(BUILD_SMACOS      OFF)
  set(BUILD_SMACOS_DVR  OFF)
  set(BUILD_GMI         OFF)
  set(USE_GIZA          OFF)

  set(CMAKE_BUILD_TYPE Debug)  # Release
  message(STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")
endif()

# ==============================================================================
# 2. Compiler Flag Logic
# ==============================================================================
if(CMAKE_Fortran_COMPILER_ID MATCHES "IntelLLVM")
  if(WIN32)
    # Windows Intel oneAPI (ifx/icx)
    if(BUILD_SMACOS_LEAN)
      set(BASE_FFLAGS "/4L132 /align /O2 /fpp /nologo /warn:nointerfaces")  # /names:uppercase /assume:nounderscore
      set(BASE_FFLAGS "${BASE_FFLAGS} /fp:strict /Qfp-speculation:strict /Qopt-mem-layout-trans")
      set(BASE_FFLAGS "${BASE_FFLAGS} /Qmkl /c /gen-interfaces /Qm64")

      set(CMAKE_Fortran_FLAGS_RELEASE "/O2 /QxHost ")
      set(CMAKE_Fortran_FLAGS_DEBUG   "/debug:full /traceback /check:bounds /warn:interfaces /warn:unused")
      set(CMAKE_Fortran_FLAGS         "${BASE_FFLAGS}")
      set(CMAKE_C_FLAGS "/nologo /D_CRT_SECURE_NO_WARNINGS /DMSWIN")

      # Fortran pre-processing options
      set(SMACOS_preprocessing_Definition CSMACOS NR_FFT MSWIN)

    else()
      set(BASE_FFLAGS "/fpp /nologo /extend-source:132 /align /names:uppercase /assume:underscore /warn:nointerfaces /DPGPLOT")
      set(CMAKE_Fortran_FLAGS_RELEASE "/O2 /fp:strict /Qmkl")
      set(CMAKE_Fortran_FLAGS_DEBUG   "/Od /debug:full /traceback /check:bounds")
      set(CMAKE_Fortran_FLAGS         "${BASE_FFLAGS}")
      set(CMAKE_C_FLAGS "/nologo /D_CRT_SECURE_NO_WARNINGS /DMSWIN")
    endif()

  else()
    # Linux Intel oneAPI (ifx/icx)
    set(BASE_FFLAGS "-fpp -extend-source 132 -align -warn nointerfaces -DPGPLOT -DLNXOS -fPIC")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -fp-model strict -qmkl")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -traceback -check bounds")
    set(CMAKE_Fortran_FLAGS         "${BASE_FFLAGS}")
    set(CMAKE_C_FLAGS "-D_CRT_SECURE_NO_WARNINGS -DLNXOS -fPIC")
  endif()
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  # Linux gfortran / gcc
  set(BASE_FFLAGS "-cpp -ffixed-line-length-132 -fd-lines-as-comments -fallow-argument-mismatch -std=legacy -DPGPLOT -DLNXOS -fPIC")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -fbounds-check -fcheck=all")
  set(CMAKE_Fortran_FLAGS         "${BASE_FFLAGS}")
  set(CMAKE_C_FLAGS "-DLNXOS -fPIC")
endif()

# ==============================================================================
# 3. Source Directory Configuration
# ==============================================================================
set(smacos_src_dir "${CMAKE_CURRENT_SOURCE_DIR}/macos_f90")
# npsol_src_dir intentionally removed: the rem_npsol PR (main) deleted
# the macos_f90/npsol/ tree.  Levenberg-Marquardt via nls_optim_dvr
# (inside design_optim_mod) is the replacement for the NPSOL-based
# constrained-optim path.
set(fits_src_dir   "${smacos_src_dir}/fits_build")

# ==============================================================================
# 4. Giza (optional — replaces pgplotdummy)
# ==============================================================================
if(USE_GIZA)
  add_subdirectory(macos_f90/giza/src)
endif()

# ==============================================================================
# 5. FITS Library
# ==============================================================================
set(FITS_SOURCES
    ${fits_src_dir}/macos_fitswrite.f90
    ${fits_src_dir}/cookbook_sub.f90
    ${fits_src_dir}/fitsfort.f90
    ${fits_src_dir}/fitsf90.f90
    ${fits_src_dir}/fitsio.f90
)

# Numerical library sources removed with rem_npsol PR -- npsol/blas/
# /lapack/ all deleted.  Targets blaslib, lapacklib, npsol below are
# also dropped; downstream target_link_libraries() calls that used to
# pull them in have been removed.

# ==============================================================================
# 6. Source Lists
# ==============================================================================

# Common Fortran/C sources shared by both macos and smacos targets
set(COMMON_SOURCES
  ${smacos_src_dir}/kinds.f90
  ${smacos_src_dir}/constants.f90
  ${smacos_src_dir}/dftsub.F
  ${smacos_src_dir}/lohpars_mod.F
  ${smacos_src_dir}/macos_IO.f90
  ${smacos_src_dir}/param_mod.F
  ${smacos_src_dir}/elt_mod.F
  ${smacos_src_dir}/plot_mode_mod.F
  ${smacos_src_dir}/nn_util.F
  ${smacos_src_dir}/mathsub.F
  ${smacos_src_dir}/macos_debug.F
  ${smacos_src_dir}/src_mod.F
  ${smacos_src_dir}/cfiles_mod.F
  ${smacos_src_dir}/validate_prescription.F90
  ${smacos_src_dir}/usersub.F
  ${smacos_src_dir}/lsq.F
  ${smacos_src_dir}/loh_mod.F
  ${smacos_src_dir}/traceutil_mod.F
  ${smacos_src_dir}/nls.F
  ${smacos_src_dir}/dopt_mod.F
  ${smacos_src_dir}/surfsub.F
  ${smacos_src_dir}/didesub.F
  ${smacos_src_dir}/elemsub.F
  ${smacos_src_dir}/sourcsub.F
  ${smacos_src_dir}/macos_mod.F
  ${smacos_src_dir}/tracesub.F
  ${smacos_src_dir}/propsub.F
  ${smacos_src_dir}/funcsub.F
  ${smacos_src_dir}/linsub.F
  ${smacos_src_dir}/macos_ops.F
  ${smacos_src_dir}/sunsub.F
  ${smacos_src_dir}/macos_init.F
  ${smacos_src_dir}/pixsub.F
  ${smacos_src_dir}/srtrace.F
  ${smacos_src_dir}/utilsub.F
  ${smacos_src_dir}/zern_wf.F
  ${smacos_src_dir}/design_optim.F
  # design_cons_optim.F removed by rem_npsol PR -- references npsol/
  # npoptn which no longer link.  nls_optim_dvr (in design_optim.F)
  # is the replacement.
  ${smacos_src_dir}/nls.F
  ${smacos_src_dir}/stop_set.F
  ${smacos_src_dir}/psesub.F
  ${smacos_src_dir}/utilsub_c.c
  ${smacos_src_dir}/mhist.c
)

# SMACOS-only sources (smacos-specific modules + Fortran splicubi2)
# pgplotdummy.F is always included for smacos: the library is computational only
# and needs MODULE pgplot_mod (stubs). Real graphics run in the driver (macos/smacos_dvr).
set(SMACOS_ONLY_SOURCES
  ${COMMON_SOURCES}
  ${smacos_src_dir}/smacos_vars_mod.F
  ${smacos_src_dir}/smacosio.F
  ${smacos_src_dir}/smacosutil.F
  ${smacos_src_dir}/pgplotdummy.F
  ${smacos_src_dir}/smacos.F
  ${smacos_src_dir}/splicubi2.f90
)
if(BUILD_GMI)
  list(APPEND SMACOS_ONLY_SOURCES ${smacos_src_dir}/fk.c)   # defines fkprocs_, needed by GMI mex link
endif()

# MACOS-only sources (macos-specific modules + C splicubi2)
set(MACOS_ONLY_SOURCES
  ${COMMON_SOURCES}
  ${smacos_src_dir}/macos_vars_mod.F
  ${smacos_src_dir}/macosio.F
  ${smacos_src_dir}/pgplotsub.F
  ${smacos_src_dir}/f90prt.F
  ${smacos_src_dir}/macos.F
  ${smacos_src_dir}/splicubi2.c
  ${smacos_src_dir}/macos_signal.c
)

# ==============================================================================
# 7. Platform-specific system libraries
# ==============================================================================
# Windows needs user32.lib for console/windowing; Linux needs nothing extra here
if(WIN32)
  set(PLATFORM_LIBS user32.lib)
else()
  set(PLATFORM_LIBS "")
endif()

# Pre-built PGPLOT library (Linux, non-Giza path)
# pgplotsub.F calls real PGPLOT symbols; must link the pre-built static lib.
# PGPLOT itself needs X11; search static first to avoid runtime path issues.
if(NOT WIN32 AND NOT USE_GIZA)
  find_library(PGPLOT_LIB libpgplot.a
      PATHS
        ${CMAKE_CURRENT_SOURCE_DIR}/macos_f90/pgplot
        /usr/local/pgplot
      NO_DEFAULT_PATH
  )
  if(NOT PGPLOT_LIB)
    message(FATAL_ERROR "Pre-built libpgplot.a not found. "
        "Build pgplot first or set USE_GIZA=ON.")
  endif()
  find_library(X11_LIB X11 PATHS /usr/lib64 /usr/lib/x86_64-linux-gnu)
  message(STATUS "Using pre-built PGPLOT: ${PGPLOT_LIB}")
  message(STATUS "Using X11: ${X11_LIB}")
endif()

# INTEL runtime libraries (ifx only — needed by smacos_dvr and GMI)
if(NOT WIN32 AND CMAKE_Fortran_COMPILER_ID MATCHES "IntelLLVM")
  set(INTEL_LIB_DIR "/opt/intel/oneapi/compiler/latest/lib" CACHE PATH
      "Intel compiler runtime library directory")
  find_library(LIBIMF    imf    PATHS ${INTEL_LIB_DIR} NO_DEFAULT_PATH)
  find_library(LIBIFPORT ifport PATHS ${INTEL_LIB_DIR} NO_DEFAULT_PATH)
  find_library(LIBIFCORE ifcore PATHS ${INTEL_LIB_DIR} NO_DEFAULT_PATH)
  set(INTEL_LIBS ${LIBIMF} ${LIBIFPORT} ${LIBIFCORE})
else()
  set(INTEL_LIBS "")
endif()

# ==============================================================================
# 8. Targets & Linking
# ==============================================================================

# FITS Library Target
if(NOT BUILD_SMACOS_LEAN)
  add_library(fitslib STATIC ${FITS_SOURCES})
  target_include_directories(fitslib PUBLIC ${fits_src_dir} ${CMAKE_Fortran_MODULE_DIRECTORY})
  if(NOT WIN32 AND CMAKE_Fortran_COMPILER_ID MATCHES "IntelLLVM")
    # fitsio.f90 is legacy code with rank/type mismatches — suppress with ifx flags
    target_compile_options(fitslib PRIVATE -mismatch -warn nousage)
  endif()
endif()

# Numerical Libraries removed with rem_npsol PR.  No blaslib /
# lapacklib / npsol targets here; SMACOS no longer needs them.

# SMACOS Static Library
if(BUILD_SMACOS_LEAN)
  add_library(smacos STATIC ${SMACOS_ONLY_SOURCES})
  target_compile_definitions(smacos PRIVATE $<$<BOOL:${WIN32}>:MSWIN> CSMACOS)
  set_target_properties(smacos PROPERTIES
      Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod
      POSITION_INDEPENDENT_CODE TRUE
      COMPILE_DEFINITIONS       "${SMACOS_preprocessing_Definition}"
  )
  target_include_directories(smacos PUBLIC
      ${smacos_src_dir}
      ${CMAKE_BINARY_DIR}/mod
  )

  target_link_libraries(smacos PRIVATE
      ${PLATFORM_LIBS}
  )
endif()

# SMACOS Static Library
if(BUILD_SMACOS)
  add_library(smacos STATIC ${SMACOS_ONLY_SOURCES})
  target_compile_definitions(smacos PRIVATE $<$<BOOL:${WIN32}>:MSWIN> CSMACOS)
  set_target_properties(smacos PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod_smacos)
  target_include_directories(smacos PUBLIC
      ${smacos_src_dir}
      ${fits_src_dir}
      ${CMAKE_BINARY_DIR}/mod_smacos
      ${CMAKE_BINARY_DIR}/mod
  )
  if(USE_GIZA)
    target_include_directories(smacos PUBLIC
        ${CMAKE_BINARY_DIR}/mod_giza
        ${CMAKE_CURRENT_SOURCE_DIR}/macos_f90/giza/src
        ${CMAKE_BINARY_DIR}/macos_f90/giza/src
    )
  endif()

  target_link_libraries(smacos PRIVATE
      fitslib
      ${PLATFORM_LIBS}
  )
endif()

# Readline (bundled) — used by macos interactive command line
set(READLINE_ROOT "${smacos_src_dir}/readline-8.2")
set(READLINE_LIB  "${READLINE_ROOT}/libreadline.a")
if(NOT EXISTS "${READLINE_LIB}")
  message(WARNING "Bundled libreadline.a not found at ${READLINE_LIB} — "
      "arrow-key/history support disabled. Build readline first: "
      "cd macos_f90/readline-8.2 && ./configure && make")
  set(READLINE_LIB "")
endif()

# MACOS Executable
if(BUILD_MACOS)
  add_executable(macos ${MACOS_ONLY_SOURCES})
  target_compile_definitions(macos PRIVATE $<$<BOOL:${WIN32}>:MSWIN> CMACOS
      $<$<BOOL:${READLINE_LIB}>:READLINE_LIBRARY>
  )
  set_target_properties(macos PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod_macos)
  target_include_directories(macos PUBLIC
      ${smacos_src_dir}
      ${fits_src_dir}
      ${CMAKE_BINARY_DIR}/mod_macos
  )
  if(READLINE_LIB)
    target_include_directories(macos PRIVATE ${READLINE_ROOT})
    target_compile_options(macos PRIVATE
        $<$<COMPILE_LANGUAGE:C>:-DHAVE_CONFIG_H -DRL_LIBRARY_VERSION="8.2" -I${READLINE_ROOT}>
    )
  endif()
  if(USE_GIZA)
    target_include_directories(macos PUBLIC
        ${CMAKE_BINARY_DIR}/mod_giza
        ${CMAKE_CURRENT_SOURCE_DIR}/macos_f90/giza/src
        ${CMAKE_BINARY_DIR}/macos_f90/giza/src
    )
  endif()

  target_link_libraries(macos PRIVATE
      fitslib
      $<$<BOOL:${USE_GIZA}>:pgplot>
      $<$<AND:$<NOT:$<BOOL:${WIN32}>>,$<NOT:$<BOOL:${USE_GIZA}>>>:${PGPLOT_LIB}>
      $<$<AND:$<NOT:$<BOOL:${WIN32}>>,$<NOT:$<BOOL:${USE_GIZA}>>>:${X11_LIB}>
      ${PLATFORM_LIBS}
      $<$<BOOL:${READLINE_LIB}>:${READLINE_LIB}>
      $<$<AND:$<NOT:$<BOOL:${WIN32}>>,$<BOOL:${READLINE_LIB}>>:ncurses>
  )

  # Embed Intel library paths in the binary so it runs without sourcing setvars.sh.
  # --disable-new-dtags forces DT_RPATH (not DT_RUNPATH) so the path propagates
  # transitively to libraries loaded by our shared libs (e.g. libimf→libintlc).
  if(NOT WIN32 AND CMAKE_Fortran_COMPILER_ID MATCHES "IntelLLVM")
    set_target_properties(macos PROPERTIES
        BUILD_RPATH   "${INTEL_LIB_DIR};/opt/intel/oneapi/mkl/latest/lib"
        INSTALL_RPATH "${INTEL_LIB_DIR};/opt/intel/oneapi/mkl/latest/lib"
    )
    target_link_options(macos PRIVATE -Wl,--disable-new-dtags)
  endif()
endif()

# SMACOS Driver Executable
if(BUILD_SMACOS_DVR)
  add_executable(smacos_dvr
      ${smacos_src_dir}/smacos_dvr.F
      ${smacos_src_dir}/macos_mod.F
      ${smacos_src_dir}/macos_vars_mod.F
      ${smacos_src_dir}/macosio.F
      ${smacos_src_dir}/pgplotsub.F
  )
  target_compile_definitions(smacos_dvr PRIVATE
      $<$<BOOL:${WIN32}>:MSWIN> CMACOS)
  set_target_properties(smacos_dvr PROPERTIES
      Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod_macos)
  target_include_directories(smacos_dvr PRIVATE
      ${smacos_src_dir}
      ${fits_src_dir}
      ${CMAKE_BINARY_DIR}/mod_macos
      ${CMAKE_BINARY_DIR}/mod_smacos
      ${CMAKE_BINARY_DIR}/mod
  )
  add_dependencies(smacos_dvr macos smacos)
  target_link_libraries(smacos_dvr PRIVATE
      smacos
      fitslib
      ${INTEL_LIBS}
      $<$<BOOL:${USE_GIZA}>:pgplot>
      $<$<AND:$<NOT:$<BOOL:${WIN32}>>,$<NOT:$<BOOL:${USE_GIZA}>>>:${PGPLOT_LIB}>
      $<$<AND:$<NOT:$<BOOL:${WIN32}>>,$<NOT:$<BOOL:${USE_GIZA}>>>:${X11_LIB}>
      ${PLATFORM_LIBS}
  )

  # Embed compiler runtime paths so smacos_dvr runs without sourcing
  # setvars.sh (ifx) or a module file (gfortran). Same pattern as macos/GMI.
  if(NOT WIN32 AND CMAKE_Fortran_COMPILER_ID MATCHES "IntelLLVM")
    set_target_properties(smacos_dvr PROPERTIES
        BUILD_RPATH   "${INTEL_LIB_DIR};/opt/intel/oneapi/mkl/latest/lib"
        INSTALL_RPATH "${INTEL_LIB_DIR};/opt/intel/oneapi/mkl/latest/lib"
    )
    target_link_options(smacos_dvr PRIVATE -Wl,--disable-new-dtags)
  elseif(NOT WIN32 AND CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set_target_properties(smacos_dvr PROPERTIES
        BUILD_RPATH   "${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES}"
        INSTALL_RPATH "${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES}"
    )
    target_link_options(smacos_dvr PRIVATE -Wl,--disable-new-dtags)
  endif()
endif()

# GMI Matlab mex interface
if(BUILD_GMI)
  set(GMI_SRC_DIR "$ENV{HOME}/dev/MACOS_resources/GMI" CACHE PATH
      "Directory containing GMI.F and GMIG.F")
  if(NOT EXISTS ${GMI_SRC_DIR}/GMI.F)
    message(FATAL_ERROR "GMI.F not found in ${GMI_SRC_DIR}")
  endif()
  add_library(GMI SHARED
      ${GMI_SRC_DIR}/GMI.F
      ${GMI_SRC_DIR}/GMIG.F
  )
  set_target_properties(GMI PROPERTIES
      OUTPUT_NAME GMI
      SUFFIX ".mexa64"
      PREFIX ""
      Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod_gmi
  )
  target_compile_definitions(GMI PRIVATE MX_COMPAT_32)
  target_compile_options(GMI PRIVATE
      -fPIC -fpp -132 -gen-interfaces -fp-model strict
      -fno-omit-frame-pointer -O2 -xHOST
  )
  target_include_directories(GMI PRIVATE
      ${smacos_src_dir}
      ${CMAKE_BINARY_DIR}/mod_smacos
      ${MATLAB_ROOT}/extern/include
  )
  target_link_libraries(GMI PRIVATE
      smacos
      fitslib
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

  # Embed compiler runtime paths in the MEX file so MATLAB can load it
  # without sourcing setvars.sh (ifx) or a module file (gfortran).
  # --disable-new-dtags forces DT_RPATH (transitive) over DT_RUNPATH
  # so paths propagate to indirect deps (e.g. libimf→libintlc).
  if(NOT WIN32 AND CMAKE_Fortran_COMPILER_ID MATCHES "IntelLLVM")
    set_target_properties(GMI PROPERTIES
        BUILD_RPATH   "${INTEL_LIB_DIR};/opt/intel/oneapi/mkl/latest/lib"
        INSTALL_RPATH "${INTEL_LIB_DIR};/opt/intel/oneapi/mkl/latest/lib"
    )
    target_link_options(GMI PRIVATE -Wl,--disable-new-dtags)
  elseif(NOT WIN32 AND CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    # gfortran: embed the directory containing libgfortran/libquadmath.
    # CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES is already a CMake list
    # (semicolon-separated), which is the form BUILD_RPATH expects.
    set_target_properties(GMI PROPERTIES
        BUILD_RPATH   "${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES}"
        INSTALL_RPATH "${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES}"
    )
    target_link_options(GMI PRIVATE -Wl,--disable-new-dtags)
  endif()
endif()
