cmake_minimum_required(VERSION 3.21)

# ==============================================================================
# 1. Compiler Selection & Project Setup
# ==============================================================================
if(NOT DEFINED CMAKE_Fortran_COMPILER)
  find_program(IFX_EXECUTABLE ifx)
  set(CMAKE_Fortran_COMPILER "${IFX_EXECUTABLE}" CACHE FILEPATH "Fortran compiler")
endif()

project(SMACOS LANGUAGES Fortran C)

option(BUILD_MACOS  "Build the MACOS executable"      ON)
option(BUILD_SMACOS "Build the SMACOS static library" ON)
option(USE_GIZA     "Build and link Giza instead of pgplotdummy" OFF)

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# ==============================================================================
# 2. Compiler Flag Logic (Intel ifx/icx — Windows and Linux)
# ==============================================================================
if(CMAKE_Fortran_COMPILER_ID MATCHES "IntelLLVM")
  if(WIN32)
    # Windows Intel oneAPI (ifx/icx)
    set(BASE_FFLAGS "/fpp /nologo /extend-source:132 /align /names:uppercase /assume:underscore /warn:nointerfaces /DPGPLOT")
    set(CMAKE_Fortran_FLAGS_RELEASE "/O2 /fp:strict /Qmkl")
    set(CMAKE_Fortran_FLAGS_DEBUG   "/Od /debug:full /traceback /check:bounds")
    set(CMAKE_Fortran_FLAGS         "${BASE_FFLAGS}")
    set(CMAKE_C_FLAGS "/nologo /D_CRT_SECURE_NO_WARNINGS /DMSWIN")
  else()
    # Linux Intel oneAPI (ifx/icx)
    set(BASE_FFLAGS "-fpp -extend-source 132 -align -warn nointerfaces -DPGPLOT -DLNXOS -fPIC")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -fp-model strict -qmkl")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -traceback -check bounds")
    set(CMAKE_Fortran_FLAGS         "${BASE_FFLAGS}")
    set(CMAKE_C_FLAGS "-D_CRT_SECURE_NO_WARNINGS -DLNXOS -fPIC")
  endif()
endif()

# ==============================================================================
# 3. Source Directory Configuration
# ==============================================================================
set(smacos_src_dir "${CMAKE_CURRENT_SOURCE_DIR}/macos_f90")
set(npsol_src_dir  "${smacos_src_dir}/npsol")
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

# Numerical library sources
file(GLOB BLAS_SOURCES   "${npsol_src_dir}/blas/*.f")
file(GLOB LAPACK_SOURCES "${npsol_src_dir}/lapack/*.f")
file(GLOB NPSOL_SOURCES  "${npsol_src_dir}/*.f")

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
  ${smacos_src_dir}/nn_util.F
  ${smacos_src_dir}/mathsub.F
  ${smacos_src_dir}/macos_debug.F
  ${smacos_src_dir}/src_mod.F
  ${smacos_src_dir}/cfiles_mod.F
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
  ${smacos_src_dir}/design_cons_optim.F
  ${smacos_src_dir}/stop_set.F
  ${smacos_src_dir}/psesub.F
  ${smacos_src_dir}/utilsub_c.c
  ${smacos_src_dir}/mhist.c
)

# SMACOS-only sources (smacos-specific modules + Fortran splicubi2)
set(SMACOS_ONLY_SOURCES
  ${COMMON_SOURCES}
  ${smacos_src_dir}/smacos_vars_mod.F
  ${smacos_src_dir}/smacosio.F
  ${smacos_src_dir}/smacosutil.F
  $<$<NOT:$<BOOL:${USE_GIZA}>>:${smacos_src_dir}/pgplotdummy.F>
  ${smacos_src_dir}/smacos.F
  ${smacos_src_dir}/splicubi2.f90
)

# MACOS-only sources (macos-specific modules + C splicubi2)
set(MACOS_ONLY_SOURCES
  ${COMMON_SOURCES}
  ${smacos_src_dir}/macos_vars_mod.F
  ${smacos_src_dir}/macosio.F
  ${smacos_src_dir}/pgplotsub.F
  ${smacos_src_dir}/f90prt.F
  ${smacos_src_dir}/macos.F
  ${smacos_src_dir}/splicubi2.c
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

# ==============================================================================
# 8. Targets & Linking
# ==============================================================================

# FITS Library Target
add_library(fitslib STATIC ${FITS_SOURCES})
target_include_directories(fitslib PUBLIC ${fits_src_dir} ${CMAKE_Fortran_MODULE_DIRECTORY})
if(NOT WIN32 AND CMAKE_Fortran_COMPILER_ID MATCHES "IntelLLVM")
  # fitsio.f90 is legacy code with rank/type mismatches — suppress with ifx flags
  target_compile_options(fitslib PRIVATE -mismatch -warn nousage)
endif()

# Numerical Libraries
add_library(blaslib   STATIC ${BLAS_SOURCES})
add_library(lapacklib STATIC ${LAPACK_SOURCES})
add_library(npsol     STATIC ${NPSOL_SOURCES})

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
      npsol
      lapacklib
      blaslib
      $<$<BOOL:${USE_GIZA}>:pgplot>
      ${PLATFORM_LIBS}
  )
endif()

# MACOS Executable
if(BUILD_MACOS)
  add_executable(macos ${MACOS_ONLY_SOURCES})
  target_compile_definitions(macos PRIVATE $<$<BOOL:${WIN32}>:MSWIN> CMACOS)
  set_target_properties(macos PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod_macos)
  target_include_directories(macos PUBLIC
      ${smacos_src_dir}
      ${fits_src_dir}
      ${CMAKE_BINARY_DIR}/mod_macos
  )
  if(USE_GIZA)
    target_include_directories(macos PUBLIC
        ${CMAKE_BINARY_DIR}/mod_giza
        ${CMAKE_CURRENT_SOURCE_DIR}/macos_f90/giza/src
        ${CMAKE_BINARY_DIR}/macos_f90/giza/src
    )
  endif()

  target_link_libraries(macos PRIVATE
      fitslib
      npsol
      lapacklib
      blaslib
      $<$<BOOL:${USE_GIZA}>:pgplot>
      $<$<AND:$<NOT:$<BOOL:${WIN32}>>,$<NOT:$<BOOL:${USE_GIZA}>>>:${PGPLOT_LIB}>
      $<$<AND:$<NOT:$<BOOL:${WIN32}>>,$<NOT:$<BOOL:${USE_GIZA}>>>:${X11_LIB}>
      ${PLATFORM_LIBS}
  )
endif()
