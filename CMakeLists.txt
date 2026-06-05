cmake_minimum_required(VERSION 3.21)

# ==============================================================================
# 1. Compiler Selection & Project Setup
# On Windows: ifx for Fortran, icx for C (icx supports VLAs needed by giza)
# On Linux:   ifx for Fortran, gcc for C (matches makeall.sh)
# ==============================================================================
if(NOT DEFINED CMAKE_Fortran_COMPILER)
  find_program(IFX_EXECUTABLE ifx)
  if(IFX_EXECUTABLE)
    set(CMAKE_Fortran_COMPILER "${IFX_EXECUTABLE}" CACHE FILEPATH "Fortran compiler")
  endif()
endif()

if(WIN32 AND NOT DEFINED CMAKE_C_COMPILER)
  # Use icx on Windows — supports VLAs needed by giza source files
  find_program(ICX_EXECUTABLE icx)
  if(ICX_EXECUTABLE)
    set(CMAKE_C_COMPILER "${ICX_EXECUTABLE}" CACHE FILEPATH "C compiler")
  endif()
endif()

project(SMACOS LANGUAGES Fortran C)

option(BUILD_MACOS       "Build the MACOS executable"          ON)
option(BUILD_SMACOS      "Build the SMACOS static library"     ON)
option(BUILD_SMACOS_LEAN "Build the SMACOS static library"     OFF)
option(BUILD_SMACOS_DVR  "Build the smacos_dvr test driver"    OFF)
option(BUILD_GMI         "Build the GMI Matlab mex interface"  OFF)
option(USE_NPSOL         "Build NPSOL-backed constrained optimization" OFF)
set(MATLAB_ROOT "/usr/local/MATLAB/R2025b" CACHE PATH "Matlab installation root")

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

if(BUILD_SMACOS_LEAN)
  if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
    message(FATAL_ERROR
        "In-source builds not allowed. Please make a new directory and run CMake from there.\n")
  endif()
  set(BUILD_MACOS      OFF)
  set(BUILD_SMACOS     OFF)
  set(BUILD_SMACOS_DVR OFF)
  set(BUILD_GMI        OFF)
  set(CMAKE_BUILD_TYPE Debug)
  message(STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")
endif()

# ==============================================================================
# 2. Compiler Flag Logic
# ==============================================================================
if(CMAKE_Fortran_COMPILER_ID MATCHES "IntelLLVM")
  if(WIN32)
    if(BUILD_SMACOS_LEAN)
      set(BASE_FFLAGS "/4L132 /align /O2 /fpp /nologo /warn:nointerfaces")
      set(BASE_FFLAGS "${BASE_FFLAGS} /fp:strict /Qfp-speculation:strict /Qopt-mem-layout-trans")
      set(BASE_FFLAGS "${BASE_FFLAGS} /Qmkl /c /gen-interfaces /Qm64")
      set(CMAKE_Fortran_FLAGS_RELEASE "/O2 /QxHost")
      set(CMAKE_Fortran_FLAGS_DEBUG   "/debug:full /traceback /check:bounds /warn:interfaces /warn:unused")
      set(CMAKE_Fortran_FLAGS         "${BASE_FFLAGS}")
      set(CMAKE_C_FLAGS "/nologo /D_CRT_SECURE_NO_WARNINGS /DMSWIN")
      set(SMACOS_preprocessing_Definition CSMACOS NR_FFT MSWIN)
    else()
      set(BASE_FFLAGS "/fpp /nologo /extend-source:132 /align /names:uppercase /assume:underscore /warn:nointerfaces /DPGPLOT")
      set(CMAKE_Fortran_FLAGS_RELEASE "/O2 /fp:strict /Qmkl")
      set(CMAKE_Fortran_FLAGS_DEBUG   "/Od /debug:full /traceback /check:bounds")
      set(CMAKE_Fortran_FLAGS         "${BASE_FFLAGS}")
      set(CMAKE_C_FLAGS "/nologo /D_CRT_SECURE_NO_WARNINGS /DMSWIN")
    endif()
  else()
    # Linux Intel oneAPI (ifx/gcc)
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
if(USE_NPSOL)
  set(npsol_src_dir "${smacos_src_dir}/npsol")
  if(NOT EXISTS "${npsol_src_dir}/CMakeLists.txt")
    message(FATAL_ERROR "USE_NPSOL=ON requires macos_f90/npsol/ tree.")
  endif()
endif()
set(fits_src_dir "${smacos_src_dir}/fits_build")

# ==============================================================================
# 4. Giza — always built (sole PGPLOT-API provider)
# ==============================================================================
add_subdirectory(macos_f90/giza/src)

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

if(USE_NPSOL)
  file(GLOB BLAS_SOURCES   "${npsol_src_dir}/blas/*.f")
  file(GLOB LAPACK_SOURCES "${npsol_src_dir}/lapack/*.f")
  file(GLOB NPSOL_SOURCES  "${npsol_src_dir}/*.f")
endif()

# ==============================================================================
# 6. Source Lists
# ==============================================================================
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
  ${smacos_src_dir}/nls.F
  ${smacos_src_dir}/stop_set.F
  ${smacos_src_dir}/psesub.F
  ${smacos_src_dir}/utilsub_c.c
  ${smacos_src_dir}/mhist.c
)

# SMACOS: pgplotdummy.F always included — provides pgplot_mod stubs
# needed by sunsub.F and others. Real graphics run via giza in macos/smacos_dvr.
set(SMACOS_ONLY_SOURCES
  ${COMMON_SOURCES}
  ${smacos_src_dir}/smacos_vars_mod.F
  ${smacos_src_dir}/smacosio.F
  ${smacos_src_dir}/smacosutil.F
  ${smacos_src_dir}/pgplotdummy.F
  ${smacos_src_dir}/smacos.F
  ${smacos_src_dir}/splicubi2.f90
  ${smacos_src_dir}/macos_api_mod.F90
  $<$<NOT:$<BOOL:${WIN32}>>:${smacos_src_dir}/fk.c>
)

set(MACOS_ONLY_SOURCES
  ${COMMON_SOURCES}
  ${smacos_src_dir}/macos_vars_mod.F
  ${smacos_src_dir}/macosio.F
  ${smacos_src_dir}/pgplotsub.F
  ${smacos_src_dir}/f90prt.F
  ${smacos_src_dir}/macos.F
  ${smacos_src_dir}/splicubi2.c
  $<$<NOT:$<BOOL:${WIN32}>>:${smacos_src_dir}/macos_signal.c>
)

if(USE_NPSOL)
  list(APPEND SMACOS_ONLY_SOURCES ${smacos_src_dir}/design_cons_optim.F)
  list(APPEND MACOS_ONLY_SOURCES  ${smacos_src_dir}/design_cons_optim.F)
endif()

# Windows stubs for Linux-only functions (readline + signal handling)
if(WIN32)
  file(WRITE "${CMAKE_BINARY_DIR}/windows_stubs.c"
"/* Windows stubs for Linux-only functions - generated by CMake */\n\
/* install_sigint_handler stub (macos_signal.c excluded on Windows) */\n\
void install_sigint_handler(void) {}\n\
\n\
/* Uppercase aliases required by ifx /names:uppercase /assume:underscore */\n\
/* mhist.c defines the lowercase versions; these forward to them */\n\
extern void set_sub_prompt_(const char *s, int slen);\n\
extern void clear_sub_prompt_(void);\n\
void SET_SUB_PROMPT_(const char *s, int slen) { set_sub_prompt_(s, slen); }\n\
void CLEAR_SUB_PROMPT_(void) { clear_sub_prompt_(); }\n\
"
  )
  list(APPEND MACOS_ONLY_SOURCES "${CMAKE_BINARY_DIR}/windows_stubs.c")
endif()

# ==============================================================================
# 7. Platform-specific system libraries
# ==============================================================================
if(WIN32)
  set(PLATFORM_LIBS user32.lib)
else()
  set(PLATFORM_LIBS "")
endif()

# Intel runtime libraries (Linux ifx only)
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

# FITS Library
if(NOT BUILD_SMACOS_LEAN)
  add_library(fitslib STATIC ${FITS_SOURCES})
  target_include_directories(fitslib PUBLIC ${fits_src_dir} ${CMAKE_Fortran_MODULE_DIRECTORY})
  if(NOT WIN32 AND CMAKE_Fortran_COMPILER_ID MATCHES "IntelLLVM")
    target_compile_options(fitslib PRIVATE -mismatch -warn nousage)
  endif()
endif()

# Numerical Libraries (NPSOL only)
if(USE_NPSOL)
  add_library(blaslib   STATIC ${BLAS_SOURCES})
  add_library(lapacklib STATIC ${LAPACK_SOURCES})
  add_library(npsol     STATIC ${NPSOL_SOURCES})
  set(NPSOL_LIBS npsol lapacklib blaslib)
  add_compile_definitions(USE_NPSOL)
else()
  set(NPSOL_LIBS "")
endif()

# SLSQP -- Kraft's Sequential Least Squares Programming (ACM Alg 733,
# BSD license).  Replacing NPSOL as MACOS's constrained-optim back end.
# Always built; no opt-in needed since the source is BSD.  See
# macos_f90/slsqp/LICENSE.txt.
set(SLSQP_SOURCES ${smacos_src_dir}/slsqp/slsqp_optmz.f)
add_library(slsqplib STATIC ${SLSQP_SOURCES})
# Kraft's source predates IMPLICIT NONE -- relax the strict-implicit
# warnings that the rest of macos uses, scoped to this target.
if(CMAKE_Fortran_COMPILER_ID MATCHES "IntelLLVM")
  target_compile_options(slsqplib PRIVATE "SHELL:-warn nointerfaces"
                                          "SHELL:-warn nounused")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  target_compile_options(slsqplib PRIVATE -w)
endif()

# SMACOS LEAN
if(BUILD_SMACOS_LEAN)
  add_library(smacos STATIC ${SMACOS_ONLY_SOURCES})
  target_compile_definitions(smacos PRIVATE $<$<BOOL:${WIN32}>:MSWIN> CSMACOS)
  set_target_properties(smacos PROPERTIES
      Fortran_MODULE_DIRECTORY    ${CMAKE_BINARY_DIR}/mod
      POSITION_INDEPENDENT_CODE   TRUE
      COMPILE_DEFINITIONS         "${SMACOS_preprocessing_Definition}"
  )
  target_include_directories(smacos PUBLIC
      ${smacos_src_dir}
      ${CMAKE_BINARY_DIR}/mod
  )
  target_link_libraries(smacos PRIVATE ${NPSOL_LIBS} ${PLATFORM_LIBS})
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
      ${CMAKE_BINARY_DIR}/mod_giza
      ${CMAKE_CURRENT_SOURCE_DIR}/macos_f90/giza/src
      ${CMAKE_BINARY_DIR}/macos_f90/giza/src
  )
  target_link_libraries(smacos PRIVATE
      fitslib
      ${NPSOL_LIBS}
      ${PLATFORM_LIBS}
  )
endif()

# Readline (Linux only)
if(NOT WIN32)
  set(READLINE_ROOT "${smacos_src_dir}/readline-8.2")
  set(READLINE_LIB  "${READLINE_ROOT}/libreadline.a")
  if(NOT EXISTS "${READLINE_LIB}")
    message(WARNING "Bundled libreadline.a not found — arrow-key/history support disabled.")
    set(READLINE_LIB "")
  endif()
else()
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
      ${CMAKE_BINARY_DIR}/mod_giza
      ${CMAKE_CURRENT_SOURCE_DIR}/macos_f90/giza/src
      ${CMAKE_BINARY_DIR}/macos_f90/giza/src
  )
  if(READLINE_LIB)
    target_include_directories(macos PRIVATE ${READLINE_ROOT})
    target_compile_options(macos PRIVATE
        $<$<COMPILE_LANGUAGE:C>:-DHAVE_CONFIG_H -DRL_LIBRARY_VERSION="8.2" -I${READLINE_ROOT}>
    )
  endif()
  add_dependencies(macos pgplot)
  target_link_libraries(macos PRIVATE
      fitslib
      ${NPSOL_LIBS}
      pgplot
      ${PLATFORM_LIBS}
      $<$<BOOL:${READLINE_LIB}>:${READLINE_LIB}>
      $<$<AND:$<NOT:$<BOOL:${WIN32}>>,$<BOOL:${READLINE_LIB}>>:ncurses>
  )
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
  target_compile_definitions(smacos_dvr PRIVATE $<$<BOOL:${WIN32}>:MSWIN> CMACOS)
  set_target_properties(smacos_dvr PROPERTIES
      Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod_macos)
  target_include_directories(smacos_dvr PRIVATE
      ${smacos_src_dir}
      ${fits_src_dir}
      ${CMAKE_BINARY_DIR}/mod_macos
      ${CMAKE_BINARY_DIR}/mod_smacos
      ${CMAKE_BINARY_DIR}/mod
      ${CMAKE_BINARY_DIR}/mod_giza
      ${CMAKE_CURRENT_SOURCE_DIR}/macos_f90/giza/src
      ${CMAKE_BINARY_DIR}/macos_f90/giza/src
  )
  add_dependencies(smacos_dvr macos smacos)
  target_link_libraries(smacos_dvr PRIVATE
      smacos
      fitslib
      ${NPSOL_LIBS}
      ${INTEL_LIBS}
      pgplot
      ${PLATFORM_LIBS}
  )
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

# GMI Matlab mex interface (Linux only)
if(BUILD_GMI)
  if(WIN32)
    message(WARNING "BUILD_GMI is not supported on Windows — skipping.")
  else()
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
    if(CMAKE_Fortran_COMPILER_ID MATCHES "IntelLLVM")
      target_compile_options(GMI PRIVATE
          -fPIC -fpp -132 -gen-interfaces -fp-model strict
          -fno-omit-frame-pointer -O2 -xHOST
      )
    elseif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
      target_compile_options(GMI PRIVATE
          -fPIC -cpp -ffixed-line-length-132
          -fno-omit-frame-pointer -O2 -march=native
          -fallow-argument-mismatch -std=legacy
      )
    endif()
    target_include_directories(GMI PRIVATE
        ${smacos_src_dir}
        ${CMAKE_BINARY_DIR}/mod_smacos
        ${MATLAB_ROOT}/extern/include
    )
    target_link_libraries(GMI PRIVATE
        smacos
        fitslib
        ${NPSOL_LIBS}
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
    if(CMAKE_Fortran_COMPILER_ID MATCHES "IntelLLVM")
      set_target_properties(GMI PROPERTIES
          BUILD_RPATH   "${INTEL_LIB_DIR};/opt/intel/oneapi/mkl/latest/lib"
          INSTALL_RPATH "${INTEL_LIB_DIR};/opt/intel/oneapi/mkl/latest/lib"
      )
      target_link_options(GMI PRIVATE -Wl,--disable-new-dtags)
    elseif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
      set_target_properties(GMI PROPERTIES
          BUILD_RPATH   "${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES}"
          INSTALL_RPATH "${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES}"
      )
      target_link_options(GMI PRIVATE -Wl,--disable-new-dtags)
    endif()
  endif()
endif()
