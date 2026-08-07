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

# ==============================================================================
# Workaround: CMake 3.28.x, 4.2.x, and 4.3.x fail to set
# CMAKE_Fortran_PREPROCESS_SOURCE for Intel ifx 2026.x on Windows.
# Must be set BEFORE project() is called.
# ==============================================================================
if(WIN32)
  # Marked as a CACHE variable (not a plain set()) so it survives into
  # the isolated try_compile sub-build that CMakeTestFortranCompiler.cmake
  # spawns during ABI detection -- a plain set() before project() is not
  # reliably propagated into that sub-process across CMake versions.
  set(CMAKE_Fortran_PREPROCESS_SOURCE
    "<CMAKE_Fortran_COMPILER> <DEFINES> <INCLUDES> <FLAGS> /fpp /E <SOURCE> > <PREPROCESSED_SOURCE>"
    CACHE STRING "" FORCE)
endif()

project(SMACOS LANGUAGES Fortran C)

option(BUILD_MACOS       "Build the MACOS executable"          ON)
option(BUILD_SMACOS      "Build the SMACOS static library"     ON)
option(BUILD_SMACOS_LEAN "Build the SMACOS static library"     OFF)
option(BUILD_SMACOS_DVR  "Build the smacos_dvr test driver"    OFF)
option(BUILD_GMI         "Build the GMI Matlab mex interface"  OFF)
option(BUILD_MMACOS      "Build the mmacos MATLAB MEX interface" OFF)
set(MMACOS_FC "" CACHE STRING
    "Fortran compiler for mmacos (gfortran or ifx). Empty = inherit parent compiler.")
option(USE_NPSOL         "Build NPSOL-backed constrained optimization" OFF)
# Matlab installation root.  Auto-detect the latest installed MATLAB (mirrors
# the GMI Makefile's `ls .../R[0-9]*[ab] | sort -V | tail -1`) rather than
# hard-coding a release -- a hard-coded path breaks the GMI/mmacos cmake builds
# on any machine lacking that exact version.  Override with -DMATLAB_ROOT=...
if(NOT MATLAB_ROOT)
  if(WIN32)
    file(GLOB _matlab_roots "C:/Program Files/MATLAB/R[0-9][0-9][0-9][0-9][ab]")
  else()
    file(GLOB _matlab_roots "/usr/local/MATLAB/R[0-9][0-9][0-9][0-9][ab]")
  endif()
  list(SORT _matlab_roots)
  if(_matlab_roots)
    list(GET _matlab_roots -1 _matlab_latest)
  else()
    set(_matlab_latest "/usr/local/MATLAB/R2025b")   # last-resort fallback
  endif()
  set(MATLAB_ROOT "${_matlab_latest}" CACHE PATH "Matlab installation root")
  unset(_matlab_roots)
  unset(_matlab_latest)
endif()

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
  set(BUILD_MMACOS     OFF)
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
      # /arch:SSE2 fixes a clang-cl/icx divergence from MSVC where mmintrin.h
      # intrinsics (used by Giza/Cairo) fail with "must be vectors" errors
      # unless an explicit SIMD target is set. MSVC cl.exe doesn't need this;
      # clang-cl-based icx does. See README Troubleshooting for details.
      set(CMAKE_C_FLAGS "/nologo /D_CRT_SECURE_NO_WARNINGS /DMSWIN /arch:SSE2")
      set(SMACOS_preprocessing_Definition CSMACOS NR_FFT MSWIN)
    else()
      set(BASE_FFLAGS "/fpp /nologo /extend-source:132 /align /names:uppercase /assume:underscore /warn:nointerfaces /DPGPLOT")
      set(CMAKE_Fortran_FLAGS_RELEASE "/O2 /fp:strict /Qmkl")
      set(CMAKE_Fortran_FLAGS_DEBUG   "/Od /debug:full /traceback /check:bounds")
      set(CMAKE_Fortran_FLAGS         "${BASE_FFLAGS}")
      # /arch:SSE2 fixes a clang-cl/icx divergence from MSVC where mmintrin.h
      # intrinsics (used by Giza/Cairo) fail with "must be vectors" errors
      # unless an explicit SIMD target is set. MSVC cl.exe doesn't need this;
      # clang-cl-based icx does. See README Troubleshooting for details.
      set(CMAKE_C_FLAGS "/nologo /D_CRT_SECURE_NO_WARNINGS /DMSWIN /arch:SSE2")
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
  # Linux + macOS gfortran / gcc.  macOS/Apple-Silicon: gfortran is the only
  # Fortran path (ifx has no arm64) and is MathWorks' supported Linux compiler.
  # -DLNXOS is KEPT on macOS: the LNXOS sites (sunsub.F:117, utilsub_c.c:24,189)
  # are POSIX and Mac shares them (the C ones already co-guard __linux__).
  set(BASE_FFLAGS "-cpp -ffixed-line-length-132 -fd-lines-as-comments -fallow-argument-mismatch -std=legacy -DPGPLOT -DLNXOS -fPIC")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -fbounds-check -fcheck=all")
  set(CMAKE_Fortran_FLAGS         "${BASE_FFLAGS}")
  set(CMAKE_C_FLAGS "-DLNXOS -fPIC")
  if(APPLE)
    # Apple clang ERRORS (not warns) on the legacy C's implicit function
    # declarations / implicit int -> downgrade to warnings so utilsub_c.c,
    # mhist.c, fk.c compile.  gfortran+clang is the standard Mac sci-computing mix.
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-implicit-function-declaration -Wno-implicit-int")
  endif()
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
  ${smacos_src_dir}/glass_builtin.f90
  ${smacos_src_dir}/param_mod.F
  ${smacos_src_dir}/elt_mod.F
  ${smacos_src_dir}/plot_mode_mod.F
  ${smacos_src_dir}/nn_util.F
  ${smacos_src_dir}/mathsub.F
  ${smacos_src_dir}/macos_debug.F
  ${smacos_src_dir}/src_mod.F
  ${smacos_src_dir}/cfiles_mod.F
  ${smacos_src_dir}/sysprop_mod.F
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

# SLSQP driver -- always built; the constrained-optim path that
# macos_cmd_loop.inc dispatches into for LOptCons=.TRUE. (Phase 3 of
# the SLSQP migration; until then, design_cons_optim.F still owns the
# call site under USE_NPSOL).
list(APPEND SMACOS_ONLY_SOURCES ${smacos_src_dir}/design_slsqp_optim.F)
list(APPEND MACOS_ONLY_SOURCES  ${smacos_src_dir}/design_slsqp_optim.F)

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
extern void print_sub_prompt_(void);\n\
void SET_SUB_PROMPT_(const char *s, int slen) { set_sub_prompt_(s, slen); }\n\
void CLEAR_SUB_PROMPT_(void) { clear_sub_prompt_(); }\n\
void PRINT_SUB_PROMPT_(void) { print_sub_prompt_(); }\n\
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
# BSD license). Replacing NPSOL as MACOS's constrained-optim back end.
# Always built; no opt-in needed since the source is BSD.
set(SLSQP_SOURCES ${smacos_src_dir}/slsqp/slsqp_optmz.f)
add_library(slsqplib STATIC ${SLSQP_SOURCES})
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
  target_link_libraries(smacos PRIVATE ${NPSOL_LIBS} slsqplib ${PLATFORM_LIBS})
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
      ${NPSOL_LIBS} slsqplib
      ${PLATFORM_LIBS}
  )
endif()

# Readline (Linux + macOS)
if(APPLE)
  # The bundled readline-8.2 tree isn't pre-built on macOS; use Homebrew's
  # readline (arm64 static lib + headers).  Headers live under
  # <prefix>/include/readline/, and mhist.c does #include "readline.h" when
  # READLINE_LIBRARY is defined, so the include dir must be that subdir.
  # WITHOUT readline the engine's sub-prompt path (mhist.c set_sub_prompt) is
  # inactive and interactive prompts render blank — so this is load-bearing
  # for a usable CLI, not just arrow-key history.
  execute_process(COMMAND brew --prefix readline
                  OUTPUT_VARIABLE _brew_readline
                  OUTPUT_STRIP_TRAILING_WHITESPACE
                  ERROR_QUIET)
  if(_brew_readline AND EXISTS "${_brew_readline}/lib/libreadline.a")
    set(READLINE_ROOT "${_brew_readline}/include/readline")
    set(READLINE_LIB  "${_brew_readline}/lib/libreadline.a")
  else()
    message(WARNING "Homebrew readline not found (brew install readline) — interactive CLI prompts will render blank.")
    set(READLINE_LIB "")
  endif()
elseif(NOT WIN32)
  set(READLINE_ROOT "${smacos_src_dir}/readline-8.2")
  set(READLINE_LIB  "${READLINE_ROOT}/libreadline.a")
  if(NOT EXISTS "${READLINE_LIB}")
    # Fresh clones carry the readline-8.2 SOURCES but not the built .a —
    # build it once, in-tree, so the interactive CLI keeps its prompts
    # and history without a manual pre-build step.
    message(STATUS "Bundled libreadline.a not found — building readline-8.2 (one-time)...")
    execute_process(COMMAND sh configure --disable-shared
                    WORKING_DIRECTORY "${READLINE_ROOT}"
                    RESULT_VARIABLE _rl_cfg OUTPUT_QUIET ERROR_QUIET)
    if(_rl_cfg EQUAL 0)
      execute_process(COMMAND make -j4
                      WORKING_DIRECTORY "${READLINE_ROOT}"
                      RESULT_VARIABLE _rl_make OUTPUT_QUIET ERROR_QUIET)
    endif()
    if(NOT EXISTS "${READLINE_LIB}")
      message(WARNING "readline-8.2 auto-build FAILED — macos builds WITHOUT "
        "readline: arrow-key history is disabled and interactive sub-prompts "
        "fall back to plain printed prompts.  To restore readline: "
        "cd macos_f90/readline-8.2 && ./configure --disable-shared && make, "
        "then re-run cmake.")
      set(READLINE_LIB "")
    endif()
  endif()
else()
  set(READLINE_LIB "")
endif()
# Test hook: force the non-readline fallback path on any platform
# (-DMACOS_FORCE_NO_READLINE=ON) to exercise the printed-prompt fallback.
option(MACOS_FORCE_NO_READLINE "Build without readline (test the fallback prompt path)" OFF)
if(MACOS_FORCE_NO_READLINE)
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
      ${NPSOL_LIBS} slsqplib
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
  # smacos_dvr re-compiles macos_mod.F / macos_vars_mod.F / macosio.F /
  # pgplotsub.F (the -DCMACOS variants) as its OWN objects, so it must emit
  # its .mod files into its OWN module directory.  Sharing mod_macos with the
  # `macos` target makes two build rules generate the same
  # mod_macos/macos_mod.mod; the Ninja generator's Fortran dyndep then aborts
  # with "multiple rules generate mod_macos/macos_mod.mod" (issue #56).
  # Unix Makefiles tolerated the duplicate, which is why it stayed latent.
  set_target_properties(smacos_dvr PROPERTIES
      Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod_smacos_dvr)
  target_include_directories(smacos_dvr PRIVATE
      ${smacos_src_dir}
      ${fits_src_dir}
      ${CMAKE_BINARY_DIR}/mod_smacos_dvr
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
      ${NPSOL_LIBS} slsqplib
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
    # --disable-new-dtags is GNU ld only; macOS ld64 rejects it.  INSTALL_RPATH
    # above becomes @rpath on macOS automatically, so no extra flag is needed.
    if(NOT APPLE)
      target_link_options(smacos_dvr PRIVATE -Wl,--disable-new-dtags)
    endif()
  endif()
endif()

# ==============================================================================
# GMI MATLAB mex interface (Linux + Windows)
# Build logic lives in GMI/CMakeLists.txt — cross-platform, auto-detects MATLAB.
# ==============================================================================
if(BUILD_GMI)
  # GMI/ is at the MACOS_resources repo root, a sibling of this macos repo:
  #   <parent>/macos/CMakeLists.txt  ->  <parent>/MACOS_resources/GMI
  # Same relative layout on Linux and Windows, so no platform split is needed.
  set(_gmi_default "${CMAKE_CURRENT_SOURCE_DIR}/../MACOS_resources/GMI")
  set(GMI_SRC_DIR "${_gmi_default}" CACHE PATH
      "Directory containing GMI.F, GMIG.F and CMakeLists.txt")
  if(NOT EXISTS "${GMI_SRC_DIR}/GMI.F")
    message(FATAL_ERROR "GMI.F not found in ${GMI_SRC_DIR}
Set -DGMI_SRC_DIR=<path> on the cmake command line.")
  endif()
  add_subdirectory("${GMI_SRC_DIR}" "${CMAKE_BINARY_DIR}/gmi")
endif()

# ==============================================================================
# mmacos MATLAB MEX interface (Linux + Windows)
# Build logic lives in mmacos/CMakeLists.txt — cross-platform, auto-detects
# MATLAB.  gfortran default; override with -DCMAKE_Fortran_COMPILER=ifx.
# ==============================================================================
if(BUILD_MMACOS)
  # mmacos/ is at the MACOS_resources repo root, a sibling of this macos repo:
  #   <parent>/macos/CMakeLists.txt  ->  <parent>/MACOS_resources/mmacos
  # Same relative layout on Linux and Windows, so no platform split is needed.
  set(_mmacos_default "${CMAKE_CURRENT_SOURCE_DIR}/../MACOS_resources/mmacos")
  set(MMACOS_SRC_DIR "${_mmacos_default}" CACHE PATH
      "Directory containing mmacos's CMakeLists.txt (mmacos_mex.F may be directly here or under src/)")
  # Newer mmacos checkouts (with their own mmacos/CMakeLists.txt) keep
  # Fortran sources under mmacos/src/ rather than directly in mmacos/.
  # Check both locations so this works against older and newer layouts.
  if(NOT EXISTS "${MMACOS_SRC_DIR}/src/mmacos_mex.F" AND
     NOT EXISTS "${MMACOS_SRC_DIR}/mmacos_mex.F")
    message(FATAL_ERROR "mmacos_mex.F not found in ${MMACOS_SRC_DIR} or ${MMACOS_SRC_DIR}/src
Set -DMMACOS_SRC_DIR=<path> on the cmake command line.")
  endif()

  # Point mmacos at the current build tree so it picks up the freshly
  # built libsmacos.a / mod_smacos without requiring a separate
  # MACOS_BUILD_DIR on the command line.
  set(MACOS_BUILD_DIR "${CMAKE_BINARY_DIR}" CACHE PATH
      "Pre-built macos CMake tree consumed by mmacos (defaults to this build)" FORCE)

  # Compiler override for mmacos — allows gfortran for mmacos while the
  # parent build uses ifx (or vice versa).  Empty = inherit parent compiler.
  # Usage: -DMMACOS_FC=gfortran  or  -DMMACOS_FC=ifx
  if(NOT MMACOS_FC STREQUAL "")
    find_program(_mmacos_fc_path "${MMACOS_FC}" REQUIRED)
    set(CMAKE_Fortran_COMPILER "${_mmacos_fc_path}" CACHE FILEPATH
        "Fortran compiler for mmacos (overridden by MMACOS_FC)" FORCE)
    message(STATUS "mmacos: forcing Fortran compiler to ${_mmacos_fc_path} (MMACOS_FC=${MMACOS_FC})")
  endif()

  add_subdirectory("${MMACOS_SRC_DIR}" "${CMAKE_BINARY_DIR}/mmacos")
endif()
