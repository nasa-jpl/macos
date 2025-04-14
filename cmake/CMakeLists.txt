# ------------------------------------------------------------------
# Create SMACOS Libs
# ------------------------------------------------------------------
# [ ] Status:
#     [ ] Windows:
#         Require: - Intel Fortran IFX Compiler
#                  - Visual Studio 2022
#
#         Info:    not compiling MACOS but SMACOS for pymacos
#
#         How:
#           ==> in cmd window with Fortran already set up
#                E.g.: C:\"Program Files (x86)"\Intel\oneAPI\setvars intel64 vs2022
#
#              mkdir build
#              cd build
#              cmake -G "NMake Makefiles" -S ..
#              nmake
#
#          ==> this will create the libs placed under .\build\libs
# ------------------------------------------------------------------
# ToDo:  [ ] define compiler flags for C / C++ compilation
#        [ ] OS dependent compilation
#            [ ] Linux
#            [ ] MacOS
#        [ ] pre-processor parameters
#        [ ] add MACOS compilation
#        [ ]
# ------------------------------------------------------------------

cmake_minimum_required(VERSION 3.31)

project(SMACOS
  VERSION 1.0
  DESCRIPTION "SMACOS Libs"
  LANGUAGES Fortran C CXX
)
enable_language(Fortran C CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Release or Debug: set to Release
# Set the Build type to Release
# if(NOT CMAKE_BUILD_TYPE)   # default is Debug
  # set(CMAKE_BUILD_TYPE Release CACHE STRING " " FORCE)
# endif()
# Set the Build type to Release
set(CMAKE_BUILD_TYPE Release)
message(STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")

# set(CMAKE_C_COMPILER "cl")  # Visual C compiler
# set(CMAKE_F_COMPILER "ifx")



# Safety net
if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
  message(
    FATAL_ERROR
      "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there.\n"
  )
endif()

# -------------------------------------------
# Fortran: Compiler options for all
# -------------------------------------------
# get_filename_component(Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

if (WIN32)
  #  Require: - Intel Fortran IFX Compiler
  #           - Visual Studio 2022 or Intel C/C++ Compiler

  set(CMAKE_Fortran_COMPILER "ifx")

  # update Fortran Flags: default: /nologo /fpp

  #    /4L132                    specify rightmost column for fixed form sources
  #    /fp:strict                Floating Point: enables /fp:precise, disables contractions and enables pragma stdc fenv_access
  #    /Qfp-speculation:strict   turns off floating point speculations
  #    /Qparallel                enable the auto-parallelizer to gen. multi-threaded code for loops that can be safely executed in parallel
  #    /MD                       use dynamically-loaded, multithread C runtime
  #    /O2                       optimize for maximum speed (DEFAULT)
  #    /Qm64                     generate code for Intel(R) 64 architecture
  #    /Qmkl                     link to the Intel(R) Math Kernel Library (Intel(R) MKL) and bring in the associated headers
  #                                parallel - link using the threaded Intel(R) MKL libraries. (default)
  #    /Qopt-mem-layout-trans    memory layout transformations performed by the compiler: enabled
  #    /align                    specify how data items are aligned  (commons)
  #    /c                        compile to object (.obj) only, do not link
  #    /fpp                      run Fortran preprocessor on source files prior to compilation
  #    /D<name>                  Preprocessor directives


  set(CMAKE_Fortran_FLAGS     "${CMAKE_Fortran_FLAGS} /4L132 /fp:strict /Qfp-speculation:strict /Qm64 /Qmkl /Qparallel /O2 /Qopt-mem-layout-trans /MD /align /fpp /nologo")
  set(CMAKE_Fortran_FLAGS_Lib "${CMAKE_Fortran_FLAGS} /c /gen-interfaces")

  # [ ] SMACOS target compile Options
  #
  #   for development:  (will disable )
  #     /traceback
  #     /check:bounds
  #     /check:stack
  #     /warn:interfaces
  #     /warn:unused
  #     /check:ALL
  #     /warn:ALL

  set(CMAKE_Fortran_FLAGS_Lib_SMACOS_Target /traceback /check:bounds /check:stack /warn:interfaces /warn:unused)

  # [ ] SMACOS/MACOS pre-processing options
  #       CMACOS      MACOS
  #       CSMACOS     SMACOS
  #       MSWIN       Target OS: Windows
  #       NR_FFT      Numerical Recipe Code
  #       FFTW_FFT
  #       INTEL_FFT
  #       SUN_FFT
  #       NLS_OPT
  #       DEBUG
  #       MACOS_CMD
  #       DREADLINE_LIBRARY

  set(SMACOS_preprocessing_CMAKE   -DCSMACOS -DMSWIN -DNR_FFT)
  set(SMACOS_preprocessing_Fortran /DCSMACOS /DMSWIN /DNR_FFT)

endif()  # WIN32



# Print out (debugging)
include(CMakePrintHelpers)
# cmake_print_variables(WIN32)
# cmake_print_variables(CMAKE_F_COMPILER)
# cmake_print_variables(CMAKE_Fortran_FLAGS)
# cmake_print_variables(PROJECT_NAME)



# Common variables
set(mod_dir "${CMAKE_BINARY_DIR}/mod")

set(smacos_src_dir "${CMAKE_SOURCE_DIR}/../macos_f90")
cmake_path(ABSOLUTE_PATH smacos_src_dir NORMALIZE OUTPUT_VARIABLE smacos_src_dir)  # eliminate /../

set(npsol_src_dir  "${smacos_src_dir}/npsol")
set(blas_src_dir   "${npsol_src_dir}/blas")
set(lapack_src_dir "${npsol_src_dir}/lapack")

make_directory(${mod_dir})
set(CMAKE_Fortran_MODULE_DIRECTORY ${mod_dir})


# -------------------------------------------
# Fortran Libs: Blas / LaPack / NPSol
# -------------------------------------------

if (1)

# -------------------------------------------
# [ ] Lib -- NPSOL: blas
# -------------------------------------------
  message("-------------------- Lib: npsol/blas")
  set(CMAKE_Fortran_FLAGS "/modul:${mod_dir} ${CMAKE_Fortran_FLAGS_Lib}")
  message(STATUS "CMAKE_Fortran_FLAGS: ${CMAKE_Fortran_FLAGS}\n")

  add_library(blaslib STATIC
         ${blas_src_dir}/daxpy.f
         ${blas_src_dir}/dcopy.f
         ${blas_src_dir}/ddot.f
         ${blas_src_dir}/dgemm.f
         ${blas_src_dir}/dgemv.f
         ${blas_src_dir}/dger.f
         ${blas_src_dir}/dnrm2.f
         ${blas_src_dir}/dscal.f
         ${blas_src_dir}/dswap.f
         ${blas_src_dir}/dsyr.f
         ${blas_src_dir}/dtrmm.f
         ${blas_src_dir}/dtrmv.f
         ${blas_src_dir}/dtrsv.f
         ${blas_src_dir}/idamax.f
         ${blas_src_dir}/lsame.f
         ${blas_src_dir}/xerbla.f
  )

  set_target_properties(blaslib PROPERTIES
      ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/libs
  )

  target_include_directories(blaslib
          PUBLIC
          ${CMAKE_SOURCE_DIR}
          ${blas_src_dir}
  )


# -------------------------------------------
# [ ] Lib -- NPSOL: lapack
# -------------------------------------------
  message("-------------------- Lib: npsol/lapack")
  set(CMAKE_Fortran_FLAGS "/modul:${mod_dir} ${CMAKE_Fortran_FLAGS_Lib}")
  message(STATUS "CMAKE_Fortran_FLAGS: ${CMAKE_Fortran_FLAGS}\n")

  add_library(lapacklib STATIC
        ${lapack_src_dir}/dbdsqr.f
        ${lapack_src_dir}/dgebd2.f
        ${lapack_src_dir}/dgebrd.f
        ${lapack_src_dir}/dgelq2.f
        ${lapack_src_dir}/dgelqf.f
        ${lapack_src_dir}/dgeqr2.f
        ${lapack_src_dir}/dgeqrf.f
        ${lapack_src_dir}/dgesvd.f
        ${lapack_src_dir}/dlabrd.f
        ${lapack_src_dir}/dlacpy.f
        ${lapack_src_dir}/dlamch.f
        ${lapack_src_dir}/dlange.f
        ${lapack_src_dir}/dlapy2.f
        ${lapack_src_dir}/dlarf.f
        ${lapack_src_dir}/dlarfb.f
        ${lapack_src_dir}/dlarfg.f
        ${lapack_src_dir}/dlarft.f
        ${lapack_src_dir}/dlartg.f
        ${lapack_src_dir}/dlas2.f
        ${lapack_src_dir}/dlascl.f
        ${lapack_src_dir}/dlaset.f
        ${lapack_src_dir}/dlasr.f
        ${lapack_src_dir}/dlassq.f
        ${lapack_src_dir}/dlasv2.f
        ${lapack_src_dir}/dorg2r.f
        ${lapack_src_dir}/dorgbr.f
        ${lapack_src_dir}/dorgl2.f
        ${lapack_src_dir}/dorglq.f
        ${lapack_src_dir}/dorgqr.f
        ${lapack_src_dir}/dorm2r.f
        ${lapack_src_dir}/dormbr.f
        ${lapack_src_dir}/dorml2.f
        ${lapack_src_dir}/dormlq.f
        ${lapack_src_dir}/dormqr.f
        ${lapack_src_dir}/ilaenv.f
  )

  set_target_properties(lapacklib PROPERTIES
      ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/libs
  )

  target_include_directories(lapacklib
          PUBLIC
          ${CMAKE_SOURCE_DIR}
          ${lapack_src_dir}
  )


# -------------------------------------------
# [ ] Lib -- NPSOL:
# -------------------------------------------
  message("-------------------- npsol")
  set(CMAKE_Fortran_FLAGS "/modul:${mod_dir} ${CMAKE_Fortran_FLAGS_Lib}")
  message(STATUS "CMAKE_Fortran_FLAGS: ${CMAKE_Fortran_FLAGS}\n")

  add_library(npsol STATIC
        ${npsol_src_dir}/blaso.f
        ${npsol_src_dir}/chsubs.f
        ${npsol_src_dir}/cmsubs.f
        ${npsol_src_dir}/lsmain.f
        ${npsol_src_dir}/lssubs.f
        ${npsol_src_dir}/mcsubs.f
        ${npsol_src_dir}/npmain.f
        ${npsol_src_dir}/npsubs.f
        ${npsol_src_dir}/opsubs.f
        ${npsol_src_dir}/srsubs.f
        )

  set_target_properties(npsol PROPERTIES
      ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/libs
  )

  target_include_directories(npsol
          PUBLIC
          ${CMAKE_SOURCE_DIR}
          ${npsol_src_dir}
  )

endif()


# -------------------------------------------
# [ ] Lib -- SMACOS:
# -------------------------------------------
message("-------------------- Lib: smacos")

set(CMAKE_Fortran_FLAGS ${CMAKE_Fortran_FLAGS_Lib})
add_library(smacos STATIC
      ${smacos_src_dir}/macos_definitions.h
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
      ${smacos_src_dir}/macos_mod.F
      ${smacos_src_dir}/dopt_mod.F
      ${smacos_src_dir}/surfsub.F
      ${smacos_src_dir}/didesub.F
      ${smacos_src_dir}/elemsub.F

      ${smacos_src_dir}/smacos_vars_mod.F
      ${smacos_src_dir}/smacosio.F          # prop_defs.inc msmacosio.inc iosub.inc
      ${smacos_src_dir}/smacos.F

      ${smacos_src_dir}/sourcsub.F   # ssrcray.inc
      ${smacos_src_dir}/smacosutil.F
      ${smacos_src_dir}/macos_init.F
      ${smacos_src_dir}/funcsub.F   # funcsub.F $(MOD_SRCS) lnk_pert.inc
      ${smacos_src_dir}/pixsub.F
      ${smacos_src_dir}/tracesub.F   #  fpp.inc  requires smacosio.F
      ${smacos_src_dir}/srtrace.F
      ${smacos_src_dir}/utilsub.F   # lensarr_indexes.inc collim.inc
      ${smacos_src_dir}/linsub.F
      ${smacos_src_dir}/pgplotdummy.F
      ${smacos_src_dir}/sunsub.F    # fits1.inc
      ${smacos_src_dir}/propsub.F
      ${smacos_src_dir}/macos_ops.F
      ${smacos_src_dir}/zern_wf.F
      ${smacos_src_dir}/design_optim.F      # smacos_compute.inc
      ${smacos_src_dir}/design_cons_optim.F # smacos_compute.inc setbeam_vars.inc setbeam.inc
      ${smacos_src_dir}/stop_set.F          # $(MOD_SRCS) design_optim.F
      ${smacos_src_dir}/psesub.F

      ${smacos_src_dir}/mhist.c
      # ${smacos_src_dir}/utilsub_c.c  # unistd.h missing

      ${smacos_src_dir}/splicubi2.c    # needed in mathsub.F
      ${smacos_src_dir}/splicubi2.f90  # DUMMY !!!  ==> ToDo: Need to find a way to incorporate splicubi2.c with f2py !!!
      )

# needed for cmake dependencies determination (together with "Fortran_PREPROCESS ON")
target_compile_definitions(smacos PUBLIC ${SMACOS_preprocessing_CMAKE})

target_compile_options(smacos
    PUBLIC
      ${SMACOS_preprocessing_Fortran}
      ${CMAKE_Fortran_FLAGS_Lib_SMACOS_Target})

set_source_files_properties(PROPERTIES Fortran_PREPROCESS ON)

# set_target_properties(smacos PROPERTIES Fortran_PREPROCESS ON)
set_target_properties(smacos PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/libs)

target_include_directories(smacos
      PUBLIC
        ${CMAKE_SOURCE_DIR}
        ${smacos_src_dir}
        ${mod_dir}
        )
 message(STATUS "CMAKE_Fortran_FLAGS: ${CMAKE_Fortran_FLAGS}\n")


