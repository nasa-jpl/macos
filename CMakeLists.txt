# ------------------------------------------------------------------
# Note:
#  -- only a rough draft without using sub-project CMakeList
#     definitions
#  -- Code needs to be simplified
#
# ------------------------------------------------------------------
# Create SMACOS Libs
# ------------------------------------------------------------------
# [ ] Status:
#     [ ] Windows:
#         Require: - Intel Fortran IFX Compiler (2025.x)
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
cmake_minimum_required(VERSION 3.31)
include(CMakePrintHelpers)

project(SMACOS
  VERSION 1.0
  DESCRIPTION "SMACOS Libs"
  LANGUAGES Fortran C CXX
)
enable_language(Fortran C CXX)


# ---------------------------------------------
# Safety net
if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
  message(FATAL_ERROR
      "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there.\n"
  )
endif()


# --------------------------------------------- Module dir
set(mod_dir "${CMAKE_CURRENT_BINARY_DIR}/mod")
make_directory(${mod_dir})
set(CMAKE_Fortran_MODULE_DIRECTORY ${mod_dir})

# --------------------------------------------- OS Definition
if(WIN32)
  message("-------------------------------\nWindows x64 Settings\n-------------------------------")
  set(WINx64 TRUE)
  set(LINUX64 FALSE)
else() # if(CMAKE_SYSTEM_NAME STREQUAL "Linux")  # (CMAKE_HOST_UNIX)
  message("-------------------------------\nLinux Settings\n-------------------------------")
  set(WINx64 FALSE)
  set(LINUX64 TRUE)
endif()

# ------------------------------------------------------- Compiler
if (LINUX64)
  set(CMAKE_C_COMPILER icx)          # Intel oneAPI
  set(CMAKE_CXX_COMPILER icpx)       # Intel oneAPI
  set(CMAKE_Fortran_COMPILER ifx)    # Intel oneAPI
else()
  set(CMAKE_Fortran_COMPILER ifx)    # Intel oneAPI
  set(CMAKE_C_COMPILER cl)           # MS VC
  set(CMAKE_CXX_COMPILER cl)         # MS VC
endif()


set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS TRUE)      # chk lnx

# ------------------------------------------------------- Build Type
# Set the Build type to Release  (Default is Debug)
set(CMAKE_BUILD_TYPE Release)
message(STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")



# -------------------------------------------
# Fortran: Compiler options for all
# -------------------------------------------
set_source_files_properties(PROPERTIES Fortran_PREPROCESS ON)

# get_filename_component(Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

# if (WINx64)

  # update Fortran Flags: default: /nologo /fpp

  #    /4L132                    specify rightmost column for fixed form sources
  #    /fp:strict                Floating Point: enables /fp:precise, disables contractions and enables pragma stdc fenv_access
  #    /Qfp-speculation:strict   turns off floating point speculations
  #    /Qparallel                enable the auto-parallelizer to gen. multi-threaded code for loops that can be safely executed in parallel
  #    /MD                       use dynamically-loaded, multithread C runtime
  #    /O2                       optimize for maximum speed (DEFAULT)
  #    /Qm64                     generate code for Intel(R) 64 architecture
  #    /QxHost                   generate instructions for the highest instruction set and processor
  #                                available on the compilation host machine
  #    /Qmkl                     link to the Intel(R) Math Kernel Library (Intel(R) MKL) and bring in the associated headers
  #                                parallel - link using the threaded Intel(R) MKL libraries. (default)
  #    /Qopt-mem-layout-trans    memory layout transformations performed by the compiler: enabled
  #    /align                    specify how data items are aligned  (commons)
  #    /c                        compile to object (.obj) only, do not link
  #    /fpp                      run Fortran preprocessor on source files prior to compilation
  #    /D<name>                  Preprocessor directives
  #    /RTCu                     report use of variable that was not initialized

  # [ ] SMACOS target compile Options
  #
  #   for development:  (will disable opt)
  #     /traceback
  #     /check:bounds
  #     /check:stack
  #     /warn:interfaces
  #     /warn:unused
  #     /check:ALL
  #     /warn:ALL

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

if (WINx64)
  # Fortran Flags
  set(FFlags  "/4L132 /align /O2 /fpp /nologo")
  set(FFlags  "${FFlags} /fp:strict /Qfp-speculation:strict")
  set(FFlags  "${FFlags} /Qm64 /Qmkl /Qparallel /Qopt-mem-layout-trans /QxHost")
  set(FFlags  "${FFlags} /c /gen-interfaces") # /RTCu")

  # Fortran Flags: SMACOS  (Debug / Development)
  set(FFlags_SMACOS "${FFlags}")
  set(FFlags_SMACOS "${FFlags_SMACOS} /traceback /check:bounds /check:stack")
  set(FFlags_SMACOS "${FFlags_SMACOS} /warn:interfaces /warn:unused") # /warn:all

  # Fortran pre-processing options
  set(SMACOS_preprocessing_Definition CSMACOS NR_FFT MSWIN)

else()
  # Fortran Flags
  set(FFlags  "-c -132 -align -O2 -fpp -nologo -fPIC -gen-interfaces ")
  set(FFlags  "${FFlags} -fp-model=strict -fp-speculation=strict")
  set(FFlags  "${FFlags} -m64 -qmkl -qopt-matmul -qopt-mem-layout-trans")

  # Fortran Flags: SMACOS  (Debug)
  set(FFlags_SMACOS "${FFlags}")
  set(FFlags_SMACOS "${FFlags_SMACOS} -traceback -check bounds -check stack -fstack-protector-all")

  # Fortran pre-processing options
  set(SMACOS_preprocessing_Definition CSMACOS NR_FFT LNXOS)

endif()

set(Fortran_FLAGS_Lib        "${FFlags}")
set(Fortran_FLAGS_Lib_SMACOS "${FFlags_SMACOS}")

# cmake_print_variables(Fortran_FLAGS_Lib)
# cmake_print_variables(Fortran_FLAGS_Lib_SMACOS)


# Common variables
set(smacos_src_dir "${CMAKE_CURRENT_SOURCE_DIR}/macos_f90")
cmake_path(ABSOLUTE_PATH smacos_src_dir NORMALIZE OUTPUT_VARIABLE smacos_src_dir)  # eliminate /../

set(lib_bin_dir    "${CMAKE_CURRENT_BINARY_DIR}/libs")
set(npsol_src_dir  "${smacos_src_dir}/npsol")
set(blas_src_dir   "${npsol_src_dir}/blas")
set(lapack_src_dir "${npsol_src_dir}/lapack")


# -------------------------------------------
# Fortran Libs: Blas / LaPack / NPSol
# -------------------------------------------


# -------------------------------------------
# [ ] Lib -- NPSOL: blas
# -------------------------------------------
  message("-------------------- Lib: npsol/blas")

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
      ARCHIVE_OUTPUT_DIRECTORY ${lib_bin_dir}
      POSITION_INDEPENDENT_CODE TRUE
      COMPILE_FLAGS ${Fortran_FLAGS_Lib}
  )

  target_include_directories(blaslib
          PUBLIC
          ${blas_src_dir}
          ${lib_bin_dir}
  )


# -------------------------------------------
# [ ] Lib -- NPSOL: lapack
# -------------------------------------------
  message("-------------------- Lib: npsol/lapack")

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
      ARCHIVE_OUTPUT_DIRECTORY ${lib_bin_dir}
      POSITION_INDEPENDENT_CODE TRUE
      COMPILE_FLAGS ${Fortran_FLAGS_Lib}
  )

  target_include_directories(lapacklib
          PUBLIC
          ${lapack_src_dir}
          ${lib_bin_dir}
  )


# -------------------------------------------
# [ ] Lib -- NPSOL:
# -------------------------------------------
  message("-------------------- npsol")

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
      ARCHIVE_OUTPUT_DIRECTORY ${lib_bin_dir}
      POSITION_INDEPENDENT_CODE TRUE
      COMPILE_FLAGS ${Fortran_FLAGS_Lib}
  )

  target_include_directories(npsol
          PUBLIC
          ${npsol_src_dir}
          ${lib_bin_dir}
  )



# -------------------------------------------
# [ ] Lib -- SMACOS:
# -------------------------------------------
message("-------------------- Lib: smacos")

# set(CMAKE_Fortran_FLAGS "/modul:${mod_dir} ${CMAKE_Fortran_FLAGS_Lib}")
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

      # ${smacos_src_dir}/mhist.c
      # ${smacos_src_dir}/utilsub_c.c  # unistd.h missing

      # ${smacos_src_dir}/splicubi2.c    # needed in mathsub.F
      ${smacos_src_dir}/splicubi2.f90  # DUMMY !!!  ==> ToDo: Need to find a way to incorporate splicubi2.c with f2py !!!
      )

# needed for cmake dependencies determination (together with "Fortran_PREPROCESS ON")
# target_compile_definitions(smacos PUBLIC ${SMACOS_preprocessing_Fortran})

set_target_properties(smacos PROPERTIES
      ARCHIVE_OUTPUT_DIRECTORY ${lib_bin_dir}
      POSITION_INDEPENDENT_CODE TRUE
      COMPILE_FLAGS             "${Fortran_FLAGS_Lib_SMACOS}"
      COMPILE_DEFINITIONS       "${SMACOS_preprocessing_Definition}"
    )

target_include_directories(smacos
      PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${smacos_src_dir}
        ${mod_dir}
        ${lib_bin_dir}
        )

