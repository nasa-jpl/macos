# ------------------------------------------------------------------
#   => set up Intel oneAPI Env
#   => activate Python >3.13 with Numpy
#
# Windows:
#   cd build
#   cmake -G "NMake Makefiles" -S ..
#   nmake
#
# Linux:
#   cd build
#   rm -r *; cmake -DCMAKE_C_COMPILER=icx -DCMAKE_CXX_COMPILER=icpx -DCMAKE_Fortran_COMPILER=ifx -S ..
#   make
#   
# ------------------------------------------------------------------
# Test:  python -c "import pymacosf90"
#
# Two ways in:
#   standalone  -- cmake -S app/pymacos/src/cmake (above): links a macos tree
#                  built earlier, found through MACOS_BUILD_DIR.
#   included    -- the top-level CMakeLists.txt with -DBUILD_PYMACOS=ON: links
#                  the smacos target of the same build, and inherits its
#                  compilers and build type.
# ------------------------------------------------------------------

set(_pymacos_standalone OFF)
if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  set(_pymacos_standalone ON)
endif()

if(_pymacos_standalone)
cmake_minimum_required(VERSION 3.31)
include(CMakePrintHelpers)

project(pymacosf90
  VERSION 1.00.00
  DESCRIPTION "Python Interface"
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
endif()

set(pymacos_target pymacosf90)

# --------------------------------------------- Module dir
set(mod_dir "${CMAKE_CURRENT_BINARY_DIR}/mod")
make_directory(${mod_dir})
set(CMAKE_Fortran_MODULE_DIRECTORY ${mod_dir})

# --------------------------------------------- OS Definition
if(WIN32)
  message("-------------------------------\nWindows x64 Settings\n-------------------------------")
  set(WINx64 TRUE)
  set(LINUX64 FALSE)

elseif(APPLE)
  message("-------------------------------\nmacOS Settings\n-------------------------------")
  set(WINx64 FALSE)
  set(LINUX64 FALSE)
  set(MACOS64 TRUE)
else() # if(CMAKE_SYSTEM_NAME STREQUAL "Linux")  # (CMAKE_HOST_UNIX)
  message("-------------------------------\nLinux Settings\n-------------------------------")
  set(WINx64 FALSE)
  set(LINUX64 TRUE)
endif()

# ------------------------------------------------------- Compiler
# Included: the compilers are the top-level build's (set before its project()).
if (NOT _pymacos_standalone)
elseif (MACOS64)
  # Apple Silicon: gfortran only (ifx has no arm64); C stays Apple clang.
  # Use the versioned Homebrew gfortran if a bare `gfortran` isn't on PATH;
  # override with -DCMAKE_Fortran_COMPILER=... on the command line.
  if(NOT CMAKE_Fortran_COMPILER)
    find_program(_GFORTRAN NAMES gfortran gfortran-16 gfortran-15 gfortran-14 gfortran-13 REQUIRED)
    set(CMAKE_Fortran_COMPILER "${_GFORTRAN}")
  endif()
elseif (LINUX64)
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
set(CMAKE_EXPORT_COMPILE_COMMANDS TRUE)


# ------------------------------------------------------- Export
# note: it will fail to create when pymacos DLL is still in use
# The .pyd goes into the pymacos package (beside macos.py), where the tests
# and `import pymacos` look for it -- in both modes.
set(pymacos_pkg_dir "${CMAKE_CURRENT_SOURCE_DIR}/../pymacos")


# ------------------------------------------------------- Build Type
# Standalone: Release, forced (Default is Debug).  Included: the top-level
# build's CMAKE_BUILD_TYPE -- forcing the cache there would turn a debug
# build into release on the next configure.
if(_pymacos_standalone)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_BUILD_TYPE Release CACHE STRING " " FORCE)     # cache may require upd.
endif()
message(STATUS "pymacos: CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")


# ------------------------------------------------------- Python > 3.13  (requires numpy > 2.2)
# Grab Python, 3.13 or newer
find_package(Python 3.13 REQUIRED
  COMPONENTS Interpreter Development.Module NumPy)

# Grab the variables from a local Python installation
# F2PY headers
execute_process(
  COMMAND "${Python_EXECUTABLE}"
  -c "import numpy.f2py; print(numpy.f2py.get_include())"
  OUTPUT_VARIABLE F2PY_INCLUDE_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
string(REPLACE "\\" "/" F2PY_INCLUDE_DIR  "${F2PY_INCLUDE_DIR}")  # Fix path formatting for Win


# ------------------------------------------------------- Fortran Pre-Processing
set_source_files_properties(PROPERTIES Fortran_PREPROCESS ON)


# ------------------------------------------------------- Common definitions
set(f2py_module_name      ${pymacos_target})
set(pymacos_src_dir      "${CMAKE_CURRENT_SOURCE_DIR}/source")

# --- Two-file split ---------------------------------------------------
# macos_api_mod.F90 : pure-Fortran language-neutral backbone (no !f2py).
#                    NOW LIVES IN libsmacos.a (macos/src/).  Only
#                    its .mod is needed here for the f2py wrapper layer.
# pymacos_f2py.f90 : module `api` -- f2py-facing wrappers that forward
#                     into macos_api_mod via NAME_impl => NAME use-rename.
# f2py parses only pymacos_f2py.f90 (so it sees the Python signatures /
# !f2py annotations).
set(fortran_src_file      ${pymacos_src_dir}/pymacos_f2py.f90)
set(fortran_src_file_fpp  ${CMAKE_CURRENT_BINARY_DIR}/pymacos_f2py.i90)
set(f2py_module_c        "${f2py_module_name}module.c")
set(f2py_wrapper_pyf     "${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_name}.pyf")
set(f2py_wrapper_f90     "${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_name}-f2pywrappers2.f90")

if(_pymacos_standalone)
# Point at the macos build tree.  Default to ~/dev/macos/build_release;
# override via -DMACOS_BUILD_DIR=... on the cmake command line.
if(NOT DEFINED MACOS_BUILD_DIR)
  set(MACOS_BUILD_DIR "$ENV{HOME}/dev/macos/build_release" CACHE PATH
      "Path to a built macos tree (must contain mod_smacos/ and lib/libsmacos.a)")
endif()
set(macos_mod_dir "${MACOS_BUILD_DIR}/mod_smacos")
set(macos_lib_dir "${MACOS_BUILD_DIR}/lib")


# ------------------------------------------------------- SMACOS dependent libs
# NPSOL / blaslib / lapacklib removed by macos rem_npsol PR
# (compatible with macos release-candidate branch).  libsmacos.a is
# now self-contained for the symbols pymacos needs.
if(WINx64)
  set(SMACOS_Libs ${macos_lib_dir}/smacos.lib)
  if(EXISTS ${macos_lib_dir}/slsqplib.lib)
    list(APPEND SMACOS_Libs ${macos_lib_dir}/slsqplib.lib)
  endif()
else()
  set(SMACOS_Libs ${macos_lib_dir}/libsmacos.a)
  # libsmacos references slsqp_ (the SLSQP constrained-optim back end);
  # link libslsqplib.a after it when present (opt-dev/main default).
  if(EXISTS ${macos_lib_dir}/libslsqplib.a)
    list(APPEND SMACOS_Libs ${macos_lib_dir}/libslsqplib.a)
  endif()
endif()
else()
# Included: the smacos and slsqplib targets of this same build (built first,
# because they are link dependencies) and smacos's module directory.
set(macos_mod_dir "${CMAKE_BINARY_DIR}/mod_smacos")
set(macos_lib_dir "${CMAKE_BINARY_DIR}/lib")
set(SMACOS_Libs smacos slsqplib)
endif()

# ------------------------------------------------------- Fortran Flags
if(WINx64)
  set(FFlags            "/nologo /fpp /Qparallel /free /O2 /gen-interfaces")
  set(FFlags  "${FFlags} /fp:strict /Qfp-speculation:strict")
  set(FFlags  "${FFlags} /Qm64 /Qmkl /Qparallel /Qopt-mem-layout-trans /QxHost")
  set(FFlags  "${FFlags} /traceback")
  # Symbol names must match smacos.lib's: the top-level CMakeLists.txt builds
  # it with /names:uppercase /assume:underscore (MACOS_API_MOD_mp_FOO_), and
  # the GMI and mmacos stubs depend on that.  /assume:nounderscore here left
  # all 169 macos_api_mod references unresolved at link time.
  set(FFlags  "${FFlags} /names:uppercase /assume:underscore")
elseif(MACOS64)
  # gfortran flag equivalents (the ifx flags below are Intel-only).
  # pymacos_f2py.f90 is free-form; -cpp for the .F preprocessing; f2py's
  # default symbol mangling is lowercase-with-underscore, which gfortran
  # produces natively (no -names/-assume knobs needed).  -mcpu=native on arm64.
  set(FFlags "-ffree-form -O2 -fPIC -cpp -fallow-argument-mismatch -std=legacy -mcpu=native")
else()
  set(FFlags            "-nologo -free -align -O2 -fpp -fPIC -gen-interfaces ")
  set(FFlags  "${FFlags} -fp-model=strict -fp-speculation=strict")
  set(FFlags  "${FFlags} -xHost -m64 -qmkl -qopt-matmul -qopt-mem-layout-trans")
  # set(FFlags  "${FFlags} -traceback -check bounds -check stack")
  set(FFlags  "${FFlags} -names lowercase -assume underscore")
endif()
# set(Fortran_FLAGS_Lib "${FFlags}")
set(CMAKE_Fortran_FLAGS "${FFlags}")


set_source_files_properties(
  ${fortran_src_file}
  PROPERTIES Fortran_PREPROCESS ON)

# ------------------------------------------------------- Fortran Pre-processing
# Intel ifx `-P` preprocesses to a `.i90` file WITHOUT compiling.  gfortran's
# bare `-P` instead does a full compile (which would need kinds.mod etc. on the
# include path and defeats the purpose); use `-E -cpp -P -o <out>` for a pure
# preprocess-to-file that f2py can then parse.
if(MACOS64)
  add_custom_command(
    OUTPUT   ${fortran_src_file_fpp}
    COMMAND ${CMAKE_Fortran_COMPILER} -E -cpp -P ${fortran_src_file} -o ${fortran_src_file_fpp}
    DEPENDS ${fortran_src_file}
  )
else()
  add_custom_command(
    OUTPUT   ${fortran_src_file_fpp}
    COMMAND ${CMAKE_Fortran_COMPILER} -P  ${fortran_src_file}
    DEPENDS ${fortran_src_file}
  )
endif()
# -------------------------------------------------------
# --lower on every OS: the generated stubs (pymacos/api_gen.py) and macos.py
# use lowercase routine and argument names.  --no-lower on Windows kept the
# Fortran spelling (elt_srf_zrn_FreeForm, nElt_out, N), which they do not find.
set(f2py_flag "--lower")


# Generate sources
add_custom_target(genpyf
  DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_c}"
  ${f2py_wrapper_pyf}
)


add_custom_command(
  OUTPUT ${f2py_wrapper_pyf}

  COMMAND ${Python_EXECUTABLE}  -m "numpy.f2py" -h ${f2py_wrapper_pyf}
            -m ${f2py_module_name}
               ${fortran_src_file_fpp}
               ${f2py_flag}  # --no-/lower
               --overwrite-signature
               --verbose
               --skip-empty-wrappers

  DEPENDS ${fortran_src_file_fpp}
  VERBATIM
)

add_custom_command(
  OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_c}"
         ${f2py_wrapper_f90}

  COMMAND ${Python_EXECUTABLE}  -m "numpy.f2py"
                -m ${f2py_module_name}
                   ${fortran_src_file_fpp}
                   ${f2py_flag}  # --no-/lower
                   --verbose
                   --skip-empty-wrappers
  DEPENDS ${fortran_src_file_fpp} ${f2py_wrapper_pyf}  # Fortran source
  VERBATIM
)

# Set up target.  macos_api_mod's compiled object now lives inside
# libsmacos.a; only the .mod is needed at compile time (taken from
# the macos build's mod_smacos dir via target_include_directories).
Python_add_library(${pymacos_target} MODULE  WITH_SOABI  # SOABI (Shared Object Application Interface)
  ${f2py_module_c}                        # Generated
  ${f2py_wrapper_pyf}                     # Generated
  ${f2py_wrapper_f90}                     # Generated -- needed when modules are present in fortran file
  "${F2PY_INCLUDE_DIR}/fortranobject.c"   # From NumPy
  "${fortran_src_file}"                   # f2py-facing wrapper
)

if(WINx64)
  # f2py's C side names the Fortran routines with F_FUNC(f,F): with
  # UPPERCASE_FORTRAN and no NO_APPEND_FORTRAN that is F##_, matching the
  # /names:uppercase /assume:underscore Fortran flags above.
  target_compile_options(${pymacos_target} PUBLIC -DUPPERCASE_FORTRAN)
endif()

target_include_directories(${pymacos_target}
        PRIVATE
        ${F2PY_INCLUDE_DIR}
        ${Python_NumPy_INCLUDE_DIRS}
        ${CMAKE_CURRENT_BINARY_DIR}
        ${macos_lib_dir}
        ${macos_mod_dir}
        ${pymacos_src_dir}
      )

# Depend on sources
target_link_libraries(${pymacos_target}
        PRIVATE
        Python::NumPy
        ${SMACOS_Libs}
)

add_dependencies(${pymacos_target} genpyf)

# The .pyd/.so goes into the pymacos package.  The Windows debug database
# (.pdb, not gitignored) stays in the build directory.
set_target_properties(${pymacos_target} PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${pymacos_pkg_dir}"
    PDB_OUTPUT_DIRECTORY     "${CMAKE_CURRENT_BINARY_DIR}"
)
foreach(_cfg Debug Release RelWithDebInfo MinSizeRel)
  string(TOUPPER "${_cfg}" _CFG)
  set_target_properties(${pymacos_target} PROPERTIES
      LIBRARY_OUTPUT_DIRECTORY_${_CFG} "${pymacos_pkg_dir}"
      PDB_OUTPUT_DIRECTORY_${_CFG}     "${CMAKE_CURRENT_BINARY_DIR}"
  )
endforeach()
