#########################################################################
#  CMake build for NPSOL (includes bundled BLAS + LAPACK subsets)
#########################################################################
cmake_minimum_required(VERSION 3.20)

# --- BLAS subset ---
set(BLAS_SRCS
    blas/ddot.f   blas/dswap.f  blas/dtrmv.f  blas/dtrsv.f
    blas/dscal.f  blas/daxpy.f  blas/dgemv.f  blas/dnrm2.f
    blas/dcopy.f  blas/idamax.f blas/dger.f   blas/dsyr.f
    blas/xerbla.f blas/dgemm.f  blas/lsame.f  blas/dtrmm.f
)

# --- LAPACK subset ---
set(LAPACK_SRCS
    lapack/dbdsqr.f  lapack/dgebrd.f  lapack/dgelqf.f  lapack/dgeqrf.f
    lapack/dgesvd.f  lapack/dlacpy.f  lapack/dlamch.f  lapack/dlange.f
    lapack/dlascl.f  lapack/dlaset.f  lapack/dorgbr.f  lapack/dorglq.f
    lapack/dorgqr.f  lapack/dormbr.f  lapack/ilaenv.f  lapack/dlassq.f
    lapack/dlarfb.f  lapack/dlarft.f  lapack/dorgl2.f  lapack/dorg2r.f
    lapack/dormlq.f  lapack/dormqr.f  lapack/dlas2.f   lapack/dlasr.f
    lapack/dlasv2.f  lapack/dlartg.f  lapack/dgebd2.f  lapack/dlabrd.f
    lapack/dgelq2.f  lapack/dgeqr2.f  lapack/dlarf.f   lapack/dorml2.f
    lapack/dorm2r.f  lapack/dlarfg.f  lapack/dlapy2.f
)

# --- NPSOL ---
set(NPSOL_SRCS
    blaso.f  chsubs.f  cmsubs.f  lssubs.f
    mcsubs.f npsubs.f  opsubs.f  srsubs.f
)

add_library(npsol STATIC ${BLAS_SRCS} ${LAPACK_SRCS} ${NPSOL_SRCS})

# Plain .f files — no preprocessing needed
set_target_properties(npsol PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    Fortran_FORMAT FIXED
)

# Per-compiler flags for npsol (these are simple .f files, minimal flags needed)
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    target_compile_options(npsol PRIVATE -nologo -fPIC -fpp -132
        -gen-interfaces -fp-model strict -O2)
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    target_compile_options(npsol PRIVATE -fPIC -O2
        -std=legacy -fno-range-check -w)
endif()
