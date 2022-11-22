        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 21 18:43:00 2022
        MODULE DGEAP__genmod
          INTERFACE 
            SUBROUTINE DGEAP(SIDE,TRANS,M,N,PERM,K,B,LDB)
              INTEGER(KIND=4) :: LDB
              CHARACTER(LEN=1) :: SIDE
              CHARACTER(LEN=1) :: TRANS
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: PERM(*)
              INTEGER(KIND=4) :: K
              REAL(KIND=8) :: B(LDB,*)
            END SUBROUTINE DGEAP
          END INTERFACE 
        END MODULE DGEAP__genmod
