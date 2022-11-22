        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 21 18:40:01 2022
        MODULE DTRSV__genmod
          INTERFACE 
            SUBROUTINE DTRSV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
              INTEGER(KIND=4) :: LDA
              CHARACTER(LEN=1) :: UPLO
              CHARACTER(LEN=1) :: TRANS
              CHARACTER(LEN=1) :: DIAG
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: A(LDA,*)
              REAL(KIND=8) :: X(*)
              INTEGER(KIND=4) :: INCX
            END SUBROUTINE DTRSV
          END INTERFACE 
        END MODULE DTRSV__genmod
