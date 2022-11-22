        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 21 18:40:07 2022
        MODULE DSYR__genmod
          INTERFACE 
            SUBROUTINE DSYR(UPLO,N,ALPHA,X,INCX,A,LDA)
              INTEGER(KIND=4) :: LDA
              CHARACTER(LEN=1) :: UPLO
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: ALPHA
              REAL(KIND=8) :: X(*)
              INTEGER(KIND=4) :: INCX
              REAL(KIND=8) :: A(LDA,*)
            END SUBROUTINE DSYR
          END INTERFACE 
        END MODULE DSYR__genmod
