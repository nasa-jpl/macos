        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 21 18:41:37 2022
        MODULE DGEQR2__genmod
          INTERFACE 
            SUBROUTINE DGEQR2(M,N,A,LDA,TAU,WORK,INFO)
              INTEGER(KIND=4) :: LDA
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: A(LDA,*)
              REAL(KIND=8) :: TAU(*)
              REAL(KIND=8) :: WORK(*)
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE DGEQR2
          END INTERFACE 
        END MODULE DGEQR2__genmod
