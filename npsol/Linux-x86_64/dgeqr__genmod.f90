        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 21 18:42:59 2022
        MODULE DGEQR__genmod
          INTERFACE 
            SUBROUTINE DGEQR(M,N,A,LDA,ZETA,INFORM)
              INTEGER(KIND=4) :: LDA
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: A(LDA,*)
              REAL(KIND=8) :: ZETA(*)
              INTEGER(KIND=4) :: INFORM
            END SUBROUTINE DGEQR
          END INTERFACE 
        END MODULE DGEQR__genmod
