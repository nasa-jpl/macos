        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 21 18:43:00 2022
        MODULE DGEQRP__genmod
          INTERFACE 
            SUBROUTINE DGEQRP(PIVOT,M,N,A,LDA,ZETA,PERM,WORK,INFORM)
              INTEGER(KIND=4) :: LDA
              CHARACTER(LEN=1) :: PIVOT
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: A(LDA,*)
              REAL(KIND=8) :: ZETA(*)
              INTEGER(KIND=4) :: PERM(*)
              REAL(KIND=8) :: WORK(*)
              INTEGER(KIND=4) :: INFORM
            END SUBROUTINE DGEQRP
          END INTERFACE 
        END MODULE DGEQRP__genmod
