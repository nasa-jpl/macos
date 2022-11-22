        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 21 18:41:29 2022
        MODULE DORMQR__genmod
          INTERFACE 
            SUBROUTINE DORMQR(SIDE,TRANS,M,N,K,A,LDA,TAU,C,LDC,WORK,    &
     &LWORK,INFO)
              INTEGER(KIND=4) :: LWORK
              INTEGER(KIND=4) :: LDC
              INTEGER(KIND=4) :: LDA
              CHARACTER(LEN=1) :: SIDE
              CHARACTER(LEN=1) :: TRANS
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: K
              REAL(KIND=8) :: A(LDA,*)
              REAL(KIND=8) :: TAU(*)
              REAL(KIND=8) :: C(LDC,*)
              REAL(KIND=8) :: WORK(LWORK)
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE DORMQR
          END INTERFACE 
        END MODULE DORMQR__genmod
