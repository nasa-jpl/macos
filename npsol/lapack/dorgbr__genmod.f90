        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 21 18:41:18 2022
        MODULE DORGBR__genmod
          INTERFACE 
            SUBROUTINE DORGBR(VECT,M,N,K,A,LDA,TAU,WORK,LWORK,INFO)
              INTEGER(KIND=4) :: LWORK
              INTEGER(KIND=4) :: LDA
              CHARACTER(LEN=1) :: VECT
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: K
              REAL(KIND=8) :: A(LDA,*)
              REAL(KIND=8) :: TAU(*)
              REAL(KIND=8) :: WORK(LWORK)
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE DORGBR
          END INTERFACE 
        END MODULE DORGBR__genmod
