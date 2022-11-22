        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 21 18:41:21 2022
        MODULE DORGQR__genmod
          INTERFACE 
            SUBROUTINE DORGQR(M,N,K,A,LDA,TAU,WORK,LWORK,INFO)
              INTEGER(KIND=4) :: LWORK
              INTEGER(KIND=4) :: LDA
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: K
              REAL(KIND=8) :: A(LDA,*)
              REAL(KIND=8) :: TAU(*)
              REAL(KIND=8) :: WORK(LWORK)
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE DORGQR
          END INTERFACE 
        END MODULE DORGQR__genmod
