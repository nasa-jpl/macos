        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 21 18:41:27 2022
        MODULE DORGL2__genmod
          INTERFACE 
            SUBROUTINE DORGL2(M,N,K,A,LDA,TAU,WORK,INFO)
              INTEGER(KIND=4) :: LDA
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: K
              REAL(KIND=8) :: A(LDA,*)
              REAL(KIND=8) :: TAU(*)
              REAL(KIND=8) :: WORK(*)
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE DORGL2
          END INTERFACE 
        END MODULE DORGL2__genmod
