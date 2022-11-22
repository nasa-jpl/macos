        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 21 18:41:13 2022
        MODULE DGEBRD__genmod
          INTERFACE 
            SUBROUTINE DGEBRD(M,N,A,LDA,D,E,TAUQ,TAUP,WORK,LWORK,INFO)
              INTEGER(KIND=4) :: LWORK
              INTEGER(KIND=4) :: LDA
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: A(LDA,*)
              REAL(KIND=8) :: D(*)
              REAL(KIND=8) :: E(*)
              REAL(KIND=8) :: TAUQ(*)
              REAL(KIND=8) :: TAUP(*)
              REAL(KIND=8) :: WORK(LWORK)
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE DGEBRD
          END INTERFACE 
        END MODULE DGEBRD__genmod
