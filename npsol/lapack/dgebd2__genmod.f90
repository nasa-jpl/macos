        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 21 18:41:34 2022
        MODULE DGEBD2__genmod
          INTERFACE 
            SUBROUTINE DGEBD2(M,N,A,LDA,D,E,TAUQ,TAUP,WORK,INFO)
              INTEGER(KIND=4) :: LDA
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: A(LDA,*)
              REAL(KIND=8) :: D(*)
              REAL(KIND=8) :: E(*)
              REAL(KIND=8) :: TAUQ(*)
              REAL(KIND=8) :: TAUP(*)
              REAL(KIND=8) :: WORK(*)
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE DGEBD2
          END INTERFACE 
        END MODULE DGEBD2__genmod
