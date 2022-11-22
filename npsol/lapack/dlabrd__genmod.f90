        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 21 18:41:35 2022
        MODULE DLABRD__genmod
          INTERFACE 
            SUBROUTINE DLABRD(M,N,NB,A,LDA,D,E,TAUQ,TAUP,X,LDX,Y,LDY)
              INTEGER(KIND=4) :: LDY
              INTEGER(KIND=4) :: LDX
              INTEGER(KIND=4) :: LDA
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: NB
              REAL(KIND=8) :: A(LDA,*)
              REAL(KIND=8) :: D(*)
              REAL(KIND=8) :: E(*)
              REAL(KIND=8) :: TAUQ(*)
              REAL(KIND=8) :: TAUP(*)
              REAL(KIND=8) :: X(LDX,*)
              REAL(KIND=8) :: Y(LDY,*)
            END SUBROUTINE DLABRD
          END INTERFACE 
        END MODULE DLABRD__genmod
