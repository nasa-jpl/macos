        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 21 18:43:00 2022
        MODULE CMRSWP__genmod
          INTERFACE 
            SUBROUTINE CMRSWP(N,NU,NRANK,NROWR,I,J,R,U,V)
              INTEGER(KIND=4) :: NROWR
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: NU
              INTEGER(KIND=4) :: NRANK
              INTEGER(KIND=4) :: I
              INTEGER(KIND=4) :: J
              REAL(KIND=8) :: R(NROWR,*)
              REAL(KIND=8) :: U(N,*)
              REAL(KIND=8) :: V(N)
            END SUBROUTINE CMRSWP
          END INTERFACE 
        END MODULE CMRSWP__genmod
