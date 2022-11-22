        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 21 18:43:01 2022
        MODULE LSCHOL__genmod
          INTERFACE 
            SUBROUTINE LSCHOL(NROWH,N,NRANK,TOLRNK,KX,H,INFORM)
              INTEGER(KIND=4) :: NROWH
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: NRANK
              REAL(KIND=8) :: TOLRNK
              INTEGER(KIND=4) :: KX(*)
              REAL(KIND=8) :: H(NROWH,*)
              INTEGER(KIND=4) :: INFORM
            END SUBROUTINE LSCHOL
          END INTERFACE 
        END MODULE LSCHOL__genmod
