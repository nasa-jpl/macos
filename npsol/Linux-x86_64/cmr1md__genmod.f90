        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 21 18:43:00 2022
        MODULE CMR1MD__genmod
          INTERFACE 
            SUBROUTINE CMR1MD(N,NU,NRANK,NROWR,LENV,LENW,R,U,V,W)
              INTEGER(KIND=4) :: NROWR
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: NU
              INTEGER(KIND=4) :: NRANK
              INTEGER(KIND=4) :: LENV
              INTEGER(KIND=4) :: LENW
              REAL(KIND=8) :: R(NROWR,*)
              REAL(KIND=8) :: U(N,*)
              REAL(KIND=8) :: V(N)
              REAL(KIND=8) :: W(N)
            END SUBROUTINE CMR1MD
          END INTERFACE 
        END MODULE CMR1MD__genmod
