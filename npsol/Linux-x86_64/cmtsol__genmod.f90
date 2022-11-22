        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 21 18:43:00 2022
        MODULE CMTSOL__genmod
          INTERFACE 
            SUBROUTINE CMTSOL(MODE,NROWT,N,T,Y)
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: NROWT
              INTEGER(KIND=4) :: MODE
              REAL(KIND=8) :: T(NROWT,*)
              REAL(KIND=8) :: Y(N)
            END SUBROUTINE CMTSOL
          END INTERFACE 
        END MODULE CMTSOL__genmod
