        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 21 18:42:59 2022
        MODULE DCOND__genmod
          INTERFACE 
            SUBROUTINE DCOND(N,X,INCX,AXMAX,AXMIN)
              INTEGER(KIND=4) :: INCX
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: X((N-1)*INCX+1)
              REAL(KIND=8) :: AXMAX
              REAL(KIND=8) :: AXMIN
            END SUBROUTINE DCOND
          END INTERFACE 
        END MODULE DCOND__genmod
