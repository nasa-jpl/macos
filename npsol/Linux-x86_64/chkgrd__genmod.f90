        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 21 18:43:00 2022
        MODULE CHKGRD__genmod
          INTERFACE 
            SUBROUTINE CHKGRD(INFORM,MSGLVL,N,BIGBND,EPSRF,OKTOL,FDCHK, &
     &OBJF,XNORM,OBJFUN,BL,BU,GRAD,UGRAD,DX,X,Y,W,LENW)
              INTEGER(KIND=4) :: LENW
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: INFORM
              INTEGER(KIND=4) :: MSGLVL
              REAL(KIND=8) :: BIGBND
              REAL(KIND=8) :: EPSRF
              REAL(KIND=8) :: OKTOL
              REAL(KIND=8) :: FDCHK
              REAL(KIND=8) :: OBJF
              REAL(KIND=8) :: XNORM
              EXTERNAL OBJFUN
              REAL(KIND=8) :: BL(N)
              REAL(KIND=8) :: BU(N)
              REAL(KIND=8) :: GRAD(N)
              REAL(KIND=8) :: UGRAD(N)
              REAL(KIND=8) :: DX(N)
              REAL(KIND=8) :: X(N)
              REAL(KIND=8) :: Y(N)
              REAL(KIND=8) :: W(LENW)
            END SUBROUTINE CHKGRD
          END INTERFACE 
        END MODULE CHKGRD__genmod
