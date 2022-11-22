        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 21 18:43:02 2022
        MODULE NPALF__genmod
          INTERFACE 
            SUBROUTINE NPALF(INFORM,N,NCLIN,NCNLN,ALFA,ALFMIN,ALFMAX,   &
     &BIGBND,DXNORM,ANORM,ADX,AX,BL,BU,DSLK,DX,SLK,X)
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: INFORM
              INTEGER(KIND=4) :: NCLIN
              INTEGER(KIND=4) :: NCNLN
              REAL(KIND=8) :: ALFA
              REAL(KIND=8) :: ALFMIN
              REAL(KIND=8) :: ALFMAX
              REAL(KIND=8) :: BIGBND
              REAL(KIND=8) :: DXNORM
              REAL(KIND=8) :: ANORM(*)
              REAL(KIND=8) :: ADX(*)
              REAL(KIND=8) :: AX(*)
              REAL(KIND=8) :: BL(*)
              REAL(KIND=8) :: BU(*)
              REAL(KIND=8) :: DSLK(*)
              REAL(KIND=8) :: DX(N)
              REAL(KIND=8) :: SLK(*)
              REAL(KIND=8) :: X(N)
            END SUBROUTINE NPALF
          END INTERFACE 
        END MODULE NPALF__genmod
