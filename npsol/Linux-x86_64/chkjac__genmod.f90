        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 21 18:43:00 2022
        MODULE CHKJAC__genmod
          INTERFACE 
            SUBROUTINE CHKJAC(INFORM,LVLDER,MSGLVL,NCSET,N,NCNLN,NROWJ, &
     &NROWUJ,BIGBND,EPSRF,OKTOL,FDCHK,XNORM,CONFUN,NEEDC,BL,BU,C,C1,CJAC&
     &,UJAC,CJDX,DX,ERR,X,Y,W,LENW)
              INTEGER(KIND=4) :: LENW
              INTEGER(KIND=4) :: NROWUJ
              INTEGER(KIND=4) :: NROWJ
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: INFORM
              INTEGER(KIND=4) :: LVLDER
              INTEGER(KIND=4) :: MSGLVL
              INTEGER(KIND=4) :: NCSET
              INTEGER(KIND=4) :: NCNLN
              REAL(KIND=8) :: BIGBND
              REAL(KIND=8) :: EPSRF
              REAL(KIND=8) :: OKTOL
              REAL(KIND=8) :: FDCHK
              REAL(KIND=8) :: XNORM
              EXTERNAL CONFUN
              INTEGER(KIND=4) :: NEEDC(*)
              REAL(KIND=8) :: BL(N)
              REAL(KIND=8) :: BU(N)
              REAL(KIND=8) :: C(*)
              REAL(KIND=8) :: C1(*)
              REAL(KIND=8) :: CJAC(NROWJ,*)
              REAL(KIND=8) :: UJAC(NROWUJ,*)
              REAL(KIND=8) :: CJDX(*)
              REAL(KIND=8) :: DX(N)
              REAL(KIND=8) :: ERR(*)
              REAL(KIND=8) :: X(N)
              REAL(KIND=8) :: Y(N)
              REAL(KIND=8) :: W(LENW)
            END SUBROUTINE CHKJAC
          END INTERFACE 
        END MODULE CHKJAC__genmod
