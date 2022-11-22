        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 21 18:43:02 2022
        MODULE NPIQP__genmod
          INTERFACE 
            SUBROUTINE NPIQP(FEASQP,UNITQ,NQPERR,MINITS,N,NCLIN,NCNLN,  &
     &NROWA,NROWJ,NROWQP,NROWR,LINACT,NLNACT,NACTIV,NFREE,NZ,NUMINF,    &
     &ISTATE,KACTIV,KX,DXNORM,GDX,QPCURV,AQP,ADX,ANORM,AX,BL,BU,C,CJAC, &
     &CLAMDA,CMUL,CS,DLAM,DSLK,DX,QPBL,QPBU,QPTOL,R,RHO,SLK,VIOLN,X,    &
     &WTINF,IW,W)
              INTEGER(KIND=4) :: NROWR
              INTEGER(KIND=4) :: NROWQP
              INTEGER(KIND=4) :: NROWJ
              INTEGER(KIND=4) :: N
              LOGICAL(KIND=4) :: FEASQP
              LOGICAL(KIND=4) :: UNITQ
              INTEGER(KIND=4) :: NQPERR
              INTEGER(KIND=4) :: MINITS
              INTEGER(KIND=4) :: NCLIN
              INTEGER(KIND=4) :: NCNLN
              INTEGER(KIND=4) :: NROWA
              INTEGER(KIND=4) :: LINACT
              INTEGER(KIND=4) :: NLNACT
              INTEGER(KIND=4) :: NACTIV
              INTEGER(KIND=4) :: NFREE
              INTEGER(KIND=4) :: NZ
              INTEGER(KIND=4) :: NUMINF
              INTEGER(KIND=4) :: ISTATE(*)
              INTEGER(KIND=4) :: KACTIV(N)
              INTEGER(KIND=4) :: KX(N)
              REAL(KIND=8) :: DXNORM
              REAL(KIND=8) :: GDX
              REAL(KIND=8) :: QPCURV
              REAL(KIND=8) :: AQP(NROWQP,*)
              REAL(KIND=8) :: ADX(*)
              REAL(KIND=8) :: ANORM(*)
              REAL(KIND=8) :: AX(*)
              REAL(KIND=8) :: BL(*)
              REAL(KIND=8) :: BU(*)
              REAL(KIND=8) :: C(*)
              REAL(KIND=8) :: CJAC(NROWJ,*)
              REAL(KIND=8) :: CLAMDA(*)
              REAL(KIND=8) :: CMUL(*)
              REAL(KIND=8) :: CS(*)
              REAL(KIND=8) :: DLAM(*)
              REAL(KIND=8) :: DSLK(*)
              REAL(KIND=8) :: DX(N)
              REAL(KIND=8) :: QPBL(*)
              REAL(KIND=8) :: QPBU(*)
              REAL(KIND=8) :: QPTOL(*)
              REAL(KIND=8) :: R(NROWR,*)
              REAL(KIND=8) :: RHO(*)
              REAL(KIND=8) :: SLK(*)
              REAL(KIND=8) :: VIOLN(*)
              REAL(KIND=8) :: X(N)
              REAL(KIND=8) :: WTINF(*)
              INTEGER(KIND=4) :: IW(*)
              REAL(KIND=8) :: W(*)
            END SUBROUTINE NPIQP
          END INTERFACE 
        END MODULE NPIQP__genmod
