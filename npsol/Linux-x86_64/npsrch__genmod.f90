        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 21 18:43:02 2022
        MODULE NPSRCH__genmod
          INTERFACE 
            SUBROUTINE NPSRCH(NEEDFD,INFORM,N,NCNLN,NROWJ,NROWUJ,NFUN,  &
     &NGRAD,NEEDC,CONFUN,OBJFUN,ALFA,ALFBND,ALFMAX,ALFSML,DXNORM,EPSRF, &
     &ETA,GDX,GRDALF,GLF1,GLF,OBJF,OBJALF,QPCURV,XNORM,C,CJAC,UJAC,CJDX,&
     &CMUL1,CMUL,CS1,CS,DX,DLAM,DSLK,GRAD,UGRAD,QPMUL,RHO,SLK1,SLK,X1,X,&
     &W,LENW)
              INTEGER(KIND=4) :: LENW
              INTEGER(KIND=4) :: NROWUJ
              INTEGER(KIND=4) :: NROWJ
              INTEGER(KIND=4) :: N
              LOGICAL(KIND=4) :: NEEDFD
              INTEGER(KIND=4) :: INFORM
              INTEGER(KIND=4) :: NCNLN
              INTEGER(KIND=4) :: NFUN
              INTEGER(KIND=4) :: NGRAD
              INTEGER(KIND=4) :: NEEDC(*)
              EXTERNAL CONFUN
              EXTERNAL OBJFUN
              REAL(KIND=8) :: ALFA
              REAL(KIND=8) :: ALFBND
              REAL(KIND=8) :: ALFMAX
              REAL(KIND=8) :: ALFSML
              REAL(KIND=8) :: DXNORM
              REAL(KIND=8) :: EPSRF
              REAL(KIND=8) :: ETA
              REAL(KIND=8) :: GDX
              REAL(KIND=8) :: GRDALF
              REAL(KIND=8) :: GLF1
              REAL(KIND=8) :: GLF
              REAL(KIND=8) :: OBJF
              REAL(KIND=8) :: OBJALF
              REAL(KIND=8) :: QPCURV
              REAL(KIND=8) :: XNORM
              REAL(KIND=8) :: C(*)
              REAL(KIND=8) :: CJAC(NROWJ,*)
              REAL(KIND=8) :: UJAC(NROWUJ,*)
              REAL(KIND=8) :: CJDX(*)
              REAL(KIND=8) :: CMUL1(*)
              REAL(KIND=8) :: CMUL(*)
              REAL(KIND=8) :: CS1(*)
              REAL(KIND=8) :: CS(*)
              REAL(KIND=8) :: DX(N)
              REAL(KIND=8) :: DLAM(*)
              REAL(KIND=8) :: DSLK(*)
              REAL(KIND=8) :: GRAD(N)
              REAL(KIND=8) :: UGRAD(N)
              REAL(KIND=8) :: QPMUL(*)
              REAL(KIND=8) :: RHO(*)
              REAL(KIND=8) :: SLK1(*)
              REAL(KIND=8) :: SLK(*)
              REAL(KIND=8) :: X1(N)
              REAL(KIND=8) :: X(N)
              REAL(KIND=8) :: W(LENW)
            END SUBROUTINE NPSRCH
          END INTERFACE 
        END MODULE NPSRCH__genmod
