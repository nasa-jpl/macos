        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 21 18:43:02 2022
        MODULE NPCHKD__genmod
          INTERFACE 
            SUBROUTINE NPCHKD(INFORM,MSGNP,NSTATE,LVLDER,NFUN,NGRAD,    &
     &NROWJ,NROWUJ,N,NCNLN,CONFUN,OBJFUN,NEEDC,BIGBND,EPSRF,CDINT,FDINT,&
     &FDCHK,FDNORM,OBJF,XNORM,BL,BU,C,C1,CJAC,UJAC,CJDX,DX,GRAD,UGRAD,  &
     &HFORWD,HCNTRL,X,WRK1,WRK2,W,LENW)
              INTEGER(KIND=4) :: LENW
              INTEGER(KIND=4) :: NCNLN
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: NROWUJ
              INTEGER(KIND=4) :: NROWJ
              INTEGER(KIND=4) :: INFORM
              INTEGER(KIND=4) :: MSGNP
              INTEGER(KIND=4) :: NSTATE
              INTEGER(KIND=4) :: LVLDER
              INTEGER(KIND=4) :: NFUN
              INTEGER(KIND=4) :: NGRAD
              EXTERNAL CONFUN
              EXTERNAL OBJFUN
              INTEGER(KIND=4) :: NEEDC(*)
              REAL(KIND=8) :: BIGBND
              REAL(KIND=8) :: EPSRF
              REAL(KIND=8) :: CDINT
              REAL(KIND=8) :: FDINT
              REAL(KIND=8) :: FDCHK
              REAL(KIND=8) :: FDNORM
              REAL(KIND=8) :: OBJF
              REAL(KIND=8) :: XNORM
              REAL(KIND=8) :: BL(N)
              REAL(KIND=8) :: BU(N)
              REAL(KIND=8) :: C(*)
              REAL(KIND=8) :: C1(*)
              REAL(KIND=8) :: CJAC(NROWJ,*)
              REAL(KIND=8) :: UJAC(NROWUJ,*)
              REAL(KIND=8) :: CJDX(*)
              REAL(KIND=8) :: DX(N)
              REAL(KIND=8) :: GRAD(N)
              REAL(KIND=8) :: UGRAD(N)
              REAL(KIND=8) :: HFORWD(*)
              REAL(KIND=8) :: HCNTRL(*)
              REAL(KIND=8) :: X(N)
              REAL(KIND=8) :: WRK1(N+NCNLN)
              REAL(KIND=8) :: WRK2(N+NCNLN)
              REAL(KIND=8) :: W(LENW)
            END SUBROUTINE NPCHKD
          END INTERFACE 
        END MODULE NPCHKD__genmod
