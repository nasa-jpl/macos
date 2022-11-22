        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 21 18:43:02 2022
        MODULE NPFD__genmod
          INTERFACE 
            SUBROUTINE NPFD(CENTRL,INFORM,NROWJ,NROWUJ,N,NCNLN,BIGBND,  &
     &CDINT,FDINT,FDNORM,OBJF,CONFUN,OBJFUN,NEEDC,BL,BU,C,C1,C2,CJAC,   &
     &UJAC,GRAD,UGRAD,HFORWD,HCNTRL,X,W,LENW)
              INTEGER(KIND=4) :: LENW
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: NROWUJ
              INTEGER(KIND=4) :: NROWJ
              LOGICAL(KIND=4) :: CENTRL
              INTEGER(KIND=4) :: INFORM
              INTEGER(KIND=4) :: NCNLN
              REAL(KIND=8) :: BIGBND
              REAL(KIND=8) :: CDINT
              REAL(KIND=8) :: FDINT
              REAL(KIND=8) :: FDNORM
              REAL(KIND=8) :: OBJF
              EXTERNAL CONFUN
              EXTERNAL OBJFUN
              INTEGER(KIND=4) :: NEEDC(*)
              REAL(KIND=8) :: BL(N)
              REAL(KIND=8) :: BU(N)
              REAL(KIND=8) :: C(*)
              REAL(KIND=8) :: C1(*)
              REAL(KIND=8) :: C2(*)
              REAL(KIND=8) :: CJAC(NROWJ,*)
              REAL(KIND=8) :: UJAC(NROWUJ,*)
              REAL(KIND=8) :: GRAD(N)
              REAL(KIND=8) :: UGRAD(N)
              REAL(KIND=8) :: HFORWD(N)
              REAL(KIND=8) :: HCNTRL(N)
              REAL(KIND=8) :: X(N)
              REAL(KIND=8) :: W(LENW)
            END SUBROUTINE NPFD
          END INTERFACE 
        END MODULE NPFD__genmod
