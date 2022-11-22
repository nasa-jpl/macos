        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 21 18:43:02 2022
        MODULE NPCORE__genmod
          INTERFACE 
            SUBROUTINE NPCORE(NAMED,NAMES,UNITQ,INFORM,MAJITS,N,NCLIN,  &
     &NCNLN,NCTOTL,NACTIV,NFREE,NZ,NROWA,NROWJ,NROWUJ,NROWQP,NROWR,NFUN,&
     &NGRAD,ISTATE,KACTIV,KX,OBJF,FDNORM,XNORM,OBJFUN,CONFUN,AQP,AX,BL, &
     &BU,C,CJAC,UJAC,CLAMDA,FEATOL,GRAD,UGRAD,R,X,IW,W,LENW)
              INTEGER(KIND=4) :: LENW
              INTEGER(KIND=4) :: NROWR
              INTEGER(KIND=4) :: NROWQP
              INTEGER(KIND=4) :: NROWUJ
              INTEGER(KIND=4) :: NROWJ
              INTEGER(KIND=4) :: NCTOTL
              INTEGER(KIND=4) :: N
              LOGICAL(KIND=4) :: NAMED
              CHARACTER(LEN=8) :: NAMES(*)
              LOGICAL(KIND=4) :: UNITQ
              INTEGER(KIND=4) :: INFORM
              INTEGER(KIND=4) :: MAJITS
              INTEGER(KIND=4) :: NCLIN
              INTEGER(KIND=4) :: NCNLN
              INTEGER(KIND=4) :: NACTIV
              INTEGER(KIND=4) :: NFREE
              INTEGER(KIND=4) :: NZ
              INTEGER(KIND=4) :: NROWA
              INTEGER(KIND=4) :: NFUN
              INTEGER(KIND=4) :: NGRAD
              INTEGER(KIND=4) :: ISTATE(*)
              INTEGER(KIND=4) :: KACTIV(N)
              INTEGER(KIND=4) :: KX(N)
              REAL(KIND=8) :: OBJF
              REAL(KIND=8) :: FDNORM
              REAL(KIND=8) :: XNORM
              EXTERNAL OBJFUN
              EXTERNAL CONFUN
              REAL(KIND=8) :: AQP(NROWQP,*)
              REAL(KIND=8) :: AX(*)
              REAL(KIND=8) :: BL(NCTOTL)
              REAL(KIND=8) :: BU(NCTOTL)
              REAL(KIND=8) :: C(*)
              REAL(KIND=8) :: CJAC(NROWJ,*)
              REAL(KIND=8) :: UJAC(NROWUJ,*)
              REAL(KIND=8) :: CLAMDA(NCTOTL)
              REAL(KIND=8) :: FEATOL(NCTOTL)
              REAL(KIND=8) :: GRAD(N)
              REAL(KIND=8) :: UGRAD(N)
              REAL(KIND=8) :: R(NROWR,*)
              REAL(KIND=8) :: X(N)
              INTEGER(KIND=4) :: IW(*)
              REAL(KIND=8) :: W(LENW)
            END SUBROUTINE NPCORE
          END INTERFACE 
        END MODULE NPCORE__genmod
