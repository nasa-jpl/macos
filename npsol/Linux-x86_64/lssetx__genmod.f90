        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 21 18:43:01 2022
        MODULE LSSETX__genmod
          INTERFACE 
            SUBROUTINE LSSETX(LINOBJ,ROWERR,UNITQ,NCLIN,NACTIV,NFREE,   &
     &NRANK,NZ,N,NCTOTL,NQ,NROWA,NROWR,NROWT,ISTATE,KACTIV,KX,JMAX,     &
     &ERRMAX,CTX,XNORM,A,AX,BL,BU,CQ,RES,RES0,FEATOL,R,T,X,ZY,P,WORK)
              INTEGER(KIND=4) :: NROWT
              INTEGER(KIND=4) :: NROWR
              INTEGER(KIND=4) :: NROWA
              INTEGER(KIND=4) :: NQ
              INTEGER(KIND=4) :: NCTOTL
              INTEGER(KIND=4) :: N
              LOGICAL(KIND=4) :: LINOBJ
              LOGICAL(KIND=4) :: ROWERR
              LOGICAL(KIND=4) :: UNITQ
              INTEGER(KIND=4) :: NCLIN
              INTEGER(KIND=4) :: NACTIV
              INTEGER(KIND=4) :: NFREE
              INTEGER(KIND=4) :: NRANK
              INTEGER(KIND=4) :: NZ
              INTEGER(KIND=4) :: ISTATE(NCTOTL)
              INTEGER(KIND=4) :: KACTIV(N)
              INTEGER(KIND=4) :: KX(N)
              INTEGER(KIND=4) :: JMAX
              REAL(KIND=8) :: ERRMAX
              REAL(KIND=8) :: CTX
              REAL(KIND=8) :: XNORM
              REAL(KIND=8) :: A(NROWA,*)
              REAL(KIND=8) :: AX(*)
              REAL(KIND=8) :: BL(NCTOTL)
              REAL(KIND=8) :: BU(NCTOTL)
              REAL(KIND=8) :: CQ(*)
              REAL(KIND=8) :: RES(*)
              REAL(KIND=8) :: RES0(*)
              REAL(KIND=8) :: FEATOL(NCTOTL)
              REAL(KIND=8) :: R(NROWR,*)
              REAL(KIND=8) :: T(NROWT,*)
              REAL(KIND=8) :: X(N)
              REAL(KIND=8) :: ZY(NQ,*)
              REAL(KIND=8) :: P(N)
              REAL(KIND=8) :: WORK(NCTOTL)
            END SUBROUTINE LSSETX
          END INTERFACE 
        END MODULE LSSETX__genmod
