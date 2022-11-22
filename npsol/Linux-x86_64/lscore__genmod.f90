        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 21 18:43:01 2022
        MODULE LSCORE__genmod
          INTERFACE 
            SUBROUTINE LSCORE(PRBTYP,NAMED,NAMES,LINOBJ,UNITQ,INFORM,   &
     &ITER,JINF,NCLIN,NCTOTL,NACTIV,NFREE,NRANK,NZ,NZ1,N,NROWA,NROWR,   &
     &ISTATE,KACTIV,KX,CTX,SSQ,SSQ1,SUMINF,NUMINF,XNORM,BL,BU,A,CLAMDA, &
     &AX,FEATOL,R,X,IW,W)
              INTEGER(KIND=4) :: NROWR
              INTEGER(KIND=4) :: NROWA
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: NCTOTL
              CHARACTER(LEN=2) :: PRBTYP
              LOGICAL(KIND=4) :: NAMED
              CHARACTER(LEN=8) :: NAMES(*)
              LOGICAL(KIND=4) :: LINOBJ
              LOGICAL(KIND=4) :: UNITQ
              INTEGER(KIND=4) :: INFORM
              INTEGER(KIND=4) :: ITER
              INTEGER(KIND=4) :: JINF
              INTEGER(KIND=4) :: NCLIN
              INTEGER(KIND=4) :: NACTIV
              INTEGER(KIND=4) :: NFREE
              INTEGER(KIND=4) :: NRANK
              INTEGER(KIND=4) :: NZ
              INTEGER(KIND=4) :: NZ1
              INTEGER(KIND=4) :: ISTATE(NCTOTL)
              INTEGER(KIND=4) :: KACTIV(N)
              INTEGER(KIND=4) :: KX(N)
              REAL(KIND=8) :: CTX
              REAL(KIND=8) :: SSQ
              REAL(KIND=8) :: SSQ1
              REAL(KIND=8) :: SUMINF
              INTEGER(KIND=4) :: NUMINF
              REAL(KIND=8) :: XNORM
              REAL(KIND=8) :: BL(NCTOTL)
              REAL(KIND=8) :: BU(NCTOTL)
              REAL(KIND=8) :: A(NROWA,*)
              REAL(KIND=8) :: CLAMDA(NCTOTL)
              REAL(KIND=8) :: AX(*)
              REAL(KIND=8) :: FEATOL(NCTOTL)
              REAL(KIND=8) :: R(NROWR,*)
              REAL(KIND=8) :: X(N)
              INTEGER(KIND=4) :: IW(*)
              REAL(KIND=8) :: W(*)
            END SUBROUTINE LSCORE
          END INTERFACE 
        END MODULE LSCORE__genmod
