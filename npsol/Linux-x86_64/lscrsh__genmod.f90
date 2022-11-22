        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 21 18:43:01 2022
        MODULE LSCRSH__genmod
          INTERFACE 
            SUBROUTINE LSCRSH(COLD,VERTEX,NCLIN,NCTOTL,NACTIV,NARTIF,   &
     &NFREE,N,NROWA,ISTATE,KACTIV,BIGBND,TOLACT,A,AX,BL,BU,X,WX,WORK)
              INTEGER(KIND=4) :: NROWA
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: NCTOTL
              LOGICAL(KIND=4) :: COLD
              LOGICAL(KIND=4) :: VERTEX
              INTEGER(KIND=4) :: NCLIN
              INTEGER(KIND=4) :: NACTIV
              INTEGER(KIND=4) :: NARTIF
              INTEGER(KIND=4) :: NFREE
              INTEGER(KIND=4) :: ISTATE(NCTOTL)
              INTEGER(KIND=4) :: KACTIV(N)
              REAL(KIND=8) :: BIGBND
              REAL(KIND=8) :: TOLACT
              REAL(KIND=8) :: A(NROWA,*)
              REAL(KIND=8) :: AX(*)
              REAL(KIND=8) :: BL(NCTOTL)
              REAL(KIND=8) :: BU(NCTOTL)
              REAL(KIND=8) :: X(N)
              REAL(KIND=8) :: WX(N)
              REAL(KIND=8) :: WORK(N)
            END SUBROUTINE LSCRSH
          END INTERFACE 
        END MODULE LSCRSH__genmod
