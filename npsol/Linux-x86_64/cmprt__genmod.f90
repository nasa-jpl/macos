        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 21 18:43:00 2022
        MODULE CMPRT__genmod
          INTERFACE 
            SUBROUTINE CMPRT(MSGLVL,NFREE,NROWA,N,NCLIN,NCNLN,NCTOTL,   &
     &BIGBND,NAMED,NAMES,LENNAM,NACTIV,ISTATE,KACTIV,KX,A,BL,BU,C,CLAMDA&
     &,RLAMDA,X)
              INTEGER(KIND=4) :: NCTOTL
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: NROWA
              INTEGER(KIND=4) :: MSGLVL
              INTEGER(KIND=4) :: NFREE
              INTEGER(KIND=4) :: NCLIN
              INTEGER(KIND=4) :: NCNLN
              REAL(KIND=8) :: BIGBND
              LOGICAL(KIND=4) :: NAMED
              CHARACTER(LEN=8) :: NAMES(*)
              INTEGER(KIND=4) :: LENNAM
              INTEGER(KIND=4) :: NACTIV
              INTEGER(KIND=4) :: ISTATE(NCTOTL)
              INTEGER(KIND=4) :: KACTIV(N)
              INTEGER(KIND=4) :: KX(N)
              REAL(KIND=8) :: A(NROWA,*)
              REAL(KIND=8) :: BL(NCTOTL)
              REAL(KIND=8) :: BU(NCTOTL)
              REAL(KIND=8) :: C(*)
              REAL(KIND=8) :: CLAMDA(NCTOTL)
              REAL(KIND=8) :: RLAMDA(N)
              REAL(KIND=8) :: X(N)
            END SUBROUTINE CMPRT
          END INTERFACE 
        END MODULE CMPRT__genmod
