        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 21 18:43:02 2022
        MODULE NPCRSH__genmod
          INTERFACE 
            SUBROUTINE NPCRSH(COLD,N,NCLIN,NCNLN,NCTOTL,NACTIV,NFREE,NZ,&
     &ISTATE,KACTIV,BIGBND,TOLACT,BL,BU,C)
              INTEGER(KIND=4) :: NCTOTL
              INTEGER(KIND=4) :: N
              LOGICAL(KIND=4) :: COLD
              INTEGER(KIND=4) :: NCLIN
              INTEGER(KIND=4) :: NCNLN
              INTEGER(KIND=4) :: NACTIV
              INTEGER(KIND=4) :: NFREE
              INTEGER(KIND=4) :: NZ
              INTEGER(KIND=4) :: ISTATE(NCTOTL)
              INTEGER(KIND=4) :: KACTIV(N)
              REAL(KIND=8) :: BIGBND
              REAL(KIND=8) :: TOLACT
              REAL(KIND=8) :: BL(NCTOTL)
              REAL(KIND=8) :: BU(NCTOTL)
              REAL(KIND=8) :: C(*)
            END SUBROUTINE NPCRSH
          END INTERFACE 
        END MODULE NPCRSH__genmod
