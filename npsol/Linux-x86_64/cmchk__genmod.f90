        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 21 18:43:00 2022
        MODULE CMCHK__genmod
          INTERFACE 
            SUBROUTINE CMCHK(NERROR,MSGLVL,COLD,USERKX,LIWORK,LWORK,    &
     &LITOTL,LWTOTL,N,NCLIN,NCNLN,ISTATE,KX,NAMED,NAMES,LENNAM,BL,BU,X)
              INTEGER(KIND=4) :: NCNLN
              INTEGER(KIND=4) :: NCLIN
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: NERROR
              INTEGER(KIND=4) :: MSGLVL
              LOGICAL(KIND=4) :: COLD
              LOGICAL(KIND=4) :: USERKX
              INTEGER(KIND=4) :: LIWORK
              INTEGER(KIND=4) :: LWORK
              INTEGER(KIND=4) :: LITOTL
              INTEGER(KIND=4) :: LWTOTL
              INTEGER(KIND=4) :: ISTATE(N+NCLIN+NCNLN)
              INTEGER(KIND=4) :: KX(N)
              LOGICAL(KIND=4) :: NAMED
              CHARACTER(LEN=8) :: NAMES(*)
              INTEGER(KIND=4) :: LENNAM
              REAL(KIND=8) :: BL(N+NCLIN+NCNLN)
              REAL(KIND=8) :: BU(N+NCLIN+NCNLN)
              REAL(KIND=8) :: X(N)
            END SUBROUTINE CMCHK
          END INTERFACE 
        END MODULE CMCHK__genmod
