        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 21 18:43:02 2022
        MODULE NPPRT__genmod
          INTERFACE 
            SUBROUTINE NPPRT(KTCOND,CONVRG,LSUMRY,MSGNP,MSGQP,NROWR,    &
     &NROWT,N,NCLIN,NCNLN,NCTOTL,NACTIV,LINACT,NLNACT,NZ,NFREE,MAJITS,  &
     &MINITS,ISTATE,ALFA,NFUN,CONDHZ,CONDH,CONDT,OBJALF,OBJF,GFNORM,    &
     &GZNORM,CVNORM,AX,C,R,T,VIOLN,X,WORK)
              INTEGER(KIND=4) :: NCTOTL
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: NROWT
              INTEGER(KIND=4) :: NROWR
              LOGICAL(KIND=4) :: KTCOND(2)
              LOGICAL(KIND=4) :: CONVRG
              CHARACTER(LEN=4) :: LSUMRY
              INTEGER(KIND=4) :: MSGNP
              INTEGER(KIND=4) :: MSGQP
              INTEGER(KIND=4) :: NCLIN
              INTEGER(KIND=4) :: NCNLN
              INTEGER(KIND=4) :: NACTIV
              INTEGER(KIND=4) :: LINACT
              INTEGER(KIND=4) :: NLNACT
              INTEGER(KIND=4) :: NZ
              INTEGER(KIND=4) :: NFREE
              INTEGER(KIND=4) :: MAJITS
              INTEGER(KIND=4) :: MINITS
              INTEGER(KIND=4) :: ISTATE(NCTOTL)
              REAL(KIND=8) :: ALFA
              INTEGER(KIND=4) :: NFUN
              REAL(KIND=8) :: CONDHZ
              REAL(KIND=8) :: CONDH
              REAL(KIND=8) :: CONDT
              REAL(KIND=8) :: OBJALF
              REAL(KIND=8) :: OBJF
              REAL(KIND=8) :: GFNORM
              REAL(KIND=8) :: GZNORM
              REAL(KIND=8) :: CVNORM
              REAL(KIND=8) :: AX(*)
              REAL(KIND=8) :: C(*)
              REAL(KIND=8) :: R(NROWR,*)
              REAL(KIND=8) :: T(NROWT,*)
              REAL(KIND=8) :: VIOLN(*)
              REAL(KIND=8) :: X(N)
              REAL(KIND=8) :: WORK(N)
            END SUBROUTINE NPPRT
          END INTERFACE 
        END MODULE NPPRT__genmod
