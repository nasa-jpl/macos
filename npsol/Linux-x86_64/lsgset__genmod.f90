        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 21 18:43:01 2022
        MODULE LSGSET__genmod
          INTERFACE 
            SUBROUTINE LSGSET(PRBTYP,LINOBJ,SINGLR,UNITGZ,UNITQ,N,NCLIN,&
     &NFREE,NROWA,NQ,NROWR,NRANK,NZ,NZ1,ISTATE,KX,BIGBND,TOLRNK,NUMINF, &
     &SUMINF,BL,BU,A,RES,FEATOL,GQ,CQ,R,X,WTINF,ZY,WRK)
              INTEGER(KIND=4) :: NROWR
              INTEGER(KIND=4) :: NQ
              INTEGER(KIND=4) :: NROWA
              INTEGER(KIND=4) :: N
              CHARACTER(LEN=2) :: PRBTYP
              LOGICAL(KIND=4) :: LINOBJ
              LOGICAL(KIND=4) :: SINGLR
              LOGICAL(KIND=4) :: UNITGZ
              LOGICAL(KIND=4) :: UNITQ
              INTEGER(KIND=4) :: NCLIN
              INTEGER(KIND=4) :: NFREE
              INTEGER(KIND=4) :: NRANK
              INTEGER(KIND=4) :: NZ
              INTEGER(KIND=4) :: NZ1
              INTEGER(KIND=4) :: ISTATE(*)
              INTEGER(KIND=4) :: KX(N)
              REAL(KIND=8) :: BIGBND
              REAL(KIND=8) :: TOLRNK
              INTEGER(KIND=4) :: NUMINF
              REAL(KIND=8) :: SUMINF
              REAL(KIND=8) :: BL(*)
              REAL(KIND=8) :: BU(*)
              REAL(KIND=8) :: A(NROWA,*)
              REAL(KIND=8) :: RES(*)
              REAL(KIND=8) :: FEATOL(*)
              REAL(KIND=8) :: GQ(N)
              REAL(KIND=8) :: CQ(*)
              REAL(KIND=8) :: R(NROWR,*)
              REAL(KIND=8) :: X(N)
              REAL(KIND=8) :: WTINF(*)
              REAL(KIND=8) :: ZY(NQ,*)
              REAL(KIND=8) :: WRK(N)
            END SUBROUTINE LSGSET
          END INTERFACE 
        END MODULE LSGSET__genmod
