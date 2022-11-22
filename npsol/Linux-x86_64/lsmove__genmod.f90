        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 21 18:43:01 2022
        MODULE LSMOVE__genmod
          INTERFACE 
            SUBROUTINE LSMOVE(HITCON,HITLOW,LINOBJ,UNITGZ,NCLIN,NRANK,  &
     &NZ1,N,NROWR,JADD,NUMINF,ALFA,CTP,CTX,XNORM,AP,AX,BL,BU,GQ,HZ,P,RES&
     &,R,X,WORK)
              INTEGER(KIND=4) :: NROWR
              INTEGER(KIND=4) :: N
              LOGICAL(KIND=4) :: HITCON
              LOGICAL(KIND=4) :: HITLOW
              LOGICAL(KIND=4) :: LINOBJ
              LOGICAL(KIND=4) :: UNITGZ
              INTEGER(KIND=4) :: NCLIN
              INTEGER(KIND=4) :: NRANK
              INTEGER(KIND=4) :: NZ1
              INTEGER(KIND=4) :: JADD
              INTEGER(KIND=4) :: NUMINF
              REAL(KIND=8) :: ALFA
              REAL(KIND=8) :: CTP
              REAL(KIND=8) :: CTX
              REAL(KIND=8) :: XNORM
              REAL(KIND=8) :: AP(*)
              REAL(KIND=8) :: AX(*)
              REAL(KIND=8) :: BL(*)
              REAL(KIND=8) :: BU(*)
              REAL(KIND=8) :: GQ(*)
              REAL(KIND=8) :: HZ(*)
              REAL(KIND=8) :: P(N)
              REAL(KIND=8) :: RES(*)
              REAL(KIND=8) :: R(NROWR,*)
              REAL(KIND=8) :: X(N)
              REAL(KIND=8) :: WORK(*)
            END SUBROUTINE LSMOVE
          END INTERFACE 
        END MODULE LSMOVE__genmod
