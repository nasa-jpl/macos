        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 21 18:43:01 2022
        MODULE LSGETP__genmod
          INTERFACE 
            SUBROUTINE LSGETP(LINOBJ,SINGLR,UNITGZ,UNITQ,N,NCLIN,NFREE, &
     &NROWA,NQ,NROWR,NRANK,NUMINF,NZ1,ISTATE,KX,CTP,PNORM,A,AP,RES,HZ,P,&
     &GQ,CQ,R,ZY,WORK)
              INTEGER(KIND=4) :: NROWR
              INTEGER(KIND=4) :: NQ
              INTEGER(KIND=4) :: NROWA
              INTEGER(KIND=4) :: NCLIN
              INTEGER(KIND=4) :: N
              LOGICAL(KIND=4) :: LINOBJ
              LOGICAL(KIND=4) :: SINGLR
              LOGICAL(KIND=4) :: UNITGZ
              LOGICAL(KIND=4) :: UNITQ
              INTEGER(KIND=4) :: NFREE
              INTEGER(KIND=4) :: NRANK
              INTEGER(KIND=4) :: NUMINF
              INTEGER(KIND=4) :: NZ1
              INTEGER(KIND=4) :: ISTATE(N+NCLIN)
              INTEGER(KIND=4) :: KX(N)
              REAL(KIND=8) :: CTP
              REAL(KIND=8) :: PNORM
              REAL(KIND=8) :: A(NROWA,*)
              REAL(KIND=8) :: AP(*)
              REAL(KIND=8) :: RES(*)
              REAL(KIND=8) :: HZ(*)
              REAL(KIND=8) :: P(N)
              REAL(KIND=8) :: GQ(N)
              REAL(KIND=8) :: CQ(*)
              REAL(KIND=8) :: R(NROWR,*)
              REAL(KIND=8) :: ZY(NQ,*)
              REAL(KIND=8) :: WORK(N)
            END SUBROUTINE LSGETP
          END INTERFACE 
        END MODULE LSGETP__genmod
