        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 21 18:43:01 2022
        MODULE LSDEL__genmod
          INTERFACE 
            SUBROUTINE LSDEL(UNITQ,N,NACTIV,NFREE,NRES,NGQ,NZ,NZ1,NROWA,&
     &NQ,NROWR,NROWT,NRANK,JDEL,KDEL,KACTIV,KX,A,RES,R,T,GQ,ZY,WORK)
              INTEGER(KIND=4) :: NROWT
              INTEGER(KIND=4) :: NROWR
              INTEGER(KIND=4) :: NQ
              INTEGER(KIND=4) :: NROWA
              INTEGER(KIND=4) :: N
              LOGICAL(KIND=4) :: UNITQ
              INTEGER(KIND=4) :: NACTIV
              INTEGER(KIND=4) :: NFREE
              INTEGER(KIND=4) :: NRES
              INTEGER(KIND=4) :: NGQ
              INTEGER(KIND=4) :: NZ
              INTEGER(KIND=4) :: NZ1
              INTEGER(KIND=4) :: NRANK
              INTEGER(KIND=4) :: JDEL
              INTEGER(KIND=4) :: KDEL
              INTEGER(KIND=4) :: KACTIV(N)
              INTEGER(KIND=4) :: KX(N)
              REAL(KIND=8) :: A(NROWA,*)
              REAL(KIND=8) :: RES(N,*)
              REAL(KIND=8) :: R(NROWR,*)
              REAL(KIND=8) :: T(NROWT,*)
              REAL(KIND=8) :: GQ(N,*)
              REAL(KIND=8) :: ZY(NQ,*)
              REAL(KIND=8) :: WORK(N)
            END SUBROUTINE LSDEL
          END INTERFACE 
        END MODULE LSDEL__genmod
