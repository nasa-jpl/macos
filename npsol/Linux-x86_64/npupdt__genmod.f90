        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 21 18:43:02 2022
        MODULE NPUPDT__genmod
          INTERFACE 
            SUBROUTINE NPUPDT(LSUMRY,UNITQ,N,NCNLN,NFREE,NZ,NROWJ1,     &
     &NROWJ2,NQ,NROWR,KX,ALFA,GLF1,GLF2,QPCURV,CJAC1,CJAC2,CJDX1,CJDX2, &
     &CS1,CS2,GQ1,GQ2,HPQ,RPQ,QPMUL,R,OMEGA,ZY,WRK1,WRK2)
              INTEGER(KIND=4) :: NROWR
              INTEGER(KIND=4) :: NQ
              INTEGER(KIND=4) :: NROWJ2
              INTEGER(KIND=4) :: NROWJ1
              INTEGER(KIND=4) :: NCNLN
              INTEGER(KIND=4) :: N
              CHARACTER(LEN=4) :: LSUMRY
              LOGICAL(KIND=4) :: UNITQ
              INTEGER(KIND=4) :: NFREE
              INTEGER(KIND=4) :: NZ
              INTEGER(KIND=4) :: KX(N)
              REAL(KIND=8) :: ALFA
              REAL(KIND=8) :: GLF1
              REAL(KIND=8) :: GLF2
              REAL(KIND=8) :: QPCURV
              REAL(KIND=8) :: CJAC1(NROWJ1,*)
              REAL(KIND=8) :: CJAC2(NROWJ2,*)
              REAL(KIND=8) :: CJDX1(*)
              REAL(KIND=8) :: CJDX2(*)
              REAL(KIND=8) :: CS1(*)
              REAL(KIND=8) :: CS2(*)
              REAL(KIND=8) :: GQ1(N)
              REAL(KIND=8) :: GQ2(N)
              REAL(KIND=8) :: HPQ(N)
              REAL(KIND=8) :: RPQ(N)
              REAL(KIND=8) :: QPMUL(*)
              REAL(KIND=8) :: R(NROWR,*)
              REAL(KIND=8) :: OMEGA(*)
              REAL(KIND=8) :: ZY(NQ,*)
              REAL(KIND=8) :: WRK1(N+NCNLN)
              REAL(KIND=8) :: WRK2(N)
            END SUBROUTINE NPUPDT
          END INTERFACE 
        END MODULE NPUPDT__genmod
