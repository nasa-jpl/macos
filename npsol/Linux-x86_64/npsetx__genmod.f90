        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 21 18:43:02 2022
        MODULE NPSETX__genmod
          INTERFACE 
            SUBROUTINE NPSETX(UNITQ,NCQP,NACTIV,NFREE,NZ,N,NLNX,NCTOTL, &
     &NQ,NROWQP,NROWR,NROWT,ISTATE,KACTIV,KX,DXNORM,GDX,AQP,ADX,BL,BU,  &
     &RPQ,RPQ0,DX,GQ,R,T,ZY,WORK)
              INTEGER(KIND=4) :: NROWT
              INTEGER(KIND=4) :: NROWR
              INTEGER(KIND=4) :: NROWQP
              INTEGER(KIND=4) :: NQ
              INTEGER(KIND=4) :: NCTOTL
              INTEGER(KIND=4) :: NLNX
              INTEGER(KIND=4) :: N
              LOGICAL(KIND=4) :: UNITQ
              INTEGER(KIND=4) :: NCQP
              INTEGER(KIND=4) :: NACTIV
              INTEGER(KIND=4) :: NFREE
              INTEGER(KIND=4) :: NZ
              INTEGER(KIND=4) :: ISTATE(NCTOTL)
              INTEGER(KIND=4) :: KACTIV(N)
              INTEGER(KIND=4) :: KX(N)
              REAL(KIND=8) :: DXNORM
              REAL(KIND=8) :: GDX
              REAL(KIND=8) :: AQP(NROWQP,*)
              REAL(KIND=8) :: ADX(*)
              REAL(KIND=8) :: BL(NCTOTL)
              REAL(KIND=8) :: BU(NCTOTL)
              REAL(KIND=8) :: RPQ(NLNX)
              REAL(KIND=8) :: RPQ0(NLNX)
              REAL(KIND=8) :: DX(N)
              REAL(KIND=8) :: GQ(N)
              REAL(KIND=8) :: R(NROWR,*)
              REAL(KIND=8) :: T(NROWT,*)
              REAL(KIND=8) :: ZY(NQ,*)
              REAL(KIND=8) :: WORK(N)
            END SUBROUTINE NPSETX
          END INTERFACE 
        END MODULE NPSETX__genmod
