C***********************************************************************
C       Begin file dftsub.F
C    +----------------------------------------------------------------+
C    |  Copyright (C) 1995-2007, California Institute of Technology.  |
C    |  U.S. Government Sponsorship Is Acknowledged.                  |
C    +----------------------------------------------------------------+
C***********************************************************************

	SUBROUTINE NFPropDFT(u1,n1,dx1,u2,n2,dx2,lambda,z)

C  Near field propagator using DFT 
C       u1 = complex field at first surface
C       n1 = dimension of field at first surface
C      dx1 = grid spacing at first surface
C       u2 = complex field at second surface
C       n2 = dimension of field at second surface
C      dx2 = grid spacing at second surface
C   lambda = wavelength
C	 z = propagation distance

	IMPLICIT NONE

C  These are local variables
	INTEGER n1,n2,ii
     	REAL*8 lambda,z,dx1,dx2,x1,c2,dfx,pi,twopi, C11
	COMPLEX*16 c1,u1(n1,n1),u2(n2,n2)

	DATA pi/3.141592653589793d0/
	DATA twopi/6.283185307179586476925286766559005768394d0/

C Apply factors to input array
	c2 = pi/lambda/z
	CALL applyfac1(u1,n1,n1,c2,dx1)

C Compute DFT
	dfx=dx2/lambda/z
	x1=-dx1*INT(n1/2)

	CALL DFT2D(u1,n1,n1,dx1,dx1,x1,x1,1,
     &	u2,n2,n2,dfx,dfx)

C Apply factors to output array
	C11=DBLE(N2*N2)
C	c1 = CDEXP(DCMPLX(0d0,twopi*z/lambda))/DCMPLX(0d0,lambda*z)
	c1 = 1d0/DCMPLX(0d0,lambda*z)
	CALL applyfac2(u2,n2,n2,c1*dx1*dx1,c2,dx2)

	RETURN
	END SUBROUTINE NFPropDFT 

C***********************************************************************

	SUBROUTINE FFPropDFT(u1,n1,dx1,u2,n2,dx2,lambda,z)

C  Far field propagator using DFT 
C       u1 = complex field at first surface
C       n1 = dimension of field at first surface
C      dx1 = grid spacing at first surface
C       u2 = complex field at second surface
C       n2 = dimension of field at second surface
C      dx2 = grid spacing at second surface
C   lambda = wavelength
C	 z = propagation distance

	IMPLICIT NONE

C  These are local variables
	INTEGER n1,n2
     	REAL*8 lambda,z,dx1,dx2,x1,c2,pi,twopi,dfx, C11
	COMPLEX*16 c1,u1(n1,n1),u2(n2,n2)

	DATA pi/3.141592653589793d0/
	DATA twopi/6.283185307179586476925286766559005768394d0/

C Compute DFT
	dfx=ABS(dx2/lambda/z)
	x1=-dx1*INT(n1/2)

	CALL DFT2D(u1,n1,n1,dx1,dx1,x1,x1,1,
     &	u2,n2,n2,dfx,dfx)

C Apply factors to output array
        C11=DBLE(N2*N2)
C	c1 = CDEXP(DCMPLX(0d0,twopi*z/lambda))/DCMPLX(0d0,lambda*z)
	c1 = 1d0/DCMPLX(0d0,lambda*z)
	c2 = pi/lambda/z
	CALL applyfac2(u2,n2,n2,c1*dx1*dx1,c2,dx2)

	RETURN
	END SUBROUTINE FFPropDFT

C***********************************************************************

	SUBROUTINE Applyfac1(u,mi,mj,c2,dx)
	
	IMPLICIT NONE

	INTEGER	i,j,mi,mj,i1,j1
	REAL*8 c2,dx,x,y,ysq
	COMPLEX*16 u(mi,mj)

	i1=mi/2+1
	j1=mj/2+1
	y =-dx*DBLE(j1)
	DO j=1,mj
	  y = y+dx
	  ysq=y*y
	  x =-dx*DBLE(i1)
	  DO i=1,mi
	    x =	x+dx 
	    u(i,j)=u(i,j)*CDEXP(DCMPLX(0d0,c2*(x*x+ysq))) 
	  END DO
	END DO

	RETURN
	END SUBROUTINE Applyfac1
	
C***********************************************************************

	SUBROUTINE Applyfac2(u,mi,mj,c1dxdy,c2,dx)
	IMPLICIT NONE

	INTEGER	i,j,mi,mj,i1,j1
	REAL*8 c2,dx,x,y,ysq
	COMPLEX*16 u(mi,mj),c1dxdy

	i1=mi/2+1
	j1=mj/2+1
	y =-dx*DBLE(j1)
	DO j=1,mj
	  y = y+dx
	  ysq=y*y
	  x =-dx*DBLE(i1)
	  DO i=1,mi
	    x =	x+dx 
	    u(i,j)=c1dxdy*u(i,j)*CDEXP(DCMPLX(0d0,c2*(x*x+ysq)))
	  END DO
	END DO

	RETURN
	END SUBROUTINE Applyfac2
	
C***********************************************************************

       SUBROUTINE RevApplyfac2(u,mi,mj,c1dxdy,c2,dx)
        IMPLICIT NONE

        INTEGER i,j,mi,mj,i1,j1
        REAL*8 c2,dx,x,y,ysq
        COMPLEX*16 u(mi,mj),c1dxdy

        i1=mi/2+1
        j1=mj/2+1
        y =-dx*DBLE(j1)
        DO j=1,mj
          y = y+dx
          ysq=y*y
          x =-dx*DBLE(i1)
          DO i=1,mi
            x = x+dx
            u(i,j)=u(i,j)/(c1dxdy*CDEXP(DCMPLX(0d0,c2*(x*x+ysq))))
          END DO
        END DO

        RETURN
        END SUBROUTINE RevApplyfac2

C***********************************************************************
c______________________________ DFT ROUTINES ___________________________
c
      SUBROUTINE DFT2D(UIN,NXIN,NYIN,DX,DY,X1,Y1,IDR,
     &                 UTR,NFXOUT,NFYOUT,DFX,DFY) 
c
c     FINITE FOURIER TRANSFORM OF COMPLEX-VALUED FUNCTION, UIN(X,Y).    
c     Input:                                                              
c      UIN(,)         complex-valued array to be transformed            
c      NXIN,NYIN      no. samples in x and in y, respectively           
c      DX,DY          sample intervals in x and y                       
c      X1,Y1          coordinates of first sample point (I=J=1)         
c      IDR            +1 for Fourier transform (Bracewell's definition) 
c                     -1  "  inverse     "                              
c      NFXOUT,NFYOUT  total number of transform values to be calculated 
c      DFX,DFY        sample spacing in transform space                 
c     Output:                                                            
c      UTR(i,j)       Fourier transform/(DX*DY); i=1,NFXOUT, j=1,NFYOUT 
c                     Zero frequency at i = (NFXOUT/2)+1                
c                                       j = (NFYOUT/2)+1                
c    Note that if NYIN is greater than 1, DFT does a TWO-DIMENSIONAL    
c    transform; otherwise, it does a ONE-DIMENSIONAL transform.         
c    In the one-dimensional case, the Y-parameters should be given      
c    dummy values (only NYIN is used).                                  
c    All input quantities are unchanged upon exiting the routine.       
c    The arrays UIN and UTR(,) may be the same array (to save storage in
c    the calling routine) when NFXOUT (NFYOUT) is no larger than NXIN 
c    (NXOUT).
c    Subroutine DSCTRN calculates sine and cosine transforms -- divided 
c    DX or DY -- of single-variable input function (array TEMP);        
c    cosine transform is returned in array A, i*(sine transform) in B.  
c
c                      Robert Benson
c                      Lockheed Martin Missiles & Space
c                      Advanced Technology Center
c                      3251 Hanover Street, B/201, O/L9-23
c                      Palo Alto  CA 94304
c                                                                       
c                      21 June 1989                                     
c                      21 Dec '95  Modified for use in Gort:            
c                                  double precision                     
c                      10 Sept'96  Added UTR as separate output and     
c                                  rearranged argument list             
c
      IMPLICIT REAL*8(A - H, O - Z) 
      IMPLICIT INTEGER (I - N) 
      INTEGER IDR 
      COMPLEX*16 UIN(*) 
      COMPLEX*16 UTR(*) 
c
c     Specify storage for work arrays - - - - - - - - - - - -           
      PARAMETER (NMAX = 2049, NHALF = (NMAX/2) + 1) 
c
c     Local variables:                                                  
      COMPLEX*16 A(NHALF), B(NHALF), TEMP(NMAX) 
      LOGICAL EVEN 
c
      NFXHAF = NFXOUT/2 
      NFYHAF = NFYOUT/2 
      NFX = NFXHAF + 1 
      NFY = NFYHAF + 1 

      IF ((NFX.GT.NHALF) .OR. (NFY.GT.NHALF) .OR.
     &    (NXIN.GT.NMAX) .OR. (NYIN.GT.NMAX))
     &    WRITE(*,*)'DFT2D: insufficient storage'
c                       
      IF (NYIN.GT.1) THEN 
c . . . Do transform in y-direction . . .                               

        IF ( 2*NFYHAF .EQ. NFYOUT) THEN 
          EVEN  = .TRUE. 
          JLAST = NFYHAF 
        ELSE 
          EVEN  = .FALSE. 
          JLAST = NFY 
        ENDIF 

        DO I = 1, NXIN 
           DO J = 1, NYIN 
              TEMP(J) = UIN(I+(J-1)*NXIN) 
           END DO 
           IF ( IDR.LT.0 ) THEN 
              DO J = 1, NYIN 
                 TEMP(J) = CONJG( TEMP(J) ) 
              END DO 
           ENDIF 

          CALL DSCTRN( TEMP, NYIN, Y1, DY, NFY, DFY, A, B ) 

c         Arrange transform with zero frequency                         
c         at J = NFY = NFYHAF + 1                                       

          UTR(I+NFYHAF*NFXOUT) = A(1) 
          DO J = 2, JLAST 
            UTR(I+(NFYHAF+(J-1))*NFXOUT) = A(J) - B(J) 
            UTR(I+(NFYHAF-(J-1))*NFXOUT) = A(J) + B(J) 
          END DO 
          IF ( EVEN ) UTR(I) = A(NFY) + B(NFY) 
c
        END DO 
        JMAX = NFYOUT 
      ELSE 
        JMAX = 1 
        DO I = 1, NXIN
           UTR(I) = UIN(I)
        END DO
      ENDIF 
c . . . Do transform in x-direction . . .                               

      IF ( 2*NFXHAF .EQ. NFXOUT) THEN 
        EVEN  = .TRUE. 
        ILAST = NFXHAF 
      ELSE 
        EVEN  = .FALSE. 
        ILAST = NFX 
      ENDIF 

      DO J = 1, JMAX 
        DO I = 1, NXIN 
          TEMP(I) = UTR(I+(J-1)*NFXOUT) 
        END DO 
        IF ( IDR.LT.0 ) THEN 
          DO I = 1, NXIN 
            TEMP(I) = CONJG( TEMP(I) ) 
          END DO 
        ENDIF 

        CALL DSCTRN( TEMP, NXIN, X1, DX, NFX, DFX, A, B ) 

c       Arrange transform with zero frequency                           
c       at I = NFX = NFXHAF + 1                                         

        UTR(NFX+(J-1)*NFXOUT) = A(1) 
        DO I = 2, ILAST 
          UTR(NFX+(I-1)+(J-1)*NFXOUT) = A(I) - B(I) 
          UTR(NFX-(I-1)+(J-1)*NFXOUT) = A(I) + B(I) 
        END DO 
        IF ( EVEN ) UTR(1+(J-1)*NFXOUT) = A(NFX) + B(NFX) 
      END DO 
c
      IF ( IDR.LT.0 ) THEN 
        IF ( NYIN.GT.1 ) THEN 
          NTOTAL = NFXOUT*NFYOUT 
        ELSE 
          NTOTAL = NFXOUT 
        ENDIF 
c
        DO IJ = 1, NTOTAL 
          UTR(IJ) = CONJG( UTR(IJ) ) 
        END DO 
      ENDIF 
c
      RETURN 
      END SUBROUTINE DFT2D
 
c---------------------------------------------------------------------

      SUBROUTINE TRIGS( FDT,NOLD,NFRQ,COSINE,SINE ) 
c
c     Called by DSCTRN to calculate trigonometric                       
c     functions for shift factors.                                      
c        ... 21 June 1989 ...                                           
c
      IMPLICIT REAL*8(A - H, O - Z) 
      IMPLICIT INTEGER (I - N) 
      REAL*8 COSINE(*), SINE(*) 
c
      twopi = 6.283185307179586476925286766559005768394D0
      COSINE(1) = 1.0D+0 
      SINE(1)   = 0.0D+0 
c
      TEMP = TWOPI*FDT 
      DO L = NOLD+1, NFRQ 
        ELm1=L-1
        THETA = TEMP*ELm1 
        COSINE(L) = COS( THETA ) 
        SINE(L)   = SIN( THETA ) 
      END DO 
c
      RETURN 
      END SUBROUTINE TRIGS

c-------------------------------------------------------------------

      SUBROUTINE TRIGM( FDT,NOLD,NFRQ,NCASE,COSPM1,SINE ) 
c
c     Called by DSCTRN to calculate trigonometric                       
c     functions for main frequencies:                                   
c                       2*[cos(theta)-1], L <= NCASE                    
c           COSPM1(L) =                                                 
c                       2*[cos(theta)+1], L > NCASE                    
c     while                                                             
c           SINE(L)   =     sin(theta), all L-values                    

c           ... 21 June 1989 ...                                        

      IMPLICIT REAL*8(A - H, O - Z) 
      IMPLICIT INTEGER (I - N) 
      REAL*8 COSPM1(*), SINE(*) 
c
      pi = dacos(-1d0)
      twopi = 2d0*pi
      TEMP = ABS( 0.25D+0 / FDT ) 
      NCASE = MIN( NFRQ, INT( TEMP ) ) 
c
      COSPM1(1) = 0.0D+0 
      SINE(1)   = 0.0D+0 
c
      TEMP = TWOPI*FDT 
      DO L = NOLD+1, NCASE 
        ELm1 = L-1
        THETA = TEMP*ELm1 
        COSPM1(L) = - 4.0D+0*( SIN( 0.5D+0*THETA ) )**2 
        SINE(L)   = SIN( THETA ) 
      END DO 
      DO L = NCASE+1, NFRQ 
        ELm1 = L-1
        THETA = TEMP*ELm1 
        COSPM1(L) = 4.0D+0*( COS( 0.5D+0*THETA ) )**2 
        SINE(L)   = SIN( THETA ) 
      END DO 
c
      RETURN 
      END SUBROUTINE TRIGM                         
 
c---------------------------------------------------------------------

      SUBROUTINE DSCTRN( DATA,NT,T1,DT,NFRQS,FZRO,AA,BB ) 
c
c     Discrete sine and cosine transforms of complex input array.       
c                     Robert Benson                                     
c                 ... 21 June 1989 ...                                  
c     Fixed bug in initializing shift factors ... 18 Dec '90            

c  Input:                                                               
c     DATA   Array containing function values (complex-valued)          
c     NT     Number of sample points                                    
c     T1     Coordinate value for first sample point                    
c     DT     Sample interval                                            
c     NFRQS  Number of frequencies at which transform to be evaluated   
c     FZRO   Frequency interval at which integrals are evaluated        
c  Output (for L = 1,2,...,NFRQS):                                      
c     AA(L) =     [SUM OF DATA(J)*COS(2*PI*(L-1)*FZRO*T)]               
c     BB(L) = i * [ "  "    "    *SIN(.................)]               
c       Summation limits correspond to spatial coordinate values        
c         TLWR = T1 ................(J=1)                               
c         TUPR = TLWR + (NT-1)*DT...(J=NT)                              
c       (Resembles trapezoid rule as approximation to integrals for
c       sine and cosine transforms of data)                            
c     * DSCTRN uses Goertzel-Reinsch algorithm -- see especially        
c       Stoer & Bulirsch, "Intro. to Numerical Analysis,"               
c                       Springer-Verlag, N.Y., 1980, sec. 2.3.3.        
c     * For Goertzel algorithm, see:                                    
c       Burrus & Parks, "DFT/FFT and Convolution Algorithms,"           
c                       John Wiley & Sons, N.Y., 1985.                  
c       Dick & Wertz, "Analog and Digital Computation of Fourier        
c                     Series and Integrals," IEEE Trans. Electronic     
c                     Computers, vol. EC-16, 8-13, 1967                 
c       Acton,F.S., "Numerical Methods that Work,"                      
c                    Harper & Row, N.Y., 1970, chapter 9                
c
      IMPLICIT REAL*8(A - H, O - Z) 
      IMPLICIT INTEGER (I - N) 
      COMPLEX*16 AA(*), BB(*), DATA(*) 

c     Specify storage for sines and cosines 
      PARAMETER (NMAX = 2049, NHALF = (NMAX/2) + 1) 
      PARAMETER (CRIT = 1.0D-6) 
c
      REAL*8 AR, AI, BR, BI 
      REAL*8 C(NHALF), S(NHALF), CS(NHALF), SS(NHALF) 
      SAVE FZRODT, FZROT1, NCASE1, NMAIN, NSHIFT, C, S, CS, SS 
      DATA FZRODT, FZROT1, NCASE1, NMAIN, NSHIFT / 2 * 0.0, 3 * 0 / 
c
c     If necessary, initialize trigonometric functions:                 
c . . . .                 Main frequencies                              

      IF ( ABS( FZRODT - (FZRO*DT) ).GT.CRIT*ABS( FZRODT ) ) THEN 
        FZRODT = FZRO*DT 
        CALL TRIGM( FZRODT, 1, NFRQS, NCASE1, C, S ) 
        NMAIN = NFRQS 
      ELSEIF ( NMAIN.LT.NFRQS ) THEN 
        CALL TRIGM( FZRODT, NMAIN, NFRQS, NCASE1, C, S ) 
        NMAIN = NFRQS 
      ENDIF 
c . . . .                 Shift factors                                 

      IF ( ABS( FZROT1 - (FZRO*T1) ).GT.CRIT*ABS( FZROT1 ) ) THEN 
        FZROT1 = FZRO*T1 
        CALL TRIGS( FZROT1, 1, NFRQS, CS, SS ) 
        NSHIFT = NFRQS 
      ELSEIF ( NSHIFT.LT.NMAIN ) THEN 
        CALL TRIGS( FZROT1, NSHIFT, NMAIN, CS, SS ) 
        NSHIFT = NMAIN 
      ENDIF 

c     In case data are non-zero at edge of grid, invoke                     
c     trapezoid rule for approximating sine & cosine transforms.        

      DATA(1) = DCMPLX(0.5D+0*DREAL(DATA(1)), 
     &                0.5D+0*DIMAG(DATA(1)) )
      DATA(NT) = DCMPLX(0.5D+0*DREAL(DATA(NT)), 
     &                 0.5D+0*DIMAG( DATA(NT) ) )
c
c     Calculate sine and cosine transforms using                        
c     Goertzel-Reinsch algorithm.                                       
c     Main loops: zero frequency is special (simple) case               

      AA(1) = DCMPLX( 0.0D+0, 0.0D+0 ) 
      BB(1) = DCMPLX( 0.0D+0, 0.0D+0 ) 
      DO K = 1, NT 
        AA(1) = AA(1) + DATA(K) 
      END DO 

c     Goertzel-Reinsch recursion for non-zero frequencies               

      DO L = 2, NCASE1 
        HALFCL = 0.5D+0 * C(L) 
        U2R = 0.0D+0 
        U2I = 0.0D+0 
        DUR = DREAL( DATA(NT) ) 
        DUI = DIMAG( DATA(NT) ) 
        DO K = NT - 1, 1, - 1 
          U1R = DUR + U2R 
          U1I = DUI + U2I 
          DUR = DREAL( DATA(K) )  + C(L) * U1R + DUR 
          DUI = DIMAG( DATA(K) ) + C(L) * U1I + DUI 
          U2R = U1R 
          U2I = U1I 
        END DO 

        AR = DUR - HALFCL * U1R 
        AI = DUI - HALFCL * U1I 
        BR = S(L) * U1R 
        BI = S(L) * U1I 

        AA(L) = DCMPLX( AR, AI ) 
        BB(L) = DCMPLX( BR, BI ) 
      END DO 
c
      DO L = NCASE1 + 1, NFRQS 
        HALFCL = 0.5D+0 * C(L) 
        U2R = 0.0D+0 
        U2I = 0.0D+0 
        DUR = DREAL( DATA(NT) ) 
        DUI = DIMAG( DATA(NT) ) 

        DO K = NT - 1, 1, - 1 
          U1R = DUR - U2R 
          U1I = DUI - U2I 
          DUR = DREAL( DATA(K) )  + C(L) * U1R - DUR 
          DUI = DIMAG( DATA(K) ) + C(L) * U1I - DUI 
          U2R = U1R 
          U2I = U1I 
        END DO 

        AR = DUR - HALFCL * U1R 
        AI = DUI - HALFCL * U1I 
        BR = S(L) * U1R 
        BI = S(L) * U1I 

        AA(L) = DCMPLX( AR, AI ) 
        BB(L) = DCMPLX( BR, BI ) 
      END DO 
c
c     Shift As and Bs back to correspond to TLWR to TUPR interval       
c     and multiply BB by SQRT(-1).  L=1 values are unchanged.           
c
      DO L = 2, NFRQS 
        AR = DREAL( AA(L) ) 
        AI = DIMAG( AA(L) ) 
        BR = REAL( BB(L) ) 
        BI = DIMAG( BB(L) ) 
        AA(L)=DCMPLX((CS(L)*AR-SS(L)*BR),
     &              (CS(L)*AI-SS(L)*BI))
        BB(L)=DCMPLX(-(CS(L)*BI+SS(L)*AI),
     &              (CS(L)*BR+SS(L)*AR))
      END DO 
c
      RETURN 
      END SUBROUTINE DSCTRN
