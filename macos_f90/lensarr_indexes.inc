C
C lensarr_indexes.inc
C Included by utilsub.F::InitLensletIndexAndCtr
C
        !HexLensletIdx(1:183)=(/ 
        HexLensletIdx(1:273)=(/ 
     &    0,0,0, 1,2,1, -1,1,2, -2,-1,1, -1,-2,-1, 1,-1,-2,
     &    2,1,-1, 2,4,2, 0,3,3, -2,2,4, -3,0,3, -4,-2,2, 
     &    -3,-3,0, -2,-4,-2, 0,-3,-3, 2,-2,-4, 3,0,-3, 4,2,-2,
     &    3,3,0,  ! end of 3rd ring
     &    3,6,3, 1,5,4, -1,4,5, -3,3,6, -4,1,5, -5,-1,4,
     &    -6,-3,3, -5,-4,1, -4,-5,-1, -3,-6,-3, -1,-5,-4,
     &    1,-4,-5, 3,-3,-6, 4,-1,-5, 5,1,-4,6,3,-3, 5,4,-1,
     &    4,5,1, ! end of 4th ring
     &    4,8,4, 2,7,5, 0,6,6, -2,5,7, -4,4,8, -5,2,7,
     &    -6,0,6, -7,-2,5, -8,-4,4, -7,-5,2, -6,-6,0, -5,-7,-2,
     &    -4,-8,-4, -2,-7,-5, 0,-6,-6, 2,-5,-7, 4,-4,-8, 5,-2,-7,
     &    6,0,-6, 7,2,-5, 8,4,-4, 7,5,-2, 6,6,0, 5,7,2, ! end of 5th ring
     &    5,10,5, 3,9,6, 1,8,7, -1,7,8, -3,6,9, -5,5,10, -6,3,9,
     &    -7,1,8, -8,-1,7, -9,-3,6, -10,-5,5, -9,-6,3, -8,-7,1, 
     &    -7,-8,-1, -6,-9,-3, -5,-10,-5, -3,-9,-6, -1,-8,-7, 1,-7,-8,
     &    3,-6,-9, 5,-5,-10, 6,-3,-9, 7,-1,-8, 8,1,-7, 9,3,-6, 10,5,-5,
     &    9,6,-3, 8,7,-1, 7,8,1, 6,9,3  ! end of 6th ring
     &    /)

C

#if 0
        RecLensletIdx(1:50)=(/ ! 5x5 square lens array
     &    0,0, 
     &    2,0, 2,2, 0,2, -2,2, -2,0, -2,-2, 0,-2, 2,-2,
     &    4,0, 4,2, 4,4, 2,4, 0,4, -2,4, -4,4, -4,2, 
     &     -4,0, -4,-2, -4,-4, -2,-4, 0,-4, 2,-4, 4,-4, 4,-2
     &    /)
	nRecLenlets=25
#endif


#if 1
!--> square lens array sized from mLenslet (SAVE-audit fix,
!    2026-07-04).  This block used to hardwire ndim=107 -> 11449
!    lenslets, but RecLensletIdx is allocated mLenslet*2 (=500) and
!    LensletStrArray mLenslet (=250) -- the fill loop below and the
!    copy loop in InitLensletIndexAndCtr stomped memory the moment
!    any Rx used LensArrayIndRef= or XYIndRefFile= (zero corpus
!    coverage; broken since the 107x107 table replaced the 5x5).
!    ndim = largest odd integer with ndim**2 <= mLenslet.
	ndim=INT(SQRT(DBLE(mLenslet)))
	If (MOD(ndim,2).EQ.0) ndim=ndim-1
	nn=(ndim-1)/2; k=1
	Do r=-nn,nn
	  Do c=-nn,nn
	    RecLensletIdx(k)=2*r;
	    RecLensletIdx(k+1)=2*c;
	    k=k+2
	  End Do
	End Do
	nRecLenlets=ndim*ndim
!<--
#endif

      
