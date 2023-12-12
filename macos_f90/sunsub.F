C***********************************************************************
C	Begin file sunsub.F
C    +----------------------------------------------------------------+
C    |  Copyright (C) 1995-2007, California Institute of Technology.  |
C    |  U.S. Government Sponsorship Is Acknowledged.                  |
C    +----------------------------------------------------------------+
C***********************************************************************


#include "realtype.h"

C***********************************************************************

        SUBROUTINE RTIME(RSECS)
        SREAL RSECS

        ! SECNDS is a Fortran library function
        RSECS=SECNDS(0e0)
        RETURN
        END SUBROUTINE RTIME


C***********************************************************************
!	DSWAP2: exchange grid blocks between 1 and 3 quarters,
!               2 and 4 quaters. It's called before and after
!		DFOURN - the FFT-based diffraction. -jzlou
C***********************************************************************

	SUBROUTINE DSWAP2(OUT,IN,N)
	COMPLEX*16 IN(N,N),OUT(N,N),CD1,CD2,CD3,CD4
	INTEGER N,ND2,J,JP,I,IP
	ND2=N/2
	DO 1 J=1,ND2
	  JP=J+ND2
	  DO 1 I=1,ND2
	    IP=I+ND2
	    CD1=IN(IP,JP)
	    CD2=IN(I,JP)
	    CD3=IN(IP,J)
	    CD4=IN(I,J)
	    OUT(I,J)=CD1
	    OUT(IP,J)=CD2
	    OUT(I,JP)=CD3
	    OUT(IP,JP)=CD4
    1	CONTINUE
	RETURN
	END SUBROUTINE DSWAP2


C***********************************************************************
!	DSWAP3: "In-Place Swap" version of DSWAP2
!       -jzlou
C***********************************************************************

	SUBROUTINE DSWAP3(IN,N)
	COMPLEX*16 IN(N,N),CD1,CD2,CD3,CD4
	INTEGER N,ND2,J,JP,I,IP
	ND2=N/2
	DO 1 J=1,ND2
	  JP=J+ND2
	  DO 1 I=1,ND2
	    IP=I+ND2
	    CD1=IN(IP,JP)
	    CD2=IN(I,JP)
	    CD3=IN(IP,J)
	    CD4=IN(I,J)
	    IN(I,J)=CD1
	    IN(IP,J)=CD2
	    IN(I,JP)=CD3
	    IN(IP,JP)=CD4
    1	CONTINUE
	RETURN
	END SUBROUTINE DSWAP3

C***********************************************************************

	SUBROUTINE CRSWAP2(OUT,IN,N)
	SREAL OUT(N,N)
	COMPLEX*16 IN(N,N),CD1,CD2,CD3,CD4
	INTEGER N,ND2,J,JP,I,IP,N2
	N2=N*N
	ND2=N/2
	DO J=1,ND2
	  JP=J+ND2
	  DO I=1,ND2
	    IP=I+ND2
	    CD1=IN(IP,JP)
	    CD2=IN(I,JP)
	    CD3=IN(IP,J)
	    CD4=IN(I,J)
	    OUT(I,J)=CD1/N2
	    OUT(IP,J)=CD2/N2
	    OUT(I,JP)=CD3/N2
	    OUT(IP,JP)=CD4/N2
    	  END DO
    	END DO
	RETURN
	END

C***********************************************************************

#undef NR_FFT
#undef SUN_FFT
#undef FFTW_FFT

! If SUNOS is defined, use SUN_FFT by default.
! When SUN_FFT and USE_FFTW are both defined, the latter
! is used.
#ifdef SUNOS
# if defined FFTW_LIB
#   define FFTW_FFT
# else
#    define SUN_FFT
# endif
#endif

#ifdef LNXOS
# if defined FFTW_LIB
#   define FFTW_FFT
# else
#   if defined INTEL_FFT_LIB
#     define INTEL_FFT
#   else
#     define NR_FFT
#   endif
# endif
#endif

#ifdef MSWIN
#   if defined INTEL_FFT_LIB
#     define INTEL_FFT
#   else
#     define NR_FFT
#   endif
#endif


! Don't use FFTW_FFT with Matlab 6.x, not compatible!  jzlou (8-26-03)


#ifdef NR_FFT
C
C  A FFT routine from Numerical Recipes, adapted from
C  /home/comp/macos-v2.86/sunsubDFOURN.f
C
      SUBROUTINE DFOURN(CDATA,NN,NDIM,ISIGN)
      Implicit none
      INTEGER :: NDIM,ISIGN,NN(NDIM)
      COMPLEX*16 :: CDATA(NN(1),NN(2))
      !Real*8 :: CDATA(1)
      REAL*8 :: WR,WI,WPR,WPI,WTEMP,THETA,TEMPI,TEMPR
      Real*8, dimension(:), allocatable, save :: DATA
      Integer :: NTOT,N,NREM,IP1,IP2,IP3,I2REV,I3REV,
     &           I1,I2,I3,IBIT,IFP1,IFP2,K1,K2,IDIM,NPREV
      Integer :: SzData,m_err_rtn,i,j,k
      Logical, save :: first_entry=.true.

      If (first_entry) Then
	SzData=2
	Do i=1,NDIM
	  SzData=SzData*NN(i)
	End Do
	if (allocated(DATA)) then
	  deallocate(DATA,stat=m_err_rtn)
	  if (m_err_rtn /= 0)
     &    call macos_memory_failure
     &               ('DFOURN: deallocate failed!')
	end if
	allocate(DATA(SzData),stat=m_err_rtn)
	if (m_err_rtn /= 0)
     &    call macos_memory_failure
     &               ('DFOURN: allocate failed!')
	first_entry=.false.
      End If

      k=1; DATA=0d0
      Do j=1,NN(2)
        Do i=1,NN(1)
	  DATA(k)=DREAL(CDATA(i,j))
	  DATA(k+1)=DIMAG(CDATA(i,j))
          k=k+2
        End Do
      End Do

#if 0
      If (NDIM==2) Then
	k=1
	Do j=1,NN(2)
          DATA(k:2:k+(2*NN(1)-1)-1)=DREAL(CDATA(1:NN(1),j))
          DATA(k+1:2:k+1+(2*NN(1))-1)=DIMAG(CDATA(1:NN(1),j))
	  k=k+NN(1)*2
	End Do
      Else
	WRITE(*,*)'**DFOURN: NDIM MUST BE 2, quiting ...'
	stop
      End If
#endif

      NTOT=1
      DO 11 IDIM=1,NDIM
        NTOT=NTOT*NN(IDIM)
11    CONTINUE
      NPREV=1
      DO 18 IDIM=1,NDIM
        N=NN(IDIM)
        NREM=NTOT/(N*NPREV)
        IP1=2*NPREV
        IP2=IP1*N
        IP3=IP2*NREM
        I2REV=1
        DO 14 I2=1,IP2,IP1
          IF(I2.LT.I2REV)THEN
            DO 13 I1=I2,I2+IP1-2,2
              DO 12 I3=I1,IP3,IP2
                I3REV=I2REV+I3-I2
                TEMPR=DATA(I3)
                TEMPI=DATA(I3+1)
                DATA(I3)=DATA(I3REV)
                DATA(I3+1)=DATA(I3REV+1)
                DATA(I3REV)=TEMPR
                DATA(I3REV+1)=TEMPI
12            CONTINUE
13          CONTINUE
          ENDIF
          IBIT=IP2/2
1         IF ((IBIT.GE.IP1).AND.(I2REV.GT.IBIT)) THEN
            I2REV=I2REV-IBIT
            IBIT=IBIT/2
          GO TO 1
          ENDIF
          I2REV=I2REV+IBIT
14      CONTINUE
        IFP1=IP1
2       IF(IFP1.LT.IP2)THEN
          IFP2=2*IFP1
c           THETA=QFLOAT(ISIGN)*
c      &            6.283185307179586476925286766559005768394Q0
c      &            /QFLOAT(IFP2/IP1)
c           WPR=-2.Q0*QSIN(0.5Q0*THETA)**2
c           WPI=QSIN(THETA)
c           WR=1.Q0
c           WI=0.Q0
          THETA=ISIGN*6.28318530717959D0/(IFP2/IP1)
          WPR=-2.D0*SIN(0.5D0*THETA)**2
          WPI=SIN(THETA)
          WR=1.D0
          WI=0.D0
          DO 17 I3=1,IFP1,IP1
            DO 16 I1=I3,I3+IP1-2,2
              DO 15 I2=I1,IP3,IFP2
                K1=I2
                K2=K1+IFP1
                TEMPR=DBLE(WR)*DATA(K2)-DBLE(WI)*DATA(K2+1)
                TEMPI=DBLE(WR)*DATA(K2+1)+DBLE(WI)*DATA(K2)
                DATA(K2)=DATA(K1)-TEMPR
                DATA(K2+1)=DATA(K1+1)-TEMPI
                DATA(K1)=DATA(K1)+TEMPR
                DATA(K1+1)=DATA(K1+1)+TEMPI
15            CONTINUE
16          CONTINUE
            WTEMP=WR
            WR=WR*WPR-WI*WPI+WR
            WI=WI*WPR+WTEMP*WPI+WI
17        CONTINUE
          IFP1=IFP2
        GO TO 2
        ENDIF
        NPREV=N*NPREV
18    CONTINUE

      k=1
      Do j=1,NN(2)
        Do i=1,NN(1)
          CDATA(i,j)=DCMPLX(DATA(k),DATA(k+1))
	  k=k+2
	End Do
      End Do
      RETURN
      END SUBROUTINE DFOURN
#endif


#ifdef FFTW_FFT
!
! This is the wrapper routine for orininal FFTW
! To use this routine, must decomment the relevant FFT calls starting
! with '!!'
!
      SUBROUTINE DFOURN(DATA,NN,NDIM,ISIGN)
      ! Need to access 'macos_model_size' variable
      ! defined in param_mod.
      use param_mod

      IMPLICIT NONE
      integer, save :: FPLAN = 0, IPLAN = 0
      INTEGER NDIM, ISIGN ,NN(NDIM)
      COMPLEX*16 DATA(NN(1),NN(2))
      Integer, save :: crnt_model_size = -1

      if (crnt_model_size < 0 .or.
     &    macos_model_size /= crnt_model_size) then
        ! first entry or model size has changed
        crnt_model_size = macos_model_size
        FPLAN = 0; IPLAN = 0
      end if

      IF (NDIM.NE.2) THEN
        PRINT *,'FFTW-FFT ERROR:  not two dimensions'
        STOP
      ENDIF

      IF (ISIGN.EQ.1) THEN
        IF (FPLAN.EQ.0) THEN
          print *,'********** INITIALIZING FORWARD FFT *************'
          call fftw2d_f77_create_plan(FPLAN,NN(1),NN(2),-1,1+8)
        ENDIF
         call fftwnd_f77(FPLAN,1,DATA,1,0,0,0,0)
      ELSEIF (ISIGN.EQ.-1) THEN
        IF (IPLAN.EQ.0) THEN
          print *,'********** INITIALIZING INVERSE FFT *************'
          call fftw2d_f77_create_plan(IPLAN,NN(1),NN(2),1,1+8)
        ENDIF
        call fftwnd_f77(IPLAN,1,DATA,1,0,0,0,0)
      ELSE
        PRINT *,'FFTW-FFT ERROR:  ISIGN not correct'
        STOP
      ENDIF

      RETURN
      END SUBROUTINE DFOURN
#endif


#ifdef INTEL_FFT
      SUBROUTINE DFOURN(DATA,NN,NDIM,ISIGN)
      IMPLICIT NONE

      INTEGER NDIM, ISIGN ,NN(NDIM)
      COMPLEX*16 DATA(NN(1),NN(2))

      IF ((ISIGN /= 1 .and. ISIGN /= -1)
     &    .or. NN(1) .le. 0 .or. NN(2) .le. 0) THEN
        PRINT *,'INTEL-FFT ERROR:  invalid ISIGN or NN value'
        STOP
      ENDIF

      print*,'=== Calling Intel FFT: ZFFT2D ==='
      CALL ZFFT2D(DATA,NN(1),NN(2),ISIGN)
      IF (ISIGN.EQ.1) DATA = DATA*NN(1)*NN(2)
      END SUBROUTINE DFOURN
#endif


#ifdef SUN_FFT
!
! This version uses Sun's Performance Library
!
      SUBROUTINE DFOURN(DATA,NN,NDIM,ISIGN)
      use sunperf

      ! Need to access 'macos_model_size' variable
      ! defined in param_mod.
      use param_mod

      IMPLICIT NONE

      INTEGER NDIM, ISIGN ,NN(NDIM)
      COMPLEX*16 DATA(NN(1),NN(2))

      Integer, save :: FPLAN = 0
      Integer, save :: crnt_model_size = -1
      Real*8, dimension(:), allocatable, save :: WS
      Integer, SAVE :: LWS, m_err_fft

      IF (NDIM.NE.2) THEN
        PRINT *,'FFT ERROR:  not two dimensions'
        STOP
      ENDIF

      if (crnt_model_size < 0 .or.
     &    macos_model_size /= crnt_model_size) then
        ! first entry or model size has changed
        crnt_model_size = macos_model_size
        FPLAN = 0
	if (allocated(WS)) then
	  deallocate(WS, stat=m_err_fft)
	  if (m_err_fft /= 0)
     &	    call macos_memory_failure
     &                ('dfourn: deallocate failed!')
	end if
	LWS = 4*(NN(1)+NN(2))+30
	allocate(WS(LWS), stat=m_err_fft)
	if (m_err_fft /= 0)
     &      call macos_memory_failure
     &                ('dfourn: allocate failed!')
	print*,'***** WS for FFT allocated *****'
      end if

      if (FPLAN==0) then
	print *,'********** INITIALIZING SunPerf FFT *************'
	Call ZFFT2I(NN(1), NN(2), WS)
	FPLAN = 1
      end if

      IF (ISIGN.EQ.1) THEN
	Call ZFFT2F(NN(1), NN(2), DATA, NN(1), WS, LWS)
      ELSE IF (ISIGN.EQ.-1) THEN
	Call ZFFT2B(NN(1), NN(2), DATA, NN(1), WS, LWS)
      END IF
      END SUBROUTINE DFOURN
#endif


C***********************************************************************

	SUBROUTINE RSWAP2(OUT,IN,N)
	COMPLEX*8 IN(N,N),OUT(N,N)
	ND2=N/2
	DO 1 J=1,ND2
	  JP=J+ND2
	  DO 1 I=1,ND2
	    IP=I+ND2
	    OUT(I,J)=IN(IP,JP)
	    OUT(IP,J)=IN(I,JP)
	    OUT(I,JP)=IN(IP,J)
	    OUT(IP,JP)=IN(I,J)
    1	CONTINUE
	RETURN
	END SUBROUTINE RSWAP2

C***********************************************************************

	SUBROUTINE RSWAP3(IN,N)
	COMPLEX*8 IN(N,N),CD1,CD2,CD3,CD4
	INTEGER N,ND2,J,JP,I,IP
	ND2=N/2
	DO 1 J=1,ND2
	  JP=J+ND2
	  DO 1 I=1,ND2
	    IP=I+ND2
	    CD1=IN(IP,JP)
	    CD2=IN(I,JP)
	    CD3=IN(IP,J)
	    CD4=IN(I,J)
	    IN(I,J)=CD1
	    IN(IP,J)=CD2
	    IN(I,JP)=CD3
	    IN(IP,JP)=CD4
    1	CONTINUE
	RETURN
	END SUBROUTINE RSWAP3

C***********************************************************************

      SUBROUTINE RFOURN(DATA,NN,NDIM,ISIGN)
      USE Constants, ONLY: TWOPI
      REAL*8 WR,WI,WPR,WPI,WTEMP,THETA
      COMPLEX*8 DATA
      DIMENSION NN(NDIM),DATA(1)
      NTOT=1
      DO 11 IDIM=1,NDIM
        NTOT=NTOT*NN(IDIM)
11    CONTINUE
      NPREV=1
      DO 18 IDIM=1,NDIM
        N=NN(IDIM)
        NREM=NTOT/(N*NPREV)
        IP1=2*NPREV
        IP2=IP1*N
        IP3=IP2*NREM
        I2REV=1
        DO 14 I2=1,IP2,IP1
          IF(I2.LT.I2REV)THEN
            DO 13 I1=I2,I2+IP1-2,2
              DO 12 I3=I1,IP3,IP2
                I3REV=I2REV+I3-I2
                TEMPR=DATA(I3)
                TEMPI=DATA(I3+1)
                DATA(I3)=DATA(I3REV)
                DATA(I3+1)=DATA(I3REV+1)
                DATA(I3REV)=TEMPR
                DATA(I3REV+1)=TEMPI
12            CONTINUE
13          CONTINUE
          ENDIF
          IBIT=IP2/2
1         IF ((IBIT.GE.IP1).AND.(I2REV.GT.IBIT)) THEN
            I2REV=I2REV-IBIT
            IBIT=IBIT/2
          GO TO 1
          ENDIF
          I2REV=I2REV+IBIT
14      CONTINUE
        IFP1=IP1
2       IF(IFP1.LT.IP2)THEN
          IFP2=2*IFP1
          THETA=ISIGN*TWOPI/(IFP2/IP1)
          WPR=-2.D0*SIN(0.5D0*THETA)**2
          WPI=SIN(THETA)
          WR=1.D0
          WI=0.D0
          DO 17 I3=1,IFP1,IP1
            DO 16 I1=I3,I3+IP1-2,2
              DO 15 I2=I1,IP3,IFP2
                K1=I2
                K2=K1+IFP1
                TEMPR=SNGL(WR)*DATA(K2)-SNGL(WI)*DATA(K2+1)
                TEMPI=SNGL(WR)*DATA(K2+1)+SNGL(WI)*DATA(K2)
                DATA(K2)=DATA(K1)-TEMPR
                DATA(K2+1)=DATA(K1+1)-TEMPI
                DATA(K1)=DATA(K1)+TEMPR
                DATA(K1+1)=DATA(K1+1)+TEMPI
15            CONTINUE
16          CONTINUE
            WTEMP=WR
            WR=WR*WPR-WI*WPI+WR
            WI=WI*WPR+WTEMP*WPI+WI
17        CONTINUE
          IFP1=IFP2
        GO TO 2
        ENDIF
        NPREV=N*NPREV
18    CONTINUE
      RETURN
      END SUBROUTINE RFOURN

C***********************************************************************

	SUBROUTINE SPOTOUT(XYMAT,NA,NC,NDA,NDC,iBeam,
     &	CXLABL,CYLABL,CTITLE,iPloTyp,image_file,iElt)

	use param_mod
	use cfiles_mod
#ifdef CSMACOS
        use smacosio_mod
#else
        use macosio_mod
#endif

	INTEGER ArchFlag,Sun,Mac,DEC
	PARAMETER(Sun=0,Mac=1,DEC=2,ArchFlag=Sun)

	LOGICAL exist
	CHARACTER(len=MacosCharLen) :: file,ANS
	CHARACTER*80 imag_file
	CHARACTER*82 SCR
	CHARACTER*(*) CXLABL,CYLABL,CTITLE,image_file
	INTEGER iElt,NA,NC,NDA,NDC,iPloTyp,L,ICLEN,RecordLength,iBeam
	REAL*8 XYMAT(NDA,NDC)

 102	FORMAT(' Text array is of dimension',i6,' by',i6)
 103	FORMAT(' Binary array is of dimension',i6,' by',i6)
 104	FORMAT(1P,512(1x,e13.6))
 105	FORMAT(1P,I4,A1,I4/E11.4,A1,E11.4)
 106	FORMAT(1P,512(E11.4,A1))
 108	FORMAT(' Matlab array is of dimension',i6,' by',i6)
 109	FORMAT(1P,2(1x,e15.8))

	print*,'** Beginning of SPOTOUT ...'
        print*,'ArchFlag =', ArchFlag
	print*,'NA, NC = ', NA,NC

C   Set record length depending on target machine

	IF (ArchFlag.EQ.Sun) THEN
	  !RecordLength=4*NC
	  RecordLength=8*NA
	ELSE IF (ArchFlag.EQ.Mac) THEN
	  RecordLength=NC
	ELSE IF (ArchFlag.EQ.DEC) THEN
	  RecordLength=NC
	END IF

C   Get output file name
	L=ICLEN(image_file)
	IF (iPloTyp.EQ.5) THEN
	  file=image_file(1:L)//'.txt'
	  image_file=file
	ELSE IF (iPloTyp.EQ.8) THEN
	  file=image_file(1:L)//'.fit'
	  image_file=file
	ELSE IF (iPloTyp.EQ.9) THEN
	  file=image_file(1:L)//'.mat'
	  image_file=file
	END IF

 1	CONTINUE
	  INQUIRE (FILE=image_file,EXIST=exist)
	  IF (exist) THEN
	    L=ICLEN(image_file)
	    SCR='Output file '//image_file(1:L)//' exists.  Overwrite?'
	    CALL CACCEPT(ANS,'YES',SCR)
	    IF ((ans(1:1).EQ.'n').OR.(ans(1:1).EQ.'N')) THEN
	      CALL CACCEPT(file,image_file,'Enter new file name:')
	      image_file=file
	      GO TO 1
	    ELSE
	      OPEN (UNIT=13,FILE=image_file,STATUS='OLD')
	      CLOSE (UNIT=13,STATUS='DELETE')
	    END IF
	  END IF

C  Text format file
	IF (iPloTyp.EQ.5) THEN
	  OPEN (UNIT=13,FILE=image_file,STATUS='NEW',FORM='FORMATTED')
	  WRITE(*,*) 'Writing ',CTITLE
	  WRITE(*,*) 'FORMATTED File=',image_file
	  WRITE(*,102)NA,NC
	  DO i=1,NA
	    WRITE(13,109)(XYMAT(i,j),j=1,NC)
	  END DO
	  CLOSE(13)

C  Binary format file
	ELSE IF (iPloTyp.EQ.6) THEN
	  OPEN (UNIT=13,FILE=image_file,ACCESS='DIRECT',
     &	        RECL=RecordLength,STATUS='NEW')
	  WRITE(*,*) 'Writing ',CTITLE
	  WRITE(*,*) 'FORMATTED File=',image_file
	  WRITE(*,103) NA,NC
	  DO 3 j=1,NC
	    WRITE(13,REC=j)(XYMAT(i,j),i=1,NA)
 3	  CONTINUE
	  CLOSE(13)

C  Fits format file
	ELSE IF (iPloTyp.EQ.8) THEN
	  WRITE(*,*)'Sorry, FITS not supported for SPOT or DRAW'

! C  Matlab format file
! 	ELSE IF (iPloTyp.EQ.9) THEN
! 	  WRITE(*,*) 'Writing ',CTITLE
! 	  WRITE(*,*) 'MATLAB File=',image_file
! 	  WRITE(*,108)NA,NC
!  	  CALL RWRITE_MAT_LOH(XYMAT,NA,NC,NDA,NDC,image_file,imag_file)

C  Unknown format
	ELSE
	  WRITE(*,*)' Invalid plot type specified: iPloTyp=',iPloTyp
	  RETURN
	END IF

C  End
 	RETURN
	END SUBROUTINE SPOTOUT

C***********************************************************************

	SUBROUTINE DRAWOUT(iPloTyp,CTITLE,nDrawElt,nDrawRay,DrawEltVec,
     &			   DrawRayVec,image_file,R1)
	use param_mod
#ifdef CSMACOS
        use smacosio_mod
#else
        use macosio_mod
#endif
	IMPLICIT NONE


	INTEGER ArchFlag,Sun,Mac,DEC
	PARAMETER(Sun=0,Mac=1,DEC=2,ArchFlag=Sun)

	LOGICAL exist

	CHARACTER(len=MacosCharLen) :: file,ANS
	CHARACTER*80 imag_file
	CHARACTER*82 SCR
	CHARACTER*(*) CTITLE,image_file

	INTEGER iPloTyp,i,j,k,L,ICLEN,RecordLength,nDrawRay,
     &	        nDrawElt(mDrawRay),DrawEltVec(mDrawElt,mDrawRay)

	SREAL DrawRayVec(2,mDrawElt,mDrawRay),R1(*)

 102	FORMAT(100(I4,','))
 103	FORMAT(1P,100(E11.4,','))
 104	FORMAT(1P,1E11.4)
 108	FORMAT(' Matlab array is of dimension',i6,' by',i6)

	WRITE(4,*)CTITLE

C   Set record length depending on target machine

	IF (ArchFlag.EQ.Sun) THEN
	  RecordLength=4
	ELSE IF (ArchFlag.EQ.Mac) THEN
	  RecordLength=1
	ELSE IF (ArchFlag.EQ.DEC) THEN
	  RecordLength=1
	END IF

C   Get output file name
	L=ICLEN(image_file)
	IF (iPloTyp.EQ.5) THEN
	  file=image_file(1:L)//'.txt'
	  image_file=file
	ELSE IF (iPloTyp.EQ.9) THEN
	  file=image_file(1:L)//'.mat'
	  image_file=file
	END IF

 1	CONTINUE
	  INQUIRE (FILE=image_file,EXIST=exist)
	  IF (exist) THEN
	    L=ICLEN(image_file)
	    SCR='Output file '//image_file(1:L)//' exists.  Overwrite?'
	    CALL CACCEPT(ANS,'YES',SCR)
	    IF ((ans(1:1).EQ.'n').OR.(ans(1:1).EQ.'N')) THEN
	      CALL CACCEPT(file,image_file(1:L),'Enter new file name:')
	      image_file=file
	      GO TO 1
	    ELSE
	      OPEN (UNIT=13,FILE=image_file,STATUS='OLD')
	      CLOSE (UNIT=13,STATUS='DELETE')
	    END IF
	  END IF

c	DO j=1,nDrawRay
c	  WRITE(4,102)((DrawEltVec(i,j)),i=1,nDrawElt(j))
c	  WRITE(4,103)((DrawRayVec(1,i,j)),i=1,nDrawElt(j))
c	  WRITE(4,103)((DrawRayVec(2,i,j)),i=1,nDrawElt(j))
c	END DO

	k=1
	R1(k)=nDrawRay
	DO j=1,nDrawRay
	  k=k+1
	  R1(k)=nDrawElt(j)
	END DO

	DO j=1,nDrawRay
	  DO i=1,nDrawElt(j)
	    k=k+1
	    R1(k)=DrawEltVec(i,j)
	  END DO
	END DO

	DO j=1,nDrawRay
	  DO i=1,nDrawElt(j)
	    k=k+1
	    R1(k)=DrawRayVec(1,i,j)
	    k=k+1
	    R1(k)=DrawRayVec(2,i,j)
	  END DO
	END DO

C  Text format file
	IF (iPloTyp.EQ.5) THEN
	  OPEN (UNIT=13,FILE=image_file,STATUS='NEW',FORM='FORMATTED')
	  WRITE(*,*) 'Writing ',CTITLE
	  WRITE(*,*) 'Text File=',image_file
	  DO i=1,k
	    WRITE(13,104)(R1(i))
	  END DO
	  CLOSE(13)

C  Binary format file
	ELSE IF (iPloTyp.EQ.6) THEN
	  OPEN (UNIT=13,FILE=image_file,ACCESS='DIRECT',
     &	        RECL=RecordLength,STATUS='NEW')
	  WRITE(*,*) 'Writing ',CTITLE
	  WRITE(*,*) 'Binary File=',image_file
	  DO i=1,k
	    WRITE(13,REC=i)(R1(i))
	  END DO
	  CLOSE(13)

C  Fits format file
	ELSE IF (iPloTyp.EQ.8) THEN
	  WRITE(*,*)'Sorry, FITS not supported for SPOT or DRAW'

! C  Matlab format file
! 	ELSE IF (iPloTyp.EQ.9) THEN
! 	  WRITE(*,*) 'Writing ',CTITLE
! 	  WRITE(*,*) 'MATLAB File=',image_file
! 	  WRITE(*,108) k,1
!  	  CALL RWRITE_MAT_LOH(R1,k,1,k,1,image_file,imag_file)

C  Unknown format
	ELSE
	  WRITE(*,*)' Invalid plot type specified: iPloTyp=',iPloTyp
	  RETURN
	END IF

C  End
	RETURN
	END SUBROUTINE DRAWOUT

C***********************************************************************

	SUBROUTINE ArchInit(iPloTyp,ifGraphics)

	INTEGER ArchFlag,Sun,Mac,DEC
	PARAMETER(Sun=0,Mac=1,DEC=2,ArchFlag=Sun)

	LOGICAL ifGraphics
	INTEGER iPloTyp

C  Set architecture-dependent variables
	IF (ArchFlag.EQ.0) THEN
	  iPloTyp=2
	  ifGraphics=.TRUE.
	ELSE IF (ArchFlag.EQ.1) THEN
	  iPloTyp=5
	  ifGraphics=.FALSE.
	ELSE
	  iPloTyp=2
	  ifGraphics=.TRUE.
	END IF

C  End
 	RETURN
	END SUBROUTINE ArchInit

C***********************************************************************
C On-screen display

	SUBROUTINE SRFPLT(iPloTyp,X,Y,Z,NA,NC,NDA,NDC,
     &	                  CXLABL,CYLABL,CTITLE,R1,R2,StretchType)
	use param_mod
#ifdef CSMACOS
        use smacosio_mod
#else
        use macosio_mod
	use pgplot_mod
#endif
	CHARACTER*(*) CXLABL,CYLABL,CTITLE
	CHARACTER*1 TAB
	CHARACTER*82 SCR
	INTEGER I,J,NA,NC,NDA,NDC,iPloTyp,nContours,StretchType,NAC
	SREAL X(NDA),Y(NDC),Z(NDA,NDC),R1(mdttl,mdttl),R2(NDA,NDC),S

 102	FORMAT(12X,',',64(E11.4,','))
 103	FORMAT(1X,256(E11.4,A1))

	TAB=','
	nContours=16

	NAC=NDA*NC
	IF (StretchType.EQ.2) THEN
	  DO J=1,NC
	    DO I=1,NA
	      S=Z(I,J)
	      IF (S.GT.1e-24) THEN
	        R2(I,J)=LOG10(S)
	      ELSE
	        R2(I,J)=-24e0
	      END IF
	    END DO
	  END DO
	ELSE IF (StretchType.EQ.3) THEN
	  DO J=1,NC
	    DO I=1,NA
	      S=Z(I,J)
	      IF (S.GE.0e0) THEN
	        R2(I,J)=SQRT(S)
	      ELSE
	        R2(I,J)=-SQRT(-S)
	      END IF
	    END DO
	  END DO
	ELSE
	  DO J=1,NC
	    DO I=1,NA
	      R2(I,J)=Z(I,J)
	    END DO
	  END DO
	END IF

	IF (iPloTyp.EQ.1) THEN
          !print*,'**B4 calling SLICE in SRFPLT in sunsub.F'
	  CALL SLICE(X,Y,R2,NA,NC,NDA,NDC,
     &	             R1(1,1),R1(1,2),CXLABL,CYLABL,CTITLE)

	ELSE IF (iPloTyp.EQ.2) THEN
	  CALL GRAY(X,Y,R2,NA,NC,NDA,NDC,CXLABL,CYLABL,CTITLE)

	ELSE IF (iPloTyp.EQ.3) THEN
	  CALL CONTOUR(X,Y,R2,NA,NC,NDA,NDC,nContours,
     &	  	       CXLABL,CYLABL,CTITLE)

	ELSE IF (iPloTyp.EQ.4) THEN
	  CALL WIRE(X,Y,R2,R1,NA,NC,NDA,NDC,md2,CXLABL,CYLABL,CTITLE)

	ELSE IF (iPloTyp.EQ.5) THEN
	  WRITE(4,*)CTITLE
	  DO 1 J=NC,1,-1
	    WRITE(4,103)(R2(I,J),TAB,I=1,NA)
   1	  CONTINUE
	  SCR = 'Writing '//CTITLE
	  WRITE(*,*) SCR

	ELSE IF (iPloTyp.EQ.7) THEN
	  I=NC/2+1
 	  CALL IACCEPT_S(J,I,1,
     &	  'Enter number of column to be plotted:')
	  IF ((J.LE.0).OR.(J.GT.NC)) RETURN
	  DO 2 I=1,NA
	    R1(I,1)=R2(J,I)
   2	  CONTINUE
	  CALL PLOTCOL(X,Y,R1,NA,NC,NDA,NDC,CXLABL,CYLABL,CTITLE)

	ELSE IF (iPloTyp.EQ.10) THEN
	  I=NC/2+1
 	  CALL IACCEPT_S(J,I,1,
     &	  'Enter number of row to be plotted:')
	  IF ((J.LE.0).OR.(J.GT.NC)) RETURN
	  CALL PLOTCOL(X,Y,R2(1,J),NA,NC,NDA,NDC,CXLABL,CYLABL,CTITLE)

	END IF
	RETURN
	END SUBROUTINE SRFPLT

C***********************************************************************

C Render in pgplot window
	SUBROUTINE DSRFPLT(iPloTyp,X,Y,Z,NA,NC,NDA,NDC,
     &		     CXLABL,CYLABL,CTITLE,R1,R2,StretchType)
	use param_mod
#ifdef CSMACOS
        use smacosio_mod
#else
        use macosio_mod
	use pgplot_mod
#endif

	!CHARACTER*(*) CXLABL,CYLABL,CTITLE
	CHARACTER(len=*) :: CXLABL,CYLABL,CTITLE
	CHARACTER*1 TAB
	CHARACTER*82 SCR
	INTEGER I,J,NA,NC,NDA,NDC,iPloTyp,nContours,StretchType,NAC
	SREAL X(NDA),Y(NDC),R1(mdttl,mdttl),R2(NDA,NDC),S
	REAL*8 Z(NDA,NDC)

 102	FORMAT(12X,',',64(E11.4,','))
 103	FORMAT(1X,256(E11.4,A1))

	TAB=','
	nContours=16

	NAC=NDA*NC
	IF (StretchType.EQ.2) THEN
	  DO J=1,NC
	    DO I=1,NA
	      S=Z(I,J)
	      IF (S.GT.1e-24) THEN
	        R2(I,J)=LOG10(S)
	      ELSE
	        R2(I,J)=-24e0
	      END IF
	    END DO
	  END DO
	ELSE IF (StretchType.EQ.3) THEN
	  DO J=1,NC
	    DO I=1,NA
	      S=Z(I,J)
	      IF (S.GE.0e0) THEN
	        R2(I,J)=SQRT(S)
	      ELSE
	        R2(I,J)=-SQRT(-S)
	      END IF
	    END DO
	  END DO
	ELSE
	  DO J=1,NC
	    DO I=1,NA
	      R2(I,J)=Z(I,J)
	    END DO
	  END DO
	END IF

	IF (iPloTyp.EQ.1) THEN
	  CALL SLICE(X,Y,R2,NA,NC,NDA,NDC,
     &	    R1(1,1),R1(1,2),CXLABL,CYLABL,CTITLE)

	ELSE IF (iPloTyp.EQ.2) THEN
	  ! here it is used when 'OPD' command is run
	  !print*,'*** b4 call GRAY: NA,NC,NDA,NDC =',NA,NC,NDA,NDC
 	  CALL GRAY(X,Y,R2,NA,NC,NDA,NDC,
     &	    CXLABL,CYLABL,CTITLE)

	ELSE IF (iPloTyp.EQ.3) THEN
	  CALL CONTOUR(X,Y,R2,NA,NC,NDA,NDC,nContours,
     &	    CXLABL,CYLABL,CTITLE)

	ELSE IF (iPloTyp.EQ.4) THEN
	  CALL WIRE(X,Y,R2,R1,NA,NC,NDA,NDC,md2,
     &	    CXLABL,CYLABL,CTITLE)

	ELSE IF (iPloTyp.EQ.5) THEN
	  WRITE(4,*)CTITLE
	  DO 1 J=NC,1,-1
	    WRITE(4,103)(R2(I,J),TAB,I=1,NA)
   1	  CONTINUE
	  SCR = 'Writing '//CTITLE
	  WRITE(*,*) SCR

	ELSE IF (iPloTyp.EQ.7) THEN
	  I=NC/2+1
 	  CALL IACCEPT_S(J,I,1,
     &	    'Enter number of column to be plotted:')
	  IF ((J.LE.0).OR.(J.GT.NC)) RETURN
	  DO 2 I=1,NA
	    R1(I,1)=R2(J,I)
   2	  CONTINUE
	  CALL PLOTCOL(X,Y,R1,NA,NC,NDA,NDC,CXLABL,CYLABL,CTITLE)

	ELSE IF (iPloTyp.EQ.10) THEN
	  I=NC/2+1
 	  CALL IACCEPT_S(J,I,1,
     &	    'Enter number of row to be plotted:')
	  IF ((J.LE.0).OR.(J.GT.NC)) RETURN
	  CALL PLOTCOL(X,Y,R2(1,J),NA,NC,NDA,NDC,CXLABL,CYLABL,CTITLE)

	END IF
	RETURN
	END SUBROUTINE DSRFPLT

C***********************************************************************

	SUBROUTINE BINOUT(iPloTyp,X,Y,Z,NA,NC,NDA,NDC,
     &	  CXLABL,CYLABL,CTITLE,image_file,dx_elt,dy_elt,iElt)
	use param_mod
	use cfiles_mod
#ifdef CSMACOS
        use smacosio_mod
#else
        use macosio_mod
#endif
	LOGICAL exist
	CHARACTER*1 TAB
	CHARACTER(len=MacosCharLen) :: file,ANS
	CHARACTER*80 imag_file
	CHARACTER*82 SCR
	CHARACTER*(*) CXLABL,CYLABL,CTITLE,image_file
	INTEGER NA,NC,NDA,NDC,iPloTyp,L,ICLEN
	INTEGER*4  iElt
	SREAL     X(NDA),Y(NDC),Z(NDA,NDC)
	REAL*8     dx_elt, dy_elt


C***********************************************************************
C    CONSTANT DEFINITIONS
C***********************************************************************

	INTEGER*4  BYTES_PER_PIXEL, BITS_PER_BYTE, BITS_PER_PIXEL
	PARAMETER(BYTES_PER_PIXEL = 4)
	PARAMETER(BITS_PER_BYTE   = 8)
	PARAMETER(BITS_PER_PIXEL  = BITS_PER_BYTE * BYTES_PER_PIXEL)

 102	FORMAT(' Text array is of dimension',i6,' by',i6,
     &	'. Magnitude of central pixel is ',D17.10)
 103	FORMAT(' Binary array is of dimension',i6,' by',i6,
     &	'. Magnitude of central pixel is ',D17.10)
 104	FORMAT(1P,512(1x,e13.6))
 105	FORMAT(1P,I4,A1,I4/E11.4,A1,E11.4)
 106	FORMAT(1P,512(E11.4,A1))

	L=ICLEN(image_file)

	IF (iPloTyp.EQ.5) THEN
	  file=image_file(1:L)//'.txt'
	  image_file=file
	ELSE IF (iPloTyp.EQ.8) THEN
	  file=image_file(1:L)//'.fit'
	  image_file=file
	ELSE IF (iPloTyp.EQ.9) THEN
	  file=image_file(1:L)//'.mat'
	  image_file=file

C Alternative code...
C   ask for file name
C	  file=image_file(1:L)//'.mat'
C	  image_file=file
C	  CALL CACCEPT(file,image_file,'Enter new file name:')
C	  imag_file=file
C
C  check to make sure it ends in .mat, if not, add it
C
C	  L=ICLEN(imag_file)
C          IF (imag_file((L-3):L) .NE. '.mat') THEN
C	    imag_file((L+1):(L+4))='.mat'
C	  END IF

	END IF

 1	CONTINUE
	  INQUIRE (FILE=image_file,EXIST=exist)
	  IF (exist) THEN
	    L=ICLEN(image_file)
C alt	    IF (iPloTyp.EQ.9) THEN
C alt	      WRITE(*,*) 'Variables being added to .mat file'
C alt	    ELSE
	      SCR='Image file '//image_file(1:L)//' exists.  Overwrite?'
	      CALL CACCEPT(ANS,'YES',SCR)
	      IF ((ans(1:1).EQ.'n').OR.(ans(1:1).EQ.'N')) THEN
	        CALL CACCEPT(file,image_file,'Enter new file name:')
	        image_file=file
	        GO TO 1
	      ELSE
	        OPEN (UNIT=13,FILE=image_file,STATUS='OLD')
	        CLOSE (UNIT=13,STATUS='DELETE')
	      END IF
C alt	    END IF
	  END IF

	IF (iPloTyp.EQ.5) THEN
	  OPEN (UNIT=13,FILE=image_file,STATUS='NEW',FORM='FORMATTED')
	  SCR = 'Writing '//CTITLE
	  WRITE(*,*) SCR
	  SCR = 'FORMATTED File='//image_file
	  WRITE(*,*) SCR
	  WRITE(*,102)NA,NC,Z(NA/2,NC/2)
	  DO 2 j=1,NC
	    WRITE(13,104)(Z(i,j),i=1,NA)
 2	  CONTINUE
	  CLOSE(13)

C "16" really means 'ASCII Special' SpyGlass format in macsub
	ELSE IF (iPloTyp.EQ.16) THEN
	  OPEN (13,FILE=imgfil,STATUS='NEW',FORM='FORMATTED')
	  WRITE(*,*)'Writing ',CTITLE
	  WRITE(*,*)'File=',imgfil
	  WRITE(*,103)NA,NC,Z(NA/2,NC/2)
	  TAB=','
	  ZMIN=1E32
	  ZMAX=-1E32
	  DO J=1,NC
	    DO I=1,NA
	     IF (ZMAX.LT.Z(I,J)) ZMAX=Z(I,J)
	     IF (ZMIN.GT.Z(I,J)) ZMIN=Z(I,J)
	    END DO
	  END DO
	  WRITE(13,105)NA,TAB,NC,ZMAX,TAB,ZMIN
	  WRITE(13,106)(X(I),TAB,I=1,NA)
	  WRITE(13,106)(Y(I),TAB,I=NC,1,-1)
	  DO J=NC,1,-1
	    WRITE(13,106)(Z(I,J),TAB,I=1,NA)
	  END DO
	  CLOSE(13)

	ELSE IF (iPloTyp.EQ.6) THEN
C -- use this for sun:
	  OPEN (UNIT=13,FILE=image_file,ACCESS='DIRECT',RECL=(4*NC),
     &	  STATUS='NEW')
C -- use this for alpha or mac:
C	  OPEN (UNIT=13,FILE=image_file,ACCESS='DIRECT',RECL=(NC),
C     &	  STATUS='NEW')
	  SCR = 'Writing '//CTITLE
	  WRITE(*,*) SCR
	  SCR = 'DIRECT ACCESS File='//image_file
	  WRITE(*,*) SCR
	  WRITE(*,103)NA,NC,Z(NA/2,NC/2)
	  DO j=1,NC
	    WRITE(13,REC=j)(Z(i,j),i=1,NA)
	  END DO
	  CLOSE(13)

	ELSE IF (iPloTyp.EQ.8) THEN
	  SCR = 'Writing '//CTITLE
	  WRITE(*,*) SCR
	  SCR = 'FORMATTED File='//image_file
	  WRITE(*,*) SCR
	  WRITE(*,102)NA,NC,Z(NA/2,NC/2)
	  CALL WRITE_FITS_FILE(Z, NDA, NDC, NA, NC, BITS_PER_PIXEL,
     &	                  filnam, image_file, dx_elt,dy_elt, iElt)

! 	ELSE IF (iPloTyp.EQ.9) THEN
! c	  OPEN (UNIT=13,FILE=image_file,STATUS='NEW',FORM='FORMATTED')
! 	  SCR = 'Writing '//CTITLE
! 	  WRITE(*,*) SCR
! 	  SCR = 'FORMATTED File='//image_file
! 	  WRITE(*,*) SCR
! 	  WRITE(*,102)NA,NC,Z(NA/2,NC/2)

!  	  CALL RWRITE_MAT_LOH(Z,NA,NC,NDA,NDC,image_file,imag_file)
	ELSE
	  WRITE(*, *) ' Unknown plot type iPloTyp = ', iPloTyp
	END IF
 	RETURN
	END SUBROUTINE BINOUT

C***********************************************************************

	SUBROUTINE DBINOUT(iPloTyp,X,Y,DZ,NA,NC,NDA,NDC,
     &	   CXLABL,CYLABL,CTITLE,image_file,dx_elt,dy_elt,iElt,Z)
	use param_mod
	use cfiles_mod
#ifdef CSMACOS
        use smacosio_mod
#else
        use macosio_mod
#endif
	LOGICAL exist
	CHARACTER*1 TAB
	CHARACTER(len=MacosCharLen) :: file,ANS
	CHARACTER*80 imag_file
	CHARACTER*82 SCR
	CHARACTER*(*) CXLABL,CYLABL,CTITLE,image_file
	INTEGER NA,NC,NDA,NDC,iPloTyp,L,ICLEN
	INTEGER*4  iElt
	SREAL     X(NDA),Y(NDC),Z(NDA,NDC)
	REAL*8    dx_elt, dy_elt, DZ(NDA,NDC)
C
C***********************************************************************
C    CONSTANT DEFINITIONS
C***********************************************************************
C
	INTEGER*4 BYTES_PER_PIXEL,BITS_PER_BYTE,BITS_PER_PIXEL
	PARAMETER(BYTES_PER_PIXEL =4)
	PARAMETER(BITS_PER_BYTE   =8)
	PARAMETER(BITS_PER_PIXEL
     &            =BITS_PER_BYTE*BYTES_PER_PIXEL)

 102	FORMAT(' Text array is of dimension',i6,' by',i6,
     &	  '. Magnitude of central pixel is ',D17.10)
 103	FORMAT(' Binary array is of dimension',i6,' by',i6,
     &	  '. Magnitude of central pixel is ',D17.10)
 104	FORMAT(1P,512(1x,e13.6))
 105	FORMAT(1P,I4,A1,I4/E11.4,A1,E11.4)
 106	FORMAT(1P,512(E11.4,A1))
	print*,'** In DBINOUT: NA,NC =',NA,NC

	L=ICLEN(image_file)

	DO j=1,NDC
	  DO i=1,NDA
	    Z(i,j)=DZ(i,j)
	  END DO
	END DO

	IF (iPloTyp.EQ.5) THEN
	  file=image_file(1:L)//'.txt'
	  image_file=file
	ELSE IF (iPloTyp.EQ.8) THEN
	  file=image_file(1:L)//'.fit'
	  image_file=file
	ELSE IF (iPloTyp.EQ.9) THEN
	  file=image_file(1:L)//'.mat'
	  image_file=file

C Alternative code...
C   ask for file name
C	  file=image_file(1:L)//'.mat'
C	  image_file=file
C	  CALL CACCEPT(file,image_file,'Enter new file name:')
C	  imag_file=file
C
C  check to make sure it ends in .mat, if not, add it
C
C	  L=ICLEN(imag_file)
C          IF (imag_file((L-3):L) .NE. '.mat') THEN
C	    imag_file((L+1):(L+4))='.mat'
C	  END IF
	END IF

 1	CONTINUE
	INQUIRE (FILE=image_file,EXIST=exist)
	IF (exist) THEN
	  L=ICLEN(image_file)
C alt	  IF (iPloTyp.EQ.9) THEN
C alt	      WRITE(*,*) 'Variables being added to .mat file'
C alt	  ELSE
C alt	  END IF
	  SCR='Image file '//image_file(1:L)//' exists.  Overwrite?'
	  CALL CACCEPT(ANS,'YES',SCR)
	  IF ((ans(1:1).EQ.'n').OR.(ans(1:1).EQ.'N')) THEN
	    CALL CACCEPT(file,image_file,'Enter new file name:')
	    image_file=file
	    GO TO 1
	  ELSE
	    OPEN (UNIT=13,FILE=image_file,STATUS='OLD')
	    CLOSE (UNIT=13,STATUS='DELETE')
	  END IF
	END IF

	IF (iPloTyp.EQ.5) THEN
	  OPEN (UNIT=13,FILE=image_file,STATUS='NEW',FORM='FORMATTED')
	  SCR = 'Writing '//CTITLE
	  WRITE(*,*) SCR
	  SCR = 'FORMATTED File='//image_file
	  WRITE(*,*) SCR
	  WRITE(*,102)NA,NC,Z(NA/2,NC/2)
	  DO j=1,NC
	    WRITE(13,104)(Z(i,j),i=1,NA)
  	  END DO
	  CLOSE(13)

C "16" really means 'ASCII Special' SpyGlass format in macsub
	ELSE IF (iPloTyp.EQ.16) THEN
	  OPEN (13,FILE=imgfil,STATUS='NEW',FORM='FORMATTED')
	  WRITE(*,*)'Writing ',CTITLE
	  WRITE(*,*)'File=',imgfil
	  WRITE(*,103)NA,NC,Z(NA/2,NC/2)
	  TAB=','
	  ZMIN=1E32
	  ZMAX=-1E32
	  DO J=1,NC
	    DO I=1,NA
	     IF (ZMAX.LT.Z(I,J)) ZMAX=Z(I,J)
	     IF (ZMIN.GT.Z(I,J)) ZMIN=Z(I,J)
	    END DO
	  END DO
	  WRITE(13,105)NA,TAB,NC,ZMAX,TAB,ZMIN
	  WRITE(13,106)(X(I),TAB,I=1,NA)
	  WRITE(13,106)(Y(I),TAB,I=NC,1,-1)
	  DO J=NC,1,-1
	    WRITE(13,106)(Z(I,J),TAB,I=1,NA)
	  END DO
	  CLOSE(13)

	ELSE IF (iPloTyp.EQ.6) THEN
	  OPEN (UNIT=13,FILE=image_file,ACCESS='DIRECT',
     &          RECL=(realRecLen*NA),STATUS='NEW')
C -- use this for alpha or mac:
C	  OPEN (UNIT=13,FILE=image_file,ACCESS='DIRECT',RECL=(NC),
C     &	        STATUS='NEW')
	  SCR = 'Writing '//CTITLE
	  WRITE(*,*) SCR
	  SCR = 'DIRECT ACCESS File='//image_file
	  WRITE(*,*) SCR
	  WRITE(*,103) NA,NC,Z(NA/2,NC/2)
	  DO j=1,NC
	    WRITE(13,REC=j)(Z(i,j),i=1,NA)
	  END DO
	  CLOSE(13)

	ELSE IF (iPloTyp.EQ.8) THEN
	  SCR = 'Writing '//CTITLE
	  WRITE(*,*) SCR
	  SCR = 'FORMATTED File='//image_file
	  WRITE(*,*) SCR
	  WRITE(*,102)NA,NC,Z(NA/2,NC/2)
	  CALL WRITE_FITS_FILE(Z, NDA, NDC, NA, NC, BITS_PER_PIXEL,
     &	                  filnam, image_file, dx_elt,dy_elt, iElt)

! 	ELSE IF (iPloTyp.EQ.9) THEN
! c	  OPEN (UNIT=13,FILE=image_file,STATUS='NEW',FORM='FORMATTED')
! 	  SCR = 'Writing '//CTITLE
! 	  WRITE(*,*) SCR
! 	  SCR = 'FORMATTED File='//image_file
! 	  WRITE(*,*) SCR
! 	  WRITE(*,102)NA,NC,Z(NA/2,NC/2)

!  	  CALL RWRITE_MAT_LOH(Z,NA,NC,NDA,NDC,image_file,imag_file)
	ELSE
	  WRITE(*, *) ' Unknown plot type iPloTyp = ', iPloTyp
	END IF
 	RETURN
	END SUBROUTINE DBINOUT


C***********************************************************************

	SUBROUTINE WRITE_FITS_FILE(image, array_row_size,
     &	  array_col_size, image_row_size, image_col_size,
     &	  bits_per_pixel, root_file_name, output_file_name,
     &	  dx_elt, dy_elt, i_elt)

	Implicit       None

	INTEGER*4      image_row_size, image_col_size
	INTEGER*4      bits_per_pixel
	INTEGER*4      array_row_size, array_col_size
	INTEGER        i, num_elm
	INTEGER        lout, ier, root_file_name_len
	INTEGER        output_file_name_len
	INTEGER        bit_representation
	INTEGER        i_elt
	SREAL          image(array_row_size, array_col_size)
	SREAL          max_image_value, min_image_value
	REAL*8         dx_elt, dy_elt
	SREAL          scale(3), zero(3), error(3)
	INTEGER*4      strLen
	CHARACTER(len=MacosCharLen) :: root_file_name,
     &                                 output_file_name
	CHARACTER*40   new_file_name, actual_file_name
	CHARACTER*24   date
	CHARACTER*16   int_str
	LOGICAL        dont_check_for_rounding_errors

	dont_check_for_rounding_errors = .FALSE.
	output_file_name_len           = strLen(output_file_name)
	root_file_name_len             = strLen(root_file_name)
c
c     see comments of fits_scalc SUBROUTINE (included in this file)
c     for explanation of 'bit_representation'
c
	IF (bits_per_pixel .EQ. 32) THEN
	   bit_representation = 3
	Elseif (bits_per_pixel .EQ. 16) THEN
	   bit_representation = 2
	Elseif (bits_per_pixel .EQ. 8) THEN
	   bit_representation = 1
	ELSE
	   WRITE(*,*) 'Cannot WRITE FITS file for given image.'
	   WRITE(*,*) 'No bit representation exists for a'
	   WRITE(*,*)  bits_per_pixel, ' bits per pixel image.'
	END IF

	num_elm = array_row_size * array_col_size

        new_file_name = output_file_name(1:output_file_name_len)

	IF (root_file_name_len .GT. 18) THEN
	   WRITE(*,*) 'File name ' //
     &	             root_file_name(1:root_file_name_len) //
     &	            ' too long for FITS FILENAME data field.' //
     &	            ' Enter a shorter one' //
     &	            '(< 19 characters).'
	   GO TO 1
	   END IF
	CALL fits_init(new_file_name, actual_file_name, lout, ier)
c sws debugging
c	print *, 'fits_init called'
c	print *,'file= [',new_file_name,'], result =[',
c	1    actual_file_name,']'
	CALL fdate(date)
c
c     fits_scalc, in addition to computing the scaling and zero values,
c     it also computes the min and max values of image array
c
	CALL fits_scalc(image, num_elm, dont_check_for_rounding_errors,
     &	min_image_value,max_image_value, scale, zero, error,
     &	image_row_size, image_col_size,array_col_size)

	CALL  fits_wlog('SIMPLE', .TRUE., 'FITS standard ', ier)
	CALL  fits_wint('BITPIX', bits_per_pixel,
     &	              'Number of bits per pixel', ier)
	CALL  fits_wint('NAXIS',  2, 'Number of axes', ier)
	CALL  fits_wint('NAXIS1', image_row_size,
     &	              'Pixel length of axis 1', ier)
	CALL  fits_wint('NAXIS2', image_col_size,
     &	              'Pixel length of axis 2', ier)
	CALL  fits_wlog('EXTEND', .TRUE.,
     &	              'There maybe standard extensions',
     &	               ier)
	CALL  fits_wreal('BSCALE', scale(bit_representation),
     &	               'REAL = array_value*bscale + bzero',
     &	               ier)
	CALL  fits_wreal('BZERO', zero(bit_representation),
     &	               ' ',
     &	               ier)
	CALL  fits_wstr('ORIGIN', 'COMP',
     &	              'Name of PROGRAM which generated image', ier)
c
c     date(5:11) represents the month and date portion of date string,
c     and date(21:24) represents the year
c
	CALL  fits_wstr('FITSDATE', date(5:11) // date(21:24),
     &	              'Date FITS file was created',
     &	               ier)
	CALL  fits_wstr('FILENAME',
     &	                root_file_name(1:root_file_name_len),
     &	               'Name of original file name',
     &	               ier)
	CALL  fits_wreal('CRVAL1', 0d0,
     &	               'Center detector position along axis 1',
     &	                ier)
	CALL  fits_wreal('CRVAL2', 0d0,
     &	               'Center detector position along axis 2',
     &	                ier)
	CALL  fits_wreal('CRPIX1', image_row_size / 2d0,
     &	               'Center reference pixel along axis 1',
     &	                ier)
	CALL  fits_wreal('CRPIX2', image_col_size / 2d0,
     &	               'Center reference pixel along axis 2',
     &	                ier)
	CALL  fits_wdble('CDELT1', dx_elt, 'dxElt(iElt)',ier)
	CALL  fits_wdble('CDELT2', dy_elt, 'dyElt(iElt)',ier)
	CALL  fits_wreal('DATAMIN', min_image_value,
     &	               'Minimum data value in image',
     &	                ier)
	CALL  fits_wreal('DATAMAX', max_image_value,
     &	               'Maximum data value in image',
     &	                ier)
	CALL  fits_wstr('CTYPE1', 'COLUMN','Coord. name of axis 1',ier)
	CALL  fits_wstr('CTYPE2', 'ROW', 'Coord. name of axis 2',ier)
	CALL  IntToChar(int_str, i_elt, i)
	CALL  fits_wstr('HISTORY', 'ELT' // int_str(1:i), ' ', ier)
	CALL  fits_wstr('FILTNAM1', 'WAVELEN=XX', 'First filter name',
     &	               ier)
	CALL  fits_wstr('IMAGETYP', 'SIMULATED',
     &	              'Image type characteristics',
     &	               ier)
c
c     creates "END" field
c
	CALL  fits_wend(ier)
c grw -- 1/18/94
	CALL  fits_wrayf(image, num_elm, bits_per_pixel,
     &	scale(bit_representation),
     &	zero(bit_representation), ier,
     &	image_row_size, image_col_size,array_col_size)
	CALL  fits_close(ier)
 1      RETURN
	END SUBROUTINE WRITE_FITS_FILE
c
c
	FUNCTION StrLen(stringName)
	CHARACTER   stringName*(*)
	INTEGER*4   index
	INTEGER*4   StrLen

	index=1
 10     IF (stringName(index:index).NE.' ') THEN
          index=index+1
	  GO TO 10
	END IF
	StrLen=index-1
	END
C
C----------------------------------------------------------------------
	SUBROUTINE FITS_WLOG (NAME,VALUE,COMMENT,STATUS)
C-----------------------------------------------------------------------
C     Writes a LOGICAL keyword into a FITS header.  Assumes
C     FITS_INIT has been called to initialise the header routines.
C
C     Parameters -  (">" input, "<" output)
C
C     (>) NAME      (CHARACTER) The name of the keyword.
C     (>) VALUE     (LOGICAL) The value of the PARAMETER.
C     (>) COMMENT   (CHARACTER) Any comment to be associated with the
C                   keyword.  Should be no longer than 46 characters.
C     (>) STATUS    (INTEGER) RETURN status code. 0 => OK, non-zero
C                   values will be error codes caused by an
C                   I/O error.  These can be decoded by FITS_ERROR.
C
C     Note: The FITS standard says keyword names and values must be
C     output in upper case.  NAME can be passed in lower case, but
C     it will be output in upper.  COMMENT will be output as passed.
C
C     COMMON variables used -
C
C     (<) BUFF      (CHARACTER) Buffer for one header line
C
C     Functions / subroutines used -
C
C     FITS_BUFFO     (FITS_ package) Outputs BUFF to main I/O buffer.
C     FITS_BUFP      ( "     "    ) Copies NAME and COMMENTS into BUFF.
C
C                                                KS / CIT  9th OCT 1983
C-----------------------------------------------------------------------
C
	LOGICAL VALUE
	INTEGER STATUS
	CHARACTER*(*) NAME,COMMENT
C
	INCLUDE 'fits1.inc'
C
	CALL FITS_BUFP(NAME,COMMENT)
	IF (VALUE) THEN
	   BUFF(30:30)='T'
	ELSE
	   BUFF(30:30)='F'
	END IF
	CALL FITS_BUFFO(STATUS)
C
	END
C
C----------------------------------------------------------------------
	SUBROUTINE FITS_WINT (NAME,VALUE,COMMENT,STATUS)
C-----------------------------------------------------------------------
C     Writes an INTEGER keyword into a FITS header.  Assumes
C     FITS_INIT has been called to initialise the header routines.
C
C     Parameters -  (">" input, "<" output)
C
C     (>) NAME      (CHARACTER) The name of the keyword.
C     (>) VALUE     (INTEGER) The value of the PARAMETER.
C     (>) COMMENT   (CHARACTER) Any comment to be associated with the
C                   keyword.  Should be no longer than 46 characters.
C     (>) STATUS    (INTEGER) RETURN status code. 0 => OK, non-zero
C                   values will be error codes caused by an
C                   I/O error.  These can be decoded by FITS_ERROR.
C
C     Note: The FITS standard says keyword names and values must be
C     output in upper case.  NAME can be passed in lower case, but
C     it will be output in upper.  COMMENT will be output as passed.
C
C     COMMON variables used -
C
C     (<) BUFF      (CHARACTER) Buffer for one header line
C
C     Functions / subroutines used -
C
C     FITS_BUFFO     (FITS_ package) Outputs BUFF to main I/O buffer.
C     FITS_BUFP      ( "     "    ) Copies NAME and COMMENTS into BUFF.
C
C                                                KS / CIT  9th OCT 1983
C-----------------------------------------------------------------------
C
	INTEGER STATUS,VALUE
	CHARACTER*(*) NAME,COMMENT
C
	INCLUDE 'fits1.inc'
C
	CALL FITS_BUFP(NAME,COMMENT)
	WRITE (BUFF(11:30),'(I20)',IOSTAT=STATUS) VALUE
	CALL FITS_BUFFO(STATUS)
C
	END
C
C----------------------------------------------------------------------
	SUBROUTINE FITS_WREAL (NAME,VALUE,COMMENT,STATUS)
C-----------------------------------------------------------------------
C     Writes a REAL keyword into a FITS header.  Assumes
C     FITS_INIT has been called to initialise the header routines.
C
C     Parameters -  (">" input, "<" output)
C
C     (>) NAME      (CHARACTER) The name of the keyword.
C     (>) VALUE     (REAL) The value of the PARAMETER.
C     (>) COMMENT   (CHARACTER) Any comment to be associated with the
C                   keyword.  Should be no longer than 46 characters.
C     (>) STATUS    (INTEGER) RETURN status code. 0 => OK, non-zero
C                   values will be error codes caused by an
C                   I/O error.  These can be decoded by FITS_ERROR.
C
C     Note: The FITS standard says keyword names and values must be
C     output in upper case.  NAME can be passed in lower case, but
C     it will be output in upper.  COMMENT will be output as passed.
C
C     COMMON variables used -
C
C     (<) BUFF      (CHARACTER) Buffer for one header line
C
C     Functions / subroutines used -
C
C     FITS_BUFFO     (FITS_ package) Outputs BUFF to main I/O buffer.
C     FITS_BUFP      ( "     "    ) Copies NAME and COMMENTS into BUFF.
C
C                                                KS / CIT  9th OCT 1983
C-----------------------------------------------------------------------
C
	INTEGER STATUS
	SREAL VALUE
	CHARACTER*(*) NAME,COMMENT
C
	INCLUDE 'fits1.inc'
C
	CALL FITS_BUFP(NAME,COMMENT)
	WRITE (BUFF(11:30),'(1PE20.6)',IOSTAT=STATUS) VALUE
	CALL FITS_BUFFO(STATUS)
C
	END
C
C----------------------------------------------------------------------
	SUBROUTINE FITS_WSTR (NAME,STRING,COMMENT,STATUS)
C-----------------------------------------------------------------------
C     Writes a string keyword into a FITS header.  Assumes
C     FITS_INIT has been called to initialise the header routines.
C
C     Parameters -  (">" input, "<" output)
C
C     (>) NAME      (CHARACTER) The name of the keyword.
C     (>) STRING    (CHARACTER) The value of the keyword - this routine
C                   appends the quote marks required by the standard, so
C                   these should not be included in STRING.
C     (>) COMMENT   (CHARACTER) Any comment to be associated with the
C                   keyword.  Normally the comment is started at col 34
C                   and so should be limited to 46 characters.  IF the
C                   keyword itself runs over past col 30, the comment
C                   will have to be shorter.
C     (>) STATUS    (INTEGER) RETURN status code. 0 => OK, non-zero
C                   values will be error codes caused by an
C                   I/O error.  These can be decoded by FITS_ERROR.
C
C     Note: The FITS standard says keyword names and values must be
C     output in upper case.  NAME & STRING can be passed in lower case,
C     but will be output in upper.  COMMENT will be output as passed.
C
C     COMMON variables used -
C
C     (<) BUFF      (CHARACTER) Buffer for one header line
C
C     Functions / subroutines used -
C
C     FITS_BUFFO     (FITS_ package) Outputs BUFF to main I/O buffer.
C     UPCASE         Convert string to upper case.
C
C                                                KS / CIT  15th Nov 1983
C-----------------------------------------------------------------------
C
	INTEGER STATUS
	CHARACTER*(*) NAME,COMMENT,STRING
C
	INTEGER IBUF,ISTR,LPTR
	CHARACTER CHR
C
	INCLUDE 'fits1.inc'
C
	BUFF=NAME
	CALL UPCASE(BUFF(1:8))
	BUFF(9:11)='= '''
	ISTR=1
	IBUF=12
	LPTR=12
	DO WHILE ((ISTR.LE.LEN(STRING)).AND.(IBUF.LT.LEN(BUFF)))
	   CHR=STRING(ISTR:ISTR)
	   CALL UPCASE(CHR)
	   BUFF(IBUF:IBUF)=CHR
	   IF (CHR.NE.' ') LPTR=IBUF
	   IBUF=IBUF+1
	   ISTR=ISTR+1
	END DO
	BUFF(MAX(20,LPTR+1):)=''''
	BUFF(MAX(32,LPTR+1):)='/ '//COMMENT
	CALL FITS_BUFFO(STATUS)
C
	END
C
C-----------------------------------------------------------------------
	SUBROUTINE FITS_BUFP (NAME,COMMENT)
C-----------------------------------------------------------------------
C     Writes NAME and COMMENT into the header line buffer BUFF
C     and sets the rest of BUFF blank.
C
C     Parameters -  (">" input, "<" output)
C
C     (>) NAME       (CHARACTER) The keyword name
C     (>) COMMENT    (CHARACTER) The associated comment
C
C     Note: NAME will be written into BUFF converted to upper case, but
C     COMMENT will be copied over as passed.
C
C     COMMON variables used -
C
C     (<) BUFF       (CHARACTER) Buffer for one header line
C
C     Subroutines / functions used -
C
C     UPCASE         Convert string to upper case.
C
C                                    KS / CIT 9th Oct 1983
C-----------------------------------------------------------------------
C
	CHARACTER*(*) NAME,COMMENT
C
	INCLUDE 'fits1.inc'
C
	BUFF=NAME
	CALL UPCASE(BUFF(1:8))
	BUFF(9:9)='='
	BUFF(32:32)='/'
	BUFF(34:)=COMMENT
C
	END
C
C-----------------------------------------------------------------------
	SUBROUTINE FITS_WEND (STATUS)
C-----------------------------------------------------------------------
C     Terminates a FITS header, by writing an END record and filling
C     up the rest of the FITS buffer with blanks.
C
C     Parameters -    (">" input, "!" modified, "<" output)
C
C     (<) STATUS     (INTEGER) RETURN code.  0 => OK.  Non-zero
C                    values will be error codes caused by an
C                    I/O error.  These can be decoded by FITS_ERROR.
C
C     COMMON variables used -
C
C     (<) BUFF       (CHARACTER) Buffer for one header line
C     (>) FPTR       (INTEGER) Pointer to the main FITS buffer FBUFF.
C
C     Subroutines / functions used -
C
C     FITS_BUFFO      (FITS_ package) Outputs BUFF to main I/O buffer
C
C                                      KS / CIT 9th Oct 1983
C-----------------------------------------------------------------------
C
	INTEGER STATUS
C
	INCLUDE 'fits1.inc'
C
C     WRITE the END line into the header
C
	BUFF='END'
	CALL FITS_BUFFO(STATUS)
C
C     Keep writing blank lines into the header until the pointer
C     is reset - this indicates that the buffer has been output.
C
	BUFF=' '
	DO WHILE ((FPTR.GT.1).AND.(STATUS.EQ.0))
	   CALL FITS_BUFFO(STATUS)
	END DO
C
	END
C
C-----------------------------------------------------------------------
	SUBROUTINE FITS_INIT (FILE,FILEXP,LF,STATUS)
C-----------------------------------------------------------------------
C     Initialises the FITS_ routines.  This should be the first
C     FITS_ routine called.
C
C     Parameters -   (">" input, "<" output)
C
C     (>) FILE     (CHARACTER) The name of the disk file to be
C                  created.
C     (<) FILEXP   (CHARACTER) The name of the file that was actually
C                  created (expanded file name)>
C     (<) STATUS   (INTEGER) A returned status code.  0 => OK, non-
C                  zero values are Fortran I/O error codes.
C
C     COMMON variables used -
C
C     (<) MTUNIT   (INTEGER) The Fortran LOGICAL unit used for output.
C     (<) FPTR     (INTEGER) The pointer to the main output buffer.
C
C-----------------------------------------------------------------------
C
	INTEGER STATUS, LF
	CHARACTER*(*) FILE, FILEXP
	INTEGER VLBOPE
C
C     grw 10/31/93
C
	CHARACTER  response
C
	INCLUDE 'fits1.inc'
C
Cvms  INTEGER LIB$GET_LUN
C
C     Assign a LOGICAL unit
C
	MTUNIT = 77
Cvms  STATUS = LIB$GET_LUN(MTUNIT)
Cvms  IF (STATUS.NE.1) RETURN
C
C     OPEN the file (2880-byte records)
C
	STATUS = VLBOPE(MTUNIT, FILE, 'FITS', 'NEW', FILEXP)
C
C     modified by grw 10/31/93
C
	IF (FILEXP .EQ. '''new'' file exists') THEN

	   WRITE(*,*) 'Would you like to overwrite existing file?'
	   READ(*,*) response
	   IF (response .EQ. 'Y' .OR. response .EQ. 'y') THEN

	     STATUS = VLBOPE(MTUNIT, FILE, 'FITS', 'OLD', FILEXP)
	     END IF
	   END IF

	IF (STATUS.NE.0) THEN
	    STATUS = 0
	ELSE
	    STATUS = 1
	END IF
	LF = LEN(FILEXP)
	DO WHILE (FILEXP(LF:LF).EQ.' ')
	    LF = LF-1
	END DO
C
C     Initialise pointers
C
	BLKNO = 1
	FPTR=1
C
	END
C
C-----------------------------------------------------------------------
	SUBROUTINE FITS_BUFFO (STATUS)
C-----------------------------------------------------------------------
C     Writes the FITS header buffer (BUFF, in COMMON) into the
C     main FITS output buffer (FBUFF).  IF necessary, the buffer
C     will be written to the output file.
C
C     PARAMETER -   (">" input,  "!" modified,  "<" output)
C
C     (<) STATUS    (INTEGER) RETURN code.  0 => OK.  Non-zero codes
C                   represent errors.  This can be decoded by FITS_ERROR.
C
C     COMMON variables used -
C
C     (>) BUFF      (CHARACTER) The single line CHARACTER buffer used
C                   for header keyword lines.
C     (<) FBUFF     (Byte array) The main FITS output buffer.
C     (!) FPTR      (INTEGER) The pointer to the next free byte in FBUFF.
C                   Should have been initialised by FITS_INIT.
C
C     Subroutines / functions used -
C
C     FITS_MWRIT      (FITS_ package) Writes FBUFF to output file.
C
C-----------------------------------------------------------------------
C
	INTEGER STATUS
C
	INCLUDE 'fits1.inc'
C
	INTEGER I
C
C     Copy BUFF into FBUFF, THEN WRITE out IF full.  (Note assumption
C     that FBUFF has been filled 80 bytes at a time, so will not fill in
C     the middle of the copy.)
C
	DO I=1,80
	   FBUFF(FPTR)=ICHAR(BUFF(I:I))
	   FPTR=FPTR+1
	END DO
	IF (FPTR.GT.2880) THEN
	   CALL FITS_MWRIT(STATUS)
	   FPTR=1
	ELSE
	   STATUS=0
	END IF
C
	END
C
C-----------------------------------------------------------------------
	SUBROUTINE FITS_CLOSE (STATUS)
C-----------------------------------------------------------------------
C     Should be the last FITS_ routine called during the writing
C     of a single image.  Flushes the main I/O buffer to disk,
C     and closes the output file.
C
C
C     (<) STATUS    (INTEGER) Returned status code. 0 => OK, non-zero
C                   codes are I/O error codes, which may be decoded
C                   by FITS_ERROR.
C
C     Subroutines / functions used -
C
C     FITS_FLUSH    (FITS_ package) Writes FBUFF to output file.
C
C                                       KS/TJP / CIT 11th Oct 1983
C-----------------------------------------------------------------------
C
	INTEGER STATUS
C
	INCLUDE 'fits1.inc'
C
C     Flush IF necessary.
C
	CALL FITS_FLUSH(STATUS)
	IF (STATUS.NE.0)  GO TO 600
C
C     CLOSE file.
C
	CLOSE (UNIT=MTUNIT, IOSTAT=STATUS)
C
  600 CONTINUE
	END
C
C-----------------------------------------------------------------------
	SUBROUTINE FITS_WRAYS (DATA,NELM,STATUS)
C-----------------------------------------------------------------------
C     Writes an INTEGER*2 array to output file as part of a FITS image.
C
C     This routine simply copies the array element by element into a
C     buffer, writing the buffer out whenever it fills.  This means
C     that this routine may be called as many times as necessary - for
C     example, if all the data is not available in memory at one time.
C     This routine should be called for the first time immediately
C     following the FITS_WEND CALL that terminated the FITS header.
C     Once all the data is output, FITS_CLOSE should be called to flush
C     the data buffer.
C
C     Parameters - (">" input, "W" workspace, "!" modified, "<" output)
C
C     (>) DATA     (INTEGER*2 array DATA(NELM)) The array to be output.
C                  DATA may be multi-dimensional, but it is treated
C                  as linear here for generality.
C     (>) NELM     (INTEGER) The number of elements in DATA.
C     (<) STATUS   (INTEGER) Returned status code.  0 => OK, non-zero
C                  values indicate an I/O error and can be
C                  decoded using FITS_ERROR.
C
C     COMMON variables used -
C
C     (W) FBUFF    (Byte array) The main FITS I/O buffer
C     (W) FBUFFS   (INTEGER*2 array) The same buffer as FBUFF.
C     (!) FPTR     (INTEGER) Byte level pointer to FBUFF.
C
C     Subroutines / functions used -
C
C     FITS_MWRIT    (FITS_    "   ) WRITE FBUFF buffer to output
C
C                                       KS / CIT 19th Oct 1983
C-----------------------------------------------------------------------
C
	INTEGER NELM,STATUS
	INTEGER*2 DATA(NELM)
C
	INCLUDE 'fits1.inc'
C
	INTEGER FPTRO,FPTRS,I
C
C     WRITE out as signed 32 bit reals.
C
	STATUS=0
	FPTRS=FPTR/2+1
	FPTRO=FPTRS
	DO I=1,NELM
	   FBUFFS(FPTRS)=DATA(I)
	   FPTRS=FPTRS+1
	   IF (FPTRS.GT.1440) THEN
csun	      CALL FITS_BSWAP(FBUFFS(FPTRO),FPTRS-FPTRO)
	      CALL FITS_MWRIT(STATUS)
	      FPTR=1
	      FPTRS=1
	      FPTRO=1
	      IF (STATUS.NE.0) GO TO 600
	   END IF
	END DO
csun	IF (FPTRS.GT.1) CALL FITS_BSWAP(FBUFFS(FPTRO),FPTRS-FPTRO)
	FPTR=(FPTRS-1)*2+1
C
  600 CONTINUE
	END

C*UPCASE -- convert CHARACTER string to upper case
C+
	SUBROUTINE UPCASE (SRC)
	CHARACTER*(*) SRC
C
C Convert CHARACTER string to upper case; all lower case letters
C (a-z) are converted to upper case equivalents (A-Z); other characters
C are unchanged. This version assumes the ASCII CHARACTER set.
C
C Arguments:
C  SRC    (in/out) : input string to be converted.
C
C History:
C  1987 nov 12 - TJP.
C-----------------------------------------------------------------------
	INTEGER I, K
C
	DO 10 I=1,LEN(SRC)
	    K = ICHAR(SRC(I:I))
	    IF (K.GE.ICHAR('a') .AND. K.LE.ICHAR('z'))
     &	        SRC(I:I) = CHAR(K-32)
   10 CONTINUE
C
	END


C*VLBOPE -- OPEN a VLB data file [Sun-UNIX]
C+
	INTEGER FUNCTION VLBOPE (UNIT, NAME, TYPE, STATUS, RESULT)
	INTEGER UNIT
	CHARACTER*(*) NAME
	CHARACTER*(*) TYPE
	CHARACTER*(*) STATUS
	CHARACTER*(*) RESULT
C
C Input parameters:
C   UNIT : Fortran unit number to be used
C   NAME : file name
C   TYPE : file type, one of the strings
C      'MERGE'  merge FORMAT
C      'FITS'   disk FITS FORMAT (image or uvfits)
C      'MODEL'  model file (text)
C      'TEXT'   arbitrary text file
C   STATUS :
C      'OLD' to OPEN a file for reading (file must already exist)
C      'NEW' to OPEN a new file for writing (file must not already
C            exist, unless the operating system can create a new
C            version)
C
C Output parameters:
C   VLBOPE : receives 1 IF the file was opened successfully,
C      0 IF the OPEN failed
C   RESULT : receives the full file name, IF the file
C      was opened successfully, or an error message (text) IF the
C      OPEN failed; the declared length of this variable should
C      be big enough to accommodate the result; 255 bytes should be
C      enough.
C
C Discussion: this encapsulates all system-dependent OPEN parameters
C like READONLY, CARRIAGECONTROL in a single SUBROUTINE.  Moving to
C a different operating system should require only this routine to be
C changed, assuming the operating system can cope with Fortran formatted
C and unformatted READ and WRITE.
C-----------------------------------------------------------------------
	INTEGER ITYPE, ISTAT, IER, LEN1, BLKNO, ICLEN
	CHARACTER*11 FMT,CTMP
	CHARACTER*4  CC
	CHARACTER*32 CNAME
	CHARACTER*80 STRING
	COMMON /FITS2/ BLKNO
C
C Verify the TYPE argument.
C
	IF (TYPE.EQ.'MERGE') THEN
	    ITYPE = 1
	    FMT = 'UNFORMATTED'
	    CC = 'NONE'
	ELSE IF (TYPE.EQ.'FITS') THEN
	    ITYPE = 2
	    FMT = 'UNFORMATTED'
	    CC = 'NONE'
	ELSE IF (TYPE.EQ.'MODEL' .OR. TYPE.EQ.'TEXT') THEN
	    ITYPE = 3
	    FMT = 'FORMATTED'
	    CC = 'LIST'
	ELSE
	    CTMP=TYPE
	    CALL ERROR('VLBOPE: invalid argument TYPE='//CTMP)
	END IF
C
C Verify the STATUS argument.
C
	IF (STATUS.EQ.'OLD') THEN
	    ISTAT = 1
	ELSE IF (STATUS.EQ.'NEW') THEN
	    ISTAT = 2
	ELSE
	    CTMP=STATUS
	    CALL ERROR('VLBOPE: invalid argument STATUS='//CTMP)
	END IF
C
C Attempt to OPEN the file.
C  Note -- RECL=720 for mac, 2880 for sun
	IF (ISTAT.EQ.1) THEN
C         -- 'OLD' file
	    IF (ITYPE.EQ.2) THEN
	        OPEN (UNIT=UNIT, FILE=NAME, STATUS='OLD',
C    1              RECORDTYPE='FIXED',
     &	             RECL=2880,
csun     2              RECL=720,
     &	             ACCESS='DIRECT',
C    4              READONLY,
     &	             FORM='UNFORMATTED',
     &	             IOSTAT=IER)
	        BLKNO = 0
	    ELSE
	        OPEN (UNIT=UNIT, FILE=NAME, STATUS='OLD',
C    1              READONLY,
     &	             ACCESS='SEQUENTIAL', FORM=FMT, IOSTAT=IER)
	    END IF
	ELSE
C         -- 'NEW' file
C  Note -- RECL=720 for mac, 2880 for sun
	    IF (ITYPE.EQ.2) THEN
	        OPEN (UNIT=UNIT,
     &              FILE=NAME(1:ICLEN(NAME)),
     &              STATUS='NEW',
C    &              RECORDTYPE='FIXED',
     &              ACCESS='DIRECT',
     &	             RECL=2880,
csun     &              RECL=720,
     &              FORM='UNFORMATTED',
     &              IOSTAT=IER)
	    ELSE
	        OPEN (UNIT=UNIT, FILE=NAME(1:ICLEN(NAME)), STATUS='NEW',
     &	             ACCESS='SEQUENTIAL', FORM=FMT,
C    2              CARRIAGECONTROL=CC,
     &	             IOSTAT=IER)
	    END IF
	END IF
C
C Success: find the complete file name.
C
	IF (IER.EQ.0) THEN
	    INQUIRE (UNIT=UNIT, NAME=RESULT)
	    VLBOPE = 1
C
C Failure: determine the error, and issue a message.
C
	ELSE
	    CNAME=NAME
	    IF (ISTAT.EQ.1) THEN
	        CTMP=TYPE
	        CALL PUTOUT('++ Cannot find '//CTMP//' file: '//CNAME)
	    ELSE
	        CTMP=TYPE
	        CALL PUTOUT('++ Cannot create '//CTMP//' file: '//CNAME)
	    END IF
	    !CALL GERROR(STRING)
	    CALL PUTOUT('++ '//STRING(1:LEN1(STRING)))
	    RESULT = STRING
	    VLBOPE = 0
	END IF
C-----------------------------------------------------------------------

c sws debugging
c	print *,'vlbope, file= [',file,'], result =[',result,']'

	END

	SUBROUTINE FITS_BSWAP (BUFFER,WORDS)
C
C     2 = NUMBER OF BYTES PER ELEMENT IN ARRAY PASSED TO BUFFER
C
	INTEGER WORDS
	integer*1 BUFFER(WORDS*2)  ! JLOU
C
C Exchanges the byte order in an array of 16 bit words
C so that it matches the convention used by non-DEC
C machines.  This is based on an MTPCKG routine written
C by Bill Sebok.
C
C Parameters -    (">" input, "!" modified, "<" output)
C;
C (!) BUFFER    (BYTE array BUFFER(WORDS*2)) The array
C 	     to be swapped.
C (>) WORDS    (INTEGER) The number of words in ARRAY
C
C Subroutines / functions used - None
C
C KS / CIT 11th Oct 1983
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
C
C     2 = NUMBER OF BYTES PER ARRAY ELEMENT PASSED TO BUFFER
C
	integer*1    TEMP  ! JLOU
	INTEGER*4    I

	DO 10, I = 2, WORDS*2, 2

	   TEMP = BUFFER(I - 1)
	   BUFFER(I - 1) = BUFFER(I)
	   BUFFER(I) = TEMP
 10      CONTINUE
	RETURN
	END


	SUBROUTINE FITS_FLUSH (STATUS)
C-----------------------------------------------------------------------
C     Flushes the main I/O buffer to disk.
C
C     Parameters -  (">" input, "W" workspace, "<" output, "!" modified)
C
C     (<) STATUS    (INTEGER) Returned status code. 0 => OK, non-zero
C                   codes are I/O error codes, which may be decoded
C                   by FITS_ERROR.
C
C     COMMON variables used -
C
C                   be filled with zeros and flushed to the output file.
C     (!) FPTR      (INTEGER) Pointer to next free byte in FBUFF.
C
C     Subroutines / functions used -
C
C     FIT_MWRIT     (FITS_ package) Writes FBUFF to output file.
C
C                                       KS / CIT 11th Oct 1983
C-----------------------------------------------------------------------
C
	INTEGER STATUS
C
	INCLUDE 'fits1.inc'
C
	INTEGER I
C
C     DO we have to flush FBUFF?
C
	IF (FPTR.GT.1) THEN
	   DO I=FPTR,2880
	      FBUFF(I)=0
	   END DO
	   CALL FITS_MWRIT(STATUS)
	   FPTR=1
	   IF (STATUS.NE.0)  GO TO 600
	END IF
C
  600 CONTINUE
	END


	SUBROUTINE FITS_MWRIT (STATUS)
C-----------------------------------------------------------------------
C     Writes the main FITS buffer, FBUFF, into the disk file.
C
C     Parameters -   (">" input, "<" output)
C
C     (<) STATUS     (INTEGER) Returned status code.  0 => OK,
C                    non-zero values are I/O error codes which can be
C                    decoded by FITS_ERROR.
C
C     COMMON variables used -
C
C     (>) MTUNIT     (INTEGER) The Fortran unit used for output.
C     (>) FBUFF      (Byte array) The main FITS buffer.
C
C                                             KS / CIT 9th Oct 1983
C-----------------------------------------------------------------------
C
	INTEGER STATUS
C
	INCLUDE 'fits1.inc'
C
	WRITE  (UNIT=MTUNIT, REC=BLKNO, IOSTAT=STATUS) FBUFF
	BLKNO = BLKNO + 1
C
	END SUBROUTINE FITS_MWRIT


C*PUTOUT -- WRITE line on standard output [Convex-UNIX]
C+
	SUBROUTINE PUTOUT(TEXT)
	CHARACTER*(*) TEXT
C
C This SUBROUTINE writes one line on the standard output; the text to be
C written is supplied as a CHARACTER-string argument.
C
C Argument:
C  TEXT   (input)  : CHARACTER string for output.
C
C Subroutines required:
C  Fortran formatted I/O.
C
C History:
C  1987 Nov 11 - TJP
C-----------------------------------------------------------------------
	WRITE (6,'(A)') TEXT
	END

C*LEN1 -- length of string excluding trailing blanks
C+
	INTEGER FUNCTION LEN1(S)
	CHARACTER*(*) S
C
C Find the length of a CHARACTER string excluding trailing blanks.
C A blank string returns a value of 0.
C
C Argument:
C  S      (input)  : CHARACTER string.
C
C Returns:
C  LEN1            : number of characters in S, excluding trailing
C                    blanks, in range 0...LEN(S). A blank string
C                    returns a value of 0.
C
C Subroutines required:
C  None
C
C Fortran 77 extensions:
C  None
C
C History:
C  1987 Nov 12 - TJP.
C-----------------------------------------------------------------------
	INTEGER  I
C
	IF (S.EQ.' ') THEN
	    LEN1 = 0
	ELSE
	    DO 10 I=LEN(S),1,-1
	        LEN1 = I
	        IF (S(I:I).NE.' ') GO TO 20
   10     CONTINUE
	    LEN1 = 0
   20     CONTINUE
	END IF
	END

C*ERROR -- PRINT error message on stdout, and abort [Convex-UNIX]
C+
	SUBROUTINE ERROR(STRING)
	CHARACTER*(*) STRING
	CHARACTER*24 CTMP
C
C This SUBROUTINE prints its argument (a CHARACTER string) on the
C standard output, prefixed with '++ERROR++', and THEN terminates
C execution of the PROGRAM, returning an error status to its parent
C (shell or command interpreter).
C
C Argument:
C  STRING (input)  : text of message
C
C Subroutines required:
C  EXIT
C  PUTOUT
C
C History:
C  1987 Apr 8 - TJP.
C-----------------------------------------------------------------------
	CTMP=STRING
	CALL PUTOUT('+++ERROR+++ '//CTMP)
	CALL EXIT(1)
C
	END

	SUBROUTINE FITS_SCALC (DATA,NELM,CHECK,VMIN,VMAX,SCALES,
     &                         ZEROS,ERRORS,
     &	                       image_row_size, image_col_size,
     &	                       max_image_col_size)
C-----------------------------------------------------------------------
C     Given a REAL array, this routine calculates the scale
C     and zero values needed to WRITE it out as a) 8 bit
C     unsigned integers, b) 16 bit signed integers, c) 32 bit
C     signed integers.  The user has to choose which of these
C     representations should be used, trading off precision
C     against compactness of DATA.
C
C     Parameters -   (">" input, "<" output)
C
C     (>) DATA     (REAL array DATA(NELM)) The DATA array.  This
C                  can have any number of dimensions, but is
C                  treated here as 1-dimensional for generality.
C     (>) NELM     (INTEGER) The number of elements in DATA.
C     (>) CHECK    (LOGICAL) IF true, the routine will pass through
C                  the DATA applying the scaling for each INTEGER
C                  representaion and THEN applying the inverse
C                  operation to the result.  The maximum errors
C                  found will be returned in ERRORS.
C     (<) VMIN     (REAL) Minimum value in DATA.
C     (<) VMAX     (REAL) Maximum value in DATA.
C     (<) SCALES   (REAL array SCALES(3)) Gives the
C                  scale factors for the three INTEGER types.
C                  SCALES(1) is for 8 bit, (2) for 16 bit, and
C                  (3) for 32 bit.
C     (<) ZEROS    (REAL array ZEROS(3)) Gives the zero values for
C                  the three INTEGER types.
C     (<) ERRORS   (REAL array ERRORS(3)) IF CHECK is specified,
C                  ERRORS will be returned with the maximum errors
C                  for each of the three representaions.  This may
C                  help in deciding that the loss in precision that
C                  be acceptable.
C     (>) IMAGE_ROW_SIZE
C                  (INTEGER) The actual row length of the image stored
C                  in DATA
C     (>) IMAGE_COL_SIZE
C                  (INTEGER) The actual column length of the image
C                  stored in DATA
C     (>) MAX_IMAGE_COL_SIZE
C                  (INTEGER) The maximum possible column length of the
C                  DATA array.
C
C     SCALES and ZEROS should be applied as -
C
C     INTEGER = ( REAL - ZEROS(N) ) / SCALES(N) + .5
C
C     COMMON variables used - None
C
C     Subroutines / functions used - None
C
C                                     KS / CIT 10th Oct 1983
C-----------------------------------------------------------------------
C
	LOGICAL CHECK
	INTEGER NELM
	SREAL DATA(NELM),VMIN,VMAX,SCALES(3),ZEROS(3),ERRORS(3)
c grw -- 1/18/94
	INTEGER image_row_size, image_col_size, max_image_col_size
C
	INTEGER i,INTVAL,j, k, l
	SREAL RANGE,RSCALE(3),RVALUE,VALUE
C
C     Find range of data array
C
	VMIN=DATA(1)
	VMAX=VMIN
c grw -- 1/18/94
	DO i = 0,
     &	    (image_row_size - 1) * max_image_col_size,
     &	     max_image_col_size
	   DO j=1, image_col_size
	      k = i + j
	      IF (DATA(k).GT.VMAX) VMAX=DATA(k)
	      IF (DATA(k).LT.VMIN) VMIN=DATA(k)
	   END DO
	END DO
C
C     Now get the three scales and zeros.  Note that we don't
C     try to make use of the VERY ends of the INTEGER ranges, to
C     allow a little slop for rounding errors etc.
C
	RANGE=VMAX-VMIN
	IF (RANGE.LT.1.E-35) THEN
	   SCALES(1)=1.
	ELSE
	   SCALES(1)=RANGE/252.
	END IF
	IF (RANGE.LT.1.E-30) THEN
	   SCALES(2)=1.
	ELSE
	   SCALES(2)=RANGE/65530.
	END IF
	IF (RANGE.LT.1.E-25) THEN
	   SCALES(3)=1.
	ELSE
	   SCALES(3)=RANGE/4.E9
	END IF
	ZEROS(1)=VMIN
	ZEROS(2)=(VMAX+VMIN)*.5
	ZEROS(3)=ZEROS(2)
C
C     Are we to check on the errors?
C
	IF (CHECK) THEN
	   DO J=1,3
	      ERRORS(J)=0.
	      RSCALE(J)=1./SCALES(J)
	   END DO
c grw -- 1/18/94
	   DO i = 0,
     &	       (image_row_size - 1) * max_image_col_size,
     &	        max_image_col_size
	      DO j=1, image_col_size
	         VALUE=DATA(i + j)
	         DO l=1,3
	            INTVAL=(VALUE-ZEROS(l))*RSCALE(l)+.5
	            RVALUE=(FLOAT(INTVAL)*SCALES(l))+ZEROS(l)
	            ERRORS(l)=MAX(ERRORS(l),ABS(VALUE-RVALUE))
	         END DO
	      END DO
	   END DO
	END IF
	END SUBROUTINE FITS_SCALC

C***********************************************************************

c grw -- 1/18/94
	SUBROUTINE FITS_WRAYF (DATA,NELM,BITPIX,BSCALE,BZERO,STATUS,
     &	                     IMAGE_ROW_SIZE, IMAGE_COL_SIZE,
     &	                     MAX_IMAGE_COL_SIZE)
C-----------------------------------------------------------------------
C     Writes a REAL array to output file as part of a FITS image,
C     converting it to one of a) 8 bit unsigned integers, b) 16 bit
C     signed integers,32 bit signed integers, by application of a scale
C     and zero value.
C
C     This routine simply copies the array element by element into a
C     buffer, writing the buffer out whenever it fills.  This means
C     that this routine may be called as many times as necessary - for
C     example, IF all the data is not available in memory at one time.
C     This routine should be called for the first time immediately
C     following the FITS_WEND CALL that terminated the FITS header.
C     Once all the data is output, FITS_CLOSE should be called to flush
C     the data buffer.
C
C     Parameters -   (">" input, "W" workspace, "!" modified, "<"
C     output)
C
C     (>) DATA     (REAL array DATA(NELM)) The DATA array to be output.
C                  DATA may be multi-dimensional, but it is treated
C                  as linear here for generality.
C     (>) NELM     (INTEGER) The number of elements in DATA.
C     (>) BITPIX   (INTEGER) Determines the number of bits per pixel
C                  to be used for the DATA.  Allowed values are
C                  8 (unsigned), 16 (signed), and 32 (signed).  IF
C                  BITPIX is none of these, 32 will be assumed.
C     (>) BSCALE   (REAL) The scale factor to apply to the DATA to
C                  convert it to integers.
C     (>) BZERO    (REAL) The offset to apply to the DATA to convert
C                  it to integers.  The integers are calculated as
C                  INTEGER = ( REAL - bzero ) / bscale + .5
C                  FITS_SCALC provides a means of calculating BSCALE
C                  and BZERO.
C     (<) STATUS   (INTEGER) Returned status code.  0 => OK, non-zero
C                  values indicate an I/O error and can be
C                  decoded using FITS_ERROR.
C grw -- 1/18/94
C     (>) IMAGE_ROW_SIZE
C                  (INTEGER) The actual row length of the image stored
C                  in DATA
C     (>) IMAGE_COL_SIZE
C                  (INTEGER) The actual column length of the image
C                  stored in DATA
C     (>) MAX_IMAGE_COL_SIZE
C                  (INTEGER) The maximum possible column length of the
C                  DATA array.
C
C     COMMON variables used -
C
C     (W) FBUFF    (Byte array) The main FITS I/O buffer
C     (W) FBUFFS   (INTEGER*2 array) The same buffer as FBUFF.
C     (W) FBUFFI   (INTEGER array) The same buffer as FBUFF.
C     (!) FPTR     (INTEGER) Byte level pointer to FBUFF.
C
C     Subroutines / functions used -
C
C     FITS_BSWAP    (FITS_ package) Swap order of bytes in words
C     FITS_WBSWAP   ( "      "   ) Swap order of bytes in longwords
C     FITS_MWRIT    (FITS_    "   ) WRITE FBUFF buffer to output file
C
C     Note: This routine assumes that the values for BSCALE and
C     BZERO will not generate DATA outside the INTEGER range
C     used.  IF this is not true, arithmetic faults will occur.
C
C                                       KS / CIT 19th Oct 1983
C-----------------------------------------------------------------------
C
	INTEGER NELM,BITPIX,STATUS
	SREAL DATA(NELM),BSCALE,BZERO
c grw -- 1/18/94
	INTEGER image_row_size, image_col_size, max_image_col_size
C
	INCLUDE 'fits1.inc'
C
	integer*1  IBYTE   ! JLOU
	INTEGER FPTRI,FPTRO,FPTRS,I,IVAL
	SREAL SCALE
C
C     The equivalence is to get round the problem that a VAX byte is
C     signed, and so cannot be assigned to an INTEGER in the range
C     0..255
C
	EQUIVALENCE (IBYTE,IVAL)
C
	SCALE=1./BSCALE
C
C     Operation depends on BITPIX
C
	IF (BITPIX.EQ.8) THEN
C
C        WRITE out as 8 bit unsigned integers.  Note the sneaky
C        use of IVAL and IBYTE.
C
c grw -- 1/18/94 modifying algorithm to handle variable dimensioned
c                fits images
C
	   DO i = 0,
     &	       (image_row_size - 1) * max_image_col_size,
     &	        max_image_col_size
	      DO j = 1, image_col_size

	         IVAL=(DATA(i + j)-BZERO)*SCALE+0.5
	         FBUFF(FPTR)=IBYTE
	         FPTR=FPTR+1
	         IF (FPTR.GT.2880) THEN
	            CALL FITS_MWRIT(STATUS)
	            FPTR=1
	            IF (STATUS.NE.0)  GO TO 600
	         END IF
	      END DO
	   END DO
C
	ELSE IF (BITPIX.EQ.16) THEN
C
C        WRITE out as signed 16 bit integers.  The CALL to FITS_BSWAP
C        is needed for DEC machines as byte order is not that required
C        by the FITS standard. Comment out for sun machines!!!
C
	   FPTRS=FPTR/2+1
	   FPTRO=FPTRS
c grw -- 1/18/94 modifying algorithm to handle variable dimensioned
c                fits images
C
	   DO i = 0,
     &	       (image_row_size - 1) * max_image_col_size,
     &	        max_image_col_size
	      DO j = 1, image_col_size

	         FBUFFS(FPTRS)=(DATA(i + j)-BZERO) *
     &	                      SCALE+0.5
	         FPTRS=FPTRS+1
	         IF (FPTRS.GT.1440) THEN
csun	            CALL FITS_BSWAP(FBUFFS(FPTRO),FPTRS-FPTRO)
	            CALL FITS_MWRIT(STATUS)
	            FPTR=1
	            FPTRS=1
	            FPTRO=1
	            IF (STATUS.NE.0) GO TO 600
	         END IF
	      END DO
	   END DO
csun	   IF (FPTRS.GT.1) CALL FITS_BSWAP(FBUFFS(FPTRO),FPTRS-FPTRO)
	   FPTR=(FPTRS-1)*2+1
C
	ELSE
C
C        WRITE out as signed 32 bit integers.  As for 16 bit integers
C        byte order is reversed for DEC machines using FITS_WBSWAP.
C
	   FPTRI=FPTR/4+1
	   FPTRO=FPTRI

c grw -- 1/18/94 modifying algorithm to handle variable dimensioned
c                fits images
C
	   DO i = 0,
     &	        (image_row_size - 1) * max_image_col_size,
     &	        max_image_col_size
	      DO j = 1, image_col_size

	         FBUFFI(FPTRI)=(DATA(i + j)-BZERO) *
     &	                      SCALE+0.5
	         FPTRI=FPTRI+1
	         IF (FPTRI.GT.720) THEN
csun	            CALL FITS_WBSWAP(FBUFFI(FPTRO),FPTRI-FPTRO)
	            CALL FITS_MWRIT(STATUS)
	            FPTR=1
	            FPTRI=1
	            FPTRO=1
	            IF (STATUS.NE.0)  GO TO 600
	         END IF
	      END DO
	   END DO
csun	   IF (FPTRI.GT.1) CALL FITS_WBSWAP(FBUFFI(FPTRO),FPTRI-FPTRO)
	   FPTR=(FPTRI-1)*4+1
	END IF
C
  600 CONTINUE
	END


	SUBROUTINE FITS_WBSWAP (BUFFER,LWORDS)
C
C     4 = NUMBER OF BYTES PER ELEMENT IN ARRAY PASSED TO BUFFER
C
	integer*1 BUFFER(LWORDS*4)  ! JLOU
	INTEGER LWORDS
C
C 	Exchanges the byte order in an array of 32 bit longwords
C       so that it matches the convention used by non-DEC
C 	machines.  This is based on an MTPCKG routine written
C 	by Bill Sebok.
C
C 	Parameters -    (">" input, "!" modified, "<" output)
C
C 	(!) BUFFER  (BYTE array BUFFER(LWORDS*4)) The array
C 		     to be swapped.
C 	(>) LWORDS   (INTEGER) The number of longwords in ARRAY
C
C 	Subroutines / functions used - None
C
C                                    KS / CIT 11th Oct 1983
C-----------------------------------------------------------------------
C
C     4 = NUMBER OF BYTES IN ARRAY ELEMENTS PASSED TO BUFFER
C
	integer*1    TEMP  ! JLOU
	INTEGER*4    I

	DO 10, I = 4, LWORDS*4, 4

c         WRITE(*,*) buffer(i-3),buffer(i-2),buffer(i-1),buffer(i)
	   TEMP = BUFFER(I - 3)
	   BUFFER(I - 3) = BUFFER(I)
	   BUFFER(I) = TEMP
	   TEMP = BUFFER(I - 2)
	   BUFFER(I - 2) = BUFFER(I - 1)
	   BUFFER(I - 1) = TEMP
c         WRITE(*,*) buffer(i-3),buffer(i-2),buffer(i-1),buffer(i)
 10      CONTINUE
	RETURN
	END

	SUBROUTINE FITS_WDBLE (NAME,VALUE,COMMENT,STATUS)
C-----------------------------------------------------------------------
C     Writes a DOUBLE PRECISION keyword into a FITS header.  Assumes
C     FITS_INIT has been called to initialise the header routines.
C
C     Parameters -  (">" input, "<" output)
C
C     (>) NAME      (CHARACTER) The name of the keyword.
C     (>) VALUE     (DOUBLE PRECISION) The value of the PARAMETER.
C     (>) COMMENT   (CHARACTER) Any comment to be associated with the
C                   keyword.  Should be no longer than 46 characters.
C     (>) STATUS    (INTEGER) RETURN status code. 0 => OK, non-zero
C                   values will be error codes caused by an
C                   I/O error.  These can be decoded by FITS_ERROR.
C
C     Note: The FITS standard says keyword names and values must be
C     output in upper case.  NAME can be passed in lower case, but
C     it will be output in upper.  COMMENT will be output as passed.
C
C     COMMON variables used -
C
C     (<) BUFF      (CHARACTER) Buffer for one header line
C
C     Functions / subroutines used -
C
C     FITS_BUFFO     (FITS_ package) Outputs BUFF to main I/O buffer.
C     FITS_BUFP      ( "     "    ) Copies NAME and COMMENTS into BUFF.
C
C                                                KS / CIT  9th OCT 1983
C-----------------------------------------------------------------------
C
	INTEGER STATUS
	DOUBLE PRECISION VALUE
	CHARACTER*(*) NAME,COMMENT
C
	INCLUDE 'fits1.inc'
C
	CALL FITS_BUFP(NAME,COMMENT)
	WRITE (BUFF(11:30),'(1PE20.13)',IOSTAT=STATUS) VALUE
	CALL FITS_BUFFO(STATUS)
C
	END


C***********************************************************************
C***********************************************************************

! 	SUBROUTINE RBININ(iPloTyp,var,row_sz,col_sz)

! C  Gets a matrix from a .mat file
! 	use param_mod
! #ifdef CSMACOS
!         use smacosio_mod
! #else
!         use macosio_mod
! #endif

! 	IMPLICIT NONE

! 	LOGICAL exist

! 	CHARACTER(len=MacosCharLen) :: file,fname,varname
! 	CHARACTER*82 SCR
! 	INTEGER iPloTyp,row_sz,col_sz,strLen,fnamelen
! 	SREAL var(mdttl,mdttl)

! C  get file name and variable name

!     1	CONTINUE

! 	CALL CACCEPT(file,'no default','Enter file name:')
! 	fnamelen=strLen(file)
! 	fname=file(1:fnamelen)

! 	INQUIRE (FILE=fname,EXIST=exist)
! 	IF (exist) THEN
! 	  CALL CACCEPT(file,'no default','Enter variable name:')
! 	  fnamelen=strLen(file)
! 	  varname=file(1:fnamelen)
! 	ELSE
! 	  SCR= '.mat file does not exist'
! 	  WRITE(*,*) SCR
! 	  GO TO 1
! 	END IF

! C	IF (iPloTyp.EQ.9) THEN

!  	    CALL RREAD_MAT_LOH(fname,varname,var,row_sz,col_sz)

! C	END IF

!  	RETURN
! 	END SUBROUTINE RBININ


C***********************************************************************

	Subroutine SysImging(mode,ca,obj,img,size)
	use macos_mod, only : CohImg,InCohImg
	Implicit none

        Integer :: mode,size,i,j,szArr(2)
	LREAL :: obj(size,size),obj2(size,size),
     &           objt(size,size),zmax,zmin
        Complex*16 :: ca(size,size),
     &    ca1(size,size),ca2(size,size),ca3(size,size)
        SREAL :: img(size,size),psf(size,size)

        ! Output image is stored in 'img' array

	szArr(1)=size; szArr(2)=size

	! normalize object image
	zmax=-1d-22; zmin=1d22
	!Do j=1,size/2
	  !Do i=1,size/2
	Do j=1,size
	  Do i=1,size
            If (obj(i,j)>zmax) zmax=obj(i,j)
            If (obj(i,j)<zmin) zmin=obj(i,j)
	  End Do
	End Do
	print*,'* obj: zmax, zmin =',zmax,zmin
	obj=(obj-zmin)/(zmax-zmin)

	!obj2(size/4+1:size*3/4,size/4+1:size*3/4)=64d0-obj
	!obj2(size/4+1:size*3/4,size/4+1:size*3/4)=1d0-obj
	!obj2=1d0-obj
	obj2=obj;

	 !print*,'obj2(1,1) =',obj2(1,1)
	 !print*,'obj2(128,128) =',obj2(128,128)

#if 1
	! invert image (top-down)
	Do j=1,size
	  objt(:,j)=obj2(:,size-j+1)
	End Do
	obj2=objt
#endif

	ca1=ca
	psf=ca1*DCONJG(ca1)

#if 0
	zmax=1d-22; zmin=1d22
        Do j=1,size
          Do i=1,size
            If (psf(i,j)>zmax) zmax=psf(i,j)
            If (psf(i,j)<zmin) zmin=psf(i,j)
          End Do
        End Do
	zmax=sqrt(zmax); zmin=sqrt(zmin)
	print*,'* psf: zmax, zmin =',zmax,zmin
	psf=(psf-zmin)/(zmax-zmin)
#endif

        If (mode==CohImg) Then
          ! img = |ca<*>obj|**2   where <*> is convolution opearator

	  ! fft ca1
          CALL DFOURN(ca1,szArr,2,1)

          ! generate ca of obj intensity and fft ca
          ca2=DCMPLX(sqrt(obj2))
          CALL DFOURN(ca2,szArr,2,1)

	  ca3=ca1*ca2
	  !ca3=ca1

	  CALL DFOURN(ca3,szArr,2,-1)
          CALL DSWAP2(ca2,ca3,size)
	  ca3=ca2/size**2
	  ca3=ca2

	  img=ca3*DCONJG(ca3)
	  !img=obj2;
	  print*,'- End of coherent imaging'

        Else If (mode==InCohImg) Then
          ! incoherent
          !psf = img
          ! img = psf<*>obj

	  ! generate ca of psf and fft the ca
	  ca1=DCMPLX(psf)
	  CALL DFOURN(ca1,szArr,2,1)

	  ! generate ca of obj intensity and fft the ca
          ca2=DCMPLX(obj2)
          CALL DFOURN(ca2,szArr,2,1)

	  ca3=ca1*ca2

	  CALL DFOURN(ca3,szArr,2,-1)
	  CALL DSWAP2(ca2,ca3,size)
	  ca3=ca2/size**2

	  img=sqrt(ca3*DCONJG(ca3))
	  print*,'- End of incoherent imaging'
        End If

        End Subroutine SysImging

C***********************************************************************
C*******END file sunsub.F***********************************************
C***********************************************************************
