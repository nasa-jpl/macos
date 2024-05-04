C***********************************************************************
C	Begin file matloh.F
C    +----------------------------------------------------------------+
C    |  Copyright (C) 1995-2007, California Institute of Technology.  |
C    |  U.S. Government Sponsorship Is Acknowledged.                  |
C    +----------------------------------------------------------------+
C***********************************************************************

C***********************************************************************
C
C    This is maintained as a separate .F, so that areas w/o access
C    to Matlab and it's library can use COMP.
C
C***********************************************************************

C***********************************************************************
C Routines in this file have been made dummy to avoid the need for
C linking matlab libraries - jzlou
C***********************************************************************

#include "realtype.h"

	! Subroutine RWRITE_MAT_LOH(image,row_sz,col_sz,brow_sz,bcol_sz,
    !  &  		          image_file,imag_file)
	! IMPLICIT NONE
	! CHARACTER*80 :: imag_file
    !     CHARACTER*(*) :: image_file
	! integer :: row_sz,col_sz, brow_sz,bcol_sz
	! SREAL :: image(brow_sz,bcol_sz)
	! return
	! End Subroutine RWRITE_MAT_LOH


	! SUBROUTINE RREAD_MAT_LOH(fname,varname,var,row_sz,col_sz)
	! use param_mod
    !     IMPLICIT NONE
    !     INTEGER row_sz,col_sz
    !     SREAL var(mdttl,mdttl)
    !     CHARACTER fname*(*)
    !     CHARACTER varname*(*)
	! return
	! End Subroutine RREAD_MAT_LOH


#if 0
	SUBROUTINE RWRITE_MAT_LOH(image,row_sz,col_sz,brow_sz,bcol_sz,
     &  image_file,imag_file)

C  matPutMatrix was used, instead of matPutFull because of problems with
C  the imaginary/real "flag" in Matlab.  Although
C  max(max(abs(imag(tr800))))=0,and mesh(tr800') works, mesh(tr800)
C  gives you the following error:
C  ??? Error using ==> surface, Argument must be real.

	use param_mod
	use cfiles_mod

	IMPLICIT NONE


	INTEGER matOpen,matClose,strLen
	INTEGER mxCreateFull,mxGetPr,matPutMatrix

	CHARACTER*80 mat_name,imag_file
	CHARACTER*(*) image_file
	INTEGER row_sz,col_sz,file_name_len,brow_sz,bcol_sz
	INTEGER fp,stat,i,j,k,a,sz
	SREAL image(brow_sz,bcol_sz)
	REAL*8 dimage(md2)

	k=0
	DO 99 j=1,col_sz
	  DO 98 i=1,row_sz
	    k=k+1
 	    dimage(k)=image(i,j)
   98	  CONTINUE
   99	CONTINUE

C  Open .mat file

	file_name_len = strLen(imag_file)
	fp=matOpen(imag_file(1:file_name_len),"u")
	file_name_len = strLen(image_file)

C  create location to put stuff...
C  0 means its a real matrix

	a=mxCreateFull(row_sz,col_sz,0)

C  set names and copy data appropriately
C	imgfil=filnam(1:j)//'.OPD'//Cinteger(1:i)

	sz=row_sz*col_sz
	CALL mxCopyReal8ToPtr(dimage,mxGetPr(a),sz)
	i=strLen(filnam)+2
	j=file_name_len-4
	IF (i.LT.j) THEN
	  mat_name=imgfil(i:j)
	ELSE
	  mat_name=filnam
	END IF
	CALL mxSetName(a,mat_name)
	stat=matPutMatrix(fp,a)

C  Close .mat file & free space

	stat=matClose(fp)
        CALL mxFreeMatrix(a)

	RETURN
	END

C***********************************************************************
C***********************************************************************

	SUBROUTINE RREAD_MAT_LOH(fname,varname,var,row_sz,col_sz)

C  fname=name of .mat file, .mat extension included
C  varname=name of variable to be gotten
	use param_mod
	IMPLICIT NONE


	INTEGER matOpen,matClose,matGetMatrix,mxGetPr,mxGetM,mxGetN
	INTEGER row_sz,col_sz

	INTEGER fp,stat,i,j,a,sz
C	SREAL var(row_sz,col_sz)
	SREAL var(mdttl,mdttl)
	REAL*8 dimage(md2)
	CHARACTER fname*(*)
	CHARACTER varname*(*)

C  Open .mat file

	fp=matOpen(fname,"r")

C  Get matrix

	a=matGetMatrix(fp,varname)
	row_sz=mxGetM(a)
	col_sz=mxGetN(a)
	sz=row_sz*col_sz
	CALL mxCopyPtrToReal8(mxGetPr(a),dimage,sz)

C  Close .mat file & free space

	stat=matClose(fp)
        CALL mxFreeMatrix(a)

C  Fortran stores as (1,1)(2,1)(3,1)(1,2)(2,2)(3,2)
C  Matlab stores as (1,1)(2,1)(3,1)(1,2)(2,2)(3,2)

	DO 98 j=1,col_sz
	  DO 99 i=1,row_sz
 	    var(i,j)=dimage((j-1)*row_sz+i)
   99	  CONTINUE
   98	CONTINUE


	RETURN
	END
#endif
C***********************************************************************
