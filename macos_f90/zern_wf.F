C***********************************************************************
C       Begin file zern_wf.F
C    +----------------------------------------------------------------+
C    |  Copyright (C) 1995-2007, California Institute of Technology.  |
C    |  U.S. Government Sponsorship Is Acknowledged.                  |
C    +----------------------------------------------------------------+
C***********************************************************************

C Zernike utility routines, adapted from Green's Matlab routines
C (Noll Circular Zernike Ordering)
C Added support for using NormHex and NormAnnular Zernike types
C
C zern_mode :      Extract a single Zernike mode from an OPD map using circular
C                  Zernike type NormNoll.
C zern_mode_nc :   Extract a single Zernike mode from an OPD map using hex or
C                  annular Zernike types.
C zern_fit:        Fit OPD with a set of given Zernike modes.
C zern_remove:     Remove given Zernike modes from an OPD map.
C zern_modes_coef: Compute and display Zernike modes coefficients.
C zern_aberrate:   Given an OPD and a set of user selected Zernike modes,
C                  return an OPD map containing those Zernike modes in
C                  the original OPD.
C
C Motivation:
C   1) Zernike mode based WF analysis, such as mode(s) fit and mode(s) removal
C   2) Supporting system optimization based on minimizing Zernike modes in WF
C
C John Z. Lou, Jet Propulsion Laboratory
C Last modified: 07/2008


	Module zern_wf_mod
	use param_mod
	use math_mod
	use lsq_mod

	use src_mod, only : nGridPts,RayIndex
   	use elt_mod, only : nRay,LRayPass,RayPos,ZernAnnuRatio,
     & 	                    ZernCent,ZernXGrid,ZernYGrid,ZernRad
	use macos_mod, only : zcoType,ZernCircular,ZernAnnular,ZernHex

#ifdef CMACOS
	use macosio_mod
#else
	use smacosio_mod
#endif

	Implicit none

	Integer, parameter :: max_zern_modes=90,
     &                        mzname=9  ! total zmode names supported
	Integer, parameter :: ZRMcmd=1,ZCOcmd=2,ZABcmd=3
        Integer :: ie_m,nzrm_modes_m,zrm_modes_m(max_zern_modes)
	Character(len=14), parameter :: zmode_name(9)=
     &    (/' Piston  ',' X-tilt  ',' Y-tilt  ','  Focus  ',' 45-Astig',
     &      ' XY-Astig',' X-coma  ',' Y-coma  ','Spherical'/)
        Real*8 :: zern_coefs_g(max_zern_modes),zRad_m
	Logical, dimension(:,:), allocatable :: mask_m
   	Real*8, dimension(:,:,:), allocatable :: OpdCoord_m
   	Logical :: usrOpt
   	! usrOpt is .true. if using user-specified parameters


	Contains


	! Extract Zernike modes from a given OPD map
	! Input: OPD, Zernike coefficients for each mode
	! Output Linear combination of Zernike modes of OPD: OPDzern
	!
	Subroutine zern_aberrate(OPD,OPDzern,mode,act)
        Real*8, dimension(:,:) :: OPD  ! nGridPts x nGridPts
        Real*8, dimension(:,:) :: OPDzern  ! output OPD w/ relevant Zernike terms
	Integer, dimension(:) :: mode
        Logical, dimension(size(OPD,1),size(OPD,2)) :: pmask
	Integer :: act,n_mode,iz

	OPDzern=0d0; pmask=(OPD/=0d0)
	Call zern_fit(OPD,pmask,mode,OPDzern,act)
	End Subroutine zern_aberrate


C  ***************************************************************************

C   Remove one or more Zernike modes from an OPD map
	Subroutine zern_remove(OPDin,OPDout,z_modes,act)
	Real*8, dimension(:,:) :: OPDin,OPDout
 	Integer, dimension(:) :: z_modes
        Integer, dimension(max_zern_modes) :: z_fit_modes
	Integer :: act,max_mode,m
	Logical, dimension(size(OPDin,1),size(OPDin,2)) :: pmask

	max_mode=maxval(z_modes); z_fit_modes=0
        Do m=1,max_mode
          z_fit_modes(m)=m
        End Do

        nzrm_modes_m=size(z_modes)
        zrm_modes_m(1:nzrm_modes_m)=z_modes

        if (.not.usrOpt) then
	  !print*,'In zern_remove(): size(z_modes)=',size(z_modes)
	  pmask=(OPDin/=0d0)
	else
          pmask(1:nGridPts,1:nGridPts)=mask_m(1:nGridPts,1:nGridPts)
        end if

	! OPDout below returns as least-square Zernike fit of all
        ! given Zernike modes
	Call zern_fit(OPDin,pmask,z_fit_modes(1:max_mode),OPDout,act)

	OPDout=OPDin-OPDout
	End Subroutine zern_remove

C  ***************************************************************************

	Subroutine zern_modes_coef(OPDin,z_modes,z_coefs,act)
        Real*8, dimension(:,:) :: OPDin
        Integer, dimension(:) :: z_modes
        Real*8, dimension(size(z_modes)) :: z_coefs
	Integer :: act
	Real*8, dimension(size(OPDin,1),size(OPDin,2)) :: OPDout
        Logical, dimension(size(OPDin,1),size(OPDin,2)) :: pmask

        if (.not.usrOpt) then
          pmask=(OPDin/=0d0)
	else
          pmask(1:nGridPts,1:nGridPts)=mask_m(1:nGridPts,1:nGridPts)
        end if
#if 0
        print*,'*** B4 calling zern_fit, size(OPDin,1) =',size(OPDin,1)
        print*,'*** B4 calling zern_fit, size(OPDin,2) =',size(OPDin,2)
        print*,'*** B4 calling zern_fit, z_modes =', z_modes
        print*,'*** B4 calling zern_fit, act =', act
#endif
        Call zern_fit(OPDin,pmask,z_modes,OPDout,act)
	z_coefs(1:size(z_modes))=zern_coefs_g(1:size(z_modes))
        End Subroutine zern_modes_coef

C  ***************************************************************************

C  Computes the Noll Zernike polynomial normalized for the support
C  of the specified pupil

	Function zern_mode(OPDin,z_mode)
	Real*8, dimension(:,:) :: OPDin
        Integer :: z_mode  ! single Zernike mode
	Real*8, dimension(size(OPDin,1),size(OPDin,2))
     &          :: zern_mode
!    &          :: zern_mode,r,a,R1
	Real*8 :: cent(2),x,y,r_tmp,rmax,q(4),Rk
	Integer, dimension(z_mode+3) :: zm, zn
	Real*8 :: p,v,pi_loc

	!Logical, dimension(size(OPDin,1),size(OPDin,2)) :: mask
	!Real*8, dimension(size(OPDin,1),size(OPDin,2)) :: dmask
        Logical, dimension(:,:), allocatable :: mask
        Real*8, dimension(:,:), allocatable :: a,r,R1,dmask
	Integer :: nrow,ncol, ir,ic, j,n,m,m1,s1,iRay,N1,N2

        If (.NOT. allocated(r)) Then
          N1=size(OPDin,1); N2=size(OPDin,2)
          allocate(r(N1,N2),a(N1,N2),R1(N1,N2),
     &             mask(N1,N2),dmask(N1,N2))
        End If

	pi_loc = 4d0*datan(1d0)

        nrow=nGridPts; ncol=nGridPts
	!nrow=size(OPDin,1); ncol=size(OPDin,2)
	!mask=(OPDin/=0d0)  ! mask out zero region

	If (.not. usrOpt) Then
	  mask=.false.
          Do iRay=2,nRay
            if (LRayPass(iRay)) then
              mask(RayIndex(1,iRay),RayIndex(2,iRay))=.true.
            end if
          End Do
	Else
          ! use user-specified Zernike center and axes
          mask(1:nrow,1:ncol)=mask_m(1:nrow,1:ncol)
	End If

	dmask=0d0
        Do ic=1,ncol
	  Do ir=1,nrow
	    if (mask(ir,ic)) dmask(ir,ic)=1d0
	  End Do
	End Do

	If (.not.usrOpt) Then
	  If (mod(nGridpts,2) /= 0) Then
            ! odd grid
            Call Centroid(dmask(1:nrow,1:ncol),nrow,ncol,
     &                    cent(1),cent(2))
          Else
            ! even grid, hardwire
            cent(1:2)=nGridpts/2+0.5d0
          End If
	  !print*,'cent =',cent(1:2)
	End If

	! Find maximum radial excusrion (rmax) from pupil centroid
	rmax=0d0
	If (.not.usrOpt) Then
          Do ic=1,ncol
	    y=ic; y=y-cent(2)
  	    Do ir=1,nrow
	      if (mask(ir,ic)) then
	        x=ir; x=x-cent(1)
	        r_tmp=dsqrt(x*x+y*y)
	        if (r_tmp>rmax) rmax=r_tmp
	      end if
	    End Do
	  End Do
	  !
          If (mod(nGridPts,2)==0) Then
            rmax=nGridPts/2-0.5  ! even grid, hardwire
          End If
        Else
	  ! use user-specified
          rmax = zRad_m
	End If
	!print*,'zern_mode: rmax=',rmax

	! Find corresponding n,m for specified Zernike mode
	j=1; zm=0; zn=0
	Do n=0,nint(dble(z_mode/2d0))
	  Do m1=0,n,2
	    if (mod(n,2) /= 0) then
	      m=m1+1
	    else
 	      m=m1
	    end if
	    !
	    zm(j)=m; zn(j)=n; j=j+1
	    if (m>0) then
	      zm(j)=m; zn(j)=n; j=j+1
	    end if
	    !
	    if (j .gt. z_mode) go to 20
	  End Do
	End Do
 20	m=zm(z_mode); n=zn(z_mode)

 	!print*,'zern_mode: j,m,n =',j,m,n

	! Compute Zernike mode over pupil support
	Do ic=1,ncol
	  if (.not.usrOpt) y=ic-cent(2)
	  Do ir=1,nrow
	    if (.not.usrOpt) then
	      x=ir-cent(1)
	    else
	      y = OpdCoord_m(2,ir,ic)
              x = OpdCoord_m(1,ir,ic)
	    end if
	    r(ir,ic) = dsqrt(x*x+y*y)/rmax
	    !a(ir,ic)=atan(y/(x+1d-9))
	    a(ir,ic)=atan(y/(x+1d-14))
	    q(1:4)=0d0
	    if (x>=0 .and. y>=0) q(1)=1
	    if (x<0 .and. y>=0)  q(2)=1
	    if (x<0 .and. y<0)   q(3)=1
	    if (x>=0 .and. y<0)  q(4)=1
	    a(ir,ic)= -1d0*(q(1)*a(ir,ic)+q(2)*(pi_loc+a(ir,ic))
     &                      +q(3)*(a(ir,ic)-pi_loc)+q(4)*a(ir,ic))
	  End Do
	End Do

	R1=0d0
	zern_mode=0d0  ! output OPD with the Zernike modes

	Do s1=0,(n-m)/2
	  Rk=PowerF(-1d0,s1)*dble(Factorial(n-s1))/dble(Factorial(s1))
     &       /dble(Factorial((n+m)/2-s1))/dble(Factorial((n-m)/2-s1))
	  Do ic=1,ncol
	    Do ir=1,nrow
	      R1(ir,ic)=R1(ir,ic)+Rk*PowerF(r(ir,ic),n-2*s1)
	    End Do
	  End Do
	End Do

	zern_mode=dsqrt(dble(n+1))*R1

	If (m>0) Then
	  If (mod(z_mode,2)==0) Then
            Do ic=1,ncol
	      Do ir=1,nrow
	        zern_mode(ir,ic)
     &            =zern_mode(ir,ic)*dcos(m*a(ir,ic))*dsqrt(2d0);
	      End Do
	    End Do
	  Else
            Do ic=1,ncol
	      Do ir=1,nrow
	    	zern_mode(ir,ic)
     &            =zern_mode(ir,ic)*dsin(m*a(ir,ic))*dsqrt(2d0);
	      End Do
            End Do
	  End If
	End If

        Do ic=1,ncol
	  Do ir=1,nrow
	    if (.not.mask(ir,ic)) zern_mode(ir,ic)=0d0
	  End Do
	End Do

	Call FindPV(zern_mode,nrow,ncol,p,v)
	!print*,'zern_mode, P,V=',p,v

        !deallocate(r,a,R1,mask,dmask)
	End Function zern_mode


C  ***************************************************************************

        ! For non-circular types
        FUNCTION zern_mode_nc(OPDin,z_mode,zt,ratio)
          use surfsub, only : ZerntoMon4, ZerntoMon7
          use elt_mod, only: mZernModes, mMonCoef

        Real*8, DIMENSION(:,:) :: OPDin
        Integer :: z_mode  ! single Zernike mode
        Real*8, DIMENSION(size(OPDin,1),size(OPDin,2)) :: zern_mode_nc
        Real*8, DIMENSION(size(OPDin,1),size(OPDin,2)) :: dmask
        Logical, DIMENSION(size(OPDin,1),size(OPDin,2)) :: mask
        !
        Real*8 :: ratio,zc(mZernModes),mc(mMonCoef),x,y,r,cent(2),r_tmp,rmax,
     &            xj,yi,yimj,p,v
        Integer :: zt,nrow,ncol, ir,ic,n,m,m1,s1,iRay,i,j,k
        Logical, save :: first_entry=.true.

        ncol=nGridPts; nrow=nGridPts
        zc=0d0; zc(z_mode)=1d0

        If (zt==ZernHex) Then
          Call ZerntoMon4(zc,mc)
        Else If (zt==ZernAnnular) Then
          Call ZerntoMon7(zc,mc,ratio)
        End If

	mask=.false.
        Do iRay=2,nRay
          if (LRayPass(iRay)) then
            mask(RayIndex(1,iRay),RayIndex(2,iRay))=.true.
          end if
        End Do
        !
        dmask=0d0
        Do ic=1,ncol
          Do ir=1,nrow
            IF (mask(ir,ic)) dmask(ir,ic)=1d0
          End Do
        End Do

        If (mod(nGridpts,2) /= 0) Then
          ! odd grid
          Call Centroid(dmask(1:nrow,1:ncol),nrow,ncol,cent(1),cent(2))
        Else
          ! even grid
          cent(1:2)=nGridpts/2+0.5d0
        End If

	DO ic=1,ncol
          y=ic; y=y-cent(2)
          DO ir=1,nrow
            x=ir; x=x-cent(1)
            IF (mask(ir,ic)) THEN
              r_tmp=dsqrt(x*x+y*y)
              IF (r_tmp>rmax) rmax=r_tmp
            END IF
          END DO
        END DO
        !
        If (mod(nGridPts,2)==0) Then
          rmax=nGridPts/2-0.5  ! even grid, test only
        End If

	Do ic=1,ncol
          y=ic; y=y-cent(2); y=y/rmax
          Do ir=1,nrow
            x=ir; x=x-cent(1); x=x/rmax
            k=1; yi=y
            zern_mode_nc(ir,ic)=mc(1)
            If ((x/=0d0).OR.(y/=0d0)) Then
              If (x==0d0) Then
                yi=y
                Do i=1,14
                  yimj=yi
                  Do j=0,i
                    k=k+1
                    If (mc(k)/=0d0) Then
                      If (j==0) Then
                        zern_mode_nc(ir,ic)
     &                       =zern_mode_nc(ir,ic)+mc(k)*yi
                      End If
                    End If
                    yimj=yimj/y
                  End Do
                  yi=yi*y
                End Do
              Else If (y==0d0) Then
	        Do i=1,14
                  xj=1
                  Do j=0,i
                    k=k+1
                    IF ((mc(k)/=0d0).AND.(i==j)) Then
                      zern_mode_nc(ir,ic)
     &                      =zern_mode_nc(ir,ic)+mc(k)*xj
                    End If
                    xj=xj*x
                  End Do
                End Do
              Else
	        Do i=1,14
                  xj=1
                  yimj=yi
                  Do j=0,i
                    k=k+1
                    If (mc(k).NE.0d0) Then
                      zern_mode_nc(ir,ic)
     &                =zern_mode_nc(ir,ic)+mc(k)*xj*yimj
                    End If
                    yimj=yimj/y
                    xj=xj*x
                  End Do
                  yi=yi*y
                End Do
              End If
            End If
          End Do
        End Do

	Do ic=1,ncol
          Do ir=1,nrow
            If (.not.mask(ir,ic)) zern_mode_nc(ir,ic)=0d0
          End Do
        End Do

#if 0
        Open(Unit=12,File='az.txt')
        Do ic=1,ncol
          Write(12,*) zern_mode_nc(1:nrow,ic)
        End Do
        Close(12)
#endif

        !Call FindPV(zern_mode_nc,nrow,ncol,p,v)
        !print*,'zern_mode_nc, P,V =',p,v

        first_entry=.false.
	End FUNCTION zern_mode_nc


C  ***************************************************************************
C Changed local arrays to dynamic arrays to avoid segmentation problem when
C compiled by Intel compiler ifort (on malfoy). -jzlou, 08/2014

	Subroutine zern_fit(OPDin,pupil_mask,modes,z_fit,act)
	Real*8, dimension(:,:) :: OPDin,z_fit
	Logical, dimension(:,:) :: pupil_mask
	Integer, dimension(:) :: modes
	Integer :: act

	Integer, parameter :: nmzline=6  ! modes displayed per line
	Real*8 :: coefs(size(modes))

!    &            ovec(size(OPDin,1)*size(OPDin,2)),
!    &            zn(size(OPDin,1),size(OPDin,2),size(modes))
        Real*8, dimension(:), allocatable :: ovec
        Real*8, dimension(:,:,:), allocatable :: zn

	Real*8, dimension(:,:), allocatable :: Cz,pCz

        !Integer, dimension(size(OPDin,1),size(OPDin,2)) :: yidx
        Integer, dimension(:), allocatable :: yidx
	Integer :: nrow,ncol,ir,ic,inz,k,strt_iz,end_iz,iz

        !print*,'** zern_fit: size(OPDin,1) =',size(OPDin,1)
        !print*,'** zern_fit: size(OPDin,2) =',size(OPDin,2)

        If (.NOT. allocated(yidx)) Then
          allocate(yidx(size(OPDin,1)*size(OPDin,2)),
     &             ovec(size(OPDin,1)*size(OPDin,2)),
     &             zn(size(OPDin,1),size(OPDin,2),size(modes)))
        End If

	nrow=size(OPDin,1); ncol=size(OPDin,2)

	! Extract non-zeros from OPDin
	inz=0; k=0
	Do ic=1,ncol
	  Do ir=1,nrow
	    k=k+1
	    if (OPDin(ir,ic)/=0) then
	      inz=inz+1; yidx(inz)=k
	      ovec(inz)=OPDin(ir,ic)
	    end if
	  End Do
	End Do

	if (inz==0) return  ! OPDin all zeros

        allocate(Cz(inz,size(modes)),pCz(size(modes),inz))
	!print*,'zern_fit: modes = ',modes(1:size(modes))

	! Assemble Cz matrix, each column of which is a Zernike mode
	Do iz=1,size(modes)
	  If (zcoType==ZernCircular) Then
            zn(:,:,iz)=zern_mode(OPDin,modes(iz))
          Else
            zn(:,:,iz)=zern_mode_nc(OPDin,modes(iz),zcoType,
     &                              ZernAnnuRatio(ie_m))
          End If
	  Do k=1,inz
	    ic=floor((yidx(k)-1)/DBLE(nrow))+1
	    ir=mod(yidx(k)-1,nrow)+1
	    Cz(k,iz)=zn(ir,ic,iz)
	  End Do
	End Do

#if 1
	! Pseudo inverse of Cz
	Call pinv(Cz,pCz)

	! Matrix arguments must have dims explicit for matmul intrinsic
	coefs=matmul(pCz(1:size(modes),1:inz),ovec(1:inz))
	zern_coefs_g(1:size(modes))=coefs(1:size(modes))

	If (act==ZCOcmd) Then
          if (act==2) WRITE(*,*)' '
          WRITE(*,*)' Zernike mode coeffs (OPD best fit):'

	  strt_iz=1; end_iz=min(strt_iz+nmzline-1,size(modes))
	  Do while (strt_iz.LE.size(modes))
	    ! Print mode names
	    WRITE(*,21)
 21         FORMAT('  ',$)
	    Do iz=strt_iz,end_iz
	      if (modes(iz)<=mzname-1) then
                WRITE(*,22) zmode_name(modes(iz))
	      else if (modes(iz)==11) then
                WRITE(*,22) zmode_name(mzname)
	      else
	        WRITE(*,23) modes(iz)
	      end if
 22           FORMAT(A,$)
!22           FORMAT(A,' ',$)
 23           FORMAT(I5,'         ',$)
	    End Do
            print*,' '
	    ! Print mode values
            Do iz=strt_iz,end_iz
              WRITE(*,24) coefs(iz)
 24           FORMAT(E14.6,$)
!24           FORMAT(E18.10,$)
            End Do
            print*,' '; print*,' '
	    strt_iz=end_iz+1
	    end_iz=min(strt_iz+nmzline-1,size(modes))
	  End Do  ! Do while loop
        End If

        z_fit=0d0
        If (act==1) Then
          ! ZRM
          Do iz=1,nzrm_modes_m
            z_fit=z_fit+coefs(zrm_modes_m(iz))*zn(:,:,zrm_modes_m(iz))
          End Do
        Else
          Do iz=1,size(modes)
            z_fit=z_fit+coefs(iz)*zn(:,:,iz)
            !print*,' **zern_fit: iz =',iz,', coef(iz) =',coefs(iz)
          End Do
        End If

	Do ic=1,ncol
          Do ir=1,nrow
	    if (.not.pupil_mask(ir,ic)) z_fit(ir,ic)=0d0
	  End Do
	End Do
#endif
	deallocate(Cz,pCz)
	!deallocate(yidx,ovec,zn,Cz,pCz)
	End Subroutine zern_fit

C ***************************************************************************

	Subroutine zern_beam_proc(rPos,zernCent,zernX,ZernY,zRad)
     	use math_mod, only : DDOTC

     	Real*8 :: rPos(:,:),zernCent(3),zernX(3),ZernY(3),zRad,dr(3)
     	Integer :: iRay
        Logical, save :: first_entry=.true.

	If (first_entry) Then
          if (allocated(mask_m)) deallocate(mask_m,OpdCoord_m)
          allocate(mask_m(mpts,mpts),OpdCoord_m(2,mpts,mpts))
          first_entry=.false.
        End If

	mask_m=.false.
        OpdCoord_m=1d0  ! avoid divide by zero
        Do iRay=2,nRay
          If (LRayPass(iRay)) Then
            mask_m(RayIndex(1,iRay),RayIndex(2,iRay))=.true.
            dr=rPos(1:3,iRay)-zernCent
            OpdCoord_m(1,RayIndex(1,iRay),
     &                 RayIndex(2,iRay))=DDOTC(dr,zernX)
            OpdCoord_m(2,RayIndex(1,iRay),
     &                 RayIndex(2,iRay))=DDOTC(dr,zernY)
          End If
        End Do
        zRad_m=zRad
   	End Subroutine zern_beam_proc

C ***************************************************************************

	Subroutine zern_interface(ie,OPDin,OPDout,wfe,usrOpt_in,
     &                            action)
	Real*8, dimension(:,:) :: OPDin,OPDout
	Real*8 :: wfe,pval,vval
	Integer :: ie,action, iz,jz,nzern,zerns(max_zern_modes)
	Logical :: usrOpt_in
	Real*8 :: zcoefs(max_zern_modes)

	ie_m=ie; usrOpt=usrOpt_in
        If (usrOpt) Then
          Call zern_beam_proc(RayPos,zernCent(1:3,ie),
     &          zernXGrid(1:3,ie),ZernYGrid(1:3,ie),zernRad(ie))
        End If

	! Get ZRM command arguments
 	CALL IACCEPT_S(nzern,1,1,
     &                 'Enter number of Zernike modes:')

	If (nzern>max_zern_modes) Then
	  WRITE(*,*)
     &    ' * Warning: Maximum Zernike modes exceeded!'
	  OPDout=OPDin
	  return
	Else If (nzern>0) Then
 	  CALL IACCEPT(zerns(1:nzern),zerns(1:nzern),-nzern,
     &                 'Enter Zernike indexes: ')
	Else If (nzern<0) Then
	  ! when nzern is negative, next integer value is starting
	  ! index of Zernike modes
	  nzern=-nzern
	  CALL IACCEPT(zerns(1:1),zerns(1:1),-1,
     &                 'Enter starting Zernike mode index: ')
	  WRITE(*,*)
     &      'Fitting OPD with Zernike modes ',zerns(1),' to ',
     &      zerns(1)+nzern-1
	  If (nzern>1) Then
	    jz=zerns(1)+1
	    Do iz=2,nzern
	      zerns(iz)=jz; jz=jz+1
	    End Do
	  End If
	Else
	  WRITE(*,*)' * Warning: Number of Zernike modes = 0'
	  return
	End If

	! Check input mode indices
	Do iz=1,nzern
	  If (zerns(iz)<1 .OR. zerns(iz)>max_zern_modes) Then
            WRITE(*,*)'Invalid Zernike mode: ',zerns(iz)
	    WRITE(*,*)'Zernike mode value must be between 1 and ',
     &                max_zern_modes
	    return
	  End If
	End Do

	If (action==ZRMcmd) Then
          Call zern_remove(OPDin,OPDout,zerns(1:nzern),action)
	Else If (action==ZCOcmd) Then
	  Call zern_modes_coef(OPDin,zerns(1:nzern),
     &                         zcoefs(1:nzern),action)
        Else If (action==ZABcmd) Then
	  Call zern_aberrate(OPDin,OPDout,zerns(1:nzern),action)
	End If
	wfe = RMS2(OPDout,size(OPDin,1),size(OPDin,2))

	!Call  FindPV(OPDout,size(OPDin,1),size(OPDin,2),pval,vval)
	!Call  FindPV(OPDin,size(OPDin,1),size(OPDin,2),pval,vval)
	!print*,'zern_interface: OPD P,V =',pval,vval
	End Subroutine zern_interface

	End Module zern_wf_mod
