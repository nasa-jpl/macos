C***********************************************************************
C       Begin file resamp.F
C    +----------------------------------------------------------------+
C    |  Copyright (C) 1995-2007, California Institute of Technology.  |
C    |  U.S. Government Sponsorship Is Acknowledged.                  |
C    +----------------------------------------------------------------+
C***********************************************************************
C
C utility routines for
C 1) Image plane intensity field resampling and combining 
C 2) Image plane complex electric field resampling and combining
C 3) Multi-path, interferometric systems modeling 
C Last updated: 02/2008
C

	    Module resamp_mod
	    use elt_mod

	    Implicit none
	    Save

	    Integer, parameter :: ToIF=1,SvRx=2,ToArm1=3,ToArm2=4,
     &                            RestoreRx=5
	    Integer :: IF_pass,iEltRx


	    Contains

#if 0
	    Subroutine IF_state_copy(dir)
	    Integer :: dir, i,j,k,iElt

	    If (dir==ToIF) Then 
              ! Copy current elements state to IF storages

	      ! Save Rx configuration
              Do iElt=1,nElt
                EltID_IF(iElt,3)=EltID(iElt)
                SrfType_IF(iElt,3)=SrfType(iElt)
                KrElt_IF(iElt,3)=KrElt(iElt)
                KcElt_IF(iElt,3)=KcElt(iElt)
                psiElt_IF(1:3,iElt,3)=psiElt(1:3,iElt)
                vptElt_IF(1:3,iElt,3)=vptElt(1:3,iElt)
                rptElt_IF(1:3,iElt,3)=rptElt(1:3,iElt)
                pMon_IF(1:3,iElt,3)=pMon(1:3,iElt)
                xMon_IF(1:3,iElt,3)=xMon(1:3,iElt)
                yMon_IF(1:3,iElt,3)=yMon(1:3,iElt)
                zMon_IF(1:3,iElt,3)=zMon(1:3,iElt)
                ApVec_IF(1:3,iElt,3)=ApVec(1:3,iElt)
                ObsVec_IF(1:3,1,iElt,3)=ObsVec(1:3,1,iElt)
                xObs_IF(1:3,iElt,3)=xObs(1:3,iElt)
                EltCoat_IF(iElt,3)=EltCoat(iElt)
                nObs_IF(iElt,3)=nObs(iElt)
                PropType_IF(iElt,3)=PropType(iElt)
	        zElt_IF(iElt,3)=zElt(iElt)
	        !GlassCoef_IF(1:6,3)=GlassCoef(1:6,iElt)
                !ZernCoef(45,mElt)=
              End Do
	      Do iElt=0,nElt
	        IndRef_IF(iElt,3)=IndRef(iElt) 
	        Extinc_IF(iElt,3)=Extinc(iElt) 
	        zElt_IF(iElt,3)=zElt(iElt)
	      End Do

	      ! Construct first-arm configuration
              Do iElt=1,EltBranch_IF(2,1)
	        Rx_iElt(iElt,1)=iElt
                EltID_IF(iElt,1)=EltID(iElt)
                SrfType_IF(iElt,1)=SrfType(iElt)
                KrElt_IF(iElt,1)=KrElt(iElt)
                KcElt_IF(iElt,1)=KcElt(iElt)
                psiElt_IF(1:3,iElt,1)=psiElt(1:3,iElt)
                vptElt_IF(1:3,iElt,1)=vptElt(1:3,iElt)
                rptElt_IF(1:3,iElt,1)=rptElt(1:3,iElt)
                pMon_IF(1:3,iElt,1)=pMon(1:3,iElt)
                xMon_IF(1:3,iElt,1)=xMon(1:3,iElt)
                yMon_IF(1:3,iElt,1)=yMon(1:3,iElt)
                zMon_IF(1:3,iElt,1)=zMon(1:3,iElt)
                ApVec_IF(1:3,iElt,1)=ApVec(1:3,iElt)
                ObsVec_IF(1:3,1,iElt,1)=ObsVec(1:3,1,iElt)
                xObs_IF(1:3,iElt,1)=xObs(1:3,iElt)
                EltCoat_IF(iElt,1)=EltCoat(iElt)
                nObs_IF(iElt,1)=nObs(iElt)
                PropType_IF(iElt,1)=PropType(iElt)
              End Do
	      Do iElt=0,EltBranch_IF(2,1)
	        IndRef_IF(iElt,1)=IndRef(iElt)
	        Extinc_IF(iElt,1)=Extinc(iElt)
	        zElt_IF(iElt,1)=zElt(iElt)
	      End Do
	      j=EltBranch_IF(2,1)+1
              Do iElt=EltBranch_IF(2,2)+1,nElt
	        Rx_iElt(j,1)=iElt
                EltID_IF(j,1)=EltID(iElt)
                SrfType_IF(j,1)=SrfType(iElt)
                KrElt_IF(j,1)=KrElt(iElt)
                KcElt_IF(j,1)=KcElt(iElt)
                psiElt_IF(1:3,j,1)=psiElt(1:3,iElt)
                vptElt_IF(1:3,j,1)=vptElt(1:3,iElt)
                rptElt_IF(1:3,j,1)=rptElt(1:3,iElt)
	        IndRef_IF(j,1)=IndRef(iElt)
	        Extinc_IF(j,1)=Extinc(iElt)
                pMon_IF(1:3,j,1)=pMon(1:3,iElt)
                xMon_IF(1:3,j,1)=xMon(1:3,iElt)
                yMon_IF(1:3,j,1)=yMon(1:3,iElt)
                zMon_IF(1:3,j,1)=zMon(1:3,iElt)
                ApVec_IF(1:3,j,1)=ApVec(1:3,iElt)
                ObsVec_IF(1:3,1,j,1)=ObsVec(1:3,1,iElt)
                xObs_IF(1:3,j,1)=xObs(1:3,iElt)
                EltCoat_IF(j,1)=EltCoat(iElt)
                nObs_IF(j,1)=nObs(iElt)
                PropType_IF(j,1)=PropType(iElt)
                j=j+1
              End Do

              ! Construct second-arm configuration
              Do iElt=1,EltBranch_IF(1,1)-1
	        Rx_iElt(iElt,2)=iElt
                EltID_IF(iElt,2)=EltID(iElt)
                SrfType_IF(iElt,2)=SrfType(iElt)
                KrElt_IF(iElt,2)=KrElt(iElt)
                KcElt_IF(iElt,2)=KcElt(iElt)
                psiElt_IF(1:3,iElt,2)=psiElt(1:3,iElt)
                vptElt_IF(1:3,iElt,2)=vptElt(1:3,iElt)
                rptElt_IF(1:3,iElt,2)=rptElt(1:3,iElt)
                pMon_IF(1:3,iElt,2)=pMon(1:3,iElt)
                xMon_IF(1:3,iElt,2)=xMon(1:3,iElt)
                yMon_IF(1:3,iElt,2)=yMon(1:3,iElt)
                zMon_IF(1:3,iElt,2)=zMon(1:3,iElt)
                ApVec_IF(1:3,iElt,2)=ApVec(1:3,iElt)
                ObsVec_IF(1:3,1,iElt,2)=ObsVec(1:3,1,iElt)
                xObs_IF(1:3,iElt,2)=xObs(1:3,iElt)
                EltCoat_IF(iElt,2)=EltCoat(iElt)
                nObs_IF(iElt,2)=nObs(iElt)
                PropType_IF(iElt,2)=PropType(iElt)
              End Do
	      Do iElt=0,EltBranch_IF(1,1)-1
	        IndRef_IF(iElt,2)=IndRef(iElt)
                Extinc_IF(iElt,2)=Extinc(iElt)
	        zElt_IF(iElt,2)=zElt(iElt)
	      End Do
              j=EltBranch_IF(1,1)
              Do iElt=EltBranch_IF(1,2),nElt
	        Rx_iElt(j,2)=iElt
                EltID_IF(j,2)=EltID(iElt)
                SrfType_IF(j,2)=SrfType(iElt)
                KrElt_IF(j,2)=KrElt(iElt)
                KcElt_IF(j,2)=KcElt(iElt)
                psiElt_IF(1:3,j,2)=psiElt(1:3,iElt)
                vptElt_IF(1:3,j,2)=vptElt(1:3,iElt)
                rptElt_IF(1:3,j,2)=rptElt(1:3,iElt)
	        IndRef_IF(j,2)=IndRef(iElt)
	        Extinc_IF(j,2)=Extinc(iElt)
                pMon_IF(1:3,j,2)=pMon(1:3,iElt)
                xMon_IF(1:3,j,2)=xMon(1:3,iElt)
                yMon_IF(1:3,j,2)=yMon(1:3,iElt)
                zMon_IF(1:3,j,2)=zMon(1:3,iElt)
                ApVec_IF(1:3,j,2)=ApVec(1:3,iElt)
                ObsVec_IF(1:3,1,j,2)=ObsVec(1:3,1,iElt)
                xObs_IF(1:3,j,2)=xObs(1:3,iElt)
                EltCoat_IF(j,2)=EltCoat(iElt)
                nObs_IF(j,2)=nObs(iElt)
                PropType_IF(j,2)=PropType(iElt)
                j=j+1
              End Do

	    Else If (dir==SvRx) Then
	      k=3
  	      Do iElt=1,nElt_IF(3)
                EltID(iElt)=EltID_IF(iElt,k)
                SrfType(iElt)=SrfType_IF(iElt,k)
                KrElt(iElt)=KrElt_IF(iElt,k)
          	KcElt(iElt)=KcElt_IF(iElt,k)
          	psiElt(1:3,iElt)=psiElt_IF(1:3,iElt,k)
          	vptElt(1:3,iElt)=vptElt_IF(1:3,iElt,k)
          	rptElt(1:3,iElt)=rptElt_IF(1:3,iElt,k)
          	pMon(1:3,iElt)=pMon_IF(1:3,iElt,k)
          	xMon(1:3,iElt)=xMon_IF(1:3,iElt,k)
          	yMon(1:3,iElt)=yMon_IF(1:3,iElt,k)
          	zMon(1:3,iElt)=zMon_IF(1:3,iElt,k)
          	ApVec(1:3,iElt)=ApVec_IF(1:3,iElt,k)
          	ObsVec(1:3,1,iElt)=ObsVec_IF(1:3,1,iElt,k)
          	xObs(1:3,iElt)=xObs_IF(1:3,iElt,k)
          	EltCoat(iElt)=EltCoat_IF(iElt,k)
          	nObs(iElt)=nObs_IF(iElt,k)
          	PropType(iElt)=PropType_IF(iElt,k)
              End Do
	      Do iElt=0,nElt_IF(3)
          	IndRef(iElt)=IndRef_IF(iElt,k)
	        Extinc(iElt)=Extinc_IF(iElt,k)
	        zElt(iElt)=zElt_IF(iElt,k)
	      End Do
	      
	    Else If (dir==ToArm1) Then 
	      ! Copy first-arm into configuration
              k=1
              Do iElt=1,nElt_IF(1)
                EltID(iElt)=EltID_IF(iElt,k)
                SrfType(iElt)=SrfType_IF(iElt,k)
                KrElt(iElt)=KrElt_IF(iElt,k)
                KcElt(iElt)=KcElt_IF(iElt,k)
                psiElt(1:3,iElt)=psiElt_IF(1:3,iElt,k)
                vptElt(1:3,iElt)=vptElt_IF(1:3,iElt,k)
                rptElt(1:3,iElt)=rptElt_IF(1:3,iElt,k)
                pMon(1:3,iElt)=pMon_IF(1:3,iElt,k)
                xMon(1:3,iElt)=xMon_IF(1:3,iElt,k)
                yMon(1:3,iElt)=yMon_IF(1:3,iElt,k)
                zMon(1:3,iElt)=zMon_IF(1:3,iElt,k)
                ApVec(1:3,iElt)=ApVec_IF(1:3,iElt,k)
                ObsVec(1:3,1,iElt)=ObsVec_IF(1:3,1,iElt,k)
                xObs(1:3,iElt)=xObs_IF(1:3,iElt,k)
                EltCoat(iElt)=EltCoat_IF(iElt,k)
                nObs(iElt)=nObs_IF(iElt,k)
                PropType(iElt)=PropType_IF(iElt,k)
              End Do
	      Do iElt=0,nElt_IF(1)
                IndRef(iElt)=IndRef_IF(iElt,k)
	        Extinc(iElt)=Extinc_IF(iElt,k)
	        zElt(iElt)=zElt_IF(iElt,k)
	      End Do
              nElt=nElt_IF(1)

            Else If (dir==ToArm2) Then
	      ! Copy second-arm into configuration
              k=2
              Do iElt=1,nElt_IF(2)
                EltID(iElt)=EltID_IF(iElt,k)
                SrfType(iElt)=SrfType_IF(iElt,k)
                KrElt(iElt)=KrElt_IF(iElt,k)
                KcElt(iElt)=KcElt_IF(iElt,k)
                psiElt(1:3,iElt)=psiElt_IF(1:3,iElt,k)
                vptElt(1:3,iElt)=vptElt_IF(1:3,iElt,k)
                rptElt(1:3,iElt)=rptElt_IF(1:3,iElt,k)
                IndRef(iElt)=IndRef_IF(iElt,k)
	        Extinc(iElt)=Extinc_IF(iElt,k)
                pMon(1:3,iElt)=pMon_IF(1:3,iElt,k)
                xMon(1:3,iElt)=xMon_IF(1:3,iElt,k)
                yMon(1:3,iElt)=yMon_IF(1:3,iElt,k)
                zMon(1:3,iElt)=zMon_IF(1:3,iElt,k)
                ApVec(1:3,iElt)=ApVec_IF(1:3,iElt,k)
                ObsVec(1:3,1,iElt)=ObsVec_IF(1:3,1,iElt,k)
                xObs(1:3,iElt)=xObs_IF(1:3,iElt,k)
                EltCoat(iElt)=EltCoat_IF(iElt,k)
                nObs(iElt)=nObs_IF(iElt,k)
                PropType(iElt)=PropType_IF(iElt,k)
              End Do
	      Do iElt=0,nElt_IF(2)
                IndRef(iElt)=IndRef_IF(iElt,k)
                Extinc(iElt)=Extinc_IF(iElt,k)
                zElt(iElt)=zElt_IF(iElt,k)
              End Do
              nElt=nElt_IF(2)

	    Else if (dir==RestoreRx) Then
	      ! Restore Rx elements order
              k=3
              Do iElt=1,nElt_IF(3)
          	EltID(iElt)=EltID_IF(iElt,k)
          	SrfType(iElt)=SrfType_IF(iElt,k)
          	KrElt(iElt)=KrElt_IF(iElt,k)
          	KcElt(iElt)=KcElt_IF(iElt,k)
          	psiElt(1:3,iElt)=psiElt_IF(1:3,iElt,k)
          	vptElt(1:3,iElt)=vptElt_IF(1:3,iElt,k)
          	rptElt(1:3,iElt)=rptElt_IF(1:3,iElt,k)
          	pMon(1:3,iElt)=pMon_IF(1:3,iElt,k)
          	xMon(1:3,iElt)=xMon_IF(1:3,iElt,k)
          	yMon(1:3,iElt)=yMon_IF(1:3,iElt,k)
          	zMon(1:3,iElt)=zMon_IF(1:3,iElt,k)
          	ApVec(1:3,iElt)=ApVec_IF(1:3,iElt,k)
          	ObsVec(1:3,1,iElt)=ObsVec_IF(1:3,1,iElt,k)
          	xObs(1:3,iElt)=xObs_IF(1:3,iElt,k)
          	EltCoat(iElt)=EltCoat_IF(iElt,k)
          	nObs(iElt)=nObs_IF(iElt,k)
          	PropType(iElt)=PropType_IF(iElt,k)
              End Do
	      Do iElt=0,nElt_IF(3)
	        IndRef(iElt)=IndRef_IF(iElt,k)
	        Extinc(iElt)=Extinc_IF(iElt,k)
	        zElt(iElt)=zElt_IF(iElt,k)
	      End Do
              nElt=nElt_IF(3)
	    End If
	    End Subroutine IF_state_copy

C  ***********************************************************************

	    ! iElt argument is 'global' from Rx order
	    Function get_IF_section(iElt)
	    Integer :: iElt, get_IF_section

	    if (iElt<EltBranch_IF(1,1)) then
	      get_IF_section=1
	    else if (iElt .ge. EltBranch_IF(1,1) .and.
     &               iElt .le. EltBranch_IF(2,1)) then
	      get_IF_section=2
	    else if (iElt .ge. EltBranch_IF(1,2) .and.
     &               iElt .le. EltBranch_IF(2,2)) then 
	      get_IF_section=3
	    else
	      ! last section, where interference occurs
	      get_IF_section=4
	    end if
	    End Function get_IF_section


	    Function get_branch_iElt(iEltg,br)
	    Integer :: get_branch_iElt,iEltg,br,ie

	    Do ie=1,nElt_IF(br)
	      if (Rx_iElt(ie,br)==iEltg) then
	        get_branch_iElt=ie
	        return
	      end if
	    End Do
	    End Function get_branch_iElt

#endif

C  ***********************************************************************

C   Test_mGridData: for testing mGridData routine only
          Subroutine Test_mGridData
	  Implicit none

          Integer, parameter :: ndim=128  ! dim of reg. grid
	  Integer :: n_mesh_points,i
          Real*8, dimension(:), allocatable :: mesh_node_coords,
     &                                         irrdata
          Real*8 :: xcoord(ndim),ycoord(ndim),regdata(ndim,ndim)
          Real*8 :: x0,y0

	  ! Get irregular mesh grid point coords
          Open(unit=10,file='SurfMap011.txt',status='old')
          read(10,*) n_mesh_points

          allocate(mesh_node_coords(2*n_mesh_points),
     &             irrdata(n_mesh_points))

          Do i=1,n_mesh_points
            read(10,*)mesh_node_coords(2*i-1),mesh_node_coords(2*i),x0,
     &                irrdata(i)
          End Do
          close(10)

	  ! On input, xcoord and ycoord provide dimensions of regular
	  ! grid in X and Y; on output, they contain coordinates of
	  ! regular grid points, and regdata contains resampled data.
	  Call mGridData(n_mesh_points,mesh_node_coords,irrdata,
     &                   xcoord,ycoord,regdata)

	  ! Write out data on regular grid
          Open(unit=11,file='regdata.txt',status='replace')
          Do i=1,ndim
            write(11,*) regdata(i,1:ndim)
          End Do
          close(11);

	  deallocate(mesh_node_coords,irrdata)
	  End Subroutine Test_mGridData

C  ***********************************************************************

C
C  Subroutine mGridData fits data defined on a set of random points
C  to a regular grid. 
C  Arguments 'n_mesh_points,mesh_node_coords,idata' define the data set
C  on the random points, and 'xdim,ydim,x0,y0,dx,dy,rdata' define the 
C  data on a regular grid.
C

	  Subroutine mGridData(n_mesh_points,mesh_node_coords,idata,
     &                         xcrd,ycrd,rdata)
	  Integer :: n_mesh_points
          Real*8, dimension(:) :: mesh_node_coords,xcrd,ycrd,idata
	  Real*8, dimension(:,:) :: rdata 

	  Real*8 :: x0,y0,dx,dy
          Integer :: max_mesh_elts,n_mesh_elts,xdim,ydim,i
	  Integer, dimension(:), allocatable :: mesh_elt_nodes

	  max_mesh_elts=3*n_mesh_points
	  allocate(mesh_elt_nodes(3*max_mesh_elts))

	  ! Define a rectangular grid covering the irregular grid
	  xdim=size(xcrd); ydim=size(ycrd)
          Call Define_rect_grid(n_mesh_points,mesh_node_coords,
     &                          xdim,ydim,x0,y0,dx,dy)

	  ! Define xcrd and ycrd, coordinates of regular grid 
	  Do i=1,xdim
	    xcrd(i)=x0+(i-1)*dx
	  End Do
	  Do i=1,ydim
            ycrd(i)=y0+(i-1)*dy
          End Do

	  ! Create a triangular mesh, calling qhull C library
          Call meshgen(n_mesh_points,mesh_node_coords,
     &                 n_mesh_elts,mesh_elt_nodes)

	  ! meshgen start node id from 0, buf in f90 it starts from 1
	  mesh_elt_nodes=mesh_elt_nodes+1

	  ! Interpolation from irr grid to reg grid
          Call Interp_over_trimesh(n_mesh_points,
     &                mesh_node_coords(1:2*n_mesh_points),
     &                n_mesh_elts,mesh_elt_nodes(1:3*n_mesh_elts),
     &                xdim,ydim,x0,y0,dx,dy,idata,rdata)

	  deallocate(mesh_elt_nodes)
          End Subroutine mGridData  

C  ***********************************************************************

C  Generate a regular grid based on a given set of random points
            Subroutine Define_rect_grid(n_nodes,mesh_coords,ndx,ndy,
     &                        xstrt,ystrt,dx,dy)
            Integer :: n_nodes,ndx,ndy,n
            Real*8, dimension(:) :: mesh_coords
            Real*8 :: xstrt,ystrt,dx,dy,xmin,xmax,ymin,ymax,span

            xmin=1d14; xmax=-1d14; ymin=1d14; ymax=-1d14
            Do n=1,n_nodes
              If (mesh_coords(2*n-1)<xmin) xmin=mesh_coords(2*n-1)
              If (mesh_coords(2*n-1)>xmax) xmax=mesh_coords(2*n-1)
              If (mesh_coords(2*n)<ymin) ymin=mesh_coords(2*n)
              If (mesh_coords(2*n)>ymax) ymax=mesh_coords(2*n)
            End Do
            xstrt=xmin; ystrt=ymin

            span=xmax-xmin; dx=span/(ndx-1)
            span=ymax-ymin; dy=span/(ndy-1)
            End Subroutine Define_rect_grid

C  ***********************************************************************

C    This is the optical beam combination routine.
C    It takes two complex fields defined on irregular grids, resample
C    them to a regular grid on their overlaping region, then combine
C    (add) them up to produce a single complex field stored in 'regdata'.
C
	    Subroutine Combining_cmplx_fields(n_fd1_points,fd1_node_coords,
     &                                        n_fd2_points,fd2_node_coords,
     &                                        irr_fd1_data,irr_fd2_data,
     &                         		      xcrd,ycrd,regdata)
	    Integer :: n_fd1_points,fd1,n_fd2_points
	    Real*8, dimension(:) :: fd1_node_coords,fd2_node_coords,
     &                              xcrd,ycrd
	    Complex*16, dimension(:) :: irr_fd1_data,irr_fd2_data
	    Complex*16, dimension(:,:) :: regdata

	    Integer, dimension(:,:), allocatable :: mesh_elt_nodes
	    Integer :: max_mesh_elts,xdim,ydim,n_mesh_elts(2),i
	    Real*8 :: x0,y0,dx,dy

	    max_mesh_elts=max(3*n_fd1_points,3*n_fd2_points)
	    allocate(mesh_elt_nodes(3*max_mesh_elts,2))

	    ! Define a rectangular grid covering the intersection of 
            ! the two irregular grids
            xdim=size(xcrd); ydim=size(ycrd)
	    Call Define_overlap_rect_grid(n_fd1_points,fd1_node_coords,
     &                                    n_fd2_points,fd2_node_coords,
     &					  xdim,ydim,x0,y0,dx,dy)

	    ! Define xcrd and ycrd, coordinates of regular grid
            Do i=1,xdim
              xcrd(i)=x0+(i-1)*dx
            End Do
            Do i=1,ydim
              ycrd(i)=y0+(i-1)*dy
            End Do

	    ! Create first triangular mesh, calling qhull C library
            Call meshgen(n_fd1_points,fd1_node_coords,
     &                   n_mesh_elts(1),mesh_elt_nodes(:,1))
	    !
	    ! meshgen start node id from 0, buf in f90 it starts from 1
            mesh_elt_nodes(:,1)=mesh_elt_nodes(:,1)+1

	    ! Create second triangular mesh, calling qhull C library
            Call meshgen(n_fd2_points,fd2_node_coords,
     &                   n_mesh_elts(2),mesh_elt_nodes(:,2))
            !
            ! meshgen start node id from 0, buf in f90 it starts from 1
            mesh_elt_nodes(:,2)=mesh_elt_nodes(:,2)+1


	    ! Combine two fields onto regular grid
	    Call Interp_add_cmplx_fields(n_fd1_points,fd1_node_coords,
     &			    n_fd2_points,fd2_node_coords,
     &                      n_mesh_elts,
     &                      mesh_elt_nodes,xdim,ydim,x0,y0,dx,dy,
     &                      irr_fd1_data,irr_fd2_data,regdata)
	    End Subroutine Combining_cmplx_fields

C  ***********************************************************************

	    Subroutine Define_overlap_rect_grid(n_fd1_pts,fd1_node_crds,
     &				 	        n_fd2_pts,fd2_node_crds,
     &	    					ndx,ndy,xstrt,ystrt,dx,dy)
	    Integer :: n_fd1_pts,n_fd2_pts,ndx,ndy,n
	    Real*8, dimension(:) :: fd1_node_crds,fd2_node_crds
	    Real*8 :: dx,dy,xmin(2),xmax(2),ymin(2),ymax(2),
     &                xstrt,ystrt,xend,yend
 
	    xmin(1:2)=1d14; xmax(1:2)=-1d14  
            ymin(1:2)=1d14; ymax(1:2)=-1d14
	    !
	    Do n=1,n_fd1_pts
              If (fd1_node_crds(2*n-1)<xmin(1)) 
     &          xmin(1)=fd1_node_crds(2*n-1)
              If (fd1_node_crds(2*n-1)>xmax(1)) x
     &          max(1)=fd1_node_crds(2*n-1)
              If (fd2_node_crds(2*n)<ymin(1)) 
     &	        ymin(1)=fd1_node_crds(2*n)
              If (fd2_node_crds(2*n)>ymax(1)) 
     &          ymax(1)=fd1_node_crds(2*n)
            End Do
	    !
	    Do n=1,n_fd2_pts
              If (fd2_node_crds(2*n-1)<xmin(2)) 
     &          xmin(2)=fd2_node_crds(2*n-1)
              If (fd2_node_crds(2*n-1)>xmax(2)) x
     &          max(2)=fd2_node_crds(2*n-1)
              If (fd2_node_crds(2*n)<ymin(2)) 
     &          ymin(2)=fd2_node_crds(2*n)
              If (fd2_node_crds(2*n)>ymax(2)) 
     &          ymax(2)=fd2_node_crds(2*n)
            End Do 
	   
	    xstrt=maxval(xmin); ystrt=maxval(ymin) 
	    xend=minval(xmax); yend=minval(ymax)
	    dx=(xend-xstrt)/ndx; dy=(yend-ystrt)/ndy
	    End Subroutine Define_overlap_rect_grid

C  ***********************************************************************

	    Subroutine Interp_over_trimesh(n_nodes,mesh_coords,n_mesh_elts,
     &                          mesh_elt_nodes,xdim,ydim,x0,y0,dx,dy,
     &                          idata,rdata)
	    Integer :: n_nodes,n_mesh_elts,xdim,ydim
	    Real*8 :: mesh_coords(2*n_nodes)
	    Integer :: mesh_elt_nodes(3*n_mesh_elts),i,j,itri,tri_nodes(3)
	    Real*8 :: x0,y0,dx,dy,x1,y1,xx(3),yy(3),v(3),A(3),AA
	    Real*8 :: idata(n_nodes), rdata(ydim,xdim)

	    Do i=1,ydim
	      y1=y0+(i-1)*dy
	      Do j=1,xdim
	        x1=x0+(j-1)*dx
	        itri=is_in_mesh(x1,y1,n_nodes,mesh_coords,
     &                          n_mesh_elts,mesh_elt_nodes)
	        if (itri>0) then
	          ! (x1,y1) is in tri-mesh, itri is the triangle id
	          tri_nodes(1:3)=mesh_elt_nodes(3*itri-2:3*itri)
	          xx(1)=mesh_coords(2*tri_nodes(1)-1)
	          yy(1)=mesh_coords(2*tri_nodes(1))
	          xx(2)=mesh_coords(2*tri_nodes(2)-1)
	          yy(2)=mesh_coords(2*tri_nodes(2))
	          xx(3)=mesh_coords(2*tri_nodes(3)-1)
	          yy(3)=mesh_coords(2*tri_nodes(3))
	          v(1:3)=idata(tri_nodes(1:3))
	          AA=det(xx(1),xx(2),xx(3),yy(1),yy(2),yy(3))
	          A(1)=det(x1,xx(2),xx(3),y1,yy(2),yy(3))
	          A(2)=det(x1,xx(3),xx(1),y1,yy(3),yy(1))
	          A(3)=det(x1,xx(1),xx(2),y1,yy(1),yy(2))
	          if (AA<1d-14) then
	             print*,'** triangle area close to zero in tri-mesh ...'
	             print*,'itri =',itri
	             print*,'AA =',AA
	             print*,'xx(1:3) =',xx(1:3)
	             print*,'yy(1:3) =',yy(1:3)
	             print*,'tri_nodes(1:3) =', tri_nodes(1:3)
	             stop
	          end if
	          rdata(i,j)=(v(1)*A(1)+v(2)*A(2)+v(3)*A(3))/AA
	        else
	          rdata(i,j)=0d0
	        end if
	        !print*,'rdata(i,j)=',rdata(i,j)
	      End Do
	    End Do
	    End Subroutine Interp_over_trimesh

C  ***********************************************************************

	    Subroutine Interp_add_cmplx_fields(n_fd1_nodes,fd1_node_crds,
     &                        n_fd2_nodes,fd2_node_crds,
     &                        n_mesh_elts,
     &                        mesh_elt_nodes,xdim,ydim,x0,y0,dx,dy,
     &                        fd1_idata,fd2_idata,rdata)
	    Integer :: n_fd1_nodes,n_fd2_nodes,n_mesh_elts(2),
     &		       xdim,ydim
            Real*8, dimension(:) :: fd1_node_crds,fd2_node_crds
	    Integer, dimension(:,:) :: mesh_elt_nodes
	    Real*8 :: x0,y0,dx,dy
	    Complex*16, dimension(:) :: fd1_idata,fd2_idata
	    Complex*16, dimension(:,:) :: rdata

	    Integer :: i,j,itri(2),tri_nodes(3)
	    Real*8 :: x1,y1,xx(3),yy(3),A(3),AA 
            Complex*16 :: v(3)

	    Do i=1,ydim
              y1=y0+(i-1)*dy
              Do j=1,xdim
                x1=x0+(j-1)*dx
                itri(1)=is_in_mesh(x1,y1,n_fd1_nodes,fd1_node_crds,
     &                          n_mesh_elts(1),mesh_elt_nodes(:,1))
                itri(2)=is_in_mesh(x1,y1,n_fd2_nodes,fd2_node_crds,
     &                          n_mesh_elts(2),mesh_elt_nodes(:,2))
	        If (itri(1)>0 .and. itri(2)>0) Then
	          ! (x1,y1) is in both tri-meshes, itri(1:2) are 
                  ! the containing triangle ids

		  tri_nodes(1:3)=mesh_elt_nodes(3*itri(1)-2:3*itri(1),1)
                  xx(1)=fd1_node_crds(2*tri_nodes(1)-1)
                  yy(1)=fd1_node_crds(2*tri_nodes(1))
                  xx(2)=fd1_node_crds(2*tri_nodes(2)-1)
                  yy(2)=fd1_node_crds(2*tri_nodes(2))
                  xx(3)=fd1_node_crds(2*tri_nodes(3)-1)
                  yy(3)=fd1_node_crds(2*tri_nodes(3))
	          v(1:3)=fd1_idata(tri_nodes(1:3))
                  AA=det(xx(1),xx(2),xx(3),yy(1),yy(2),yy(3))
                  A(1)=det(x1,xx(2),xx(3),y1,yy(2),yy(3))
                  A(2)=det(x1,xx(3),xx(1),y1,yy(3),yy(1))
                  A(3)=det(x1,xx(1),xx(2),y1,yy(1),yy(2))
                  if (AA<1d-14) then
	            print*,'** triangle area close to zero in tri-mesh ...'
	            stop
	          end if
	          rdata(i,j)=(v(1)*DCMPLX(A(1))+v(2)*DCMPLX(A(2))+
     &                        v(3)*DCMPLX(A(3)))/DCMPLX(AA)
		  !
	          tri_nodes(1:3)=mesh_elt_nodes(3*itri(2)-2:3*itri(2),2)
                  xx(1)=fd1_node_crds(2*tri_nodes(1)-1)
                  yy(1)=fd1_node_crds(2*tri_nodes(1))
                  xx(2)=fd1_node_crds(2*tri_nodes(2)-1)
                  yy(2)=fd1_node_crds(2*tri_nodes(2))
                  xx(3)=fd1_node_crds(2*tri_nodes(3)-1)
                  yy(3)=fd1_node_crds(2*tri_nodes(3))
                  v(1:3)=fd2_idata(tri_nodes(1:3))
                  AA=det(xx(1),xx(2),xx(3),yy(1),yy(2),yy(3))
                  A(1)=det(x1,xx(2),xx(3),y1,yy(2),yy(3))
                  A(2)=det(x1,xx(3),xx(1),y1,yy(3),yy(1))
                  A(3)=det(x1,xx(1),xx(2),y1,yy(1),yy(2))
                  if (AA<1d-14) then
                    print*,'** triangle area close to zero in tri-mesh ...'
                    stop
                  end if
                  rdata(i,j)=rdata(i,j)+(v(1)*DCMPLX(A(1))+
     &                       v(2)*DCMPLX(A(2))+
     &                       v(3)*DCMPLX(A(3)))/DCMPLX(AA)
	        Else
	          rdata(i,j)=DCMPLX(0d0)
	        End If
	      End Do
	    End Do
	    End Subroutine Interp_add_cmplx_fields

C  ***********************************************************************

	    Function det(x1,x2,x3,y1,y2,y3)
	    Real*8 :: x1,x2,x3,y1,y2,y3,det
	    det = abs((x2*y3-x3*y2)-(x1*y3-x3*y1)+(x1*y2-x2*y1)) 
	    End Function det

	    Function is_in_mesh(x,y,n_nodes,mesh_coords,
     &                          n_mesh_elts,mesh_elt_nodes)
	    Real*8 :: x,y
	    Integer :: n_nodes, n_mesh_elts
            Real*8, dimension(:) :: mesh_coords
            Integer, dimension(:) :: mesh_elt_nodes
	    Integer :: ie,is_in_mesh
 	    Real*8 :: xt(3),yt(3),xin,yin,xm,ym
	  
	    is_in_mesh = -1
            Do ie=1,n_mesh_elts
              xt(1)=mesh_coords(2*mesh_elt_nodes(3*ie-2)-1)
              yt(1)=mesh_coords(2*mesh_elt_nodes(3*ie-2))
              xt(2)=mesh_coords(2*mesh_elt_nodes(3*ie-1)-1)
              yt(2)=mesh_coords(2*mesh_elt_nodes(3*ie-1))
              xt(3)=mesh_coords(2*mesh_elt_nodes(3*ie)-1)
              yt(3)=mesh_coords(2*mesh_elt_nodes(3*ie))

	      ! find a point interior to triangle
              xm=(xt(2)+xt(3))/2d0  
              ym=(yt(2)+yt(3))/2d0
              xin=(xt(1)+xm)/2d0  ! (xin,yin) lies in triangle
              yin=(yt(1)+ym)/2d0

              if (is_in_triangle(x,y,xin,yin,xt,yt)) then
                is_in_mesh = ie
                return
              end if
            End Do 
	    End Function is_in_mesh

C  ***********************************************************************

            Function is_in_triangle(xr,yr,xi,yi,xarr,yarr)
	    Real*8 :: xr,yr,xi,yi,xarr(3),yarr(3)
	    Real*8 :: x(2),y(2),a,b,c,t,s
	    Integer :: e
	    Logical :: is_in_triangle

	    Do e=1,3
	      if (e==1 .or. e==2) then
	        x(1)=xarr(e); x(2)=xarr(e+1)
	        y(1)=yarr(e); y(2)=yarr(e+1)
	      else 
	        x(1)=xarr(e); x(2)=xarr(1)
                y(1)=yarr(e); y(2)=yarr(1)
	      end if
	      a=x(2)-x(1); b=-(y(2)-y(1))
	      c=x(1)*(y(2)-y(1))-y(1)*(x(2)-x(1))
              t=a*yr+b*xr+c; s=a*yi+b*xi+c
	      if (t*s .lt. 0d0) then
	        is_in_triangle=.false.
	        return 
	      end if
	    End Do
	    is_in_triangle=.true.
	    End Function is_in_triangle

	    End module resamp_mod
