! Dummy Routines for splicubi2.c
! ==> not yet found a way to incorporate the file with f2Py
!
! They are used in mathsub.F

! Call SPLIE2F2C(ya(:,:),xa(:,:),GridMat,nGridMat,nGridMat,srfIdx)
  SUBROUTINE SPLIE2F2C(a,b,mat,n,m,idx)
    REAL*8,  dimension(n), intent(inout):: a
    REAL*8,  dimension(n), intent(inout):: b
    REAL*8,  dimension(n,m), intent(inout):: mat
    INTEGER,                 intent(inout):: n,m,idx

  END SUBROUTINE SPLIE2F2C

!  Call SPLIN2F2C(ya(:,srfIdx),xa(:,srfIdx),nGridMat,nGridMat,y,x,
! &               fh,dfdy,dfdx,j0,i0,srfIdx)
  SUBROUTINE SPLIN2F2C(a,b,n,m,x,y,fh,dfdx,dfdy,i,j,idx)
    REAL*8,  dimension(n), intent(inout):: a
    REAL*8,  dimension(n), intent(inout):: b
    INTEGER,                 intent(inout):: n,m
    REAL*8,                  intent(inout):: x,y,fh,dfdx,dfdy
    INTEGER,                 intent(inout):: i,j,idx

  END SUBROUTINE SPLIN2F2C


