C***********************************************************************
C       Begin file nn_util.F
C    +----------------------------------------------------------------+
C    |  Copyright (C) 1995-2007, California Institute of Technology.  |
C    |  U.S. Government Sponsorship Is Acknowledged.                  |
C    +----------------------------------------------------------------+
C***********************************************************************

#include "realtype.h"


	Module nn_type
	implicit none

	LREAL, parameter :: pi=3.1415926535897973
	LREAL, parameter :: eps=2.2204d-18

        End Module nn_type


C ------

	Module nn_util
	use nn_type

	implicit none

	interface swap
	  module procedure swap_i, swap_r, swap_rv
	end interface

	interface assert_eq	
	  module procedure assert_eq2, assert_eq3, assert_eq4, assert_eqn
	end interface

        interface diagmult
	  module procedure diagmult_rv, diagmult_r
	end interface

	interface outerprod
	  module procedure outerprod_d
	end interface

	Contains

	Function assert_eq2(n1,n2,string)
	character(len=*), intent(in) :: string
	integer, intent(in) :: n1, n2
	integer :: assert_eq2

	if (n1==n2) then
	  assert_eq2=n1
        else
	  print*,'assert_eq2(): failed with tag: ', string 
	  stop 'program terminated by assert_eq2'
	end if
	End Function assert_eq2

	Function assert_eq3(n1,n2,n3,string)
        character(len=*), intent(in) :: string
        integer, intent(in) :: n1, n2,n3
	integer :: assert_eq3

        if (n1==n2 .and. n2==n3) then
          assert_eq3=n1 
        else
          print*,'assert_eq3(): failed with tag: ', string 
          stop 'program terminated by assert_eq3'
        end if
        End Function assert_eq3

	Function assert_eq4(n1,n2,n3,n4,string)
        character(len=*), intent(in) :: string
        integer, intent(in) :: n1, n2,n3,n4
	integer :: assert_eq4

        if (n1==n2 .and. n2==n3 .and. n3==n4) then
          assert_eq4=n1 
        else
          print*,'assert_eq4(): failed with tag: ', string 
          stop 'program terminated by assert_eq4'
        end if
        End Function assert_eq4

	Function assert_eqn(nn,string)
        character(len=*), intent(in) :: string
        integer, dimension(:), intent(in) :: nn
	integer :: assert_eqn

        if (all(nn(2:)==nn(1))) then
          assert_eqn=nn(1) 
        else
          print*,'assert_eqn(): failed with tag: ', string 
          stop 'program terminated by assert_eqn'
        end if
        End Function assert_eqn


	Subroutine diagmult_rv(mat, diag)
	LREAL, dimension(:,:), intent(inout) :: mat
	LREAL, dimension(:), intent(in) :: diag
	integer :: j, n

	n = assert_eq2(size(diag),
     &		       min(size(mat,1),size(mat,2)),'diagmult_rv')
	do j=1,n
	  mat(j,j) = mat(j,j)*diag(j)
	end do
	End Subroutine diagmult_rv

	Subroutine diagmult_r(mat, diag)
        LREAL, dimension(:,:), intent(inout) :: mat
        LREAL, intent(in) :: diag
        integer :: j, n

        n = min(size(mat,1),size(mat,2))
        do j=1,n
          mat(j,j) = mat(j,j)*diag
        end do
        End Subroutine diagmult_r


	Subroutine swap_i(a,b)
	integer, intent(inout) :: a, b
	integer :: dum
	dum=a; a=b; b=dum
	End Subroutine swap_i

        Subroutine swap_r(a,b)
        LREAL, intent(inout) :: a, b
        LREAL :: dum
        dum=a; a=b; b=dum
        End Subroutine swap_r

	Subroutine swap_rv(a,b)
        LREAL, dimension(:), intent(inout) :: a, b
        LREAL, dimension(size(a)) :: dum
        dum=a; a=b; b=dum
        End Subroutine swap_rv

	Function outerprod_d(a,b)
        LREAL, dimension(:), intent(in) :: a, b
        LREAL, dimension(size(a),size(b)) :: outerprod_d
        outerprod_d = spread(a,dim=2,ncopies=size(b)) *
     &    spread(b,dim=1,ncopies=size(a))
        End Function outerprod_d

	Function outerand(a,b)
	logical, dimension(:), intent(in) :: a, b
	logical, dimension(size(a),size(b)) :: outerand
 	outerand = spread(a,dim=2,ncopies=size(b)) .and.
     &    spread(b,dim=1,ncopies=size(a))
	End Function outerand

	Subroutine nerror(string)
	character(len=*), intent(in) :: string
	print*,'nerror: ', string
	stop 'program terminated by nerror'
	End Subroutine nerror

	end module nn_util


C -----

	Module nn_mod
	use nn_type
	
	implicit none

	Contains

	Subroutine covsrt(covar, maska)
	use nn_type
	use nn_util, only : assert_eq, swap

	LREAL, dimension(:,:), intent(inout) :: covar
	logical, dimension(:), intent(in) :: maska
  	integer :: ma, mfit, j, k
	
	ma = assert_eq(size(covar,1), size(covar,2),size(maska),'covsrt')
	mfit = count(maska)
	covar(mfit+1:ma,1:ma) = 0.0
	covar(1:ma,mfit+1:ma) = 0.0
	k = mfit	
	do j=ma, 1, -1
	  if (maska(j)) then
	    call swap(covar(1:ma,k), covar(1:ma,j))
	    call swap(covar(k,1:ma),covar(j,1:ma))
	    k = k-1
	  end if
	end do
	End Subroutine covsrt


	Subroutine gaussj(a,b,success)
        use nn_type
	use nn_util, only : assert_eq, nerror, outerand, outerprod, swap

	LREAL, dimension(:,:), intent(inout) :: a, b
	logical :: success
	integer, dimension(size(a,1)) :: ipiv, indxr, indxc
	logical, dimension(size(a,1)) :: lpiv
	LREAL :: pivinv
	LREAL, dimension(size(a,1)) :: dumc
	integer, target :: irc(2)
	integer :: i, l, n
	integer, pointer :: irow, icol

	success=.true.

	if (.false.) then
	! Debugging: check matrix 'a'
	print*,'gaussj: matrix a: '
	do i=1, size(a,1)	
	  print*, a(i,:)
	end do	
	print*,' '
	end if

	n = assert_eq(size(a,1),size(a,2),size(b,1),'gaussj')
	irow => irc(1); icol => irc(2)
	ipiv=0
	
	Do i=1, n
	  lpiv = (ipiv == 0)
	  irc = maxloc(abs(a), outerand(lpiv,lpiv))
	  ipiv(icol) = ipiv(icol)+1 
	  if (ipiv(icol) > 1) then
            print*,'gaussj: singular matrix (1)'
	    success=.false.
	    return
	  end if
	  if (irow /= icol) then
	    call swap(a(irow,:),a(icol,:))
	    call swap(b(irow,:),b(icol,:))
	  end if
	  indxr(i)=irow; indxc(i)=icol
	  if (a(icol,icol)==0.0) then
	    print*,'gaussj: singular matrix (2)'
            success=.false.
            return
          end if
!    &	    call nerror('gaussj: singular matrix (2)')
	  pivinv = 1.0/a(icol,icol)
	  a(icol,icol) = 1.0
	  a(icol,:) = a(icol,:)*pivinv
	  b(icol,:) = b(icol,:)*pivinv
	  dumc = a(:,icol)
	  a(:,icol) = 0.0
	  a(icol,icol) = pivinv
	  a(1:icol-1,:) =
     &      a(1:icol-1,:)-outerprod(dumc(1:icol-1),a(icol,:))
	  b(1:icol-1,:) = 
     &      b(1:icol-1,:)-outerprod(dumc(1:icol-1),b(icol,:))
	  a(icol+1:,:) = a(icol+1:,:)-outerprod(dumc(icol+1:),a(icol,:))
	  b(icol+1:,:) = b(icol+1:,:)-outerprod(dumc(icol+1:),b(icol,:))
	End Do

	Do l=n,1,-1
          call swap(a(:,indxr(l)),a(:,indxc(l)))
	End Do
	End Subroutine gaussj


	Function pythag(a,b)
	implicit none

	LREAL, intent(in) :: a, b
	LREAL :: pythag, absa, absb

	absa=abs(a); absb=abs(b)
	if (absa > absb) then
	  pythag = absa*sqrt(1.0d0+(absb/absa)**2)
	else
	  if (absb==0) then
	    pythag = 0d0
	  else
	    pythag = absb*sqrt(1.0d0+(absa/absb)**2)
	  end if
	end if
	End Function pythag

	End Module nn_mod
