C***********************************************************************
C       Begin file lsq.F
C    +----------------------------------------------------------------+
C    |  Copyright (C) 1995-2007, California Institute of Technology.  |
C    |  U.S. Government Sponsorship Is Acknowledged.                  |
C    +----------------------------------------------------------------+
C***********************************************************************

C  This file provides subroutines for computing SVD and pseudo inverse,
C  useful for linear least-squre based system optimizations.
C  Subroutines:
C  pinv - calculate pseudo-inverse of a matrix
C  svd_decomp - calculate SVD decomposition of a matrix
C  John Z. Lou, Jet Propulsion Laboratory


#include "realtype.h"

	Module lsq_mod
	use nn_type
        use nn_util, only : assert_eq,nerror,outerprod

	Implicit none

	Contains

C
C Given an M x N matrix a, this routine computes its SVD,
C A = U*WVt. The matrix U replaces matrix 'a' on output. The 
C diagonal matrix of singular values W is output as the N-dim 
C vector w. The N x N matrix V is output as v.
C
	Subroutine svd_decomp(a,w,v)
        use nn_mod, only : pythag

	LREAL, dimension(:,:), intent(inout) :: a
	LREAL, dimension(:), intent(out) :: w
	LREAL, dimension(:,:) :: v
	Integer :: i, its, j, k, l, m, n, nm
	LREAL :: anorm, c, f, g, h, s, scale, x, y, z
	LREAL, dimension(size(a,1)) :: tempm
	LREAL, dimension(size(a,2)) :: rv1, tempn

        m=size(a,1)
	n=assert_eq(size(a,2),size(v,1),size(v,2),size(w),
     &              'svd_decomp')
	g=0.0; scale=0.0
	Do i=1,n
	  l=i+1; rv1(i)=scale*g
          g=0.0; scale=0.0
	  if (i <= m) then
	    scale = sum(abs(a(i:m,i)))
	    if (scale /= 0.0) then
	      a(i:m,i) = a(i:m,i)/scale
	      s=dot_product(a(i:m,i),a(i:m,i))
	      f=a(i,i)
	      g=-sign(sqrt(s),f)
	      h=f*g-s; a(i,i)=f-g
	      tempn(l:n)=matmul(a(i:m,i),a(i:m,l:n))/h
	      a(i:m,l:n)=a(i:m,l:n)
     &                   +outerprod(a(i:m,i),tempn(l:n))
	      a(i:m,i)=scale*a(i:m,i)
	    end if
	  end if
	  w(i)=scale*g; g=0.0; scale=0.0
	  if ((i <= m) .and. (i /= n)) then
	    scale=sum(abs(a(i,l:n)))
	    if (scale /= 0.0) then
	      a(i,l:n) = a(i,l:n)/scale
	      s=dot_product(a(i,l:n),a(i,l:n))
	      f=a(i,l); g=-sign(sqrt(s),f); h=f*g-s;
	      a(i,l)=f-g; rv1(l:n)=a(i,l:n)/h
	      tempm(l:m) = matmul(a(l:m,l:n),a(i,l:n))
	      a(l:m,l:n) = a(l:m,l:n)
     &                     +outerprod(tempm(l:m),rv1(l:n)) 
	      a(i,l:n)=scale*a(i,l:n)
	    end if
	  end if
	End Do
	anorm = maxval(abs(w)+abs(rv1))

	Do i=n,1,-1
	  if (i < n) then
	    if (g /= 0.0) then
	      v(l:n,i) = (a(i,l:n)/a(i,l))/g
	      tempn(l:n) = matmul(a(i,l:n),v(l:n,l:n))
	      v(l:n,l:n) = v(l:n,l:n)
     &			   +outerprod(v(l:n,i),tempn(l:n))
	    end if
	    v(i,l:n)=0.0; v(l:n,i)=0.0
	  end if
	  v(i,i)=1.0; g=rv1(i); l=i 
	End Do
	
	Do i=min(m,n),1,-1
	  l=i+1; g=w(i); a(i,l:n)=0.0
	  if (g /= 0.0) then
	    g=1.0d0/g
	    tempn(l:n) = (matmul(a(l:m,i),a(l:m,l:n))/a(i,i))*g
	    a(i:m,l:n) = a(i:m,l:n)+outerprod(a(i:m,i),tempn(l:n))
	    a(i:m,i)=a(i:m,i)*g
	  else
	    a(i:m,i)=0.0
	  end if
	  a(i,i)=a(i,i)+1.0d0
	End Do

	Do k=n,1,-1
	  do its=1,30
	    do l=k,1,-1
	      nm=l-1
	      if ((abs(rv1(l))+anorm) == anorm) exit
	      if ((abs(w(nm))+anorm) == anorm) then
	        c=0.0; s=1.0
	        do i=l,k
	          f=s*rv1(i); rv1(i)=c*rv1(i)
	          if ((abs(f)+anorm) == anorm) exit
	          g=w(i); h=pythag(f,g); w(i)=h
	          h=1.0d0/h; c=(g*h); s=-(f*h)
	          tempm(1:m)=a(1:m,nm)
	          a(1:m,nm)=a(1:m,nm)*c+a(1:m,i)*s
	          a(1:m,i)=-tempm(1:m)*s+a(1:m,i)*c 
	        end do ! i
	        exit
	      end if
	    end do ! l
	    z=w(k)

	    if (l==k) then
	      if (z < 0.0) then
	        w(k)=-z; v(1:n,k)=-v(1:n,k)
	      end if
	      exit
	    end if
	    if (its == 30) 
     &        call nerror('svd_decomp: no convergence in svd_decomp')
	    x=w(l); nm=k-1; y=w(nm)
            g=rv1(nm); h=rv1(k)
            f = ((y-z)*(y+z)+(g-h)*(g+h))/(2.0d0*h*y)
	    g = pythag(f,1.0d0)
	    f = ((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
	   
	    c=1.0; s=1.0
	    do j=l,nm
	      i=j+1; g=rv1(i); y=w(i); h=s*g; g=c*g
              z=pythag(f,h)
	      rv1(j)=z; c=f/z; s=h/z  
	      f=(x*c)+(g*s); g=-(x*s)+(g*c)
              h=y*s; y=y*c
	      tempn(1:n) = v(1:n,j)
	      v(1:n,j)=v(1:n,j)*c+v(1:n,i)*s
	      v(1:n,i)=-tempn(1:n)*s+v(1:n,i)*c
	      z=pythag(f,h); w(j)=z
	      if (z /= 0.0) then
	        z=1.0d0/z; c=f*z; s=h*z
	      end if
	      f=(c*g)+(s*y); x=-(s*g)+(c*y)
	      tempm(1:m)=a(1:m,j)
	      a(1:m,j)=a(1:m,j)*c+a(1:m,i)*s
	      a(1:m,i)=-tempm(1:m)*s+a(1:m,i)*c
	    end do ! l
	    rv1(l)=0.0; rv1(k)=f; w(k)=x
	  end do ! its
	End Do ! k

	! Sort singular values in 'w' and re-arrange rows
	! in 'a' and columns in v.
        ! *** Disabled ***, it could mess up the pseudo-inverse
	! based least-square solution !!!  
	Call svd_sort(a,w,v,size(a,1),size(v,1))
	End Subroutine svd_decomp


	! Bubble sort
	Subroutine svd_sort(Am,wv,Vm,m,n)
	Implicit none

        Real*8, dimension(m,m) :: Am
        Real*8, dimension(n,n) :: Vm
        Real*8, dimension(n) :: wv
        Real*8 :: tmpa(m), tmpv(n), wt
        Integer :: m, n, i, j

        Do i=1,n-1
          Do j=i,n
            if (abs(wv(i))<abs(wv(j))) then
              wt=wv(i); wv(i)=wv(j); wv(j)=wt
              tmpa=Am(:,i); Am(:,i)=Am(:,j); Am(:,j)=tmpa
              !tmpv=Vm(i,:); Vm(i,:)=Vm(j,:); Vm(j,:)=tmpv
              tmpv=Vm(:,i); Vm(:,i)=Vm(:,j); Vm(:,j)=tmpv
            end if
          End Do
	End Do	
	End Subroutine svd_sort

	
	!
	! Compute pseudo-inverse of given matrix 'a',
	! output is 'pa' - pseudo-inverse of 'a' 
	!
	Subroutine pinv(a,pa,tol)
	LREAL, dimension(:,:) :: a, pa
	LREAL, intent(in), optional :: tol
	LREAL, dimension(:), allocatable, save :: w
	LREAL, dimension(:,:), allocatable, save :: v
	Integer, save :: m, n, i, pinv_err
	LREAL :: cutoff
	Logical, save :: first_entry= .true.

        m=assert_eq(size(a,1),size(pa,2),'pinv')
        n=assert_eq(size(a,2),size(pa,1),'pinv')

	if (first_entry .or. n /= size(w)) then
	  if (.not. first_entry) then
	    deallocate(w,v,stat=pinv_err)
 	    if (pinv_err /= 0)
     &      call macos_memory_failure
     &             ('pinv: deallocate failed!') 
	  else
	    first_entry = .false.
	  end if
	  allocate(w(n),v(n,n),stat=pinv_err)
 	  if (pinv_err /= 0)
     &      call macos_memory_failure
     &             ('pinv: allocate failed!')
	end if

	! SVD, and singular values are sorted
	Call svd_decomp(a,w,v)

	! Compute pinv 
 	If (present(tol)) Then
 	  cutoff=max(size(a,1),size(a,2))*abs(w(1))*tol
 	Else
 	  cutoff=max(size(a,1),size(a,2))*abs(w(1))*eps
        End If
	!
	pa=0.d0
	Do i=1,n
	  if (abs(w(i)) < cutoff) w(i)=0.d0 
	  if (w(i) /= 0) then
	    Call RankOneMultiply(a(:,i),v(:,i),w(i),pa)
	  end if 
	End Do
	End Subroutine pinv


	Subroutine RankOneMultiply(u,v,w,p)
	LREAL, dimension(:) :: u, v
	LREAL :: w
	LREAL, dimension(:,:) :: p
	Integer :: m, n, i, j

	m=size(u); n=size(v)
	!print*,'size u=', m
	!print*,'size v=', n

	Do i=1,n ! row 
	  Do j=1,m ! col
	    p(i,j)=p(i,j)+v(i)*u(j)/w
	  End Do
	End Do
	End Subroutine RankOneMultiply


#if 0
! 'matmul' intrinsic function is supported on Linux Lahey f95 compiler
!
! Dummy one since can't make Sun f90 compiler to support them
        Function matmul(m1,m2)
        LREAL, dimension(:,:) :: m1, m2
        LREAL, dimension(size(m1,1),size(m2,2)) :: matmul
        ! dummy ...
        End Function matmul
#endif
        End Module lsq_mod


