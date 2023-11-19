C***********************************************************************
C       Begin file nls.F
C    +----------------------------------------------------------------+
C    |  Copyright (C) 1995-2007, California Institute of Technology.  |
C    |  U.S. Government Sponsorship Is Acknowledged.                  |
C    +----------------------------------------------------------------+
C***********************************************************************

C
C This implemenation is based on the Levenberg-Marquardt
C nonlinear least-square model-fitting method.
C John Z. Lou, Jet Propulsion Laboratory
C

#include "realtype.h"


	Module nls
	use nn_util, only : assert_eq, diagmult
	use nn_mod, only : covsrt, gaussj

	Implicit none

        LREAL, dimension(:,:), allocatable :: dyda
        LREAL, dimension(:), allocatable :: dy,sig2i,wt,ymod


	Contains

 	Subroutine lmlsq(y,sig,a,maska,covar,
     &                   alpha,chisq,funcs_app,alamda,success)
	LREAL, dimension(:), intent(in) :: y,sig
	LREAL, dimension(:), intent(inout) :: a
	LREAL, dimension(:,:), intent(out) :: covar,alpha
	LREAL, intent(out) :: chisq
	LREAL, intent(inout) :: alamda
	Logical, dimension(:), intent(in) :: maska
        Logical :: success

 	Interface
	  Subroutine funcs_app(a_p, yfit_p, dyda_p)
          use traceutil_mod, only : LUseChfRayIfOK
	  LREAL, dimension(:), intent(in) :: a_p
	  LREAL, dimension(:), intent(out) :: yfit_p
	  LREAL, dimension(:,:), intent(out) :: dyda_p
	  End Subroutine funcs_app
	End Interface

	Integer :: ma, ndata
	Integer, save :: mfit


	Call mrqmin_private


	Contains

	 ! An inner subroutine
	 Subroutine mrqmin_private
	 LREAL, save :: ochisq
	 LREAL, dimension(:), allocatable, save :: atry, beta
	 LREAL, dimension(:,:), allocatable, save :: da

	 !integer, dimension(5) :: tmp
	 !tmp(1:5) = (/1,2,2,3,4/)

	 ndata = assert_eq(size(y),size(y),size(sig),'lmlsq: ndata')
	 ma=assert_eq((/size(a),size(maska),size(covar,1),size(covar,2),
     &                size(alpha,1),size(alpha,2)/),'lmlsq: ma')

	 mfit = count(maska)

	 if (alamda < 0.0) then  ! initialization
	   allocate(atry(ma),beta(ma),da(ma,1));
	   alamda = 0.001

           print*,'**mrqmin_private: size(a) =',size(a)
           print*,'**mrqmin_private: size(alpha,1) =',size(alpha,1)
           print*,'**mrqmin_private: size(alpha,2) =',size(alpha,2)
           print*,'**mrqmin_private: size(beta) =',size(beta)

           if (allocated(dyda)) then
              deallocate(dyda,dy,sig2i,wt,ymod)
           end if
           allocate(dyda(size(y),size(a)),dy(size(y)),
     &       sig2i(size(y)),wt(size(y)),ymod(size(y)))

	   Call mrqcof(a, alpha, beta)
	   ochisq = chisq
	   atry = a
	 end if

	 covar(1:mfit,1:mfit) = alpha(1:mfit,1:mfit)
	 call diagmult(covar(1:mfit,1:mfit),1.0+alamda)

	 da(1:mfit,1) = beta(1:mfit)
	 Call gaussj(covar(1:mfit,1:mfit),da(1:mfit,1:1),success)
	 if (.not. success) return
	 if (alamda == 0.0) then
	   Call covsrt(covar,maska)
	   Call covsrt(alpha,maska)
	   deallocate(atry, beta, da)
	   return
	 end if

	 atry = a + unpack(da(1:mfit,1),maska,0d0)
	 Call mrqcof(atry, covar, da(1:mfit,1))

	 if (chisq < ochisq) then
	   alamda = 0.1*alamda
	   ochisq = chisq
	   alpha(1:mfit,1:mfit) = covar(1:mfit,1:mfit)
	   beta(1:mfit) = da(1:mfit,1)
	   a = atry
	 else
	   alamda = 10.0*alamda
	   chisq = ochisq
	 end if
	 End Subroutine mrqmin_private


#if 0
         Subroutine mrqcof(a, alpha, beta)
         LREAL, dimension(:), intent(in) :: a
         LREAL, dimension(:), intent(out) :: beta
         LREAL, dimension(:,:), intent(out) :: alpha

         integer :: j, k, l, m
         LREAL, dimension(size(y),size(a)) :: dyda
         LREAL, dimension(size(y)) :: dy, sig2i
         !LREAL, dimension(size(y)) :: wt, ymod
         LREAL, dimension(261121) :: wt, ymod

         print*,'**mrqcof: size(y) =',size(y)
         Stop
         End Subroutine mrqcof
#endif

#if 1
	 Subroutine mrqcof(a, alpha, beta)
	 LREAL, dimension(:), intent(in) :: a
	 LREAL, dimension(:), intent(out) :: beta
	 LREAL, dimension(:,:), intent(out) :: alpha

	 integer :: j, k, l, m
	 !LREAL, dimension(size(y),size(a)) :: dyda
	 !LREAL, dimension(size(y)) :: dy, sig2i, wt, ymod

         ! Variables below are now defined at module level -05/2013
	 !LREAL, dimension(261121,3) :: dyda
         !LREAL, dimension(261121) :: dy, sig2i
         !LREAL, dimension(261121) :: wt, ymod

         !print*,'**mrqcof: size(y),size(a) =',size(y),size(a)
	 Call funcs_app(a, ymod, dyda)

	 sig2i = 1.0/(sig**2)
	 dy = y - ymod
	 j = 0

	 do l=1, ma
	   if (maska(l)) then
	     j = j+1
	     wt = dyda(:,l)*sig2i
	     k = 0
	     do m=1, l
	       if (maska(m)) then
	         k = k+1
	         alpha(j,k) = dot_product_jl(wt,dyda(:,m))
	         alpha(k,j) = alpha(j,k)
	       end if
	     end do
	     beta(j) = dot_product_jl(dy,wt)
	   end if
	 end do
	 chisq = dot_product_jl(dy**2,sig2i)
	 End Subroutine mrqcof
#endif


	 Function dot_product_jl(xx,yy)
	 LREAL, dimension(:) :: xx, yy
	 LREAL :: dot_product_jl
	 Integer :: i

	 dot_product_jl = 0d0
	 do i=1, size(xx)
	   dot_product_jl = dot_product_jl + xx(i)*yy(i)
	 end do
	 End Function dot_product_jl

  	End Subroutine lmlsq
	End module nls
