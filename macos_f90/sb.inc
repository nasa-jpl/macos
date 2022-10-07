C
C  ApType= SunBlades
C  ApVec=  xc yc R1 R2 N  % R1=innerRadius,R2=outerRadius,N=numberOfBlades
C 
	ELSE IF (ApType.EQ.SunBlades) THEN
	  px2=px-ApVec(1); py2=py-ApVec(2)
	  r1=ApVec(3); r2=ApVec(4)
          nBlades=ApVec(5)	  

	  dth=2d0*pi/nBlades
	  th1=0d0; th2=th1+dth

	  If (firstEntry_sunBlades) Then
	    ! Precompute x(3,0:nBlades),y(3,0:nBlades),z(3,0:nBlades)
	    Do ib=1,nBlades
	      xsb(1)=r1*cos(th1); ysb(1)=r1*sin(th1)    
	      xsb(2)=r1*cos(th2); ysb(2)=r1*sin(th2)
	      xsb(3)=r2*cos((th1+th2)/2d0)  
	      ysb(3)=r2*sin((th1+th2)/2d0)
	      Do j=1,3
	        if (j==1) then
	          i1=1; i2=2
	        else if (j==2) then
	          i1=2; i2=3
	        else
		  i1=3; i2=1
	        end if
	        if (abs(xsb(i1)-xsb(i2))<1d-20) then
	          ! horizontal line
	          Asb(1,ib)=1;Bsb(1,ib)=0d0;Csb(1,ib)=-xsb(i1)
	        else
	          Asb(1,ib)=(ysb(i2)-ysb(i1))/(xsb(i2)-xsb(i1))
	          Bsb(1,ib)=-1d0
	          Csb(1,ib)=ysb(i1)-xsb(i1)*Asb(1,ib) 
	        end if
	      End Do
	      ! xcm,ycm store center of mass of triangle
	      xcm(ib)=(xsb(1)+xsb(2)+xsb(3))/3d0 
	      ycm(ib)=(ysb(1)+ysb(2)+ysb(3))/3d0 
	    End Do
	    firstEntry_sunBlades=.false.
	  End If

	  inSB=.false.
	  Do ib=1,nBlades
	    r2=px2**2+py2**2 
	    If ((r2.LE.r1**2) .OR. 
     &          inTrig(Asb(1:3,ib),Bsb(1:3,ib),Csb(1:3,ib),
     &                 xcm(ib),ycm(ib),px2,py2)) Then
	      inSB=.true.
	      exit
	    End If 
	  End Do 

C *************************************************************************

	Logical Function inTrig(A,B,C,xin,yin,xp,yp)
	Real*8 :: A(3),B(3),C(3),xin,yin,xp,yp
	Integer :: i

	Do i=1,3
  	  vin(i)=A(i)*xin+B(i)*xin+C(i)
	  vp(i)=A(i)*xp+B(i)*yp+C(i)
 	End Do
	inTrig=.false.
	If ((vin(1)*vp(1).GE.0) .AND. (vin(2)*vp(2).GE.0)
     &      .AND. (vin(3)*vp(3).GE.0)) Then	
          inTrig=.true.
	End If
	End Function inTrig

