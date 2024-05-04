C
C propsub.inc
C

#if 0
 115    FORMAT(1P,' Ray ',i4,' segment from Element ',i4,' (',a8,
     &  ') to Element ',i4,' (',a8,'):'/
     &  '   Starting point=',3D17.9/'        End point=',3D17.9/
     &  '        Direction=',3D17.9/'   Segment Length=',D17.9/
     &  '     Total Length=',D17.9)
#endif
 115    FORMAT(1P,' Ray ',i4,' segment from Element ',i4,' (',a8,
     &  ') to Element ',i4,' (',a8,'):'/
     &  '   Starting point=',3D22.14/'        End point=',3D22.14/
     &  '        Direction=',3D22.14/'   Segment Length=',D22.14/
     &  '     Total Length=',D22.14)

 116    FORMAT(' Ray',i6,' has aperture coordinates',i5,' and',i5,
     &  '.'/' WF coordinates are ',i4,' and ',i4,'.')
 117    FORMAT(' Ray',i6,' is the chief ray.')
 118    FORMAT(1P,'       Index(i-1)=',2d17.9/
     &  '         Index(i)=',2d17.9)
 119    FORMAT(1P,15X,'Ex=',2d17.9/15X,'Ey=',2d17.9/
     &  15X,'Ez=',2d17.9/13X,'Ixyz=',d17.9)
 122    FORMAT(/' Ray',i6,' becomes undefined after element',i3)
 123    FORMAT(/' Ray',i6,' is obscured by element',i3)
 124    FORMAT(' Interpolated surface',i4,' traced using base ',
     &  'conicoid only.'/' Use SINT command to set interpolation data.')
 502    FORMAT(1P,'  NF/GS Prop between Elt',i3,' and Elt',i3,
     &  ' to WF',i3,':'/
     &  '  z1=',d11.4,' dx1=',d11.4/'  z2=',d11.4,' dx2=',d11.4,
     &  ' min=',d11.4,' max=',d11.4,' dev=',d11.4,0P,' lin=',f6.2,'%')
 503    FORMAT(1P,'  Scalar FF Prop between Elt',i3,' and Elt',i3,
     &  ' to WF',i3,':'/
     &  '  z1=',d11.4,' dx1=',d11.4,' min=',d11.4,' max=',d11.4,
     &  ' dev=',d11.4,0P,' lin=',f6.2,'%'/1P,'  z2=',d11.4,
     &  ' dx2=',d11.4)
 504    FORMAT(1P,'  NF/PP Prop between Elt',i3,' and Elt',i3,
     &  ' to WF',i3,':'/
     &  '  z1=',d11.4,' dx1=',d11.4/'  z2=',d11.4,' dx2=',d11.4,
     &  ' min=',d11.4,' max=',d11.4,' dev=',d11.4,0P,' lin=',f6.2,'%')
 505    FORMAT(1P,'  NF/SS Prop between Elt',i3,' and Elt',i3,
     &  ' to WF',i3,':'/
     &  '  z1=',d11.4,' dx1=',d11.4,' min=',d11.4,' max=',d11.4,
     &  ' dev=',d11.4,0P,' lin=',f6.2,'%'/1P,'  z2=',d11.4,
     &  ' dx2=',d11.4,
     &  ' min=',d11.4,' max=',d11.4,' dev=',d11.4,0P,' lin=',f6.2,'%')
 506    FORMAT('  Geometric Prop between Elt',i3,' and Elt',i3,
     &  ' to WF',i3)
 507    FORMAT(1P,'  SF/SS Prop between Elt',i3,' and Elt',i3,
     &  ' to WF',i3,':'/
     &  'Pin hole radius=',d11.4,' Far-field dx=',d11.4/
     &  '  z1=',d11.4,' dx1=',d11.4,' min=',d11.4,' max=',d11.4,
     &  ' dev=',d11.4,0P,' lin=',f6.2,'%'/1P,'  z2=',d11.4,
     &  ' dx2=',d11.4,
     &  ' min=',d11.4,' max=',d11.4,' dev=',d11.4,0P,' lin=',f6.2,'%')
 508    FORMAT(1P,'  SF/SS Prop between Elt',i3,' and Elt',i3,
     &  ' to WF',i3,':'/
     &  'Pin hole radius=',d11.4,' Far-field dx=',d11.4/
     &  '  z1=',d11.4,' dx1=',d11.4,/'  z2=',d11.4,
     &  ' dx2=',d11.4,
     &  ' min=',d11.4,' max=',d11.4,' dev=',d11.4,0P,' lin=',f6.2,'%')
 509    FORMAT(1P,'  Vector FF Prop between Elt',i3,' and Elt',i3,
     &  ' to WF',i3,':'/
     &  '  z1=',d11.4,' dx1=',d11.4,' min=',d11.4,' max=',d11.4,
     &  ' dev=',d11.4,0P,' lin=',f6.2,'%'/1P,'  z2=',d11.4,
     &  ' dx2=',d11.4)
 510    FORMAT(1P,'  Prop to reference point Elt',i3,' and Elt',i3,
     &  ' to WF',i3,':'/
     &  '  z1=',d11.4,' dx1=',d11.4,' min=',d11.4,' max=',d11.4,
     &  ' dev=',d11.4,0P,' lin=',f6.2,'%'/1P,'  z2=',d11.4,
     &  ' dx2=',d11.4)
 511    FORMAT(1P,'  Prop from reference point Elt',i3,' and Elt',i3,
     &  ' to WF',i3,':'/
     &  '  z1=',d11.4,' dx1=',d11.4/1P,'  z2=',d11.4,' dx2=',d11.4,
     &  ' min=',d11.4,' max=',d11.4,' dev=',d11.4,0P,' lin=',f6.2,'%')
 512    FORMAT(1P,'  Fresnel NF Prop between Elt',i3,' and Elt',i3,
     &  ' to WF',i3,':'/
     &  '  z1=',d11.4,' dx1=',d11.4,' min=',d11.4,' max=',d11.4,
     &  ' dev=',d11.4,0P,' lin=',f6.2,'%'/1P,'  z2=',d11.4,
     &  ' dx2=',d11.4)
 513    FORMAT(1P,'  NF DFT PP Prop between Elt',i3,' and Elt',i3,
     &  ' to WF',i3,':'/
     &  '  z1=',d11.4,' dx1=',d11.4,' min=',d11.4,' max=',d11.4,
     &  ' dev=',d11.4,0P,' lin=',f6.2,'%'/1P,'  z2=',d11.4,
     &  ' dx2=',d11.4)
 514    FORMAT(1P,'  NF DFT SS Prop between Elt',i3,' and Elt',i3,
     &  ' to WF',i3,':'/
     &  '  z1=',d11.4,' dx1=',d11.4,' min=',d11.4,' max=',d11.4,
     &  ' dev=',d11.4,0P,' lin=',f6.2,'%'/1P,'  z2=',d11.4,
     &  ' dx2=',d11.4)
 515    FORMAT(1P,'  FF DFT Prop between Elt',i3,' and Elt',i3,
     &  ' to WF',i3,':'/
     &  '  z1=',d11.4,' dx1=',d11.4,' min=',d11.4,' max=',d11.4,
     &  ' dev=',d11.4,0P,' lin=',f6.2,'%'/1P,'  z2=',d11.4,
     &  ' dx2=',d11.4)
