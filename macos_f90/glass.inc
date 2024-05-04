C***********************************************************************
C	Begin file glass.inc
C    +-------------------------------------------------------------+
C    |  Copyright (C) 1995-9, California Institute of Technology.  |  
C    |  U.S. Government Sponsorship under NASA Contract NAS7-918   | 
C    |  is acknowledged.					   |
C    +-------------------------------------------------------------+
C***********************************************************************


c	DATA GlassName(1)/'Air'/
c	DATA (GlassTable(i,1),i=1,6)/0d0,0d0,0d0,0d0,0d0,0d0/

c	DATA GlassName(2)/'BK7'/
c	DATA (GlassTable(i,2),i=1,6)/1.03961212,2.31792344e-1,
c     &	1.01046945,6.00069867e-3,2.00179144e-2,1.0350653e+2/



	GlassName(1)='Air'
	CALL DZERO(GlassTable(1,1),6)

	GlassName(2)='BK7'
	GlassTable(1,2)=1.03961212
	GlassTable(2,2)=2.31792344e-1
	GlassTable(3,2)=1.01046945
	GlassTable(4,2)=6.00069867e-3
	GlassTable(5,2)=2.00179144e-2
	GlassTable(6,2)=1.0350653e+2

	GlassName(3)='LAK9'
	GlassTable(1,3)=1.46231905d0
	GlassTable(2,3)=3.44399589d-1
	GlassTable(3,3)=1.15508372d0
	GlassTable(4,3)=7.24270156d-3
	GlassTable(5,3)=2.43353131d-2
	GlassTable(6,3)=8.54686868d1

C	End Glass Include File - 'glass.inc'
