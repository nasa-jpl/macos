C***********************************************************************
C       Begin file pgplotdummy.F
C    +----------------------------------------------------------------+
C    |  Copyright (C) 1995-2007, California Institute of Technology.  |
C    |  U.S. Government Sponsorship Is Acknowledged.                  |
C    +----------------------------------------------------------------+
C***********************************************************************

C This file contains dummy versions of corresponding PGPLOT routines.

C***********************************************************************

#include "realtype.h"

        SUBROUTINE GRAEND
        END SUBROUTINE GRAEND

C***********************************************************************

        SUBROUTINE GRAINI
        END SUBROUTINE GRAINI

C***********************************************************************

        SUBROUTINE CONTOUR(X,Y,Z,NROW,NCOL,NDROW,NDCOL,NCONTOURS,
     &  CXLABL,CZLABL,CTITLE)

        CHARACTER*(*) CXLABL,CZLABL,CTITLE
        INTEGER NROW,NCOL,NDROW,NDCOL,I,J,NCONTOURS
        SREAL X(NROW),Y(NCOL),Z(NDROW,NDCOL),TRCONT(6),
     &  DX,DY,ZMIN,ZMAX,CONTOURS(16)

        END SUBROUTINE CONTOUR

C***********************************************************************

        SUBROUTINE SPOTDIAG(DXYMAT,NROW,NCOL,NDROW,NDCOL,NBEAMS,
     &  		    CXLABL,CZLABL,CTITLE)

        IMPLICIT NONE
        CHARACTER*(*) CXLABL,CZLABL,CTITLE
        INTEGER NROW,NCOL,NDROW,NDCOL,NBEAMS
        SREAL XYMAT(NDROW,NDCOL)
        REAL*8 DXYMAT(NDROW,NDCOL)

        END SUBROUTINE SPOTDIAG

C***********************************************************************

        SUBROUTINE LINSPOTDIAG(XYMAT,NROW,NCOL,NDROW,NDCOL,NBEAMS,
     &  CXLABL,CZLABL,CTITLE)

        CHARACTER*(*) CXLABL,CZLABL,CTITLE
        INTEGER NROW,NCOL,NDROW,NDCOL,I,J,K,NBEAMS,ISYM(4),CI
        SREAL XYMAT(NDROW,NDCOL),XMIN,XMAX,YMIN,YMAX,PAD,DX,DY
        DATA ISYM/20,5,8,16/

        END SUBROUTINE LINSPOTDIAG

C***********************************************************************

        SUBROUTINE SLICE(X,Y,Z,NROW,NCOL,NDROW,NDCOL,PX,PZ,
     &  CXLABL,CYLABL,CTITLE)

        CHARACTER*(*) CXLABL,CYLABL,CTITLE
        INTEGER NROW,NCOL,NDROW,NDCOL,I,J
        SREAL X(NROW),Y(NCOL),Z(NDROW,NDCOL),PX(NROW),PZ(NROW),
     &  DX,DZ,ZMIN,ZMAX,XSKEW,ZSKEW,XBMIN,XBMAX,ZBMIN,ZBMAX

        END SUBROUTINE SLICE

C***********************************************************************

        SUBROUTINE PGSLICE(X,Y,Z,PX,PZ,NROW,NCOL,NDROW,NDCOL,DX,DZ,ZMIN)

        INTEGER NROW,NCOL,NDROW,NDCOL,I,J,IP1,CI
        SREAL X(NROW),Y(NROW),Z(NDROW,NDCOL),PX(NROW),PZ(NROW),DX,DZ,
     &  XSTEP,ZSTEP,ZMIN,POLYX(4),POLYZ(4)

        END SUBROUTINE PGSLICE

C***********************************************************************

        SUBROUTINE GRAY(X,Y,Z,NROW,NCOL,NDROW,NDCOL,
     &  CXLABL,CZLABL,CTITLE)

        CHARACTER*(*) CXLABL,CZLABL,CTITLE
        INTEGER NROW,NCOL,NDROW,NDCOL,I,J
        SREAL X(NROW),Y(NCOL),Z(NDROW,NDCOL),TRCONT(6),
     &  DX,DY,ZMIN,ZMAX

        END SUBROUTINE GRAY

C***********************************************************************

        SUBROUTINE FREDDY(ARRAY,KX,NY,SIZE,ANGLE)
        INTEGER KX, NY
        REAL ARRAY(KX,NY), SIZE, ANGLE

        END SUBROUTINE FREDDY

C***********************************************************************

        SUBROUTINE FREDGO(ARRAY,MN)
        INTEGER MN
        REAL ARRAY(MN)

        END SUBROUTINE FREDGO

C***********************************************************************

        SUBROUTINE WIRE(X,Y,Z,ZVEC,NROW,NCOL,NDROW,NDCOL,ND2,
     &  CXLABL,CZLABL,CTITLE)

        CHARACTER*(*) CXLABL,CZLABL,CTITLE
        INTEGER NROW,NCOL,NDROW,NDCOL,ND2,I,J,K
        SREAL X(NROW),Y(NCOL),Z(NDROW,NDCOL),ZVEC(ND2),DX,DZ,SIZE,
     &  ZMIN,ZMAX

        END SUBROUTINE WIRE

C***********************************************************************

        SUBROUTINE PLOTCOL(X,Y,Z,NROW,NCOL,NDROW,NDCOL,
     &  CXLABL,CYLABL,CTITLE)

        CHARACTER*(*) CXLABL,CYLABL,CTITLE
        INTEGER NROW,NCOL,NDROW,NDCOL,I,J,IMAX,JMAX,CI
        SREAL X(NROW),Y(NCOL),Z(NDROW),ZMIN,ZMAX

        END SUBROUTINE PLOTCOL

C***********************************************************************

        SUBROUTINE DRAW(CTITLE,nDrawElt,nDrawRay,DrawEltVec,
     &                  DrawRayVec_in,XLabel,YLabel)
           use param_mod
           use elt_mod

           IMPLICIT NONE

           INTEGER :: nDrawRay,nDrawElt(mDrawRay),
     &                DrawEltVec(mDrawElt,mDrawRay)
           SREAL :: DrawRayVec_in(2,mDrawElt,mDrawRay)
           CHARACTER*(*) CTITLE
           CHARACTER*8 XLabel,YLabel

           END SUBROUTINE DRAW

C***********************************************************************a


	SUBROUTINE GERROR(STRING)
	CHARACTER*80 STRING
	END SUBROUTINE GERROR
