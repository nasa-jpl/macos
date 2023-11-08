C***********************************************************************
C	Begin file fits1.inc
C    +-------------------------------------------------------------+
C    |  Copyright (C) 1995-9, California Institute of Technology.  |  
C    |  U.S. Government Sponsorship under NASA Contract NAS7-918   | 
C    |  is acknowledged.                                           |
C    +-------------------------------------------------------------+
C***********************************************************************
      CHARACTER BUFF*80,MTERR*64
      COMMON /FIG_COMB/ BUFF,MTERR
      integer*1 FBUFF(2880)
      INTEGER*2 FBUFFS(1440)
      INTEGER FBUFFI(720)
      INTEGER FPTR,MTUNIT,MTSTAT,BLKNO
      COMMON /FIT_COMF/ FPTR,MTUNIT,MTSTAT,BLKNO,FBUFF
      EQUIVALENCE (FBUFFS(1),FBUFF(1)),(FBUFFI(1),FBUFF(1))
C***********************************************************************
 

