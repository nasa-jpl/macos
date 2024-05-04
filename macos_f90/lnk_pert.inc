
	SUBROUTINE LnkEltCPERTURB(jElt,th,del,LLocCrd,LUpdateTElt)
	use param_mod
	use elt_mod
        use math_mod
	use traceutil_mod, only : CRIncidPosNom

	Integer :: jElt,j
	Real*8 :: th(3),del(3),Q(3,3),dQ(3,3),D1(9),D2(9)
	Logical :: LLocCrd,LUpdateTElt

        CALL Qform(Q,dQ,th)
        CALL DMPROD(D1,Q,psiElt(1,jElt),3,3,1)
        CALL DEQUATE(psiElt(1,jElt),D1,3)

        CALL DSUB(D1,VptElt(1,jElt),RptElt(1,jElt),3)
        CALL DMPROD(D2,dQ,D1,3,3,1)
        CALL DADD(D1,VptElt(1,jElt),D2,3)
        CALL DADD(VptElt(1,jElt),D1,del,3)

        ! Since RptElt is not affected by the rotation (about itself),
        ! its update is only affected by the translation part of
        ! perturbation.
        CALL DADD(D1,RptElt(1,jElt),del,3)
        CALL DEQUATE(RptElt(1,jElt),D1,3)

        ! Update CRIncidPosNom, added by jzlou
        CALL DSUB(D1,CRIncidPosNom(1,jElt),RptElt(1,jElt),3)
        CALL DMPROD(D2,dQ,D1,3,3,1)
        CALL DADD(D1,CRIncidPosNom(1,jElt),D2,3)
        CALL DADD(CRIncidPosNom(1,jElt),D1,del,3)

	!IF (nObs(jElt).NE.0) THEN
        ! This IF test is removed by jzlou, since an aperture needs
        ! to be perturbed with the element.
            CALL DMPROD(D1,Q,xObs(1,jElt),3,3,1)
            CALL DEQUATE(xObs(1,jElt),D1,3)
        !END IF

        IF ((EltID(jElt).EQ.4).OR.(EltID(jElt).EQ.5).OR.
     &       (EltID(jElt).EQ.13)) THEN
          CALL DSUB(D1,h1HOE(1,jElt),RptElt(1,jElt),3)
          CALL DMPROD(D2,dQ,D1,3,3,1)
          CALL DADD(D1,h1HOE(1,jElt),D2,3)
          CALL DADD(h1HOE(1,jElt),D1,del,3)
          CALL DSUB(D1,h2HOE(1,jElt),RptElt(1,jElt),3)
          CALL DMPROD(D2,dQ,D1,3,3,1)
          CALL DADD(D1,h2HOE(1,jElt),D2,3)
          CALL DADD(h2HOE(1,jElt),D1,del,3)
        END IF

	IF ((SrfType(jElt).GE.4).AND.(SrfType(jElt).LE.9)) THEN
          CALL DSUB(D1,pMon(1,jElt),RptElt(1,jElt),3)
          CALL DMPROD(D2,dQ,D1,3,3,1)
          CALL DADD(D1,pMon(1,jElt),D2,3)
          CALL DADD(pMon(1,jElt),D1,del,3)
          CALL DMPROD(D1,Q,xMon(1,jElt),3,3,1)
          CALL DEQUATE(xMon(1,jElt),D1,3)
          CALL DMPROD(D1,Q,yMon(1,jElt),3,3,1)
          CALL DEQUATE(yMon(1,jElt),D1,3)
          CALL DMPROD(D1,Q,zMon(1,jElt),3,3,1)
          CALL DEQUATE(zMon(1,jElt),D1,3)
        END IF

	If (LLocCrd .and. LUpdateTElt) Then
          ! TElt(1:6,1:6,jElt) is updated here
          Do j=1,3
            CALL DMPROD(D1,Q,TElt(1:3,j,jElt),3,3,1)
            CALL DEQUATE(TElt(1:3,j,jElt),D1,3)
            CALL DMPROD(D1,Q,TElt(4:6,3+j,jElt),3,3,1)
            CALL DEQUATE(TElt(4:6,3+j,jElt),D1,3)
          End Do
        End If
	END SUBROUTINE LnkEltCPERTURB
