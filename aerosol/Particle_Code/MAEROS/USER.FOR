      FUNCTION BETA(Y,X,TGAS,PGAS,NBTYPE)
C................................................................
C	CONSTANT BETA
C................................................................
      COMMON/PHYSPT/ACELOV,AFLROV,AWALOV,CHI,DELDIF,DELSAT,DENSTY,FSLIP,
     $FTHERM,GAMMA,PSAT,STICK,TGRADC,TGRADF,TGRADW,TKGOP,TURBDS,VOLUME,
     $WTCONM,WTMOL
      BETA=4.0E-16
      V=EXP(X)
      U=EXP(Y)
C...............................................................
C	CONVERT FOR SECTIONALIZATION BY MASS
C...............................................................
      GO TO (2,1,2,3,3,1),NBTYPE
1     BETA=BETA/V
      RETURN
2     BETA=BETA/U
      RETURN
3     BETA=BETA*(1./V+1./U)
      RETURN
      END
C.......................................................   
      FUNCTION GROWTH(X,DUMMY,TGAS,PGAS,NDUM)
C.......................................................
C    VOLUME GROWTH; 8613 HOMEWORK PROBLEM
C.......................................................
      V=EXP(X)
      GROWTH=0.02
      RETURN
      END
	
	
