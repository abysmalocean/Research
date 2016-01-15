*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE LSGETP( LINOBJ, SINGLR, UNITGZ, UNITQ,
     $                   N, NCLIN, NFREE,
     $                   NROWA, NQ, NROWR, NRANK, NUMINF, NZ1,
     $                   ISTATE, KX, CTP, PNORM,
     $                   A, AP, RES, HZ, P,
     $                   GQ, CQ, R, ZY, WORK )

*     IMPLICIT           DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION  ONE,    ZERO
      DOUBLE PRECISION  CTP,    PNORM,  GTP
      DOUBLE PRECISION  DDOT,   DNRM2
      INTEGER           LDBG,   NOUT,  ILSDBG, N,      NCLIN,  NFREE,  
     &                  NROWA,  NQ,    NROWR,  NRANK,  NUMINF, NZ1,    
     &                  J,      I
      LOGICAL            LINOBJ, SINGLR, UNITGZ, UNITQ
      INTEGER            ISTATE(N+NCLIN), KX(N)
      DOUBLE PRECISION   A(NROWA,*), AP(*), RES(*), HZ(*), P(N),
     $                   GQ(N), CQ(*), R(NROWR,*), ZY(NQ,*)
      DOUBLE PRECISION   WORK(N)

************************************************************************
*     LSGETP  computes the following quantities for  LSCORE.
*     (1) The vector  (hz1) = (Rz1)(pz1).
*         If X is not yet feasible,  the product is computed directly.
*         If  Rz1 is singular,  hz1  is zero.  Otherwise  hz1  satisfies
*         the equations
*                        Rz1'hz1 = -gz1,
*         where  g  is the total gradient.  If there is no linear term
*         in the objective,  hz1  is set to  dz1  directly.
*     (2) The search direction P (and its 2-norm).  The vector P is
*         defined as  Z*(pz1), where  (pz1)  depends upon whether or
*         not X is feasible and the nonsingularity of  (Rz1).
*         If  NUMINF .GT. 0,  (pz1)  is the steepest-descent direction.
*         Otherwise,  x  is the solution of the  NZ1*NZ1  triangular
*         system   (Rz1)*(pz1) = (hz1).
*     (3) The vector Ap,  where A is the matrix of linear constraints.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original version written 31-October-1984.
*     Level 2 Blas added 11-June-1986.
*     This version of LSGETP dated 11-June-1986.
************************************************************************
      COMMON    /SOL1CM/ NOUT

      LOGICAL            LSDBG
      PARAMETER         (LDBG = 5)
      COMMON    /LSDEBG/ ILSDBG(LDBG), LSDBG

      EXTERNAL           DDOT  , DNRM2
      INTRINSIC          MIN
      PARAMETER        ( ZERO = 0.0D+0, ONE  = 1.0D+0 )

      IF (SINGLR) THEN
*        ---------------------------------------------------------------
*        The triangular factor for the current objective function is
*        singular,  i.e., the objective is linear along the last column
*        of Z1.  This can only occur when UNITGZ is TRUE.
*        ---------------------------------------------------------------
         IF (NZ1 .GT. 1) THEN
            CALL DCOPY ( NZ1-1, R(1,NZ1), 1, P, 1 )
            CALL DTRSV ( 'U', 'N', 'N', NZ1-1, R, NROWR, P, 1 )
         END IF
         P(NZ1) = - ONE

         GTP = DDOT  ( NZ1, GQ, 1, P, 1 )
         IF (GTP .GT. ZERO) CALL DSCAL ( NZ1, (-ONE), P, 1 )

         IF (NZ1 .LE. NRANK) THEN
            IF (NUMINF .EQ. 0) THEN
               IF (UNITGZ) THEN
                  HZ(NZ1) = R(NZ1,NZ1)*P(NZ1)
               ELSE
                  CALL DLOAD ( NZ1, (ZERO), HZ, 1 )
               END IF
            ELSE
               HZ(1)   = R(1,1)*P(1)
            END IF
         END IF
      ELSE
*        ---------------------------------------------------------------
*        The objective is quadratic in the space spanned by Z1.
*        ---------------------------------------------------------------
         IF (LINOBJ) THEN
            IF (UNITGZ) THEN
               IF (NZ1 .GT. 1)
     $            CALL DLOAD ( NZ1-1, (ZERO), HZ, 1 )
               HZ(NZ1) = - GQ(NZ1)/R(NZ1,NZ1)
            ELSE
               CALL DCOPY ( NZ1, GQ  , 1, HZ, 1 )
               CALL DSCAL ( NZ1, (-ONE), HZ, 1 )
               CALL DTRSV ( 'U', 'T', 'N', NZ1, R, NROWR, HZ, 1 )
            END IF
         ELSE
            CALL DCOPY ( NZ1, RES, 1, HZ, 1 )
         END IF

*        Solve  Rz1*pz1 = hz1.

         CALL DCOPY ( NZ1, HZ, 1, P, 1 )
         CALL DTRSV ( 'U', 'N', 'N', NZ1, R, NROWR, P, 1 )
      END IF

*     Compute  p = Z1*pz1  and its norm.

      IF (LINOBJ)
     $   CTP = DDOT  ( NZ1, CQ, 1, P, 1 )
      PNORM  = DNRM2 ( NZ1, P, 1 )

      CALL CMQMUL( 1, N, NZ1, NFREE, NQ, UNITQ, KX, P, ZY, WORK )

      IF (LSDBG  .AND.  ILSDBG(2) .GT. 0)
     $   WRITE (NOUT, 1000) (P(J), J = 1, N)

*     Compute  Ap.

      IF (NCLIN .GT. 0) THEN
         CALL DLOAD ( NCLIN, ZERO, AP, 1 )
         DO 410 J = 1, N
            IF (ISTATE(J) .LE. 0)
     $         CALL DAXPY( NCLIN, P(J), A(1,J), 1, AP, 1 )
  410    CONTINUE
         IF (LSDBG  .AND.  ILSDBG(2) .GT. 0)
     $   WRITE (NOUT, 1100) (AP(I), I = 1, NCLIN)
      END IF

      RETURN

 1000 FORMAT(/ ' //LSGETP//   P ... ' / (1P5E15.5))
 1100 FORMAT(/ ' //LSGETP//  AP ... ' / (1P5E15.5))

*     End of  LSGETP.

      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE LSGSET( PRBTYP, LINOBJ, SINGLR, UNITGZ, UNITQ,
     $                   N, NCLIN, NFREE,
     $                   NROWA, NQ, NROWR, NRANK, NZ, NZ1,
     $                   ISTATE, KX,
     $                   BIGBND, TOLRNK, NUMINF, SUMINF,
     $                   BL, BU, A, RES, FEATOL,
     $                   GQ, CQ, R, X, WTINF, ZY, WRK )

*     IMPLICIT           DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION  HALF,   ONE,    ZERO
      INTEGER           N,      NCLIN,  NFREE,  NROWA,  NQ,     NROWR,  
     &                  NRANK,  NZ,     NZ1,    J,      NUMINF, K
      DOUBLE PRECISION  BIGBND, TOLRNK, SUMINF, BIGUPP, BIGLOW, FEASJ,  
     &                  CTX,    WEIGHT, S
      DOUBLE PRECISION  DDOT
      INTEGER           IDRANK
      CHARACTER*2        PRBTYP
      LOGICAL            LINOBJ, SINGLR, UNITGZ, UNITQ
      INTEGER            ISTATE(*), KX(N)
      DOUBLE PRECISION   BL(*), BU(*), A(NROWA,*),
     $                   RES(*), FEATOL(*)
      DOUBLE PRECISION   GQ(N), CQ(*), R(NROWR,*), X(N), WTINF(*),
     $                   ZY(NQ,*)
      DOUBLE PRECISION   WRK(N)

************************************************************************
*     LSGSET  finds the number and weighted sum of infeasibilities for
*     the bounds and linear constraints.   An appropriate transformed
*     gradient vector is returned in  GQ.
*
*     Positive values of  ISTATE(j)  will not be altered.  These mean
*     the following...
*
*               1             2           3
*           a'x = bl      a'x = bu     bl = bu
*
*     Other values of  ISTATE(j)  will be reset as follows...
*           a'x lt bl     a'x gt bu     a'x free
*              - 2           - 1           0
*
*     If  x  is feasible,  LSGSET computes the vector Q(free)'g(free),
*     where  g  is the gradient of the the sum of squares plus the
*     linear term.  The matrix Q is of the form
*                    ( Q(free)  0       ),
*                    (   0      I(fixed))
*     where  Q(free)  is the orthogonal factor of  A(free)  and  A  is
*     the matrix of constraints in the working set.  The transformed
*     gradients are stored in GQ.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original version written 31-October-1984.
*     Level 2 Blas added 11-June-1986.
*     This version of LSGSET dated 24-June-1986.
************************************************************************
      EXTERNAL           DDOT  , IDRANK
      INTRINSIC          DABS   , DMAX1   , MIN
      PARAMETER        ( ZERO = 0.0D+0, HALF = 0.5D+0, ONE = 1.0D+0 )

      BIGUPP =   BIGBND
      BIGLOW = - BIGBND

      NUMINF =   0
      SUMINF =   ZERO
      CALL DLOAD ( N, ZERO, GQ, 1 )

      DO 200 J = 1, N+NCLIN
         IF (ISTATE(J) .LE. 0) THEN
            FEASJ  = FEATOL(J)
            IF (J .LE. N) THEN
               CTX = X(J)
            ELSE
               K   = J - N
               CTX = DDOT  ( N, A(K,1), NROWA, X, 1 )
            END IF
            ISTATE(J) = 0

*           See if the lower bound is violated.

            IF (BL(J) .GT. BIGLOW) THEN
               S = BL(J) - CTX
               IF (S     .GT. FEASJ ) THEN
                  ISTATE(J) = - 2
                  WEIGHT    = - WTINF(J)
                  GO TO 160
               END IF
            END IF

*           See if the upper bound is violated.

            IF (BU(J) .GE. BIGUPP) GO TO 200
            S = CTX - BU(J)
            IF (S     .LE. FEASJ ) GO TO 200
            ISTATE(J) = - 1
            WEIGHT    =   WTINF(J)

*           Add the infeasibility.

  160       NUMINF = NUMINF + 1
            SUMINF = SUMINF + DABS( WEIGHT ) * S
            IF (J .LE. N) THEN
               GQ(J) = WEIGHT
            ELSE
               CALL DAXPY ( N, WEIGHT, A(K,1), NROWA, GQ, 1 )
            END IF
         END IF
  200 CONTINUE

*     ------------------------------------------------------------------
*     Install  GQ,  the transformed gradient.
*     ------------------------------------------------------------------
      SINGLR = .FALSE.
      UNITGZ = .TRUE.

      IF (NUMINF .GT. 0) THEN
         CALL CMQMUL( 6, N, NZ, NFREE, NQ, UNITQ, KX, GQ, ZY, WRK )
      ELSE IF (NUMINF .EQ. 0  .AND.  PRBTYP .EQ. 'FP') THEN
         CALL DLOAD ( N, ZERO, GQ, 1 )
      ELSE

*        Ready for the Optimality Phase.
*        Set NZ1 so that Rz1 is nonsingular.

         IF (NRANK .EQ. 0) THEN
            IF (LINOBJ) THEN
               CALL DCOPY ( N, CQ, 1, GQ, 1 )
            ELSE
               CALL DLOAD ( N, ZERO, GQ, 1 )
            END IF
            NZ1    = 0
         ELSE

*           Compute  GQ = - R' * (transformed residual)

            CALL DCOPY ( NRANK, RES, 1, GQ, 1 )
            CALL DSCAL ( NRANK, (-ONE), GQ, 1 )
            CALL DTRMV ( 'U', 'T', 'N', NRANK, R, NROWR, GQ, 1 )
            IF (NRANK .LT. N)
     $         CALL DGEMV( 'T', NRANK, N-NRANK, -ONE,R(1,NRANK+1),NROWR,
     $                      RES, 1, ZERO, GQ(NRANK+1), 1 )

            IF (LINOBJ) CALL DXPY ( N, CQ, 1, GQ, 1 )
            UNITGZ = .FALSE.
            NZ1    = IDRANK( MIN(NRANK, NZ), R, NROWR+1, TOLRNK )
         END IF
      END IF

      RETURN

*     End of  LSGSET.

      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE LSKEY ( NOUT, BUFFER, KEY )

*     IMPLICIT           DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION  RDUMMY, ZERO
      INTEGER           MXPARM, MAXKEY, MAXTIE, MAXTOK, MAXTYP, IDUMMY
      INTEGER           IDBGLS, ITMAX1, ITMAX2, LCRASH, LDBGLS, LPROB, 
     &                  MSGLS,  NN,     NNCLIN, NPROB,  IPADLS, NOUT,   
     &                  I,      LENBUF, NTOKEN, LOC1,   LOC2,   IVALUE, 
     &                  LOC3
      DOUBLE PRECISION  BIGBND, BIGDX,  BNDLOW, BNDUPP, TOLACT, TOLFEA, 
     &                  TOLRNK, RPADLS, RVALUE
      CHARACTER*(*)      BUFFER

************************************************************************
*     LSKEY   decodes the option contained in  BUFFER  in order to set
*     a parameter value in the relevant element of  IPRMLS  or  RPRMLS.
*
*
*     Input:
*
*     NOUT   A unit number for printing error messages.
*            NOUT  must be a valid unit.
*
*     Output:
*
*     KEY    The first keyword contained in BUFFER.
*
*
*     LSKEY  calls OPNUMB and the subprograms
*                 LOOKUP, SCANNR, TOKENS, UPCASE
*     (now called OPLOOK, OPSCAN, OPTOKN, OPUPPR)
*     supplied by Informatics General, Inc., Palo Alto, California.
*
*     Systems Optimization Laboratory, Stanford University.
*     This version dated Jan 22, 1986.
************************************************************************
*-----------------------------------------------------------------------
      PARAMETER         (MXPARM = 30)
      INTEGER            IPRMLS(MXPARM), IPSVLS
      DOUBLE PRECISION   RPRMLS(MXPARM), RPSVLS

      COMMON    /LSPAR1/ IPSVLS(MXPARM),
     $                   IDBGLS, ITMAX1, ITMAX2, LCRASH, LDBGLS, LPROB ,
     $                   MSGLS , NN    , NNCLIN, NPROB , IPADLS(20)

      COMMON    /LSPAR2/ RPSVLS(MXPARM),
     $                   BIGBND, BIGDX , BNDLOW, BNDUPP, TOLACT, TOLFEA,
     $                   TOLRNK, RPADLS(23)

      EQUIVALENCE       (IPRMLS(1), IDBGLS), (RPRMLS(1), BIGBND)

      SAVE      /LSPAR1/, /LSPAR2/
*-----------------------------------------------------------------------

      EXTERNAL           OPNUMB
      LOGICAL            FIRST , MORE  , NUMBER, OPNUMB, SORTED
      SAVE               FIRST

      PARAMETER         (     MAXKEY = 27,  MAXTIE = 10,   MAXTOK = 10,
     $                        MAXTYP = 16)
      CHARACTER*16       KEYS(MAXKEY), TIES(MAXTIE), TOKEN(MAXTOK),
     $                   TYPE(MAXTYP)
      CHARACTER*16       KEY, KEY2, KEY3, VALUE

      PARAMETER         (IDUMMY = -11111,  RDUMMY = -11111.0,
     $                   SORTED = .TRUE.,  ZERO   =  0.0     )

      DATA                FIRST
     $                  /.TRUE./
      DATA   KEYS
     $ / 'BEGIN           ',
     $   'COLD            ', 'CONSTRAINTS     ', 'CRASH           ',
     $   'DEBUG           ', 'DEFAULTS        ', 'END             ',
     $   'FEASIBILITY     ', 'HOT             ', 'INFINITE        ',
     $   'IPRMLS          ', 'ITERATIONS      ', 'ITERS:ITERATIONS',
     $   'ITNS :ITERATIONS', 'LINEAR          ', 'LIST            ',
     $   'LOWER           ', 'NOLIST          ', 'OPTIMALITY      ',
     $   'PRINT           ', 'PROBLEM         ', 'RANK            ',
     $   'RPRMLS          ', 'START           ', 'UPPER           ',
     $   'VARIABLES       ', 'WARM            '/

      DATA   TIES
     $ / 'BOUND           ', 'CONSTRAINTS     ',
     $   'NO              ', 'NO.      :NUMBER', 'NUMBER          ',
     $   'PHASE           ', 'STEP            ',
     $   'TOLERANCE       ', 'TYPE            ', 'YES             '/

      DATA   TYPE
     $ / 'FP              ',
     $   'LEAST       :LS1', 'LINEAR       :LP', 'LP              ',
     $   'LS          :LS1', 'LS1             ', 'LS2             ',
     $   'LS3             ', 'LS4             ', 'LSQ         :LS1',
     $   'QP          :QP2', 'QP1             ', 'QP2             ',
     $   'QP3             ', 'QP4             ', 'QUADRATIC   :QP2'/
*-----------------------------------------------------------------------

      IF (FIRST) THEN
         FIRST  = .FALSE.
         DO 10 I = 1, MXPARM
            IPRMLS(I) = IDUMMY
            RPRMLS(I) = RDUMMY
   10    CONTINUE
      END IF

*     Eliminate comments and empty lines.
*     A '*' appearing anywhere in BUFFER terminates the string.

      I      = INDEX( BUFFER, '*' )
      IF (I .EQ. 0) THEN
         LENBUF = LEN( BUFFER )
      ELSE
         LENBUF = I - 1
      END IF
      IF (LENBUF .LE. 0) THEN
         KEY = '*'
         GO TO 900
      END IF

*     ------------------------------------------------------------------
*     Extract up to MAXTOK tokens from the record.
*     NTOKEN returns how many were actually found.
*     KEY, KEY2, KEY3 are the first tokens if any, otherwise blank.
*     ------------------------------------------------------------------
      NTOKEN = MAXTOK
      CALL OPTOKN( BUFFER(1:LENBUF), NTOKEN, TOKEN )
      KEY    = TOKEN(1)
      KEY2   = TOKEN(2)
      KEY3   = TOKEN(3)

*     Certain keywords require no action.

      IF (KEY .EQ. ' '     .OR.  KEY .EQ. 'BEGIN' ) GO TO 900
      IF (KEY .EQ. 'LIST'  .OR.  KEY .EQ. 'NOLIST') GO TO 900
      IF (KEY .EQ. 'END'                          ) GO TO 900

*     Most keywords will have an associated integer or real value,
*     so look for it no matter what the keyword.

      I      = 1
      NUMBER = .FALSE.

   50 IF (I .LT. NTOKEN  .AND.  .NOT. NUMBER) THEN
         I      = I + 1
         VALUE  = TOKEN(I)
         NUMBER = OPNUMB( VALUE )
         GO TO 50
      END IF

      IF (NUMBER) THEN
         READ (VALUE, '(BN, E16.0)') RVALUE
      ELSE
         RVALUE = ZERO
      END IF

*     Convert the keywords to their most fundamental form
*     (upper case, no abbreviations).
*     SORTED says whether the dictionaries are in alphabetic order.
*     LOCi   says where the keywords are in the dictionaries.
*     LOCi = 0 signals that the keyword wasn't there.

      CALL OPLOOK( MAXKEY, KEYS, 1, SORTED, KEY , LOC1 )
      CALL OPLOOK( MAXTIE, TIES, 1, SORTED, KEY2, LOC2 )

*     ------------------------------------------------------------------
*     Decide what to do about each keyword.
*     The second keyword (if any) might be needed to break ties.
*     Some seemingly redundant testing of MORE is used
*     to avoid compiler limits on the number of consecutive ELSE IFs.
*     ------------------------------------------------------------------
      MORE   = .TRUE.
      IF (MORE) THEN
         MORE   = .FALSE.
         IF (KEY .EQ. 'COLD        ') THEN
            LCRASH = 0
         ELSE IF (KEY .EQ. 'CONSTRAINTS ') THEN
            NNCLIN = RVALUE
         ELSE IF (KEY .EQ. 'CRASH       ') THEN
            TOLACT = RVALUE
         ELSE IF (KEY .EQ. 'DEBUG       ') THEN
            LDBGLS = RVALUE
         ELSE IF (KEY .EQ. 'DEFAULTS    ') THEN
            DO 20 I = 1, MXPARM
               IPRMLS(I) = IDUMMY
               RPRMLS(I) = RDUMMY
   20       CONTINUE
         ELSE IF (KEY .EQ. 'FEASIBILITY ') THEN
              IF (KEY2.EQ. 'PHASE       ') ITMAX1 = RVALUE
              IF (KEY2.EQ. 'TOLERANCE   ') TOLFEA = RVALUE
              IF (LOC2.EQ.  0            ) WRITE(NOUT, 2320) KEY2
         ELSE
            MORE   = .TRUE.
         END IF
      END IF

      IF (MORE) THEN
         MORE   = .FALSE.
         IF (KEY .EQ. 'HOT         ') THEN
            LCRASH = 2
         ELSE IF (KEY .EQ. 'INFINITE    ') THEN
              IF (KEY2.EQ. 'BOUND       ') BIGBND = RVALUE * 0.99999
              IF (KEY2.EQ. 'STEP        ') BIGDX  = RVALUE
              IF (LOC2.EQ.  0            ) WRITE(NOUT, 2320) KEY2
         ELSE IF (KEY .EQ. 'IPRMLS      ') THEN
*           Allow things like  IPRMLS 21 = 100  to set IPRMLS(21) = 100
            IVALUE = RVALUE
            IF (IVALUE .GE. 1  .AND. IVALUE .LE. MXPARM) THEN
               READ (KEY3, '(BN, I16)') IPRMLS(IVALUE)
            ELSE
               WRITE(NOUT, 2400) IVALUE
            END IF
         ELSE IF (KEY .EQ. 'ITERATIONS  ') THEN
            ITMAX2 = RVALUE
         ELSE IF (KEY .EQ. 'LINEAR      ') THEN
            NNCLIN = RVALUE
         ELSE IF (KEY .EQ. 'LOWER       ') THEN
            BNDLOW = RVALUE
         ELSE
            MORE   = .TRUE.
         END IF
      END IF

      IF (MORE) THEN
         MORE   = .FALSE.
         IF      (KEY .EQ. 'OPTIMALITY  ') THEN
            ITMAX2 = RVALUE
         ELSE IF (KEY .EQ. 'PROBLEM     ') THEN
            IF      (KEY2 .EQ. 'NUMBER') THEN
               NPROB  = RVALUE
            ELSE IF (KEY2 .EQ. 'TYPE  ') THEN

*              Recognize     Problem type = LP     etc.

               CALL OPLOOK( MAXTYP, TYPE, 1, SORTED, KEY3, LOC3 )
               IF (KEY3 .EQ. 'FP' ) LPROB = 1
               IF (KEY3 .EQ. 'LP' ) LPROB = 2
               IF (KEY3 .EQ. 'QP1') LPROB = 3
               IF (KEY3 .EQ. 'QP2') LPROB = 4
               IF (KEY3 .EQ. 'QP3') LPROB = 5
               IF (KEY3 .EQ. 'QP4') LPROB = 6
               IF (KEY3 .EQ. 'LS1') LPROB = 7
               IF (KEY3 .EQ. 'LS2') LPROB = 8
               IF (KEY3 .EQ. 'LS3') LPROB = 9
               IF (KEY3 .EQ. 'LS4') LPROB = 10
               IF (LOC3 .EQ.  0  ) WRITE(NOUT, 2330) KEY3
            ELSE
               WRITE(NOUT, 2320) KEY2
            END IF
         ELSE
            MORE   = .TRUE.
         END IF
      END IF

      IF (MORE) THEN
         MORE   = .FALSE.
         IF      (KEY .EQ. 'PRINT       ') THEN
            MSGLS  = RVALUE
         ELSE IF (KEY .EQ. 'RANK        ') THEN
            TOLRNK = RVALUE
         ELSE IF (KEY .EQ. 'RPRMLS      ') THEN
*           Allow things like  RPRMLS 21 = 2  to set RPRMLS(21) = 2.0
            IVALUE = RVALUE
            IF (IVALUE .GE. 1  .AND. IVALUE .LE. MXPARM) THEN
               READ (KEY3, '(BN, E16.0)') RPRMLS(IVALUE)
            ELSE
               WRITE(NOUT, 2400) IVALUE
            END IF
         ELSE IF (KEY .EQ. 'START       ') THEN
            IDBGLS = RVALUE
         ELSE IF (KEY .EQ. 'UPPER       ') THEN
            BNDUPP = RVALUE
         ELSE IF (KEY .EQ. 'VARIABLES   ') THEN
            NN     = RVALUE
         ELSE IF (KEY .EQ. 'WARM        ') THEN
            LCRASH = 1
         ELSE
            WRITE(NOUT, 2300) KEY
         END IF
      END IF

  900 RETURN

 2300 FORMAT(' XXX  Keyword not recognized:         ', A)
 2320 FORMAT(' XXX  Second keyword not recognized:  ', A)
 2330 FORMAT(' XXX  Third  keyword not recognized:  ', A)
 2400 FORMAT(' XXX  The PARM subscript is out of range:', I10)

*     End of LSKEY

      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE LSLOC ( LPROB, N, NCLIN, LITOTL, LWTOTL )

*     IMPLICIT           DOUBLE PRECISION(A-H,O-Z)
      INTEGER            LENLS,  LDBG,   NOUT,   ILSDBG, LENNAM, NROWT,  
     &                   NCOLT,  NQ,     LOCLS,  LPROB,  N,      NCLIN,  
     &                   LITOTL, LWTOTL, MINIW,  LENT,   LENZY,  LENCQ,  
     &                   LENRES, LKACTV, MINW,   LANORM, LAP,    LPX,    
     &                   LRES,   LHZ,    LGQ,    LCQ,    LRLAM,  LT,     
     &                   LZY,    LWTINF, LWRK,   LFEATL, LRES0

************************************************************************
*     LSLOC   allocates the addresses of the work arrays for  LSCORE.
*
*     Note that the arrays  ( GQ, CQ )  and  ( RES, RES0, HZ )  lie in
*     contiguous areas of workspace.
*     RES, RES0 and HZ are not needed for LP.
*     CQ is defined when the objective has an explicit linear term.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original version written  29-October-1984.
*     This version of LSLOC dated 16-February-1986.
************************************************************************
      COMMON    /SOL1CM/ NOUT
      COMMON    /SOL3CM/ LENNAM, NROWT, NCOLT, NQ

      PARAMETER        ( LENLS = 20 )
      COMMON    /SOL1LS/ LOCLS(LENLS)

      LOGICAL            LSDBG
      PARAMETER        ( LDBG = 5 )
      COMMON    /LSDEBG/ ILSDBG(LDBG), LSDBG

      MINIW     = LITOTL + 1
      MINW      = LWTOTL + 1


*     Assign array lengths that depend upon the problem dimensions.

      IF (NCLIN .EQ. 0) THEN
         LENT      = 0
         LENZY     = 0
      ELSE
         LENT  = NROWT*NCOLT
         LENZY = NQ   *NQ
      END IF

      LENCQ  = 0
      IF (LPROB .EQ. 2*(LPROB/2)) LENCQ  = N
      LENRES = 0
      IF (LPROB .GT. 2          ) LENRES = N

      LKACTV    = MINIW
      MINIW     = LKACTV + N

      LANORM    = MINW
      LAP       = LANORM + NCLIN
      LPX       = LAP    + NCLIN
      LGQ       = LPX    + N
      LCQ       = LGQ    + N
      LRES      = LCQ    + LENCQ
      LRES0     = LRES   + LENRES
      LHZ       = LRES0  + LENRES
      LRLAM     = LHZ    + LENRES
      LT        = LRLAM  + N
      LZY       = LT     + LENT
      LWTINF    = LZY    + LENZY
      LWRK      = LWTINF + N  + NCLIN
      LFEATL    = LWRK   + N  + NCLIN
      MINW      = LFEATL + N  + NCLIN

      LOCLS( 1) = LKACTV
      LOCLS( 2) = LANORM
      LOCLS( 3) = LAP
      LOCLS( 4) = LPX
      LOCLS( 5) = LRES
      LOCLS( 6) = LRES0
      LOCLS( 7) = LHZ
      LOCLS( 8) = LGQ
      LOCLS( 9) = LCQ
      LOCLS(10) = LRLAM
      LOCLS(11) = LT
      LOCLS(12) = LZY
      LOCLS(13) = LWTINF
      LOCLS(14) = LWRK
      LOCLS(15) = LFEATL

      LITOTL    = MINIW - 1
      LWTOTL    = MINW  - 1

      RETURN

*     End of  LSLOC .

      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE LSMOVE( HITCON, HITLOW, LINOBJ, UNITGZ,
     $                   NCLIN, NRANK, NZ1,
     $                   N, NROWR, JADD, NUMINF,
     $                   ALFA, CTP, CTX, XNORM,
     $                   AP, AX, BL, BU, GQ, HZ, P, RES,
     $                   R, X, WORK )

*     IMPLICIT           DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION  ONE,    ZERO
      INTEGER           LDBG
      INTEGER           NOUT,   ILSDBG, NCLIN,  NRANK,  NZ1,    N,      
     &                  NROWR,  JADD,   NUMINF
      DOUBLE PRECISION  ALFA,   CTP,    CTX,    XNORM,  BND,    DNRM2
      LOGICAL            HITCON, HITLOW, LINOBJ, UNITGZ
      DOUBLE PRECISION   AP(*), AX(*), BL(*), BU(*), GQ(*), HZ(*),
     $                   P(N), RES(*), R(NROWR,*), X(N)
      DOUBLE PRECISION   WORK(*)

************************************************************************
*     LSMOVE  changes X to X + ALFA*P and updates CTX, AX, RES and GQ
*     accordingly.
*
*     If a bound was added to the working set,  move X exactly on to it,
*     except when a negative step was taken (CMALF may have had to move
*     to some other closer constraint.)
*
*     Systems Optimization Laboratory, Stanford University.
*     Original version written 27-December-1985.
*     Level 2 BLAS added 11-June-1986.
*     This version of LSMOVE dated 11-June-1986.
************************************************************************
      COMMON    /SOL1CM/ NOUT

      LOGICAL            LSDBG
      PARAMETER         (LDBG = 5)
      COMMON    /LSDEBG/ ILSDBG(LDBG), LSDBG

      EXTERNAL           DDOT  , DNRM2
      INTRINSIC          DABS   , MIN
      PARAMETER        ( ZERO  = 0.0D+0, ONE = 1.0D+0 )

      CALL DAXPY ( N, ALFA, P, 1, X, 1 )
      IF (LINOBJ) CTX = CTX + ALFA*CTP

      IF (HITCON  .AND.  JADD .LE. N) THEN
         BND = BU(JADD)
         IF (HITLOW) BND = BL(JADD)
         IF (ALFA .GE. ZERO) X(JADD) = BND
      END IF
      XNORM  = DNRM2 ( N, X, 1 )

      IF (NCLIN .GT. 0)
     $   CALL DAXPY ( NCLIN, ALFA, AP, 1, AX, 1 )

      IF (NZ1 .LE. NRANK) THEN
         IF (UNITGZ) THEN
            RES(NZ1) = RES(NZ1) - ALFA*HZ(NZ1)
         ELSE
            CALL DAXPY ( NZ1, (-ALFA), HZ, 1, RES, 1  )
         END IF

         IF (NUMINF .EQ. 0) THEN

*           Update the transformed gradient GQ so that
*           GQ = GQ + ALFA*R'( HZ ).
*                            ( 0  )

            IF (UNITGZ) THEN
               CALL DAXPY ( N-NZ1+1, ALFA*HZ(NZ1), R(NZ1,NZ1), NROWR,
     $                                             GQ(NZ1)   , 1      )
            ELSE
               CALL DCOPY ( NZ1, HZ, 1, WORK, 1 )
               CALL DTRMV ( 'U', 'T', 'N', NZ1, R, NROWR, WORK, 1 )
               IF (NZ1 .LT. N)
     $            CALL DGEMV ( 'T', NZ1, N-NZ1, ONE, R(1,NZ1+1), NROWR,
     $                         HZ, 1, ZERO, WORK(NZ1+1), 1 )
               CALL DAXPY ( N, ALFA, WORK, 1, GQ, 1 )
            END IF
         END IF
      END IF

      RETURN

*     End of  LSMOVE.

      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE LSMULS( PRBTYP,
     $                   MSGLVL, N, NACTIV, NFREE,
     $                   NROWA, NROWT, NUMINF, NZ, NZ1,
     $                   ISTATE, KACTIV, KX, DINKY,
     $                   JSMLST, KSMLST, JINF, JTINY,
     $                   JBIGST, KBIGST, TRULAM,
     $                   A, ANORMS, GQ, RLAMDA, T, WTINF )

*     IMPLICIT           DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION  ONE,    ZERO
      INTEGER           MSGLVL, NACTIV, NFREE,  NROWA,  NROWT,  NUMINF, 
     &                  NZ,     NZ1,    JSMLST, KSMLST, JINF,   JTINY,  
     &                  JBIGST, KBIGST, NFIXED, J,      K
      INTEGER           NOUT,   ILSDBG, N,      L,      I,      
     &                  IS
      DOUBLE PRECISION  DINKY,  TRULAM, TINYLM, SMLLST, BIGGST, RLAM,   
     &                  BLAM,   ANORMJ, SCDLAM
      INTEGER           LDBG
      CHARACTER*2        PRBTYP
      INTEGER            ISTATE(*), KACTIV(N), KX(N)
      DOUBLE PRECISION   A(NROWA,*), ANORMS(*),
     $                   GQ(N), RLAMDA(N), T(NROWT,*), WTINF(*)

************************************************************************
*     LSMULS  first computes the Lagrange multiplier estimates for the
*     given working set.  It then determines the values and indices of
*     certain significant multipliers.  In this process, the multipliers
*     for inequalities at their upper bounds are adjusted so that a
*     negative multiplier for an inequality constraint indicates non-
*     optimality.  All adjusted multipliers are scaled by the 2-norm
*     of the associated constraint row.  In the following, the term
*     minimum refers to the ordering of numbers on the real line,  and
*     not to their magnitude.
*
*     JSMLST  is the index of the minimum of the set of adjusted
*             multipliers with values less than  - DINKY.  A negative
*             JSMLST defines the index in Q'g of the artificial
*             constraint to be deleted.
*     KSMLST  marks the position of general constraint JSMLST in KACTIV.
*
*     JBIGST  is the index of the largest of the set of adjusted
*             multipliers with values greater than (1 + DINKY).
*     KBIGST  marks its position in KACTIV.
*
*     On exit,  elements 1 thru NACTIV of RLAMDA contain the unadjusted
*     multipliers for the general constraints.  Elements NACTIV onwards
*     of RLAMDA contain the unadjusted multipliers for the bounds.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original version written 31-October-1984.
*     This version of LSMULS dated  30-June-1986.
************************************************************************
      COMMON    /SOL1CM/ NOUT

      LOGICAL            LSDBG
      PARAMETER         (LDBG = 5)
      COMMON    /LSDEBG/ ILSDBG(LDBG), LSDBG

      INTRINSIC          DABS, MIN
      PARAMETER        ( ZERO   =0.0D+0,ONE    =1.0D+0 )

      NFIXED =   N - NFREE

      JSMLST =   0
      KSMLST =   0
      SMLLST = - DINKY

      TINYLM =   DINKY
      JTINY  =   0

      JBIGST =   0
      KBIGST =   0
      BIGGST =   ONE + DINKY

      IF (NZ1 .LT. NZ) THEN
*        ---------------------------------------------------------------
*        Compute JSMLST for the artificial constraints.
*        ---------------------------------------------------------------
         DO 100 J = NZ1+1, NZ
            RLAM = - DABS( GQ(J) )
            IF (RLAM .LT. SMLLST) THEN
               SMLLST =   RLAM
               JSMLST = - J
            ELSE IF (RLAM .LT. TINYLM) THEN
               TINYLM =   RLAM
               JTINY  =   J
            END IF
  100    CONTINUE

         IF (MSGLVL .GE. 20)
     $      WRITE (NOUT, 1000) (GQ(K), K=NZ1+1,NZ)

      END IF

*     ------------------------------------------------------------------
*     Compute JSMLST for regular constraints and temporary bounds.
*     ------------------------------------------------------------------
*     First, compute the Lagrange multipliers for the general
*     constraints in the working set, by solving  T'*lamda = Y'g.

      IF (N .GT. NZ)
     $   CALL DCOPY ( N-NZ, GQ(NZ+1), 1, RLAMDA, 1 )
      IF (NACTIV .GT. 0)
     $   CALL CMTSOL( 2, NROWT, NACTIV, T(1,NZ+1), RLAMDA )

*     -----------------------------------------------------------------
*     Now set elements NACTIV, NACTIV+1,... of  RLAMDA  equal to
*     the multipliers for the bound constraints.
*     -----------------------------------------------------------------
      DO 190 L = 1, NFIXED
         J     = KX(NFREE+L)
         BLAM  = RLAMDA(NACTIV+L)
         DO 170 K = 1, NACTIV
            I    = KACTIV(K)
            BLAM = BLAM - A(I,J)*RLAMDA(K)
  170    CONTINUE
         RLAMDA(NACTIV+L) = BLAM
  190 CONTINUE

*     -----------------------------------------------------------------
*     Find JSMLST and KSMLST.
*     -----------------------------------------------------------------
      DO 330 K = 1, N - NZ
         IF (K .GT. NACTIV) THEN
            J = KX(NZ+K)
         ELSE
            J = KACTIV(K) + N
         END IF

         IS   = ISTATE(J)

         I    = J - N
         IF (J .LE. N) ANORMJ = ONE
         IF (J .GT. N) ANORMJ = ANORMS(I)

         RLAM = RLAMDA(K)

*        Change the sign of the estimate if the constraint is in
*        the working set at its upper bound.

         IF (IS .EQ. 2) RLAM =      - RLAM
         IF (IS .EQ. 3) RLAM =   DABS( RLAM )
         IF (IS .EQ. 4) RLAM = - DABS( RLAM )

         IF (IS .NE. 3) THEN
            SCDLAM = RLAM * ANORMJ
            IF      (SCDLAM .LT. SMLLST) THEN
               SMLLST = SCDLAM
               JSMLST = J
               KSMLST = K
            ELSE IF (SCDLAM .LT. TINYLM) THEN
               TINYLM = SCDLAM
               JTINY  = J
            END IF
         END IF

         IF (NUMINF .GT. 0  .AND.  J .GT. JINF) THEN
            SCDLAM = RLAM/WTINF(J)
            IF (SCDLAM .GT. BIGGST) THEN
               BIGGST = SCDLAM
               TRULAM = RLAMDA(K)
               JBIGST = J
               KBIGST = K
            END IF
         END IF
  330 CONTINUE

*     -----------------------------------------------------------------
*     If required, print the multipliers.
*     -----------------------------------------------------------------
      IF (MSGLVL .GE. 20) THEN
         IF (NFIXED .GT. 0)
     $      WRITE (NOUT, 1100) PRBTYP, (KX(NFREE+K),
     $                         RLAMDA(NACTIV+K), K=1,NFIXED)
         IF (NACTIV .GT. 0)
     $      WRITE (NOUT, 1200) PRBTYP, (KACTIV(K),
     $                         RLAMDA(K), K=1,NACTIV)
      END IF

      IF (LSDBG  .AND.  ILSDBG(1) .GT. 0) THEN
         WRITE (NOUT, 9000) JSMLST, SMLLST, KSMLST
         WRITE (NOUT, 9100) JBIGST, BIGGST, KBIGST
         WRITE (NOUT, 9200) JTINY , TINYLM
      END IF

      RETURN

 1000 FORMAT(/ ' Multipliers for the artificial constraints        '
     $       / 4(5X, 1PE11.2))
 1100 FORMAT(/ ' Multipliers for the ', A2, ' bound  constraints   '
     $       / 4(I5, 1PE11.2))
 1200 FORMAT(/ ' Multipliers for the ', A2, ' linear constraints   '
     $       / 4(I5, 1PE11.2))
 9000 FORMAT(/ ' //LSMULS//  JSMLST     SMLLST     KSMLST (Scaled) '
     $       / ' //LSMULS//  ', I6, 1PE11.2, 5X, I6 )
 9100 FORMAT(  ' //LSMULS//  JBIGST     BIGGST     KBIGST (Scaled) '
     $       / ' //LSMULS//  ', I6, 1PE11.2, 5X, I6 )
 9200 FORMAT(  ' //LSMULS//   JTINY     TINYLM                     '
     $       / ' //LSMULS//  ', I6, 1PE11.2)

*     End of  LSMULS.

      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE LSOPTN( STRING )
      INTEGER           NOUT
      CHARACTER*(*)      STRING

************************************************************************
*     LSOPTN  loads the option supplied in  STRING  into the relevant
*     element of  IPRMLS  or  RPRMLS.
************************************************************************

      LOGICAL             NEWOPT
      COMMON     /SOL3LS/ NEWOPT
      SAVE       /SOL3LS/

      DOUBLE PRECISION    WMACH(15)
      COMMON     /SOLMCH/ WMACH
      SAVE       /SOLMCH/

      EXTERNAL            MCHPAR
      CHARACTER*16        KEY
      CHARACTER*72        BUFFER
      LOGICAL             FIRST , PRNT
      SAVE                FIRST , NOUT  , PRNT
      DATA                FIRST /.TRUE./

*     If first time in, set  NOUT.
*     NEWOPT  is true first time into  LSOPTN
*     and just after a call to  LSSOL.
*     PRNT    is set to true whenever  NEWOPT  is true.

      IF (FIRST) THEN
         FIRST  = .FALSE.
         NEWOPT = .TRUE.
         CALL MCHPAR()
         NOUT   =  WMACH(11)
      END IF
      BUFFER = STRING

*     Call  LSKEY   to decode the option and set the parameter value.
*     If NEWOPT is true, reset PRNT and test specially for NOLIST.

      IF (NEWOPT) THEN
         NEWOPT = .FALSE.
         PRNT   = .TRUE.
         CALL LSKEY ( NOUT, BUFFER, KEY )

         IF (KEY .EQ. 'NOLIST') THEN
            PRNT   = .FALSE.
         END IF
      ELSE
         CALL LSKEY ( NOUT, BUFFER, KEY )

         IF (KEY .EQ.   'LIST') PRNT = .TRUE.
         IF (KEY .EQ. 'NOLIST') PRNT = .FALSE.
      END IF

      RETURN

*     End of  LSOPTN.

      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE LSPRT ( PRBTYP, PRNT1, ISDEL, ITER, JADD, JDEL,
     $                   MSGLVL, NACTIV, NFREE, N, NCLIN,
     $                   NRANK, NROWR, NROWT, NZ, NZ1, ISTATE,
     $                   ALFA, CONDRZ, CONDT, GFNORM, GZNORM, GZ1NRM,
     $                   NUMINF, SUMINF, CTX, SSQ,
     $                   AX, R, T, X, WORK )

*     IMPLICIT           DOUBLE PRECISION(A-H,O-Z)
      INTEGER           LDBG,   NOUT,   ILSDBG, ISDEL,  ITER,   JADD,   
     &                  JDEL,   MSGLVL, NACTIV, NFREE,  N,      NCLIN,  
     &                  NRANK,  NROWR,  NROWT,  NZ,     NUMINF, KDEL,   
     &                  KADD,   J,      NZ1,    K
      DOUBLE PRECISION  ALFA,   CONDRZ, CONDT,  GFNORM, GZ1NRM, SUMINF, 
     &                  CTX,    SSQ,    OBJ,    GZNORM
      CHARACTER*2        PRBTYP
      LOGICAL            PRNT1
      INTEGER            ISTATE(*)
      DOUBLE PRECISION   AX(*), R(NROWR,*), T(NROWT,*), X(N)
      DOUBLE PRECISION   WORK(N)

************************************************************************
*  LSPRT  prints various levels of output for  LSCORE.
*
*           Msg        Cumulative result
*           ---        -----------------
*
*        le   0        no output.
*
*        eq   1        nothing now (but full output later).
*
*        eq   5        one terse line of output.
*
*        ge  10        same as 5 (but full output later).
*
*        ge  20        constraint status,  x  and  Ax.
*
*        ge  30        diagonals of  T  and  R.
*
*
*  Debug printing is performed depending on the logical variable  LSDBG.
*  LSDBG  is set true when  IDBG  major iterations have been performed.
*  At this point,  printing is done according to a string of binary
*  digits of the form  SVT  (stored in the integer array  ILSDBG).
*
*  S  set 'on'  gives information from the maximum step routine  CMALF.
*  V  set 'on'  gives various vectors in  LSCORE  and its auxiliaries.
*  T  set 'on'  gives a trace of which routine was called and an
*               indication of the progress of the run.
*
*  Systems Optimization Laboratory, Stanford University.
*  Original version written 31-October-1984.
*  This version of LSPRT dated 14-January-1985.
************************************************************************
      COMMON    /SOL1CM/ NOUT

      LOGICAL            LSDBG
      PARAMETER         (LDBG = 5)
      COMMON    /LSDEBG/ ILSDBG(LDBG), LSDBG

      CHARACTER*2        LADD, LDEL
      CHARACTER*2        LSTATE(0:5)
      DATA               LSTATE(0), LSTATE(1), LSTATE(2)
     $                  /'  '     , 'L '     , 'U '     /
      DATA               LSTATE(3), LSTATE(4), LSTATE(5)
     $                  /'E '     , 'T '     , 'Z '     /

      IF (MSGLVL .GE. 15) WRITE (NOUT, 1000) PRBTYP, ITER

      IF (MSGLVL .GE. 5) THEN
         IF      (JDEL .GT. 0) THEN
            KDEL =   ISDEL
         ELSE IF (JDEL .LT. 0) THEN
            JDEL = - JDEL
            KDEL =   5
         ELSE
            KDEL =   0
         END IF

         IF (JADD .GT. 0) THEN
            KADD = ISTATE(JADD)
         ELSE
            KADD = 0
         END IF

         LDEL   = LSTATE(KDEL)
         LADD   = LSTATE(KADD)

         IF (NUMINF .GT. 0) THEN
            OBJ    = SUMINF
         ELSE
            OBJ    = SSQ + CTX
         END IF

*        ---------------------------------------------------------------
*        Print the terse line.
*        ---------------------------------------------------------------
         IF (NRANK .EQ. 0) THEN
            IF (PRNT1  .OR.  MSGLVL .GE. 15) WRITE (NOUT, 1100)
            WRITE (NOUT, 1200) ITER, JDEL, LDEL, JADD, LADD,
     $                         ALFA, NUMINF, OBJ, N-NFREE, NACTIV,
     $                         NZ, NZ1, GFNORM, GZ1NRM, CONDT
         ELSE
            IF (PRNT1  .OR.  MSGLVL .GE. 15) WRITE (NOUT, 1110)
            WRITE (NOUT, 1200) ITER, JDEL, LDEL, JADD, LADD,
     $                         ALFA, NUMINF, OBJ, N-NFREE, NACTIV,
     $                         NZ, NZ1, GFNORM, GZ1NRM, CONDT, CONDRZ
         END IF

         IF (MSGLVL .GE. 20) THEN
            WRITE (NOUT, 2000) PRBTYP
            WRITE (NOUT, 2100) (X(J) , ISTATE(J)  ,  J=1,N)
            IF (NCLIN .GT. 0)
     $      WRITE (NOUT, 2200) (AX(K), ISTATE(N+K), K=1,NCLIN )

            IF (MSGLVL .GE. 30) THEN
*              ---------------------------------------------------------
*              Print the diagonals of  T  and  R.
*              ---------------------------------------------------------
               IF (NACTIV .GT. 0) THEN
                  CALL DCOPY ( NACTIV, T(NACTIV,NZ+1), NROWT-1, WORK,1 )
                  WRITE (NOUT, 3000) PRBTYP, (WORK(J), J=1,NACTIV)
               END IF
               IF (NRANK  .GT. 0)
     $            WRITE (NOUT, 3100) PRBTYP, (R(J,J) , J=1,NRANK )
            END IF
            WRITE (NOUT, 5000)
         END IF
      END IF

      PRNT1 = .FALSE.

      RETURN

 1000 FORMAT(/// ' ', A2, ' iteration', I5
     $         / ' =================' )
 1100 FORMAT(// '  Itn Jdel  Jadd      Step',
     $          ' Ninf  Sinf/Objective', '  Bnd', '  Lin', '    Nz',
     $          '   Nz1   Norm Gf  Norm Gz1   Cond T' )
 1110 FORMAT(// '  Itn Jdel  Jadd      Step',
     $          ' Ninf  Sinf/Objective', '  Bnd', '  Lin', '    Nz',
     $          '   Nz1   Norm Gf  Norm Gz1   Cond T Cond Rz1' )
 1200 FORMAT(I5, I5, A1, I5, A1, 1PE9.1, I5, 1X, 1PE15.6, 2I5,
     $       2I6, 1P2E10.2, 1P2E9.1 )
 2000 FORMAT(/ ' Values and status of the ', A2, ' constraints'
     $       / ' ---------------------------------------' )
 2100 FORMAT(/ ' Variables...'                 /   (1X, 5(1PE15.6, I5)))
 2200 FORMAT(/ ' General linear constraints...'/   (1X, 5(1PE15.6, I5)))
 3000 FORMAT(/ ' Diagonals of ' , A2,' working set factor T'/(1P5E15.6))
 3100 FORMAT(/ ' Diagonals of ' , A2, ' triangle R         '/(1P5E15.6))
 5000 FORMAT(/// ' ---------------------------------------------------',
     $           '--------------------------------------------' )

*     End of  LSPRT .

      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE LSSETX( LINOBJ, ROWERR, UNITQ,
     $                   NCLIN, NACTIV, NFREE, NRANK, NZ,
     $                   N, NCTOTL, NQ, NROWA, NROWR, NROWT,
     $                   ISTATE, KACTIV, KX,
     $                   JMAX, ERRMAX, CTX, XNORM,
     $                   A, AX, BL, BU, CQ, RES, RES0, FEATOL,
     $                   R, T, X, ZY, P, WORK )

*     IMPLICIT           DOUBLE PRECISION (A-H,O-Z)
      INTEGER           LDBG,   NTRY,   NOUT,   ILSDBG, NCLIN,  NACTIV, 
     &                  NFREE,  NRANK,  NZ,     N,      NCTOTL, NQ,     
     &                  NROWA,  NROWR,  NROWT,  JMAX,   K,      J,      
     &                  IS,     I,      KTRY,   IDAMAX
      DOUBLE PRECISION  ONE,    ZERO,   ERRMAX, CTX,    XNORM,  
     &                  BND,    DDOT,   DNRM2
      LOGICAL            LINOBJ, ROWERR, UNITQ
      INTEGER            ISTATE(NCTOTL), KACTIV(N), KX(N)
      DOUBLE PRECISION   A(NROWA,*), AX(*), BL(NCTOTL), BU(NCTOTL),
     $                   CQ(*), RES(*), RES0(*), FEATOL(NCTOTL), P(N),
     $                   R(NROWR,*), T(NROWT,*), ZY(NQ,*), X(N)
      DOUBLE PRECISION   WORK(NCTOTL)

************************************************************************
*  LSSETX  computes the point on a working set that is closest to the
*  input vector  x  (in the least-squares sense).  The norm of  x, the
*  transformed residual vector  Pr - RQ'x,  and the constraint values
*  Ax  are also initialized.
*
*  If the computed point gives a row error of more than the feasibility
*  tolerance, an extra step of iterative refinement is used.  If  x  is
*  still infeasible,  the logical variable  ROWERR  is set.
*
*  Systems Optimization Laboratory, Stanford University.
*  Original version written 31-October-1984.
*  This version dated 29-December-1985.
************************************************************************
      COMMON    /SOL1CM/ NOUT

      LOGICAL            LSDBG
      PARAMETER         (LDBG = 5)
      COMMON    /LSDEBG/ ILSDBG(LDBG), LSDBG

      EXTERNAL           IDAMAX, DDOT
      INTRINSIC          DABS, MIN
      PARAMETER        ( NTRY  = 2 )
      PARAMETER        ( ZERO  = 0.0D+0, ONE = 1.0D+0 )

*     ------------------------------------------------------------------
*     Move  x  onto the simple bounds in the working set.
*     ------------------------------------------------------------------
      DO 100 K = NFREE+1, N
          J   = KX(K)
          IS  = ISTATE(J)
          BND = BL(J)
          IF (IS .GE. 2) BND  = BU(J)
          IF (IS .NE. 4) X(J) = BND
  100 CONTINUE

*     ------------------------------------------------------------------
*     Move  x  onto the general constraints in the working set.
*     We shall make  ntry  tries at getting acceptable row errors.
*     ------------------------------------------------------------------
      KTRY   = 1
      JMAX   = 1
      ERRMAX = ZERO

*     REPEAT
  200    IF (NACTIV .GT. 0) THEN

*           Set  work = residuals for constraints in the working set.
*           Solve for p, the smallest correction to x that gives a point
*           on the constraints in the working set.  Define  p = Y*(py),
*           where  py  solves the triangular system  T*(py) = residuals.

            DO 220 I = 1, NACTIV
               K   = KACTIV(I)
               J   = N + K
               BND = BL(J)
               IF (ISTATE(J) .EQ. 2) BND = BU(J)
               WORK(I) = BND - DDOT  ( N, A(K,1), NROWA, X, 1 )
  220       CONTINUE

            CALL CMTSOL( 1, NROWT, NACTIV, T(1,NZ+1), WORK )
            CALL DLOAD ( N, ZERO, P, 1 )
            CALL DCOPY ( NACTIV, WORK, 1, P(NZ+1), 1 )

            CALL CMQMUL( 2, N, NZ, NFREE, NQ, UNITQ, KX, P, ZY, WORK )
            CALL DXPY ( N, P, 1, X, 1 )
         END IF

*        ---------------------------------------------------------------
*        Compute the 2-norm of  x.
*        Initialize  Ax  for all the general constraints.
*        ---------------------------------------------------------------
         XNORM  = DNRM2 ( N, X, 1 )
         IF (NCLIN .GT. 0)
     $      CALL DGEMV ( 'N', NCLIN, N, ONE, A, NROWA,
     $                   X, 1, ZERO, AX, 1 )

*        ---------------------------------------------------------------
*        Check the row residuals.
*        ---------------------------------------------------------------
         IF (NACTIV .GT. 0) THEN
            DO 300 K = 1, NACTIV
               I   = KACTIV(K)
               J   = N + I
               IS  = ISTATE(J)
               IF (IS .EQ. 1) WORK(K) = BL(J) - AX(I)
               IF (IS .GE. 2) WORK(K) = BU(J) - AX(I)
  300       CONTINUE

            JMAX   = IDAMAX( NACTIV, WORK, 1 )
            ERRMAX = DABS( WORK(JMAX) )
         END IF

         KTRY = KTRY + 1
*     UNTIL    (ERRMAX .LE. FEATOL(JMAX) .OR. KTRY .GT. NTRY
      IF (.NOT.(ERRMAX .LE. FEATOL(JMAX) .OR. KTRY .GT. NTRY)) GO TO 200

      ROWERR = ERRMAX .GT. FEATOL(JMAX)

*     ==================================================================
*     Compute the linear objective value  c'x  and the transformed
*     residual  Pr  -  RQ'x = RES0  -  RQ'x.
*     ==================================================================
      IF (NRANK .GT. 0  .OR.  LINOBJ) THEN
         CALL DCOPY ( N, X, 1, P, 1 )
         CALL CMQMUL( 6, N, NZ, NFREE, NQ, UNITQ, KX, P, ZY, WORK )
      END IF

      CTX = ZERO
      IF (LINOBJ)
     $   CTX = DDOT  ( N, CQ, 1, P, 1 )

      IF (NRANK .GT. 0) THEN

         CALL DTRMV ( 'U', 'N', 'N', NRANK, R, NROWR, P, 1 )
         IF (NRANK .LT. N)
     $      CALL DGEMV ( 'N', NRANK, N-NRANK, ONE, R(1,NRANK+1), NROWR,
     $                   P(NRANK+1), 1, ONE, P, 1 )

         CALL DCOPY ( NRANK,       RES0, 1, RES, 1 )
         CALL DXMY ( NRANK, P   , 1, RES, 1 )

      END IF

      IF (LSDBG  .AND.  ILSDBG(2) .GT. 0)
     $   WRITE (NOUT, 2200) (X(J), J = 1, N)

      RETURN

 2200 FORMAT(/ ' //LSSETX// Variables after refinement ... '/ (5G12.3))

*     End of  LSSETX.

      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE LSSOL ( fail, MM, N,
     $                   NCLIN, NROWA, NROWR,
     $                   A, BL, BU, CVEC,
     $                   ISTATE, KX, X, R, B,
     $                   INFORM, ITER, OBJ, CLAMDA,
     $                   IW, LENIW, W, LENW, SSQ1 )

*     IMPLICIT           DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION  POINT1, POINT3, POINT8, POINT9, ONE,    ZERO
      INTEGER           LENLS,  LDBG,   MXPARM, NOUT,   ILSDBG, ITMAX1, 
     &                  ITMAX2, LCRASH, LDBGLS, LPROB,  MSGLS,  NN,     
     &                  NNCLIN, NPROB,  IPADLS, LENNAM, NROWT,  NCOLT,  
     &                  NQ,     LOCLS,  MSGLVL, IDBG,   MSGDBG, MM,     
     &                  N,      NCLIN,  NROWA,  NROWR,  INFORM, ITER,   
     &                  M,      MINACT, MINFXD, MAXACT, MXFREE, MAXNZ,  
     &                  NCNLN,  LITOTL, LAX,    LWTOTL, NERROR, LKACTV
      INTEGER           LANORM, LPX,    LRES,   LRES0,  LGQ,    LCQ,
     &                  LRLAM,  LT,     LZY,    LWTINF, LWRK,   LFEATL, 
     &                  IANRMJ, J,      IDBGLS, I,      NCTOTL, NRANK,  
     &                  INFO,   LJ,     JMAX,   JSAVE,  NTEMP,  IDRANK, 
     &                  NACTIV, NARTIF, NFREE,  NGQ,    NRES,   NZ,     
     &                  NACT1,  NREJTD, JINF,   NZ1,    NUMINF
      DOUBLE PRECISION  BIGBND, BIGDX,  BNDLOW, BNDUPP, TOLACT, TOLFEA, 
     &                  RPADLS, EPSPT3, EPSPT5, EPSPT8, ASIZE,  DTMAX,  
     &                  DTMIN,  OBJ,    SSQ1,   EPSMCH, RTEPS,  TOLRNK, 
     &                  EPSPT9, CONDMX, DNRM2,  AMIN,   FEAMAX, FEAMIN, 
     &                  ERRMAX, CTX,    XNORM,  SUMINF
      INTEGER            fail, LENIW, LENW
      INTEGER            ISTATE(N+NCLIN), KX(N)
      INTEGER            IW(LENIW)
      DOUBLE PRECISION   BL(N+NCLIN), BU(N+NCLIN), A(NROWA,*)
      DOUBLE PRECISION   CLAMDA(N+NCLIN), CVEC(*)
      DOUBLE PRECISION   R(NROWR,*), X(N), B(*)
      DOUBLE PRECISION   W(LENW)
      dimension          fail(*)

************************************************************************
*  LSSOL  solves problems of the form
*
*           Minimize               F(x)
*              x
*                                 (  x )
*           subject to    bl  .le.(    ).ge.  bu,
*                                 ( Ax )
*
*  where  '  denotes the transpose of a column vector,  x  denotes the
*  n-vector of parameters and  F(x) is one of the following functions..
*
*  FP =  None                         (find a feasible point).
*  LP =  c'x
*  QP1=        1/2 x'Rx                R  n times n, symmetric pos. def.
*  QP2=  c'x + 1/2 x'Rx                .  .   ..        ..       ..  ..
*  QP3=        1/2 x'R'Rx              R  m times n, upper triangular.
*  QP4=  c'x + 1/2 x'R'Rx              .  .   ..  .   ..      ...
*  LS1=        1/2 (b - Rx)'(b - Rx)   R  m times n, rectangular.
*  LS2=  c'x + 1/2 (b - Rx)'(b - Rx)   .  .   ..  .     ...
*  LS3=        1/2 (b - Rx)'(b - Rx)   R  m times n, upper triangular.
*  LS4=  c'x + 1/2 (b - Rx)'(b - Rx)   .  .   ..  .   ..      ...
*
*  The matrix  R  is entered as the two-dimensional array  R  (of row
*  dimension  NROWR).  If  NROWR = 0,  R  is not accessed.
*
*  The vector  c  is entered in the one-dimensional array  CVEC.
*
*  NCLIN  is the number of general linear constraints (rows of  A).
*  (NCLIN may be zero.)
*
*  The first  N  components of  BL  and   BU  are lower and upper
*  bounds on the variables.  The next  NCLIN  components are
*  lower and upper bounds on the general linear constraints.
*
*  The matrix  A  of coefficients in the general linear constraints
*  is entered as the two-dimensional array  A  (of dimension
*  NROWA by N).  If NCLIN = 0, A is not accessed.
*
*  The vector  x  must contain an initial estimate of the solution,
*  and will contain the computed solution on output.
*
*
*  Complete documentation for  LSSOL  is contained in Report SOL 86-1,
*  Users Guide for LSSOL (Version 1.0), by P.E. Gill, S. J. Hammarling,
*  W. Murray, M.A. Saunders and M.H. Wright, Department of
*  Operations Research, Stanford University, Stanford, California 94305.
*
*  Systems Optimization Laboratory, Stanford University.
*  Version 1.01 Dated  30-June-1986.
*
*  Copyright  1984  Stanford University.
*
*  This material may be reproduced by or for the U.S. Government pursu-
*  ant to the copyright license under DAR clause 7-104.9(a) (1979 Mar).
*
*  This material is based upon work partially supported by the National
*  Science Foundation under Grants MCS-7926009 and ECS-8312142; the
*  Department of Energy Contract AM03-76SF00326, PA No. DE-AT03-
*  76ER72018; the Army Research Office Contract DAA29-84-K-0156;
*  and the Office of Naval Research Grant N00014-75-C-0267.
************************************************************************
      DOUBLE PRECISION   WMACH
      COMMON    /SOLMCH/ WMACH(15)
      SAVE      /SOLMCH/
      COMMON    /SOL1CM/ NOUT
      COMMON    /SOL3CM/ LENNAM, NROWT, NCOLT, NQ
      COMMON    /SOL4CM/ EPSPT3, EPSPT5, EPSPT8, EPSPT9
      COMMON    /SOL5CM/ ASIZE, DTMAX, DTMIN

      PARAMETER         (LENLS = 20)
      COMMON    /SOL1LS/ LOCLS(LENLS)

      LOGICAL            LSDBG
      PARAMETER         (LDBG = 5)
      COMMON    /LSDEBG/ ILSDBG(LDBG), LSDBG
*-----------------------------------------------------------------------
      PARAMETER         (MXPARM = 30)
      INTEGER            IPRMLS(MXPARM), IPSVLS
      DOUBLE PRECISION   RPRMLS(MXPARM), RPSVLS

      COMMON    /LSPAR1/ IPSVLS(MXPARM),
     $                   IDBGLS, ITMAX1, ITMAX2, LCRASH, LDBGLS, LPROB ,
     $                   MSGLS , NN    , NNCLIN, NPROB , IPADLS(20)

      COMMON    /LSPAR2/ RPSVLS(MXPARM),
     $                   BIGBND, BIGDX , BNDLOW, BNDUPP, TOLACT, TOLFEA,
     $                   TOLRNK, RPADLS(23)

      EQUIVALENCE       (IPRMLS(1), IDBGLS), (RPRMLS(1), BIGBND)

      SAVE      /LSPAR1/, /LSPAR2/
*-----------------------------------------------------------------------
      EQUIVALENCE   (MSGLS , MSGLVL), (IDBGLS, IDBG), (LDBGLS, MSGDBG)

      INTRINSIC          MAX, MIN

*     Local variables.

      LOGICAL            COLD  , FACTRZ, LINOBJ, NAMED , ROWERR,
     $                   UNITQ , VERTEX
      CHARACTER*2        PRBTYP
      CHARACTER*8        NAMES(1)
      PARAMETER        ( ZERO   =0.0D+0, POINT1 =0.1D+0, POINT3 =3.3D-1)
      PARAMETER        ( POINT8 =0.8D+0, POINT9 =0.9D+0, ONE    =1.0D+0)

      CHARACTER*40       TITLE
      DATA               TITLE
     $                 / 'SOL/LSSOL  ---  Version 1.01   June 1986' /

*     Set the machine-dependent constants.
      if (fail(1) .gt. 0) return

      CALL MCHPAR()

      EPSMCH = WMACH( 3)
      RTEPS  = WMACH( 4)
      NOUT   = WMACH(11)

      EPSPT3 = EPSMCH**POINT3
      EPSPT5 = RTEPS
      EPSPT8 = EPSMCH**POINT8
      EPSPT9 = EPSMCH**POINT9

      NAMED  = .FALSE.

      INFORM = 0
      ITER   = 0

      CONDMX = ONE / EPSPT5

      NCTOTL = N + NCLIN

*     Set the default values of the parameters.

      CALL LSDFLT( MM, N, NCLIN, TITLE )

*     Set all parameters determined by the problem type.

      IF      (LPROB .EQ. 1 ) THEN
         PRBTYP    = 'FP'
         M      = 0
         LINOBJ = .FALSE.
         FACTRZ = .TRUE.
      ELSE IF (LPROB .EQ. 2 ) THEN
         PRBTYP    = 'LP'
         M      = 0
         LINOBJ = .TRUE.
         FACTRZ = .TRUE.
      ELSE IF (LPROB .EQ. 3 ) THEN
         PRBTYP    = 'QP'
         M      = MM
         LINOBJ = .FALSE.
         FACTRZ = .TRUE.
      ELSE IF (LPROB .EQ. 4 ) THEN
         PRBTYP    = 'QP'
         M      = MM
         LINOBJ = .TRUE.
         FACTRZ = .TRUE.
      ELSE IF (LPROB .EQ. 5 ) THEN
         PRBTYP    = 'QP'
         M      = MM
         LINOBJ = .FALSE.
         FACTRZ = .FALSE.
      ELSE IF (LPROB .EQ. 6 ) THEN
         PRBTYP    = 'QP'
         M      = MM
         LINOBJ = .TRUE.
         FACTRZ = .FALSE.
      ELSE IF (LPROB .EQ. 7 ) THEN
         PRBTYP    = 'LS'
         M      = MM
         LINOBJ = .FALSE.
         FACTRZ = .TRUE.
      ELSE IF (LPROB .EQ. 8 ) THEN
         PRBTYP    = 'LS'
         M      = MM
         LINOBJ = .TRUE.
         FACTRZ = .TRUE.
      ELSE IF (LPROB .EQ. 9 ) THEN
         PRBTYP    = 'LS'
         M      = MM
         LINOBJ = .FALSE.
         FACTRZ = .FALSE.
      ELSE IF (LPROB .EQ. 10) THEN
         PRBTYP    = 'LS'
         M      = MM
         LINOBJ = .TRUE.
         FACTRZ = .FALSE.
      END IF

*     Assign the dimensions of arrays in the parameter list of LSCORE.
*     Economies of storage are possible if the minimum number of active
*     constraints and the minimum number of fixed variables are known in
*     advance.  The expert user should alter MINACT and MINFXD
*     accordingly.
*     If a linear program is being solved and the matrix of general
*     constraints is fat,  i.e.,  NCLIN .LT. N,  a non-zero value is
*     known for MINFXD.  Note that in this case, VERTEX must be
*     set  .TRUE..

      MINACT = 0
      MINFXD = 0

      VERTEX = .FALSE.
      IF (      (PRBTYP .EQ. 'LP'  .OR.  PRBTYP .EQ. 'FP')
     $    .AND.  NCLIN  .LT. N   ) THEN
         MINFXD = N - NCLIN - 1
         VERTEX = .TRUE.
      END IF

      MXFREE = N - MINFXD
      MAXACT = MAX( 1, MIN( N, NCLIN ) )
      MAXNZ  = N - ( MINFXD + MINACT )

      IF (NCLIN .EQ. 0) THEN
         NQ     = 1
         NROWT  = 1
         NCOLT  = 1
         VERTEX = .FALSE.
      ELSE
         NQ     = MAX( 1, MXFREE )
         NROWT  = MAX( MAXNZ, MAXACT )
         NCOLT  = MXFREE
      END IF

      NCNLN  = 0
      LENNAM = 1

*     Allocate certain arrays that are not done in LSLOC.

      LITOTL = 0

      LAX    = 1
      LWTOTL = LAX + NCLIN  - 1

*     Allocate remaining work arrays.

      CALL LSLOC ( LPROB, N, NCLIN, LITOTL, LWTOTL )

      COLD  = LCRASH .EQ. 0

*     Check input parameters and storage limits.

      CALL CMCHK ( NERROR, MSGLVL, COLD, (.NOT.FACTRZ),
     $             LENIW, LENW, LITOTL, LWTOTL,
     $             N, NCLIN, NCNLN,
     $             ISTATE, KX, NAMED, NAMES, LENNAM,
     $             BL, BU, X )

      IF (NERROR .GT. 0) THEN
         INFORM = 6
         GO TO 800
      END IF

      LKACTV = LOCLS( 1)

      LANORM = LOCLS( 2)
      LPX    = LOCLS( 4)
      LRES   = LOCLS( 5)
      LRES0  = LOCLS( 6)
      LGQ    = LOCLS( 8)
      LCQ    = LOCLS( 9)
      LRLAM  = LOCLS(10)
      LT     = LOCLS(11)
      LZY    = LOCLS(12)
      LWTINF = LOCLS(13)
      LWRK   = LOCLS(14)
      LFEATL = LOCLS(15)

      IF (TOLFEA .GT. ZERO)
     $   CALL DLOAD ( N+NCLIN, (TOLFEA), W(LFEATL), 1 )

      IANRMJ = LANORM
      DO 200 J = 1, NCLIN
         W(IANRMJ) = DNRM2 ( N, A(J,1), NROWA )
         IANRMJ    = IANRMJ + 1
  200 CONTINUE
      IF (NCLIN .GT. 0)
     $   CALL DCOND ( NCLIN, W(LANORM), 1, ASIZE, AMIN )

      CALL DCOND ( NCTOTL, W(LFEATL), 1, FEAMAX, FEAMIN )
      CALL DCOPY ( NCTOTL, W(LFEATL), 1, W(LWTINF), 1 )
      CALL DSCAL ( NCTOTL, (ONE/FEAMIN), W(LWTINF), 1 )

*     SSQ1   = ZERO

      IF (FACTRZ) THEN
*        ===============================================================
*        Factorize R using QR or Cholesky.  KX must be initialized.
*        ===============================================================
         DO 210 I = 1, N
            KX(I) = I
  210    CONTINUE

         IF      (PRBTYP .EQ. 'LP'  .OR.  PRBTYP .EQ. 'FP') THEN
            NRANK = 0
         ELSE IF (PRBTYP .EQ. 'QP') THEN
*           ------------------------------------------------------------
*           Compute the Cholesky factorization of R.  The Hessian is
*           M times M and resides in the upper left-hand corner of R.
*           ------------------------------------------------------------
            DO 220 J = M+1, N
               CALL DLOAD ( M, (ZERO), R(1,J), 1 )
  220       CONTINUE

            CALL LSCHOL( NROWR, M, NRANK, TOLRNK, KX, R, INFO, NOUT )

            IF (NRANK .GT. 0)
     $         CALL DLOAD ( NRANK, (ZERO), W(LRES0), 1 )

         ELSE IF (PRBTYP .EQ. 'LS') THEN
*           ------------------------------------------------------------
*           Compute the orthogonal factorization PRQ = ( U ),  where P
*                                                      ( 0 )
*           is an orthogonal matrix and Q is a permutation matrix.
*           Overwrite R with the upper-triangle U.  The orthogonal
*           matrix P is applied to the residual and discarded.  The
*           permutation is stored in the array KX.  Once U has been
*           computed we need only work with vectors of length N within
*           LSCORE.  However, it is necessary to store the sum of
*           squares of the terms  B(NRANK+1),...,B(M),  where B = Pr.
*           ------------------------------------------------------------
            CALL DGEQRP( 'Column iterchanges', M, N, R, NROWR,
     $                   W(LWRK), IW(LKACTV), W(LGQ), INFO )

            LJ  = LKACTV
            DO 230 J = 1, N
               JMAX = IW(LJ)
               IF (JMAX .GT. J) THEN
                  JSAVE    = KX(JMAX)
                  KX(JMAX) = KX(J)
                  KX(J)    = JSAVE
               END IF
               LJ = LJ + 1
  230       CONTINUE

            NTEMP = MIN(N,M)
            CALL DGEAPQ( 'Transpose', 'Separate', M, NTEMP, R, NROWR,
     $                   W(LWRK), 1, B, M, W(LGQ), INFO )

            NRANK = IDRANK( MIN(N, M), R, NROWR+1, TOLRNK )

            IF (M .GT. NRANK) SSQ1 = DNRM2 ( M-NRANK, B(NRANK+1), 1 )

            IF (NRANK .GT. 0)
     $         CALL DCOPY ( NRANK, B, 1, W(LRES0), 1 )
         END IF
      ELSE
*        ===============================================================
*        R is input as an upper-triangular matrix with M rows.
*        ===============================================================
         NRANK = M
         IF (NRANK .GT. 0) THEN
            IF      (PRBTYP .EQ. 'QP') THEN
               CALL DLOAD ( NRANK, (ZERO), W(LRES0), 1 )
            ELSE IF (PRBTYP .EQ. 'LS') THEN
               CALL DCOPY ( NRANK, B, 1, W(LRES0), 1 )
            END IF
         END IF
      END IF

      IF (       MSGLVL .GT. 0     .AND.  NRANK  .LT. N
     $    .AND.  PRBTYP .NE. 'LP'  .AND.  PRBTYP .NE. 'FP')
     $   WRITE (NOUT, 9000) NRANK

*     ------------------------------------------------------------------
*     Find an initial working set.
*     ------------------------------------------------------------------
      CALL LSCRSH( COLD, VERTEX,
     $             NCLIN, NCTOTL, NACTIV, NARTIF,
     $             NFREE, N, NROWA,
     $             ISTATE, IW(LKACTV),
     $             BIGBND, TOLACT,
     $             A, W(LAX), BL, BU, X, W(LGQ), W(LWRK) )

*     ------------------------------------------------------------------
*     Compute the TQ factorization of the constraints while keeping R in
*     upper-triangular form.  Transformations associated with Q are
*     applied to CQ.  Transformations associated with P are applied to
*     RES0.  If some simple bounds are in the working set,  KX is
*     re-ordered so that the free variables come first.
*     ------------------------------------------------------------------
*     First, add the bounds. To save a bit of work, CQ is not loaded
*     until after KX has been re-ordered.

      NGQ   = 0
      NRES  = 0
      IF (NRANK .GT. 0) NRES = 1
      UNITQ = .TRUE.

      CALL LSBNDS( UNITQ,
     $             INFORM, NZ, NFREE, NRANK, NRES, NGQ,
     $             N, NQ, NROWA, NROWR, NROWT,
     $             ISTATE, KX,
     $             CONDMX,
     $             A, R, W(LT), W(LRES0), W(LCQ),
     $             W(LZY), W(LGQ), W(LWRK) )

      IF (LINOBJ) THEN

*        Install the transformed linear term in CQ.
*        CMQMUL applies the permutations in KX to CVEC.

         NGQ = 1
         CALL DCOPY ( N, CVEC, 1, W(LCQ), 1 )
         CALL CMQMUL( 6, N, NZ, NFREE, NQ, UNITQ,
     $                KX, W(LCQ), W(LZY), W(LWRK) )
      END IF

      IF (NACTIV .GT. 0) THEN
         NACT1  = NACTIV
         NACTIV = 0

         CALL LSADDS( UNITQ, VERTEX,
     $                INFORM, 1, NACT1, NACTIV, NARTIF, NZ, NFREE,
     $                NRANK, NREJTD, NRES, NGQ,
     $                N, NQ, NROWA, NROWR, NROWT,
     $                ISTATE, IW(LKACTV), KX,
     $                CONDMX,
     $                A, R, W(LT), W(LRES0), W(LCQ),
     $                W(LZY), W(LGQ), W(LWRK) )
      END IF

*     ------------------------------------------------------------------
*     Move the initial  x  onto the constraints in the working set.
*     Compute the transformed residual vector  Pr = Pb - RQ'x.
*     ------------------------------------------------------------------
      CALL LSSETX( LINOBJ, ROWERR, UNITQ,
     $             NCLIN, NACTIV, NFREE, NRANK, NZ,
     $             N, NCTOTL, NQ, NROWA, NROWR, NROWT,
     $             ISTATE, IW(LKACTV), KX,
     $             JMAX, ERRMAX, CTX, XNORM,
     $             A, W(LAX), BL, BU, W(LCQ), W(LRES), W(LRES0),
     $             W(LFEATL), R, W(LT), X, W(LZY), W(LPX), W(LWRK) )

      JINF = 0

      CALL LSCORE( PRBTYP, NAMED, NAMES, LINOBJ, UNITQ,
     $             INFORM, ITER, JINF, NCLIN, NCTOTL,
     $             NACTIV, NFREE, NRANK, NZ, NZ1,
     $             N, NROWA, NROWR,
     $             ISTATE, IW(LKACTV), KX,
     $             CTX, OBJ, SSQ1,
     $             SUMINF, NUMINF, XNORM,
     $             BL, BU, A, CLAMDA, W(LAX),
     $             W(LFEATL), R, X, IW, W )

      OBJ    = OBJ    + CTX
      IF (PRBTYP .EQ. 'LS'  .AND.  NRANK .GT. 0)
     $   CALL DCOPY ( NRANK, W(LRES), 1, B, 1 )

*     ==================================================================
*     Print messages if required.
*     ==================================================================
  800 IF (MSGLVL .GT.   0) THEN
         IF (INFORM .EQ.   0) THEN
            IF (PRBTYP .EQ. 'FP') THEN
               WRITE (NOUT, 2001)
            ELSE
               WRITE (NOUT, 2002) PRBTYP
            END IF
         END IF
         IF (INFORM .EQ.   1) WRITE (NOUT, 2010) PRBTYP
         IF (INFORM .EQ.   2) WRITE (NOUT, 2020) PRBTYP
         IF (INFORM .EQ.   3) WRITE (NOUT, 2030)
         IF (INFORM .EQ.   4) WRITE (NOUT, 2040)
         IF (INFORM .EQ.   5) WRITE (NOUT, 2050)
         IF (INFORM .EQ.   6) WRITE (NOUT, 2060) NERROR

         IF (INFORM .LT.   6) THEN
            IF      (NUMINF .EQ. 0) THEN
                IF (PRBTYP .NE. 'FP') WRITE (NOUT, 3000) PRBTYP, OBJ
            ELSE IF (INFORM .EQ. 3) THEN
               WRITE (NOUT, 3010) SUMINF
            ELSE
               WRITE (NOUT, 3020) SUMINF
            END IF
            IF (NUMINF .GT. 0) OBJ = SUMINF
         END IF
      END IF

*     Recover the optional parameters set by the user.

      CALL ICOPY ( MXPARM, IPSVLS, 1, IPRMLS, 1 )
      CALL DCOPY ( MXPARM, RPSVLS, 1, RPRMLS, 1 )

      RETURN

 2001 FORMAT(/ ' Exit LSSOL - Feasible point found.     ')
 2002 FORMAT(/ ' Exit LSSOL - Optimal ', A2, ' solution.')
 2010 FORMAT(/ ' Exit LSSOL - Weak ',    A2, ' solution.')
 2020 FORMAT(/ ' Exit LSSOL - ', A2,         ' solution is unbounded.' )
 2030 FORMAT(/ ' Exit LSSOL - Cannot satisfy the linear constraints. ' )
 2040 FORMAT(/ ' Exit LSSOL - Too many iterations.')
 2050 FORMAT(/ ' Exit LSSOL - Too many iterations without changing X.' )
 2060 FORMAT(/ ' Exit LSSOL - ', I10, ' errors found in the input',
     $         ' parameters.  Problem abandoned.'         )
 3000 FORMAT(/ ' Final ', A2, ' objective value =', G16.7 )
 3010 FORMAT(/ ' Minimum sum of infeasibilities =', G16.7 )
 3020 FORMAT(/ ' Final sum of infeasibilities =',   G16.7 )

 9000 FORMAT(/ ' Rank of the objective function data matrix = ', I5 )

*     End of  LSSOL .

      END
