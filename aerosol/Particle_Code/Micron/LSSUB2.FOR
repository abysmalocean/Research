*     File  LSSUBS FORTRAN
*
*     LSADD    LSADDS   LSBNDS   LSCHOL   LSCORE   LSCRSH   LSDEL
*     LSDFLT   LSFEAS   LSGETP   LSGSET   LSKEY    LSLOC
*     LSMOVE   LSMULS   LSOPTN   LSPRT    LSSETX   LSSOL
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE LSADD ( UNITQ,
     $                   INFORM, IFIX, IADD, JADD,
     $                   NACTIV, NZ, NFREE, NRANK, NRES, NGQ,
     $                   N, NROWA, NQ, NROWR, NROWT,
     $                   KX, CONDMX,
     $                   A, R, T, RES, GQ, ZY, WRK1, WRK2 )

*     IMPLICIT           DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION  CONDBD, DNRM2,  DTNEW,  DDIV,   COND,   TDTMAX, 
     &                  TDTMIN, CS,     SN
      LOGICAL            UNITQ
      INTEGER            INFORM, IFIX,   IADD,   JADD,   NACTIV, NZ,
     &                   NFREE,  NRANK,  NRES,   NGQ,    N,      NROWA,
     &                   NQ,     NROWR,  NROWT
      INTEGER            KX(N)
      DOUBLE PRECISION   A(NROWA,*), R(NROWR,*), T(NROWT,*),
     $                   RES(N,*), GQ(N,*), ZY(NQ,*)
      DOUBLE PRECISION   WRK1(N), WRK2(N)
      DOUBLE PRECISION   CONDMX
************************************************************************
*  LSADD   updates the factorization,  A(free) * (Z Y) = (0 T),  when a
*  constraint is added to the working set.  If  NRANK .gt. 0, the
*  factorization  ( R ) = PWQ  is also updated,  where  W  is the
*                 ( 0 )
*  least squares matrix,  R  is upper-triangular,  and  P  is an
*  orthogonal matrix.  The matrices  W  and  P  are not stored.
*
*  There are three separate cases to consider (although each case
*  shares code with another)...
*
*  (1) A free variable becomes fixed on one of its bounds when there
*      are already some general constraints in the working set.
*
*  (2) A free variable becomes fixed on one of its bounds when there
*      are only bound constraints in the working set.
*
*  (3) A general constraint (corresponding to row  IADD  of  A) is
*      added to the working set.
*
*  In cases (1) and (2), we assume that  KX(IFIX) = JADD.
*  In all cases,  JADD  is the index of the constraint being added.
*
*  matrix  Q = (Z Y)  is the identity and will not be touched.
*
*  If  NRES .GT. 0,  the row transformations are applied to the rows of
*  the  (N by NRES)  matrix  RES.
*  If  NGQ .GT. 0,  the column transformations are applied to the
*  columns of the  (NGQ by N)  matrix  GQ'.
*
*  Systems Optimization Laboratory, Stanford University.
*  Original version written 31-October--1984.
*  This version of LSADD dated 29-December-1985.
************************************************************************
      INTEGER            I,      J,      K,      NOUT,   LDBG,   ILSDBG,
     &                   NANEW,  NZERO,  LROWR,  IROWT,  NELM,   NFMIN,
     &                   LCOL
      DOUBLE PRECISION   EPSPT3, EPSPT5, EPSPT8, EPSPT9, ASIZE,  DTMAX,
     &                   DTMIN,  ONE,    ZERO
      COMMON    /SOL1CM/ NOUT
      COMMON    /SOL4CM/ EPSPT3, EPSPT5, EPSPT8, EPSPT9
      COMMON    /SOL5CM/ ASIZE, DTMAX, DTMIN

      LOGICAL            LSDBG
      PARAMETER         (LDBG = 5)
      COMMON    /LSDEBG/ ILSDBG(LDBG), LSDBG

      LOGICAL            BOUND , OVERFL
      LOGICAL            WARNED
      EXTERNAL           DDOT  , DDIV  , DNRM2
      INTRINSIC          DMAX1   , DMIN1
      PARAMETER         (ZERO = 0.0D+0, ONE = 1.0D+0)
      SAVE               WARNED
      DATA               WARNED  /.FALSE./

*     If the condition estimator of the updated factors is greater than
*     CONDBD,  a warning message is printed.

      CONDBD = ONE / EPSPT9

      OVERFL = .FALSE.
      BOUND  = JADD .LE. N
      IF (BOUND) THEN
*        ===============================================================
*        A simple bound has entered the working set.  IADD  is not used.
*        ===============================================================
         IF (LSDBG  .AND.  ILSDBG(1) .GT. 0)
     $      WRITE (NOUT, 1010) NACTIV, NZ, NFREE, IFIX, JADD, UNITQ
         NANEW = NACTIV

         IF (UNITQ) THEN

*           Q  is not stored, but KX defines an ordering of the columns
*           of the identity matrix that implicitly define  Q.
*           Reorder KX so that variable IFIX is moved to position
*           NFREE+1 and variables IFIX+1,...,NFREE+1 are moved one
*           position to the left.

            CALL DLOAD ( NFREE, (ZERO), WRK1, 1 )
            WRK1(IFIX) = ONE

            DO 100 I = IFIX, NFREE-1
               KX(I) = KX(I+1)
  100       CONTINUE
         ELSE
*           ------------------------------------------------------------
*           Q  is stored explicitly.
*           ------------------------------------------------------------
*           Set  WRK1 = the  (IFIX)-th  row of  Q.
*           Move the  (NFREE)-th  row of  Q  to position  IFIX.

            CALL DCOPY ( NFREE, ZY(IFIX,1), NQ, WRK1, 1 )
            IF (IFIX .LT. NFREE) THEN
               CALL DCOPY ( NFREE, ZY(NFREE,1), NQ, ZY(IFIX,1), NQ )
               KX(IFIX) = KX(NFREE)
            END IF
         END IF
         KX(NFREE) = JADD
      ELSE
*        ===============================================================
*        A general constraint has entered the working set.
*        IFIX  is not used.
*        ===============================================================
         IF (LSDBG  .AND.  ILSDBG(1) .GT. 0)
     $      WRITE (NOUT, 1020) NACTIV, NZ, NFREE, IADD, JADD, UNITQ

         NANEW  = NACTIV + 1

*        Transform the incoming row of  A  by  Q'.

         CALL DCOPY ( N, A(IADD,1), NROWA, WRK1, 1 )
         CALL CMQMUL( 8, N, NZ, NFREE, NQ, UNITQ, KX, WRK1, ZY, WRK2)

*        Check that the incoming row is not dependent upon those
*        already in the working set.

         DTNEW  = DNRM2 ( NZ, WRK1, 1 )
         IF (NACTIV .EQ. 0) THEN

*           This is the only general constraint in the working set.

            COND   = DDIV  ( ASIZE, DTNEW, OVERFL )
            TDTMAX = DTNEW
            TDTMIN = DTNEW
         ELSE

*           There are already some general constraints in the working
*           set. Update the estimate of the condition number.

            TDTMAX = DMAX1( DTNEW, DTMAX )
            TDTMIN = DMIN1( DTNEW, DTMIN )
            COND   = DDIV  ( TDTMAX, TDTMIN, OVERFL )
         END IF

         IF (COND .GT. CONDMX  .OR.  OVERFL) GO TO 900

         IF (UNITQ) THEN

*           This is the first general constraint to be added.
*           Set  Q = I.

            DO 200 J = 1, NFREE
               CALL DLOAD ( NFREE, (ZERO), ZY(1,J), 1 )
               ZY(J,J) = ONE
  200       CONTINUE
            UNITQ  = .FALSE.
         END IF
      END IF

      NZERO  = NZ - 1
      IF (BOUND) NZERO = NFREE - 1

*     ------------------------------------------------------------------
*     Use a sequence of 2*2 column transformations to reduce the
*     first NZERO elements of WRK1 to zero.  This affects ZY, except
*     when UNITQ is true.  The transformations may also be applied
*     to R, T and GQ'.
*     ------------------------------------------------------------------
      LROWR  = N
      NELM   = 1
      IROWT  = NACTIV

      DO 300 K = 1, NZERO

*        Compute the transformation that reduces WRK1(K) to zero,
*        then apply it to the relevant columns of  Z  and  GQ'.

         CALL DROT3G( WRK1(K+1), WRK1(K), CS, SN )
         IF (.NOT. UNITQ)
     $      CALL DROT3 ( NFREE, ZY(1,K+1), 1, ZY(1,K), 1, CS, SN )
         IF (NGQ .GT. 0)
     $      CALL DROT3 ( NGQ  , GQ(K+1,1), N, GQ(K,1), N, CS, SN )

         IF (K .GE. NZ  .AND.  NACTIV .GT. 0) THEN

*           Apply the rotation to  T.

            T(IROWT,K) = ZERO
            CALL DROT3 ( NELM, T(IROWT,K+1), 1, T(IROWT,K), 1, CS, SN )
            NELM  = NELM  + 1
            IROWT = IROWT - 1
         END IF

         IF (NRANK .GT. 0) THEN

*           Apply the same transformation to the columns of R.
*           This generates a subdiagonal element in R that must be
*           eliminated by a row rotation.

            IF (K .LT. NRANK) R(K+1,K) = ZERO
            LCOL   = MIN( K+1, NRANK )

            CALL DROT3 ( LCOL, R(1,K+1), 1, R(1,K), 1, CS, SN )
            IF (K .LT. NRANK) THEN
               CALL DROT3G( R(K,K), R(K+1,K), CS, SN )
               LROWR  = LROWR - 1
               CALL DROT3 ( LROWR,   R(K,K+1)  , NROWR,
     $                               R(K+1,K+1), NROWR, CS, SN )

               IF (NRES .GT. 0)
     $            CALL DROT3 ( NRES, RES(K,1)  , N    ,
     $                               RES(K+1,1), N    , CS, SN )
            END IF
         END IF
  300 CONTINUE

      IF (BOUND) THEN

*        The last row and column of ZY has been transformed to plus
*        or minus the unit vector E(NFREE).  We can reconstitute the
*        columns of GQ and R corresponding to the new fixed variable.

         IF (WRK1(NFREE) .LT. ZERO) THEN
            NFMIN = MIN( NRANK, NFREE )
            IF (NFMIN .GT. 0) CALL DSCAL ( NFMIN, -ONE, R(1,NFREE) , 1 )
            IF (NGQ   .GT. 0) CALL DSCAL ( NGQ  , -ONE, GQ(NFREE,1), N )
         END IF

*        ---------------------------------------------------------------
*        The diagonals of T have been altered.  Recompute the
*        largest and smallest values.
*        ---------------------------------------------------------------
         IF (NACTIV .GT. 0) THEN
            CALL DCOND( NACTIV, T(NACTIV,NZ), NROWT-1, TDTMAX, TDTMIN )
            COND   = DDIV  ( TDTMAX, TDTMIN, OVERFL )
         END IF
      ELSE
*        ---------------------------------------------------------------
*        General constraint.  Install the new row of T.
*        ---------------------------------------------------------------
         CALL DCOPY ( NANEW, WRK1(NZ), 1, T(NANEW,NZ), NROWT )
      END IF

*     ==================================================================
*     Prepare to exit.  Check the magnitude of the condition estimator.
*     ==================================================================
  900 IF (NANEW .GT. 0) THEN
         IF (COND .LT. CONDMX  .AND.  .NOT. OVERFL) THEN

*           The factorization has been successfully updated.

            INFORM = 0
            DTMAX  = TDTMAX
            DTMIN  = TDTMIN
CJKW        IF (COND .GE. CONDBD) WRITE (NOUT, 2000) JADD
            IF (COND .GE. CONDBD) THEN
              WARNED = .TRUE.
              WRITE(NOUT, 2000)
            END IF
         ELSE

*           The proposed working set appears to be linearly dependent.

            INFORM = 1
            IF (LSDBG  .AND.  ILSDBG(1) .GT. 0) THEN
               WRITE( NOUT, 3000 )
               IF (BOUND) THEN
                  WRITE (NOUT, 3010) ASIZE, DTMAX, DTMIN
               ELSE
                  IF (NACTIV .GT. 0) THEN
                     WRITE (NOUT, 3020) ASIZE, DTMAX, DTMIN, DTNEW
                  ELSE
                     WRITE (NOUT, 3030) ASIZE, DTNEW
                  END IF
               END IF
            END IF
         END IF
      END IF

      RETURN

 1010 FORMAT(/ ' //LSADD //  Simple bound added.'
     $       / ' //LSADD //  NACTIV    NZ NFREE  IFIX  JADD UNITQ'
     $       / ' //LSADD //  ', 5I6, L6 )
 1020 FORMAT(/ ' //LSADD //  General constraint added.           '
     $       / ' //LSADD //  NACTIV    NZ NFREE  IADD  JADD UNITQ'
     $       / ' //LSADD //  ', 5I6, L6 )
CJKW 2000 FORMAT(/ ' XXX  Serious ill-conditioning in the working set',
CJKW     $         ' after adding constraint ',  I5
CJKW     $       / ' XXX  Overflow may occur in subsequent iterations.'//)
 2000 FORMAT(/ ' Constraint addition has caused ill-conditioning in',
     &         ' the working set.'//)
 3000 FORMAT(/ ' //LSADD //  Dependent constraint rejected.' )
 3010 FORMAT(/ ' //LSADD //     ASIZE     DTMAX     DTMIN        '
     $       / ' //LSADD //', 1P3E10.2 )
 3020 FORMAT(/ ' //LSADD //     ASIZE     DTMAX     DTMIN     DTNEW'
     $       / ' //LSADD //', 1P4E10.2 )
 3030 FORMAT(/ ' //LSADD //     ASIZE     DTNEW'
     $       / ' //LSADD //', 1P2E10.2 )

*     End of  LSADD .

      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE LSADDS( UNITQ, VERTEX,
     $                   INFORM, K1, K2, NACTIV, NARTIF, NZ, NFREE,
     $                   NRANK, NREJTD, NRES, NGQ,
     $                   N, NQ, NROWA, NROWR, NROWT,
     $                   ISTATE, KACTIV, KX,
     $                   CONDMX,
     $                   A, R, T, RES, GQ,
     $                   ZY, WRK1, WRK2 )

*     IMPLICIT           DOUBLE PRECISION(A-H,O-Z)
      INTEGER           INFORM, K1,     K2,     NACTIV, NARTIF, NZ,     
     &                  NFREE,  NRANK,  NREJTD, NRES,   NGQ,    N,      
     &                  NQ,     NROWA,  NROWR,  NROWT,  K,      L,      
     &                  I,      IADD,   IFIX,   JADD,   ISWAP,  NZADD,  
     &                  IARTIF
      DOUBLE PRECISION  ASIZE,  DTMAX,  DTMIN,  RTMAX,  CNDMAX, ONE,    
     &                  ZERO,   ROWMAX, DNRM2,  RNORM
      LOGICAL            UNITQ, VERTEX
      INTEGER            ISTATE(*), KACTIV(N), KX(N)
      DOUBLE PRECISION   CONDMX
      DOUBLE PRECISION   A(NROWA,*), R(NROWR,*),
     $                   T(NROWT,*), RES(N,*), GQ(N,*), ZY(NQ,*)
      DOUBLE PRECISION   WRK1(N), WRK2(N)

************************************************************************
*     LSADDS  includes general constraints K1 thru K2 as new rows of
*     the TQ factorization stored in T, ZY.  If NRANK is nonzero, the
*     changes in Q are reflected in NRANK by N triangular factor R such
*     that
*                         W  =  P ( R ) Q,
*                                 ( 0 )
*     where  P  is orthogonal.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original version written  October-31-1984.
*     This version of LSADDS dated 30-December-1985.
************************************************************************
      DOUBLE PRECISION   WMACH
      COMMON    /SOLMCH/ WMACH(15)
      SAVE      /SOLMCH/
      COMMON    /SOL5CM/ ASIZE, DTMAX, DTMIN

      EXTERNAL           DNRM2
      PARAMETER        ( ZERO = 0.0D+0, ONE = 1.0D+0 )

      RTMAX  = WMACH(8)

*     Estimate the condition number of the constraints that are not
*     to be refactorized.

      IF (NACTIV .EQ. 0) THEN
         DTMAX = ZERO
         DTMIN = ONE
      ELSE
         CALL DCOND ( NACTIV, T(NACTIV,NZ+1), NROWT-1, DTMAX, DTMIN )
      END IF

      DO 200 K = K1, K2
         IADD = KACTIV(K)
         JADD = N + IADD
         IF (NACTIV .LT. NFREE) THEN

            CALL LSADD ( UNITQ,
     $                   INFORM, IFIX, IADD, JADD,
     $                   NACTIV, NZ, NFREE, NRANK, NRES, NGQ,
     $                   N, NROWA, NQ, NROWR, NROWT,
     $                   KX, CONDMX,
     $                   A, R, T, RES, GQ, ZY,
     $                   WRK1, WRK2 )

            IF (INFORM .EQ. 0) THEN
               NACTIV = NACTIV + 1
               NZ     = NZ     - 1
            ELSE
               ISTATE(JADD) =   0
               KACTIV(K)    = - KACTIV(K)
            END IF
         END IF
  200 CONTINUE

      IF (NACTIV .LT. K2) THEN

*        Some of the constraints were classed as dependent and not
*        included in the factorization.  Re-order the part of  KACTIV
*        that holds the indices of the general constraints in the
*        working set.  Move accepted indices to the front and shift
*        rejected indices (with negative values) to the end.

         L      = K1 - 1
         DO 300 K = K1, K2
            I         = KACTIV(K)
            IF (I .GE. 0) THEN
               L      = L + 1
               IF (L .NE. K) THEN
                  ISWAP     = KACTIV(L)
                  KACTIV(L) = I
                  KACTIV(K) = ISWAP
               END IF
            END IF
  300    CONTINUE

*        If a vertex is required, add some temporary bounds.
*        We must accept the resulting condition number of the working
*        set.

         IF (VERTEX) THEN
            CNDMAX = RTMAX
            NZADD  = NZ
            DO 320 IARTIF = 1, NZADD
               ROWMAX = ZERO
               DO 310 I = 1, NFREE
                  RNORM = DNRM2 ( NZ, ZY(I,1), NQ )
                  IF (ROWMAX .LT. RNORM) THEN
                     ROWMAX = RNORM
                     IFIX   = I
                  END IF
  310          CONTINUE
               JADD = KX(IFIX)

               CALL LSADD ( UNITQ,
     $                      INFORM, IFIX, IADD, JADD,
     $                      NACTIV, NZ, NFREE, NRANK, NRES, NGQ,
     $                      N, NROWA, NQ, NROWR, NROWT,
     $                      KX, CNDMAX,
     $                      A, R, T, RES, GQ, ZY,
     $                      WRK1, WRK2 )

               NFREE  = NFREE  - 1
               NZ     = NZ     - 1
               NARTIF = NARTIF + 1
               ISTATE(JADD) = 4
  320       CONTINUE
         END IF
      END IF

      NREJTD = K2 - NACTIV

      RETURN

*     End of  LSADDS.

      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE LSBNDS( UNITQ,
     $                   INFORM, NZ, NFREE, NRANK, NRES, NGQ,
     $                   N, NQ, NROWA, NROWR, NROWT,
     $                   ISTATE, KX,
     $                   CONDMX,
     $                   A, R, T, RES, GQ,
     $                   ZY, WRK1, WRK2 )

*     IMPLICIT           DOUBLE PRECISION(A-H,O-Z)
      INTEGER           INFORM, NZ,     NFREE,  NRANK,  NRES,   NGQ,    
     &                  N,      NQ,     NROWA,  NROWT,  NFIXED, NACTV,  
     &                  J,      JADD,   IFIX,   J2,     L,      NROWR,  
     &                  IADD,   LSTART, K
      LOGICAL            UNITQ
      INTEGER            ISTATE(*), KX(N)
      DOUBLE PRECISION   CONDMX
      DOUBLE PRECISION   A(NROWA,*), R(NROWR,*),
     $                   T(NROWT,*), RES(N,*), GQ(N,*), ZY(NQ,*)
      DOUBLE PRECISION   WRK1(N), WRK2(N)

************************************************************************
*     LSBNDS updates the factor R as KX is reordered to reflect the
*     status of the bound constraints given by ISTATE.  KX is reordered
*     so that the fixed variables come last.  One of two alternative
*     are used to reorder KX. One method needs fewer accesses to KX, the
*     other gives a matrix Rz with more rows and columns.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original version written  30-December-1985.
*     This version dated 30-December-1985.
************************************************************************

      NFIXED = N - NFREE

      IF (NRANK .LT. N  .AND.  NRANK .GT. 0) THEN
*        ---------------------------------------------------------------
*        R is specified but singular.  Try and keep the dimension of Rz
*        as large as possible.
*        ---------------------------------------------------------------
         NACTV = 0
         NFREE = N
         NZ    = N

         J     = N
*+       WHILE (J .GT. 0  .AND.  N-NFREE .LT. NFIXED) DO
  100    IF    (J .GT. 0  .AND.  N-NFREE .LT. NFIXED) THEN
            IF (ISTATE(J) .GT. 0) THEN
               JADD = J
               DO 110 IFIX = NFREE, 1, -1
                  IF (KX(IFIX) .EQ. JADD) GO TO 120
  110          CONTINUE

*              Add bound JADD.

  120          CALL LSADD ( UNITQ,
     $                      INFORM, IFIX, IADD, JADD,
     $                      NACTV, NZ, NFREE, NRANK, NRES, NGQ,
     $                      N, NROWA, NQ, NROWR, NROWT,
     $                      KX, CONDMX,
     $                      A, R, T, RES, GQ, ZY,
     $                      WRK1, WRK2 )

               NFREE = NFREE - 1
               NZ    = NZ    - 1
            END IF
            J = J - 1
            GO TO 100
*+       END WHILE
         END IF
      ELSE
*        ---------------------------------------------------------------
*        R is of full rank,  or is not specified.
*        ---------------------------------------------------------------
         IF (NFIXED .GT. 0) THEN

*           Order KX so that the free variables come first.

            LSTART = NFREE + 1
            DO 250 K = 1, NFREE
               J = KX(K)
               IF (ISTATE(J) .GT. 0) THEN
                  DO 220 L = LSTART, N
                     J2 = KX(L)
                     IF (ISTATE(J2) .EQ. 0) GO TO 230
  220             CONTINUE

  230             KX(K)  = J2
                  KX(L)  = J
                  LSTART = L + 1

                  IF (NRANK .GT. 0)
     $               CALL CMRSWP( N, NRES, NRANK, NROWR, K, L,
     $                            R, RES, WRK1 )
               END IF
  250       CONTINUE

         END IF
         NZ = NFREE
      END IF

      RETURN

*     End of  LSBNDS.

      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE LSCHOL( NROWH, N, NRANK, TOLRNK, KX, H, INFORM, NOUT)

*     IMPLICIT           DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION  ONE,    ZERO,   TOLRNK, DMAX,   D,      SUPMAX
      INTEGER           NROWH,  N,      NRANK,  INFORM, NOUT,   J,      
     &                  IDAMAX, KMAX,   K,      I
      INTEGER            KX(*)
      DOUBLE PRECISION   H(NROWH,*)

************************************************************************
*     LSCHOL  forms the Cholesky factorization of the positive
*     semi-definite matrix H such that
*                   PHP'  =  R'R
*     where  P  is a permutation matrix and  R  is upper triangular.
*     The permutation P is chosen to maximize the diagonal of R at each
*     stage.  Only the diagonal and super-diagonal elements of H are
*     used.
*
*     Output:
*
*         INFORM = 0   the factorization was computed successfully,
*                      with the Cholesky factor written in the upper
*                      triangular part of H and P stored in KX.
*                  1   the matrix H was indefinite.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original version of LSCHOL dated  2-February-1981.
*     Level 2 Blas added 29-June-1986.
*     This version of LSCHOL dated  30-June-1986.
************************************************************************

      DOUBLE PRECISION   TEMP
C     COMMON    /SOL1CM/ NOUT
      INTRINSIC          DABS   , DMAX1   , DSQRT
      EXTERNAL           IDAMAX
      PARAMETER        ( ZERO = 0.0D+0, ONE = 1.0D+0 )

      INFORM = 0
      NRANK  = 0

*     Main loop for computing rows of  R.

      DO 200 J = 1, N

*        Find maximum available diagonal.

         KMAX = J - 1 + IDAMAX( N-J+1, H(J,J), NROWH+1 )
         DMAX = H(KMAX,KMAX)

         IF (DMAX .LE. TOLRNK*DABS(H(1,1))) GO TO 300

*        Perform a symmetric interchange if necessary.

         IF (KMAX .NE. J) THEN
            K        = KX(KMAX)
            KX(KMAX) = KX(J)
            KX(J)    = K

            CALL DSWAP ( J       , H(1,J)   , 1, H(1,KMAX), 1     )
            CALL DSWAP ( KMAX-J,   H(J,KMAX), 1, H(J,J)   , NROWH )
            TEMP         = H(KMAX,KMAX)
            H(KMAX,KMAX) = H(J,KMAX)
            H(J,KMAX)    = TEMP
            CALL DSWAP ( N-KMAX+1, H(KMAX,KMAX), NROWH,
     $                             H(J,KMAX)   , NROWH )

            H(KMAX,KMAX) = H(J,J)
         END IF

*        Set the diagonal of  R.

         D      = DSQRT( DMAX )
         H(J,J) = D
         NRANK  = NRANK + 1

         IF (J .LT. N) THEN

*           Set the super-diagonal elements of this row of R and update
*           the elements of the block that is yet to be factorized.

            CALL DSCAL ( N-J,   (ONE/D), H(J  ,J+1), NROWH )
            CALL DSYRM1( 'U', N-J, H(J  ,J+1), NROWH,
     $                        H(J+1,J+1), NROWH )
         END IF

  200 CONTINUE
*     ------------------------------------------------------------------
*     Check for the semi-definite case.
*     ------------------------------------------------------------------
  300 IF (NRANK .LT. N) THEN

*        Find the largest element in the unfactorized block.

         SUPMAX = ZERO
         DO 310 I = J, N-1
            K      = I + IDAMAX( N-I, H(I,I+1), NROWH )
            SUPMAX = MAX( SUPMAX, DABS(H(I,K)) )
  310    CONTINUE

         IF (SUPMAX .GT. TOLRNK*DABS(H(1,1))) THEN
            WRITE (NOUT, 1000) DMAX, SUPMAX
            INFORM = 1
         END IF
      END IF

      RETURN

 1000 FORMAT(' XXX  Hessian appears to be indefinite.'
     $      /' XXX  Maximum diagonal and off-diagonal ignored',
     $             ' in the Cholesky factorization:', 1P2E22.14 )

*     End of LSCHOL.

      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE LSCORE( PRBTYP, NAMED, NAMES, LINOBJ, UNITQ,
     $                   INFORM, ITER, JINF, NCLIN, NCTOTL,
     $                   NACTIV, NFREE, NRANK, NZ, NZ1,
     $                   N, NROWA, NROWR,
     $                   ISTATE, KACTIV, KX,
     $                   CTX, SSQ, SSQ1, SUMINF, NUMINF, XNORM,
     $                   BL, BU, A, CLAMDA, AX,
     $                   FEATOL, R, X, IW, W )

*     IMPLICIT           DOUBLE PRECISION(A-H,O-Z)
      INTEGER           LENLS,  LDBG,   MXPARM, MREFN,  MSTALL, NOUT,   
     &                  ILSDBG, LENNAM, NROWT,  NCOLT,  NQ,     ICMDBG, 
     &                  IDBGLS, LPROB,  MSGLS,  NN,     NNCLIN, NPROB,  
     &                  IPADLS, ITMAX1, ITMAX2, LCRASH, LDBGLS, MSGLVL, 
     &                  IDBG,   MSGDBG, INFORM, ITER,   JINF,   NCLIN,  
     &                  NCTOTL, NACTIV, NFREE,  NRANK,  NZ,     NZ1
      INTEGER           N,      NROWA,  NROWR,  NUMINF, LANORM, LAP,    
     &                  LPX,    LRES,   LRES0,  LHZ,    LGQ,    LCQ,    
     &                  LRLAM,  LT,     LZY,    LWTINF, LWRK,   NRES,   
     &                  NGQ,    IREFN,  ITMAX,  JADD,   JDEL,   NCNLN, 
     &                  NPHASE, NSTALL, ISDEL,  IFIX,   IADD,   JMAX1,  
     &                  NVIOL,  JTINY,  JSMLST, JBIGST, KSMLST, KDEL
      INTEGER           IS,     MSGSVD, KBIGST
      DOUBLE PRECISION  HALF,   ONE,    ZERO,   EPSPT3, EPSPT5, EPSPT8, 
     &                  EPSPT9, ASIZE,  DTMAX,  DTMIN,  BIGBND, BIGDX,  
     &                  BNDLOW, BNDUPP, TOLACT, TOLFEA, TOLRNK, RPADLS, 
     &                  CTX,    SSQ,    SSQ1,   SUMINF, XNORM,  EPSMCH, 
     &                  FLMAX,  RTMAX,  ALFA,   CONDMX, DRZMAX, DRZMIN, 
     &                  DNRM2,  RESNRM, GZNORM, GFNORM, ABSRZZ, 
     &                  DDIV,   CONDRZ, CONDT,  DINKY,  OBJSIZ, WSSIZE
      DOUBLE PRECISION  CTP,    PNORM,  BIGALF, ALFHIT, PALFA,  ATPHIT, 
     &                  GZ1NRM, CNORM,  ERR1,   ERR2,   TRULAM
      CHARACTER*2        PRBTYP
      CHARACTER*8        NAMES(*)
      INTEGER            ISTATE(NCTOTL), KACTIV(N), KX(N)
      INTEGER            IW(*)
      DOUBLE PRECISION   BL(NCTOTL), BU(NCTOTL), A(NROWA,*),
     $                   CLAMDA(NCTOTL), AX(*),
     $                   FEATOL(NCTOTL), R(NROWR,*), X(N)
      DOUBLE PRECISION   W(*)
      LOGICAL            NAMED, LINOBJ, UNITQ

************************************************************************
*     LSCORE  is a subroutine for linearly constrained linear-least
*     squares.  On entry, it is assumed that an initial working set of
*     linear constraints and bounds is available.
*     The arrays ISTATE, KACTIV and KX will have been set accordingly
*     and the arrays T and ZY will contain the TQ factorization of
*     the matrix whose rows are the gradients of the active linear
*     constraints with the columns corresponding to the active bounds
*     removed.  the TQ factorization of the resulting (NACTIV by NFREE)
*     matrix is  A(free)*Q = (0 T),  where Q is (NFREE by NFREE) and T
*     is reverse-triangular.
*
*     Values of ISTATE(J) for the linear constraints.......
*
*     ISTATE(J)
*     ---------
*          0    constraint J is not in the working set.
*          1    constraint J is in the working set at its lower bound.
*          2    constraint J is in the working set at its upper bound.
*          3    constraint J is in the working set as an equality.
*
*     Constraint J may be violated by as much as FEATOL(J).
*
*     Systems Optimization Laboratory, Stanford University.
*     This version of  LSCORE  dated  1-August-1986.
*
*     Copyright  1984  Stanford University.
*
*  This material may be reproduced by or for the U.S. Government pursu-
*  ant to the copyright license under DAR clause 7-104.9(a) (1979 Mar).
*
*  This material is based upon work partially supported by the National
*  Science Foundation under grants MCS-7926009 and ECS-8012974; the
*  Department of Energy Contract AM03-76SF00326, PA No. DE-AT03-
*  76ER72018; and the Army Research Office Contract DAA29-79-C-0110.
************************************************************************
      DOUBLE PRECISION   WMACH
      COMMON    /SOLMCH/ WMACH(15)
      SAVE      /SOLMCH/
      COMMON    /SOL1CM/ NOUT
      COMMON    /SOL3CM/ LENNAM, NROWT, NCOLT, NQ
      COMMON    /SOL4CM/ EPSPT3, EPSPT5, EPSPT8, EPSPT9
      COMMON    /SOL5CM/ ASIZE, DTMAX, DTMIN

      INTEGER            LOCLS
      PARAMETER         (LENLS = 20)
      COMMON    /SOL1LS/ LOCLS(LENLS)

      LOGICAL            CMDBG, LSDBG
      PARAMETER         (LDBG = 5)
      COMMON    /LSDEBG/ ILSDBG(LDBG), LSDBG
      COMMON    /CMDEBG/ ICMDBG(LDBG), CMDBG
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

      EXTERNAL           DDIV  , DDOT  , DNRM2
      INTRINSIC          DABS   , MAX   , DSQRT
      LOGICAL            CONVRG, CYCLIN, ERROR , FIRSTV, HITCON,
     $                   HITLOW, NEEDFG, OVERFL, PRNT  , PRNT1 , ROWERR
      LOGICAL            SINGLR, STALL , STATPT, UNBNDD, UNCON , UNITGZ,
     $                   WEAK
      PARAMETER        ( ZERO   =0.0D+0, HALF   =0.5D+0, ONE   =1.0D+0 )
      PARAMETER        ( MREFN  =1     , MSTALL =50                    )

*     Specify the machine-dependent parameters.

      EPSMCH = WMACH(3)
      FLMAX  = WMACH(7)
      RTMAX  = WMACH(8)

      LANORM = LOCLS( 2)
      LAP    = LOCLS( 3)
      LPX    = LOCLS( 4)
      LRES   = LOCLS( 5)
      LRES0  = LOCLS( 6)
      LHZ    = LOCLS( 7)
      LGQ    = LOCLS( 8)
      LCQ    = LOCLS( 9)
      LRLAM  = LOCLS(10)
      LT     = LOCLS(11)
      LZY    = LOCLS(12)
      LWTINF = LOCLS(13)
      LWRK   = LOCLS(14)

*     Set up the adresses of the contiguous arrays  ( RES0, RES )
*     and  ( GQ, CQ ).

      NRES   = 0
      IF (NRANK .GT. 0) NRES = 2
      NGQ    = 1
      IF (LINOBJ) NGQ = 2

*     Initialize.

      IREFN  =   0
      ITER   =   0
      ITMAX  =   ITMAX1
      JADD   =   0
      JDEL   =   0
      NCNLN  =   0
      NPHASE =   1
      NSTALL =   0
      NUMINF = - 1
      NZ1    =   0

      ALFA   = ZERO
      CONDMX = FLMAX
      DRZMAX = ONE
      DRZMIN = ONE
      SSQ    = ZERO

      CYCLIN = .FALSE.
      ERROR  = .FALSE.
      FIRSTV = .FALSE.
      PRNT   = .TRUE.
      PRNT1  = .TRUE.
      NEEDFG = .TRUE.
      STALL  = .TRUE.
      UNCON  = .FALSE.
      UNBNDD = .FALSE.

*     If debug output is required,  print nothing until iteration IDBG.

      MSGSVD = MSGLVL
      IF (IDBG .GT. 0  .AND.  IDBG .LE. ITMAX) THEN
         MSGLVL = 0
      END IF

*======================== start of the main loop =======================
*
*      cyclin = false
*      unbndd = false
*      error  = false
*      k      = 0
*
*      repeat
*            repeat
*                  compute Z'g,  print details of this iteration
*                  stat pt = (Z'g .eq. 0)
*                  if (not stat pt) then
*                     error =  k .ge. itmax
*                     if (not error) then
*                        compute p, alfa
*                        error = unbndd  or  cyclin
*                        if (not error) then
*                           k = k + 1
*                           x = x + alfa p
*                           if (feasible) update Z'g
*                           if necessary, add a constraint
*                        end if
*                     end if
*                  end if
*            until  stat pt  or  error
*
*            compute lam1, lam2, smllst
*            optmul =  smllst .gt. 0
*            if ( not (optmul .or. error) ) then
*                  delete an artificial or regular constraint
*            end if
*      until optmul  or  error
*
*=======================================================================

*     REPEAT
*        REPEAT
  100       IF (NEEDFG) THEN
               IF (NRANK .GT. 0) THEN
                  RESNRM = DNRM2 ( NRANK, W(LRES), 1 )
                  SSQ    = HALF*(SSQ1**2 + RESNRM**2 )
               END IF

               IF (NUMINF .NE. 0) THEN

*                 Compute the transformed gradient of either the sum of
*                 of infeasibilities or the objective.  Initialize
*                 SINGLR and UNITGZ.

                  CALL LSGSET( PRBTYP, LINOBJ, SINGLR, UNITGZ, UNITQ,
     $                         N, NCLIN, NFREE,
     $                         NROWA, NQ, NROWR, NRANK, NZ, NZ1,
     $                         ISTATE, KX,
     $                         BIGBND, TOLRNK, NUMINF, SUMINF,
     $                         BL, BU, A, W(LRES), FEATOL,
     $                         W(LGQ), W(LCQ), R, X, W(LWTINF),
     $                         W(LZY), W(LWRK) )

                  IF (PRBTYP .NE. 'FP'  .AND.  NUMINF .EQ. 0
     $                                  .AND.  NPHASE .EQ. 1) THEN
                     ITMAX  = ITER + ITMAX2
                     NPHASE = 2
                  END IF
               END IF
            END IF

            GZNORM = ZERO
            IF (NZ  .GT. 0 ) GZNORM = DNRM2 ( NZ, W(LGQ), 1 )

            IF (NZ1 .EQ. NZ) THEN
               GZ1NRM = GZNORM
            ELSE
               GZ1NRM = ZERO
               IF (NZ1 .GT. 0) GZ1NRM = DNRM2 ( NZ1, W(LGQ), 1 )
            END IF

            GFNORM = GZNORM
            IF (NFREE .GT. 0  .AND.  NACTIV .GT. 0)
     $         GFNORM = DNRM2 ( NFREE, W(LGQ), 1 )

*           ------------------------------------------------------------
*           Print the details of this iteration.
*           ------------------------------------------------------------
*           Define small quantities that reflect the size of X, R and
*           the constraints in the working set.  If feasible,  estimate
*           the rank and condition number of Rz1.
*           Note that NZ1 .LE. NRANK + 1.

            IF (NZ1 .EQ. 0) THEN
               SINGLR = .FALSE.
            ELSE
               IF (NUMINF .GT. 0  .OR.  NZ1 .GT. NRANK) THEN
                  ABSRZZ = ZERO
               ELSE
                  CALL DCOND ( NZ1, R, NROWR+1, DRZMAX, DRZMIN )
                  ABSRZZ = DABS( R(NZ1,NZ1) )
               END IF
               SINGLR = ABSRZZ .LE. DRZMAX*TOLRNK

               IF (LSDBG  .AND.  ILSDBG(1) .GT. 0)
     $            WRITE (NOUT, 9100) SINGLR, ABSRZZ, DRZMAX, DRZMIN

            END IF

            CONDRZ = DDIV  ( DRZMAX, DRZMIN, OVERFL )
            CONDT  = ONE
            IF (NACTIV .GT. 0)
     $         CONDT  = DDIV  ( DTMAX , DTMIN , OVERFL )

            IF (PRNT) THEN
               CALL LSPRT ( PRBTYP, PRNT1, ISDEL, ITER, JADD, JDEL,
     $                      MSGLVL, NACTIV, NFREE, N, NCLIN,
     $                      NRANK, NROWR, NROWT, NZ, NZ1,
     $                      ISTATE,
     $                      ALFA, CONDRZ, CONDT, GFNORM, GZNORM, GZ1NRM,
     $                      NUMINF, SUMINF, CTX, SSQ,
     $                      AX, R, W(LT), X, W(LWRK) )

               JDEL  = 0
               JADD  = 0
               ALFA  = ZERO
            END IF

            IF (NUMINF .GT. 0) THEN
               DINKY  = ZERO
            ELSE
               OBJSIZ = ONE  + DABS( SSQ + CTX )
               WSSIZE = ZERO
               IF (NACTIV .GT. 0) WSSIZE = DTMAX
               DINKY  = EPSPT8 * MAX( WSSIZE, OBJSIZ, GFNORM )
               IF (UNCON) THEN
                  UNITGZ = GZ1NRM .LE. DINKY
               END IF
            END IF

            IF (LSDBG  .AND.  ILSDBG(1) .GT. 0)
     $         WRITE (NOUT, 9000) UNITGZ, IREFN, GZ1NRM, DINKY

*           If the projected gradient  Z'g  is small and Rz is of full
*           rank, X is a minimum on the working set.  An additional
*           refinement step is allowed to take care of an inaccurate
*           value of DINKY.

            STATPT = .NOT. SINGLR  .AND.  GZ1NRM .LE. DINKY
     $                             .OR.   IREFN  .GT. MREFN

            IF (.NOT. STATPT) THEN
*              ---------------------------------------------------------
*              Compute a search direction.
*              ---------------------------------------------------------
               PRNT  = .TRUE.

               ERROR = ITER .GE. ITMAX
               IF (.NOT. ERROR) THEN

                  IREFN = IREFN + 1
                  ITER  = ITER  + 1

                  IF (ITER .EQ. IDBG) THEN
                     LSDBG  = .TRUE.
                     CMDBG  =  LSDBG
                     MSGLVL =  MSGSVD
                  END IF

                  CALL LSGETP( LINOBJ, SINGLR, UNITGZ, UNITQ,
     $                         N, NCLIN, NFREE,
     $                         NROWA, NQ, NROWR, NRANK, NUMINF, NZ1,
     $                         ISTATE, KX, CTP, PNORM,
     $                         A, W(LAP), W(LRES), W(LHZ), W(LPX),
     $                         W(LGQ), W(LCQ), R, W(LZY), W(LWRK) )

*                 ------------------------------------------------------
*                 Find the constraint we bump into along P.
*                 Update X and AX if the step ALFA is nonzero.
*                 ------------------------------------------------------
*                 ALFHIT is initialized to BIGALF.  If it remains
*                 that way after the call to CMALF, it will be
*                 regarded as infinite.

                  BIGALF = DDIV  ( BIGDX, PNORM, OVERFL )

                  CALL CMALF ( FIRSTV, HITLOW,
     $                         ISTATE, INFORM, JADD, N, NROWA,
     $                         NCLIN, NCTOTL, NUMINF,
     $                         ALFHIT, PALFA, ATPHIT,
     $                         BIGALF, BIGBND, PNORM,
     $                         W(LANORM), W(LAP), AX,
     $                         BL, BU, FEATOL, W(LPX), X )

*                 If  Rz1  is nonsingular,  ALFA = 1.0  will be the
*                 step to the least-squares minimizer on the
*                 current subspace. If the unit step does not violate
*                 the nearest constraint by more than FEATOL,  the
*                 constraint is not added to the working set.

                  HITCON = SINGLR  .OR.  PALFA  .LE. ONE
                  UNCON  = .NOT. HITCON

                  IF (HITCON) THEN
                     ALFA = ALFHIT
                  ELSE
                     JADD   = 0
                     ALFA   = ONE
                  END IF

*                 Check for an unbounded solution or negligible step.

                  UNBNDD =  ALFA .GE. BIGALF
                  STALL  = DABS( ALFA*PNORM ) .LE. EPSPT9*XNORM
                  IF (STALL) THEN
                     NSTALL = NSTALL + 1
                     CYCLIN = NSTALL .GT. MSTALL
                  ELSE
                     NSTALL = 0
                  END IF

                  ERROR = UNBNDD  .OR.  CYCLIN
                  IF (.NOT.  ERROR) THEN
*                    ---------------------------------------------------
*                    Set X = X + ALFA*P.  Update AX, GQ, RES and CTX.
*                    ---------------------------------------------------
                     IF (ALFA .NE. ZERO)
     $                  CALL LSMOVE( HITCON, HITLOW, LINOBJ, UNITGZ,
     $                               NCLIN, NRANK, NZ1,
     $                               N, NROWR, JADD, NUMINF,
     $                               ALFA, CTP, CTX, XNORM,
     $                               W(LAP), AX, BL, BU, W(LGQ),
     $                               W(LHZ), W(LPX), W(LRES),
     $                               R, X, W(LWRK) )

                     IF (HITCON) THEN
*                       ------------------------------------------------
*                       Add a constraint to the working set.
*                       Update the TQ factors of the working set.
*                       Use P as temporary work space.
*                       ------------------------------------------------
*                       Update  ISTATE.

                        IF (BL(JADD) .EQ. BU(JADD)) THEN
                           ISTATE(JADD) = 3
                        ELSE IF (HITLOW) THEN
                           ISTATE(JADD) = 1
                        ELSE
                           ISTATE(JADD) = 2
                        END IF
                        IADD = JADD - N
                        IF (JADD .LE. N) THEN

                           DO 510 IFIX = 1, NFREE
                              IF (KX(IFIX) .EQ. JADD) GO TO 520
  510                      CONTINUE
  520                   END IF

                        CALL LSADD ( UNITQ,
     $                               INFORM, IFIX, IADD, JADD,
     $                               NACTIV, NZ, NFREE, NRANK, NRES,NGQ,
     $                               N, NROWA, NQ, NROWR, NROWT,
     $                               KX, CONDMX,
     $                               A, R, W(LT), W(LRES), W(LGQ),
     $                               W(LZY), W(LWRK), W(LRLAM) )

                        NZ1    = NZ1 - 1
                        NZ     = NZ  - 1

                        IF (JADD .LE. N) THEN

*                          A simple bound has been added.

                           NFREE  = NFREE  - 1
                        ELSE

*                          A general constraint has been added.

                           NACTIV = NACTIV + 1
                           KACTIV(NACTIV) = IADD
                        END IF

                        IREFN  = 0

                     END IF

*                    ---------------------------------------------------
*                    Check the feasibility of constraints with non-
*                    negative ISTATE values.  If some violations have
*                    occurred.  Refine the current X and set INFORM so
*                    that feasibility is checked in LSGSET.
*                    ---------------------------------------------------
                     CALL LSFEAS( N, NCLIN, ISTATE,
     $                            BIGBND, CNORM, ERR1, JMAX1, NVIOL,
     $                            AX, BL, BU, FEATOL, X, W(LWRK) )

                     IF (ERR1 .GT. FEATOL(JMAX1)) THEN
                        CALL LSSETX( LINOBJ, ROWERR, UNITQ,
     $                               NCLIN, NACTIV, NFREE, NRANK, NZ,
     $                               N, NCTOTL, NQ, NROWA, NROWR, NROWT,
     $                               ISTATE, KACTIV, KX,
     $                               JMAX1, ERR2, CTX, XNORM,
     $                               A, AX, BL, BU, W(LCQ),
     $                               W(LRES), W(LRES0), FEATOL, R,
     $                               W(LT), X, W(LZY), W(LPX), W(LWRK) )

                        IF (LSDBG  .AND.  ILSDBG(1) .GT. 0)
     $                     WRITE (NOUT, 2100) ERR1, ERR2
                        IF (ROWERR)       WRITE (NOUT, 2200)

                        UNCON  =   .FALSE.
                        IREFN  =   0
                        NUMINF = - 1
                     END IF
                     NEEDFG = ALFA .NE. ZERO
                  END IF
               END IF
            END IF

*        UNTIL      STATPT  .OR.  ERROR
         IF (.NOT. (STATPT  .OR.  ERROR) ) GO TO 100

*        ===============================================================
*        Try and find the index JDEL of a constraint to drop from
*        the working set.
*        ===============================================================
         JDEL   = 0

         IF (NUMINF .EQ. 0  .AND.  PRBTYP .EQ. 'FP') THEN
            IF (N .GT. NZ)
     $         CALL DLOAD ( N-NZ, (ZERO), W(LRLAM), 1 )
            JTINY  = 0
            JSMLST = 0
            JBIGST = 0
         ELSE

            CALL LSMULS( PRBTYP,
     $                   MSGLVL, N, NACTIV, NFREE,
     $                   NROWA, NROWT, NUMINF, NZ, NZ1,
     $                   ISTATE, KACTIV, KX, DINKY,
     $                   JSMLST, KSMLST, JINF, JTINY,
     $                   JBIGST, KBIGST, TRULAM,
     $                   A, W(LANORM), W(LGQ), W(LRLAM),
     $                   W(LT), W(LWTINF) )

         END IF

         IF (.NOT. ERROR) THEN
            IF (     JSMLST .GT. 0) THEN

*              LSMULS found a regular constraint with multiplier less
*              than (-DINKY).

               JDEL   = JSMLST
               KDEL   = KSMLST
               ISDEL  = ISTATE(JDEL)
               ISTATE(JDEL) = 0

            ELSE IF (JSMLST .LT. 0) THEN

               JDEL   = JSMLST

            ELSE IF (NUMINF .GT. 0  .AND.  JBIGST .GT. 0) THEN

*              No feasible point exists for the constraints but the
*              sum of the constraint violations may be reduced by
*              moving off constraints with multipliers greater than 1.

               JDEL   = JBIGST
               KDEL   = KBIGST
               ISDEL  = ISTATE(JDEL)
               IF (TRULAM .LE. ZERO) IS = - 1
               IF (TRULAM .GT. ZERO) IS = - 2
               ISTATE(JDEL) = IS
               FIRSTV = .TRUE.
               NUMINF = NUMINF + 1
            END IF

            IF      (JDEL .NE. 0  .AND.  SINGLR) THEN

*              Cannot delete a constraint when Rz is singular.
*              Probably a weak minimum.

               JDEL = 0
            ELSE IF (JDEL .NE. 0               ) THEN

*              Constraint JDEL has been deleted.
*              Update the matrix factorizations.

               CALL LSDEL ( UNITQ,
     $                      N, NACTIV, NFREE, NRES, NGQ, NZ, NZ1,
     $                      NROWA, NQ, NROWR, NROWT, NRANK,
     $                      JDEL, KDEL, KACTIV, KX,
     $                      A, W(LRES), R, W(LT), W(LGQ),W(LZY),W(LWRK))

            END IF
         END IF

         IREFN  =  0
         CONVRG =  JDEL .EQ. 0
         PRNT   = .FALSE.
         UNCON  = .FALSE.
         NEEDFG = .FALSE.

*     until       convrg  .or.  error
      IF (.NOT.  (CONVRG  .OR.  ERROR)) GO TO 100

*  .........................End of main loop............................

      WEAK = JTINY .GT. 0  .OR.  SINGLR

      IF (ERROR) THEN
         IF (UNBNDD) THEN
            INFORM = 2
            IF (NUMINF .GT. 0) INFORM = 3
         ELSE IF (ITER .GE. ITMAX) THEN
            INFORM = 4
         ELSE IF (CYCLIN) THEN
            INFORM = 5
         END IF
      ELSE IF (CONVRG) THEN
         INFORM = 0
         IF (NUMINF .GT. 0) THEN
            INFORM = 3
         ELSE IF (PRBTYP .NE. 'FP'  .AND.  WEAK) THEN
            INFORM = 1
         END IF
      END IF

*     ------------------------------------------------------------------
*     Set   CLAMDA.  Print the full solution.
*     ------------------------------------------------------------------
      MSGLVL = MSGSVD
      IF (MSGLVL .GT. 0) WRITE (NOUT, 2000) PRBTYP, ITER, INFORM

      CALL CMPRT ( MSGLVL, NFREE, NROWA,
     $             N, NCLIN, NCNLN, NCTOTL, BIGBND,
     $             NAMED, NAMES, LENNAM,
     $             NACTIV, ISTATE, KACTIV, KX,
     $             A, BL, BU, X, CLAMDA, W(LRLAM), X )

      RETURN

 2000 FORMAT(/ ' Exit from ', A2, ' problem after ', I4, ' iterations.',
     $         '  INFORM =', I3 )
 2100 FORMAT(  ' XXX  Iterative refinement.  Maximum errors before and',
     $         ' after refinement are ',  1P2E14.2 )
 2200 FORMAT(  ' XXX  Warning.  Cannot satisfy the constraints to the',
     $         ' accuracy requested.')
 9000 FORMAT(/ ' //LSCORE//  UNITGZ IREFN     GZ1NRM      DINKY'
     $       / ' //LSCORE//  ', L6, I6, 1P2E11.2 )
 9100 FORMAT(/ ' //LSCORE//  SINGLR   DABS(RZZ1)      DRZMAX     DRZMIN'
     $       / ' //LSCORE//  ', L6,     1P3E12.4 )

*     End of  LSCORE.

      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE LSCRSH( COLD, VERTEX,
     $                   NCLIN, NCTOTL, NACTIV, NARTIF,
     $                   NFREE, N, NROWA,
     $                   ISTATE, KACTIV,
     $                   BIGBND, TOLACT,
     $                   A, AX, BL, BU, X, WX, WORK )

*     IMPLICIT           DOUBLE PRECISION(A-H,O-Z)
      INTEGER           LDBG,   NOUT,   ILSDBG, NCLIN,  NCTOTL, NACTIV, 
     &                  NARTIF, NFREE,  N,      NROWA,  J,      NFIXED,  
     &                  IS,     I,      IMIN,   K,      JMIN
      DOUBLE PRECISION  ONE,    ZERO,   TOLACT, FLMAX,  BIGBND, BIGLOW, 
     &                  BIGUPP, B1,     B2,     DDOT,   TOOBIG, RESMIN, 
     &                  RESL,   RESU,   RESIDL, COLSIZ, COLMIN
      LOGICAL            COLD, VERTEX
      INTEGER            ISTATE(NCTOTL), KACTIV(N)
      DOUBLE PRECISION   A(NROWA,*), AX(*), BL(NCTOTL), BU(NCTOTL),
     $                   X(N), WX(N), WORK(N)

************************************************************************
*     LSCRSH  computes the quantities  ISTATE (optionally), KACTIV,
*     NACTIV, NZ and NFREE  associated with the working set at X.
*     The computation depends upon the value of the input parameter
*     COLD,  as follows...
*
*     COLD = TRUE.  An initial working set will be selected. First,
*                   nearly-satisfied or violated bounds are added.
*                   Next,  general linear constraints are added that
*                   have small residuals.
*
*     COLD = FALSE. The quantities KACTIV, NACTIV, NZ and NFREE are
*                   computed from ISTATE,  specified by the user.
*
*     Values of ISTATE(j)....
*
*        - 2         - 1         0           1          2         3
*     a'x lt bl   a'x gt bu   a'x free   a'x = bl   a'x = bu   bl = bu
*
*     Systems Optimization Laboratory, Stanford University.
*     Original version written 31-October-1984.
*     This version of LSCRSH dated 27-December-1985.
************************************************************************
      DOUBLE PRECISION   WMACH
      COMMON    /SOLMCH/ WMACH(15)
      SAVE      /SOLMCH/
      COMMON    /SOL1CM/ NOUT

      LOGICAL            LSDBG
      PARAMETER         (LDBG = 5)
      COMMON    /LSDEBG/ ILSDBG(LDBG), LSDBG

      EXTERNAL           DDOT
      INTRINSIC          DABS, MIN
      PARAMETER        ( ZERO = 0.0D+0, ONE = 1.0D+0 )

      FLMAX  = WMACH(7)
      CALL DCOPY ( N, X, 1, WX, 1 )

      IF (LSDBG) THEN
         IF (ILSDBG(1) .GT. 0)
     $      WRITE (NOUT, 1000) COLD, NCLIN, NCTOTL
         IF (ILSDBG(2) .GT. 0)
     $      WRITE (NOUT, 1100) (WX(J), J = 1, N)
      END IF

      NFIXED = 0
      NACTIV = 0
      NARTIF = 0

*     If a cold start is being made, initialize  ISTATE.
*     If  BL(j) = BU(j),  set  ISTATE(j)=3  for all variables and linear
*     constraints.

      IF (COLD) THEN
         DO 100 J = 1, NCTOTL
            ISTATE(J) = 0
            IF (BL(J) .EQ. BU(J)) ISTATE(J) = 3
  100    CONTINUE
      ELSE
         DO 110 J = 1, NCTOTL
            IF (ISTATE(J) .GT. 3  .OR.  ISTATE(J) .LT. 0) ISTATE(J) = 0
  110    CONTINUE
      END IF

*     Initialize NFIXED, NFREE and KACTIV.
*     Ensure that the number of bounds and general constraints in the
*     working set does not exceed N.

      DO 200 J = 1, NCTOTL
         IF (NFIXED + NACTIV .EQ. N) ISTATE(J) = 0
         IF (ISTATE(J) .GT. 0) THEN
            IF (J .LE. N) THEN
               NFIXED = NFIXED + 1
               IF (ISTATE(J) .EQ. 1) WX(J) = BL(J)
               IF (ISTATE(J) .GE. 2) WX(J) = BU(J)
            ELSE
               NACTIV = NACTIV + 1
               KACTIV(NACTIV) = J - N
            END IF
         END IF
  200 CONTINUE

*     ------------------------------------------------------------------
*     If a cold start is required,  attempt to add as many
*     constraints as possible to the working set.
*     ------------------------------------------------------------------
      IF (COLD) THEN
         BIGLOW = - BIGBND
         BIGUPP =   BIGBND

*        See if any bounds are violated or nearly satisfied.
*        If so,  add these bounds to the working set and set the
*        variables exactly on their bounds.

         J = N
*+       WHILE (J .GE. 1  .AND.  NFIXED + NACTIV .LT. N) DO
  300    IF    (J .GE. 1  .AND.  NFIXED + NACTIV .LT. N) THEN
            IF (ISTATE(J) .EQ. 0) THEN
               B1     = BL(J)
               B2     = BU(J)
               IS     = 0
               IF (B1 .GT. BIGLOW) THEN
                  IF (WX(J) - B1 .LE. (ONE + DABS( B1 ))*TOLACT) IS = 1
               END IF
               IF (B2 .LT. BIGUPP) THEN
                  IF (B2 - WX(J) .LE. (ONE + DABS( B2 ))*TOLACT) IS = 2
               END IF
               IF (IS .GT. 0) THEN
                  ISTATE(J) = IS
                  IF (IS .EQ. 1) WX(J) = B1
                  IF (IS .EQ. 2) WX(J) = B2
                  NFIXED = NFIXED + 1
               END IF
            END IF
            J = J - 1
            GO TO 300
*+       END WHILE
         END IF

*        ---------------------------------------------------------------
*        The following loop finds the linear constraint (if any) with
*        smallest residual less than or equal to TOLACT  and adds it
*        to the working set.  This is repeated until the working set
*        is complete or all the remaining residuals are too large.
*        ---------------------------------------------------------------
*        First, compute the residuals for all the constraints not in the
*        working set.

         IF (NCLIN .GT. 0  .AND.  NACTIV+NFIXED .LT. N) THEN
            DO 410 I = 1, NCLIN
               IF (ISTATE(N+I) .LE. 0)
     $         AX(I) = DDOT  (N, A(I,1), NROWA, WX, 1 )
  410       CONTINUE

            IS     = 1
            TOOBIG = TOLACT + TOLACT

*+          WHILE (IS .GT. 0  .AND.  NFIXED + NACTIV .LT. N) DO
  500       IF    (IS .GT. 0  .AND.  NFIXED + NACTIV .LT. N) THEN
               IS     = 0
               RESMIN = TOLACT

               DO 520 I = 1, NCLIN
                  J      = N + I
                  IF (ISTATE(J) .EQ. 0) THEN
                     B1     = BL(J)
                     B2     = BU(J)
                     RESL   = TOOBIG
                     RESU   = TOOBIG
                     IF (B1 .GT. BIGLOW)
     $                  RESL  = DABS( AX(I) - B1 ) / (ONE + DABS( B1 ))
                     IF (B2 .LT. BIGUPP)
     $                  RESU  = DABS( AX(I) - B2 ) / (ONE + DABS( B2 ))
                     RESIDL   = MIN( RESL, RESU )
                     IF(RESIDL .LT. RESMIN) THEN
                        RESMIN = RESIDL
                        IMIN   = I
                        IS     = 1
                        IF (RESL .GT. RESU) IS = 2
                     END IF
                  END IF
  520          CONTINUE

               IF (IS .GT. 0) THEN
                  NACTIV = NACTIV + 1
                  KACTIV(NACTIV) = IMIN
                  J         = N + IMIN
                  ISTATE(J) = IS
               END IF
               GO TO 500
*+          END WHILE
            END IF
         END IF

*        ---------------------------------------------------------------
*        If required, add temporary bounds to make a vertex.
*        ---------------------------------------------------------------
         IF (VERTEX  .AND.  NACTIV+NFIXED .LT. N) THEN

*           Compute lengths of columns of selected linear constraints
*           (just the ones corresponding to free variables).

            DO 630 J = 1, N
               IF (ISTATE(J) .EQ. 0) THEN
                  COLSIZ = ZERO
                  DO 620 K = 1, NCLIN
                     IF (ISTATE(N+K) .GT. 0)
     $               COLSIZ = COLSIZ + DABS( A(K,J) )
  620             CONTINUE
                  WORK(J) = COLSIZ
               END IF
  630       CONTINUE

*           Find the  NARTIF  smallest such columns.
*           This is an expensive loop.  Later we can replace it by a
*           4-pass process (say), accepting the first col that is within
*           T  of  COLMIN, where  T = 0.0, 0.001, 0.01, 0.1 (say).
*           (This comment written in 1980).

*+          WHILE (NFIXED + NACTIV .LT. N) DO
  640       IF    (NFIXED + NACTIV .LT. N) THEN
               COLMIN = FLMAX
               DO 650 J = 1, N
                  IF (ISTATE(J) .EQ. 0) THEN
                     IF (NCLIN .EQ. 0) GO TO 660
                     COLSIZ = WORK(J)
                     IF (COLMIN .GT. COLSIZ) THEN
                        COLMIN = COLSIZ
                        JMIN   = J
                     END IF
                  END IF
  650          CONTINUE
               J      = JMIN
  660          ISTATE(J) = 4
               NARTIF = NARTIF + 1
               NFIXED = NFIXED + 1
               GO TO 640
*+          END WHILE
            END IF
         END IF
      END IF

      NFREE = N - NFIXED

      IF (LSDBG) THEN
         IF (ILSDBG(1) .GT. 0)
     $       WRITE (NOUT, 1300) NFIXED, NACTIV, NARTIF
         IF (ILSDBG(2) .GT. 0)
     $       WRITE (NOUT, 1200) (WX(J), J = 1, N)
      END IF

      RETURN

 1000 FORMAT(/ ' //LSCRSH// COLD NCLIN NCTOTL'
     $       / ' //LSCRSH// ', L4, I6, I7 )
 1100 FORMAT(/ ' //LSCRSH// Variables before crash... '/ (5G12.3))
 1200 FORMAT(/ ' //LSCRSH// Variables after  crash... '/ (5G12.3))
 1300 FORMAT(/ ' //LSCRSH// Working set selected ...             '
     $       / ' //LSCRSH// NFIXED NACTIV NARTIF      '
     $       / ' //LSCRSH// ', I6, 2I7 )

*     End of  LSCRSH.

      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE LSDEL ( UNITQ,
     $                   N, NACTIV, NFREE, NRES, NGQ, NZ, NZ1,
     $                   NROWA, NQ, NROWR, NROWT, NRANK,
     $                   JDEL, KDEL, KACTIV, KX,
     $                   A, RES, R, T, GQ, ZY, WORK )

*     IMPLICIT           DOUBLE PRECISION(A-H,O-Z)
      INTEGER           LDBG,   NOUT,   ILSDBG, N,      NACTIV, NFREE,  
     &                  NRES,   NGQ,    NZ,     NZ1,    NROWA,  NQ,     
     &                  NROWR,  NROWT,  NRANK,  JDEL,   KDEL,   IR,     
     &                  IBEGIN, KA,     I,      K,      L,      LD,     
     &                  LROWR,  LCOL,   IDAMAX, JART
      DOUBLE PRECISION  ONE,    ZERO,   ASIZE,  DTMAX,  DTMIN,  CS,     
     &                  SN
      LOGICAL            UNITQ
      INTEGER            KACTIV(N), KX(N)
      DOUBLE PRECISION   A(NROWA,*), RES(N,*), R(NROWR,*), T(NROWT,*),
     $                   GQ(N,*), ZY(NQ,*)
      DOUBLE PRECISION   WORK(N)

************************************************************************
*     LSDEL   updates the least-squares factor R and the factorization
*     A(free) (Z Y) = (0 T) when a regular, temporary or artificial
*     constraint is deleted from the working set.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original version written 31-October-1984.
*     This version of LSDEL dated 10-June-1986.
************************************************************************
      COMMON    /SOL1CM/ NOUT
      COMMON    /SOL5CM/ ASIZE, DTMAX, DTMIN

      LOGICAL            LSDBG
      PARAMETER         (LDBG = 5)
      COMMON    /LSDEBG/ ILSDBG(LDBG), LSDBG

      EXTERNAL           IDAMAX
      PARAMETER        ( ZERO = 0.0D+0, ONE = 1.0D+0 )

      IF (JDEL .GT. 0) THEN

*        Regular constraint or temporary bound deleted.

         IF (JDEL .LE. N) THEN

*           Case 1.  A simple bound has been deleted.
*           =======  Columns NFREE+1 and IR of R must be swapped.

            IR     = NZ    + KDEL
            IF (LSDBG  .AND.  ILSDBG(1) .GT. 0)
     $         WRITE (NOUT, 1100) NACTIV, NZ, NFREE, IR, JDEL, UNITQ

            IBEGIN = 1
            NFREE  = NFREE + 1
            IF (NFREE .LT. IR) THEN
               KX(IR)    = KX(NFREE)
               KX(NFREE) = JDEL
               IF (NRANK .GT. 0)
     $            CALL CMRSWP( N, NRES, NRANK, NROWR, NFREE, IR,
     $                         R, RES, WORK )
               CALL DSWAP ( NGQ, GQ(NFREE,1), N, GQ(IR,1), N )
            END IF

            IF (.NOT. UNITQ) THEN

*              Copy the incoming column of  A(free)  into the end of T.

               DO 130 KA = 1, NACTIV
                  I = KACTIV(KA)
                  T(KA,NFREE) = A(I,JDEL)
  130          CONTINUE

*              Expand Q by adding a unit row and column.

               IF (NFREE .GT. 1) THEN
                  CALL DLOAD ( NFREE-1, ZERO, ZY(NFREE,1), NQ )
                  CALL DLOAD ( NFREE-1, ZERO, ZY(1,NFREE), 1  )
               END IF
               ZY(NFREE,NFREE) = ONE
            END IF
         ELSE

*           Case 2.  A general constraint has been deleted.
*           =======

            IF (LSDBG  .AND.  ILSDBG(1) .GT. 0)
     $         WRITE (NOUT, 1200) NACTIV, NZ, NFREE, KDEL, JDEL, UNITQ

            IBEGIN = KDEL
            NACTIV = NACTIV - 1

*           Delete a row of T and move the ones below it up.

            DO 220 I = KDEL, NACTIV
               KACTIV(I) = KACTIV(I+1)
               LD        = NFREE - I
               CALL DCOPY ( I+1, T(I+1,LD), NROWT, T(I,LD), NROWT )
  220       CONTINUE
         END IF

*        ---------------------------------------------------------------
*        Eliminate the super-diagonal elements of  T,
*        using a backward sweep of 2*2 transformations.
*        ---------------------------------------------------------------
         K     = NFREE  - IBEGIN
         L     = NACTIV - IBEGIN
         LROWR = N      - K

         DO 420 I = IBEGIN, NACTIV
            CALL DROT3G( T(I,K+1), T(I,K), CS, SN )

            IF (L .GT. 0)
     $      CALL DROT3 ( L    , T(I+1,K+1), 1, T(I+1,K ), 1, CS, SN )
            CALL DROT3 ( NFREE, ZY(1,K+1) , 1, ZY(1,K  ), 1, CS, SN )
            CALL DROT3 ( NGQ  , GQ(K+1,1) , N, GQ(K,1)  , N, CS, SN )

*           Apply the column transformations to  R.  The non-zero
*           sub-diagonal that is generated must be eliminated by a row
*           rotation.

            IF (K .LT. NRANK) R(K+1,K) = ZERO
            LCOL   = MIN( K+1, NRANK )
            IF (LCOL .GT. 0)
     $         CALL DROT3 ( LCOL, R(1,K+1), 1, R(1,K), 1, CS, SN )

            IF (K .LT. NRANK) THEN
               CALL DROT3G( R(K,K), R(K+1,K), CS, SN )

               CALL DROT3 ( LROWR, R(K,K+1)    , NROWR,
     $                             R(K+1,K+1)  , NROWR, CS, SN )
               CALL DROT3 ( NRES , RES(K,1)    , N    ,
     $                             RES(K+1,1)  , N    , CS, SN )
            END IF
            K     = K     - 1
            L     = L     - 1
            LROWR = LROWR + 1
  420    CONTINUE

         NZ  = NZ  + 1

*        ---------------------------------------------------------------
*        Estimate the condition number of  T.
*        ---------------------------------------------------------------
         IF (NACTIV .EQ. 0) THEN
            DTMAX = ONE
            DTMIN = ONE
         ELSE
            CALL DCOND ( NACTIV, T(NACTIV,NZ+1), NROWT-1, DTMAX, DTMIN )
         END IF

      END IF

      NZ1 = NZ1 + 1

      IF (NZ .GT. NZ1) THEN
         IF (JDEL .GT. 0) THEN
            JART =   NZ1 - 1 + IDAMAX( NZ-NZ1+1, GQ(NZ1,1), 1 )
         ELSE
            JART = - JDEL
         END IF

         IF (LSDBG  .AND.  ILSDBG(1) .GT. 0)
     $      WRITE( NOUT, 1000 ) NZ, NZ1, JART

         IF (JART .GT. NZ1) THEN

*           Swap columns NZ1 and JART of R.

            IF (UNITQ) THEN
               K        = KX(NZ1)
               KX(NZ1)  = KX(JART)
               KX(JART) = K
            ELSE
               CALL DSWAP ( NFREE, ZY(1,NZ1), 1, ZY(1,JART), 1 )
            END IF

            CALL DSWAP ( NGQ, GQ(NZ1,1), N, GQ(JART,1), N )
            IF (NRANK .GT. 0)
     $         CALL CMRSWP( N, NRES, NRANK, NROWR, NZ1, JART,
     $                      R, RES, WORK )
         END IF
      END IF

      RETURN

 1000 FORMAT(/ ' //LSDEL //  Artificial constraint deleted.      '
     $       / ' //LSDEL //      NZ   NZ1   JART                 '
     $       / ' //LSDEL //  ', 3I6 )
 1100 FORMAT(/ ' //LSDEL //  Simple bound deleted.               '
     $       / ' //LSDEL //  NACTIV    NZ NFREE    IR  JDEL UNITQ'
     $       / ' //LSDEL //  ', 5I6, L6 )
 1200 FORMAT(/ ' //LSDEL //  General constraint deleted.         '
     $       / ' //LSDEL //  NACTIV    NZ NFREE  KDEL  JDEL UNITQ'
     $       / ' //LSDEL //  ', 5I6, L6 )

*     End of  LSDEL .

      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE LSDFLT( M, N, NCLIN, TITLE )

*     IMPLICIT           DOUBLE PRECISION(A-H,O-Z)
      INTEGER           LDBG,   IDUMMY, NOUT,   ILSDBG, ICMDBG, 
     &                  ITMAX1, ITMAX2, LCRASH, LDBGLS, LPROB,  MSGLS,  
     &                  NN,     NNCLIN, NPROB,  IPADLS, MSGLVL, IDBG,   
     &                  MSGDBG, M,      N,      NCLIN,  K,      MSG,    
     &                  I,      MXPARM, IDBGLS, LENT,   J,      NSPACE
      DOUBLE PRECISION  TEN,    ZERO,   RDUMMY, GIGANT, WRKTOL, EPSPT3, 
     &                  EPSPT5, EPSPT8, EPSPT9, BIGBND, BIGDX,  BNDLOW, 
     &                  BNDUPP, TOLACT, TOLFEA, TOLRNK, RPADLS, EPSMCH

      CHARACTER*(*)      TITLE

************************************************************************
*  LSDFLT  loads the default values of parameters not set by the user.
*
*  Systems Optimization Laboratory, Stanford University.
*  Original Fortran 77 version written 17-September-1985.
*  This version of LSDFLT dated   9-September-1986.
************************************************************************
      DOUBLE PRECISION   WMACH
      COMMON    /SOLMCH/ WMACH(15)
      SAVE      /SOLMCH/
      COMMON    /SOL1CM/ NOUT
      COMMON    /SOL4CM/ EPSPT3, EPSPT5, EPSPT8, EPSPT9

      LOGICAL            CMDBG, LSDBG
      PARAMETER         (LDBG = 5)
      COMMON    /LSDEBG/ ILSDBG(LDBG), LSDBG
      COMMON    /CMDEBG/ ICMDBG(LDBG), CMDBG

      LOGICAL            NEWOPT
      COMMON    /SOL3LS/ NEWOPT
      SAVE      /SOL3LS/

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

      LOGICAL            CDEFND
      CHARACTER*4        ICRSH(0:2)
      CHARACTER*3        LSTYPE(1:10)
      CHARACTER*16       KEY
      INTRINSIC          LEN    ,  MAX   , MOD
      PARAMETER        ( ZERO   =  0.0D+0, TEN    = 10.0D+0)
      PARAMETER        ( RDUMMY = -11111., IDUMMY = -11111 )
      PARAMETER        ( GIGANT = 1.0D+10*.99999           )
      PARAMETER        ( WRKTOL = 1.0D-2                   )
      DATA               ICRSH(0), ICRSH(1), ICRSH(2)
     $                 /'COLD'   ,'WARM'   ,'HOT '   /
      DATA               LSTYPE(1), LSTYPE(2)
     $                 /' FP'     ,' LP'     /
      DATA               LSTYPE(3), LSTYPE(4), LSTYPE(5), LSTYPE(6)
     $                 /'QP1'     ,'QP2'     ,'QP3'     ,'QP4'     /
      DATA               LSTYPE(7), LSTYPE(8), LSTYPE(9), LSTYPE(10)
     $                 /'LS1'     ,'LS2'     ,'LS3'     ,'LS4'     /

      EPSMCH = WMACH( 3)

*     Make a dummy call to LSKEY to ensure that the defaults are set.

      CALL LSKEY ( NOUT, '*', KEY )
      NEWOPT = .TRUE.

*     Save the optional parameters set by the user.  The values in
*     RPRMLS and IPRMLS may be changed to their default values.

      CALL ICOPY ( MXPARM, IPRMLS, 1, IPSVLS, 1 )
      CALL DCOPY ( MXPARM, RPRMLS, 1, RPSVLS, 1 )

      IF (       LPROB  .LT. 0      )  LPROB   = 4
                                       CDEFND  = LPROB .EQ. 2*(LPROB/2)
      IF (       LCRASH .LT. 0
     $    .OR.   LCRASH .GT. 2      )  LCRASH  = 0
      IF (       ITMAX1 .LT. 0      )  ITMAX1  = MAX(50, 5*(N+NCLIN))
      IF (       ITMAX2 .LT. 0      )  ITMAX2  = MAX(50, 5*(N+NCLIN))
      IF (       MSGLVL .EQ. IDUMMY )  MSGLVL  = 0
      IF (       IDBG   .LT. 0
     $    .OR.   IDBG   .GT. ITMAX1 + ITMAX2
     $                              )  IDBG    = 0
      IF (       MSGDBG .LT. 0      )  MSGDBG  = 0
      IF (       MSGDBG .EQ. 0      )  IDBG    = ITMAX1 + ITMAX2 + 1
      IF (       TOLACT .LT. ZERO   )  TOLACT  = WRKTOL
      IF (       TOLFEA .EQ. RDUMMY
     $    .OR.  (TOLFEA .GE. ZERO
     $    .AND.  TOLFEA .LT. EPSMCH))  TOLFEA  = EPSPT5
      IF (       TOLRNK .LE. ZERO
     $    .AND.  CDEFND             )  TOLRNK  = EPSPT5
      IF (       TOLRNK .LE. ZERO   )  TOLRNK  = TEN*EPSMCH
      IF (       BIGBND .LE. ZERO   )  BIGBND  = WMACH(8)
      IF (       BIGDX  .LE. ZERO   )  BIGDX   = WMACH(8) * .99D+0

      LSDBG = IDBG .EQ. 0
      CMDBG = LSDBG
      K     = 1
      MSG   = MSGDBG
      DO 200 I = 1, LDBG
         ILSDBG(I) = MOD( MSG/K, 10 )
         ICMDBG(I) = ILSDBG(I)
         K = K*10
  200 CONTINUE

      IF (MSGLVL .GT. 0) THEN

*        Print the title.

         LENT = LEN( TITLE )
         IF (LENT .GT. 0) THEN
            NSPACE = (81 - LENT)/2 + 1
            WRITE (NOUT, '(///// (80A1) )')
     $         (' ', J=1, NSPACE), (TITLE(J:J), J=1,LENT)
            WRITE (NOUT, '(80A1 //)')
     $         (' ', J=1, NSPACE), ('='       , J=1,LENT)
         END IF

         WRITE (NOUT, 2000)
         WRITE (NOUT, 2100) LSTYPE(LPROB),
     $                      NCLIN , TOLFEA, ICRSH(LCRASH),
     $                      N     , BIGBND, TOLACT,
     $                      M     , BIGDX , TOLRNK
         WRITE (NOUT, 2200) EPSMCH, ITMAX1, MSGLVL,
     $                              ITMAX2
      END IF

      RETURN

 2000 FORMAT(
     $//' Parameters'
     $/ ' ----------' )
 2100 FORMAT(
     $/ ' Problem type...........', 7X, A3
     $/ ' Linear constraints.....', I10,     6X,
     $  ' Feasibility tolerance..', 1PE10.2, 6X,
     $  1X, A4, ' start.............'
     $/ ' Variables..............', I10,     6X,
     $  ' Infinite bound size....', 1PE10.2, 6X,
     $  ' Crash tolerance........', 1PE10.2
     $/ ' Objective matrix rows..', I10,     6X,
     $  ' Infinite step size.....', 1PE10.2, 6X,
     $  ' Rank tolerance.........', 1PE10.2 )
 2200 FORMAT(
     $/ ' EPS (machine precision)', 1PE10.2, 6X,
     $  ' Feasibility phase itns.', I10, 6X,
     $  ' Print level............', I10
     $/ 40X,
     $  ' Optimality  phase itns.', I10 )

*     End of  LSDFLT.

      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE LSFEAS( N, NCLIN, ISTATE,
     $                   BIGBND, CVNORM, ERRMAX, JMAX, NVIOL,
     $                   AX, BL, BU, FEATOL, X, WORK )

*     IMPLICIT           DOUBLE PRECISION(A-H,O-Z)
      INTEGER           ILSDBG, LDBG,   NOUT,   N,      NCLIN,  JMAX,  
     &                  J,      I,      NVIOL,  IS,     IDAMAX
      DOUBLE PRECISION  BIGBND, CVNORM, ERRMAX, BIGLOW, BIGUPP, ZERO,   
     &                  FEASJ,  RES,    CON,    TOLJ,   DNRM2
      INTEGER            ISTATE(N+NCLIN)
      DOUBLE PRECISION   AX(*), BL(N+NCLIN), BU(N+NCLIN)
      DOUBLE PRECISION   FEATOL(N+NCLIN), X(N)
      DOUBLE PRECISION   WORK(N+NCLIN)

************************************************************************
*  LSFEAS  computes the following...
*  (1)  The number of constraints that are violated by more
*       than  FEATOL  and the 2-norm of the constraint violations.
*
*  Systems Optimization Laboratory, Stanford University.
*  Original version      April    1984.
*  This version of  LSFEAS  dated  17-October-1985.
************************************************************************
      COMMON    /SOL1CM/ NOUT

      LOGICAL            LSDBG
      PARAMETER         (LDBG = 5)
      COMMON    /LSDEBG/ ILSDBG(LDBG), LSDBG

      EXTERNAL           IDAMAX, DNRM2
      INTRINSIC          DABS
      PARAMETER        ( ZERO = 0.0D+0 )

      BIGLOW = - BIGBND
      BIGUPP =   BIGBND

*     ==================================================================
*     Compute NVIOL,  the number of constraints violated by more than
*     FEATOL,  and CVNORM,  the 2-norm of the constraint violations and
*     residuals of the constraints in the working set.
*     ==================================================================
      NVIOL  = 0

      DO 200 J = 1, N+NCLIN
         FEASJ  = FEATOL(J)
         IS     = ISTATE(J)
         RES    = ZERO

         IF (IS .GE. 0  .AND.  IS .LT. 4) THEN
            IF (J .LE. N) THEN
               CON =  X(J)
            ELSE
               I   = J - N
               CON = AX(I)
            END IF

            TOLJ   = FEASJ

*           Check for constraint violations.

            IF (BL(J) .GT. BIGLOW) THEN
               RES    = BL(J) - CON
               IF (RES .GT.   FEASJ ) NVIOL = NVIOL + 1
               IF (RES .GT.    TOLJ ) GO TO 190
            END IF

            IF (BU(J) .LT. BIGUPP) THEN
               RES    = BU(J) - CON
               IF (RES .LT. (-FEASJ)) NVIOL = NVIOL + 1
               IF (RES .LT.  (-TOLJ)) GO TO 190
            END IF

*           This constraint is satisfied,  but count the residual as a
*           violation if the constraint is in the working set.

            IF (IS .LE. 0) RES = ZERO
            IF (IS .EQ. 1) RES = BL(J) - CON
            IF (IS .GE. 2) RES = BU(J) - CON
            IF (DABS( RES ) .GT. FEASJ) NVIOL = NVIOL + 1
         END IF
  190    WORK(J) = RES
  200 CONTINUE

      JMAX   = IDAMAX( N+NCLIN, WORK, 1 )
      ERRMAX = DABS ( WORK(JMAX) )

      IF (LSDBG  .AND.  ILSDBG(1) .GT. 0)
     $   WRITE (NOUT, 1000) ERRMAX, JMAX

      CVNORM  = DNRM2 ( N+NCLIN, WORK, 1 )

      RETURN

 1000 FORMAT(/ ' //LSFEAS//  The maximum violation is ', 1PE14.2,
     $                     ' in constraint', I5 )

*     End of  LSFEAS.

      END
