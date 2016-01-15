*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE DSYRM1( UPLO, N, X, INCX, A, LDA )
C     SUBROUTINE DSYR  ( UPLO, N, ALPHA, X, INCX, A, LDA )
*     .. Scalar Arguments ..
C     DOUBLE PRECISION   ALPHA
      INTEGER            INCX, LDA, N
      CHARACTER*1        UPLO
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * )
*     ..
*
*  Purpose
*  =======
*
*  DSYR   performs the symmetric rank 1 operation
*
*     A := alpha*x*x' + A,
*
*  where alpha is a real scalar, x is an n element vector and A is an
*  n by n symmetric matrix.
*
*  Parameters
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the upper or lower
*           triangular part of the array A is to be referenced as
*           follows:
*
*              UPLO = 'U' or 'u'   Only the upper triangular part of A
*                                  is to be referenced.
*
*              UPLO = 'L' or 'l'   Only the lower triangular part of A
*                                  is to be referenced.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the order of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  X      - DOUBLE PRECISION array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the n
*           element vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
*           Before entry with  UPLO = 'U' or 'u', the leading n by n
*           upper triangular part of the array A must contain the upper
*           triangular part of the symmetric matrix and the strictly
*           lower triangular part of A is not referenced. On exit, the
*           upper triangular part of the array A is overwritten by the
*           upper triangular part of the updated matrix.
*           Before entry with UPLO = 'L' or 'l', the leading n by n
*           lower triangular part of the array A must contain the lower
*           triangular part of the symmetric matrix and the strictly
*           upper triangular part of A is not referenced. On exit, the
*           lower triangular part of the array A is overwritten by the
*           lower triangular part of the updated matrix.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           max( 1, n ).
*           Unchanged on exit.
*
*
*  Level 2 Blas routine.
*
*  -- Written on 20-July-1986.
*     Sven Hammarling, Nag Central Office.
*
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
*     .. Local Scalars ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, INFO, IX, J, JX, KX
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF      ( .NOT.LSAME( UPLO, 'U' ).AND.
     $          .NOT.LSAME( UPLO, 'L' )      ) THEN
         INFO = 1
      ELSE IF ( N.LT.0 ) THEN
         INFO = 2
      ELSE IF ( INCX.EQ.0 ) THEN
         INFO = 5
      ELSE IF ( LDA.LT.MAX(1,N) ) THEN
         INFO = 7
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DSYR  ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( N.EQ.0 ) )
C     IF( ( N.EQ.0 ).OR.( ALPHA.EQ.ZERO ) )
     $   RETURN
*
*     Set the start point in X if the increment is not unity.
*
      IF( INCX.LE.0 )THEN
         KX = 1 - ( N - 1 )*INCX
      ELSE IF( INCX.NE.1 )THEN
         KX = 1
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through the triangular part
*     of A.
*
      IF( LSAME( UPLO, 'U' ) )THEN
*
*        Form  A  when A is stored in upper triangle.
*
         IF( INCX.EQ.1 )THEN
            DO 20, J = 1, N
               IF( X( J ).NE.ZERO )THEN
                  TEMP = -X( J )
C                 TEMP = ALPHA*X( J )
                  DO 10, I = 1, J
                     A( I, J ) = A( I, J ) + X( I )*TEMP
   10             CONTINUE
               END IF
   20       CONTINUE
         ELSE
            JX = KX
            DO 40, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
C                 TEMP = ALPHA*X( JX )
                  TEMP = -X( JX )
                  IX   = KX
                  DO 30, I = 1, J
                     A( I, J ) = A( I, J ) + X( IX )*TEMP
                     IX        = IX        + INCX
   30             CONTINUE
               END IF
               JX = JX + INCX
   40       CONTINUE
         END IF
      ELSE
*
*        Form  A  when A is stored in lower triangle.
*
         IF( INCX.EQ.1 )THEN
            DO 60, J = 1, N
               IF( X( J ).NE.ZERO )THEN
                  TEMP = -X( J )
C                 TEMP = ALPHA*X( J )
                  DO 50, I = J, N
                     A( I, J ) = A( I, J ) + X( I )*TEMP
   50             CONTINUE
               END IF
   60       CONTINUE
         ELSE
            JX = KX
            DO 80, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = -X( JX )
C                 TEMP = ALPHA*X( JX )
                  IX   = JX
                  DO 70, I = J, N
                     A( I, J ) = A( I, J ) + X( IX )*TEMP
                     IX        = IX        + INCX
   70             CONTINUE
               END IF
               JX = JX + INCX
   80       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of DSYRM1.
*
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE DTRMV ( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
*     .. Scalar Arguments ..
      INTEGER            INCX, LDA, N
      CHARACTER*1        DIAG, TRANS, UPLO
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * )
*     ..
*
*  Purpose
*  =======
*
*  DTRMV  performs one of the matrix-vector operations
*
*     x := A*x,   or   x := A'*x,
*
*  where x is n element vector and A is an n by n unit, or non-unit,
*  upper or lower triangular matrix.
*
*  Parameters
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the matrix is an upper or
*           lower triangular matrix as follows:
*
*              UPLO = 'U' or 'u'   A is an upper triangular matrix.
*
*              UPLO = 'L' or 'l'   A is a lower triangular matrix.
*
*           Unchanged on exit.
*
*  TRANS  - CHARACTER*1.
*           On entry, TRANS specifies the operation to be performed as
*           follows:
*
*              TRANS = 'N' or 'n'   x := A*x.
*
*              TRANS = 'T' or 't'   x := A'*x.
*
*              TRANS = 'C' or 'c'   x := A'*x.
*
*           Unchanged on exit.
*
*  DIAG   - CHARACTER*1.
*           On entry, DIAG specifies whether or not A is unit
*           triangular as follows:
*
*              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
*
*              DIAG = 'N' or 'n'   A is not assumed to be unit
*                                  triangular.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the order of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
*           Before entry with  UPLO = 'U' or 'u', the leading n by n
*           upper triangular part of the array A must contain the upper
*           triangular matrix and the strictly lower triangular part of
*           A is not referenced.
*           Before entry with UPLO = 'L' or 'l', the leading n by n
*           lower triangular part of the array A must contain the lower
*           triangular matrix and the strictly upper triangular part of
*           A is not referenced.
*           Note that when  DIAG = 'U' or 'u', the diagonal elements of
*           A are not referenced either, but are assumed to be unity.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           max( 1, n ).
*           Unchanged on exit.
*
*  X      - DOUBLE PRECISION array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the n
*           element vector x. On exit, X is overwritten with the
*           tranformed vector x.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*
*  Level 2 Blas routine.
*
*  -- Written on 20-July-1986.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
*     .. Local Scalars ..
      INTEGER            I, INFO, IX, J, JX, KX
      LOGICAL            NOUNIT
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF      ( .NOT.LSAME( UPLO, 'U' ).AND.
     $          .NOT.LSAME( UPLO, 'L' )      ) THEN
         INFO = 1
      ELSE IF ( .NOT.LSAME( TRANS, 'N' ).AND.
     $          .NOT.LSAME( TRANS, 'T' ).AND.
     $          .NOT.LSAME( TRANS, 'C' )      ) THEN
         INFO = 2
      ELSE IF ( .NOT.LSAME( DIAG, 'U' ).AND.
     $          .NOT.LSAME( DIAG, 'N' )      ) THEN
         INFO = 3
      ELSE IF ( N.LT.0 ) THEN
         INFO = 4
      ELSE IF ( LDA.LT.MAX(1,N) ) THEN
         INFO = 6
      ELSE IF ( INCX.EQ.0 ) THEN
         INFO = 8
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DTRMV ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( N.EQ.0 )
     $   RETURN
*
      NOUNIT = LSAME( DIAG, 'N' )
*
*     Set up the start point in X if the increment is not unity. This
*     will be  ( N - 1 )*INCX  too small for descending loops.
*
      IF( INCX.LE.0 )THEN
         KX = 1 - ( N - 1 )*INCX
      ELSE IF( INCX.NE.1 )THEN
         KX = 1
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
      IF( LSAME( TRANS, 'N' ) )THEN
*
*        Form  x := A*x.
*
         IF( LSAME( UPLO, 'U' ) )THEN
            IF( INCX.EQ.1 )THEN
               DO 20, J = 1, N
                  IF( X( J ).NE.ZERO )THEN
                     DO 10, I = 1, J - 1
                        X( I ) = X( I ) + X( J )*A( I, J )
   10                CONTINUE
                     IF( NOUNIT )
     $                  X( J ) = X( J )*A( J, J )
                  END IF
   20          CONTINUE
            ELSE
               JX = KX
               DO 40, J = 1, N
                  IF( X( JX ).NE.ZERO )THEN
                     IX = KX
                     DO 30, I = 1, J - 1
                        X( IX ) = X( IX ) + X( JX )*A( I, J )
                        IX      = IX      + INCX
   30                CONTINUE
                     IF( NOUNIT )
     $                  X( JX ) = X( JX )*A( J, J )
                  END IF
                  JX = JX + INCX
   40          CONTINUE
            END IF
         ELSE
            IF( INCX.EQ.1 )THEN
               DO 60, J = N, 1, -1
                  IF( X( J ).NE.ZERO )THEN
                     DO 50, I = N, J + 1, -1
                        X( I ) = X( I ) + X( J )*A( I, J )
   50                CONTINUE
                     IF( NOUNIT )
     $                  X( J ) = X( J )*A( J, J )
                  END IF
   60          CONTINUE
            ELSE
               KX = KX + ( N - 1 )*INCX
               JX = KX
               DO 80, J = N, 1, -1
                  IF( X( JX ).NE.ZERO )THEN
                     IX  = KX
                     DO 70, I = N, J + 1, -1
                        X( IX ) = X( IX ) + X( JX )*A( I, J )
                        IX      = IX      - INCX
   70                CONTINUE
                     IF( NOUNIT )
     $                  X( JX ) = X( JX )*A( J, J )
                  END IF
                  JX = JX - INCX
   80          CONTINUE
            END IF
         END IF
      ELSE
*
*        Form  x := A'*x.
*
         IF( LSAME( UPLO, 'U' ) )THEN
            IF( INCX.EQ.1 )THEN
               DO 100, J = N, 1, -1
                  IF( NOUNIT )
     $               X( J ) = X( J )*A( J, J )
                  DO 90, I = J - 1, 1, -1
                     X( J ) = X( J ) + A( I, J )*X( I )
   90             CONTINUE
  100          CONTINUE
            ELSE
               JX = KX + ( N - 1 )*INCX
               DO 120, J = N, 1, -1
                  IX = JX
                  IF( NOUNIT )
     $               X( JX ) = X( JX )*A( J, J )
                  DO 110, I = J - 1, 1, -1
                     IX      = IX      - INCX
                     X( JX ) = X( JX ) + A( I, J )*X( IX )
  110             CONTINUE
                  JX = JX - INCX
  120          CONTINUE
            END IF
         ELSE
            IF( INCX.EQ.1 )THEN
               DO 140, J = 1, N
                  IF( NOUNIT )
     $               X( J ) = X( J )*A( J, J )
                  DO 130, I = J + 1, N
                     X( J ) = X( J ) + A( I, J )*X( I )
  130             CONTINUE
  140          CONTINUE
            ELSE
               JX = KX
               DO 160, J = 1, N
                  IX = JX
                  IF( NOUNIT )
     $               X( JX ) = X( JX )*A( J, J )
                  DO 150, I = J + 1, N
                     IX      = IX      + INCX
                     X( JX ) = X( JX ) + A( I, J )*X( IX )
  150             CONTINUE
                  JX = JX + INCX
  160          CONTINUE
            END IF
         END IF
      END IF
*
      RETURN
*
*     End of DTRMV .
*
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE DTRSV ( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
*     .. Scalar Arguments ..
      INTEGER            INCX, LDA, N
      CHARACTER*1        DIAG, TRANS, UPLO
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * )
*     ..
*
*  Purpose
*  =======
*
*  DTRSV  solves one of the systems of equations
*
*     A*x = b,   or   A'*x = b,
*
*  where b and x are n element vectors and A is an n by n unit, or
*  non-unit, upper or lower triangular matrix.
*
*  No test for singularity or near-singularity is included in this
*  routine. Such tests must be performed before calling this routine.
*
*  Parameters
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the matrix is an upper or
*           lower triangular matrix as follows:
*
*              UPLO = 'U' or 'u'   A is an upper triangular matrix.
*
*              UPLO = 'L' or 'l'   A is a lower triangular matrix.
*
*           Unchanged on exit.
*
*  TRANS  - CHARACTER*1.
*           On entry, TRANS specifies the equations to be solved as
*           follows:
*
*              TRANS = 'N' or 'n'   A*x = b.
*
*              TRANS = 'T' or 't'   A'*x = b.
*
*              TRANS = 'C' or 'c'   A'*x = b.
*
*           Unchanged on exit.
*
*  DIAG   - CHARACTER*1.
*           On entry, DIAG specifies whether or not A is unit
*           triangular as follows:
*
*              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
*
*              DIAG = 'N' or 'n'   A is not assumed to be unit
*                                  triangular.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the order of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
*           Before entry with  UPLO = 'U' or 'u', the leading n by n
*           upper triangular part of the array A must contain the upper
*           triangular matrix and the strictly lower triangular part of
*           A is not referenced.
*           Before entry with UPLO = 'L' or 'l', the leading n by n
*           lower triangular part of the array A must contain the lower
*           triangular matrix and the strictly upper triangular part of
*           A is not referenced.
*           Note that when  DIAG = 'U' or 'u', the diagonal elements of
*           A are not referenced either, but are assumed to be unity.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           max( 1, n ).
*           Unchanged on exit.
*
*  X      - DOUBLE PRECISION array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the n
*           element right-hand side vector b. On exit, X is overwritten
*           with the solution vector x.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*
*  Level 2 Blas routine.
*
*  -- Written on 20-July-1986.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
*     .. Local Scalars ..
      INTEGER            I, INFO, IX, J, JX, KX
      LOGICAL            NOUNIT
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF      ( .NOT.LSAME( UPLO, 'U' ).AND.
     $          .NOT.LSAME( UPLO, 'L' )      ) THEN
         INFO = 1
      ELSE IF ( .NOT.LSAME( TRANS, 'N' ).AND.
     $          .NOT.LSAME( TRANS, 'T' ).AND.
     $          .NOT.LSAME( TRANS, 'C' )      ) THEN
         INFO = 2
      ELSE IF ( .NOT.LSAME( DIAG, 'U' ).AND.
     $          .NOT.LSAME( DIAG, 'N' )      ) THEN
         INFO = 3
      ELSE IF ( N.LT.0 ) THEN
         INFO = 4
      ELSE IF ( LDA.LT.MAX(1,N) ) THEN
         INFO = 6
      ELSE IF ( INCX.EQ.0 ) THEN
         INFO = 8
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DTRSV ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( N.EQ.0 )
     $   RETURN
*
      NOUNIT = LSAME( DIAG, 'N' )
*
*     Set up the start point in X if the increment is not unity. This
*     will be  ( N - 1 )*INCX  too small for descending loops.
*
      IF( INCX.LE.0 )THEN
         KX = 1 - ( N - 1 )*INCX
      ELSE IF( INCX.NE.1 )THEN
         KX = 1
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
      IF( LSAME( TRANS, 'N' ) )THEN
*
*        Form  x := inv( A )*x.
*
         IF( LSAME( UPLO, 'U' ) )THEN
            IF( INCX.EQ.1 )THEN
               DO 20, J = N, 1, -1
                  IF( X( J ).NE.ZERO )THEN
                     IF( NOUNIT )
     $                  X( J ) = X( J )/A( J, J )
                     DO 10, I = J - 1, 1, -1
                        X( I ) = X( I ) - X( J )*A( I, J )
   10                CONTINUE
                  END IF
   20          CONTINUE
            ELSE
               JX = KX + ( N - 1 )*INCX
               DO 40, J = N, 1, -1
                  IF( X( JX ).NE.ZERO )THEN
                     IF( NOUNIT )
     $                  X( JX ) = X( JX )/A( J, J )
                     IX = JX
                     DO 30, I = J - 1, 1, -1
                        IX      = IX      - INCX
                        X( IX ) = X( IX ) - X( JX )*A( I, J )
   30                CONTINUE
                  END IF
                  JX = JX - INCX
   40          CONTINUE
            END IF
         ELSE
            IF( INCX.EQ.1 )THEN
               DO 60, J = 1, N
                  IF( X( J ).NE.ZERO )THEN
                     IF( NOUNIT )
     $                  X( J ) = X( J )/A( J, J )
                     DO 50, I = J + 1, N
                        X( I ) = X( I ) - X( J )*A( I, J )
   50                CONTINUE
                  END IF
   60          CONTINUE
            ELSE
               JX = KX
               DO 80, J = 1, N
                  IF( X( JX ).NE.ZERO )THEN
                     IF( NOUNIT )
     $                  X( JX ) = X( JX )/A( J, J )
                     IX = JX
                     DO 70, I = J + 1, N
                        IX      = IX      + INCX
                        X( IX ) = X( IX ) - X( JX )*A( I, J )
   70                CONTINUE
                  END IF
                  JX = JX + INCX
   80          CONTINUE
            END IF
         END IF
      ELSE
*
*        Form  x := inv( A' )*x.
*
         IF( LSAME( UPLO, 'U' ) )THEN
            IF( INCX.EQ.1 )THEN
               DO 100, J = 1, N
                  DO 90, I = 1, J - 1
                     X( J ) = X( J ) - A( I, J )*X( I )
   90             CONTINUE
                  IF( NOUNIT )
     $               X( J ) = X( J )/A( J, J )
  100          CONTINUE
            ELSE
               JX = KX
               DO 120, J = 1, N
                  IX = KX
                  DO 110, I = 1, J - 1
                     X( JX ) = X( JX ) - A( I, J )*X( IX )
                     IX      = IX      + INCX
  110             CONTINUE
                  IF( NOUNIT )
     $               X( JX ) = X( JX )/A( J, J )
                  JX = JX + INCX
  120          CONTINUE
            END IF
         ELSE
            IF( INCX.EQ.1 )THEN
               DO 140, J = N, 1, -1
                  DO 130, I = N, J + 1, -1
                     X( J ) = X( J ) - A( I, J )*X( I )
  130             CONTINUE
                  IF( NOUNIT )
     $               X( J ) = X( J )/A( J, J )
  140          CONTINUE
            ELSE
               KX = KX + ( N - 1 )*INCX
               JX = KX
               DO 160, J = N, 1, -1
                  IX = KX
                  DO 150, I = N, J + 1, -1
                     X( JX ) = X( JX ) - A( I, J )*X( IX )
                     IX      = IX      - INCX
  150             CONTINUE
                  IF( NOUNIT )
     $               X( JX ) = X( JX )/A( J, J )
                  JX = JX - INCX
  160          CONTINUE
            END IF
         END IF
      END IF
*
      RETURN
*
*     End of DTRSV .
*
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      LOGICAL FUNCTION LSAME ( CA, CB )
*     .. Scalar Arguments ..
      CHARACTER*1              CA, CB
*     ..
*
*  Purpose
*  =======
*
*  LSAME  tests if CA is the same letter as CB regardless of case.
*  CB is assumed to be an upper case letter. LSAME returns .TRUE. if
*  CA is either the same as CB or the equivalent lower case letter.
*
*  N.B. This version of the routine is only correct for ASCII code.
*       Installers must modify the routine for other character-codes.
*
*       For EBCDIC systems the constant IOFF must be changed to -64.
*       For CDC systems using 6-12 bit representations, the system-
*       specific code in comments must be activated.
*
*  Parameters
*  ==========
*
*  CA     - CHARACTER*1
*  CB     - CHARACTER*1
*           On entry, CA and CB specify characters to be compared.
*           Unchanged on exit.
*
*
*  Auxiliary routine for Level 2 Blas.
*
*  -- Written on 20-July-1986
*     Richard Hanson, Sandia National Labs.
*     Jeremy Du Croz, Nag Central Office.
*
*     .. Parameters ..
      INTEGER                IOFF
      PARAMETER            ( IOFF=32 )
*     .. Intrinsic Functions ..
      INTRINSIC              ICHAR
*     .. Executable Statements ..
*
*     Test if the characters are equal
*
      LSAME = CA .EQ. CB
*
*     Now test for equivalence
*
      IF ( .NOT.LSAME ) THEN
         LSAME = ICHAR(CA) - IOFF .EQ. ICHAR(CB)
      END IF
*
      RETURN
*
*  The following comments contain code for CDC systems using 6-12 bit
*  representations.
*
*     .. Parameters ..
*     INTEGER                ICIRFX
*     PARAMETER            ( ICIRFX=62 )
*     .. Scalar Arguments ..
*     CHARACTER*1            CB
*     .. Array Arguments ..
*     CHARACTER*1            CA(*)
*     .. Local Scalars ..
*     INTEGER                IVAL
*     .. Intrinsic Functions ..
*     INTRINSIC              ICHAR, CHAR
*     .. Executable Statements ..
*
*     See if the first character in string CA equals string CB.
*
*     LSAME = CA(1) .EQ. CB .AND. CA(1) .NE. CHAR(ICIRFX)
*
*     IF (LSAME) RETURN
*
*     The characters are not identical. Now check them for equivalence.
*     Look for the 'escape' character, circumflex, followed by the
*     letter.
*
*     IVAL = ICHAR(CA(2))
*     IF (IVAL.GE.ICHAR('A') .AND. IVAL.LE.ICHAR('Z')) THEN
*        LSAME = CA(1) .EQ. CHAR(ICIRFX) .AND. CA(2) .EQ. CB
*     END IF
*
*     RETURN
*
*     End of LSAME.
*
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE XERBLA( SRNAME, INFO )
      INTEGER            INFO
      CHARACTER*6        SRNAME
*
*     This is a special version of XERBLA to be used only as part of
*     the test program for testing error exits from the Level 2 BLAS
*     routines.
*
*     XERBLA  is an error handler for the Level 2 BLAS routines.
*
*     It is called by the Level 2 BLAS routines if an input parameter is
*     invalid.
*
      INTEGER            INFOT, NOUT
      LOGICAL            OK, LERR
      CHARACTER*6        SRNAMT
      COMMON    /INFOC / INFOT, NOUT, OK, LERR
      COMMON    /SRNAMC/ SRNAMT

      LERR = .TRUE.
      IF (INFO.NE.INFOT) THEN
         WRITE (NOUT,99998) INFO, INFOT
         OK = .FALSE.
      ENDIF
      IF (SRNAME.NE.SRNAMT) THEN
         WRITE (NOUT,99997) SRNAME, SRNAMT
         OK = .FALSE.
      END IF
      RETURN
*
99998 FORMAT (' XXXX    XERBLA was called with INFO = ', I6,
     $        ' instead of ',I2,'   XXXX')
99997 FORMAT (' XXXX    XERBLA was called with SRNAME = ', A6,
     $        ' instead of ',A6,'   XXXX')
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE DCOND ( N, X, INCX, AXMAX, AXMIN )

      INTEGER            N, INCX
      DOUBLE PRECISION   AXMAX, AXMIN
      DOUBLE PRECISION   X( (N-1)*INCX+1 )
C
C     DCOND   finds the elements in  x  that are largest and smallest
C     in magnitude.
C
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
      INTEGER            I, IX
      INTRINSIC          DABS, DMAX1, DMIN1

      IF (N .EQ. 0) THEN
         AXMAX = ZERO
         AXMIN = ZERO
      ELSE
         AXMAX = DABS( X(1) )
         AXMIN = AXMAX
         IX    = 1
         DO 100 I = 2, N
            IX    = IX + INCX
            AXMAX = DMAX1( AXMAX, DABS( X(IX) ) )
            AXMIN = DMIN1( AXMIN, DABS( X(IX) ) )
  100    CONTINUE
      END IF

      RETURN

*     End of  DCOND

      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      DOUBLE PRECISION FUNCTION DDIV  ( A, B, FAIL )
      DOUBLE PRECISION                  A, B
      LOGICAL                           FAIL
C
C  DDIV   returns the value div given by
C
C     div = ( a/b                 if a/b does not overflow,
C           (
C           ( 0.0                 if a .eq. 0.0,
C           (
C           ( sign( a/b )*flmax   if a .ne. 0.0 and a/b would overflow,
C
C  where flmax is a large value, via the function name. In addition if
C  a/b would overflow then fail is returned as true, otherwise fail is
C  returned as false.
C
C  Note that when a and b are both zero, fail is returned as true,
C  but div is returned as 0.0. in all other cases of overflow div is
C  such that abs( div ) = flmax.
C
C
C  Nag Fortran 77 O( 1 ) basic linear algebra routine.
C
C  -- Written on 26-October-1982.
C     Sven Hammarling, Nag Central Office.
C
      INTRINSIC           DABS   , SIGN
      LOGICAL             FIRST
      DOUBLE PRECISION    ABSB  , FLMAX , FLMIN
      DOUBLE PRECISION    ONE   ,         ZERO
      PARAMETER         ( ONE   = 1.0D+0, ZERO  = 0.0D+0 )

      DOUBLE PRECISION   WMACH
      COMMON    /SOLMCH/ WMACH(15)
      SAVE      /SOLMCH/

      SAVE                FIRST , FLMIN , FLMAX
      DATA                FIRST / .TRUE. /

      IF( A.EQ.ZERO )THEN
         DDIV   = ZERO
         IF( B.EQ.ZERO )THEN
            FAIL = .TRUE.
         ELSE
            FAIL = .FALSE.
         END IF
         RETURN
      END IF

      IF( FIRST )THEN
         FIRST  = .FALSE.
         FLMIN  = WMACH( 5 )
         FLMAX  = WMACH( 7 )
      END IF

      IF( B.EQ.ZERO )THEN
         DDIV   = SIGN( FLMAX, A )
         FAIL   = .TRUE.
      ELSE
         ABSB   = DABS( B )
         IF( ABSB.GE.ONE )THEN
            FAIL = .FALSE.
            IF( DABS( A ).GE.ABSB*FLMIN )THEN
               DDIV   = A/B
            ELSE
               DDIV   = ZERO
            END IF
         ELSE
            IF( DABS( A ).LE.ABSB*FLMAX )THEN
               FAIL   = .FALSE.
               DDIV   = A/B
            ELSE
               FAIL   = .TRUE.
               DDIV   = FLMAX
               IF( ( ( A.LT.ZERO ).AND.( B.GT.ZERO ) ).OR.
     $             ( ( A.GT.ZERO ).AND.( B.LT.ZERO ) )     )
     $            DDIV   = -DDIV
            END IF
         END IF
      END IF

      RETURN

*     End of DDIV  .

      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE DDDIV ( N, D, INCD, X, INCX )
      INTEGER            N, INCD, INCX
      DOUBLE PRECISION   D( * ), X( * )
C
C     DDDIV  performs the operation
C
C     x := diag( d )(inverse)*x
C
      EXTERNAL           DSCAL
      INTEGER            I     , ID    , IX
      DOUBLE PRECISION   ONE
      PARAMETER        ( ONE = 1.0 )

      IF( N.GE.1 )THEN
         IF( INCD.EQ.0 )THEN

            CALL DSCAL ( N, (ONE/D( 1 )), X, INCX )

         ELSE IF( ( INCD.EQ.INCX ).AND.( INCD.GT.0 ) )THEN
            DO 10, ID = 1, 1 + ( N - 1 )*INCD, INCD
               X( ID ) = X( ID )/D( ID )
   10       CONTINUE
         ELSE
            IF( INCX.GE.0 )THEN
               IX = 1
            ELSE
               IX = 1 - ( N - 1 )*INCX
            END IF
            IF( INCD.GT.0 )THEN
               DO 20, ID = 1, 1 + ( N - 1 )*INCD, INCD
                  X( IX ) = X( IX )/D( ID )
                  IX      = IX + INCX
   20          CONTINUE
            ELSE
               ID = 1 - ( N - 1 )*INCD
               DO 30, I = 1, N
                  X( IX ) = X( IX )/D( ID )
                  ID      = ID + INCD
                  IX      = IX + INCX
   30          CONTINUE
            END IF
         END IF
      END IF

      RETURN

*     End of DDDIV .

      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE DDSCL ( N, D, INCD, X, INCX )
      INTEGER            N, INCD, INCX
      DOUBLE PRECISION   D( * ), X( * )
C
C  DDSCL  performs the operation
C
C     x := diag( d )*x
C
C
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 22-September-1983.
C     Sven Hammarling, Nag Central Office.
C
      EXTERNAL           DSCAL
      INTEGER            I     , ID    , IX

      IF( N.GE.1 )THEN
         IF( INCD.EQ.0 )THEN

            CALL DSCAL ( N, D( 1 ), X, INCX )

         ELSE IF( ( INCD.EQ.INCX ).AND.( INCD.GT.0 ) )THEN
            DO 10, ID = 1, 1 + ( N - 1 )*INCD, INCD
               X( ID ) = D( ID )*X( ID )
   10       CONTINUE
         ELSE
            IF( INCX.GE.0 )THEN
               IX = 1
            ELSE
               IX = 1 - ( N - 1 )*INCX
            END IF
            IF( INCD.GT.0 )THEN
               DO 20, ID = 1, 1 + ( N - 1 )*INCD, INCD
                  X( IX ) = D( ID )*X( IX )
                  IX      = IX + INCX
   20          CONTINUE
            ELSE
               ID = 1 - ( N - 1 )*INCD
               DO 30, I = 1, N
                  X( IX ) = D( ID )*X( IX )
                  ID      = ID + INCD
                  IX      = IX + INCX
   30          CONTINUE
            END IF
         END IF
      END IF

      RETURN

*     End of DDSCL .

      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE DGRFG ( N, ALPHA, X, INCX, TOL, ZETA )
      INTEGER            N, INCX
      DOUBLE PRECISION   ALPHA, X( * ), TOL, ZETA
C
C  DGRFG  generates details of a generalized Householder reflection such
C  that
C
C     P*( alpha ) = ( beta ),   P'*P = I.
C       (   x   )   (   0  )
C
C  P is given in the form
C
C     P = I - ( zeta )*( zeta  z' ),
C             (   z  )
C
C  where z is an n element vector and zeta is a scalar that satisfies
C
C     1.0 .le. zeta .le. sqrt( 2.0 ).
C
C  zeta is returned in ZETA unless x is such that
C
C     max( abs( x( i ) ) ) .le. max( eps*abs( alpha ), tol )
C
C  where eps is the relative machine precision and tol is the user
C  supplied value TOL, in which case ZETA is returned as 0.0 and P can
C  be taken to be the unit matrix.
C
C  beta is overwritten on alpha and z is overwritten on x.
C  the routine may be called with  n = 0  and advantage is taken of the
C  case where  n = 1.
C
C
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 30-August-1984.
C     Sven Hammarling, Nag Central Office.
C     This version dated 28-September-1984.
C
      EXTERNAL           DSSQ  , DSCAL
      INTRINSIC          DABS   , DMAX1   , SIGN  , DSQRT
      DOUBLE PRECISION   BETA  , EPS   , SCALE , SSQ
      DOUBLE PRECISION   ONE   ,         ZERO
      PARAMETER        ( ONE   = 1.0D+0, ZERO  = 0.0D+0 )

      DOUBLE PRECISION   WMACH
      COMMON    /SOLMCH/ WMACH(15)
      SAVE      /SOLMCH/

      IF( N.LT.1 )THEN
         ZETA = ZERO
      ELSE IF( ( N.EQ.1 ).AND.( X( 1 ).EQ.ZERO ) )THEN
         ZETA = ZERO
      ELSE

         EPS    =  WMACH( 3 )

*        Treat case where P is a 2 by 2 matrix specially.

         IF( N.EQ.1 )THEN

*           Deal with cases where  ALPHA = zero  and
*           abs( X( 1 ) ) .le. max( EPS*abs( ALPHA ), TOL )  first.

            IF( ALPHA.EQ.ZERO )THEN
               ZETA   =  ONE
               ALPHA  =  DABS( X( 1 ) )
               X( 1 ) = -SIGN( ONE, X( 1 ) )
            ELSE IF( DABS( X( 1 ) ).LE.DMAX1( EPS*DABS( ALPHA ),
     $                                     TOL ) )THEN
               ZETA   =  ZERO
            ELSE
               IF( DABS( ALPHA ).GE.DABS( X( 1 ) ) )THEN
                  BETA = DABS ( ALPHA  )*
     $                   DSQRT( ONE + ( X( 1 )/ALPHA )**2 )
               ELSE
                  BETA = DABS ( X( 1 ) )*
     $                   DSQRT( ONE + ( ALPHA/X( 1 ) )**2 )
               END IF
               ZETA   =  DSQRT( ( DABS( ALPHA ) + BETA )/BETA )
               IF( ALPHA.GE.ZERO )BETA = -BETA
               X( 1 ) = -X( 1 )/( ZETA*BETA )
               ALPHA  =  BETA
            END IF
         ELSE

*           Now P is larger than 2 by 2.

            SSQ   = ONE
            SCALE = ZERO

            CALL DSSQ  ( N, X, INCX, SCALE, SSQ )

*           Treat cases where  SCALE = zero,
*           SCALE .le. max( EPS*abs( ALPHA ), TOL )  and
*           ALPHA = zero  specially.
*           Note that  SCALE = max( abs( X( i ) ) ).

            IF( ( SCALE.EQ.ZERO ).OR.
     $          ( SCALE.LE.DMAX1( EPS*DABS( ALPHA ), TOL ) ) )THEN
               ZETA  = ZERO
            ELSE IF( ALPHA.EQ.ZERO )THEN
               ZETA  = ONE
               ALPHA = SCALE*DSQRT( SSQ )

               CALL DSCAL ( N, -ONE/ALPHA, X, INCX )

            ELSE
               IF( SCALE.LT.DABS( ALPHA ) )THEN
                  BETA = DABS ( ALPHA )*
     $                   DSQRT( ONE + SSQ*( SCALE/ALPHA )**2 )
               ELSE
                  BETA = SCALE*
     $                   DSQRT( SSQ +     ( ALPHA/SCALE )**2 )
               END IF
               ZETA = DSQRT( ( BETA + DABS( ALPHA ) )/BETA )
               IF( ALPHA.GT.ZERO )BETA = -BETA

               CALL DSCAL( N, -ONE/( ZETA*BETA ), X, INCX )

               ALPHA = BETA
            END IF
         END IF
      END IF
      RETURN

*     End of DGRFG .

      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE DLOAD ( N, CONST, X, INCX )
      INTEGER            N, INCX
      DOUBLE PRECISION   CONST
      DOUBLE PRECISION   X( * )
C
C  DLOAD  performs the operation
C
C     x = const*e,   e' = ( 1  1 ... 1 ).
C
C
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 22-September-1983.
C     Sven Hammarling, Nag Central Office.
C
      INTEGER            IX
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )

      IF( N.LT.1 )RETURN

      IF( CONST.NE.ZERO )THEN
         DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
            X( IX ) = CONST
   10    CONTINUE
      ELSE
         DO 20, IX = 1, 1 + ( N - 1 )*INCX, INCX
            X( IX ) = ZERO
   20    CONTINUE
      END IF

      RETURN

*     End of DLOAD .

      END
