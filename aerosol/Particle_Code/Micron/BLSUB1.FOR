*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*     File  BLAS FORTRAN
*
*                         Level 1  BLAS
*                         -------  ----
*     DAXPY    DCOPY    DDOT     DNRM2    DSCAL    DSWAP    IDAMAX
*
*                         Level 2  BLAS
*                         -------  ----
*     DGEMV    DGER     DSYMV +  DSYR     DTRMV    DTRSV    LSAME
*     XERBLA
*                         Others
*                         ------
*     DCOND*   DDDIV*   DDIV     DDSCL    DGRFG    DLOAD +  DNORM
*     DROT3*   DROT3G*  DSSQ     ICOPY*   ILOAD    IDRANK+
*
*    *Not in the Nag Blas.
*    +Differs from the Nag Blas.
*
*                         QR Routines
*                         -- --------
*     DGEQR    DGEQRP   DGEAP    DGEAPQ
*
*                         ADDED BY KEN
*                         ----- -- ---
*     IDAMIN   DXPY     DXMY     DGERM1   DSYRM1
*
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE DAXPY ( N, ALPHA, X, INCX, Y, INCY )
      INTEGER            N, INCX, INCY
      DOUBLE PRECISION   ALPHA
      DOUBLE PRECISION   X( * ), Y( * )

C  DAXPY  performs the operation
C
C     y := alpha*x + y
C
C
C  Nag Fortran 77 version of the Blas routine DAXPY .
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 3-September-1982.
C     Sven Hammarling, Nag Central Office.

      INTEGER            I     , IX    , IY
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO  = 0.0D+0 )

      IF( N    .LT.1    )RETURN
      IF( ALPHA.EQ.ZERO )RETURN

      IF( ( INCX.EQ.INCY ).AND.( INCX.GT.0 ) )THEN
         DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
            Y( IX ) = ALPHA*X( IX ) + Y( IX )
   10    CONTINUE
      ELSE
         IF( INCY.GE.0 )THEN
            IY = 1
         ELSE
            IY = 1 - ( N - 1 )*INCY
         END IF
         IF( INCX.GT.0 )THEN
            DO 20, IX = 1, 1 + ( N - 1 )*INCX, INCX
               Y( IY ) = ALPHA*X( IX ) + Y( IY )
               IY      = IY + INCY
   20       CONTINUE
         ELSE
            IX = 1 - ( N - 1 )*INCX
            DO 30, I = 1, N
               Y( IY ) = ALPHA*X( IX ) + Y( IY )
               IX      = IX + INCX
               IY      = IY + INCY
   30       CONTINUE
         END IF
      END IF
      RETURN

*     End of DAXPY .

      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE DXPY  ( N, X, INCX, Y, INCY )
      INTEGER            N, INCX, INCY
      DOUBLE PRECISION   X( * ), Y( * )

C  DXPY  performs the operation
C
C     y := x + y
C
C
C  Nag Fortran 77 version of the Blas routine DAXPY .
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 3-September-1982.
C     Sven Hammarling, Nag Central Office.

      INTEGER            I     , IX    , IY
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO  = 0.0D+0 )

      IF( N    .LT.1    )RETURN

      IF( ( INCX.EQ.INCY ).AND.( INCX.GT.0 ) )THEN
         DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
            Y( IX ) = X( IX ) + Y( IX )
   10    CONTINUE
      ELSE
         IF( INCY.GE.0 )THEN
            IY = 1
         ELSE
            IY = 1 - ( N - 1 )*INCY
         END IF
         IF( INCX.GT.0 )THEN
            DO 20, IX = 1, 1 + ( N - 1 )*INCX, INCX
               Y( IY ) = X( IX ) + Y( IY )
               IY      = IY + INCY
   20       CONTINUE
         ELSE
            IX = 1 - ( N - 1 )*INCX
            DO 30, I = 1, N
               Y( IY ) = X( IX ) + Y( IY )
               IX      = IX + INCX
               IY      = IY + INCY
   30       CONTINUE
         END IF
      END IF
      RETURN

*     End of DXPY .

      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE DXMY  ( N, X, INCX, Y, INCY )
      INTEGER            N, INCX, INCY
      DOUBLE PRECISION   X( * ), Y( * )

C  DXMY  performs the operation
C
C     y := y - x
C
C
C  Nag Fortran 77 version of the Blas routine DAXPY .
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 3-September-1982.
C     Sven Hammarling, Nag Central Office.

      INTEGER            I     , IX    , IY
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO  = 0.0D+0 )

      IF( N    .LT.1    )RETURN

      IF( ( INCX.EQ.INCY ).AND.( INCX.GT.0 ) )THEN
         DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
            Y( IX ) = -X( IX ) + Y( IX )
   10    CONTINUE
      ELSE
         IF( INCY.GE.0 )THEN
            IY = 1
         ELSE
            IY = 1 - ( N - 1 )*INCY
         END IF
         IF( INCX.GT.0 )THEN
            DO 20, IX = 1, 1 + ( N - 1 )*INCX, INCX
               Y( IY ) = -X( IX ) + Y( IY )
               IY      = IY + INCY
   20       CONTINUE
         ELSE
            IX = 1 - ( N - 1 )*INCX
            DO 30, I = 1, N
               Y( IY ) = -X( IX ) + Y( IY )
               IX      = IX + INCX
               IY      = IY + INCY
   30       CONTINUE
         END IF
      END IF
      RETURN

*     End of DXMY .

      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE DCOPY ( N, X, INCX, Y, INCY )
      INTEGER            N, INCX, INCY
      DOUBLE PRECISION   X( * ), Y( * )

C  DCOPY  performs the operation
C
C     y := x
C
C
C  Nag Fortran 77 version of the Blas routine DCOPY .
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 26-November-1982.
C     Sven Hammarling, Nag Central Office.

      INTEGER            I     , IX    , IY

      IF( N.LT.1 )RETURN

      IF( ( INCX.EQ.INCY ).AND.( INCY.GT.0 ) )THEN
         DO 10, IY = 1, 1 + ( N - 1 )*INCY, INCY
            Y( IY ) = X( IY )
   10    CONTINUE
      ELSE
         IF( INCX.GE.0 )THEN
            IX = 1
         ELSE
            IX = 1 - ( N - 1 )*INCX
         END IF
         IF( INCY.GT.0 )THEN
            DO 20, IY = 1, 1 + ( N - 1 )*INCY, INCY
               Y( IY ) = X( IX )
               IX      = IX + INCX
   20       CONTINUE
         ELSE
            IY = 1 - ( N - 1 )*INCY
            DO 30, I = 1, N
               Y( IY ) = X( IX )
               IY      = IY + INCY
               IX      = IX + INCX
   30       CONTINUE
         END IF
      END IF

      RETURN

*     End of DCOPY .

      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      DOUBLE PRECISION FUNCTION DDOT  ( N, X, INCX, Y, INCY )
      INTEGER                           N, INCX, INCY
      DOUBLE PRECISION                  X( * ), Y( * )

C  DDOT   returns the value
C
C     DDOT   = x'y
C
C
C  Nag Fortran 77 version of the Blas routine DDOT  .
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 21-September-1982.
C     Sven Hammarling, Nag Central Office.

      INTEGER             I     , IX    , IY
      DOUBLE PRECISION    SUM
      DOUBLE PRECISION    ZERO
      PARAMETER         ( ZERO  = 0.0D+0 )

      SUM = ZERO
      IF( N.GE.1 )THEN
         IF( ( INCX.EQ.INCY ).AND.( INCX.GT.0 ) )THEN
            DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
               SUM = SUM + X( IX )*Y( IX )
   10       CONTINUE
         ELSE
            IF( INCY.GE.0 )THEN
               IY = 1
            ELSE
               IY = 1 - ( N - 1 )*INCY
            END IF
            IF( INCX.GT.0 )THEN
               DO 20, IX = 1, 1 + ( N - 1 )*INCX, INCX
                  SUM = SUM + X( IX )*Y( IY )
                  IY  = IY  + INCY
   20          CONTINUE
            ELSE
               IX = 1 - ( N - 1 )*INCX
               DO 30, I = 1, N
                  SUM = SUM + X( IX )*Y( IY )
                  IX  = IX  + INCX
                  IY  = IY  + INCY
   30          CONTINUE
            END IF
         END IF
      END IF

      DDOT   = SUM
      RETURN

*     End of DDOT  .

      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      DOUBLE PRECISION FUNCTION DNRM2 ( N, X, INCX )
      INTEGER                           N, INCX
      DOUBLE PRECISION                  X( * )

C  DNRM2  returns the Euclidean norm of a vector via the function
C  name, so that
C
C     DNRM2  := sqrt( x'*x )
C
C
C  Nag Fortran 77 version of the Blas routine DNRM2 .
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 25-October-1982.
C     Sven Hammarling, Nag Central Office.

      EXTERNAL            DNORM , DSSQ
      INTRINSIC           DABS
      DOUBLE PRECISION    SCALE , DNORM , SSQ
      DOUBLE PRECISION    ONE   ,         ZERO
      PARAMETER         ( ONE   = 1.0D+0, ZERO  = 0.0D+0 )

      IF( N.LT.1 )THEN
         DNRM2  = ZERO
      ELSE IF( N.EQ.1 )THEN
         DNRM2  = DABS( X( 1 ) )
      ELSE
         SCALE  = ZERO
         SSQ    = ONE

         CALL DSSQ  ( N, X, INCX, SCALE, SSQ )

         DNRM2  = DNORM ( SCALE, SSQ )

      END IF
      RETURN

*     End of DNRM2 .

      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE DSCAL ( N, ALPHA, X, INCX )
      INTEGER            N, INCX
      DOUBLE PRECISION   ALPHA
      DOUBLE PRECISION   X( * )

C  DSCAL  performs the operation
C
C     x := alpha*x
C
C
C  Nag Fortran 77 version of the Blas routine DSCAL .
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 26-November-1982.
C     Sven Hammarling, Nag Central Office.

      INTEGER            IX
      DOUBLE PRECISION   ONE   ,         ZERO
      PARAMETER        ( ONE   = 1.0D+0, ZERO  = 0.0D+0 )

      IF( N.GE.1 )THEN
         IF( ALPHA.EQ.ZERO )THEN
            DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
               X( IX ) = ZERO
   10       CONTINUE
         ELSE IF( ALPHA.EQ.( -ONE ) )THEN
            DO 20, IX = 1, 1 + ( N - 1 )*INCX, INCX
               X( IX ) = -X( IX )
   20       CONTINUE
         ELSE IF( ALPHA.NE.ONE )THEN
            DO 30, IX = 1, 1 + ( N - 1 )*INCX, INCX
               X( IX ) = ALPHA*X( IX )
   30       CONTINUE
         END IF
      END IF

      RETURN

*     End of DSCAL .

      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE DSWAP ( N, X, INCX, Y, INCY )
      INTEGER            N, INCX, INCY
      DOUBLE PRECISION   X( * ), Y( * )

C  DSWAP  performs the operations
C
C     temp := x,   x := y,   y := temp.
C
C
C  Nag Fortran 77 version of the Blas routine DSWAP .
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 26-November-1982.
C     Sven Hammarling, Nag Central Office.

      INTEGER            I     , IX    , IY
      DOUBLE PRECISION   TEMP

      IF( N.LT.1 )RETURN

      IF( ( INCX.EQ.INCY ).AND.( INCY.GT.0 ) )THEN
         DO 10, IY = 1, 1 + ( N - 1 )*INCY, INCY
            TEMP    = X( IY )
            X( IY ) = Y( IY )
            Y( IY ) = TEMP
   10    CONTINUE
      ELSE
         IF( INCX.GE.0 )THEN
            IX = 1
         ELSE
            IX = 1 - ( N - 1 )*INCX
         END IF
         IF( INCY.GT.0 )THEN
            DO 20, IY = 1, 1 + ( N - 1 )*INCY, INCY
               TEMP    = X( IX )
               X( IX ) = Y( IY )
               Y( IY ) = TEMP
               IX      = IX + INCX
   20       CONTINUE
         ELSE
            IY = 1 - ( N - 1 )*INCY
            DO 30, I = 1, N
               TEMP    = X( IX )
               X( IX ) = Y( IY )
               Y( IY ) = TEMP
               IY      = IY + INCY
               IX      = IX + INCX
   30       CONTINUE
         END IF
      END IF

      RETURN

*     End of DSWAP .

      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      INTEGER FUNCTION IDAMAX( N, X, INCX )
      INTEGER                  N, INCX
      DOUBLE PRECISION         X( * )

C  IDAMAX returns the smallest value of i such that
C
C     abs( x( i ) ) = max( abs( x( j ) ) )
C                      j
C
C  Nag Fortran 77 version of the Blas routine IDAMAX.
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 31-May-1983.
C     Sven Hammarling, Nag Central Office.

      INTRINSIC           DABS
      INTEGER             I     , IMAX  , IX
      DOUBLE PRECISION    XMAX

      IF( N.LT.1 )THEN
         IDAMAX = 0
         RETURN
      END IF

      IMAX = 1
      IF( N.GT.1 )THEN
         XMAX = DABS( X( 1 ) )
         IX   = 1
         DO 10, I = 2, N
            IX = IX + INCX
            IF( XMAX.LT.DABS( X( IX ) ) )THEN
               XMAX = DABS( X( IX ) )
               IMAX = I
            END IF
   10    CONTINUE
      END IF

      IDAMAX = IMAX
      RETURN

*     End of IDAMAX.

      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE DGEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
*     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA, BETA
      INTEGER            INCX, INCY, LDA, M, N
      CHARACTER*1        TRANS
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * ), Y( * )
*     ..
*
*  Purpose
*  =======
*
*  DGEMV  performs one of the matrix-vector operations
*
*     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
*
*  where alpha and beta are scalars, x and y are vectors and A is an
*  m by n matrix.
*
*  Parameters
*  ==========
*
*  TRANS  - CHARACTER*1.
*           On entry, TRANS specifies the operation to be performed as
*           follows:
*
*              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
*
*              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
*
*              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of the matrix A.
*           M must be at least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
*           Before entry, the leading m by n part of the array A must
*           contain the matrix of coefficients.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           max( 1, m ).
*           Unchanged on exit.
*
*  X      - DOUBLE PRECISION array of DIMENSION at least
*           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
*           and at least
*           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
*           Before entry, the incremented array X must contain the
*           vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*  BETA   - DOUBLE PRECISION.
*           On entry, BETA specifies the scalar beta. When BETA is
*           supplied as zero then Y need not be set on input.
*           Unchanged on exit.
*
*  Y      - DOUBLE PRECISION array of DIMENSION at least
*           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
*           and at least
*           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
*           Before entry with BETA non-zero, the incremented array Y
*           must contain the vector y. On exit, Y is overwritten by the
*           updated vector y.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y. INCY must not be zero.
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
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     .. Local Scalars ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, INFO, IX, IY, J, JX, JY, KX, KY, LENX, LENY
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
      IF      ( .NOT.LSAME( TRANS, 'N' ).AND.
     $          .NOT.LSAME( TRANS, 'T' ).AND.
     $          .NOT.LSAME( TRANS, 'C' )      ) THEN
         INFO = 1
      ELSE IF ( M.LT.0 ) THEN
         INFO = 2
      ELSE IF ( N.LT.0 ) THEN
         INFO = 3
      ELSE IF ( LDA.LT.MAX(1,M) ) THEN
         INFO = 6
      ELSE IF ( INCX.EQ.0 ) THEN
         INFO = 8
      ELSE IF ( INCY.EQ.0 ) THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DGEMV ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.
     $    ( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
*
*     Set LENX and LENY, the lengths of the vectors x and y.
*
      IF( LSAME( TRANS, 'N' ) )THEN
         LENX = N
         LENY = M
      ELSE
         LENX = M
         LENY = N
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
*     First form  y := beta*y  and set up the start points in X and Y if
*     the increments are not both unity.
*
      IF( ( INCX.EQ.1 ).AND.( INCY.EQ.1 ) )THEN
         IF( BETA.NE.ONE )THEN
            IF( BETA.EQ.ZERO )THEN
               DO 10, I = 1, LENY
                  Y( I ) = ZERO
   10          CONTINUE
            ELSE
               DO 20, I = 1, LENY
                  Y( I ) = BETA*Y( I )
   20          CONTINUE
            END IF
         END IF
      ELSE
         IF( INCX.GT.0 )THEN
            KX = 1
         ELSE
            KX = 1 - ( LENX - 1 )*INCX
         END IF
         IF( INCY.GT.0 )THEN
            KY = 1
         ELSE
            KY = 1 - ( LENY - 1 )*INCY
         END IF
         IF( BETA.NE.ONE )THEN
            IY = KY
            IF( BETA.EQ.ZERO )THEN
               DO 30, I = 1, LENY
                  Y( IY ) = ZERO
                  IY      = IY   + INCY
   30          CONTINUE
            ELSE
               DO 40, I = 1, LENY
                  Y( IY ) = BETA*Y( IY )
                  IY      = IY           + INCY
   40          CONTINUE
            END IF
         END IF
      END IF
      IF( ALPHA.EQ.ZERO )
     $   RETURN
      IF( LSAME( TRANS, 'N' ) )THEN
*
*        Form  y := alpha*A*x + y.
*
         IF( ( INCX.EQ.1 ).AND.( INCY.EQ.1 ) )THEN
            DO 60, J = 1, N
               IF( X( J ).NE.ZERO )THEN
                  TEMP = ALPHA*X( J )
                  DO 50, I = 1, M
                     Y( I ) = Y( I ) + TEMP*A( I, J )
   50             CONTINUE
               END IF
   60       CONTINUE
         ELSE
            JX = KX
            DO 80, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  IY   = KY
                  DO 70, I = 1, M
                     Y( IY ) = Y( IY ) + TEMP*A( I, J )
                     IY      = IY      + INCY
   70             CONTINUE
               END IF
               JX = JX + INCX
   80       CONTINUE
         END IF
      ELSE
*
*        Form  y := alpha*A'*x + y.
*
         IF( ( INCX.EQ.1 ).AND.( INCY.EQ.1 ) )THEN
            DO 100, J = 1, N
               TEMP = ZERO
               DO 90, I = 1, M
                  TEMP = TEMP + A( I, J )*X( I )
   90          CONTINUE
               Y( J ) = Y( J ) + ALPHA*TEMP
  100       CONTINUE
         ELSE
            JY = KY
            DO 120, J = 1, N
               TEMP = ZERO
               IX   = KX
               DO 110, I = 1, M
                  TEMP = TEMP + A( I, J )*X( IX )
                  IX   = IX   + INCX
  110          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  120       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of DGEMV .
*
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE DGERM1( M, N, X, INCX, Y, INCY, A, LDA )
C     SUBROUTINE DGER  ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
*     .. Scalar Arguments ..
C     DOUBLE PRECISION   ALPHA
      INTEGER            INCX, INCY, LDA, M, N
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * ), Y( * )
*     ..
*
*  Purpose
*  =======
*
*  DGER   performs the rank 1 operation
*
*     A := alpha*x*y' + A,
*
*  where alpha is a scalar, x is an m element vector, y is an n element
*  vector and A is an m by n matrix.
*
*  Parameters
*  ==========
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of the matrix A.
*           M must be at least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  X      - DOUBLE PRECISION array of dimension at least
*           ( 1 + ( m - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the m
*           element vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*  Y      - DOUBLE PRECISION array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCY ) ).
*           Before entry, the incremented array Y must contain the n
*           element vector y.
*           Unchanged on exit.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y. INCY must not be zero.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
*           Before entry, the leading m by n part of the array A must
*           contain the matrix of coefficients. On exit, A is
*           overwritten by the updated matrix.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           max( 1, m ).
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
      DOUBLE PRECISION   TEMP
      INTEGER            I, INFO, IX, J, JY, KX
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
      IF      ( M.LT.0 ) THEN
         INFO = 1
      ELSE IF ( N.LT.0 ) THEN
         INFO = 2
      ELSE IF ( INCX.EQ.0 ) THEN
         INFO = 5
      ELSE IF ( INCY.EQ.0 ) THEN
         INFO = 7
      ELSE IF ( LDA.LT.MAX(1,M) ) THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DGER  ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
C     IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.( ALPHA.EQ.ZERO ) )
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ) )
     $   RETURN
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
      IF( ( INCX.EQ.1 ).AND.( INCY.EQ.1 ) )THEN
         DO 20, J = 1, N
            IF( Y( J ).NE.ZERO )THEN
               TEMP = -Y( J )
C              TEMP = ALPHA*Y( J )
               DO 10, I = 1, M
                  A( I, J ) = A( I, J ) + X( I )*TEMP
   10          CONTINUE
            END IF
   20    CONTINUE
      ELSE
         IF( INCX.GT.0 )THEN
            KX = 1
         ELSE
            KX = 1 - ( M - 1 )*INCX
         END IF
         IF( INCY.GT.0 )THEN
            JY = 1
         ELSE
            JY = 1 - ( N - 1 )*INCY
         END IF
         DO 40, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = -Y( JY )
C              TEMP = ALPHA*Y( JY )
               IX   = KX
               DO 30, I = 1, M
                  A( I, J ) = A( I, J ) + X( IX )*TEMP
                  IX        = IX        + INCX
   30          CONTINUE
            END IF
            JY = JY + INCY
   40    CONTINUE
      END IF
*
      RETURN
*
*     End of DGERM1.
*
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE DSYMV ( UPLO, N, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
*     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA, BETA
      INTEGER            INCX, INCY, LDA, N
      CHARACTER*1        UPLO
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * ), Y( * )
*     ..
*
*  Purpose
*  =======
*
*  DSYMV  performs the matrix-vector  operation
*
*     y := alpha*A*x + beta*y,
*
*  where alpha and beta are scalars, x and y are n element vectors and
*  A is an n by n symmetric matrix.
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
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
*           Before entry with  UPLO = 'U' or 'u', the leading n by n
*           upper triangular part of the array A must contain the upper
*           triangular part of the symmetric matrix and the strictly
*           lower triangular part of A is not referenced.
*           Before entry with UPLO = 'L' or 'l', the leading n by n
*           lower triangular part of the array A must contain the lower
*           triangular part of the symmetric matrix and the strictly
*           upper triangular part of A is not referenced.
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
*           element vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*  BETA   - DOUBLE PRECISION.
*           On entry, BETA specifies the scalar beta. When BETA is
*           supplied as zero then Y need not be set on input.
*           Unchanged on exit.
*
*  Y      - DOUBLE PRECISION array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCY ) ).
*           Before entry, the incremented array Y must contain the n
*           element vector y. On exit, Y is overwritten by the updated
*           vector y.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y. INCY must not be zero.
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
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     .. Local Scalars ..
      DOUBLE PRECISION   TEMP1, TEMP2
      INTEGER            I, INFO, IX, IY, J, JX, JY, KX, KY
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
      ELSE IF ( LDA.LT.MAX(1,N) ) THEN
         INFO = 5
      ELSE IF ( INCX.EQ.0 ) THEN
         INFO = 7
      ELSE IF ( INCY.EQ.0 ) THEN
         INFO = 10
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DSYMV ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( N.EQ.0 ).OR.( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through the triangular part
*     of A.
*
*     First form  y := beta*y  and set up the start points in X and Y if
*     the increments are not both unity.
*
      IF( ( INCX.EQ.1 ).AND.( INCY.EQ.1 ) )THEN
         IF( BETA.NE.ONE )THEN
            IF( BETA.EQ.ZERO )THEN
               DO 10, I = 1, N
                  Y( I ) = ZERO
   10          CONTINUE
            ELSE
               DO 20, I = 1, N
                  Y( I ) = BETA*Y( I )
   20          CONTINUE
            END IF
         END IF
      ELSE
         IF( INCX.GT.0 )THEN
            KX = 1
         ELSE
            KX = 1 - ( N - 1 )*INCX
         END IF
         IF( INCY.GT.0 )THEN
            KY = 1
         ELSE
            KY = 1 - ( N - 1 )*INCY
         END IF
         IF( BETA.NE.ONE )THEN
            IY = KY
            IF( BETA.EQ.ZERO )THEN
               DO 30, I = 1, N
                  Y( IY ) = ZERO
                  IY      = IY   + INCY
   30          CONTINUE
            ELSE
               DO 40, I = 1, N
                  Y( IY ) = BETA*Y( IY )
                  IY      = IY           + INCY
   40          CONTINUE
            END IF
         END IF
      END IF
      IF( ALPHA.EQ.ZERO )
     $   RETURN
      IF( LSAME( UPLO, 'U' ) )THEN
*
*        Form  y  when A is stored in upper triangle.
*
         IF( ( INCX.EQ.1 ).AND.( INCY.EQ.1 ) )THEN
            DO 60, J = 1, N
               TEMP1 = ALPHA*X( J )
               TEMP2 = ZERO
               DO 50, I = 1, J - 1
                  Y( I ) = Y( I ) + TEMP1*A( I, J )
                  TEMP2  = TEMP2  + A( I, J )*X( I )
   50          CONTINUE
               Y( J ) = Y( J ) + TEMP1*A( J, J ) + ALPHA*TEMP2
   60       CONTINUE
         ELSE
            IX = KX - INCX
            DO 80, J = 1, N
               TEMP1 = ALPHA*X( IX + INCX )
               TEMP2 = ZERO
               IX    = KX
               IY    = KY
               DO 70, I = 1, J - 1
                  Y( IY ) = Y( IY ) + TEMP1*A( I, J )
                  TEMP2   = TEMP2   + A( I, J )*X( IX )
                  IX      = IX      + INCX
                  IY      = IY      + INCY
   70          CONTINUE
               Y( IY ) = Y( IY ) + TEMP1*A( J, J ) + ALPHA*TEMP2
   80       CONTINUE
         END IF
      ELSE
*
*        Form  y  when A is stored in lower triangle.
*
         IF( ( INCX.EQ.1 ).AND.( INCY.EQ.1 ) )THEN
            DO 100, J = 1, N
               TEMP1  = ALPHA*X( J )
               TEMP2  = ZERO
               Y( J ) = Y( J )       + TEMP1*A( J, J )
               DO 90, I = J + 1, N
                  Y( I ) = Y( I ) + TEMP1*A( I, J )
                  TEMP2  = TEMP2  + A( I, J )*X( I )
   90          CONTINUE
               Y( J ) = Y( J ) + ALPHA*TEMP2
  100       CONTINUE
         ELSE
            JX = KX
            JY = KY
            DO 120, J = 1, N
               TEMP1   = ALPHA*X( JX )
               TEMP2   = ZERO
               Y( JY ) = Y( JY )       + TEMP1*A( J, J )
               IX      = JX
               IY      = JY
               DO 110, I = J + 1, N
                  IX      = IX      + INCX
                  IY      = IY      + INCY
                  Y( IY ) = Y( IY ) + TEMP1*A( I, J )
                  TEMP2   = TEMP2   + A( I, J )*X( IX )
  110          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP2
               JX      = JX      + INCX
               JY      = JY      + INCY
  120       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of DSYMV .
*
      END
