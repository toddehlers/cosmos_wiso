      SUBROUTINE DLAED2( K, N, D, Q, LDQ, INDXQ, RHO, CUTPNT, Z, DLAMDA,
     $                   Q2, LDQ2, INDXC, W, INDXP, INDX, COLTYP, INFO )
*
*  -- LAPACK routine (version 2.0) --
*     Univ. of Tennessee, Oak Ridge National Lab, Argonne National Lab,
*     Courant Institute, NAG Ltd., and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      INTEGER            CUTPNT, INFO, K, LDQ, LDQ2, N
      DOUBLE PRECISION   RHO
*     ..
*     .. Array Arguments ..
      INTEGER            COLTYP( * ), INDX( * ), INDXC( * ), INDXP( * ),
     $                   INDXQ( * )
      DOUBLE PRECISION   D( * ), DLAMDA( * ), Q( LDQ, * ),
     $                   Q2( LDQ2, * ), W( * ), Z( * )
*     ..
*
*  Purpose
*  =======
*
*  DLAED2 merges the two sets of eigenvalues together into a single
*  sorted set.  Then it tries to deflate the size of the problem.
*  There are two ways in which deflation can occur:  when two or more
*  eigenvalues are close together or if there is a tiny entry in the
*  Z vector.  For each such occurrence the order of the related secular
*  equation problem is reduced by one.
*
*  Arguments
*  =========
*
*  K      (output) INTEGER
*         The number of non-deflated eigenvalues, and the order of the
*         related secular equation. 0 <= K <=N.
*
*  N      (input) INTEGER
*         The dimension of the symmetric tridiagonal matrix.  N >= 0.
*
*  D      (input/output) DOUBLE PRECISION array, dimension (N)
*         On entry, D contains the eigenvalues of the two submatrices to
*         be combined.
*         On exit, D contains the trailing (N-K) updated eigenvalues
*         (those which were deflated) sorted into increasing order.
*
*  Q      (input/output) DOUBLE PRECISION array, dimension (LDQ, N)
*         On entry, Q contains the eigenvectors of two submatrices in
*         the two square blocks with corners at (1,1), (CUTPNT,CUTPNT)
*         and (CUTPNT+1, CUTPNT+1), (N,N).
*         On exit, Q contains the trailing (N-K) updated eigenvectors
*         (those which were deflated) in its last N-K columns.
*
*  LDQ    (input) INTEGER
*         The leading dimension of the array Q.  LDQ >= max(1,N).
*
*  INDXQ  (input/output) INTEGER array, dimension (N)
*         The permutation which separately sorts the two sub-problems
*         in D into ascending order.  Note that elements in the second
*         half of this permutation must first have CUTPNT added to their
*         values. Destroyed on exit.
*
*  RHO    (input/output) DOUBLE PRECISION
*         On entry, the off-diagonal element associated with the rank-1
*         cut which originally split the two submatrices which are now
*         being recombined.
*         On exit, RHO has been modified to the value required by
*         DLAED3.
*
*  CUTPNT (input) INTEGER
*         The location of the last eigenvalue in the leading sub-matrix.
*         min(1,N) <= CUTPNT <= N.
*
*  Z      (input) DOUBLE PRECISION array, dimension (N)
*         On entry, Z contains the updating vector (the last
*         row of the first sub-eigenvector matrix and the first row of
*         the second sub-eigenvector matrix).
*         On exit, the contents of Z have been destroyed by the updating
*         process.
*
*  DLAMDA (output) DOUBLE PRECISION array, dimension (N)
*         A copy of the first K eigenvalues which will be used by
*         DLAED3 to form the secular equation.
*
*  Q2     (output) DOUBLE PRECISION array, dimension (LDQ2, N)
*         A copy of the first K eigenvectors which will be used by
*         DLAED3 in a matrix multiply (DGEMM) to solve for the new
*         eigenvectors.   Q2 is arranged into three blocks.  The
*         first block contains non-zero elements only at and above
*         CUTPNT, the second contains non-zero elements only below
*         CUTPNT, and the third is dense.
*
*  LDQ2   (input) INTEGER
*         The leading dimension of the array Q2.  LDQ2 >= max(1,N).
*
*  INDXC  (output) INTEGER array, dimension (N)
*         The permutation used to arrange the columns of the deflated
*         Q matrix into three groups:  the first group contains non-zero
*         elements only at and above CUTPNT, the second contains
*         non-zero elements only below CUTPNT, and the third is dense.
*
*  W      (output) DOUBLE PRECISION array, dimension (N)
*         The first k values of the final deflation-altered z-vector
*         which will be passed to DLAED3.
*
*  INDXP  (workspace) INTEGER array, dimension (N)
*         The permutation used to place deflated values of D at the end
*         of the array.  INDXP(1:K) points to the nondeflated D-values
*         and INDXP(K+1:N) points to the deflated eigenvalues.
*
*  INDX   (workspace) INTEGER array, dimension (N)
*         The permutation used to sort the contents of D into ascending
*         order.
*
*  COLTYP (workspace/output) INTEGER array, dimension (N)
*         During execution, a label which will indicate which of the
*         following types a column in the Q2 matrix is:
*         1 : non-zero in the upper half only;
*         2 : non-zero in the lower half only;
*         3 : dense;
*         4 : deflated.
*         On exit, COLTYP(i) is the number of columns of type i,
*         for i=1 to 4 only.
*
*  INFO   (output) INTEGER
*          = 0:  successful exit.
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   MONE, ZERO, ONE, TWO, EIGHT
      PARAMETER          ( MONE = -1.0D0, ZERO = 0.0D0, ONE = 1.0D0,
     $                   TWO = 2.0D0, EIGHT = 8.0D0 )
*     ..
*     .. Local Arrays ..
      INTEGER            CTOT( 4 ), PSM( 4 )
*     ..
*     .. Local Scalars ..
      INTEGER            CT, I, IMAX, J, JLAM, JMAX, JP, K2, N1, N1P1,
     $                   N2
      DOUBLE PRECISION   C, EPS, S, T, TAU, TOL
*     ..
*     .. External Functions ..
      INTEGER            IDAMAX
      DOUBLE PRECISION   DLAMCH, DLAPY2
      EXTERNAL           IDAMAX, DLAMCH, DLAPY2
*     ..
*     .. External Subroutines ..
      EXTERNAL           DCOPY, DLACPY, DLAMRG, DROT, DSCAL, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
*
      IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDQ.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( MIN( 1, N ).GT.CUTPNT .OR. N.LT.CUTPNT ) THEN
         INFO = -8
      ELSE IF( LDQ2.LT.MAX( 1, N ) ) THEN
         INFO = -12
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLAED2', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      N1 = CUTPNT
      N2 = N - N1
      N1P1 = N1 + 1
*
      IF( RHO.LT.ZERO ) THEN
         CALL DSCAL( N2, MONE, Z( N1P1 ), 1 )
      END IF
*
*     Normalize z so that norm(z) = 1.  Since z is the concatenation of
*     two normalized vectors, norm2(z) = sqrt(2).
*
      T = ONE / SQRT( TWO )
      DO 10 J = 1, N
         INDX( J ) = J
   10 CONTINUE
      CALL DSCAL( N, T, Z, 1 )
*
*     RHO = ABS( norm(z)**2 * RHO )
*
      RHO = ABS( TWO*RHO )
*
      DO 20 I = 1, CUTPNT
         COLTYP( I ) = 1
   20 CONTINUE
      DO 30 I = CUTPNT + 1, N
         COLTYP( I ) = 2
   30 CONTINUE
*
*     Sort the eigenvalues into increasing order
*
      DO 40 I = CUTPNT + 1, N
         INDXQ( I ) = INDXQ( I ) + CUTPNT
   40 CONTINUE
*
*     re-integrate the deflated parts from the last pass
*
      DO 50 I = 1, N
         DLAMDA( I ) = D( INDXQ( I ) )
         W( I ) = Z( INDXQ( I ) )
         INDXC( I ) = COLTYP( INDXQ( I ) )
   50 CONTINUE
      CALL DLAMRG( N1, N2, DLAMDA, 1, 1, INDX )
      DO 60 I = 1, N
         D( I ) = DLAMDA( INDX( I ) )
         Z( I ) = W( INDX( I ) )
         COLTYP( I ) = INDXC( INDX( I ) )
   60 CONTINUE
*
*     Calculate the allowable deflation tolerance
*
      IMAX = IDAMAX( N, Z, 1 )
      JMAX = IDAMAX( N, D, 1 )
      EPS = DLAMCH( 'Epsilon' )
      TOL = EIGHT*EPS*MAX( ABS( D( JMAX ) ), ABS( Z( IMAX ) ) )
*
*     If the rank-1 modifier is small enough, no more needs to be done
*     except to reorganize Q so that its columns correspond with the
*     elements in D.
*
      IF( RHO*ABS( Z( IMAX ) ).LE.TOL ) THEN
         K = 0
         DO 70 J = 1, N
            CALL DCOPY( N, Q( 1, INDXQ( INDX( J ) ) ), 1, Q2( 1, J ),
     $                  1 )
   70    CONTINUE
         CALL DLACPY( 'A', N, N, Q2, LDQ2, Q, LDQ )
         GO TO 180
      END IF
*
*     If there are multiple eigenvalues then the problem deflates.  Here
*     the number of equal eigenvalues are found.  As each equal
*     eigenvalue is found, an elementary reflector is computed to rotate
*     the corresponding eigensubspace so that the corresponding
*     components of Z are zero in this new basis.
*
      K = 0
      K2 = N + 1
      DO 80 J = 1, N
         IF( RHO*ABS( Z( J ) ).LE.TOL ) THEN
*
*           Deflate due to small z component.
*
            K2 = K2 - 1
            INDXP( K2 ) = J
            COLTYP( J ) = 4
            IF( J.EQ.N )
     $         GO TO 120
         ELSE
            JLAM = J
            GO TO 90
         END IF
   80 CONTINUE
   90 CONTINUE
      J = J + 1
      IF( J.GT.N )
     $   GO TO 110
      IF( RHO*ABS( Z( J ) ).LE.TOL ) THEN
*
*        Deflate due to small z component.
*
         K2 = K2 - 1
         INDXP( K2 ) = J
         COLTYP( J ) = 4
      ELSE
*
*        Check if eigenvalues are close enough to allow deflation.
*
         S = Z( JLAM )
         C = Z( J )
*
*        Find sqrt(a**2+b**2) without overflow or
*        destructive underflow.
*
         TAU = DLAPY2( C, S )
         T = D( J ) - D( JLAM )
         C = C / TAU
         S = -S / TAU
         IF( ABS( T*C*S ).LE.TOL ) THEN
*
*           Deflation is possible.
*
            Z( J ) = TAU
            Z( JLAM ) = ZERO
            IF( COLTYP( J ).NE.COLTYP( JLAM ) )
     $         COLTYP( J ) = 3
            COLTYP( JLAM ) = 4
            CALL DROT( N, Q( 1, INDXQ( INDX( JLAM ) ) ), 1,
     $                 Q( 1, INDXQ( INDX( J ) ) ), 1, C, S )
            T = D( JLAM )*C**2 + D( J )*S**2
            D( J ) = D( JLAM )*S**2 + D( J )*C**2
            D( JLAM ) = T
            K2 = K2 - 1
            I = 1
  100       CONTINUE
            IF( K2+I.LE.N ) THEN
               IF( D( JLAM ).LT.D( INDXP( K2+I ) ) ) THEN
                  INDXP( K2+I-1 ) = INDXP( K2+I )
                  INDXP( K2+I ) = JLAM
                  I = I + 1
                  GO TO 100
               ELSE
                  INDXP( K2+I-1 ) = JLAM
               END IF
            ELSE
               INDXP( K2+I-1 ) = JLAM
            END IF
            JLAM = J
         ELSE
            K = K + 1
            W( K ) = Z( JLAM )
            DLAMDA( K ) = D( JLAM )
            INDXP( K ) = JLAM
            JLAM = J
         END IF
      END IF
      GO TO 90
  110 CONTINUE
*
*     Record the last eigenvalue.
*
      K = K + 1
      W( K ) = Z( JLAM )
      DLAMDA( K ) = D( JLAM )
      INDXP( K ) = JLAM
*
  120 CONTINUE
*
*     Count up the total number of the various types of columns, then
*     form a permutation which positions the four column types into
*     four uniform groups (although one or more of these groups may be
*     empty).
*
      DO 130 J = 1, 4
         CTOT( J ) = 0
  130 CONTINUE
      DO 140 J = 1, N
         CT = COLTYP( J )
         CTOT( CT ) = CTOT( CT ) + 1
  140 CONTINUE
*
*     PSM(*) = Position in SubMatrix (of types 1 through 4)
*
      PSM( 1 ) = 1
      PSM( 2 ) = 1 + CTOT( 1 )
      PSM( 3 ) = PSM( 2 ) + CTOT( 2 )
      PSM( 4 ) = PSM( 3 ) + CTOT( 3 )
*
*     Fill out the INDXC array so that the permutation which it induces
*     will place all type-1 columns first, all type-2 columns next,
*     then all type-3's, and finally all type-4's.
*
      DO 150 J = 1, N
         JP = INDXP( J )
         CT = COLTYP( JP )
         INDXC( PSM( CT ) ) = J
         PSM( CT ) = PSM( CT ) + 1
  150 CONTINUE
*
*     Sort the eigenvalues and corresponding eigenvectors into DLAMDA
*     and Q2 respectively.  The eigenvalues/vectors which were not
*     deflated go into the first K slots of DLAMDA and Q2 respectively,
*     while those which were deflated go into the last N - K slots.
*
      DO 160 J = 1, N
         JP = INDXP( J )
         DLAMDA( J ) = D( JP )
         CALL DCOPY( N, Q( 1, INDXQ( INDX( INDXP( INDXC( J ) ) ) ) ), 1,
     $               Q2( 1, J ), 1 )
  160 CONTINUE
*
*     The deflated eigenvalues and their corresponding vectors go back
*     into the last N - K slots of D and Q respectively.
*
      CALL DCOPY( N-K, DLAMDA( K+1 ), 1, D( K+1 ), 1 )
      CALL DLACPY( 'A', N, N-K, Q2( 1, K+1 ), LDQ2, Q( 1, K+1 ), LDQ )
*
*     Copy CTOT into COLTYP for referencing in DLAED3.
*
      DO 170 J = 1, 4
         COLTYP( J ) = CTOT( J )
  170 CONTINUE
*
  180 CONTINUE
      RETURN
*
*     End of DLAED2
*
      END
