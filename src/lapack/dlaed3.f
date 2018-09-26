      SUBROUTINE DLAED3( K, KSTART, KSTOP, N, D, Q, LDQ, RHO, CUTPNT,
     $                   DLAMDA, Q2, LDQ2, INDXC, CTOT, W, S, LDS,
     $                   INFO )
*
*  -- LAPACK routine (version 2.0) --
*     Univ. of Tennessee, Oak Ridge National Lab, Argonne National Lab,
*     Courant Institute, NAG Ltd., and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      INTEGER            CUTPNT, INFO, K, KSTART, KSTOP, LDQ, LDQ2, LDS,
     $                   N
      DOUBLE PRECISION   RHO
*     ..
*     .. Array Arguments ..
      INTEGER            CTOT( * ), INDXC( * )
      DOUBLE PRECISION   D( * ), DLAMDA( * ), Q( LDQ, * ),
     $                   Q2( LDQ2, * ), S( LDS, * ), W( * )
*     ..
*
*  Purpose
*  =======
*
*  DLAED3 finds the roots of the secular equation, as defined by the
*  values in D, W, and RHO, between KSTART and KSTOP.  It makes the
*  appropriate calls to DLAED4 and then updates the eigenvectors by
*  multiplying the matrix of eigenvectors of the pair of eigensystems
*  being combined by the matrix of eigenvectors of the K-by-K system
*  which is solved here.
*
*  This code makes very mild assumptions about floating point
*  arithmetic. It will work on machines with a guard digit in
*  add/subtract, or on those binary machines without guard digits
*  which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or Cray-2.
*  It could conceivably fail on hexadecimal or decimal machines
*  without guard digits, but we know of none.
*
*  Arguments
*  =========
*
*  K       (input) INTEGER
*          The number of terms in the rational function to be solved by
*          DLAED4.  K >= 0.
*
*  KSTART  (input) INTEGER
*  KSTOP   (input) INTEGER
*          The updated eigenvalues Lambda(I), KSTART <= I <= KSTOP
*          are to be computed.  1 <= KSTART <= KSTOP <= K.
*
*  N       (input) INTEGER
*          The number of rows and columns in the Q matrix.
*          N >= K (deflation may result in N>K).
*
*  D       (output) DOUBLE PRECISION array, dimension (N)
*          D(I) contains the updated eigenvalues for
*          KSTART <= I <= KSTOP.
*
*  Q       (output) DOUBLE PRECISION array, dimension (LDQ,N)
*          Initially the first K columns are used as workspace.
*          On output the columns KSTART to KSTOP contain
*          the updated eigenvectors.
*
*  LDQ     (input) INTEGER
*          The leading dimension of the array Q.  LDQ >= max(1,N).
*
*  RHO     (input) DOUBLE PRECISION
*          The value of the parameter in the rank one update equation.
*          RHO >= 0 required.
*
*  CUTPNT  (input) INTEGER
*          The location of the last eigenvalue in the leading submatrix.
*          min(1,N) <= CUTPNT <= N.
*
*  DLAMDA  (input/output) DOUBLE PRECISION array, dimension (K)
*          The first K elements of this array contain the old roots
*          of the deflated updating problem.  These are the poles
*          of the secular equation. May be changed on output by
*          having lowest order bit set to zero on Cray X-MP, Cray Y-MP,
*          Cray-2, or Cray C-90, as described above.
*
*  Q2      (input) DOUBLE PRECISION array, dimension (LDQ2, N)
*          The first K columns of this matrix contain the non-deflated
*          eigenvectors for the split problem.
*
*  LDQ2    (input) INTEGER
*          The leading dimension of the array Q2.  LDQ2 >= max(1,N).
*
*  INDXC   (input) INTEGER array, dimension (N)
*          The permutation used to arrange the columns of the deflated
*          Q matrix into three groups:  the first group contains
*          non-zero elements only at and above CUTPNT, the second
*          contains non-zero elements only below CUTPNT, and the third
*          is dense.  The rows of the eigenvectors found by DLAED4
*          must be likewise permuted before the matrix multiply can take
*          place.
*
*  CTOT    (input) INTEGER array, dimension (4)
*          A count of the total number of the various types of columns
*          in Q, as described in INDXC.  The fourth column type is any
*          column which has been deflated.
*
*  W       (input/output) DOUBLE PRECISION array, dimension (K)
*          The first K elements of this array contain the components
*          of the deflation-adjusted updating vector. Destroyed on
*          output.
*
*  S       (workspace) DOUBLE PRECISION array, dimension (LDS, K)
*          Will contain the eigenvectors of the repaired matrix which
*          will be multiplied by the previously accumulated eigenvectors
*          to update the system.
*
*  LDS     (input) INTEGER
*          The leading dimension of S.  LDS >= max(1,K).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit.
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*          > 0:  if INFO = 1, an eigenvalue did not converge
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D0, ZERO = 0.0D0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J, JC, KTEMP, PARTS
      DOUBLE PRECISION   TEMP
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMC3, DNRM2
      EXTERNAL           DLAMC3, DNRM2
*     ..
*     .. External Subroutines ..
      EXTERNAL           DCOPY, DGEMM, DLAED4, DLASET, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, SIGN, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
*
      IF( K.LT.0 ) THEN
         INFO = -1
      ELSE IF( KSTART.LT.1 .OR. KSTART.GT.MAX( 1, K ) ) THEN
         INFO = -2
      ELSE IF( MAX( 1, KSTOP ).LT.KSTART .OR. KSTOP.GT.MAX( 1, K ) )
     $          THEN
         INFO = -3
      ELSE IF( N.LT.K ) THEN
         INFO = -4
      ELSE IF( LDQ.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDQ2.LT.MAX( 1, N ) ) THEN
         INFO = -12
      ELSE IF( LDS.LT.MAX( 1, K ) ) THEN
         INFO = -17
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLAED3', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( K.EQ.0 )
     $   RETURN
*
*     Modify values DLAMDA(i) to make sure all DLAMDA(i)-DLAMDA(j) can
*     be computed with high relative accuracy (barring over/underflow).
*     This is a problem on machines without a guard digit in
*     add/subtract (Cray XMP, Cray YMP, Cray C 90 and Cray 2).
*     The following code replaces DLAMDA(I) by 2*DLAMDA(I)-DLAMDA(I),
*     which on any of these machines zeros out the bottommost
*     bit of DLAMDA(I) if it is 1; this makes the subsequent
*     subtractions DLAMDA(I)-DLAMDA(J) unproblematic when cancellation
*     occurs. On binary machines with a guard digit (almost all
*     machines) it does not change DLAMDA(I) at all. On hexadecimal
*     and decimal machines with a guard digit, it slightly
*     changes the bottommost bits of DLAMDA(I). It does not account
*     for hexadecimal or decimal machines without guard digits
*     (we know of none). We use a subroutine call to compute
*     2*DLAMBDA(I) to prevent optimizing compilers from eliminating
*     this code.
*
      DO 10 I = 1, N
         DLAMDA( I ) = DLAMC3( DLAMDA( I ), DLAMDA( I ) ) - DLAMDA( I )
   10 CONTINUE
*
      KTEMP = KSTOP - KSTART + 1
      DO 20 J = KSTART, KSTOP
         CALL DLAED4( K, J, DLAMDA, W, Q( 1, J ), RHO, D( J ), INFO )
*
*        If the zero finder fails, the computation is terminated.
*
         IF( INFO.NE.0 )
     $      GO TO 130
   20 CONTINUE
*
      IF( K.EQ.1 .OR. K.EQ.2 ) THEN
         DO 40 I = 1, K
            DO 30 J = 1, K
               JC = INDXC( J )
               S( J, I ) = Q( JC, I )
   30       CONTINUE
   40    CONTINUE
         GO TO 120
      END IF
*
*     Compute updated W.
*
      CALL DCOPY( K, W, 1, S, 1 )
*
*     Initialize W(I) = Q(I,I)
*
      CALL DCOPY( K, Q, LDQ+1, W, 1 )
      DO 70 J = 1, K
         DO 50 I = 1, J - 1
            W( I ) = W( I )*( Q( I, J ) / ( DLAMDA( I )-DLAMDA( J ) ) )
   50    CONTINUE
         DO 60 I = J + 1, K
            W( I ) = W( I )*( Q( I, J ) / ( DLAMDA( I )-DLAMDA( J ) ) )
   60    CONTINUE
   70 CONTINUE
      DO 80 I = 1, K
         W( I ) = SIGN( SQRT( -W( I ) ), S( I, 1 ) )
   80 CONTINUE
*
*     Compute eigenvectors of the modified rank-1 modification.
*
      DO 110 J = 1, K
         DO 90 I = 1, K
            Q( I, J ) = W( I ) / Q( I, J )
   90    CONTINUE
         TEMP = DNRM2( K, Q( 1, J ), 1 )
         DO 100 I = 1, K
            JC = INDXC( I )
            S( I, J ) = Q( JC, J ) / TEMP
  100    CONTINUE
  110 CONTINUE
*
*     Compute the updated eigenvectors.
*
  120 CONTINUE
*
      PARTS = 0
      IF( CTOT( 1 ).GT.0 ) THEN
         PARTS = PARTS + 1
         CALL DGEMM( 'N', 'N', CUTPNT, KTEMP, CTOT( 1 ), ONE,
     $               Q2( 1, 1 ), LDQ2, S( 1, KSTART ), LDS, ZERO,
     $               Q( 1, KSTART ), LDQ )
      END IF
      IF( CTOT( 2 ).GT.0 ) THEN
         PARTS = PARTS + 2
         CALL DGEMM( 'N', 'N', N-CUTPNT, KTEMP, CTOT( 2 ), ONE,
     $               Q2( 1+CUTPNT, 1+CTOT( 1 ) ), LDQ2,
     $               S( 1+CTOT( 1 ), KSTART ), LDS, ZERO,
     $               Q( 1+CUTPNT, KSTART ), LDQ )
      END IF
      IF( PARTS.EQ.1 )
     $   CALL DLASET( 'A', N-CUTPNT, KTEMP, ZERO, ZERO,
     $                Q( 1+CUTPNT, KSTART ), LDQ )
      IF( PARTS.EQ.2 )
     $   CALL DLASET( 'A', CUTPNT, KTEMP, ZERO, ZERO, Q( 1, KSTART ),
     $                LDQ )
      IF( CTOT( 3 ).GT.0 ) THEN
         IF( PARTS.GT.0 ) THEN
            CALL DGEMM( 'N', 'N', N, KTEMP, CTOT( 3 ), ONE,
     $                  Q2( 1, 1+CTOT( 1 )+CTOT( 2 ) ), LDQ2,
     $                  S( 1+CTOT( 1 )+CTOT( 2 ), KSTART ), LDS, ONE,
     $                  Q( 1, KSTART ), LDQ )
         ELSE
            CALL DGEMM( 'N', 'N', N, KTEMP, CTOT( 3 ), ONE,
     $                  Q2( 1, 1+CTOT( 1 )+CTOT( 2 ) ), LDQ2,
     $                  S( 1+CTOT( 1 )+CTOT( 2 ), KSTART ), LDS, ZERO,
     $                  Q( 1, KSTART ), LDQ )
         END IF
      END IF
*
  130 CONTINUE
      RETURN
*
*     End of DLAED3
*
      END
