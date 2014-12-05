      SUBROUTINE Int2_IncFun
!
      USE MP2_Constant_Module, ONLY : Zero, One, Half, Two, Pi
      USE Int2_Gamma_Module, ONLY : FF, FTab, MaxGam, MaxtuvGam
!
      IMPLICIT NONE
!
      REAL(8), PARAMETER :: TM78 = 1.0D-78, TM29 = 1.0D-29
      INTEGER :: N, M, I, J, MaxDim
      REAL(8) :: X, XX, E, QXX, Tol, FacMin, Term, Sum, Fac, T, S, A, TMax, SqPie4
!
      SqPie4 = SQRT(Pi * 0.25D0)   ! SqPie4 = SQRT(Pie4) = SQRT(0.78539816339744831D+00)
!
!YA062312      N = MaxtuvGam - 1
      N = (MaxtuvGam + 1) - 1   ! YA062312
      Tol = TM29
      MaxDim = MaxGam
!
!     o Set the FFunc table
!
      X = Zero
      DO I = 1, MaxDim
!
         XX = X + X
         FacMin = XX
         E = TM78
         QXX = -X
         IF (FacMin < 360.0D0) E = EXP(QXX)
         IF (FacMin >  80.0D0) GO TO 100
         Term = One
         Sum = One
         Fac = DBLE(N)
         Fac = Fac + Half
   10    Fac = Fac + One
         Term = Term * X / Fac
         Sum = Sum + Term
         IF (Fac <= FacMin) GO TO 10
         T = Term
         S = Sum
         IF (T > S * Tol) GO TO 10
         Fac = DBLE(N + N + 1)
         FF(N+1) = Sum * E / Fac
         M = N - 1
         Fac = DBLE(M + M + 1)
   20    IF (M < 0) GO TO 130
         FF(M+1) = (E + XX * FF(M+2)) / Fac
         M = M - 1
         Fac = Fac - Two
         GO TO 20
!
!     o Asymptotic expansion for large arguments
!
  100    A = SqPie4 / SQRT(X)   ! A = SQRT(Pie4 / X)
         TMax = A * Tol / E
         Term = One / XX
         Sum = Term
         Fac = One
  110    Fac = Fac - Two
         Term = Fac * Term / XX
         Sum = Term + Sum
         T = Term
         IF (ABS(T) > TMax) GO TO 110
         FF(1) = A - E * Sum
         Fac = -One
         M = 0
  120    IF (M == N) GO TO 130
         M = M + 1
         Fac = Fac + Two
         FF(M+1) = (Fac * FF(M) - E) / XX
         GO TO 120
  130    CONTINUE
!
!YA062312         DO J = 1, MaxtuvGam
         DO J = 1, (MaxtuvGam + 1)   ! YA062312
            FTab(J,I) = FF(J)
         END DO
!
         X = X + 0.05D+00
!
      END DO
!
      END SUBROUTINE
