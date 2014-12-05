      SUBROUTINE Int2_FFunc(X, NNN)
!
      USE MP2_Constant_Module, ONLY : Zero, One, Half, Pi
      USE Int2_Module, ONLY : FMM_LongRange, FMM_shortRange
      USE Int2_Gamma_Module, ONLY : FF, FTab, MaxGam
!
      IMPLICIT NONE
!
!     o Generate FFunc value by seven-term Taylor expansion
!
      INTEGER :: NNN
      REAL(8) :: X
!
      INTEGER :: K0, K, J, N21, KK
      REAL(8) :: XI, TN, TAT, Sum
!
      IF (X < Zero) THEN
         STOP 'Error: X < 0 in FFunc'
      END IF
!
      IF (FMM_LongRange .AND. (.NOT. FMM_ShortRange)) GO TO 120
!
      K0 = INT((X + 0.025D0) * 20.D0)
      K = K0 + 1
      IF (K <= 0) GO TO 120
      IF (K <= MaxGam) GO TO 320
!
!     o Asymptotic expansion
!
  120 CONTINUE
      XI = One / (X + X)
      J = 1
      TN = SQRT(Pi / X) * Half
      FF(J) = TN
  100 IF (J >= NNN) GO TO 500
      J = J + 1
      N21 = J + J - 3
      TN = TN * XI * N21
      FF(J) = TN
      GO TO 100
!
!     o Seven-term Taylor expansion
!
  320 TAT = 0.05D0 * DBLE(K0) - X
!
      DO KK = 1, NNN
         Sum = FTab(KK  ,K) +                       TAT * (                               &
     &         FTab(KK+1,K) +               0.5D0 * TAT * (                               &
     &         FTab(KK+2,K) + 0.333333333333333D0 * TAT * (                               &
     &         FTab(KK+3,K) +              0.25D0 * TAT * (                               &
     &         FTab(KK+4,K) +               0.2D0 * TAT * (                               &
     &         FTab(KK+5,K) + 0.166666666666666D0 * TAT *                                 &
     &         FTab(KK+6,K)                                )))))
         FF(KK) = Sum
      END DO
!
  500 CONTINUE
!
      IF (FMM_ShortRange .AND. (.NOT. FMM_LongRange)) THEN
         XI = One / (X + X)
         J = 1
         TN = SQRT(Pi / X) * Half
         FF(J) = FF(J) - TN
  600    IF (J >= NNN) GO TO 700
         J = J + 1
         N21 = J + J - 3
         TN = TN * XI * N21
         FF(J) = FF(J) - TN
         GO TO 600
      END IF
  700 CONTINUE
!
      END SUBROUTINE
