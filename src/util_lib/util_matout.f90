      SUBROUTINE Util_MatOut(XX, N, M)
!
      IMPLICIT NONE
!
      INTEGER :: N, M
      REAL(8) :: XX(N,*)
!
      INTEGER, PARAMETER :: NColumn = 10
      INTEGER :: I, J, K, II
      INTEGER :: IAmari, IRow, IStart, IEnd
!
      IAmari = MOD(M, NColumn)
      IRow = M / NColumn
      IF (IAmari /= 0) THEN
         IRow = IRow + 1
      ELSE
         IAmari = NColumn
      END IF
!
      DO I = 1, IRow
         IStart = (I - 1) * NColumn + 1
         IF (I == IRow) THEN
            IEnd = (I - 1) * NColumn + IAmari
         ELSE
            IEnd = I * NColumn
         END IF
!         WRITE(*, '(4X, 10(5X, I3, 5X))') (II, II = IStart, IEnd)
         WRITE(*, '(4X, 10(4X, I3, 5X))') (II, II = IStart, IEnd)
         DO J = 1, N
!            WRITE(*, '(1X, I3, 10(1X, F12.6))') J, (XX(J,K), K = IStart, IEnd)
            WRITE(*, '(1X, I3, 10(1X, F11.5))') J, (XX(J,K), K = IStart, IEnd)
         END DO
         WRITE(*, *)
      END DO
!
      END SUBROUTINE
