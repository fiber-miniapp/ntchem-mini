      SUBROUTINE Util_LowTrMZero(XX2, NDim, N)
!
!     o Symmetrisation of a upper triangle matrix
!
      IMPLICIT NONE
!
      INTEGER :: N, NDim
      REAL(8) :: XX2(NDim,NDim)
!
      REAL(8), PARAMETER :: Zero = 0.0D+00
      INTEGER :: I, J
!
      DO I = 1, (N - 1)
        DO J = (I + 1), N
           XX2(J,I) = Zero
         END DO
      END DO
!
      END SUBROUTINE
