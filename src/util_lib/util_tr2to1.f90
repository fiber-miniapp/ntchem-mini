      SUBROUTINE Util_Tr2To1(XX2, XX1, N)
!
      IMPLICIT NONE
!
      INTEGER :: N
      REAL(8) :: XX1(*), XX2(N,*)
!
      INTEGER :: I, J, IJ
!
      IJ = 0
      DO I = 1, N
         DO J = 1, I
            IJ = IJ + 1
            XX1(IJ) = XX2(J,I)
         END DO
      END DO
!
      END SUBROUTINE
