      SUBROUTINE Util_LinOut(XX1, N)
!
      IMPLICIT NONE
!
      INTEGER :: N
      REAL(8) :: XX1(*)
!
      INTEGER :: I, J, IJ
      REAL(8) :: XX2(N,N)
!
      IJ = 0
      DO I = 1, N
         DO J = 1, I
            IJ = IJ + 1
            XX2(I,J) = XX1(IJ)
            XX2(J,I) = XX2(I,J)
         END DO
      END DO
!
      CALL Util_MatOut(XX2, N, N)
!
      END SUBROUTINE
