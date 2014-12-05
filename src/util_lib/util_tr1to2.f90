      SUBROUTINE Util_Tr1To2(XX1, XX2, N)
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
         DO J = 1, I - 1
            IJ = IJ + 1
            XX2(J,I) = XX1(IJ)
            XX2(I,J) = XX2(J,I)
         END DO
         IJ = IJ + 1
         XX2(I,I) = XX1(IJ)
      END DO
!
      END SUBROUTINE
