      SUBROUTINE MP2_IOff
!
      USE MP2_Module, ONLY : DoRI, IOff, NBF_Car
      USE RIMP2_Module, ONLY : NBF_RI
!
      IMPLICIT NONE
!
      INTEGER :: IBF, MBF
!
      MBF = NBF_Car
      IF (DoRI) THEN
         MBF = MAX(NBF_Car, NBF_RI)
      END IF
!
      DO IBF = 1, MBF
         IOff(IBF) = (IBF * (IBF - 1)) / 2
      END DO
!
      END SUBROUTINE
