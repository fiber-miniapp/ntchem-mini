      SUBROUTINE MP2_Allocate
!
      USE MP2_Module, ONLY : IOff, NBF_Car, DoRI
      USE RIMP2_Module, ONLY : NBF_RI
!
      IMPLICIT NONE
!
      IF (.NOT. DoRI) THEN
         ALLOCATE(IOff(NBF_Car))
      ELSE
         ALLOCATE(IOff(MAX(NBF_Car, NBF_RI)))
      END IF
!
      END SUBROUTINE
