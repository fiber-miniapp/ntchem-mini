      SUBROUTINE Int2_ECoef_Allocate
!
      USE MP2_Parameter_Module, ONLY : MaxAngl, MaxAng2, MaxSgm2
      USE Int2_Module, ONLY : DoMD2
      USE Int2_ECoef_Module, ONLY : ECoefXAB, ECoefYAB, ECoefZAB, ECoefXCD, ECoefYCD, ECoefZCD, ExpPHalf, ExpQHalf, &
     &   ECoefXAB_MD2, ECoefYAB_MD2, ECoefZAB_MD2, ECoefXCD_MD2, ECoefYCD_MD2, ECoefZCD_MD2
!
      IMPLICIT NONE
!
!     o E-coefficients
!
      ALLOCATE(ECoefXAB(MaxSgm2,0:MaxAngl,0:MaxAngl,0:MaxAng2))
      ALLOCATE(ECoefYAB(MaxSgm2,0:MaxAngl,0:MaxAngl,0:MaxAng2))
      ALLOCATE(ECoefZAB(MaxSgm2,0:MaxAngl,0:MaxAngl,0:MaxAng2))
      ALLOCATE(ECoefXCD(MaxSgm2,0:MaxAngl,0:MaxAngl,0:MaxAng2))
      ALLOCATE(ECoefYCD(MaxSgm2,0:MaxAngl,0:MaxAngl,0:MaxAng2))
      ALLOCATE(ECoefZCD(MaxSgm2,0:MaxAngl,0:MaxAngl,0:MaxAng2))
      ALLOCATE(ExpPHalf(MaxSgm2))
      ALLOCATE(ExpQHalf(MaxSgm2))
      IF (DoMD2)  THEN
         ALLOCATE(ECoefXAB_MD2(MaxSgm2,0:MaxAng2,0:MaxAng2))
         ALLOCATE(ECoefYAB_MD2(MaxSgm2,0:MaxAng2,0:MaxAng2))
         ALLOCATE(ECoefZAB_MD2(MaxSgm2,0:MaxAng2,0:MaxAng2))
         ALLOCATE(ECoefXCD_MD2(MaxSgm2,0:MaxAng2,0:MaxAng2))
         ALLOCATE(ECoefYCD_MD2(MaxSgm2,0:MaxAng2,0:MaxAng2))
         ALLOCATE(ECoefZCD_MD2(MaxSgm2,0:MaxAng2,0:MaxAng2))
      END IF
!
      END SUBROUTINE
