      SUBROUTINE Int2_ECoef_Deallocate
!
      USE Int2_Module, ONLY : DoMD2
      USE Int2_ECoef_Module, ONLY : ECoefXAB, ECoefYAB, ECoefZAB, ECoefXCD, ECoefYCD, ECoefZCD, ExpPHalf, ExpQHalf, &
     &   ECoefXAB_MD2, ECoefYAB_MD2, ECoefZAB_MD2, ECoefXCD_MD2, ECoefYCD_MD2, ECoefZCD_MD2
!
      IMPLICIT NONE
!
!     o E-coefficients
!
      DEALLOCATE(ECoefXAB)
      DEALLOCATE(ECoefYAB)
      DEALLOCATE(ECoefZAB)
      DEALLOCATE(ECoefXCD)
      DEALLOCATE(ECoefYCD)
      DEALLOCATE(ECoefZCD)
      DEALLOCATE(ExpPHalf)
      DEALLOCATE(ExpQHalf)
      IF (DoMD2)  THEN
         DEALLOCATE(ECoefXAB_MD2)
         DEALLOCATE(ECoefYAB_MD2)
         DEALLOCATE(ECoefZAB_MD2)
         DEALLOCATE(ECoefXCD_MD2)
         DEALLOCATE(ECoefYCD_MD2)
         DEALLOCATE(ECoefZCD_MD2)
      END IF
!
      END SUBROUTINE
