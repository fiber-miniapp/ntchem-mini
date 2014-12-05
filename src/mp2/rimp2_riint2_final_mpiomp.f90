      SUBROUTINE RIMP2_RIInt2_Final_MPIOMP
!
      USE MP2_Module, ONLY : PScreen, SchwInt
      USE RIMP2_Module, ONLY : SchwInt_RI
      USE Int2_Module, ONLY : DoMD4, DoMD2
      USE Int2_Gamma_Module, ONLY : FTab
!
      IMPLICIT NONE
!
!     o Memory deallocation
!
      IF (DoMD4 .OR. DoMD2) THEN
         DEALLOCATE(FTab)
      END IF
      IF (PScreen) THEN
         DEALLOCATE(SchwInt)
         DEALLOCATE(SchwInt_RI)
      END IF
!
      END SUBROUTINE
