      SUBROUTINE Int2_Gamma_Deallocate
!
      USE Int2_Gamma_Module, ONLY : FF, FTab
!
      IMPLICIT NONE
!
!     o Imcomplete gamma functions
!
      DEALLOCATE(FTab)
      DEALLOCATE(FF)
!
      END SUBROUTINE
