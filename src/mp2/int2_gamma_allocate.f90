      SUBROUTINE Int2_Gamma_Allocate
!
      USE Int2_Gamma_Module, ONLY : FF, FTab, MaxGam, MaxtuvGam
!
      IMPLICIT NONE
!
!     o Imcomplete gamma functions
!
!YA062312      ALLOCATE(FTab(MaxtuvGam,MaxGam))
!YA062312      ALLOCATE(FF(MaxtuvGam))
      ALLOCATE(FTab(MaxtuvGam+1,MaxGam))
      ALLOCATE(FF(MaxtuvGam+1))
!
      END SUBROUTINE
