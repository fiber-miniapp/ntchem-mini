      SUBROUTINE MP2_Int2_PScreen_Allocate
!
      USE MP2_Module, ONLY : SchwInt
      USE MP2_Parameter_Module, ONLY : MaxShel
!
      IMPLICIT NONE
!
      ALLOCATE(SchwInt((MaxShel*(MaxShel+1))/2))
!
      END SUBROUTINE
