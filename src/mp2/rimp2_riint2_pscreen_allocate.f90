      SUBROUTINE RIMP2_RIInt2_PScreen_Allocate
!
      USE RIMP2_Module, ONLY : SchwInt_RI
      USE RIMP2_Parameter_Module, ONLY : MaxShel_RI
!
      IMPLICIT NONE
!
      ALLOCATE(SchwInt_RI(MaxShel_RI))
!
      END SUBROUTINE
