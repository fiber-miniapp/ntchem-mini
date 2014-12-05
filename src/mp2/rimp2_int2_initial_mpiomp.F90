      SUBROUTINE RIMP2_Int2_Initial_MPIOMP
!
      USE MP2_Module, ONLY : PScreen
      USE Int2_Gamma_Module, ONLY : FF
!
      IMPLICIT NONE
!
!     o Read NAMELIST RIInt2
!
      CALL RIMP2_RIInt2_Read_Input_MPI
!
!     o Generate incomplete gamma functions
!
      CALL Int2_Gamma_Allocate
      CALL Int2_IncFun
      DEALLOCATE(FF)
!
!     o Calculate two-center integrals for prescreening
!
      IF (PScreen) THEN
         CALL MP2_Int2_PScreen_Allocate
         CALL MP2_MDInt2_PScreen_MPIOMP
         CALL RIMP2_RIInt2_PScreen_Allocate
         CALL RIMP2_RIInt2_MDInt2_PScreen_MPIOMP
      END IF
!
      END SUBROUTINE
