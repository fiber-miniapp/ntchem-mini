      SUBROUTINE Int2_Int2e_Allocate
!
      USE MP2_Parameter_Module
      USE Int2_Int2e_Module
!
      IMPLICIT NONE
!
!     o R-integrals for two-electron integrals
!
      ALLOCATE(R0tuv(MaxSgm2,MaxSgm2,0:MaxLtuv-1))
      ALLOCATE(ExpntPQ1(MaxSgm2,MaxSgm2))
      ALLOCATE(ExpntPQ2(MaxSgm2,MaxSgm2))
      ALLOCATE(PQX(MaxSgm2,MaxSgm2))
      ALLOCATE(PQY(MaxSgm2,MaxSgm2))
      ALLOCATE(PQZ(MaxSgm2,MaxSgm2))
!
      END SUBROUTINE
