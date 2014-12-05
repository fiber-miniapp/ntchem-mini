      SUBROUTINE Int2_Int2e_Deallocate
!
      USE Int2_Int2e_Module
!
      IMPLICIT NONE
!
!     o R-integrals for two-electron integrals
!
      DEALLOCATE(R0tuv)
      DEALLOCATE(ExpntPQ1)
      DEALLOCATE(ExpntPQ2)
      DEALLOCATE(PQX)
      DEALLOCATE(PQY)
      DEALLOCATE(PQZ)
!
      END SUBROUTINE
