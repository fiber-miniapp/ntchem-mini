      SUBROUTINE Int2_Array_Deallocate
!
      USE Int2_Array_Module, ONLY : ERI_Array, ERI_Array_Temp, IJKLG, JInt
!
      IMPLICIT NONE
!
      DEALLOCATE(ERI_Array)
      DEALLOCATE(ERI_Array_Temp)
      DEALLOCATE(IJKLG)

      DEALLOCATE(JInt)  ! test
!
      END SUBROUTINE
