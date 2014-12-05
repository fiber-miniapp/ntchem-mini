      SUBROUTINE Int2_Array_Allocate
!
      USE MP2_Parameter_Module, ONLY : MaxAngl, MaxSgm2, MaxLtuv
      USE Int2_Array_Module, ONLY : ERI_Array, ERI_Array_Temp, IJKLG, JInt
!
      IMPLICIT NONE
!
!     o ERI array (Cartesian, Spherical and work space)
!
      ALLOCATE(ERI_Array(1:((MaxAngl+1)*(MaxAngl+2)/2)**4))
      ALLOCATE(ERI_Array_Temp(1:((MaxAngl+1)*(MaxAngl+2)/2)**2*(MaxAngl+MaxAngl+1)**2))
!
!     o Cartesian ERI index array
!
      ALLOCATE(IJKLG(1:(MaxAngl+1)*(MaxAngl+2)/2,1:(MaxAngl+1)*(MaxAngl+2)/2, &
     &               1:(MaxAngl+1)*(MaxAngl+2)/2,1:(MaxAngl+1)*(MaxAngl+2)/2))
!
      ALLOCATE(JInt(MaxSgm2,0:MaxLtuv))  ! test
!
      END SUBROUTINE
