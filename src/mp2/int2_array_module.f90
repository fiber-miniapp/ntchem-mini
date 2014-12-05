      MODULE Int2_Array_Module
!
!    o Array of Cartesian ERI index
!
      REAL(8), ALLOCATABLE ::ERI_Array(:)
      REAL(8), ALLOCATABLE ::ERI_Array_Temp(:)
      INTEGER, ALLOCATABLE :: IJKLG(:,:,:,:)

      REAL(8), ALLOCATABLE :: JInt(:,:)  ! test
!
!$OMP THREADPRIVATE(ERI_Array, ERI_Array_Temp, IJKLG, JInt)
      END MODULE
