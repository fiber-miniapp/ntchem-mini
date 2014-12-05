      MODULE RIMP2_Basis_Module
!
      INTEGER :: NShel_RI
!
      INTEGER, ALLOCATABLE :: KStart_RI(:)
      INTEGER, ALLOCATABLE :: KAtom_RI(:)
      INTEGER, ALLOCATABLE :: KType_RI(:)
      INTEGER, ALLOCATABLE :: KontG_RI(:)
      INTEGER, ALLOCATABLE :: KLoc_RI_Car(:)
      INTEGER, ALLOCATABLE :: KLoc_RI_Sph(:)
      INTEGER, ALLOCATABLE :: KMin_RI(:)   ! For DRK
      INTEGER, ALLOCATABLE :: KMax_RI(:)   ! For DRK
      INTEGER, ALLOCATABLE :: NLtuv_RI_Car(:)
      INTEGER, ALLOCATABLE :: NLtuv_RI_Sph(:)
      REAL(8), ALLOCATABLE :: Expnt_RI(:)
      REAL(8), ALLOCATABLE :: CCoef_RI(:)
!
!     o For recurence relation
!
      INTEGER, ALLOCATABLE :: KStart_RI_HRR(:)
      INTEGER, ALLOCATABLE :: KAtom_RI_HRR(:)
      INTEGER, ALLOCATABLE :: KType_RI_HRR(:)
      INTEGER, ALLOCATABLE :: KontG_RI_HRR(:)
      INTEGER, ALLOCATABLE :: KLoc_RI_HRR(:)
      INTEGER, ALLOCATABLE :: KLoc_RI_HRR_Sph(:)
      INTEGER, ALLOCATABLE :: MapShel_RI(:)
      REAL(8), ALLOCATABLE :: Expnt_RI_HRR(:)
      REAL(8), ALLOCATABLE :: CCoef_RI_HRR(:)
!
      END MODULE
