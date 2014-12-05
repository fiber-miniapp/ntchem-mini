      MODULE RIMP2_Module
!
!     o Module for RI-MP2 calculation
!
      INTEGER :: NBF_RI
      INTEGER :: NBC_RI
      INTEGER, ALLOCATABLE :: NBF_RI_MyRank(:)
      INTEGER, ALLOCATABLE :: IdxBF_RI_MyRank(:)
      REAL(8), ALLOCATABLE :: RIInt2c(:)
      REAL(8), ALLOCATABLE :: SchwInt_RI(:)
      ! REAL(8), ALLOCATABLE, target :: RIInt3c2a(:,:), RIInt3c2b(:,:), RIInt3c3a(:,:), RIInt3c3b(:,:)
      REAL(8), ALLOCATABLE :: RIInt3c2b(:,:), RIInt3c3b(:,:)
#ifdef USE_GPU
      REAL(8), ALLOCATABLE, target, pinned :: RIInt3c3a(:,:)
      REAL(8), ALLOCATABLE, pinned :: RI2cInv(:,:)
      REAL(8), ALLOCATABLE, pinned :: RIInt3c2a(:,:)
#else
      REAL(8), ALLOCATABLE, target :: RIInt3c3a(:,:)
      REAL(8), ALLOCATABLE :: RI2cInv(:,:)
      REAL(8), ALLOCATABLE :: RIInt3c2a(:,:)
#endif
!
      END MODULE
