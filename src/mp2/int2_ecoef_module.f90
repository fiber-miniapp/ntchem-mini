      MODULE Int2_ECoef_Module
!
!     o E-coefficients
!
      REAL(8), ALLOCATABLE :: ECoefXAB(:,:,:,:)
      REAL(8), ALLOCATABLE :: ECoefYAB(:,:,:,:)
      REAL(8), ALLOCATABLE :: ECoefZAB(:,:,:,:)
      REAL(8), ALLOCATABLE :: ECoefXCD(:,:,:,:)
      REAL(8), ALLOCATABLE :: ECoefYCD(:,:,:,:)
      REAL(8), ALLOCATABLE :: ECoefZCD(:,:,:,:)
      REAL(8), ALLOCATABLE :: ExpPHalf(:)
      REAL(8), ALLOCATABLE :: ExpQHalf(:)
      REAL(8), ALLOCATABLE :: ECoefXAB_MD2(:,:,:)
      REAL(8), ALLOCATABLE :: ECoefYAB_MD2(:,:,:)
      REAL(8), ALLOCATABLE :: ECoefZAB_MD2(:,:,:)
      REAL(8), ALLOCATABLE :: ECoefXCD_MD2(:,:,:)
      REAL(8), ALLOCATABLE :: ECoefYCD_MD2(:,:,:)
      REAL(8), ALLOCATABLE :: ECoefZCD_MD2(:,:,:)
!
!$OMP THREADPRIVATE(ECoefXAB, ECoefYAB, ECoefZAB, ECoefXCD, ECoefYCD, ECoefZCD, &
!$OMP&   ExpPHalf, ExpQHalf, &
!$OMP&   ECoefXAB_MD2, ECoefYAB_MD2, ECoefZAB_MD2, &
!$OMP&   ECoefXCD_MD2, ECoefYCD_MD2, ECoefZCD_MD2)
!
      END MODULE
