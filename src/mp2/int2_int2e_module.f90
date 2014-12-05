      MODULE Int2_Int2e_Module
!
!     o R-integrals for two-electron integrals
!
      REAL(8), ALLOCATABLE :: R0tuv(:,:,:)
      REAL(8), ALLOCATABLE :: ExpntPQ1(:,:)
      REAL(8), ALLOCATABLE :: ExpntPQ2(:,:)
      REAL(8), ALLOCATABLE :: PQX(:,:)
      REAL(8), ALLOCATABLE :: PQY(:,:)
      REAL(8), ALLOCATABLE :: PQZ(:,:)
!
!$OMP THREADPRIVATE(R0tuv, ExpntPQ1, ExpntPQ2, PQX, PQY, PQZ)
!
      END MODULE
