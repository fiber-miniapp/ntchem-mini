      MODULE Int2_Gamma_Module
!
!     o Imcomplete gamma functions
!
      INTEGER, PARAMETER :: MaxGam = 1400
      INTEGER :: MaxtuvGam
      REAL(8), ALLOCATABLE :: FTab(:,:)
      REAL(8), ALLOCATABLE :: FF(:)
!
!$OMP THREADPRIVATE(FF)
!
      END MODULE
