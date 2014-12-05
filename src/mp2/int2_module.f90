      MODULE Int2_Module
!
!     o Module for McMurchie-Davidson and Dupuis-Rys-King integral schemes
!
      CHARACTER(LEN=8) :: IntType
      LOGICAL :: PHInt
      LOGICAL :: DoMD4
      LOGICAL :: DoMD2
      LOGICAL :: DoDRK2
      LOGICAL :: DoLibint
      LOGICAL :: DoPH
      LOGICAL :: Only1c
      LOGICAL :: Only2c
      LOGICAL :: NDDO
      LOGICAL :: FMM_ShortRange = .TRUE.
      LOGICAL :: FMM_LongRange = .TRUE.
      LOGICAL :: LibInit = .FALSE.
      LOGICAL :: Int2Init = .FALSE.
      LOGICAL :: Int2Fin = .FALSE.
      INTEGER :: IAnglA, IAnglB, IAnglC, IAnglD
      INTEGER :: IAnglAB, IAnglCD
      INTEGER :: IAtom1, IAtom2, IAtom3, IAtom4
      INTEGER :: LibVer
      REAL(8) :: ThrInt
      REAL(8) :: ThrPrim
      REAL(8), ALLOCATABLE :: PreFactAB(:), PreFactCD(:)
      REAL(8), ALLOCATABLE :: CCoefAB(:), CCoefCD(:)
      REAL(8), ALLOCATABLE :: PX(:), PY(:), PZ(:)
      REAL(8), ALLOCATABLE :: PAX(:), PAY(:), PAZ(:)
      REAL(8), ALLOCATABLE :: PBX(:), PBY(:), PBZ(:)
      REAL(8), ALLOCATABLE :: QX(:), QY(:), QZ(:)
      REAL(8), ALLOCATABLE :: QCX(:), QCY(:), QCZ(:)
      REAL(8), ALLOCATABLE :: QDX(:), QDY(:), QDZ(:)
      REAL(8), ALLOCATABLE :: CContAB(:,:,:,:,:,:,:)
      REAL(8), ALLOCATABLE :: CContCD(:,:,:,:,:,:,:)
      REAL(8), ALLOCATABLE :: ExpntA(:), ExpntB(:), ExpntC(:), ExpntD(:)
      REAL(8), ALLOCATABLE :: ExpntP(:), ExpntQ(:)
!
!     o OMP private variables
!
!$OMP THREADPRIVATE(IAnglA, IAnglB, IAnglC, IAnglD, IAnglAB, IAnglCD, &
!$OMP&   PreFactAB, PreFactCD, CCoefAB, CCoefCD, &
!$OMP&   PX, PY, PZ, PAX, PAY, PAZ, PBX, PBY, PBZ, &
!$OMP&   QX, QY, QZ, QCX, QCY, QCZ, QDX, QDY, QDZ, &
!$OMP&   CContAB, CContCD, ExpntA, ExpntB, ExpntC, ExpntD, ExpntP, ExpntQ)
!
      END MODULE
