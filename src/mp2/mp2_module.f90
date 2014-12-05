      MODULE MP2_Module
!
!     o Module for MO transformation and MP2 calculation
!
      CHARACTER(LEN=255) :: Name
      LOGICAL :: UHF
      LOGICAL :: DoRI
      LOGICAL :: InCore
      INTEGER :: MP2BatchLv
      INTEGER :: NOccBat
      INTEGER :: LenOccBat
      INTEGER :: NOccBat_per_Proc
      INTEGER :: NAlpBet
      INTEGER :: NBF
      INTEGER :: NBF_Car
      INTEGER :: NBC
      INTEGER :: NBC_Car
      INTEGER :: NMO
      INTEGER :: NActOA
      INTEGER :: NActOB
      INTEGER :: NActVA
      INTEGER :: NActVB
      INTEGER :: NFrzOA
      INTEGER :: NFrzOB
      INTEGER :: NFrzVA
      INTEGER :: NFrzVB
      INTEGER :: MaxContS
      INTEGER :: IPrint
      INTEGER :: MaxNActO
      INTEGER :: NActO(2)
      INTEGER :: NActV(2)
      INTEGER :: NFrzO(2)
      INTEGER :: NFrzV(2)
      REAL(8) :: EneMP2
      REAL(8) :: EMP2
      REAL(8) :: ESCSMP2
      REAL(8) :: E1
      REAL(8) :: E2T
      REAL(8) :: E2S
      REAL(8) :: E2
      REAL(8) :: E2SCS
      INTEGER, ALLOCATABLE :: IOff(:)
      INTEGER, ALLOCATABLE :: IOccBat(:,:,:)
!
!     o Module for the prescreening
!
      LOGICAL :: PScreen
      REAL(8) :: ThrPre
      REAL(8), ALLOCATABLE :: SchwInt(:)
!
      END MODULE MP2_Module
