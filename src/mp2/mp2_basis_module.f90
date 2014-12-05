      MODULE MP2_Basis_Module
!
      CHARACTER(LEN=6), ALLOCATABLE :: CAtom(:)
      REAL(8), ALLOCATABLE :: Centr(:,:)
      REAL(8), ALLOCATABLE :: ZAtom(:)
      REAL(8), ALLOCATABLE :: Expnt(:)
      REAL(8), ALLOCATABLE :: CCoef(:)
      REAL(8), ALLOCATABLE :: ZAtomNo(:)   ! Atomic numbers
      INTEGER, ALLOCATABLE :: KStart(:)
      INTEGER, ALLOCATABLE :: KAtom(:)
      INTEGER, ALLOCATABLE :: KType(:)
      INTEGER, ALLOCATABLE :: KontG(:)
      INTEGER, ALLOCATABLE :: KLoc_Car(:)
      INTEGER, ALLOCATABLE :: KLoc_Sph(:)
      INTEGER, ALLOCATABLE :: KMin(:)   ! For DRK and PH
      INTEGER, ALLOCATABLE :: KMax(:)   ! For DRK and PH
      INTEGER, ALLOCATABLE :: Lt(:)
      INTEGER, ALLOCATABLE :: Lu(:)
      INTEGER, ALLOCATABLE :: Lv(:)
      INTEGER, ALLOCATABLE :: IndLtuv(:,:,:)
      INTEGER, ALLOCATABLE :: NLtuv_Car(:)
      INTEGER, ALLOCATABLE :: NLtuv_Sph(:)
      INTEGER, ALLOCATABLE :: LtuvMin_Car(:)
      INTEGER, ALLOCATABLE :: LtuvMin_Sph(:)
      INTEGER, ALLOCATABLE :: LtuvMax_Car(:)
      INTEGER, ALLOCATABLE :: LtuvMax_Sph(:)
      REAL(8), ALLOCATABLE :: FacDen(:,:,:)
      REAL(8), ALLOCATABLE :: SphCoef(:,:,:)
      LOGICAL :: NormP
      LOGICAL :: NormF
!     LLShl(IAtom): Lower shell limit of atom IAtom.
      INTEGER, ALLOCATABLE :: LLShl(:)
!     UlShl(IAtom): Upper shell limit of atom IAtom.
      INTEGER, ALLOCATABLE :: ULShl(:)
!     LLAO(IAtom)   : Lower AO limit of atom IAtom.
      INTEGER, ALLOCATABLE :: LLAO_Car(:)
      INTEGER, ALLOCATABLE :: LLAO(:)
!     ULAO(IAtom)   : Upper AO limit of atom IAtom.
      INTEGER, ALLOCATABLE :: ULAO_Car(:)
      INTEGER, ALLOCATABLE :: ULAO(:)
!     LLChi(IShl) : Lower AO limit from IShl.
      INTEGER, ALLOCATABLE :: LLChi_Car(:)
!     ULChi(IShl) : Upper AO limit of shell IShl.
      INTEGER, ALLOCATABLE :: ULChi_Car(:)
!
      CHARACTER(LEN=9) :: GTOType
      LOGICAL :: Spherical
!
!     NAtom     : Number of atoms.
      INTEGER :: NAtom
!     NShel     : Number of shells.
      INTEGER :: NShel
!
!     o For recurence relation
!
      REAL(8), ALLOCATABLE :: Expnt_HRR(:)
      REAL(8), ALLOCATABLE :: CCoef_HRR(:)
      INTEGER, ALLOCATABLE :: KStart_HRR(:)
      INTEGER, ALLOCATABLE :: KAtom_HRR(:)
      INTEGER, ALLOCATABLE :: KType_HRR(:)
      INTEGER, ALLOCATABLE :: KontG_HRR(:)
      INTEGER, ALLOCATABLE :: KLoc_HRR_Car(:)
      INTEGER, ALLOCATABLE :: KLoc_HRR_Sph(:)
      INTEGER, ALLOCATABLE :: MapShel(:)
!
!     o For parallel DFT
!
      INTEGER :: MaxAtomBlock
!
      END MODULE
