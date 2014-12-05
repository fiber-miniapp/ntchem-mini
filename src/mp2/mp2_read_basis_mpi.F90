      SUBROUTINE MP2_Read_Basis_MPI
!
      USE MP2_Module, ONLY : Name
      USE MP2_Basis_Module, ONLY : Expnt, CCoef, KStart, KAtom, KType, KontG, &
     &   KMin, KMax, KLoc_Car, KLoc_Sph, NLtuv_Car, NLtuv_Sph, LLShl, ULShl, LLAO_Car, ULAO_Car, LLAO, ULAO, &
     &   LLChi_Car, ULChi_Car, NAtom, NShel, GTOType, NormF, NormP, Spherical
      USE MPI_Module, ONLY : MPIIO, IORank, MPI_COMM_IO
!
      IMPLICIT NONE
!
      INCLUDE 'mpif.h'
!
#ifdef MPIINT8
#define MPI_INTEGER MPI_INTEGER8
#endif
!
      INTEGER, PARAMETER :: MinF(0:4) = (/ 1, 2, 5, 11, 21/)   ! For DRK
      INTEGER, PARAMETER :: MaxF(0:4) = (/ 1, 4, 10, 20, 35/)   ! For DRK
      INTEGER, PARAMETER :: IO = 99
      CHARACTER(LEN=255) :: FBuf
      CHARACTER(LEN=1) :: CAngl
      INTEGER :: IAtom, IShel, ICont, IPrim, IAngl, ILoc_Car, ILoc_Sph, JAtom, NShel0, NCont0
      INTEGER :: IErr
!
!     o Open Basis file and read Basis information
!
      IF (MPIIO) THEN
!
         FBuf = TRIM(Name)//".Basis"
         OPEN(UNIT=IO, FILE=TRIM(FBuf), STATUS='OLD', ACCESS='SEQUENTIAL', FORM='FORMATTED')
         REWIND(IO)
!
         READ(IO, *) GTOType
         READ(IO, *) NormP
         READ(IO, *) NormF
         READ(IO, *) NAtom
!  
         NShel = 0
         IPrim = 0
         ILoc_Car = 0
         ILoc_Sph = 0
         DO IAtom = 1, NAtom
            READ(IO, *)
            READ(IO, *) NShel0
            IF (IAtom == 1) THEN
               LLShl(1) = 1
               ULShl(1) = NShel0
            ELSE
               LLShl(IAtom) = ULShl(IAtom-1) + 1
               ULShl(IAtom) = ULShl(IAtom-1) + NShel0
            END IF
            DO IShel = 1, NShel0
               NShel = NShel + 1
               READ(IO, *) CAngl, NCont0
               IF (TRIM(CAngl) == "S") THEN
                  IAngl = 0
               ELSE IF (TRIM(CAngl) == "P") THEN
                  IAngl = 1
               ELSE IF (TRIM(CAngl) == "D") THEN
                  IAngl = 2
               ELSE IF (TRIM(CAngl) == "F") THEN
                  IAngl = 3
               ELSE IF (TRIM(CAngl) == "G") THEN
                  IAngl = 4
               ELSE IF (TRIM(CAngl) == "H") THEN
                  IAngl = 5
               ELSE IF (TRIM(CAngl) == "I") THEN
                  IAngl = 6
               ELSE IF (TRIM(CAngl) == "K") THEN
                  IAngl = 7
               ELSE IF (TRIM(CAngl) == "L") THEN
                  IAngl = 8
               ELSE
                  CALL Util_AbortMPI('Error: Check IAngl')
               END IF
               KStart(NShel) = IPrim + 1
               KAtom(NShel) = IAtom
               KType(NShel) = IAngl
               KontG(NShel) = NCont0
               KLoc_Car(NShel) = ILoc_Car
               ILoc_Car = ILoc_Car + NLtuv_Car(IAngl)
               KLoc_Sph(NShel) = ILoc_Sph
               ILoc_Sph = ILoc_Sph + NLtuv_Sph(IAngl)
               KMin(NShel) = MinF(IAngl)
               KMax(NShel) = MaxF(IAngl)
               LLChi_Car(NShel) = KLoc_Car(NShel) + 1
               ULChi_Car(NShel) = ILoc_Car
               DO ICont = 1, NCont0
                  IPrim = IPrim + 1
                  READ(IO, *) Expnt(IPrim), CCoef(IPrim)
               END DO
!
!              o Normalization of basis sets
!
               CALL MP2_NormBasis(IAngl, KStart(NShel), IPrim)
!
            END DO
         END DO
!
         CLOSE(IO)
!
      END IF
!
      CALL MPI_Bcast(NormP, 1, MPI_LOGICAL, IORank, MPI_COMM_IO, IErr)
      CALL MPI_Bcast(NormF, 1, MPI_LOGICAL, IORank, MPI_COMM_IO, IErr)
      CALL MPI_Bcast(NAtom, 1, MPI_INTEGER, IORank, MPI_COMM_IO, IErr)
      CALL MPI_Bcast(NShel, 1, MPI_INTEGER, IORank, MPI_COMM_IO, IErr)
      CALL MPI_Bcast(KAtom, NShel, MPI_INTEGER, IORank, MPI_COMM_IO, IErr)
      CALL MPI_Bcast(KType, NShel, MPI_INTEGER, IORank, MPI_COMM_IO, IErr)
      CALL MPI_Bcast(KStart, NShel, MPI_INTEGER, IORank, MPI_COMM_IO, IErr)
      CALL MPI_Bcast(KType, NShel, MPI_INTEGER, IORank, MPI_COMM_IO, IErr)
      CALL MPI_Bcast(KontG, NShel, MPI_INTEGER, IORank, MPI_COMM_IO, IErr)
      CALL MPI_Bcast(KMin, NShel, MPI_INTEGER, IORank, MPI_COMM_IO, IErr)
      CALL MPI_Bcast(KMax, NShel, MPI_INTEGER, IORank, MPI_COMM_IO, IErr)
      CALL MPI_Bcast(KLoc_Car, NShel, MPI_INTEGER, IORank, MPI_COMM_IO, IErr)
      CALL MPI_Bcast(KLoc_Sph, NShel, MPI_INTEGER, IORank, MPI_COMM_IO, IErr)
      CALL MPI_Bcast(LLChi_Car, NShel, MPI_INTEGER, IORank, MPI_COMM_IO, IErr)
      CALL MPI_Bcast(ULChi_Car, NShel, MPI_INTEGER, IORank, MPI_COMM_IO, IErr)
      CALL MPI_Bcast(LLShl, NAtom, MPI_INTEGER, IORank, MPI_COMM_IO, IErr)
      CALL MPI_Bcast(ULShl, NAtom, MPI_INTEGER, IORank, MPI_COMM_IO, IErr)
      CALL MPI_Bcast(IPrim, 1, MPI_INTEGER, IORank, MPI_COMM_IO, IErr)
      CALL MPI_Bcast(Expnt, IPrim, MPI_DOUBLE_PRECISION, IORank, MPI_COMM_IO, IErr)
      CALL MPI_Bcast(CCoef, IPrim, MPI_DOUBLE_PRECISION, IORank, MPI_COMM_IO, IErr)
!
      ULAO_Car(1) = 0
      JAtom = 0
      DO IShel = 1, NShel
         IAtom = KAtom(IShel)
         IAngl = KType(IShel)
         IF (IAtom == JAtom) THEN
            ULAO_Car(IAtom) = ULAO_Car(IAtom) + NLtuv_Car(IAngl)
         ELSE IF ((IAtom /= JAtom) .AND. (IAtom /= 1)) THEN
            JAtom = JAtom + 1
            LLAO_Car(IAtom) = ULAO_Car(IAtom-1) + 1
            ULAO_Car(IAtom) = ULAO_Car(IAtom-1) + NLtuv_Car(IAngl)
         ELSE IF ((IAtom /= JAtom) .AND. (IAtom == 1)) THEN
            JAtom = JAtom + 1
            LLAO_Car(1) = 1
            ULAO_Car(1) = ULAO_Car(1) + NLtuv_Car(IAngl)
         ELSE
            CALL Util_AbortMPI('Error: Stop in MP2_Read_Basis')
         END IF
      END DO
!
      IF (Spherical) THEN
         ULAO(1) = 0
         JAtom = 0
         DO IShel = 1, NShel
            IAtom = KAtom(IShel)
            IAngl = KType(IShel)
            IF (IAtom == JAtom) THEN
               ULAO(IAtom) = ULAO(IAtom) + NLtuv_Sph(IAngl)
            ELSE IF ((IAtom /= JAtom) .AND. (IAtom /= 1)) THEN
               JAtom = JAtom + 1
               LLAO(IAtom) = ULAO(IAtom-1) + 1
               ULAO(IAtom) = ULAO(IAtom-1) + NLtuv_Sph(IAngl)
            ELSE IF ((IAtom /= JAtom) .AND. (IAtom == 1)) THEN
               JAtom = JAtom + 1
               LLAO(1) = 1
               ULAO(1) = ULAO(1) + NLtuv_Sph(IAngl)
            ELSE
               CALL Util_AbortMPI('Error: Stop in MP2_Read_Basis')
            END IF
         END DO
      ELSE
         LLAO(:) = LLAO_Car(:)
         ULAO(:) = ULAO_Car(:)
      END IF
!
      END SUBROUTINE
