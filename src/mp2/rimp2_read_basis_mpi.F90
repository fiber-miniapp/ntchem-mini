      SUBROUTINE RIMP2_Read_Basis_MPI
!
      USE MP2_Module, ONLY : Name
      USE MP2_Basis_Module, ONLY : NAtom
      USE RIMP2_Basis_Module, ONLY : KMin_RI, KMax_RI, NLtuv_RI_Car, KLoc_RI_Car, NLtuv_RI_Sph, KLoc_RI_Sph, &
     &   CCoef_RI, Expnt_RI, KAtom_RI, KStart_RI, NShel_RI, KontG_RI, KType_RI
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
      INTEGER, PARAMETER :: MinF(0:4) = (/ 1,  2,  5, 11, 21/)   ! For DRK
      INTEGER, PARAMETER :: MaxF(0:4) = (/ 1,  4, 10, 20, 35/)   ! For DRK
      INTEGER, PARAMETER :: IO = 99
      CHARACTER(LEN=255) :: FBuf
      CHARACTER(LEN=1) :: CAngl
      INTEGER :: IAtom, IShel, ICont, IPrim, IAngl, ILoc_Car, ILoc_Sph, NShel0, NCont0
      INTEGER :: IErr
!
!     o Open Basis_RIC file and read Basis_RIC information
!
      IF (MPIIO) THEN
!
         FBuf = TRIM(Name)//".Basis_RIC"
         OPEN(UNIT=IO, FILE=TRIM(FBuf), STATUS='OLD', ACCESS='SEQUENTIAL', FORM='FORMATTED')
         REWIND(IO)
!
         READ(IO, *)
         READ(IO, *)
         READ(IO, *)
         READ(IO, *)
!
         NShel_RI = 0
         IPrim = 0
         ILoc_Car = 0
         ILoc_Sph = 0
         DO IAtom = 1, NAtom
            READ(IO, *)
            READ(IO, *) NShel0
            DO IShel = 1, NShel0
               NShel_RI = NShel_RI + 1
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
               KStart_RI(NShel_RI) = IPrim + 1
               KAtom_RI(NShel_RI) = IAtom
               KType_RI(NShel_RI) = IAngl
               KontG_RI(NShel_RI) = NCont0
               KLoc_RI_Car(NShel_RI) = ILoc_Car
               ILoc_Car = ILoc_Car + NLtuv_RI_Car(IAngl)
               KLoc_RI_Sph(NShel_RI) = ILoc_Sph
               ILoc_Sph = ILoc_Sph + NLtuv_RI_Sph(IAngl)
               KMin_RI(NShel_RI) = MinF(IAngl)
               KMax_RI(NShel_RI) = MaxF(IAngl)
               DO ICont = 1, NCont0
                  IPrim = IPrim + 1
                  READ(IO, *) Expnt_RI(IPrim), CCoef_RI(IPrim)
               END DO
!
!              o Normalization of basis sets
!
               CALL RIMP2_NormBasis(IAngl, KStart_RI(NShel_RI), IPrim)
!
            END DO
         END DO
!
!        o Close Basis_RIC file
!
         CLOSE(IO)
!
      END IF
!
      CALL MPI_Bcast(NShel_RI, 1, MPI_INTEGER, IORank, MPI_COMM_IO, IErr)
      CALL MPI_Bcast(KMin_RI, NShel_RI, MPI_INTEGER, IORank, MPI_COMM_IO, IErr)
      CALL MPI_Bcast(KMax_RI, NShel_RI, MPI_INTEGER, IORank, MPI_COMM_IO, IErr)
      CALL MPI_Bcast(KLoc_RI_Car, NShel_RI, MPI_INTEGER, IORank, MPI_COMM_IO, IErr)
      CALL MPI_Bcast(KLoc_RI_Sph, NShel_RI, MPI_INTEGER, IORank, MPI_COMM_IO, IErr)
      CALL MPI_Bcast(KAtom_RI, NShel_RI, MPI_INTEGER, IORank, MPI_COMM_IO, IErr)
      CALL MPI_Bcast(KStart_RI, NShel_RI, MPI_INTEGER, IORank, MPI_COMM_IO, IErr)
      CALL MPI_Bcast(KontG_RI, NShel_RI, MPI_INTEGER, IORank, MPI_COMM_IO, IErr)
      CALL MPI_Bcast(KType_RI, NShel_RI, MPI_INTEGER, IORank, MPI_COMM_IO, IErr)
      CALL MPI_Bcast(IPrim, 1, MPI_INTEGER, IORank, MPI_COMM_IO, IErr)
      CALL MPI_Bcast(CCoef_RI, IPrim, MPI_DOUBLE_PRECISION, IORank, MPI_COMM_IO, IErr)
      CALL MPI_Bcast(Expnt_RI, IPrim, MPI_DOUBLE_PRECISION, IORank, MPI_COMM_IO, IErr)
!
      END SUBROUTINE
