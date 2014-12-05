      SUBROUTINE MP2_MaxParam_MPI
!
      USE MP2_Module, ONLY : NBF_Car, NBC_Car, MaxContS, Name
      USE MP2_Basis_Module, ONLY : Spherical, GTOType
      USE MP2_Parameter_Module, ONLY : MaxCont, MaxCont_Sph, MaxCont_Car, MaxConC, MaxConC_Sph, MaxConC_Car, MaxSgmt, MaxSgm2, &
     &   MaxPrim, MaxAngl, MaxAng2, MaxShel, MaxAtom, Maxtuv, MaxLtuv
      USE Int2_Gamma_Module, ONLY : MaxtuvGam
      USE MPI_Module, ONLY : MPIIO, IORank, MPI_COMM_IO
!
!     o Determine SCF parameters
!
      IMPLICIT NONE
!
      INCLUDE 'mpif.h'
!
#ifdef MPIINT8
#define MPI_INTEGER MPI_INTEGER8
#endif
!
      INTEGER, PARAMETER :: IO = 99
      CHARACTER(LEN=255) :: FBuf
      CHARACTER(LEN=1) :: CAngl
      INTEGER :: IAtom, IShel, ICont, IAngl, NShel0, NCont0
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
         READ(IO, *)   ! NormP
         READ(IO, *)   ! NormF
         READ(IO, *) MaxAtom
!
         MaxShel = 0
         MaxPrim = 0
         MaxAngl = 0
         MaxSgmt = 1
         MaxCont_Car = 0
         MaxCont_Sph = 0
         DO IAtom = 1, MaxAtom
            READ(IO, *)
            READ(IO, *) NShel0
            DO IShel = 1, NShel0
               MaxShel = MaxShel + 1
               READ(IO, *) CAngl, NCont0
               MaxSgmt = MAX(MaxSgmt, NCont0)
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
               MaxAngl = MAX(MaxAngl, IAngl)
               MaxCont_Car = MaxCont_Car + ((IAngl + 1) * (IAngl + 2)) / 2
               MaxCont_Sph = MaxCont_Sph + IAngl + IAngl + 1
               DO ICont = 1, NCont0
                  MaxPrim = MaxPrim + 1
                  READ(IO, *)
               END DO
            END DO
         END DO
!
         CLOSE(IO)
!
      END IF
!
      CALL MPI_Bcast(GTOType, 9, MPI_CHARACTER, IORank, MPI_COMM_IO, IErr)
      CALL MPI_Bcast(MaxCont_Sph, 1, MPI_INTEGER, IORank, MPI_COMM_IO, IErr)
      CALL MPI_Bcast(MaxCont_Car, 1, MPI_INTEGER, IORank, MPI_COMM_IO, IErr)
      CALL MPI_Bcast(MaxSgmt, 1, MPI_INTEGER, IORank, MPI_COMM_IO, IErr)
      CALL MPI_Bcast(MaxPrim, 1, MPI_INTEGER, IORank, MPI_COMM_IO, IErr)
      CALL MPI_Bcast(MaxAngl, 1, MPI_INTEGER, IORank, MPI_COMM_IO, IErr)
      CALL MPI_Bcast(MaxShel, 1, MPI_INTEGER, IORank, MPI_COMM_IO, IErr)
      CALL MPI_Bcast(MaxAtom, 1, MPI_INTEGER, IORank, MPI_COMM_IO, IErr)
!
      IF (TRIM(GTOType) == 'SPHERICAL') THEN
         Spherical = .TRUE.
      ELSE IF (TRIM(GTOType) == 'CARTESIAN') THEN
         Spherical = .FALSE.
      END IF
!
      MaxAng2 = MaxAngl * 2
!NTQC      Maxtuv = MaxAngl * 2 + 1
      Maxtuv = MaxAngl * 4   ! For 2e integrals
      MaxtuvGam = Maxtuv + 7
      MaxLtuv = ((Maxtuv + 1) * (Maxtuv + 2) * (Maxtuv + 3)) / 6
      MaxSgm2 = MaxSgmt * MaxSgmt
      MaxConC_Car = (MaxCont_Car * (MaxCont_Car + 1)) / 2
      MaxConC_Sph = (MaxCont_Sph * (MaxCont_Sph + 1)) / 2
!
      IF (Spherical) THEN
         MaxCont = MaxCont_Sph
         MaxConC = MaxConC_Sph
      ELSE
         MaxCont = MaxCont_Car
         MaxConC = MaxConC_Car
      END IF
      NBF_Car = MaxCont_Car
      NBC_Car = (NBF_Car * (NBF_Car + 1)) / 2
!
      IF (Spherical) THEN
         MaxContS = MaxAngl + MaxAngl + 1
      ELSE
         MaxContS = ((MaxAngl + 1) * (MaxAngl + 2)) / 2
      END IF
!
      END SUBROUTINE
