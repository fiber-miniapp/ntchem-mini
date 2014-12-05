      SUBROUTINE RIMP2_MaxParam_MPI
!
      USE MP2_Module, ONLY : NBF, NBF_Car, NBC, NBC_Car, MaxContS, Name
      USE RIMP2_Module, ONLY : NBF_RI, NBC_RI
      USE MP2_Basis_Module, ONLY : Spherical, GTOType
      USE RIMP2_Parameter_Module, ONLY : MaxSgmt_RI, MaxCont_RI, Maxtuv_RI, MaxShel_RI, MaxPrim_RI, MaxAngl_RI
      USE MP2_Parameter_Module, ONLY : MaxCont, MaxCont_Sph, MaxCont_Car, MaxConC, MaxConC_Sph, MaxConC_Car, MaxSgmt, MaxSgm2, &
     &   MaxPrim, MaxAngl, MaxAng2, MaxShel, MaxAtom, Maxtuv, MaxLtuv
      USE Int2_Gamma_Module, ONLY : MaxtuvGam
      USE MPI_Module, ONLY : MPIIO, IORank, MPI_COMM_IO
!
!     o Determine RIC (RIMP2) parameters
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
      INTEGER :: IAtom, IShel, ICont, IAngl, NShel0, NCont0, MCont_Car, MCont_Sph
      INTEGER :: IErr
!
!     o Open Basis_RIC file and read Basis_RIC information
!
      IF (MPIIO) THEN
!
         FBuf = TRIM(Name)//".Basis_RIC"
         OPEN(UNIT=IO, FILE=TRIM(FBuf), STATUS='OLD', ACCESS='SEQUENTIAL', FORM='FORMATTED', ERR=9000)
         REWIND(IO)
         GO TO 9010
 9000    REWIND(IO)
         CLOSE(IO)
         RETURN
 9010    CONTINUE
!
         READ(IO, *)
         READ(IO, *)
         READ(IO, *)
         READ(IO, *)
!
         MaxShel_RI = 0
         MaxPrim_RI = 0
         MaxAngl_RI = 0
         MaxSgmt_RI = 1
         MCont_Car = 0
         MCont_Sph = 0
         DO IAtom = 1, MaxAtom
            READ(IO, *)
            READ(IO, *) NShel0
            DO IShel = 1, NShel0
               MaxShel_RI = MaxShel_RI + 1
               READ(IO, *) CAngl, NCont0
               MaxSgmt_RI = MAX(MaxSgmt_RI, NCont0)
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
               MaxAngl_RI = MAX(MaxAngl_RI, IAngl)
               MCont_Car = MCont_Car + ((IAngl + 1) * (IAngl + 2)) / 2
               MCont_Sph = MCont_Sph + IAngl + IAngl + 1
               DO ICont = 1, NCont0
                  MaxPrim_RI = MaxPrim_RI + 1
                  READ(IO, *)
               END DO
!
            END DO
         END DO
!
         CLOSE(IO)
!
         IF (Spherical) THEN
            MaxCont_RI = MCont_Sph
         ELSE
            MaxCont_RI = MCont_Car
         END IF
!
      END IF
!
      CALL MPI_Bcast(MaxSgmt_RI, 1, MPI_INTEGER, IORank, MPI_COMM_IO, IErr)
      CALL MPI_Bcast(MaxCont_RI, 1, MPI_INTEGER, IORank, MPI_COMM_IO, IErr)
      CALL MPI_Bcast(MaxShel_RI, 1, MPI_INTEGER, IORank, MPI_COMM_IO, IErr)
      CALL MPI_Bcast(MaxPrim_RI, 1, MPI_INTEGER, IORank, MPI_COMM_IO, IErr)
      CALL MPI_Bcast(MaxAngl_RI, 1, MPI_INTEGER, IORank, MPI_COMM_IO, IErr)
!
      NBF_RI = MaxCont_RI
      NBC_RI = (NBF_RI * (NBF_RI + 1)) / 2
!
!     o Basis de kimeta parameter wo kokode saidai no hou ni okikaeru.
!
      MaxAngl = MAX(MaxAngl, MaxAngl_RI)
      MaxAng2 = MaxAngl * 2
      Maxtuv = MaxAngl * 4
      MaxLtuv = ((Maxtuv + 1) * (Maxtuv + 2) * (Maxtuv + 3)) / 6
      MaxtuvGam = Maxtuv + 7
      MaxSgmt = MAX(MaxSgmt, MaxSgmt_RI)
      MaxSgm2 = MaxSgmt * MaxSgmt
!
      IF (Spherical) THEN
         MaxContS = MaxAngl + MaxAngl + 1
      ELSE
         MaxContS = ((MaxAngl + 1) * (MaxAngl + 2)) / 2
      END IF
!
      Maxtuv_RI = Maxtuv
!
      END SUBROUTINE
