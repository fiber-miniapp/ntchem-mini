      SUBROUTINE Util_InitMPI(Hybrid)
!
!$    USE OMP_LIB
      USE MPI_Module, ONLY : NThreads, NProcs, MyRank, NProcsIO, MyRankIO, MainRank, IORank, &
     &   MPIMain, MPIIO, NCorePerIO, MPI_COMM_IO, MyGroupIO, NumGroupIO, MPI_COMM_Node, &
     &   NCorePerMat, MPI_COMM_MAT, NProcsMat, MyRankMat, MyGroupMat, NumGroupMat, &
     &   MPI_COMM_MO, NProcsMO, MyRankMO, MyGroupMO, NumGroupMO
!
!     o Initialize MPI
!
      IMPLICIT NONE
!
      INCLUDE "mpif.h"
!
#ifdef MPIINT8
#define MPI_INTEGER MPI_INTEGER8
#endif
!
      LOGICAL, INTENT(IN) :: Hybrid
!
      INTEGER, PARAMETER :: IO1 = 5
      INTEGER :: MyColor, MyKey
      INTEGER :: IErr
      INTEGER :: MyRankNode, NProcsNode
!
      NAMELIST /Parallel/ NCorePerIO, NCorePerMat
!
      CALL MPI_INIT(IErr)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NProcs, IErr)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD, MyRank, IErr)
!
      MPIMain = (MyRank == MainRank)
      NCorePerIO = NProcs
      NCorePerMat = 1
!
!     o Read NCorePerIO from INPUT
!
      IF (MPIMain) THEN
!
         OPEN(UNIT=IO1, FILE='INPUT', STATUS='OLD', ACCESS='SEQUENTIAL', ERR=100)
         GO TO 200
  100    CONTINUE
         CALL Util_AbortMPI('Error in opening input file')
  200    CONTINUE
!
         REWIND(IO1)
         READ(IO1, Parallel, END=300, ERR=300)
         GO TO 400
  300    CONTINUE
         CALL Util_AbortMPI('Error: Check NAMELIST PARALLEL')
  400    CONTINUE
         CLOSE(IO1)
!
      END IF
!
      CALL MPI_BCAST(NCorePerIO, 1, MPI_INTEGER, MainRank, MPI_COMM_WORLD, IErr)
      CALL MPI_BCAST(NCorePerMat, 1, MPI_INTEGER, MainRank, MPI_COMM_WORLD, IErr)
!
!     o Local Communicator for file I/O
!       Definition of MyColor and MyKey may be dependent on hardware
!
      MyColor = MyRank / NCorePerIO
      MyKey = MOD(MyRank, NCorePerIO)
!
      CALL MPI_COMM_SPLIT(MPI_COMM_WORLD, MyColor, MyKey, MPI_COMM_IO, IErr)
      CALL MPI_COMM_SIZE(MPI_COMM_IO, NProcsIO, IErr)
      CALL MPI_COMM_RANK(MPI_COMM_IO, MyRankIO, IErr)
      MyGroupIO = MyColor
      NumGroupIO = NProcs / NCorePerIO
      IF (MOD(NProcs, NCorePerIO) /= 0) NumGroupIO = NumGroupIO + 1
!
!     o Assign I/O rank
!       MainRank needs to be the IORank of its IO communicator
!
      MPIIO = (MyRankIO == IORank)
!
      IF (MPIMain .AND. (.NOT. MPIIO)) THEN
         CALL Util_AbortMPI('Error: MainRank must be IORank')
      END IF
!
!     o Global Comminicator involving only IORank
!
      MyColor = 0
      IF (MPIIO) MyColor = 1
      MyKey = MyRank / NCorePerIO   ! May depend on hardware
!
      CALL MPI_COMM_SPLIT(MPI_COMM_WORLD, MyColor, MyKey, MPI_COMM_Node, IErr)
      IF (MPIIO) THEN
         CALL MPI_COMM_SIZE(MPI_COMM_Node, NProcsNode, IErr)
         CALL MPI_COMM_RANK(MPI_COMM_Node, MyRankNode, IErr)
      END IF
!
!     o Local Communicator for matrix multiplication
!       Definition of MyColor and MyKey may be dependent on hardware
!
      MyColor = MyRank / NCorePerMat
      MyKey = MOD(MyRank, NCorePerMat)
!
      CALL MPI_COMM_SPLIT(MPI_COMM_WORLD, MyColor, MyKey, MPI_COMM_MAT, IErr)
      CALL MPI_COMM_SIZE(MPI_COMM_MAT, NProcsMat, IErr)
      CALL MPI_COMM_RANK(MPI_COMM_MAT, MyRankMat, IErr)
      MyGroupMat = MyColor
      NumGroupMat = NProcs / NCorePerMat
      IF (MOD(NProcs, NCorePerMat) /= 0) NumGroupMat = NumGroupMat + 1
!
!     o Local Comminicator for orbital distribution
!
      MyColor = MOD(MyRank, NCorePerMat)
      MyKey = MyRank / NCorePerMat
!
      CALL MPI_COMM_SPLIT(MPI_COMM_WORLD, MyColor, MyKey, MPI_COMM_MO, IErr)
      CALL MPI_COMM_SIZE(MPI_COMM_MO, NProcsMO, IErr)
      CALL MPI_COMM_RANK(MPI_COMM_MO, MyRankMO, IErr)
      MyGroupMO = MyColor
      NumGroupMO = NProcs / NCorePerMat
      IF (MOD(NProcs, NCorePerMat) /= 0) NumGroupMO = NumGroupMO + 1
!
      IF (MPIMain) THEN
         WRITE(*, '("MPI has been initialized successfully")')
         WRITE(*, '("The job is running on", I8, " processes")') NProcs
         WRITE(*, '("Number of local I/O groups: ", I7)') NumGroupIO
         WRITE(*, '("Number of local matrix groups: ", I7)') NumGroupMat
         WRITE(*, '("Number of local MO groups: ", I7)') NumGroupMO
      END IF
!
!     o MPI/OpenMP hybrid
!
      IF (Hybrid) THEN
!$       NThreads = OMP_GET_MAX_THREADS()
!$       CALL OMP_SET_DYNAMIC(.FALSE.)
!$       CALL OMP_SET_NUM_THREADS(NThreads)
!$       IF (MPIMain) THEN
!$          WRITE(*, '("Number of threads :", I3)') NThreads
!$       END IF
      END IF
!
      END SUBROUTINE
