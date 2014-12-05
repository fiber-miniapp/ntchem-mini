      MODULE MPI_Module
!
!     o MPI information
!
      LOGICAL :: MPIMain
      LOGICAL :: MPIIO
      INTEGER(4) :: NThreads
      INTEGER, PARAMETER :: MainRank = 0
      INTEGER, PARAMETER :: IORank = 0
      INTEGER :: NCorePerIO   ! # of cores sharing a local disk
      INTEGER :: NCorePerMat  ! # of cores distributing a matrix
      INTEGER :: NProcs
      INTEGER :: MyRank
      INTEGER :: MPI_COMM_IO
      INTEGER :: MPI_COMM_Node
      INTEGER :: MPI_COMM_MAT
      INTEGER :: MPI_COMM_MO
      INTEGER :: NProcsIO
      INTEGER :: MyRankIO
      INTEGER :: MyGroupIO
      INTEGER :: NumGroupIO
      INTEGER :: NProcsMat
      INTEGER :: MyRankMat
      INTEGER :: MyGroupMat
      INTEGER :: NumGroupMat
      INTEGER :: NProcsMO
      INTEGER :: MyRankMO
      INTEGER :: MyGroupMO
      INTEGER :: NumGroupMO
!YA050712      INTEGER :: MPI_COMM_Node
!YA050712      INTEGER :: NProcsNode
!YA050712      INTEGER :: MyRankNode
!
      END MODULE
