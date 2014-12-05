      SUBROUTINE Util_AbortMPI(Msg)
!
      USE MPI_Module, ONLY : MPIMain
!
!     o Error temination of MPI
!
      IMPLICIT NONE
!
      INCLUDE "mpif.h"
!
      CHARACTER(LEN=255) :: Msg
      INTEGER :: IErr, IFailure
!
      IF (MPIMain) THEN
         WRITE(*, *) TRIM(Msg)
      END IF
      CALL MPI_ABORT(MPI_COMM_WORLD, IFailure, IErr)
      STOP
!
      END SUBROUTINE
