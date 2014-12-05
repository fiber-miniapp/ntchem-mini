      SUBROUTINE Util_FinMPI
!
      USE MPI_Module, ONLY : MPIMain
!
!     o Terminate MPI
!
      IMPLICIT NONE
!
      INCLUDE "mpif.h"
!
      INTEGER :: IErr
!
      CALL MPI_FINALIZE(IErr)
!
      IF (MPIMain) THEN
         WRITE(*, '("MPI has been terminated")')
      END IF
!
      END SUBROUTINE
