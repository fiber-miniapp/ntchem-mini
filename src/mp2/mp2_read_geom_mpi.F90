      SUBROUTINE MP2_Read_Geom_MPI
!
      USE MP2_Module, ONLY : Name
      USE MP2_Basis_Module, ONLY : Centr, NAtom
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
      INTEGER, PARAMETER :: IO = 99
      CHARACTER(LEN=255) :: FBuf
      CHARACTER(LEN=6) :: CAtm
      INTEGER :: IAtom
      INTEGER :: IErr
      REAL(8) :: ZAtm, ZAtmNo
!
!     o Open Geom file and read Geom information
!
      IF (MPIIO) THEN
         FBuf = TRIM(Name)//".Geom"
         OPEN(UNIT=IO, FILE=TRIM(FBuf), STATUS='OLD', ACCESS='SEQUENTIAL', FORM='FORMATTED')
         REWIND(IO)
         READ(IO, *) NAtom
         DO IAtom = 1, NAtom
            READ(IO, *) CAtm, Centr(1:3,IAtom), ZAtm, ZAtmNo
         END DO
         CLOSE(IO)
      END IF
!
      CALL MPI_Bcast(NAtom, 1, MPI_INTEGER, IORank, MPI_COMM_IO, IErr)
      CALL MPI_Bcast(Centr, 3*NAtom, MPI_DOUBLE_PRECISION, IORank, MPI_COMM_IO, IErr)
!
      END SUBROUTINE
