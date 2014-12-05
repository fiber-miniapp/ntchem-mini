      SUBROUTINE MP2_Read_Input_MPI
!
      USE MP2_Module, ONLY : MP2BatchLv, DoRI, InCore, NBF, NMO, NActOA, NActOB, NActVA, NActVB, &
     &   NFrzOA, NFrzOB, NFrzVA, NFrzVB, E1, Name, IPrint
      USE MP2_Constant_Module, ONLY : Zero
      USE MPI_Module, ONLY : MyRank, MPIIO, IORank, MPI_COMM_IO
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
      CHARACTER(LEN=6) :: MP2Type
      CHARACTER(LEN=9) :: OrthType
      INTEGER :: IErr
!
      NAMELIST /Prefix/ Name
      NAMELIST /MP2/ MP2Type, NFrzOA, NFrzOB, NFrzVA, NFrzVB, MP2BatchLv, InCore, IPrint
!
      Name = 'ntchem'
      MP2Type = 'RIMP2'
      OrthType = 'Cholesky'   ! 'Canonical'
      NFrzOA = 0
      NFrzOB = 0
      NFrzVA = 0
      NFrzVB = 0
      MP2BatchLv = 0
      IPrint = 0
      DoRI = .TRUE.
      InCore = .TRUE.
!
      E1 = Zero
!
      IF (MPIIO) THEN
!
         OPEN(UNIT=IO, FILE='INPUT', STATUS='OLD', ACCESS='SEQUENTIAL')
!
!        o Read NAMELIST Prefix
!
         REWIND(IO)
         READ(IO, Prefix, END=8000, ERR=8000)
         GO TO 8010
 8000    CONTINUE
 8010    CONTINUE
!
!        o Read NAMELIST MP2
!
         REWIND(IO)
         READ(IO, MP2, END=9000, ERR=9000)
         GO TO 9010

 9000    CONTINUE
         CALL Util_AbortMPI('Error: Check NAMELIST MP2')
 9010    CONTINUE
!
         CLOSE(IO)
!
      END IF
!
      CALL MPI_Bcast(Name, 255, MPI_CHARACTER, IORank, MPI_COMM_IO, IErr)
      CALL MPI_Bcast(MP2Type,  6, MPI_CHARACTER, IORank, MPI_COMM_IO, IErr)
      CALL MPI_Bcast(OrthType, 9, MPI_CHARACTER, IORank, MPI_COMM_IO, IErr)
      CALL MPI_Bcast(NFrzOA, 1, MPI_INTEGER, IORank, MPI_COMM_IO, IErr)
      CALL MPI_Bcast(NFrzOB, 1, MPI_INTEGER, IORank, MPI_COMM_IO, IErr)
      CALL MPI_Bcast(NFrzVA, 1, MPI_INTEGER, IORank, MPI_COMM_IO, IErr)
      CALL MPI_Bcast(NFrzVB, 1, MPI_INTEGER, IORank, MPI_COMM_IO, IErr)
      CALL MPI_Bcast(MP2BatchLv, 1, MPI_INTEGER, IORank, MPI_COMM_IO, IErr)
      CALL MPI_Bcast(InCore, 1, MPI_LOGICAL, IORank, MPI_COMM_IO, IErr)
      CALL MPI_Bcast(IPrint, 1, MPI_INTEGER, IORank, MPI_COMM_IO, IErr)
!
!     o Choose MP2Type
!
      MP2Type = 'RIMP2'
      CALL Util_TransChar(MP2Type, 6)
      IF ((TRIM(MP2Type) == 'RI') .OR. (TRIM(MP2Type) == 'RIMP2')) DoRI = .TRUE.
!
!     o Choose OrthType
!
      OrthType = 'CHOLESKY'
      CALL Util_TransChar(OrthType, 9)
      IF (TRIM(OrthType) == '') OrthType = 'CHOLESKY'
!
      IF (MyRank == 0) THEN
         WRITE(*, *) 'o MP2Type  = ', TRIM(MP2Type)
         WRITE(*, *) 'o OrthType = ', TRIM(OrthType)
         WRITE(*, *) 'o InCore   = ', InCore
         WRITE(*, *) 'o Name     = ', TRIM(Name)
      END IF
!
      END SUBROUTINE
