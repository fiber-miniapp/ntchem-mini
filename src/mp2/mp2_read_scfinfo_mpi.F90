      SUBROUTINE MP2_Read_SCFInfo_MPI
!
      USE MP2_Module, ONLY : DoRI, UHF, NAlpBet, NBF, NBC, NMO, &
     &       NActOA, NActOB, NActVA, NActVB, NFrzOA, NFrzOB, NFrzVA, NFrzVB, &
     &       NActO, NActV, NFrzO, NFrzV, Name
      USE RIMP2_Module, ONLY : NBF_RI
      USE MP2_Parameter_Module, ONLY : MaxCont
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
      CHARACTER(LEN=255) :: FBuf
      INTEGER :: NOccA, NOccB, NVirA, NVirB
      INTEGER :: IErr
!
!     o Read SCF information
!
      IF (MPIIO) THEN
         FBuf = TRIM(Name)//".SCF_Info"
         OPEN(UNIT=IO, FILE=TRIM(FBuf), STATUS='OLD', ACCESS='SEQUENTIAL', FORM='FORMATTED')
         REWIND(IO)
         READ(IO, *) UHF, NBF, NMO, NOccA, NOccB
         CLOSE(IO)
      END IF
!
      CALL MPI_Bcast(UHF, 1, MPI_LOGICAL, IORank, MPI_COMM_IO, IErr)
      CALL MPI_Bcast(NBF, 1, MPI_INTEGER, IORank, MPI_COMM_IO, IErr)
      CALL MPI_Bcast(NMO, 1, MPI_INTEGER, IORank, MPI_COMM_IO, IErr)
      CALL MPI_Bcast(NOccA, 1, MPI_INTEGER, IORank, MPI_COMM_IO, IErr)
      CALL MPI_Bcast(NOccB, 1, MPI_INTEGER, IORank, MPI_COMM_IO, IErr)
!
!      WRITE(*, *) UHF, NBF, NMO, NOccA, NOccB
      IF (NBF /= MaxCont) THEN
         CALL Util_AbortMPI('Error: Check NBF')
      END IF
!
      NAlpBet = 1
      IF (UHF) NAlpBet = 2
!
      NBC = (NBF * (NBF + 1)) / 2
      NVirA  = NMO - NOccA
      NVirB  = NMO - NOccB
      NActOA = NOccA - NFrzOA
      NActOB = NOccB - NFrzOB
      NActVA = NVirA - NFrzVA
      NActVB = NVirB - NFrzVB
!
      NActO(1) = NActOA
      NActO(2) = NActOB
      NActV(1) = NActVA
      NActV(2) = NActVB
      NFrzO(1) = NFrzOA
      NFrzO(2) = NFrzOB
      NFrzV(1) = NFrzVA
      NFrzV(2) = NFrzVB
!
      IF (MyRank == 0) THEN
         WRITE(*, *) 'o NBF      = ', NBF
         WRITE(*, *) 'o NMO      = ', NMO
         WRITE(*, *) 'o NOccA    = ', NOccA
         WRITE(*, *) 'o NOccB    = ', NOccB
         WRITE(*, *) 'o NVirA    = ', NVirA
         WRITE(*, *) 'o NVirB    = ', NVirB
         WRITE(*, *) 'o NActOA   = ', NActOA
         WRITE(*, *) 'o NActOB   = ', NActOB
         WRITE(*, *) 'o NActVA   = ', NActVA
         WRITE(*, *) 'o NActVB   = ', NActVB
         WRITE(*, *) 'o NFrzOA   = ', NFrzOA
         WRITE(*, *) 'o NFrzOB   = ', NFrzOB
         WRITE(*, *) 'o NFrzVA   = ', NFrzVA
         WRITE(*, *) 'o NFrzVB   = ', NFrzVB
         IF (DoRI) THEN
            WRITE(*, *) 'o NBF_RI   = ', NBF_RI
         END IF
      END IF
!
      END SUBROUTINE
