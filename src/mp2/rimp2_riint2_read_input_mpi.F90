      SUBROUTINE RIMP2_RIInt2_Read_Input_MPI
!
      USE MP2_Module, ONLY : ThrPre, PScreen
      USE Int2_Module, ONLY : IntType, ThrInt, ThrPrim, DoMD4
      USE RIInt2_Module, ONLY : ThrRI
      USE MPI_Module, ONLY : MyRank, MPIIO, IORank, MPI_COMM_IO
!
      IMPLICIT NONE
!
      INCLUDE 'mpif.h'
!
      INTEGER, PARAMETER :: IO = 99
      INTEGER :: IErr
!
      NAMELIST /RIInt2/ IntType, ThrInt, ThrPrim, ThrPre, PScreen, ThrRI
!
      IntType = 'MD4'
      ThrInt = 1.0D-13
      ThrPrim = 1.0D-20
      PScreen = .FALSE.
      ThrPre = 1.0D-12
      ThrRI = 1.0D-6
!
      IF (MPIIO) THEN
!
         OPEN(UNIT=IO, FILE='INPUT', STATUS='OLD', ACCESS='SEQUENTIAL')
!
!        o Read NAMELIST RIInt2
!
         REWIND(IO)
         READ(IO, RIInt2, END=9010, ERR=9000)
         GO TO 9010
 9000    CONTINUE
 9010    CONTINUE
!
         CLOSE(IO)
!
      END IF
!
      CALL MPI_Bcast(IntType, 8, MPI_CHARACTER, IORank, MPI_COMM_IO, IErr)
      CALL MPI_Bcast(PScreen, 1, MPI_LOGICAL, IORank, MPI_COMM_IO, IErr)
      CALL MPI_Bcast(ThrInt, 1, MPI_DOUBLE_PRECISION, IORank, MPI_COMM_IO, IErr)
      CALL MPI_Bcast(ThrPre, 1, MPI_DOUBLE_PRECISION, IORank, MPI_COMM_IO, IErr)
      CALL MPI_Bcast(ThrRI,  1, MPI_DOUBLE_PRECISION, IORank, MPI_COMM_IO, IErr)
!
!     o Transform lower-case letters to upper-case ones
!
      CALL Util_TransChar(IntType, 8)
!
      IntType = 'MD4'
      IF (MyRank == 0) THEN
         WRITE(*,*) 'o IntType  = ', TRIM(IntType)
         WRITE(*,*) 'o ThrInt   = ', ThrInt
         WRITE(*,*) 'o ThrPrim  = ', ThrPrim
         WRITE(*,*) 'o PScreen  = ', PScreen
         WRITE(*,*) 'o ThrPre   = ', ThrPre
         WRITE(*,*) 'o ThrRI    = ', ThrRI
      END IF
!
!     o Type of integrals
!
      DoMD4 = .FALSE.
      IF ((TRIM(IntType) == 'MD') .OR. (TRIM(IntType) == 'MD4')) DoMD4 = .TRUE.
!
      END SUBROUTINE
