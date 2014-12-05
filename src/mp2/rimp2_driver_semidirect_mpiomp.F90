      SUBROUTINE RIMP2_Driver_SemiDirect_MPIOMP
!
      USE MPI_Module, ONLY : NProcs, MyRank
!
!     o Driver subroutine for RI-MP2 energy evaluation
!
      IMPLICIT NONE
!
      INCLUDE 'mpif.h'
!
      INTEGER :: IErr
      REAL(8) :: TimeBgn, TimeEnd, WTimeBgn, WTimeEnd
!
!     o Obtaining RI-MP2 batch infomation and memory allocation
!
      CALL RIMP2_Get_BatchInfo_MPI
!
!     o Initialization of the RI integral evaluation
!
      IF (MyRank == 0) THEN
         PRINT '(" ...   Enter    (RIInt2_Initial  )")'
      END IF
      WTimeBgn = MPI_WTIME()
      CALL CPU_TIME(TimeBgn)
      CALL RIMP2_Int2_Initial_MPIOMP
      CALL CPU_TIME(TimeEnd)
      WTimeEnd = MPI_WTIME()
      IF (MyRank == 0) THEN
         PRINT '(" ...   CPU time (RIInt2_Initial  ) :", F12.2)', TimeEnd - TimeBgn
         PRINT '(" ...   WALL time (RIMP2_Initial   ) :", F12.2)', WTimeEnd - WTimeBgn
      END IF
!
!     o Construction of (ia|P) integrals
!
      IF (MyRank == 0) THEN
         PRINT '(" ...   Enter    (RIMP2_Tran3c1   )")'
      END IF
      WTimeBgn = MPI_WTIME()
      CALL CPU_TIME(TimeBgn)
      CALL RIMP2_Tran3c1_SemiDirect_V_MPIOMP
      CALL MPI_Barrier(MPI_COMM_WORLD, IErr)
      CALL CPU_TIME(TimeEnd)
      WTimeEnd = MPI_WTIME()
      IF (MyRank == 0) THEN
         PRINT '(" ...   CPU time (RIMP2_Tran3c1   ) :", F12.2)', TimeEnd - TimeBgn
         PRINT '(" ...   WALL time (RIMP2_Tran3c1   ) :", F12.2)', WTimeEnd - WTimeBgn
      END IF
!
!     o Construction of (ia|Q) integrals
!
      IF (MyRank == 0) THEN
         PRINT '(" ...   Enter    (RIMP2_Tran3c2   )")'
      END IF
      WTimeBgn = MPI_WTIME()
      CALL CPU_TIME(TimeBgn)
      CALL RIMP2_Tran3c2_SemiDirect_V_MPIOMP
      CALL MPI_Barrier(MPI_COMM_WORLD, IErr)
      CALL CPU_TIME(TimeEnd)
      WTimeEnd = MPI_WTIME()
      IF (MyRank == 0) THEN
         PRINT '(" ...   CPU time (RIMP2_Tran3c2   ) :", F12.2)', TimeEnd - TimeBgn
         PRINT '(" ...   WALL time (RIMP2_Tran3c2   ) :", F12.2)', WTimeEnd - WTimeBgn
      END IF
!
!     o Finalization of the RI integral evaluation
!
      CALL RIMP2_RIInt2_Final_MPIOMP
!
!     o Construction of (ia|jb) integrals and evaluation of MP2 correlation energy!
!     o RMP2 case
!
      IF (MyRank == 0) THEN
         PRINT '(" ...   Enter    (RIMP2_RMP2Energy)")'
      END IF
      WTimeBgn = MPI_WTIME()
      CALL CPU_TIME(TimeBgn)
      CALL RIMP2_RMP2Energy_SemiDirect_V_MPIOMP
      CALL CPU_TIME(TimeEnd)
      WTimeEnd = MPI_WTIME()
      IF (MyRank == 0) THEN
         PRINT '(" ...   CPU time (RIMP2_RMP2Energy) :", F12.2)', TimeEnd - TimeBgn
         PRINT '(" ...   WALL time (RIMP2_RMP2Energy) :", F12.2)', WTimeEnd - WTimeBgn
      END IF
!
!     o Memory deallocation
!
      CALL RIMP2_Deallocate
!
      END SUBROUTINE
