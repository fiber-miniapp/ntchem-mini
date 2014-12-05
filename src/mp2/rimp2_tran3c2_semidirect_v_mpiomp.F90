      SUBROUTINE RIMP2_Tran3c2_SemiDirect_V_MPIOMP
!
!     o three-center RI integral transformation from original aux. basis to transformed aux. basis
!       (ia|c) -> (ia|d)
!
      USE MP2_Module, ONLY : IOccBat, NOccBat, LenOccBat, NOccBat_per_Proc, NAlpBet, NActO, NActV, Name
      USE RIMP2_Module, ONLY : NBF_RI, NBC_RI, RIInt2c, RI2cInv, NBF_RI_MyRank, IdxBF_RI_MyRank
      USE MP2_Constant_Module, ONLY : Zero, One
      USE MPI_Module, ONLY : NProcs, MyRank
!
      IMPLICIT NONE
!
      INCLUDE 'mpif.h'
!
#ifdef MPIINT8
#define MPI_INTEGER MPI_INTEGER8
#endif
!
      INTEGER, PARAMETER :: IO1 = 99, IO2 = 98
      CHARACTER(LEN=1), PARAMETER :: AlpBet(1:2) = (/ "A", "B"/) 
      CHARACTER(LEN=255) :: FBuf
      CHARACTER(LEN=10)  :: RankNo
      INTEGER :: MXNActO
      INTEGER :: IAlpBet, IaBat, Ia, Ii, I, K
      INTEGER :: IaBg, IaEd
      INTEGER :: IaBat_Proc, Irank_diff, Iranksend, Irankrecv
      INTEGER :: MXNBF_RI_MyRank
      INTEGER :: NT2BufSize
      INTEGER :: IErr
      INTEGER :: ireq(4)
      REAL(8) :: TimeBgn, TimeEnd, WTimeBgn, WTimeEnd
      INTEGER, ALLOCATABLE :: IdxBF_RI_Irank(:)
      INTEGER, ALLOCATABLE :: istat1(:), istat2(:), istat3(:), istat4(:)
      REAL(8), ALLOCATABLE :: T2Int(:,:,:), T3Int(:,:), T2BufSend(:,:,:), T2BufRecv(:,:,:)
!
      REAL(8) :: Time_T2R, Time_T2C, Time_T3, Time_T3W, WTime_T2R, WTime_T2C, WTime_T3, WTime_T3W
!
      Time_T2R = Zero
      Time_T2C = Zero
      Time_T3 = Zero
      Time_T3W = Zero
      WTime_T2R = Zero
      WTime_T2C = Zero
      WTime_T3 = Zero
      WTime_T3W = Zero
!
      ALLOCATE(RIInt2c(NBC_RI))
      CALL DCOPY(NBC_RI, Zero, 0, RIInt2c, 1)
!
!     o Evaluation of two-center RI integrals (P|Q)
!
      IF (MyRank == 0) THEN
         PRINT '(" ..... Enter    (RIInt2_MDInt2c  )")'
      END IF
      WTimeBgn = MPI_WTIME()
      CALL CPU_TIME(TimeBgn)
      CALL RIMP2_RIInt2_MDInt2_Int2c_MPIOMP
      CALL CPU_TIME(TimeEnd)
      WTimeEnd = MPI_WTIME()
      IF (MyRank == 0) THEN
         PRINT '(" ..... CPU time (RIInt2_MDInt2c  ) :", F12.2)', TimeEnd - TimeBgn
         PRINT '(" ..... WALL time (RIInt2_MDInt2c  ) :", F12.2)', WTimeEnd - WTimeBgn
      END IF
!
      ALLOCATE(RI2cInv(NBF_RI, NBF_RI))
!
!     o Calculation of the inverse matrix of two-center RI integrals (P|Q)^-1
!
      IF (MyRank == 0) THEN
         PRINT '(" ..... Enter    (RIInt2_Inv2c    )")'
      END IF
      WTimeBgn = MPI_WTIME()
      CALL CPU_TIME(TimeBgn)
      CALL RIMP2_Inv2c_MPI
      CALL CPU_TIME(TimeEnd)
      WTimeEnd = MPI_WTIME()
      IF (MyRank == 0) THEN
         PRINT '(" ..... CPU time (RIInt2_Inv2c    ) :", F12.2)', TimeEnd - TimeBgn
         PRINT '(" ..... WALL time (RIInt2_Inv2c    ) :", F12.2)', WTimeEnd - WTimeBgn
      END IF
!
      DEALLOCATE(RIInt2c)
!
!     o Allocation of memory
!
      MXNActO = MAX(NActO(1), NActO(2))
      MXNBF_RI_MyRank = MAXVAL(NBF_RI_MyRank)
      NT2BufSize = MXNActO * LenOccBat * MXNBF_RI_MyRank
      ALLOCATE(T2Int(NBF_RI,MXNActO,LenOccBat))
      ALLOCATE(IdxBF_RI_Irank(MXNBF_RI_MyRank))
      ALLOCATE(istat1(MPI_STATUS_SIZE))
      ALLOCATE(istat2(MPI_STATUS_SIZE))
      ALLOCATE(istat3(MPI_STATUS_SIZE))
      ALLOCATE(istat4(MPI_STATUS_SIZE))
!
      WRITE(RankNo, '(I0)') MyRank
!
!     o Open three-center MO integral file
!
      FBuf = TRIM(Name)//".MOInt3ci"//AlpBet(1)//"."//TRIM(RankNo)
      OPEN(UNIT=IO1, FILE=TRIM(FBuf), STATUS='UNKNOWN', ACCESS='DIRECT', FORM='UNFORMATTED', RECL=(8*NActO(1)))
      FBuf = TRIM(Name)//".MOInt3c"//AlpBet(1)//"."//TRIM(RankNo)
      OPEN(UNIT=IO2, FILE=TRIM(FBuf), STATUS='UNKNOWN', ACCESS='DIRECT', FORM='UNFORMATTED', RECL=(8*NBF_RI*NActO(1)))
!
!MPI Parallel
      DO IaBat_Proc = 1, NOccBat_per_Proc
!MPI Parallel
!
!        o Send and recieve the 2/3 transformed 3c RI integrals (ia|c)
!
         ALLOCATE(T2BufSend(MXNActO,MXNBF_RI_MyRank,LenOccBat))
         ALLOCATE(T2BufRecv(MXNActO,MXNBF_RI_MyRank,LenOccBat))
         DO Irank_diff = 0, (NProcs - 1)
            Irankrecv = MyRank + Irank_diff
            Iranksend = MyRank - Irank_diff
            IF (Irankrecv >= Nprocs) Irankrecv = Irankrecv - NProcs
            IF (Iranksend < 0)       Iranksend = Iranksend + NProcs
            IaBat = Iranksend * NOccBat_per_Proc + IaBat_Proc
            IaBg = IOccBat(1,IaBat,1) + 1
            IaEd = IOccBat(1,IaBat,1) + IOccBat(2,IaBat,1)
!
            WTimeBgn = MPI_WTIME()
            CALL CPU_TIME(TimeBgn)
            DO K = 1, NBF_RI_MyRank(MyRank)
               I = 0
               DO Ia = IaBg, IaEd
                  I = I + 1
                  READ(IO1, REC=(Ia+(K-1)*NActV(1))) T2BufSend(1:NActO(1),K,I)
               END DO
            END DO
            CALL CPU_TIME(TimeEnd)
            WTimeEnd = MPI_WTIME()
            Time_T2R = Time_T2R + TimeEnd - TimeBgn
            WTime_T2R = WTime_T2R + WTimeEnd - WTimeBgn
!
            WTimeBgn = MPI_WTIME()
            CALL CPU_TIME(TimeBgn)
            CALL MPI_ISend(IdxBF_RI_MyRank, NBF_RI_MyRank(MyRank), MPI_INTEGER, Iranksend, 0, &
     &         MPI_COMM_WORLD, ireq(1), IErr)
            CALL MPI_IRecv(IdxBF_RI_Irank,  NBF_RI_MyRank(Irankrecv), MPI_INTEGER, Irankrecv, 0, &
     &         MPI_COMM_WORLD, ireq(2), IErr)
            CALL MPI_ISend(T2BufSend, NT2BufSize, MPI_DOUBLE_PRECISION, Iranksend, 0, MPI_COMM_WORLD, ireq(3), IErr)
            CALL MPI_IRecv(T2BufRecv, NT2BufSize, MPI_DOUBLE_PRECISION, Irankrecv, 0, MPI_COMM_WORLD, ireq(4), IErr)
            CALL MPI_Wait(ireq(1), istat1, IErr)
            CALL MPI_Wait(ireq(2), istat2, IErr)
            CALL MPI_Wait(ireq(3), istat3, IErr)
            CALL MPI_Wait(ireq(4), istat4, IErr)
            CALL CPU_TIME(TimeEnd)
            WTimeEnd = MPI_WTIME()
            Time_T2C = Time_T2C + TimeEnd - TimeBgn
            WTime_T2C = WTime_T2C + WTimeEnd - WTimeBgn
!
            IaBat = MyRank * NOccBat_per_Proc + IaBat_Proc
            DO I = 1, IOccBat(2,IaBat,1)
               DO K = 1, NBF_RI_MyRank(Irankrecv)
                  DO Ii = 1, NActO(1)
                     T2Int(IdxBF_RI_Irank(K),Ii,I) = T2BufRecv(Ii,K,I)
                  END DO
               END DO
            END DO
         END DO
         DEALLOCATE(T2BufSend)
         DEALLOCATE(T2BufRecv)
!
         ALLOCATE(T3Int(NBF_RI,MXNActO))
         IaBat = MyRank * NOccBat_per_Proc + IaBat_Proc
         IaBg = IOccBat(1,IaBat,1) + 1
         IaEd = IOccBat(1,IaBat,1) + IOccBat(2,IaBat,1)
!
         I = 0
         DO Ia = IaBg, IaEd
            I = I + 1
!
!           o 3/3 integral transformation
!             (ia|P) -> (iq|Q)
!
            WTimeBgn = MPI_WTIME()
            CALL CPU_TIME(TimeBgn)
            CALL DGEMM('T', 'N', NBF_RI, NActO(1), NBF_RI, One, RI2cInv, NBF_RI, T2Int(1,1,I), NBF_RI, Zero, &
     &         T3Int, NBF_RI)
            CALL CPU_TIME(TimeEnd)
            WTimeEnd = MPI_WTIME()
            Time_T3 = Time_T3 + TimeEnd - TimeBgn
            WTime_T3 = WTime_T3 + WTimeEnd - WTimeBgn
!
!           o Storage of three-center integrals (ia|Q) to disk
!
            WTimeBgn = MPI_WTIME()
            CALL CPU_TIME(TimeBgn)
            WRITE(IO2, REC=Ia) T3Int(1:NBF_RI,1:NActO(1))
            CALL CPU_TIME(TimeEnd)
            WTimeEnd = MPI_WTIME()
            Time_T3W = Time_T3W + TimeEnd - TimeBgn
            WTime_T3W = WTime_T3W + WTimeEnd - WTimeBgn
!
         END DO
         DEALLOCATE(T3Int)
!
      END DO
!
!     o Close files
!
      CLOSE(IO1, STATUS='DELETE')
      CLOSE(IO2)

      IF (MyRank == 0) THEN
         PRINT '(" ..... CPU time (2/3 tran3c2 read) :", F12.2)', Time_T2R
         PRINT '(" ..... CPU time (2/3 tran3c2 comm) :", F12.2)', Time_T2C
         PRINT '(" ..... CPU time (3/3 tran3c2 tran) :", F12.2)', Time_T3
         PRINT '(" ..... CPU time (3/3 tran3c2 writ) :", F12.2)', Time_T3W
         PRINT '(" ..... WALL time (2/3 tran3c2 read) :", F12.2)', WTime_T2R
         PRINT '(" ..... WALL time (2/3 tran3c2 comm) :", F12.2)', WTime_T2C
         PRINT '(" ..... WALL time (3/3 tran3c2 tran) :", F12.2)', WTime_T3
         PRINT '(" ..... WALL time (3/3 tran3c2 writ) :", F12.2)', WTime_T3W
      END IF
!
!     o deallocate memory
!
      DEALLOCATE(RI2cInv)
      DEALLOCATE(T2Int)
      DEALLOCATE(IdxBF_RI_MyRank)
      DEALLOCATE(NBF_RI_MyRank)
      DEALLOCATE(IdxBF_RI_Irank)
      DEALLOCATE(istat1)
      DEALLOCATE(istat2)
      DEALLOCATE(istat3)
      DEALLOCATE(istat4)
!
      END SUBROUTINE
