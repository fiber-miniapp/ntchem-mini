      SUBROUTINE RIMP2_Driver_InCore_MPIOMP
!
      USE MP2_Module, ONLY : IOccBat, LenOccBat, NOccBat_per_Proc, NActO, NActV
      USE RIMP2_Module, ONLY : NBF_RI, RIInt3c2a, RIInt3c2b, RIInt3c3a, RIInt3c3b
      USE MP2_Basis_Module, ONLY : Spherical, LtuvMin_Car, LtuvMin_Sph, LtuvMax_Car, LtuvMax_Sph
      USE RIMP2_Basis_Module, ONLY : NShel_RI, KType_RI
      USE MPI_Module, ONLY : NProcs, MyRank, NProcsMat, MyRankMat, NProcsMO, MyRankMO
!
!     o Driver subroutine for RI-MP2 energy evaluation
!
      IMPLICIT NONE
!
      INCLUDE 'mpif.h'
!
      INTEGER :: LenOccBat_per_Proc, NBF_RI_per_ProcMat
      INTEGER :: KK, KBF_RI, NK, IAnglC
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
!
      KBF_RI = 0
!MPI parallel
      DO KK = 1, NShel_RI
         IF (MOD(KK, NProcsMO) /= MyRankMO) CYCLE
!MPI parallel
!
         IAnglC = KType_RI(KK)
         IF (Spherical) THEN
            NK = LtuvMax_Sph(IAnglC) - LtuvMin_Sph(IAnglC) + 1
         ELSE
            NK = LtuvMax_Car(IAnglC) - LtuvMin_Car(IAnglC) + 1
         END IF
         KBF_RI = KBF_RI + NK
!
      END DO
      ALLOCATE(RIInt3c2a(NActO(1)*NActV(1),KBF_RI))

      ! test
      RIInt3c2a(:,:) = 0.0
!
      CALL RIMP2_Tran3c1_InCore_V_MPIOMP
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
!
      LenOccBat_per_Proc = LenOccBat * NOccBat_per_Proc
      NBF_RI_per_ProcMat = NBF_RI / NProcsMat
      if (MOD(NBF_RI, NProcsMat) > MyRankMat) then
         NBF_RI_per_ProcMat = NBF_RI_per_ProcMat + 1
      end if
      ALLOCATE(RIInt3c3a(NBF_RI_per_ProcMat*NActO(1),LenOccBat_per_Proc))
      CALL RIMP2_Tran3c2_InCore_V_MPIOMP
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
      DEALLOCATE(RIInt3c2a)
!
!     o Construction of (ia|jb) integrals and evaluation of MP2 correlation energy
!     o RMP2 case
!
      IF (MyRank == 0) THEN
         PRINT '(" ...   Enter    (RIMP2_RMP2Energy)")'
      END IF
      WTimeBgn = MPI_WTIME()
      CALL CPU_TIME(TimeBgn)
      CALL RIMP2_RMP2Energy_InCore_V_MPIOMP
      CALL CPU_TIME(TimeEnd)
      WTimeEnd = MPI_WTIME()
      IF (MyRank == 0) THEN
         PRINT '(" ...   CPU time (RIMP2_RMP2Energy) :", F12.2)', TimeEnd - TimeBgn
         PRINT '(" ...   WALL time (RIMP2_RMP2Energy) :", F12.2)', WTimeEnd - WTimeBgn
      END IF
!
!     o Memory deallocation
!
      DEALLOCATE(RIInt3c3a)
      DEALLOCATE(IOccBat)
!
      END SUBROUTINE
