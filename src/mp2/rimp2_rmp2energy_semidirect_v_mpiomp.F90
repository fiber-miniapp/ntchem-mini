      SUBROUTINE RIMP2_RMP2Energy_SemiDirect_V_MPIOMP
!
!     o 4c integral generation and RMP2 energy accumulation
!
      USE MP2_Module, ONLY : IOccBat, NOccBat, LenOccBat, NOccBat_per_Proc, NMO, &
     &   NActO, NActV, NFrzO, EMP2, ESCSMP2, E1, E2T, E2S, E2, E2SCS, Name, IPrint
      USE RIMP2_Module, ONLY : NBF_RI
      Use MP2_Constant_Module, ONLY : Zero, One, Two, Three, P12
      USE MPI_Module, ONLY : NProcs, MyRank, MPIIO, IORank, MPI_COMM_IO
!
      IMPLICIT NONE
!
      INCLUDE 'mpif.h'
!
      INTEGER, PARAMETER :: IO = 99
      CHARACTER(LEN=255) :: FBuf
      CHARACTER(LEN=10)  :: RankNo
      REAL(8) :: E2TP, E2SP, E2Tab, E2Sab
      REAL(8) :: T2, Fac
      REAL(8) :: EigIb, EigIab, EigIjab, EigIi
      INTEGER :: IaBat, IbBat, Ii, Ij, Ia, Ib
      INTEGER :: IaBg, IaEd, IbBg, IbEd
      INTEGER :: I, J
      INTEGER :: NMOInt3BufSize
      INTEGER :: NHOMO
      INTEGER :: IaBat_Proc, IbBat_Proc, Jranksend, Jrankrecv, Jrank_diff, IaBat_Proc_End
      INTEGER :: NOccBat_per_Proc_half, NProcs_half, IaBatBg
      INTEGER :: IErr
      INTEGER :: ireq(2)
      LOGICAL :: EvenProcs, ExchIBat
      INTEGER, ALLOCATABLE :: istat1(:), istat2(:)
      REAL(8), ALLOCATABLE :: MOInt3cb(:,:), MOInt3ck(:,:), MOInt(:,:), Eig(:)
      REAL(8), ALLOCATABLE :: MOInt3ck0(:,:)
!
      REAL(8) :: TimeBgn, TimeEnd, Time_MOI, Time_EMP2, WTimeBgn, WTimeEnd, WTime_MOI, WTime_EMP2
      REAL(8) :: Time_T3RK, Time_T3C, Time_T3RB, WTime_T3RK, WTime_T3C, WTime_T3RB
!
      Time_T3RK = Zero
      Time_T3C = Zero
      Time_T3RB = Zero
      Time_MOI = Zero
      Time_EMP2 = Zero
      WTime_T3RK = Zero
      WTime_T3C = Zero
      WTime_T3RB = Zero
      WTime_MOI = Zero
      WTime_EMP2 = Zero
!
!     o Memory allocation
!
      ALLOCATE(MOInt3cb(NBF_RI*NActO(1),LenOccBat))
      ALLOCATE(MOInt3ck(NBF_RI*NActO(1),LenOccBat))
      ALLOCATE(MOInt(NActO(1),NActO(1)))
      ALLOCATE(Eig(NMO))
      ALLOCATE(MOInt3ck0(NBF_RI*NActO(1),LenOccBat))
      ALLOCATE(istat1(MPI_STATUS_SIZE))
      ALLOCATE(istat2(MPI_STATUS_SIZE))
!
      NMOInt3BufSize = NBF_RI * NActO(1) * LenOccBat
      NHOMO = NFrzO(1) + NActO(1)
!
      EvenProcs = .FALSE.
      IF (MOD(NProcs, 2) == 0) THEN
         EvenProcs = .TRUE.
      END IF
!
!     o Read orbital energies
!
      IF (MPIIO) THEN
         FBuf = TRIM(Name)//".OrbEne"
         OPEN(UNIT=IO, FILE=TRIM(FBuf), STATUS='OLD', ACCESS='SEQUENTIAL', FORM='FORMATTED')
         READ(IO, *) Eig(1:NMO)
         CLOSE(IO)
      END IF
      CALL MPI_Bcast(Eig, NMO, MPI_DOUBLE_PRECISION, IORank, MPI_COMM_IO, IErr)
!
      IF ((MyRank == 0) .AND. (IPrint >= 1)) THEN
         WRITE(*, *) '+++++ Orbital energy (Alpha) +++++'
         WRITE(*, '(10F12.6)') Eig(1:NMO)
      END IF
!
!     o Open 3c MO integral file (ia|q)
!
      WRITE(RankNo,'(I0)') MyRank
      FBuf = TRIM(Name)//".MOInt3cA."//Trim(RankNo)
      OPEN(UNIT=IO, FILE=TRIM(FBuf), STATUS='UNKNOWN', ACCESS='DIRECT', FORM='UNFORMATTED', RECL=(8*NBF_RI*NActO(1)))
!
!     o Calculation of RMP2 correlation energy
!
      NProcs_half = NProcs / 2
      NOccBat_per_Proc_half = NOccBat_per_Proc / 2
!
      E2TP = Zero
      E2SP = Zero
!
!MPI Parallel
      DO IbBat_Proc = 1, NOccBat_per_Proc
!MPI Parallel
         IbBat = MyRank * NOccBat_per_Proc + IbBat_Proc
         IbBg = IOccBat(1,IbBat,1) + 1
         IbEd = IOccBat(1,IbBat,1) + IOccBat(2,IbBat,1)
!
!        o Reading of three-center MO integral (jb|Q) from file
!
         WTimeBgn = MPI_WTIME()
         CALL CPU_TIME(TimeBgn)
         J = 0
         DO Ib = IbBg, IbEd
            J = J + 1
            READ(UNIT=IO, REC=Ib) MOInt3ck0(:,J)
         END DO
         CALL CPU_TIME(TimeEnd)
         WTimeEnd = MPI_WTIME()
         Time_T3RK = Time_T3RK + TimeEnd - TimeBgn
         WTime_T3RK = WTime_T3RK + WTimeEnd - WTimeBgn
!
         DO Jrank_diff = 0, NProcs_half
            ExchIBat = EvenProcs .AND. (jrank_diff == NProcs_half)
            Jrankrecv = MyRank + jrank_diff
            Jranksend = MyRank - jrank_diff
            IF (Jrankrecv >= Nprocs) Jrankrecv = Jrankrecv - NProcs
            IF (Jranksend < 0)       Jranksend = Jranksend + NProcs
            IbBat = Jrankrecv * NOccBat_per_Proc + IbBat_Proc
            IbBg = IOccBat(1,IbBat,1) + 1
            IbEd = IOccBat(1,IbBat,1) + IOccBat(2,IbBat,1)
!
!           o Communicate three-center MO integral (jb|Q) with each process
!
            WTimeBgn = MPI_WTIME()
            CALL CPU_TIME(TimeBgn)
            CALL MPI_ISend(MOInt3ck0, NMOInt3BufSize, MPI_DOUBLE_PRECISION, Jranksend, 0, MPI_COMM_WORLD, ireq(1), IErr)
            CALL MPI_IRecv(MOInt3ck,  NMOInt3BufSize, MPI_DOUBLE_PRECISION, Jrankrecv, 0, MPI_COMM_WORLD, ireq(2), IErr)
            CALL MPI_Wait(ireq(1), istat1, IErr)
            CALL MPI_Wait(ireq(2), istat2, IErr)
            CALL CPU_TIME(TimeEnd)
            WTimeEnd = MPI_WTIME()
            Time_T3C = Time_T3C + TimeEnd - TimeBgn
            WTime_T3C = WTime_T3C + WTimeEnd - WTimeBgn
!
            IaBatBg = MyRank * NOccBat_per_Proc
            IF (ExchIBat .AND. (MyRank <= Jrankrecv)) THEN
               IF (IbBat_Proc > NOccBat_per_Proc_half) THEN
                  IaBatBg = MyRank * NOccBat_per_Proc + NOccBat_per_Proc / 2
               END IF
            END IF
            IF (ExchIBat .AND. (MyRank > Jrankrecv)) THEN
               IF (IbBat_Proc <= NOccBat_per_Proc_half) THEN
                  IaBatBg = MyRank * NOccBat_per_Proc + NOccBat_per_Proc / 2
               END IF
            END IF
!
            IaBat_Proc_End = NOccBat_per_Proc
            IF (Jrank_diff == 0) THEN
               IaBat_Proc_End = IbBat_Proc
            ELSE IF (ExchIBat) THEN
               IaBat_Proc_End = NOccBat_per_Proc_half
            END IF
!MPI Parallel
            DO IaBat_Proc = 1, IaBat_Proc_End
!MPI Parallel
               IaBat = IaBatBg + IaBat_Proc
               IaBg = IOccBat(1,IaBat,1) + 1
               IaEd = IOccBat(1,IaBat,1) + IOccBat(2,IaBat,1)
!
!              o Reading of three-center MO integral (ia|Q) from file
!
               WTimeBgn = MPI_WTIME()
               CALL CPU_TIME(TimeBgn)
               I = 0
               DO Ia = IaBg, IaEd
                  I = I + 1
                  READ(UNIT=IO, REC=Ia) MOInt3cb(:,I)
               END DO
               CALL CPU_TIME(TimeEnd)
               WTimeEnd = MPI_WTIME()
               Time_T3RB = Time_T3RB + TimeEnd - TimeBgn
               WTime_T3RB = WTime_T3RB + WTimeEnd - WTimeBgn
!
               J = 0
               DO Ib = IbBg, IbEd
                  EigIb = - Eig(Ib+NHOMO)
                  J = J + 1
                  IF (IaBat == IbBat) IaEd = Ib
                  I = 0
                  DO Ia = IaBg, IaEd
                     EigIab = EigIb - Eig(Ia+NHOMO)
                     I = I + 1
!
!                    o Evaluation of four-center MO integrals (ia|jb) from three-center integrals
!
                     WTimeBgn = MPI_WTIME()
                     CALL CPU_TIME(TimeBgn)
                     CALL DGEMM('T', 'N', NActO(1), NActO(1), NBF_RI, One, MOInt3cb(1,I), NBF_RI, &
     &                  MOInt3ck(1,J), NBF_RI, Zero, MOInt, NActO(1))
                     CALL CPU_TIME(TimeEnd)
                     WTimeEnd = MPI_WTIME()
                     Time_MOI = Time_MOI + TimeEnd - TimeBgn
                     WTime_MOI = WTime_MOI + WTimeEnd - WTimeBgn
!
                     E2Tab = Zero
                     E2Sab = Zero
!
!                    o Evaluation of MP2 correlation energy for ij orbital pair
!
                     WTimeBgn = MPI_WTIME()
                     CALL CPU_TIME(TimeBgn)
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(Ii, Ij, EigIjab, EigIi, T2)
!$OMP DO REDUCTION(+: E2Tab, E2Sab)
                     DO Ij = 1, NActO(1)
                        EigIjab = EigIab + Eig(Ij+NFrzO(1))
                        DO Ii = 1, NActO(1)
                           EigIi = Eig(Ii+NFrzO(1))
                           T2 = MOInt(Ii,Ij) / (EigIjab + EigIi)
                           E2Tab = E2Tab + T2 * (MOInt(Ii,Ij) - MOInt(Ij,Ii))
                           E2Sab = E2Sab + T2 * MOInt(Ii,Ij)
                        END DO
                     END DO
!$OMP END DO
!$OMP END PARALLEL
!
                     IF (Ia /= Ib) THEN
                        Fac = Two
                     ELSE
                        Fac = One
                     END IF
                     E2TP = E2TP + Fac * E2Tab
                     E2SP = E2SP + Fac * E2Sab
                     CALL CPU_TIME(TimeEnd)
                     WTimeEnd = MPI_WTIME()
                     Time_EMP2 = Time_EMP2 + TimeEnd - TimeBgn
                     WTime_EMP2 = WTime_EMP2 + WTimeEnd - WTimeBgn
!
                  END DO
               END DO
!
            END DO
         END DO
      END DO
!
!     o Close 3c MO integral file
!
      CLOSE(IO, STATUS='DELETE')
!
!     o Memory deallocation
!
      DEALLOCATE(MOInt3cb)
      DEALLOCATE(MOInt3ck)
      DEALLOCATE(MOInt)
      DEALLOCATE(Eig)
      DEALLOCATE(MOInt3ck0)
      DEALLOCATE(istat1)
      DEALLOCATE(istat2)
!
!     o Correct MP2 correlation energy to master process (rank=0)
!
      CALL MPI_Reduce(E2SP, E2S, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, IErr)
      CALL MPI_Reduce(E2TP, E2T, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, IErr)
!
      IF (MyRank == 0) THEN
!
!        o Read the SCF energy
!
         FBuf = TRIM(Name)//".TotEne"
         OPEN(UNIT=IO, FILE=TRIM(FBuf), STATUS='OLD', ACCESS='SEQUENTIAL', FORM='FORMATTED')
         READ(IO, *) EMP2
         CLOSE(IO)
!
!        o Write the MP2 correlation energy
!
         E2 = E2S + E2T
         E2SCS = E2S * P12 + E2T / Three
!
         WRITE(*, *) 'SCF energy                   =', EMP2
         WRITE(*, *) 'MP1 energy                   =', E1
         WRITE(*, *) 'MP2 energy (Singlet corr )   =', E2S
         WRITE(*, *) 'MP2 energy (Triplet corr )   =', E2T
         WRITE(*, *) 'MP2 energy (Total corr )     =', E2
         WRITE(*, *) 'SCS-MP2 energy (Total corr ) =', E2SCS
!
!        o Write the total MP2 energy
!
         ESCSMP2 = EMP2
         EMP2 = EMP2 + E1 + E2
         ESCSMP2 = ESCSMP2 + E1 + E2SCS
!
         WRITE(*, *) 'Total MP2 energy         =', EMP2
         WRITE(*, *) 'Total SCS-MP2 energy     =', ESCSMP2
!
         PRINT '(" ..... CPU time (3/3k 2cints read) :", F12.2)', Time_T3RK
         PRINT '(" ..... CPU time (3/3k 2cints comm) :", F12.2)', Time_T3C
         PRINT '(" ..... CPU time (3/3b 2cints read) :", F12.2)', Time_T3RB
         PRINT '(" ..... CPU time (4c Ints         ) :", F12.2)', Time_MOI
         PRINT '(" ..... CPU time (EMP2 corr.      ) :", F12.2)', Time_EMP2
         PRINT '(" ..... WALL time (3/3k 2cints read) :", F12.2)', WTime_T3RK
         PRINT '(" ..... WALL time (3/3k 2cints comm) :", F12.2)', WTime_T3C
         PRINT '(" ..... WALL time (3/3b 2cints read) :", F12.2)', WTime_T3RB
         PRINT '(" ..... WALL time (4c Ints         ) :", F12.2)', WTime_MOI
         PRINT '(" ..... WALL time (EMP2 corr.      ) :", F12.2)', WTime_EMP2
!
      END IF
!
      END SUBROUTINE
