      SUBROUTINE RIMP2_RMP2Energy_InCore_V_MPIOMP
!
!     o 4c integral generation and RMP2 energy accumulation
!
      USE MP2_Module, ONLY : IOccBat, NOccBat, LenOccBat, NOccBat_per_Proc, NMO, &
     &   NActO, NActV, NFrzO, EMP2, ESCSMP2, E1, E2T, E2S, E2, E2SCS, Name, IPrint
      USE RIMP2_Module, ONLY : NBF_RI, RIInt3c3a
      Use MP2_Constant_Module, ONLY : Zero, One, Two, Three, P12
      USE MPI_Module, ONLY : NProcs, MyRank, MPIIO, IORank, MPI_COMM_IO, MPI_COMM_MO, NProcsMO, MyRankMO, &
     &   MPI_COMM_MAT, NProcsMat, MyRankMat
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
      REAL(8) :: E2TP, E2SP, E2Tab, E2Sab
      REAL(8) :: T2, Fac
      REAL(8) :: EigIb, EigIab, EigIjab, EigIi
      INTEGER :: IaBat, IbBat, Ii, Ij, Ia, Ib
      INTEGER :: IaBg, IaEd, IbBg, IbEd, IaBg_Proc, IbBg_Proc, Ib_Send
      INTEGER :: I, J
      INTEGER :: NMOInt3BufSize
      INTEGER :: NHOMO
      INTEGER :: IaBat_Proc, IbBat_Proc, Jranksend, Jrankrecv, Jrank_diff, IaBat_Proc_End
      integer :: Jranksend_1, Jrankrecv_1
      INTEGER :: NOccBat_per_Proc_half, NProcs_half, IaBatBg
      INTEGER :: IErr
      INTEGER :: ireq(2)
      LOGICAL :: EvenProcs, ExchIBat
      INTEGER :: NBF_RI_per_ProcMat
      INTEGER, ALLOCATABLE :: istat1(:), istat2(:), istat3(:)

      ! REAL(8), ALLOCATABLE, target :: MOInt3ck(:,:)
      REAL(8), pointer :: MOInt3ck(:,:)
#ifndef USE_GPU
      integer, parameter :: NUM_STREAM = 1
      REAL(8), allocatable, target :: RWork2_Pool(:,:,:)
      REAL(8), allocatable, target :: CommBuf(:,:,:)
#else
      ! *** for GPU ***
      integer, parameter :: NUM_STREAM = 3
      REAL(8), allocatable, target, pinned :: RWork2_Pool(:,:,:)
      REAL(8), allocatable, target, pinned :: CommBuf(:,:,:)
#endif
      REAL(8), pointer :: RWork2(:,:)
      REAL(8), pointer :: SendBuf(:,:)
      REAL(8), pointer :: RecvBuf(:,:)
      integer :: SendBufId
      integer :: RecvBufId

      REAL(8), ALLOCATABLE :: Eig(:)
      REAL(8), allocatable :: MOIntSome(:,:)
!
!
      REAL(8) :: TimeBgn, TimeEnd, Time_MOI, Time_EMP2, WTimeBgn, WTimeEnd, WTime_MOI, WTime_EMP2
      REAL(8) :: Time_T3C, Time_MOIC, Time_EMP2C, WTime_T3C, WTime_MOIC, WTime_EMP2C
!
      real(8), parameter :: maxMem = 1e9
      ! real(8), parameter :: maxMem = 5e8
      integer :: id_st
      integer, parameter :: NUM_EVENT = NUM_STREAM*2
      integer :: id_ev, id_ev_next
      integer :: id_A, id_Am
      integer :: id_B, id_Bm
      integer :: maxMN
      integer :: nGrp, nGrpMax, nBlk
      integer :: lenB, lenA
      integer :: OfstI, OfstJ
      integer :: Ia_0, Ia_1    ! Ia = Ia_0 + Ia_1
      integer :: Ib_0, Ib_1    ! Ib = Ib_0 + Ib_1
      integer :: m, n, k
      integer :: lda, ldb, ldc
      integer, allocatable :: devptr_A(:), devptr_B(:), devptr_C(:)
      integer :: count_dgemm   ! debug

      integer :: Mod_NActO_NProcsMat
      integer :: IiBgn, IiEnd, leni
      integer :: Max_LCount, LCount, LNumber, LNumber_Base
      integer, allocatable :: Ia_0_list(:), Ib_0_list(:)
      integer, allocatable :: cflag_set_mat_A(:), cflag_set_mat_B(:)

      integer, allocatable :: commIndexEach(:)
      integer, allocatable :: commSizeEach(:)
      integer :: commPhase, commCount, commSizeTotal
      integer :: chunk, myChunk

      Time_T3C = Zero
      Time_MOI = Zero
      Time_MOIC = Zero
      Time_EMP2 = Zero
      Time_EMP2C = Zero
      WTime_T3C = Zero
      WTime_MOI = Zero
      WTime_MOIC = Zero
      WTime_EMP2 = Zero
      WTime_EMP2C = Zero

      NBF_RI_per_ProcMat = NBF_RI / NProcsMat
      if (MOD(NBF_RI, NProcsMat) > MyRankMat) then
         NBF_RI_per_ProcMat = NBF_RI_per_ProcMat + 1
      end if

      leni = NActO(1) / NProcsMat
      Mod_NActO_NProcsMat = MOD(NActO(1), NProcsMat)
      if (Mod_NActO_NProcsMat > MyRankMat) then
         IiBgn = leni * MyRankMat + MyRankMat + 1
         IiEnd = IiBgn + leni
      else
         IiBgn = leni * MyRankMat + Mod_NActO_NProcsMat + 1
         IiEnd = IiBgn + leni - 1
      end if

!      write(*, *) 'MyRankMat=', MyRankMat, 'NBF_RI_per_ProcMat =', NBF_RI_per_ProcMat

      !
      m = NActO(1)
      n = NActO(1)
      k = NBF_RI_per_ProcMat
      maxMN = sqrt(maxMem/8 + k*k) - k

      nGrp = (NActO(1)*LenOccBat + maxMN-1) / maxMN
      if (nGrp <= 1) nGrp = 2    ! test

      nBlk = (LenOccBat + nGrp-1) / nGrp
#ifndef USE_MERGE
      nBlk = 1   ! When nBlk=1, this code behave as if original one
#endif

      nGrp = (LenOccBat + nBlk-1) / nBlk
      CALL MPI_Allreduce( nGrp, nGrpMax, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr )
      nGrp = nGrpMax

      ALLOCATE(devptr_A(1:NUM_STREAM), devptr_B(1:NUM_STREAM), devptr_C(1:NUM_STREAM))

      Max_LCount = (nGrp+1) * (nGrp+1)
      ALLOCATE(Ia_0_list(Max_LCount), Ib_0_list(Max_LCount))
      ALLOCATE(cflag_set_mat_A(nGrp+1), cflag_set_mat_B(nGrp+1))

      allocate(commSizeEach(Max_LCount), commIndexEach(Max_LCount))

      commCount = (nGrp) * (nGrp+1) / 2

#ifdef DEBUG
      if (MyRank == 0) then
         write(*,*) "# NBF_RI:", NBF_RI
         write(*,*) "# NActO(1):", NActO(1)
         write(*,*) "# LenOccBat:", LenOccBat
         write(*,*) "# maxMem:", maxMem
         write(*,*) "# maxMN:", maxMN
         write(*,*) "# nBlk:", nBlk
         write(*,*) "# nGrp:", nGrp
         write(*,*) "# Max_LCount:", Max_LCount
         write(*,*) "# NUM_STREAM:", NUM_STREAM
         write(*,'(a,3i8)') "# M,N,K:", m, n, k
      endif
#endif
      count_dgemm = 0  ! debug
!
!     o Memory allocation
!
      ! ALLOCATE(MOInt3ck(NBF_RI*NActO(1),LenOccBat))
      ! ALLOCATE(CommBuf(NBF_RI_per_ProcMat*NActO(1),LenOccBat,2))
      ALLOCATE(CommBuf(NBF_RI_per_ProcMat*NActO(1),LenOccBat,2))

      Jrankrecv_1 = mod( MyRankMO + 1 + NProcsMO, NProcsMO )
      Jranksend_1 = mod( MyRankMO - 1 + NProcsMO, NProcsMO )

      m = NActO(1) * nBlk
      n = NActO(1) * nBlk
      k = NBF_RI_per_ProcMat
      lda = k
      ldb = k
      ldc = m
      ALLOCATE(MOIntSome(m,n))
      ALLOCATE(RWork2_Pool(m,n,NUM_STREAM))

#ifdef USE_GPU
      !
      ! initialize cublaas
      !
      CALL CPU_TIME(TimeBgn)
      call mpi_barrier(MPI_COMM_WORLD,ierr)  ! not essential, just make measured time meaningful
      call cublas_init()
      ! allocate memory space for matrix A,B and C on GPU
      do id_st = 1, NUM_STREAM
         call cublas_alloc( devptr_A(id_st), m, k )
         call cublas_alloc( devptr_B(id_st), k, n )
         call cublas_alloc( devptr_C(id_st), m, n )
      enddo
      call mpi_barrier(MPI_COMM_WORLD,ierr)  ! not essential, just make measured time meaningful
      CALL CPU_TIME(TimeEnd)
#ifdef DEBUG
      if ( MyRank == 0 ) then
         write(*,'("# gpu init/alloc time",F12.2)') TimeEnd - TimeBgn
      end if
#endif
#endif

      ALLOCATE(Eig(NMO))
      ALLOCATE(istat1(MPI_STATUS_SIZE))
      ALLOCATE(istat2(MPI_STATUS_SIZE))
      ALLOCATE(istat3(MPI_STATUS_SIZE))
!
      NMOInt3BufSize = NBF_RI_per_ProcMat * NActO(1) * LenOccBat
      NHOMO = NFrzO(1) + NActO(1)
!
      EvenProcs = .FALSE.
      IF (MOD(NProcsMO, 2) == 0) THEN
         EvenProcs = .TRUE.
      END IF

      !
      commIndexEach(:) = 0
      commSizeEach(:) = 0
      commSizeTotal = 0
      chunk = (LenOccBat + commCount - 1) / commCount
      DO commPhase = 1, commCount
         commIndexEach(commPhase) = chunk * (commPhase-1) + 1
         myChunk = LenOccBat - commIndexEach(commPhase) + 1
         if ( myChunk > chunk ) myChunk = chunk
         if ( myChunk < 0 ) myChunk = 0

         commSizeEach(commPhase) = myChunk * (NBF_RI_per_ProcMat * NActO(1))
         commSizeTotal = commSizeTotal + commSizeEach(commPhase)
#ifdef DEBUG
         if ( myRank == 0 ) then
            write(*,'(a,3i6,i10)') &
                 "# comm phase, start index, chunk, comm size, ", commPhase, &
                 commIndexEach(commPhase), &
                 myChunk, &
                 commSizeEach(commPhase)
         endif
#endif
      End DO
      if ( commSizeTotal /= NMOInt3BufSize ) then
         write(*,'(a,2i10)') "# wrong comm size, ", commSizeTotal, NMOInt3BufSize
         stop
      endif
      !

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
!     o Calculation of RMP2 correlation energy
!
      NProcs_half = NProcsMO / 2
      NOccBat_per_Proc_half = NOccBat_per_Proc / 2
!
      E2TP = Zero
      E2SP = Zero
!
      IbBat = MyRankMO * NOccBat_per_Proc + 1
      IbBg_Proc = IOccBat(1,IbBat,1)
      IaBg_Proc = IbBg_Proc
!MPI Parallel
      DO IbBat_Proc = 1, NOccBat_per_Proc
!MPI Parallel

         IbBat = MyRankMO * NOccBat_per_Proc + IbBat_Proc
         IbBg = IOccBat(1,IbBat,1) + 1
         Ib_Send = IbBg - IbBg_Proc
!
#ifdef DEBUG
         if ( MyRank == 0 ) then
            write(*,'(a,i6,a,i6,a,i6)') '# Ib_Send=', Ib_Send,', IbBg=',IbBg,', IbBg_Proc=',IbBg_Proc
         endif
#endif
!
         RecvBuf => RIInt3c3a(:,Ib_Send:)
!
         DO Jrank_diff = 0, NProcs_half
#ifdef DEBUG
            if ( MyRank == 0 ) then
               write(*,'(a,i6,a,i6)') '# Jrank_diff =', Jrank_diff, ', NProcs_half=', NProcs_half
            endif
#endif
            ExchIBat = EvenProcs .AND. (jrank_diff == NProcs_half)
            Jrankrecv = MyRankMO + jrank_diff
            Jranksend = MyRankMO - jrank_diff
            IF (Jrankrecv >= NProcsMO) Jrankrecv = Jrankrecv - NProcsMO
            IF (Jranksend < 0)         Jranksend = Jranksend + NProcsMO
            IbBat = Jrankrecv * NOccBat_per_Proc + IbBat_Proc
            IbBg = IOccBat(1,IbBat,1) + 1
            IbEd = IOccBat(1,IbBat,1) + IOccBat(2,IbBat,1)
!
!           o Communicate three-center MO integral (jb|Q) with each process
!
!            WTimeBgn = MPI_WTIME()
!            CALL CPU_TIME(TimeBgn)

            MOInt3ck => RecvBuf

            if ( (Jrank_diff /= NProcs_half) .OR. (NProcs_half == 0) ) then
               ! call non-blocking MPI reqs to exchange data that will be used at next iteratoin (not for now)
               RecvBufId = mod(Jrank_diff, 2) + 1
               SendBuf => MOInt3ck
               RecvBuf => CommBuf(:,:,RecvBufId)

               ! at once
               ! CALL MPI_ISend(SendBuf, NMOInt3BufSize, MPI_DOUBLE_PRECISION, Jranksend_1, 0, MPI_COMM_MO, ireq(1), IErr)
               ! CALL MPI_IRecv(RecvBuf, NMOInt3BufSize, MPI_DOUBLE_PRECISION, Jrankrecv_1, 0, MPI_COMM_MO, ireq(2), IErr)
            endif

!            CALL CPU_TIME(TimeEnd)
!            WTimeEnd = MPI_WTIME()
!            Time_T3C = Time_T3C + TimeEnd - TimeBgn
!            WTime_T3C = WTime_T3C + WTimeEnd - WTimeBgn
!
            IaBatBg = MyRankMO * NOccBat_per_Proc
            IF (ExchIBat .AND. (MyRankMO <= Jrankrecv)) THEN
               IF (IbBat_Proc > NOccBat_per_Proc_half) THEN
                  IaBatBg = MyRankMO * NOccBat_per_Proc + NOccBat_per_Proc / 2
               END IF
            END IF
            IF (ExchIBat .AND. (MyRankMO > Jrankrecv)) THEN
               IF (IbBat_Proc <= NOccBat_per_Proc_half) THEN
                  IaBatBg = MyRankMO * NOccBat_per_Proc + NOccBat_per_Proc / 2
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
               LCount = 0
               DO Ib_0 = IbBg, IbEd, nBlk
                  lenB = IbEd - Ib_0 + 1
                  if ( lenB > nBlk ) lenB = nBlk
                  J = Ib_0 - IbBg + 1

                  IF (IaBat == IbBat) IaEd = Ib_0 + lenB - 1

                  DO Ia_0 = IaBg, IaEd, nBlk
                     lenA = IaEd - Ia_0 + 1
                     if ( lenA > nBlk ) lenA = nBlk
                     I = Ia_0 - IaBg_Proc

                     LCount = LCount + 1
                     Ib_0_List(LCount) = Ib_0
                     Ia_0_List(LCount) = Ia_0
                  ENDDO
               ENDDO
               cflag_set_mat_A(:) = 0
               cflag_set_mat_B(:) = 0
!
               ! if ( LCount < commCount ) then
               !    write(*,'(a,i6,a,i6,a,i6)') "[WARN] myRank=", myRank, ", LCount=", LCount, ", commCount=", commCount
               !    ! stop
               ! endif
!
#ifdef DEBUG
               if ( MyRank == 0 ) then
                  write(*,'(a,i6)') "# LCount:", LCount
                  write(*,'("# ",i6,":",i6,"-",i6,"(",i3,")",",",i6,"-",i6,"(",i3,")")') MyRank, &
                       IbBg, IbEd, IbEd - IbBg + 1, &
                       IaBg, IaEd, IaEd - IaBg + 1
               endif
#endif

               WTimeBgn = MPI_WTIME()
               CALL CPU_TIME(TimeBgn)
               commPhase = 1
               if ( commSizeEach(commPhase) > 0 ) then
                  CALL MPI_ISend(SendBuf(1,commIndexEach(commPhase)), commSizeEach(commPhase), &
                       MPI_DOUBLE_PRECISION, Jranksend_1, commPhase, MPI_COMM_MO, ireq(1), IErr)
                  CALL MPI_IRecv(RecvBuf(1,commIndexEach(commPhase)), commSizeEach(commPhase), &
                       MPI_DOUBLE_PRECISION, Jrankrecv_1, commPhase, MPI_COMM_MO, ireq(2), IErr)
               endif
               CALL CPU_TIME(TimeEnd)
               WTimeEnd = MPI_WTIME()
               Time_T3C = Time_T3C + TimeEnd - TimeBgn
               WTime_T3C = WTime_T3C + WTimeEnd - WTimeBgn

               !
               !
               DO LNumber_Base = 1, LCount + (NUM_STREAM-1)

                  !
                  ! 1st half
                  !
                  LNumber = LNumber_Base
                  if ( LNumber >= 1 .and. LNumber <= LCount ) then

                     Ib_0 = Ib_0_List(LNumber)
                     lenB = IbEd - Ib_0 + 1
                     if ( lenB > nBlk ) lenB = nBlk
                     J = Ib_0 - IbBg + 1
                     ! test
                     id_B  = (Ib_0-IbBg)/nBlk + 1
                     id_Bm = mod(id_B - 1, NUM_STREAM) + 1

                     IF (IaBat == IbBat) IaEd = Ib_0 + lenB - 1

                     Ia_0 = Ia_0_List(LNumber)
                     lenA = IaEd - Ia_0 + 1
                     if ( lenA > nBlk ) lenA = nBlk
                     I = Ia_0 - IaBg_Proc
                     ! test
                     id_A  = (Ia_0-IaBg)/nBlk + 1
                     id_Am = mod(id_A - 1, NUM_STREAM) + 1

                     m = NActO(1)*lenA
                     n = NActO(1)*lenB
                     k = NBF_RI_per_ProcMat

                     id_st = mod(LNumber-1, NUM_STREAM) + 1
                     RWork2 => RWork2_Pool(:,:,id_st)

                     id_ev = mod(LNumber-1, NUM_EVENT) + 1
                     id_ev_next = mod(id_ev, NUM_EVENT) + 1
#ifdef DEBUG
                     if ( MyRank == 0 ) then
                        write(*,'(a,i4,a,i4,a,i6,a,i4,a,i6,a,i4,a,3i6)') &
                             "# H1: LNumber:", LNumber, ", id_st:", id_st, &
                             ", Ib_0:", Ib_0, ", lenB:", lenB, &
                             ", Ia_0:", Ia_0, ", lenA:", lenA, &
                             ", (m,n,k):", m, n, k
                        write(*,'(a,i3,a,i3,a,i3,a,i3)') "# id_B=", id_B, ", id_A=", id_A, &
                             ", id_ev=", id_ev, ", id_ev_next=", id_ev_next
                     endif
#endif
!
!                    o Evaluation of four-center MO integrals (ia|jb) from three-center integrals
!
                     WTimeBgn = MPI_WTIME()
                     CALL CPU_TIME(TimeBgn)
                     ! DGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC )
#ifdef USE_GPU
                     ! wait for previous send to complete
                     call cublas_ev_wait( id_ev, id_st )

                     ! send matrix B to GPU
                     if ( cflag_set_mat_B(id_B) == 0 ) then
                        call cublas_set_matrix_async( k, n, MOInt3ck (1,J), ldb, devptr_B(id_Bm), ldb, id_st )
                        cflag_set_mat_B(id_B) = 1
                     endif
                     ! send matrix A to GPU
                     if ( nGrp <= NUM_STREAM ) then
                        if ( cflag_set_mat_A(id_A) == 0 ) then
                           call cublas_set_matrix_async( k, m, RIInt3c3a(1,I), lda, devptr_A(id_Am), lda, id_st )
                           cflag_set_mat_A(id_A) = 1
                        endif
                     else
                        id_Am = id_st
                        call cublas_set_matrix_async( k, m, RIInt3c3a(1,I), lda, devptr_A(id_Am), lda, id_st )
                     endif

                     !
                     call cublas_ev_rec( id_ev_next, id_st )

                     ! wait for previous dgemm to complete
                     call cublas_ev2_wait( id_ev, id_st )

                     ! run dgemm on GPU
                     call cublas_dgemm_async('T', 'N', m, n, k, One, &
                          devptr_A(id_Am), lda, &
                          devptr_B(id_Bm), ldb, Zero, &
                          devptr_C(id_st), ldc, id_st )

                     !
                     call cublas_ev2_rec( id_ev_next, id_st )

                     ! get matrix C from GPU
                     call cublas_get_matrix_async( m, n, devptr_C(id_st), ldc, RWork2, ldc, id_st )

                     ! debug
                     ! call cublas_st_sync( id_st )
#else
                     CALL DGEMM('T', 'N', m, n, k, One, &
                          RIInt3c3a(1,I), lda, &
                          MOInt3ck (1,J), ldb, Zero, &
                          RWork2, ldc )
#endif
                     count_dgemm = count_dgemm + 1  ! debug
                     CALL CPU_TIME(TimeEnd)
                     WTimeEnd = MPI_WTIME()
                     Time_MOI = Time_MOI + TimeEnd - TimeBgn
                     WTime_MOI = WTime_MOI + WTimeEnd - WTimeBgn
                  endif
                  !
                  ! end of 1st half
                  !

                  !
                  ! 2nd half
                  !
                  LNumber = LNumber_Base - (NUM_STREAM-1)

                  if ( LNumber >= 1 .and. LNumber <= LCount ) then

                     Ib_0 = Ib_0_List(LNumber)
                     lenB = IbEd - Ib_0 + 1
                     if ( lenB > nBlk ) lenB = nBlk
                     J = Ib_0 - IbBg + 1

                     IF (IaBat == IbBat) IaEd = Ib_0 + lenB - 1

                     Ia_0 = Ia_0_List(LNumber)
                     lenA = IaEd - Ia_0 + 1
                     if ( lenA > nBlk ) lenA = nBlk
                     I = Ia_0 - IaBg_Proc

                     m = NActO(1)*lenA
                     n = NActO(1)*lenB
                     k = NBF_RI_per_ProcMat

                     id_st = mod(LNumber-1, NUM_STREAM) + 1
                     RWork2 => RWork2_Pool(:,:,id_st)
#ifdef DEBUG
                     if ( MyRank == 0 ) then
                        write(*,'(a,i4,a,i4,a,i6,a,i4,a,i6,a,i4,a,3i6)') &
                             "# H2: LNumber:", LNumber, ", id_st:", id_st, &
                             ", Ib_0:", Ib_0, ", lenB:", lenB, &
                             ", Ia_0:", Ia_0, ", lenA:", lenA, &
                             ", (m,n,k):", m, n, k
                     endif
#endif
                     commPhase = LNumber
                     if ( commPhase <= commCount .and. commSizeEach(commPhase) > 0 ) then
                        CALL MPI_Wait(ireq(1), istat1, IErr)
                        CALL MPI_Wait(ireq(2), istat2, IErr)
                     endif

                     WTimeBgn = MPI_WTIME()
                     CALL CPU_TIME(TimeBgn)
                     commPhase = LNumber + 1
                     if ( commPhase <= commCount .and. commSizeEach(commPhase) > 0 ) then
                        CALL MPI_ISend(SendBuf(1,commIndexEach(commPhase)), commSizeEach(commPhase), &
                             MPI_DOUBLE_PRECISION, Jranksend_1, commPhase, MPI_COMM_MO, ireq(1), IErr)
                        CALL MPI_IRecv(RecvBuf(1,commIndexEach(commPhase)), commSizeEach(commPhase), &
                             MPI_DOUBLE_PRECISION, Jrankrecv_1, commPhase, MPI_COMM_MO, ireq(2), IErr)
                     endif
                     CALL CPU_TIME(TimeEnd)
                     WTimeEnd = MPI_WTIME()
                     Time_T3C = Time_T3C + TimeEnd - TimeBgn
                     WTime_T3C = WTime_T3C + WTimeEnd - WTimeBgn

#ifdef USE_GPU
                     call cublas_st_sync( id_st )
#endif
                     !
                     ! do reductoin sum for dgemm result among processes in a MPI_COMM_MAT group
                     ! typical number of process in a MPI_COMM_MAT is upto 8 (might be more in case of 10K processes)
                     !
                     WTimeBgn = MPI_WTIME()
                     CALL CPU_TIME(TimeBgn)
                     CALL MPI_Allreduce(RWork2, MOIntSome, ldc*n, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_MAT, IErr)
                     CALL CPU_TIME(TimeEnd)
                     WTimeEnd = MPI_WTIME()
                     Time_MOIC = Time_MOIC + TimeEnd - TimeBgn
                     WTime_MOIC = WTime_MOIC + WTimeEnd - WTimeBgn

!
!                    o Evaluation of MP2 correlation energy for ij orbital pair
!
                     WTimeBgn = MPI_WTIME()
                     CALL CPU_TIME(TimeBgn)
                     !
                     do Ib_1 = 1, lenB
                        do Ia_1 = 1, lenA
                           Ib = Ib_0 + ib_1 - 1
                           Ia = Ia_0 + Ia_1 - 1
                           if ((IaBat == IbBat) .and. (Ia > Ib)) cycle

                           EigIb = - Eig(Ib+NHOMO)
                           EigIab = EigIb - Eig(Ia+NHOMO)
!
                           OfstJ = (Ib_1 - 1)*NActO(1)
                           OfstI = (Ia_1 - 1)*NActO(1)
!
                           E2Tab = Zero
                           E2Sab = Zero

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(Ii, Ij, EigIjab, EigIi, T2)
!$OMP DO REDUCTION(+: E2Tab, E2Sab)
                           DO Ij = 1, NActO(1)
                              EigIjab = EigIab + Eig(Ij+NFrzO(1))
                              DO Ii = IiBgn, IiEnd
                                 EigIi = Eig(Ii+NFrzO(1))
                                 T2 = MOIntSome(Ii+OfstI,Ij+OfstJ) / (EigIjab + EigIi)
                                 E2Tab = E2Tab + T2 * (MOIntSome(Ii+OfstI,Ij+OfstJ) - MOIntSome(Ij+OfstI,Ii+OfstJ))
                                 E2Sab = E2Sab + T2 *  MOIntSome(Ii+OfstI,Ij+OfstJ)
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
!
                        end DO
                     end DO
                     !
                     CALL CPU_TIME(TimeEnd)
                     WTimeEnd = MPI_WTIME()
                     Time_EMP2 = Time_EMP2 + TimeEnd - TimeBgn
                     WTime_EMP2 = WTime_EMP2 + WTimeEnd - WTimeBgn
                     !
                  end if
                  ! end of 2nd half

               END DO

               !
               !  when commCount > LCount, following MPI reqs are necessary ...
               !
               do LNumber = LCount+1, commCount
                     commPhase = LNumber
                     if ( commPhase <= commCount .and. commSizeEach(commPhase) > 0 ) then
                        CALL MPI_Wait(ireq(1), istat1, IErr)
                        CALL MPI_Wait(ireq(2), istat2, IErr)
                     endif

                     WTimeBgn = MPI_WTIME()
                     CALL CPU_TIME(TimeBgn)
                     commPhase = LNumber + 1
                     if ( commPhase <= commCount .and. commSizeEach(commPhase) > 0 ) then
                        CALL MPI_ISend(SendBuf(1,commIndexEach(commPhase)), commSizeEach(commPhase), &
                             MPI_DOUBLE_PRECISION, Jranksend_1, commPhase, MPI_COMM_MO, ireq(1), IErr)
                        CALL MPI_IRecv(RecvBuf(1,commIndexEach(commPhase)), commSizeEach(commPhase), &
                             MPI_DOUBLE_PRECISION, Jrankrecv_1, commPhase, MPI_COMM_MO, ireq(2), IErr)
                     endif
                     CALL CPU_TIME(TimeEnd)
                     WTimeEnd = MPI_WTIME()
                     Time_T3C = Time_T3C + TimeEnd - TimeBgn
                     WTime_T3C = WTime_T3C + WTimeEnd - WTimeBgn
               enddo

!
            END DO

         END DO
      END DO
!
!     o Memory deallocation
!
#ifdef USE_GPU
      !
      ! finalize cublas
      !
      CALL CPU_TIME(TimeBgn)
      call mpi_barrier(MPI_COMM_WORLD,ierr)  ! not essential, just make measured time meaningful
      ! free buffer on GPU
      do id_st = 1, NUM_STREAM
         call cublas_free( devptr_A(id_st) )
         call cublas_free( devptr_B(id_st) )
         call cublas_free( devptr_C(id_st) )
      enddo
      call cublas_fin()
      call mpi_barrier(MPI_COMM_WORLD,ierr)  ! not essential, just make measured time meaningful
      CALL CPU_TIME(TimeEnd)
#ifdef DEBUG
      if ( MyRank == 0 ) then
         write(*,'("# gpu free/fin time",F12.2)') TimeEnd - TimeBgn
         write(*,*) "# count_dgemm:", count_dgemm
      end if
#endif
#endif

      ! DEALLOCATE(MOInt3ck)
      DEALLOCATE(CommBuf)
      DEALLOCATE(MOIntSome)
      DEALLOCATE(Eig)
      DEALLOCATE(istat1)
      DEALLOCATE(istat2)
      DEALLOCATE(istat3)
      DEALLOCATE(RWork2_Pool)
      DEALLOCATE(Ia_0_list, Ib_0_list)
      DEALLOCATE(cflag_set_mat_A, cflag_set_mat_B)
      DEALLOCATE(devptr_A, devptr_B, devptr_C)

!
!     o Correct MP2 correlation energy
!
      WTimeBgn = MPI_WTIME()
      CALL CPU_TIME(TimeBgn)
      CALL MPI_Allreduce(E2SP, E2S, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, IErr)
      CALL MPI_Allreduce(E2TP, E2T, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, IErr)
      CALL CPU_TIME(TimeEnd)
      WTimeEnd = MPI_WTIME()
      Time_EMP2C = Time_EMP2C + TimeEnd - TimeBgn
      WTime_EMP2C = WTime_EMP2C + WTimeEnd - WTimeBgn
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
         PRINT '(" ..... CPU time (3/3k 2cints comm) :", F12.2)', Time_T3C
         PRINT '(" ..... CPU time (4c Ints         ) :", F12.2)', Time_MOI
         PRINT '(" ..... CPU time (4c Ints comm    ) :", F12.2)', Time_MOIC
         PRINT '(" ..... CPU time (EMP2 corr.      ) :", F12.2)', Time_EMP2
         PRINT '(" ..... CPU time (EMP2 corr. comm ) :", F12.2)', Time_EMP2C
         PRINT '(" ..... WALL time (3/3k 2cints comm) :", F12.2)', WTime_T3C
         PRINT '(" ..... WALL time (4c Ints         ) :", F12.2)', WTime_MOI
         PRINT '(" ..... WALL time (4c Ints comm    ) :", F12.2)', WTime_MOIC
         PRINT '(" ..... WALL time (EMP2 corr.      ) :", F12.2)', WTime_EMP2
         PRINT '(" ..... WALL time (EMP2 corr. comm ) :", F12.2)', WTime_EMP2C
!
      END IF
!
      END SUBROUTINE
