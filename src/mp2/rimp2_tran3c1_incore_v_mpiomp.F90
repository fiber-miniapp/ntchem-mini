      SUBROUTINE RIMP2_Tran3c1_InCore_V_MPIOMP
!
      USE MP2_Module, ONLY : IOccBat, NOccBat, PScreen, ThrPre, SchwInt, NBF, NMO, NActO, NActV, NFrzO, &
     &   MaxContS, Name, IPrint
      USE RIMP2_Module, ONLY : SchwInt_RI, RIInt3c2a, NBF_RI, NBF_RI_MyRank, IdxBF_RI_MyRank
      USE MP2_Basis_Module, ONLY : Expnt, KontG, Spherical, CCoef, KStart, KAtom, NShel, Centr, KType, &
     &   LtuvMin_Car, LtuvMin_Sph, LtuvMax_Car, LtuvMax_Sph
      USE RIMP2_Basis_Module, ONLY : Expnt_RI, KontG_RI, CCoef_RI, KStart_RI, KAtom_RI, NShel_RI, KType_RI, KLoc_RI_Car, KLoc_RI_Sph
      USE Int2_Module, ONLY : ExpntA, ExpntB, ExpntC, ExpntD, ExpntP, ExpntQ, &
     &   PX, PY, PZ, QX, QY, QZ, PAX, PAY, PAZ, PBX, PBY, PBZ, QCX, QCY, QCZ, QDX, QDY, QDZ, &
     &   CContAB, CContCD, CCoefAB, CCoefCD, PreFactAB, PreFactCD, IAnglA, IAnglB, IAnglC, IAnglD, DoPH, ThrPrim
      USE Int2_ECoef_Module, ONLY : ECoefXAB, ECoefYAB, ECoefZAB, ECoefXCD, ECoefYCD, ECoefZCD, ExpPHalf, ExpQHalf
      USE Int2_Int2e_Module, ONLY : PQX, PQY, PQZ, ExpntPQ1, ExpntPQ2
      USE Int2_Gamma_Module, ONLY : FF, MaxtuvGam
      USE MP2_Constant_Module, ONLY : Zero, One, Half, Pi252, RLN10
      USE MPI_Module, ONLY : NProcs, MyRank, MPIIO, IORank, MPI_COMM_IO, MPI_COMM_MO, NProcsMO, MyRankMO, &
     &   MPI_COMM_Mat, NProcsMat, MyRankMat
!
!     o 3c RI integral transformation from AO to MO basis
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
      INTEGER :: II, JJ, KK, IJ, KL, I, J, K, ITemp, JTemp, KTemp
      INTEGER :: IAtomA, IAtomB, IAtomC, IAtomD
      INTEGER :: IPrim1, JPrim1, KPrim1, IPrim2, JPrim2, KPrim2, JPTemp2
      INTEGER :: IIOff
      INTEGER :: Mod_NBF_NProcsMat
      INTEGER :: lenj, jdim, jbgn
      REAL(8) :: ExpA, ExpB, ExpC, ExpD, ExpP, ExpQ, ExpPI, ExpQI, ExpAR2, ExpCR2, ExpnPQ
      REAL(8) :: ACentX, ACentY, ACentZ, BCentX, BCentY, BCentZ
      REAL(8) :: CCentX, CCentY, CCentZ, DCentX, DCentY, DCentZ
      REAL(8) :: ABCentX, ABCentY, ABCentZ, R2AB, R2CD
      REAL(8) :: ThrFac
!
!     o E-coefficients
!
      INTEGER :: IJPrim, KLPrim
      INTEGER :: NPIJ, NPKL
      REAL(8) :: ExpKAB, ExpKCD
!
      INTEGER :: KBF_RI, LabelK, NK, MAX_NK
      INTEGER :: MActO, MActV
!
!     o Integral prescreening
!
      REAL(8) :: SchwAB, SchwCD
!
#ifdef USE_GPU
      REAL(8), ALLOCATABLE, pinned :: CMO(:,:,:), ERIMat(:,:), T1Int(:,:)
      REAL(8), ALLOCATABLE, pinned :: RWork1(:,:), RWork2(:,:)
#else
      REAL(8), ALLOCATABLE :: CMO(:,:,:), ERIMat(:,:), T1Int(:,:)
      REAL(8), ALLOCATABLE :: RWork1(:,:), RWork2(:,:)
#endif
      integer :: devptr_CMO1, devptr_CMO2
      integer, allocatable :: devptr_ERIMAT(:), devptr_T1Int(:), devptr_RWork1(:)
!
      REAL(8) :: TimeBgn, TimeEnd, Time_T0, Time_T1, Time_T2, WTimeBgn, WTimeEnd, WTime_T0, WTime_T1, WTime_T2
      REAL(8) :: Time_T0C, Time_T2C, WTime_T0C, WTime_T2C
!
      INTEGER :: IRank
      INTEGER :: IErr
      integer, parameter :: Num_Stream = 3
      integer :: id_st
      integer :: lm, ln, lk
!
      Time_T0 = Zero
      Time_T1 = Zero
      Time_T2 = Zero
      Time_T0C = Zero
      Time_T2C = Zero
      WTime_T0 = Zero
      WTime_T1 = Zero
      WTime_T2 = Zero
      WTime_T0C = Zero
      WTime_T2C = Zero
!
!     o Initialization
!
      ThrFac = RLN10 * LOG10(ThrPrim)
      ! MActO = MAX(NActO(1), NActO(2))
      MActO = NActO(1)
      MActV = MAX(NActV(1), NActV(2))
!
      lenj = NBF / NProcsMat
      Mod_NBF_NProcsMat = MOD(NBF, NProcsMat)
      if (Mod_NBF_NProcsMat > MyRankMat) then
         jbgn = lenj * MyRankMat + MyRankMat + 1
         jdim = lenj + 1
      else
         jbgn = lenj * MyRankMat + Mod_NBF_NProcsMat + 1
         jdim = lenj
      end if
!debug
!      if (myrank == 0) then
!         write(*, *) 'MyRankMat=', MyRankMat
!         write(*, *) 'jbgn=', jbgn, ' jdim=', jdim
!      end if
!
!     o Allocate memory
!
      ALLOCATE(CMO(NBF,NMO,1))
!MPI
      ALLOCATE(NBF_RI_MyRank(0:NProcsMO-1))
      ALLOCATE(IdxBF_RI_MyRank(NBF_RI))
!MPI
!
!     o Get infomation of parallel distribution of aux. basis
!
      MAX_NK = 0  ! test
      DO IRank = 0, (NProcsMO - 1)
         KBF_RI = 0
         DO KK = 1, NShel_RI
            IF (MOD(KK, NProcsMO) /= IRank) CYCLE
            IAnglC = KType_RI(KK)
            IF (Spherical) THEN
               NK = LtuvMax_Sph(IAnglC) - LtuvMin_Sph(IAnglC) + 1
            ELSE
               NK = LtuvMax_Car(IAnglC) - LtuvMin_Car(IAnglC) + 1
            END IF
            KBF_RI = KBF_RI + NK
            if (MAX_NK < NK) MAX_NK = NK  ! test
         END DO
         NBF_RI_MyRank(IRank) = KBF_RI
#ifdef DEBUG
         if ( MyRank == 0 ) then
            write(*,'(a,i4,a,i6)') "# NBF_RI_MyRank(", IRank, ") = ", KBF_RI
         endif
#endif
      END DO
#ifdef DEBUG
      if ( MyRank == 0 ) then
         write(*,'(a,i4)') "# MAX_NK = ", MAX_NK
         write(*,'(a,i4)') "# NBF = ", NBF
         write(*,'(a,i4)') "# NMO = ", NMO
      endif
#endif
!
!     o Read MO coefficients
!
      IF (MPIIO) THEN
         FBuf = TRIM(Name)//".MO"
         OPEN(UNIT=IO, FILE=TRIM(FBuf), STATUS='OLD', ACCESS='SEQUENTIAL', FORM='FORMATTED')
         READ(IO, *) CMO(1:NBF,1:NMO,1)
         CLOSE(IO)
      END IF
      CALL MPI_Bcast(CMO, NBF*NMO, MPI_DOUBLE_PRECISION, IORank, MPI_COMM_IO, IErr)
!
      IF ((MyRank == 0) .AND. (IPrint >= 2)) THEN
         WRITE(*, *) '+++++ MO coefficient (Alpha) +++++'
         CALL Util_MatOut(CMO(1,1,1), NBF, NMO)
      END IF
!
! anaruse, test
      ALLOCATE(ERIMat(NBF*NBF,MAX_NK))
      ALLOCATE(T1Int(MActO*NBF,MAX_NK))
      ! ALLOCATE(RWork1(NActO(1)*NActV(1),MAX_NK))
      ALLOCATE(RWork2(NBF*NBF,MAX_NK))

#ifdef USE_GPU
      allocate(devptr_ERIMAT(MAX_NK), devptr_T1Int(MAX_NK), devptr_RWork1(MAX_NK))
      call cublas_init()
      !
      call cublas_alloc( devptr_CMO1, NBF, NMO );
      call cublas_alloc( devptr_CMO2, NBF, NMO );
      do k = 1, MAX_NK
         call cublas_alloc( devptr_ERIMAT(k), NBF, NBF );
         call cublas_alloc( devptr_T1Int(k),  MActO, NBF );
         call cublas_alloc( devptr_Rwork1(k), NActO(1), NActV(1) );
      end do
      !
      call cublas_set_matrix( NBF, NActO(1), CMO(1,(NFrzO(1)+1),1), NBF, devptr_CMO1, NBF );
      call cublas_set_matrix( NBF, NActV(1), CMO(1,(NFrzO(1)+NActO(1)+1),1), NBF, devptr_CMO2, NBF );
#endif

!
!$omp parallel
!
!     o Allocate memory for 3c AO integrals
!
      CALL Int2_Allocate
      CALL Int2_ECoef_Allocate
      CALL Int2_Array_Allocate
      CALL Int2_Int2e_Allocate
      ALLOCATE(FF(MaxtuvGam))
!$omp end parallel
!
      KBF_RI = 0
!MPI parallel
      DO KK = 1, NShel_RI
         IF (MOD(KK, NProcsMO) /= MyRankMO) CYCLE
!MPI parallel

!#ifdef DEBUG
#if 1
         if ( MyRank == 0 ) then
            write(*,'(a,i6,a,i6)') "# KK = ", KK, ", NShel_RI = ", NShel_RI
         endif
#endif

         WTimeBgn = MPI_WTIME()
         CALL CPU_TIME(TimeBgn)
!
         IAnglC = KType_RI(KK)
         IF (Spherical) THEN
            LabelK = KLoc_RI_Sph(KK) 
            NK = LtuvMax_Sph(IAnglC) - LtuvMin_Sph(IAnglC) + 1
         ELSE
            LabelK = KLoc_RI_Car(KK)
            NK = LtuvMax_Car(IAnglC) - LtuvMin_Car(IAnglC) + 1
         END IF
         ! ALLOCATE(RWork2(NBF*NBF,NK))
         CALL DCOPY(NBF*NBF*NK, Zero, 0, RWork2, 1)
!
!        o Evaluation of three-center RI integrals (PQ|C)
!
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP&PRIVATE(II, IAtomA, ACentX, ACentY, ACentZ, IPrim1, IPrim2, ITemp, I, IIOff, &
!$OMP&        JJ, IAtomB, BCentX, BCentY, BCentZ, JPrim1, JPrim2, JTemp, J, JPTemp2, &
!$OMP&            IAtomC, CCentX, CCentY, CCentZ, KPrim1, KPrim2, KTemp, K, &
!$OMP&            IAtomD, DCentX, DCentY, DCentZ, &
!$OMP&        IJPrim, KLPrim, IJ, KL, NPIJ, NPKL, &
!$OMP&        ABCentX, ABCentY, ABCentZ, R2AB, R2CD, SchwAB, SchwCD, &
!$OMP&        ExpA, ExpB, ExpC, ExpD, ExpP, ExpQ, ExpPI, ExpQI, ExpAR2, ExpCR2, ExpnPQ, &
!$OMP&        ExpKAB, ExpKCD)
!
!        o Allocate memory for 3c AO integrals
!
         ! CALL Int2_Allocate
         ! CALL Int2_ECoef_Allocate
         ! CALL Int2_Array_Allocate
         ! CALL Int2_Int2e_Allocate
         ! ALLOCATE(FF(MaxtuvGam))
!
         IAtomC = KAtom_RI(KK)
         IAnglC = KType_RI(KK)
         CCentX = Centr(1,IAtomC)
         CCentY = Centr(2,IAtomC)
         CCentZ = Centr(3,IAtomC)
         KPrim1 = KStart_RI(KK)
         KPrim2 = KPrim1 + KontG_RI(KK) - 1
         ! KTemp = 0
         ! DO K = KPrim1, KPrim2
         !    KTemp = KTemp + 1
         !    ExpntC(KTemp) = Expnt_RI(K)
         ! END DO
!
         IAtomD = IAtomC
         IAnglD = 0   ! s-type
         DCentX = Centr(1,IAtomD)
         DCentY = Centr(2,IAtomD)
         DCentZ = Centr(3,IAtomD)
         ExpntD(1) = Zero
!
         R2CD = Zero
         KLPrim = 0
         DO K = KPrim1, KPrim2
            ExpC = Expnt_RI(K)
            ExpCR2 = ExpC * R2CD
            ExpD = ExpntD(1)
            ExpQ = ExpC + ExpD
            ExpQI = One / ExpQ
            ExpKCD = - ExpCR2 * ExpD * ExpQI
            KLPrim = KLPrim + 1
            PreFactCD(KLPrim) = EXP(ExpKCD)
            ExpntQ(KLPrim) = ExpQ
            ExpQHalf(KLPrim) = Half * ExpQI
            CCoefCD(KLPrim) = CCoef_RI(K)
            QX(KLPrim) = (ExpC * CCentX + ExpD * DCentX) * ExpQI
            QY(KLPrim) = (ExpC * CCentY + ExpD * DCentY) * ExpQI
            QZ(KLPrim) = (ExpC * CCentZ + ExpD * DCentZ) * ExpQI
            QCX(KLPrim) = QX(KLPrim) - CCentX
            QCY(KLPrim) = QY(KLPrim) - CCentY
            QCZ(KLPrim) = QZ(KLPrim) - CCentZ
            QDX(KLPrim) = QX(KLPrim) - DCentX
            QDY(KLPrim) = QY(KLPrim) - DCentY
            QDZ(KLPrim) = QZ(KLPrim) - DCentZ
         END DO
         NPKL = KLPrim

         IF (PScreen) THEN
            SchwCD = SchwInt_RI(KK)
         END IF
!
!        o Normalization
!
         CALL RIInt2_MDInt2_CCont(CCoefCD, CContCD, IAnglC, IAnglD, NPKL)
!
!        o Generate E-coefficients for C and D
!
         CALL MDInt2_ECoef1(ECoefXCD, ECoefYCD, ECoefZCD, QCX, QCY, QCZ, QDX, QDY, QDZ, ExpQHalf, PreFactCD, &
     &      IAnglC, IAnglD, NPKL)
!
!        o calculate 3c AO integrals (PQ|C)
!
!$OMP DO SCHEDULE(DYNAMIC, 1)
         DO II = NShel, 1, -1
!         DO II = 1, NShel
            ! debug
            ! if ( MyRank == 0 ) then
            !    write(*,'(a,2i8)') "II:", II, NShel
            ! endif
            
            IAtomA = KAtom(II)
            IAnglA = KType(II)
            ACentX = Centr(1,IAtomA)
            ACentY = Centr(2,IAtomA)
            ACentZ = Centr(3,IAtomA)
            IPrim1 = KStart(II)
            IPrim2 = IPrim1 + KontG(II) - 1
            IIOff = (II * (II - 1)) / 2
!
!MPI parallel
            DO JJ = 1, II
               IF (MOD(JJ, NProcsMat) /= MyRankMat) CYCLE
!MPI parallel
!
!              o Check integrals by the Schwarz inequality
!
               IF (PScreen) THEN
                  SchwAB = SchwInt(IIOff+JJ)
                  IF ((SchwAB * SchwCD) < ThrPre) CYCLE
               END IF
!
               IAtomB = KAtom(JJ)
               IAnglB = KType(JJ)
               BCentX = Centr(1,IAtomB)
               BCentY = Centr(2,IAtomB)
               BCentZ = Centr(3,IAtomB)
               JPrim1 = KStart(JJ)
               JPrim2 = JPrim1 + KontG(JJ) - 1
!
               ABCentX = ACentX - BCentX
               ABCentY = ACentY - BCentY
               ABCentZ = ACentZ - BCentZ
               R2AB = ABCentX * ABCentX + ABCentY * ABCentY + ABCentZ * ABCentZ
               IJPrim = 0
               DO I = IPrim1, IPrim2
                  ExpA = Expnt(I)
                  ExpAR2 = ExpA * R2AB
                  JPTemp2 = JPrim2
                  IF (II == JJ) JPTemp2 = I
                  DO J = JPrim1, JPTemp2
                     ExpB = Expnt(J)
                     ExpP = ExpA + ExpB
                     ExpPI = One / ExpP
                     ExpKAB = - ExpAR2 * ExpB * ExpPI
                     IF (ExpKAB >= ThrFac) THEN
                        IJPrim = IJPrim + 1
                        PreFactAB(IJPrim) = EXP(ExpKAB)
                        ExpntP(IJPrim) = ExpP
                        ExpPHalf(IJPrim) = Half * ExpPI
                        CCoefAB(IJPrim) = CCoef(I) * CCoef(J)
                        IF ((II == JJ) .AND. (I /= J)) CCoefAB(IJPrim) = CCoefAB(IJPrim) + CCoefAB(IJPrim)
                        PX(IJPrim) = (ExpA * ACentX + ExpB * BCentX) * ExpPI
                        PY(IJPrim) = (ExpA * ACentY + ExpB * BCentY) * ExpPI
                        PZ(IJPrim) = (ExpA * ACentZ + ExpB * BCentZ) * ExpPI
                        PAX(IJPrim) = PX(IJPrim) - ACentX
                        PAY(IJPrim) = PY(IJPrim) - ACentY
                        PAZ(IJPrim) = PZ(IJPrim) - ACentZ
                        PBX(IJPrim) = PX(IJPrim) - BCentX
                        PBY(IJPrim) = PY(IJPrim) - BCentY
                        PBZ(IJPrim) = PZ(IJPrim) - BCentZ
                     END IF
                  END DO
               END DO
               NPIJ = IJPrim
               IF (NPIJ == 0) CYCLE
!
!              o Normalization
!
               CALL Int2_CCont(CCoefAB, CContAB, II, JJ, IAnglA, IAnglB, NPIJ)
!
!              o Generate E-coefficients for A and B
!
               IF (IAnglA >= IAnglB) THEN
                  CALL MDInt2_ECoef1(ECoefXAB, ECoefYAB, ECoefZAB, PAX, PAY, PAZ, PBX, PBY, PBZ, ExpPHalf, PreFactAB, &
     &               IAnglA, IAnglB, NPIJ)
               ELSE
                  CALL MDInt2_ECoef2(ECoefXAB, ECoefYAB, ECoefZAB, PAX, PAY, PAZ, PBX, PBY, PBZ, ExpPHalf, PreFactAB, &
     &               IAnglA, IAnglB, NPIJ)
               END IF
!
!              o Generate R-integrals
!
#if 0
               DO IJ = 1, NPIJ
                  DO KL = 1, NPKL
                     ExpnPQ = ExpntP(IJ) + ExpntQ(KL)
                     ExpntPQ1(KL,IJ) = ExpntP(IJ) * ExpntQ(KL) / ExpnPQ
                     ExpnPQ = ExpntP(IJ) * ExpntQ(KL) * SQRT(ExpnPQ)
                     ExpntPQ2(KL,IJ) = Pi252 / ExpnPQ   ! Adsorption
                     PQX(KL,IJ) = PX(IJ) - QX(KL)
                     PQY(KL,IJ) = PY(IJ) - QY(KL)
                     PQZ(KL,IJ) = PZ(IJ) - QZ(KL)
                  END DO
               END DO
               CALL MDInt2_R0tuv(NPIJ, NPKL)
#else
               CALL MDInt2_R0tuv_test(NPIJ, NPKL)
#endif
!
!              o Construct three-center two-electron ERI matrix
!
               IF (Spherical) THEN
                  CALL RIMP2_RIInt2_MDInt2_ERI3c_Sph(II, JJ, KK, NPIJ, NPKL, RWork2)
               ELSE
                  CALL RIMP2_RIInt2_MDInt2_ERI3c_Car(II, JJ, KK, NPIJ, NPKL, RWork2)
               END IF
!
            END DO
         END DO
!$OMP END DO
!
!        o Deallocate memory for 3c AO integrals
!
         ! CALL Int2_Deallocate
         ! CALL Int2_ECoef_Deallocate
         ! CALL Int2_Int2e_Deallocate
         ! CALL Int2_Array_Deallocate
         ! DEALLOCATE(FF)
!$OMP END PARALLEL
!
         CALL CPU_TIME(TimeEnd)
         WTimeEnd = MPI_WTIME()
         Time_T0 = Time_T0 + TimeEnd - TimeBgn
         WTime_T0 = WTime_T0 + WTimeEnd - WTimeBgn
!
         IF (Spherical) THEN
            LabelK = KLoc_RI_Sph(KK) 
            NK = LtuvMax_Sph(IAnglC) - LtuvMin_Sph(IanglC) + 1
         ELSE
            LabelK = KLoc_RI_Car(KK)
            NK = LtuvMax_Car(IAnglC) - LtuvMin_Car(IAnglC) + 1
         END IF
!
#ifdef USE_GPU
         DO id_st = 1, Num_Stream
            call cublas_st_sync( id_st )
         End DO
#endif

!
!        o Communicate 3c integral matrix
!
         ! WTimeBgn = MPI_WTIME()
         ! CALL CPU_TIME(TimeBgn)
         ! CALL MPI_Allreduce(RWork2(1,1), ERIMat(1,1), NBF*NBF*NK, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_MAT, IErr)
         ! CALL CPU_TIME(TimeEnd)
         ! WTimeEnd = MPI_WTIME()
         ! Time_T0C = Time_T0C + TimeEnd - TimeBgn
         ! WTime_T0C = WTime_T0C + WTimeEnd - WTimeBgn

         DO K = 1, NK
            KBF_RI = KBF_RI + 1
            LabelK = LabelK + 1
            IdxBF_RI_MyRank(KBF_RI) = LabelK
            id_st = mod(K-1, Num_Stream) + 1
!
!           o Communicate 3c integral matrix
!
            WTimeBgn = MPI_WTIME()
            CALL CPU_TIME(TimeBgn)
            CALL MPI_Allreduce(RWork2(1,K), ERIMat(1,K), NBF*NBF, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_MAT, IErr)
            CALL CPU_TIME(TimeEnd)
            WTimeEnd = MPI_WTIME()
            Time_T0C = Time_T0C + TimeEnd - TimeBgn
            WTime_T0C = WTime_T0C + WTimeEnd - WTimeBgn
!
!           o 1/3 integral transformation
!                (pq|c) -> (iq|c)
!
#ifdef DEBUG
            ! if ( MyRank == 0 ) then
            !    write(*,'(a,i4,a,3i6)') "  K=",K, ", 1/3, (m,n,k)=", NActO(1), NBF, NBF
            ! endif
#endif
            WTimeBgn = MPI_WTIME()
            CALL CPU_TIME(TimeBgn)
            lm = NActO(1)
            ln = NBF
            lk = NBF
#ifdef USE_GPU
            call cublas_set_matrix_async( lk, ln, ERIMat(1,K), lk, devptr_ERIMAT(K), lk, id_st )
            call cublas_dgemm_async('T', 'N', lm, ln, lk, One, &
                 devptr_CMO1, NBF, &
                 devptr_ERIMat(K), NBF, Zero, &
                 devptr_T1Int(K), NActO(1), id_st )
            ! call cublas_get_matrix_async( lm, ln, devptr_T1Int(K), lm, T1Int(1,K), lm, id_st )
#else
            CALL DGEMM('T', 'N', lm, ln, lk, One, &
                 CMO(1,(NFrzO(1)+1),1), NBF, &
                 ERIMat(1,K), NBF, Zero, &
                 T1Int(1,K), NActO(1))
#endif
            CALL CPU_TIME(TimeEnd)
            WTimeEnd = MPI_WTIME()
            Time_T1 = Time_T1 + TimeEnd - TimeBgn
            WTime_T1 = WTime_T1 + WTimeEnd - WTimeBgn
!
!           o 2/3 integral transformation
!               (iq|c) -> (ia|c)
!
#ifdef DEBUG
            ! if ( MyRank == 0 ) then
            !    write(*,'(a,i4,a,3i6)') "  K=",K, ", 2/3, (m,n,k)=", NActO(1), NActV(1), NBF
            ! endif
#endif
            WTimeBgn = MPI_WTIME()
            CALL CPU_TIME(TimeBgn)
            lm = NActO(1)
            ln = NActV(1)
            lk = NBF
#ifdef USE_GPU
            ! call cublas_set_matrix_async( lm, lk, T1Int(1,K), lm, devptr_T1Int(K), lm, id_st )
            call cublas_dgemm_async('N', 'N', lm, ln, lk, One, &
                 devptr_T1Int(K), NActO(1), &
                 devptr_CMO2, NBF, Zero, &
                 devptr_RWork1(K), NActO(1), id_st)
            call cublas_get_matrix_async( lm, ln, devptr_RWork1(K), lm, RIInt3c2a(1,KBF_RI), lm, id_st )
#else
            CALL DGEMM('N', 'N', lm, ln, lk, One, &
                 T1Int(1,K), NActO(1), &
                 CMO(1,(NFrzO(1)+NActO(1)+1),1), NBF, Zero, &
                 RIInt3c2a(1,KBF_RI), NActO(1))
#endif
            CALL CPU_TIME(TimeEnd)
            WTimeEnd = MPI_WTIME()
            Time_T2 = Time_T2 + TimeEnd - TimeBgn
            WTime_T2 = WTime_T2 + WTimeEnd - WTimeBgn
         END DO
!
      END DO

#ifdef USE_GPU
      DO id_st = 1, Num_Stream
         call cublas_st_sync( id_st )
      End DO
#endif

!$omp parallel
!
!     o Deallocate memory for 3c AO integrals
!
      CALL Int2_Deallocate
      CALL Int2_ECoef_Deallocate
      CALL Int2_Int2e_Deallocate
      CALL Int2_Array_Deallocate
      DEALLOCATE(FF)
!$omp end parallel

! anaruse, test
      DEALLOCATE(ERIMat)
      DEALLOCATE(T1Int)
      ! DEALLOCATE(RWork1)
      DEALLOCATE(RWork2)

#ifdef USE_GPU
      call cublas_free( devptr_CMO1 );
      call cublas_free( devptr_CMO2 );
      do k = 1, MAX_NK
         call cublas_free( devptr_ERIMAT(k) );
         call cublas_free( devptr_T1Int(k) );
         call cublas_free( devptr_Rwork1(k) );
      end do
      ! call cublas_fin()
      deallocate( devptr_ERIMAT, devptr_T1Int, devptr_RWork1 )
#endif
!
!     o Deallocate memory
!
      DEALLOCATE(CMO)
!
#ifdef DEBUG
      write(*,'(a,i6,F12.2)') "# ", MyRank, WTime_T0
#endif
      call MPI_Barrier( MPI_COMM_WORLD, IErr )
!
      IF (MyRank == 0) THEN
         PRINT '(" ..... CPU time (3c-RIInt        ) :", F12.2)', Time_T0
         PRINT '(" ..... CPU time (3c-RIInt comm   ) :", F12.2)', Time_T0C
         PRINT '(" ..... CPU time (1/3 tran3c1     ) :", F12.2)', Time_T1
         PRINT '(" ..... CPU time (2/3 tran3c1     ) :", F12.2)', Time_T2
         PRINT '(" ..... CPU time (2/3 tran3c1 comm) :", F12.2)', Time_T2C
         PRINT '(" ..... WALL time (3c-RIInt        ) :", F12.2)', WTime_T0
         PRINT '(" ..... WALL time (3c-RIInt comm   ) :", F12.2)', WTime_T0C
         PRINT '(" ..... WALL time (1/3 tran3c1     ) :", F12.2)', WTime_T1
         PRINT '(" ..... WALL time (2/3 tran3c1     ) :", F12.2)', WTime_T2
         PRINT '(" ..... WALL time (2/3 tran3c1 comm) :", F12.2)', WTime_T2C
      END IF
!
      END SUBROUTINE
