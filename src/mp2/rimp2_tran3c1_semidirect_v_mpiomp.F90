      SUBROUTINE RIMP2_Tran3c1_SemiDirect_V_MPIOMP
!
      USE MP2_Module, ONLY : IOccBat, NOccBat, NAlpBet, PScreen, ThrPre, SchwInt, NBF, NMO, NActO, NActV, NFrzO, &
     &   MaxContS, Name, IPrint
      USE RIMP2_Module, ONLY : SchwInt_RI, NBF_RI, NBF_RI_MyRank, IdxBF_RI_MyRank
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
      USE MPI_Module, ONLY : NProcs, MyRank, MPIIO, IORank, MPI_COMM_IO
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
      CHARACTER(LEN=1), PARAMETER :: AlpBet(1:2) = (/ "A", "B"/) 
      INTEGER, PARAMETER :: IO(1:2) = (/ 99, 98/)
      CHARACTER(LEN=255) :: FBuf
      INTEGER :: II, JJ, KK, IJ, KL, I, J, K, ITemp, JTemp, KTemp
      INTEGER :: IAtomA, IAtomB, IAtomC, IAtomD
      INTEGER :: IPrim1, JPrim1, KPrim1, IPrim2, JPrim2, KPrim2, JPTemp2
      INTEGER :: IIOff
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
      INTEGER :: KBF_RI, LabelK, NK
      INTEGER :: Ia
      INTEGER :: MActO, MActV
!
!     o Integral prescreening
!
      REAL(8) :: SchwAB, SchwCD
!
      REAL(8), ALLOCATABLE :: CMO(:,:,:), ERIMat(:,:), T1Int(:), T2Int(:)
!
      REAL(8) :: TimeBgn, TimeEnd, Time_T0, Time_T1, Time_T2, WTimeBgn, WTimeEnd, WTime_T0, WTime_T1, WTime_T2
      REAL(8) :: Time_T2W, WTime_T2W
!
      CHARACTER(LEN=10)  :: RankNo
      INTEGER :: IRank
      INTEGER :: IErr
!
      Time_T0 = Zero
      Time_T1 = Zero
      Time_T2 = Zero
      Time_T2W = Zero
      WTime_T0 = Zero
      WTime_T1 = Zero
      WTime_T2 = Zero
      WTime_T2W = Zero
!
!     o Initialization
!
      ThrFac = RLN10 * LOG10(ThrPrim)
      MActO = MAX(NActO(1), NActO(2))
      MActV = MAX(NActV(1), NActV(2))
!
!     o Allocate memory
!
      ALLOCATE(CMO(NBF,NMO,NAlpBet))
      ALLOCATE(ERIMat(NBF*NBF,MaxContS))
      ALLOCATE(T1Int(MActO*NBF))
      ALLOCATE(T2Int(MActO*MActV))
!MPI
      ALLOCATE(NBF_RI_MyRank(0:NProcs-1))
      ALLOCATE(IdxBF_RI_MyRank(NBF_RI))
!MPI
!
!     o Get infomation of parallel distribution of aux. basis
!
      DO IRank = 0, (NProcs - 1)
         KBF_RI = 0
         DO KK = 1, NShel_RI
            IF (MOD(KK, NProcs) /= IRank) CYCLE
            IAnglC = KType_RI(KK)
            IF (Spherical) THEN
               NK = LtuvMax_Sph(IAnglC) - LtuvMin_Sph(IAnglC) + 1
            ELSE
               NK = LtuvMax_Car(IAnglC) - LtuvMin_Car(IAnglC) + 1
            END IF
            KBF_RI = KBF_RI + NK
         END DO
         NBF_RI_MyRank(IRank) = KBF_RI
      END DO
!
!     o Read MO coefficients
!
      IF (MPIIO) THEN
         FBuf = TRIM(Name)//".MO"
         OPEN(UNIT=IO(1), FILE=TRIM(FBuf), STATUS='OLD', ACCESS='SEQUENTIAL', FORM='FORMATTED')
         READ(IO(1), *) CMO(1:NBF,1:NMO,1)
         CLOSE(IO(1))
      END IF
      CALL MPI_Bcast(CMO, NBF*NMO, MPI_DOUBLE_PRECISION, IORank, MPI_COMM_IO, IErr)
!
      IF ((MyRank == 0) .AND. (IPrint >= 2)) THEN
         WRITE(*, *) '+++++ MO coefficient (Alpha) +++++'
         CALL Util_MatOut(CMO(1,1,1), NBF, NMO)
      END IF
!
!     o Open three-center MO integral file
!
      WRITE(RankNo,'(I0)') MyRank
      FBuf = TRIM(Name)//".MOInt3ci"//AlpBet(1)//"."//TRIM(RankNo)
      OPEN(UNIT=IO(1), FILE=TRIM(FBuf), STATUS='UNKNOWN', ACCESS='DIRECT', FORM='UNFORMATTED', RECL=(8*NActO(1)))
!
      KBF_RI = 0
!MPI parallel
      DO KK = 1, NShel_RI
         IF (MOD(KK, NProcs) /= MyRank) CYCLE
!MPI parallel
!
         IAnglC = KType_RI(KK)
         IF (Spherical) THEN
            NK = LtuvMax_Sph(IAnglC) - LtuvMin_Sph(IAnglC) + 1
         ELSE
            NK = LtuvMax_Car(IAnglC) - LtuvMin_Car(IAnglC) + 1
         END IF
         CALL DCOPY(NBF*NBF*NK, Zero, 0, ERIMat, 1)
         WTimeBgn = MPI_WTIME()
         CALL CPU_TIME(TimeBgn)
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
         CALL Int2_Allocate
         CALL Int2_ECoef_Allocate
         CALL Int2_Array_Allocate
         CALL Int2_Int2e_Allocate
         ALLOCATE(FF(MaxtuvGam))
!
         IAtomC = KAtom_RI(KK)
         IAnglC = KType_RI(KK)
         CCentX = Centr(1,IAtomC)
         CCentY = Centr(2,IAtomC)
         CCentZ = Centr(3,IAtomC)
         KPrim1 = KStart_RI(KK)
         KPrim2 = KPrim1 + KontG_RI(KK) - 1
         KTemp = 0
         DO K = KPrim1, KPrim2
            KTemp = KTemp + 1
            ExpntC(KTemp) = Expnt_RI(K)
         END DO
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
            ExpC = ExpntC(K-KPrim1+1)
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
!
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
            IAtomA = KAtom(II)
            IAnglA = KType(II)
            ACentX = Centr(1,IAtomA)
            ACentY = Centr(2,IAtomA)
            ACentZ = Centr(3,IAtomA)
            IPrim1 = KStart(II)
            IPrim2 = IPrim1 + KontG(II) - 1
            ITemp = 0
            DO I = IPrim1, IPrim2
               ITemp = ITemp + 1
               ExpntA(ITemp) = Expnt(I)
            END DO
            IIOff = (II * (II - 1)) / 2
!
            DO JJ = 1, II
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
               JTemp = 0
               DO J = JPrim1, JPrim2
                  JTemp = JTemp + 1
                  ExpntB(JTemp) = Expnt(J)
               END DO
!
               ABCentX = ACentX - BCentX
               ABCentY = ACentY - BCentY
               ABCentZ = ACentZ - BCentZ
               R2AB = ABCentX * ABCentX + ABCentY * ABCentY + ABCentZ * ABCentZ
               IJPrim = 0
               DO I = IPrim1, IPrim2
                  ExpA = ExpntA(I-IPrim1+1)
                  ExpAR2 = ExpA * R2AB
                  JPTemp2 = JPrim2
                  IF (II == JJ) JPTemp2 = I
                  DO J = JPrim1, JPTemp2
                     ExpB = ExpntB(J-JPrim1+1)
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
!
!              o Construct three-center two-electron ERI matrix
!
               IF (Spherical) THEN
                  CALL RIMP2_RIInt2_MDInt2_ERI3c_Sph(II, JJ, KK, NPIJ, NPKL, ERIMat)
               ELSE
                  CALL RIMP2_RIInt2_MDInt2_ERI3c_Car(II, JJ, KK, NPIJ, NPKL, ERIMat)
               END IF
!
            END DO
         END DO
!$OMP END DO
!
!        o Deallocate memory for 3c AO integrals
!
         CALL Int2_Deallocate
         CALL Int2_ECoef_Deallocate
         CALL Int2_Int2e_Deallocate
         CALL Int2_Array_Deallocate
         DEALLOCATE(FF)
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
         DO K = 1, NK
            KBF_RI = KBF_RI + 1
            LabelK = LabelK + 1
            IdxBF_RI_MyRank(KBF_RI) = LabelK
!
!           o 1/3 integral transformation
!              (pq|c) -> (iq|c)
!
            WTimeBgn = MPI_WTIME()
            CALL CPU_TIME(TimeBgn)
            CALL DGEMM('T', 'N', NActO(1), NBF, NBF, One, CMO(1,(NFrzO(1)+1),1), NBF, &
     &         ERIMat(1,K), NBF, Zero, T1Int, NActO(1))
            CALL CPU_TIME(TimeEnd)
            WTimeEnd = MPI_WTIME()
            Time_T1 = Time_T1 + TimeEnd - TimeBgn
            WTime_T1 = WTime_T1 + WTimeEnd - WTimeBgn
!
!           o 2/3 integral transformation
!             (iq|c) -> (ia|c)
!
            WTimeBgn = MPI_WTIME()
            CALL CPU_TIME(TimeBgn)
            CALL DGEMM('N', 'N', NActO(1), NActV(1), NBF, One, T1Int, NActO(1), &
     &         CMO(1,(NFrzO(1)+NActO(1)+1),1), NBF, Zero, T2Int, NActO(1))
            CALL CPU_TIME(TimeEnd)
            WTimeEnd = MPI_WTIME()
            Time_T2 = Time_T2 + TimeEnd - TimeBgn
            WTime_T2 = WTime_T2 + WTimeEnd - WTimeBgn
!
!           o 2/3 integral storeing
!
            WTimeBgn = MPI_WTIME()
            CALL CPU_TIME(TimeBgn)
            DO Ia = 1, NActV(1)
               WRITE(IO(1), REC=(Ia+(KBF_RI-1)*NActV(1))) (T2Int(Ii+(Ia-1)*NActO(1)),Ii=1,NActO(1))
            END DO
            CALL CPU_TIME(TimeEnd)
            WTimeEnd = MPI_WTIME()
            Time_T2W = Time_T2W + TimeEnd - TimeBgn
            WTime_T2W = WTime_T2W + WTimeEnd - WTimeBgn
!
         END DO
!
      END DO
!
!     o Close file
!
      CLOSE(IO(1))
!
!     o Deallocate memory
!
      DEALLOCATE(CMO)
      DEALLOCATE(ERIMat)
      DEALLOCATE(T1Int)
      DEALLOCATE(T2Int)
!
      IF (MyRank == 0) THEN
         PRINT '(" ..... CPU time (3c-RIInt        ) :", F12.2)', Time_T0
         PRINT '(" ..... CPU time (1/3 tran3c1     ) :", F12.2)', Time_T1
         PRINT '(" ..... CPU time (2/3 tran3c1     ) :", F12.2)', Time_T2
         PRINT '(" ..... CPU time (2/3 tran3c1 writ) :", F12.2)', Time_T2W
         PRINT '(" ..... WALL time (3c-RIInt        ) :", F12.2)', WTime_T0
         PRINT '(" ..... WALL time (1/3 tran3c1     ) :", F12.2)', WTime_T1
         PRINT '(" ..... WALL time (2/3 tran3c1     ) :", F12.2)', WTime_T2
         PRINT '(" ..... WALL time (2/3 tran3c1 writ) :", F12.2)', WTime_T2W
      END IF
!
      END SUBROUTINE
