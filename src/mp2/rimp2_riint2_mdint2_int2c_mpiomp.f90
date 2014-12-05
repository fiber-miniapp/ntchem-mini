      SUBROUTINE RIMP2_RIInt2_MDInt2_Int2c_MPIOMP
!
      USE MP2_Basis_Module, ONLY : Spherical, Centr
      USE MP2_Constant_Module, ONLY : Zero, Half, One, Pi252, RLN10
      USE MP2_Parameter_Module, ONLY : MaxSgmt
      USE RIMP2_Module, ONLY : RIInt2c, NBC_RI
      USE Int2_Module, ONLY : IAnglA, IAnglB, IAnglC, IAnglD, PreFactAB, PreFactCD, CCoefAB, CCoefCD, &
     &   PX, PY, PZ, PAX, PAY, PAZ, PBX, PBY, PBZ, QX, QY, QZ, QCX, QCY, QCZ, QDX, QDY, QDZ, &
     &   CContAB, CContCD, ExpntA, ExpntB, ExpntC, ExpntD, ExpntP, ExpntQ, ThrPrim
      USE Int2_ECoef_Module, ONLY : ECoefXAB, ECoefYAB, ECoefZAB, ECoefXCD, ECoefYCD, ECoefZCD, ExpPHalf, ExpQHalf
      USE Int2_Int2e_Module, ONLY : ExpntPQ1, ExpntPQ2, PQX, PQY, PQZ
      USE Int2_Gamma_Module, ONLY : FF, MaxtuvGam
      USE RIMP2_Basis_Module, ONLY : KAtom_RI, KType_RI, KStart_RI, KontG_RI, Expnt_RI, CCoef_RI, NShel_RI
      USE MPI_Module, ONLY : NProcs, MyRank, MPIIO, IORank, MPI_COMM_IO, MPI_COMM_MO, NProcsMO, MyRankMO, &
     &   MPI_COMM_Mat, NProcsMat, MyRankMat
!
      IMPLICIT NONE
!
      INCLUDE 'mpif.h'
!
      INTEGER :: II, JJ, KK, LL, IJ, KL, I, K, ITemp, KTemp, NPIJ, NPKL
      INTEGER :: IAtomA, IAtomB, IAtomC, IAtomD
      INTEGER :: IPrim1, KPrim1, IPrim2, KPrim2, IJPrim, KLPrim
      INTEGER :: IErr
      REAL(8) :: ExpA, ExpC, ExpP, ExpQ, ExpPI, ExpQI, ExpAR2, ExpCR2, ExpnPQ, ExpKAB, ExpKCD
      REAL(8) :: ACentX, ACentY, ACentZ, BCentX, BCentY, BCentZ
      REAL(8) :: CCentX, CCentY, CCentZ, DCentX, DCentY, DCentZ
      REAL(8) :: R2AB, R2CD
      REAL(8) :: ThrFac
      REAL(8), ALLOCATABLE :: RWork1(:)
!
      REAL(8) :: TimeBgn, TimeEnd, WTimeBgn, WTimeEnd
      REAL(8) :: Time_RIInt2c, Time_RIInt2cC, WTime_RIInt2c, WTime_RIInt2cC
!
!     o Initialization
!
      ThrFac = RLN10 * LOG10(ThrPrim)
!
      ALLOCATE(RWork1(NBC_RI))
      CALL DCOPY(NBC_RI, Zero, 0, RWork1, 1)
!
      WTimeBgn = MPI_WTIME()
      CALL CPU_TIME(TimeBgn)
!
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP&PRIVATE(II, IAtomA, ACentX, ACentY, ACentZ, I, IPrim1, IPrim2, ITemp, &
!$OMP&        JJ, IAtomB, BCentX, BCentY, BCentZ, &
!$OMP&        KK, IAtomC, CCentX, CCentY, CCentZ, K, KPrim1, KPrim2, KTemp, &
!$OMP&        LL, IAtomD, DCentX, DCentY, DCentZ, &
!$OMP&        R2AB, IJPrim, ExpA, ExpAR2, ExpP, ExpPI, ExpKAB, NPIJ, IJ, &
!$OMP&        R2CD, KLPrim, ExpC, ExpCR2, ExpQ, ExpQI, ExpKCD, NPKL, KL, &
!$OMP&        ExpnPQ)
!
!     o Allocate memory
!
      CALL Int2_Allocate
      CALL Int2_ECoef_Allocate
      CALL Int2_Int2e_Allocate
      CALL Int2_Array_Allocate
      ALLOCATE(FF(MaxtuvGam))
!
!$OMP DO SCHEDULE(DYNAMIC, 1)
      DO II = NShel_RI, 1, -1
!      DO II = 1, NShel_RI
         IAtomA = KAtom_RI(II)
         IAnglA = KType_RI(II)
         ACentX = Centr(1,IAtomA)
         ACentY = Centr(2,IAtomA)
         ACentZ = Centr(3,IAtomA)
         IPrim1 = KStart_RI(II)
         IPrim2 = IPrim1 + KontG_RI(II) - 1
         ITemp = 0
         DO I = IPrim1, IPrim2
            ITemp = ITemp + 1
            ExpntA(ITemp) = Expnt_RI(I)
         END DO
!
         JJ = 1
         IAtomB = IAtomA
         IAnglB = 0   ! s-type
         BCentX = Centr(1,IAtomB)
         BCentY = Centr(2,IAtomB)
         BCentZ = Centr(3,IAtomB)
         ExpntB(1) = Zero
!
         R2AB = Zero
         IJPrim = 0
         DO I = IPrim1, IPrim2
            ExpA = ExpntA(I-IPrim1+1)
            ExpAR2 = ExpA * R2AB
            ExpP = ExpA
            ExpPI = One / ExpP
            ExpKAB = Zero
               IJPrim = IJPrim + 1
               PreFactAB(IJPrim) = One
               ExpntP(IJPrim) = ExpP
               ExpPHalf(IJPrim) = Half * ExpPI
               CCoefAB(IJPrim) = CCoef_RI(I)
               PX(IJPrim) = ACentX
               PY(IJPrim) = ACentY
               PZ(IJPrim) = ACentZ
               PAX(IJPrim) = Zero
               PAY(IJPrim) = Zero
               PAZ(IJPrim) = Zero
               PBX(IJPrim) = Zero
               PBY(IJPrim) = Zero
               PBZ(IJPrim) = Zero
         END DO
         NPIJ = IJPrim
         IF (NPIJ == 0) GO TO 100
!
!        o Normalization
!
         CALL RIInt2_MDInt2_CCont(CCoefAB, CContAB, IAnglA, IAnglB, NPIJ)
!
!        o Generate E-coefficients for A and B
!
         CALL MDInt2_ECoef1(ECoefXAB, ECoefYAB, ECoefZAB, PAX, PAY, PAZ, PBX, PBY, PBZ, ExpPHalf, PreFactAB, &
     &      IAnglA, IAnglB, NPIJ)
!
         DO KK = 1, II
!MPI parallel
            IF (MOD(KK, NProcs) /= MyRank) GO TO 200
!            IF (MOD(KK, NProcsMO) /= MyRankMO) GO TO 200
!            IF (MOD(KK, NProcsMat) /= MyRankMat) GO TO 200
!MPI parallel
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
            LL = 1
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
               ExpQ = ExpC
               ExpQI = One / ExpQ
               ExpKCD = Zero
                  KLPrim = KLPrim + 1
                  PreFactCD(KLPrim) = One
                  ExpntQ(KLPrim) = ExpQ
                  ExpQHalf(KLPrim) = Half * ExpQI
                  CCoefCD(KLPrim) = CCoef_RI(K)
                  QX(KLPrim) = CCentX
                  QY(KLPrim) = CCentY
                  QZ(KLPrim) = CCentZ
                  QCX(KLPrim) = Zero
                  QCY(KLPrim) = Zero
                  QCZ(KLPrim) = Zero
                  QDX(KLPrim) = Zero
                  QDY(KLPrim) = Zero
                  QDZ(KLPrim) = Zero
            END DO
            NPKL = KLPrim
            IF (NPKL == 0) GO TO 200
!
!           o Normalization
!
            CALL RIInt2_MDInt2_CCont(CCoefCD, CContCD, IAnglC, IAnglD, NPKL)
!
!           o Generate E-coefficients for C and D
!
            CALL MDInt2_ECoef1(ECoefXCD, ECoefYCD, ECoefZCD, QCX, QCY, QCZ, QDX, QDY, QDZ, ExpQHalf, PreFactCD, &
     &         IAnglC, IAnglD, NPKL)
!
!           o Generate R-integrals
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
!           o Construct two-electron integrals
!
            IF (Spherical) THEN
               CALL RIMP2_RIInt2_MDInt2_ERI2c_Sph(II, JJ, KK, LL, NPIJ, NPKL, RWork1)
            ELSE
               CALL RIMP2_RIInt2_MDInt2_ERI2c_Car(II, JJ, KK, LL, NPIJ, NPKL, RWork1)
            END IF
!
  200       CONTINUE
!
         END DO
!
  100    CONTINUE
!
      END DO
!$OMP END DO
!
!     o Deallocate memory
!
      CALL Int2_Deallocate
      CALL Int2_ECoef_Deallocate
      CALL Int2_Int2e_Deallocate
      CALL Int2_Array_Deallocate
      DEALLOCATE(FF)
!$OMP END PARALLEL
      CALL CPU_TIME(TimeEnd)
      WTimeEnd = MPI_WTIME()
      Time_RIInt2c = TimeEnd - TimeBgn
      WTime_RIInt2c = TimeEnd - TimeBgn
!
      WTimeBgn = MPI_WTIME()
      CALL CPU_TIME(TimeBgn)
      CALL MPI_Allreduce(RWork1, RIInt2c, NBC_RI, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, IErr)
!      CALL MPI_Allreduce(RWork1, RIInt2c, NBC_RI, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_MO, IErr)
!      CALL MPI_Allreduce(RWork1, RIInt2c, NBC_RI, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_MAT, IErr)
      DEALLOCATE(RWork1)
      CALL CPU_TIME(TimeEnd)
      WTimeEnd = MPI_WTIME()
      Time_RIInt2cC = TimeEnd - TimeBgn
      WTime_RIInt2cC = TimeEnd - TimeBgn
!
      IF (MyRank == 0) THEN
         PRINT '(" ..... CPU time (RIInt2c calc    ) :", F12.2)', Time_RIInt2c
         PRINT '(" ..... CPU time (RIInt2c comm    ) :", F12.2)', Time_RIInt2cC
         PRINT '(" ..... WALL time (RIInt2c calc    ) :", F12.2)', WTime_RIInt2c
         PRINT '(" ..... WALL time (RIInt2c comm    ) :", F12.2)', WTime_RIInt2cC
      END IF
!
      END SUBROUTINE
