      SUBROUTINE MP2_MDInt2_PScreen_MPIOMP
!
      USE MP2_Module, ONLY : SchwInt
      USE MP2_Basis_Module, ONLY : Expnt, KontG, Spherical, CCoef, KStart, KAtom, NShel, Centr, KType
      USE MP2_Constant_Module, ONLY : Zero, Half, Pi252, One, RLN10
      USE MP2_Parameter_Module, ONLY : MaxShel
      USE Int2_Module, ONLY : ExpntA, ExpntB, ExpntC, ExpntD, ExpntP, ExpntQ, &
     &   PX, PY, PZ, QX, QY, QZ, PAX, PAY, PAZ, PBX, PBY, PBZ, QCX, QCY, QCZ, QDX, QDY, QDZ, &
     &   CContAB, CContCD, CCoefAB, CCoefCD, PreFactAB, PreFactCD, IAnglA, IAnglB, IAnglC, IAnglD, ThrPrim
      USE Int2_ECoef_Module, ONLY : ECoefXAB, ECoefYAB, ECoefZAB, ECoefXCD, ECoefYCD, ECoefZCD, ExpPHalf, ExpQHalf
      USE Int2_Int2e_Module, ONLY : PQX, PQY, PQZ, ExpntPQ1, ExpntPQ2
      USE Int2_Gamma_Module, ONLY : FF, MaxtuvGam
!
      IMPLICIT NONE
!
      INTEGER :: II, JJ, KK, LL, IJ, KL, I, J, K, L
      INTEGER :: ITemp, JTemp, KTemp, LTemp, JPTemp2, LPTemp2
      INTEGER :: IAtomA, IAtomB, IAtomC, IAtomD
      INTEGER :: IPrim1, JPrim1, KPrim1, LPrim1, IPrim2, JPrim2, KPrim2, LPrim2
      REAL(8) :: ExpA, ExpB, ExpC, ExpD, ExpP, ExpQ, ExpPI, ExpQI, ExpAR2, ExpCR2, ExpnPQ
      REAL(8) :: ACentX, ACentY, ACentZ, BCentX, BCentY, BCentZ
      REAL(8) :: CCentX, CCentY, CCentZ, DCentX, DCentY, DCentZ
      REAL(8) :: ABCentX, ABCentY, ABCentZ, CDCentX, CDCentY, CDCentZ, R2AB, R2CD
      REAL(8) :: ThrFac
!
!     o E-coefficients
!
      INTEGER :: IJPrim, KLPrim, NPIJ, NPKL
      REAL(8) :: ExpKAB, ExpKCD, ExpPQ
!
!     o Initialization
!
      ThrFac = RLN10 * LOG10(ThrPrim)
!
!     o Initialization of SchwInt
!
      CALL DCOPY(MaxShel*(MaxShel+1)/2, Zero, 0, SchwInt, 1)
!
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP&PRIVATE(II, IAtomA, ACentX, ACentY, ACentZ, I, IPrim1, IPrim2, ITemp, &
!$OMP&        JJ, IAtomB, BCentX, BCentY, BCentZ, J, JPrim1, JPrim2, JTemp, JPTemp2, &
!$OMP&        KK, IAtomC, CCentX, CCentY, CCentZ, K, KPrim1, KPrim2, KTemp, &
!$OMP&        LL, IAtomD, DCentX, DCentY, DCentZ, L, LPrim1, LPrim2, LTemp, LPTemp2, &
!$OMP&        R2AB, IJPrim, ExpA, ExpB, ExpAR2, ExpP, ExpPI, ExpKAB, NPIJ, IJ, &
!$OMP&        R2CD, KLPrim, ExpC, ExpD, ExpCR2, ExpQ, ExpQI, ExpKCD, NPKL, KL, &
!$OMP&        ABCentX, ABCentY, ABCentZ, CDCentX, CDCentY, CDCentZ, &
!$OMP&        ExpnPQ, ExpPQ)
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
      DO II = 1, NShel
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
!
         DO JJ = 1, II
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
            IF (NPIJ == 0) GO TO 100
!
!           o Normalization
!
            CALL Int2_CCont(CCoefAB, CContAB, II, JJ, IAnglA, IAnglB, NPIJ)
!
!           o Generate E-coefficients for A and B
!
            IF (IAnglA >= IAnglB) THEN
               CALL MDInt2_ECoef1(ECoefXAB, ECoefYAB, ECoefZAB, PAX, PAY, PAZ, PBX, PBY, PBZ, ExpPHalf, PreFactAB, &
     &            IAnglA, IAnglB, NPIJ)
            ELSE
               CALL MDInt2_ECoef2(ECoefXAB, ECoefYAB, ECoefZAB, PAX, PAY, PAZ, PBX, PBY, PBZ, ExpPHalf, PreFactAB, &
     &            IAnglA, IAnglB, NPIJ)
            END IF
!
            KK = II
!
            IAtomC = KAtom(KK)
            IAnglC = KType(KK)
            CCentX = Centr(1,IAtomC)
            CCentY = Centr(2,IAtomC)
            CCentZ = Centr(3,IAtomC)
            KPrim1 = KStart(KK)
            KPrim2 = KPrim1 + KontG(KK) - 1
            KTemp = 0
            DO K = KPrim1, KPrim2
               KTemp = KTemp + 1
               ExpntC(KTemp) = Expnt(K)
            END DO
!
            LL = JJ
!
            IAtomD = KAtom(LL)
            IAnglD = KType(LL)
            DCentX = Centr(1,IAtomD)
            DCentY = Centr(2,IAtomD)
            DCentZ = Centr(3,IAtomD)
            LPrim1 = KStart(LL)
            LPrim2 = LPrim1 + KontG(LL) - 1
            LTemp = 0
            DO L = LPrim1, LPrim2
               LTemp = LTemp + 1
               ExpntD(LTemp) = Expnt(L)
            END DO
!
            CDCentX = CCentX - DCentX
            CDCentY = CCentY - DCentY
            CDCentZ = CCentZ - DCentZ
            R2CD = CDCentX * CDCentX + CDCentY * CDCentY + CDCentZ * CDCentZ
            KLPrim = 0
            DO K = KPrim1, KPrim2
               ExpC = ExpntC(K-KPrim1+1)
               ExpCR2 = ExpC * R2CD
               LPTemp2 = LPrim2
               IF (KK == LL) LPTemp2 = K
               DO L = LPrim1, LPTemp2
                  ExpD = ExpntD(L-LPrim1+1)
                  ExpQ = ExpC + ExpD
                  ExpQI = One / ExpQ
                  ExpKCD = - ExpCR2 * ExpD * ExpQI
                  IF (ExpKCD >= ThrFac) THEN
                     KLPrim = KLPrim + 1
                     PreFactCD(KLPrim) = EXP(ExpKCD)
                     ExpntQ(KLPrim) = ExpQ
                     ExpQHalf(KLPrim) = Half * ExpQI
                     CCoefCD(KLPrim) = CCoef(K) * CCoef(L)
                     IF ((KK == LL) .AND. (K /= L)) CCoefCD(KLPrim) = CCoefCD(KLPrim) + CCoefCD(KLPrim)
                     QX(KLPrim) = (ExpC * CCentX + ExpD * DCentX) * ExpQI
                     QY(KLPrim) = (ExpC * CCentY + ExpD * DCentY) * ExpQI
                     QZ(KLPrim) = (ExpC * CCentZ + ExpD * DCentZ) * ExpQI
                     QCX(KLPrim) = QX(KLPrim) - CCentX
                     QCY(KLPrim) = QY(KLPrim) - CCentY
                     QCZ(KLPrim) = QZ(KLPrim) - CCentZ
                     QDX(KLPrim) = QX(KLPrim) - DCentX
                     QDY(KLPrim) = QY(KLPrim) - DCentY
                     QDZ(KLPrim) = QZ(KLPrim) - DCentZ
                  END IF
               END DO
            END DO
            NPKL = KLPrim
            IF (NPKL == 0) GO TO 200
!
!           o Normalization
!
            CALL Int2_CCont(CCoefCD, CContCD, KK, LL, IAnglC, IAnglD, NPKL)
!
!           o Generate E-coefficients for C and D
!
            IF (IAnglC >= IAnglD) THEN
               CALL MDInt2_ECoef1(ECoefXCD, ECoefYCD, ECoefZCD, QCX, QCY, QCZ, QDX, QDY, QDZ, ExpQHalf, PreFactCD, &
     &            IAnglC, IAnglD, NPKL)
            ELSE
               CALL MDInt2_ECoef2(ECoefXCD, ECoefYCD, ECoefZCD, QCX, QCY, QCZ, QDX, QDY, QDZ, ExpQHalf, PreFactCD, &
     &            IAnglC, IAnglD, NPKL)
            END IF
!
!           o Generate R-integrals
!
            DO IJ = 1, NPIJ
               DO KL = 1, NPKL
                  ExpPQ = ExpntP(IJ) * ExpntQ(KL)
                  ExpnPQ = ExpntP(IJ) + ExpntQ(KL)
                  ExpntPQ1(KL,IJ) = ExpPQ / ExpnPQ
                  ExpnPQ = ExpPQ * SQRT(ExpnPQ)
                  ExpntPQ2(KL,IJ) = Pi252 / ExpnPQ   ! Adsorption
                  PQX(KL,IJ) = PX(IJ) - QX(KL)
                  PQY(KL,IJ) = PY(IJ) - QY(KL)
                  PQZ(KL,IJ) = PZ(IJ) - QZ(KL)
               END DO
            END DO
!
            CALL MDInt2_R0tuv(NPIJ, NPKL)
!
!           o Evaluate ERIs
!
            IF (Spherical) THEN
               CALL MP2_MDInt2_ERI_PScreen_Sph(II, JJ, KK, LL, NPIJ, NPKL)
            ELSE
               CALL MP2_MDInt2_ERI_PScreen_Car(II, JJ, KK, LL, NPIJ, NPKL)
            END IF
!
  200       CONTINUE
!
  100       CONTINUE
!
         END DO
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
!
      END SUBROUTINE
