      SUBROUTINE MDInt2_CarToSph
!
      USE MP2_Basis_Module, ONLY : SphCoef
      USE MP2_Constant_Module, ONLY : Zero
      USE Int2_Module, ONLY : IAnglA, IAnglB, IAnglC, IAnglD
      USE Int2_Array_Module, ONLY : ERI_Array, ERI_Array_Temp
!
      IMPLICIT NONE
!
      INTEGER:: NICar, NJCar, NKCar, NLCar, NIJCar, NKLCar
      INTEGER:: NISph, NJSph, NKSph, NLSph, NIJSph, NKLSph
      INTEGER:: L_IJ, L_KL
      INTEGER:: NIJ, NIJC, IJ, K, L, IJKL, KL
      INTEGER:: I, J
      INTEGER:: IC, JC, KC, LC, IJC, KLC
      INTEGER:: LCMax, JCMax, LMax, JMax
      REAL(8) :: Dum
!
      NICar = ((IAnglA + 1) * (IAnglA + 2)) / 2
      NJCar = ((IAnglB + 1) * (IAnglB + 2)) / 2
      NKCar = ((IAnglC + 1) * (IAnglC + 2)) / 2
      NLCar = ((IAnglD + 1) * (IAnglD + 2)) / 2
      NISph = IAnglA + IAnglA + 1
      NJSph = IAnglB + IAnglB + 1
      NKSph = IAnglC + IAnglC + 1
      NLSph = IAnglD + IAnglD + 1
      NIJCar = NICar * NJCar
      NKLCar = NKCar * NLCar
      NIJSph = NISph * NJSph
      NKLSph = NKSph * NLSph
      L_IJ = IAnglA * 10 + IAnglB
      L_KL = IAnglC * 10 + IAnglD
!
!     --- Cartesian -> Solid Harmonics transformation (Ket)
!
      NIJ = 0
      NIJC = 0
!
      SELECT CASE(L_KL)
      CASE (0)
         DO IJKL = 1, NIJCar * NKLCar
            ERI_Array_Temp(IJKL) = ERI_Array(IJKL)
         END DO
      CASE (1)
         DO IJ = 1, NIJCar
!PH            ERI_Array_Temp(NIJ+   1) =   ERI_Array(NIJC+   2)
!PH            ERI_Array_Temp(NIJ+   2) =   ERI_Array(NIJC+   3)
!PH            ERI_Array_Temp(NIJ+   3) =   ERI_Array(NIJC+   1)
            ERI_Array_Temp(NIJ+   1) =   ERI_Array(NIJC+   1)
            ERI_Array_Temp(NIJ+   2) =   ERI_Array(NIJC+   2)
            ERI_Array_Temp(NIJ+   3) =   ERI_Array(NIJC+   3)
            NIJ = NIJ + 3
            NIJC = NIJC + 3
         END DO
      CASE (10)
         DO IJ = 1, NIJCar
!PH            ERI_Array_Temp(NIJ+   1) =   ERI_Array(NIJC+   2)
!PH            ERI_Array_Temp(NIJ+   2) =   ERI_Array(NIJC+   3)
!PH            ERI_Array_Temp(NIJ+   3) =   ERI_Array(NIJC+   1)
            ERI_Array_Temp(NIJ+   1) =   ERI_Array(NIJC+   1)
            ERI_Array_Temp(NIJ+   2) =   ERI_Array(NIJC+   2)
            ERI_Array_Temp(NIJ+   3) =   ERI_Array(NIJC+   3)
            NIJ = NIJ + 3
            NIJC = NIJC + 3
         END DO
      CASE (11)
         DO IJ = 1, NIJCar
!PH            ERI_Array_Temp(NIJ+   1) =   ERI_Array(NIJC+   5)
!PH            ERI_Array_Temp(NIJ+   2) =   ERI_Array(NIJC+   6)
!PH            ERI_Array_Temp(NIJ+   3) =   ERI_Array(NIJC+   4)
!PH            ERI_Array_Temp(NIJ+   4) =   ERI_Array(NIJC+   8)
!PH            ERI_Array_Temp(NIJ+   5) =   ERI_Array(NIJC+   9)
!PH            ERI_Array_Temp(NIJ+   6) =   ERI_Array(NIJC+   7)
!PH            ERI_Array_Temp(NIJ+   7) =   ERI_Array(NIJC+   2)
!PH            ERI_Array_Temp(NIJ+   8) =   ERI_Array(NIJC+   3)
!PH            ERI_Array_Temp(NIJ+   9) =   ERI_Array(NIJC+   1)
            ERI_Array_Temp(NIJ+   1) =   ERI_Array(NIJC+   1)
            ERI_Array_Temp(NIJ+   2) =   ERI_Array(NIJC+   2)
            ERI_Array_Temp(NIJ+   3) =   ERI_Array(NIJC+   3)
            ERI_Array_Temp(NIJ+   4) =   ERI_Array(NIJC+   4)
            ERI_Array_Temp(NIJ+   5) =   ERI_Array(NIJC+   5)
            ERI_Array_Temp(NIJ+   6) =   ERI_Array(NIJC+   6)
            ERI_Array_Temp(NIJ+   7) =   ERI_Array(NIJC+   7)
            ERI_Array_Temp(NIJ+   8) =   ERI_Array(NIJC+   8)
            ERI_Array_Temp(NIJ+   9) =   ERI_Array(NIJC+   9)
            NIJ = NIJ + 9
            NIJC = NIJC + 9
         END DO
      CASE DEFAULT   ! Ket contains angular momentums larger than D
         DO IJ = 1, NIJCar
            KL = 0
            DO K = 1, NKSph
               LMax = NLSph
               DO L = 1, LMax
                  KL = KL + 1
                  Dum = Zero
                  KLC = 0
                  DO KC = 1, NKCar
                     LCMax = NLCar
                     DO LC = 1, LCMax
                        KLC = KLC + 1
                        Dum = Dum + SphCoef(K,KC,IAnglC) * SphCoef(L,LC,IAnglD) * ERI_Array(NIJC+KLC)
                     END DO
                  END DO
                  ERI_Array_Temp(NIJ+KL) = Dum
               END DO
            END DO
            NIJ  = NIJ  + NKLSph
            NIJC = NIJC + NKLCar
         END DO
      END SELECT   ! End Ket transformation
!
!     --- Cartesian -> Solid Harmonics transformation (Bra)
!
      SELECT CASE(L_IJ)
      CASE (0)
         DO IJKL = 1, NIJSph * NKLSph
            ERI_Array(IJKL) = ERI_Array_Temp(IJKL)
         END DO
      CASE (1)
         DO KL = 1, NKLSph
!PH            ERI_Array(NKLSph*(1-1)+KL) = ERI_Array_Temp(NKLSph*(2-1)+KL)
!PH            ERI_Array(NKLSph*(2-1)+KL) = ERI_Array_Temp(NKLSph*(3-1)+KL)
!PH            ERI_Array(NKLSph*(3-1)+KL) = ERI_Array_Temp(NKLSph*(1-1)+KL)
            ERI_Array(NKLSph*(1-1)+KL) = ERI_Array_Temp(NKLSph*(1-1)+KL)
            ERI_Array(NKLSph*(2-1)+KL) = ERI_Array_Temp(NKLSph*(2-1)+KL)
            ERI_Array(NKLSph*(3-1)+KL) = ERI_Array_Temp(NKLSph*(3-1)+KL)
         END DO
      CASE (10)
         DO KL = 1, NKLSph
!PH            ERI_Array(NKLSph*(1-1)+KL) = ERI_Array_Temp(NKLSph*(2-1)+KL)
!PH            ERI_Array(NKLSph*(2-1)+KL) = ERI_Array_Temp(NKLSph*(3-1)+KL)
!PH            ERI_Array(NKLSph*(3-1)+KL) = ERI_Array_Temp(NKLSph*(1-1)+KL)
            ERI_Array(NKLSph*(1-1)+KL) = ERI_Array_Temp(NKLSph*(1-1)+KL)
            ERI_Array(NKLSph*(2-1)+KL) = ERI_Array_Temp(NKLSph*(2-1)+KL)
            ERI_Array(NKLSph*(3-1)+KL) = ERI_Array_Temp(NKLSph*(3-1)+KL)
         END DO
      CASE (11)
         DO KL = 1, NKLSph
!PH            ERI_Array(NKLSph*(1-1)+KL) = ERI_Array_Temp(NKLSph*(5-1)+KL)
!PH            ERI_Array(NKLSph*(2-1)+KL) = ERI_Array_Temp(NKLSph*(6-1)+KL)
!PH            ERI_Array(NKLSph*(3-1)+KL) = ERI_Array_Temp(NKLSph*(4-1)+KL)
!PH            ERI_Array(NKLSph*(4-1)+KL) = ERI_Array_Temp(NKLSph*(8-1)+KL)
!PH            ERI_Array(NKLSph*(5-1)+KL) = ERI_Array_Temp(NKLSph*(9-1)+KL)
!PH            ERI_Array(NKLSph*(6-1)+KL) = ERI_Array_Temp(NKLSph*(7-1)+KL)
!PH            ERI_Array(NKLSph*(7-1)+KL) = ERI_Array_Temp(NKLSph*(2-1)+KL)
!PH            ERI_Array(NKLSph*(8-1)+KL) = ERI_Array_Temp(NKLSph*(3-1)+KL)
!PH            ERI_Array(NKLSph*(9-1)+KL) = ERI_Array_Temp(NKLSph*(1-1)+KL)
            ERI_Array(NKLSph*(1-1)+KL) = ERI_Array_Temp(NKLSph*(1-1)+KL)
            ERI_Array(NKLSph*(2-1)+KL) = ERI_Array_Temp(NKLSph*(2-1)+KL)
            ERI_Array(NKLSph*(3-1)+KL) = ERI_Array_Temp(NKLSph*(3-1)+KL)
            ERI_Array(NKLSph*(4-1)+KL) = ERI_Array_Temp(NKLSph*(4-1)+KL)
            ERI_Array(NKLSph*(5-1)+KL) = ERI_Array_Temp(NKLSph*(5-1)+KL)
            ERI_Array(NKLSph*(6-1)+KL) = ERI_Array_Temp(NKLSph*(6-1)+KL)
            ERI_Array(NKLSph*(7-1)+KL) = ERI_Array_Temp(NKLSph*(7-1)+KL)
            ERI_Array(NKLSph*(8-1)+KL) = ERI_Array_Temp(NKLSph*(8-1)+KL)
            ERI_Array(NKLSph*(9-1)+KL) = ERI_Array_Temp(NKLSph*(9-1)+KL)
         END DO
      CASE DEFAULT   ! Bra contains angular momentums larger than D
         DO KL = 1, NKLSph
            IJ = 0
            DO I = 1, NISph
               JMax = NJSph
               DO J = 1, JMax
                  IJ = IJ + 1
                  Dum = Zero
                  IJC = 0
                  DO IC = 1, NICar
                     JCMax = NJCar
                     DO JC = 1, JCMax
                        IJC = IJC + 1
                        Dum = Dum + SphCoef(I,IC,IAnglA) * SphCoef(J,JC,IAnglB) * ERI_Array_Temp(NKLSph*(IJC-1)+KL)
                     END DO
                  END DO
                  ERI_Array(NKLSph*(IJ-1)+KL) = Dum
               END DO
            END DO
         END DO
      END SELECT   ! End Bra transformation
!
      END SUBROUTINE
