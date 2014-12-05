      SUBROUTINE MP2_MDInt2_ERI_PScreen_Car(ISh, JSh, KSh, LSh, NPIJ, NPKL)
!
      USE MP2_Module, ONLY : SchwInt
      USE MP2_Basis_Module, ONLY : KLoc_Sph, Lt, Lu, Lv, LtuvMin_Car, LtuvMax_Car, LtuvMin_Sph, LtuvMax_Sph, IndLtuv
      USE MP2_Constant_Module, ONLY : Zero, One
      USE MP2_Parameter_Module, ONLY : MaxSgm2, MaxLtuv
      USE Int2_Module, ONLY : CContAB, CContCD, IAnglA, IAnglB, IAnglC, IAnglD, ThrInt
      USE Int2_ECoef_Module, ONLY : ECoefXAB, ECoefYAB, ECoefZAB, ECoefXCD, ECoefYCD, ECoefZCD
      USE Int2_Int2e_Module, ONLY : R0tuv
!
!     o Evaluate ERIs
!
      IMPLICIT NONE
!
      INTEGER :: ISh, JSh, KSh, LSh, NPIJ, NPKL
!
      LOGICAL :: IEqJ, KEqL, Same
      INTEGER :: IL1, IL2, IL3, IL4, IL2Temp, IL4Temp
      INTEGER :: It1, Iu1, Iv1, It2, Iu2, Iv2, It3, Iu3, Iv3, It4, Iu4, Iv4
      INTEGER :: It12, Iu12, Iv12, It34, Iu34, Iv34, It1234, Iu1234, Iv1234
      INTEGER :: It1It2, Iu1Iu2, Iv1Iv2, It3It4, Iu3Iu4, Iv3Iv4
      INTEGER :: IJSh, KLSh, IJ, KL, Ituv, IAngl34, MxInd, Ind34, Ind1234
      REAL(8) :: ERI, ERITemp, Factor, Temp, SchwIn
      REAL(8) :: JInt(MaxSgm2,0:MaxLtuv)
!
      IAngl34 = IAnglC + IAnglD
!
      IEqJ = (ISh == JSh)
      KEqL = (KSh == LSh)
      Same = (ISh == KSh) .AND. (JSh == LSh)
!
      SchwIn = Zero
      IJSh = 0
      DO IL1 = LtuvMin_Car(IAnglA), LtuvMax_Car(IAnglA)
         It1 = Lt(IL1)
         Iu1 = Lu(IL1)
         Iv1 = Lv(IL1)
!
         IL2Temp = LtuvMax_Car(IAnglB)
         IF (IEqJ) IL2Temp = IL1
         DO IL2 = LtuvMin_Car(IAnglB), IL2Temp
            IJSh = IJSh + 1
            It2 = Lt(IL2)
            Iu2 = Lu(IL2)
            Iv2 = Lv(IL2)
!
            It1It2 = It1 + It2
            Iu1Iu2 = Iu1 + Iu2
            Iv1Iv2 = Iv1 + Iv2
!
!           o Initialization of JInt
!
            MxInd = IndLtuv(0,0,IAngl34)
            DO Ituv = 0, MxInd
               DO KL = 1, NPKL
                  JInt(KL,Ituv) = Zero
               END DO
            END DO
!
!           o Evaluate J-integrals
!
            DO It12 = 0, It1It2
               DO Iu12 = 0, Iu1Iu2
                  DO Iv12 = 0, Iv1Iv2
!
                     DO IJ = 1, NPIJ
                        Temp = ECoefXAB(IJ,It1,It2,It12) * ECoefYAB(IJ,Iu1,Iu2,Iu12) * ECoefZAB(IJ,Iv1,Iv2,Iv12) &
     &                       * CContAB(IJ,It1,Iu1,Iv1,It2,Iu2,Iv2)
!
                        IF (ABS(Temp) >= ThrInt) THEN
                           DO It34 = 0, IAngl34
                              It1234 = It12 + It34
                              DO Iu34 = 0, IAngl34 - It34
                                 Iu1234 = Iu12 + Iu34
                                 DO Iv34 = 0, IAngl34 - It34 - Iu34
                                    Iv1234 = Iv12 + Iv34
!
                                    Ind34 = IndLtuv(It34,Iu34,Iv34)
                                    Ind1234  = IndLtuv(It1234,Iu1234,Iv1234)
                                    DO KL = 1, NPKL
                                       JInt(KL,Ind34) = JInt(KL,Ind34) + Temp * R0tuv(KL,IJ,Ind1234)
                                    END DO
!
                                 END DO
                              END DO
                           END DO
                        END IF
!
                     END DO
!
                  END DO
               END DO
            END DO
!
            KLSh = 0
            DO IL3 = LtuvMin_Car(IAnglC), LtuvMax_Car(IAnglC)
               It3 = Lt(IL3)
               Iu3 = Lu(IL3)
               Iv3 = Lv(IL3)
!
               IL4Temp = LtuvMax_Car(IAnglD)
               IF (KEqL) IL4Temp = IL3
               DO IL4 = LtuvMin_Car(IAnglD), IL4Temp
                  KLSh = KLSh + 1
                  IF (Same .AND. (IJSh < KLSh)) GO TO 100
                  It4 = Lt(IL4)
                  Iu4 = Lu(IL4)
                  Iv4 = Lv(IL4)
!
                  It3It4 = It3 + It4
                  Iu3Iu4 = Iu3 + Iu4
                  Iv3Iv4 = Iv3 + Iv4
!
!                 o Evaluate ERIs
!
                  ERI = Zero
                  DO It34 = 0, It3It4
                     DO Iu34 = 0, Iu3Iu4
                        DO Iv34 = 0, Iv3Iv4
!
                           Factor = (-One) ** (It34 + Iu34 + Iv34)
                           Ind34 = IndLtuv(It34,Iu34,Iv34)
                           ERITemp = Zero
                           DO KL = 1, NPKL
                              Temp = ECoefXCD(KL,It3,It4,It34) * ECoefYCD(KL,Iu3,Iu4,Iu34) * ECoefZCD(KL,Iv3,Iv4,Iv34) &
     &                             * CContCD(KL,It3,Iu3,Iv3,It4,Iu4,Iv4)
                              IF (ABS(Temp) >= ThrInt) THEN
                                 Temp = Temp * JInt(KL,Ind34)
                                 ERITemp = ERITemp + Temp
                              END IF
                           END DO
                           ERI = ERI + ERITemp * Factor
!
                        END DO
                     END DO
                  END DO
!
                  SchwIn = MAX(SchwIn, ERI)
!
               END DO
!
            END DO
!
  100       CONTINUE
!
         END DO
!
      END DO
!
      SchwInt((ISh*(ISh-1))/2+JSh) = SQRT(SchwIn)
!
      END SUBROUTINE
