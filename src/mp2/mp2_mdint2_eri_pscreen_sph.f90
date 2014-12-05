      SUBROUTINE MP2_MDInt2_ERI_PScreen_Sph(ISh, JSh, KSh, LSh, NPIJ, NPKL)
!
      USE MP2_Module, ONLY : SchwInt
      USE MP2_Basis_Module, ONLY : KLoc_Sph, Lt, Lu, Lv, LtuvMin_Car, LtuvMax_Car, LtuvMin_Sph, LtuvMax_Sph, IndLtuv
      USE MP2_Constant_Module, ONLY : Zero, One
      USE MP2_Parameter_Module, ONLY : MaxSgm2, MaxLtuv
      USE Int2_Module, ONLY : CCoefAB, CCoefCD, IAnglA, IAnglB, IAnglC, IAnglD, ThrInt
      USE Int2_Array_Module, ONLY : ERI_Array, IJKLG
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
      INTEGER :: IJSh, KLSh, IJ, KL, Ituv, N_IJKL, IAngl34, MxInd, Ind34, Ind1234
      INTEGER :: I1, I2, I3, I4, N0, N1, N2, N3, N4, N5, N6, N7
      REAL(8) :: ERI, ERITemp, Factor, Temp, SchwIn
      REAL(8) :: JInt(MaxSgm2,0:MaxLtuv)
!
      IAngl34 = IAnglC + IAnglD
!
      IEqJ = (ISh == JSh)
      KEqL = (KSh == LSh)
      Same = (ISh == KSh) .AND. (JSh == LSh)
!
!    o Prepare two-electron integral index
!
      N_IJKL = 0
      I1 = 0
      DO IL1 = LtuvMin_Car(IAnglA), LtuvMax_Car(IAnglA)
         I1 = I1 + 1
         I2 = 0
         DO IL2 = LtuvMin_Car(IAnglB), LtuvMax_Car(IAnglB)
            I2 = I2 + 1
            I3 = 0
            DO IL3 = LtuvMin_Car(IAnglC), LtuvMax_Car(IAnglC)
               I3 = I3 + 1
               I4 = 0
               DO IL4 = LtuvMin_Car(IAnglD), LtuvMax_Car(IAnglD)
                  I4 = I4 + 1
                  N_IJKL = N_IJKL + 1
                  IJKLG(I1,I2,I3,I4) = N_IJKL
               END DO
            END DO
         END DO
      END DO
!
      SchwIn = Zero
      IJSh = 0
      I1 = 0
      DO IL1 = LtuvMin_Car(IAnglA), LtuvMax_Car(IAnglA)
         It1 = Lt(IL1)
         Iu1 = Lu(IL1)
         Iv1 = Lv(IL1)
         I1 = I1 + 1
!
         IL2Temp = LtuvMax_Car(IAnglB)
         IF (IEqJ) IL2Temp = IL1
         I2 = 0
         DO IL2 = LtuvMin_Car(IAnglB), IL2Temp
            IJSh = IJSh + 1
            It2 = Lt(IL2)
            Iu2 = Lu(IL2)
            Iv2 = Lv(IL2)
            I2 = I2 + 1
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
     &                       * CCoefAB(IJ)
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
                                    Ind1234 = IndLtuv(It1234,Iu1234,Iv1234)
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
            I3 = 0
            DO IL3 = LtuvMin_Car(IAnglC), LtuvMax_Car(IAnglC)
               It3 = Lt(IL3)
               Iu3 = Lu(IL3)
               Iv3 = Lv(IL3)
               I3 = I3 + 1
!
               IL4Temp = LtuvMax_Car(IAnglD)
               IF (KEqL) IL4Temp = IL3
               I4 = 0
               DO IL4 = LtuvMin_Car(IAnglD), IL4Temp
                  KLSh = KLSh + 1
                  IF (Same .AND. (IJSh < KLSh)) GO TO 100
                  It4 = Lt(IL4)
                  Iu4 = Lu(IL4)
                  Iv4 = Lv(IL4)
                  I4 = I4 + 1
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
     &                             * CCoefCD(KL)
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
!                 o Store ERI to buffer
!
                  N0 = IJKLG(I1,I2,I3,I4)
                  N1 = IJKLG(I2,I1,I3,I4)
                  N2 = IJKLG(I1,I2,I4,I3)
                  N3 = IJKLG(I2,I1,I4,I3)
                  N4 = IJKLG(I3,I4,I1,I2)
                  N5 = IJKLG(I4,I3,I1,I2)
                  N6 = IJKLG(I3,I4,I2,I1)
                  N7 = IJKLG(I4,I3,I2,I1)
                  ERI_Array(N0) = ERI
                  IF (IEqJ) ERI_Array(N1) = ERI
                  IF (KEqL) ERI_Array(N2) = ERI
                  IF (IEqJ .AND. KEqL) ERI_Array(N3) = ERI
                  IF (Same) ERI_Array(N4) = ERI
                  IF (Same .AND. IEqJ) THEN
                     ERI_Array(N5) = ERI
                     ERI_Array(N6) = ERI
                     ERI_Array(N7) = ERI
                  END IF
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
!     o Cartesian -> Spherical transformation
!
      IF (MAX(IAnglA, IAnglB, IAnglC, IAnglD) >= 1) THEN
         CALL MDInt2_CarToSph
      END IF
!
!     o Write ERIs
!
      N_IJKL = 0
      IJSh = 0
      DO IL1 = LtuvMin_Sph(IAnglA), LtuvMax_Sph(IAnglA)
!
         IL2Temp = LtuvMax_Sph(IAnglB)
         DO IL2 = LtuvMin_Sph(IAnglB), IL2Temp
            IJSh = IJSh + 1
!
            KLSh = 0
            DO IL3 = LtuvMin_Sph(IAnglC), LtuvMax_Sph(IAnglC)
!
               IL4Temp = LtuvMax_Sph(IAnglD)
               DO IL4 = LtuvMin_Sph(IAnglD), IL4Temp
                  KLSh = KLSh + 1
!
                  N_IJKL = N_IJKL + 1
                  IF ((IEqJ .AND. (IL2 > IL1)) .OR. (KEqL .AND. (IL4 > IL3)) .OR. (Same .AND. (IJSh < KLSh))) CYCLE
!
                  ERI = ERI_Array(N_IJKL)
                  SchwIn = MAX(SchwIn, ERI)
!
               END DO
            END DO
         END DO
      END DO
!
      SchwInt((ISh*(ISh-1))/2+JSh) = SQRT(SchwIn)
!
      END SUBROUTINE
