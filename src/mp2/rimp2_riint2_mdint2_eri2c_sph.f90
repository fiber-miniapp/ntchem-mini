      SUBROUTINE RIMP2_RIInt2_MDInt2_ERI2c_Sph(ISh, JSh, KSh, LSh, NPIJ, NPKL, RIInt2c)
!
      USE MP2_Module, ONLY : IOff
      USE MP2_Basis_Module, ONLY : KLoc_Sph, Lt, Lu, Lv, LtuvMin_Car, LtuvMax_Car, LtuvMin_Sph, LtuvMax_Sph, IndLtuv
      USE MP2_Constant_Module, ONLY : Zero, One
      USE MP2_Parameter_Module, ONLY : MaxSgm2, MaxLtuv
      USE Int2_Module, ONLY : CCoefAB, CCoefCD, IAnglA, IAnglB, IAnglC, IAnglD, ThrInt
!NT101811      USE RIInt2_Array_Module, ONLY : ERI_Array, IJKLG
!NT101811      USE RIInt2_ECoef_Module, ONLY : ECoefXAB, ECoefYAB, ECoefZAB, ECoefXCD, ECoefYCD, ECoefZCD
!NT101811      USE RIInt2_Int2e_Module, ONLY : R0tuv
      USE Int2_Array_Module, ONLY : ERI_Array, IJKLG
      USE Int2_ECoef_Module, ONLY : ECoefXAB, ECoefYAB, ECoefZAB, ECoefXCD, ECoefYCD, ECoefZCD
      USE Int2_Int2e_Module, ONLY : R0tuv
      USE RIMP2_Basis_Module, ONLY : KLoc_RI_Sph
!
!     o Evaluate ERIs
!
      IMPLICIT NONE
!
      INTEGER :: NPIJ, NPKL, ISh, JSh, KSh, LSh
      REAL(8) :: RIInt2c(*)
!
      LOGICAL :: Same
      INTEGER :: IL1, IL2, IL3, IL4, I1, I2, I3, I4
      INTEGER :: It1, Iu1, Iv1, It2, Iu2, Iv2, It3, Iu3, Iv3, It4, Iu4, Iv4
      INTEGER :: It12, Iu12, Iv12, It34, Iu34, Iv34, It1234, Iu1234, Iv1234
      INTEGER :: IJSh, KLSh, IJ, KL, Ituv, IAngl34, MxInd, Ind34, Ind1234
      INTEGER :: Labelp, Labelr
      INTEGER :: It1It2, Iu1Iu2, Iv1Iv2, It3It4, Iu3Iu4, Iv3Iv4
      INTEGER :: N_IJKL, N0, N4
      REAL(8) :: ERI, ERITemp, Factor, Temp
      REAL(8) :: JInt(MaxSgm2,0:MaxLtuv)
!
      IAngl34 = IAnglC + IAnglD
!
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
      IJSh = 0
      I1 = 0
      DO IL1 = LtuvMin_Car(IAnglA), LtuvMax_Car(IAnglA)
         It1 = Lt(IL1)
         Iu1 = Lu(IL1)
         Iv1 = Lv(IL1)
         I1 = I1 + 1
!
         I2 = 0
         DO IL2 = LtuvMin_Car(IAnglB), LtuvMax_Car(IAnglB)
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
               I4 = 0
               DO IL4 = LtuvMin_Car(IAnglD), LtuvMax_Car(IAnglD)
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
                  N4 = IJKLG(I3,I4,I1,I2)
                  ERI_Array(N0) = ERI
                  IF (Same) ERI_Array(N4) = ERI
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
!     o Construct two-electron Fock matrices
!
      N_IJKL = 0
      IJSh = 0
      Labelp = KLoc_RI_Sph(ISh)
      DO IL1 = LtuvMin_Sph(IAnglA), LtuvMax_Sph(IAnglA)
         Labelp = Labelp + 1
!
         DO IL2 = LtuvMin_Sph(IAnglB), LtuvMax_Sph(IAnglB)
            IJSh = IJSh + 1
!
            KLSh = 0
            Labelr = KLoc_RI_Sph(KSh)
            DO IL3 = LtuvMin_Sph(IAnglC), LtuvMax_Sph(IAnglC)
               Labelr = Labelr + 1
!
               DO IL4 = LtuvMin_Sph(IAnglD), LtuvMax_Sph(IAnglD)
                  KLSh = KLSh + 1
!
                  N_IJKL = N_IJKL + 1
                  IF (Same .AND. (IJSh < KLSh)) CYCLE
!
                  ERI = ERI_Array(N_IJKL)
!
!                 o Write ERIs
!
                  IF (ABS(ERI) >= ThrInt) THEN
!Debug               IF (ABS(ERI) >= 1.0D-06) WRITE(*, *) Labelp, Labelr, IOff(Labelp)+Labelr, ERI
                     RIInt2c(IOff(Labelp)+Labelr) = ERI
                  END IF
!
               END DO
            END DO
         END DO
      END DO
!
      END SUBROUTINE
