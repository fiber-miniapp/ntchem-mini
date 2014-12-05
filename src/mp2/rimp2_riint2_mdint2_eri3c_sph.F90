      SUBROUTINE RIMP2_RIInt2_MDInt2_ERI3c_Sph(ISh, JSh, KSh, NPIJ, NPKL, ERI2e)
!
      USE MP2_Module, ONLY : NBF
      USE MP2_Basis_Module, ONLY : KLoc_Sph, Lt, Lu, Lv, IndLtuv, LtuvMin_Car, LtuvMin_Sph, LtuvMax_Car, LtuvMax_Sph
      USE MP2_Constant_Module, ONLY : Zero, One
      USE MP2_Parameter_Module, ONLY : MaxSgm2, MaxLtuv
      USE Int2_Module, ONLY : IAnglA, IAnglB, IAnglC, IAnglD, CCoefAB, CCoefCD, ThrInt
!NT101811      USE RIInt2_Array_Module
!NT101811      USE RIInt2_ECoef_Module
!NT101811      USE RIInt2_Int2e_Module
      USE Int2_Array_Module, ONLY : ERI_Array, IJKLG, JInt
      USE Int2_ECoef_Module, ONLY : ECoefXAB, ECoefYAB, ECoefZAB, ECoefXCD, ECoefYCD, ECoefZCD
      USE Int2_Int2e_Module, ONLY : R0tuv
      USE RIMP2_Basis_Module, ONLY : KLoc_RI_Sph
!
!     o Evaluate ERIs
!
      IMPLICIT NONE
!
      INTEGER :: NPIJ, NPKL, ISh, JSh, KSh
      REAL(8) :: ERI2e(NBF, NBF, *)
!
      LOGICAL :: IEqJ
      INTEGER :: IL1, IL2, IL3, IL4, IL2Temp, IL4Temp
      INTEGER :: It1, Iu1, Iv1, It2, Iu2, Iv2, It3, Iu3, Iv3, It4, Iu4, Iv4
      INTEGER :: It12, Iu12, Iv12, It34, Iu34, Iv34, It1234, Iu1234, Iv1234
      INTEGER :: It1It2, Iu1Iu2, Iv1Iv2, It3It4, Iu3Iu4, Iv3Iv4
      INTEGER :: IJ, KL, Ituv, Ind34, Ind1234, MxInd, IAngl34
!      INTEGER :: IJ, KL, Ituv, IJSh, KLSh, Ind34, Ind1234, MxInd, IAngl34
      INTEGER :: I1, I2, I3, I4, N_IJKL, N0, N1
      INTEGER :: Labelp, Labelq !, Labelr
      ! REAL(8) :: JInt(MaxSgm2,0:MaxLtuv)
      REAL(8) :: ERI, ERITemp, Factor, Temp
!
      INTEGER :: NR
      integer :: MinA, MinB, MinC, MinD
      integer :: MaxA, MaxB, MaxC, MaxD
      integer :: LenA, LenB, LenC, LenD
!
      IAngl34 = IAnglC + IAnglD
!
      IEqJ = (ISh == JSh)
!
!    o Prepare two-electron integral index
!
      MinA = LtuvMin_Car(IAnglA)
      MinB = LtuvMin_Car(IAnglB)
      MinC = LtuvMin_Car(IAnglC)
      MinD = LtuvMin_Car(IAnglD)
      MaxA = LtuvMax_Car(IAnglA)
      MaxB = LtuvMax_Car(IAnglB)
      MaxC = LtuvMax_Car(IAnglC)
      MaxD = LtuvMax_Car(IAnglD)
      LenA = MaxA - MinA + 1
      LenB = MaxB - MinB + 1
      LenC = MaxC - MinC + 1
      LenD = MaxD - MinD + 1

      ! write(*,'(8i4)') MinA, MaxA, MinB, MaxB, MinC, MaxC, MinD, MaxD

#if 0
      N_IJKL = 0
      DO IL1 = MinA, MaxA
         I1 = IL1 - MinA + 1
         DO IL2 = MinB, MaxB
            I2 = IL2 - MinB + 1
            DO IL3 = MinC, MaxC
               I3 = IL3 - MinC + 1
               DO IL4 = MinD, MaxD
                  I4 = IL4 - MinD + 1
                  N_IJKL = N_IJKL + 1  ! I4 + LenD*( (I3-1) + LenC*( (I2-1) + LenB*(I1-1) ) )
                  IJKLG(I1,I2,I3,I4) = N_IJKL
               END DO
            END DO
         END DO
      END DO
#endif
!
!      IJSh = 0
      DO IL1 = MinA, MaxA
         It1 = Lt(IL1)
         Iu1 = Lu(IL1)
         Iv1 = Lv(IL1)
         I1 = IL1 - MinA + 1
!
         IL2Temp = MaxB
         IF (IEqJ) IL2Temp = IL1
         DO IL2 = MinB, IL2Temp
!            IJSh = IJSh + 1
            It2 = Lt(IL2)
            Iu2 = Lu(IL2)
            Iv2 = Lv(IL2)
            I2 = IL2 - MinB + 1
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
                             * CCoefAB(IJ)
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
                     END DO

                  END DO
               END DO
            END DO
!
!            KLSh = 0
            DO IL3 = MinC, MaxC
               It3 = Lt(IL3)
               Iu3 = Lu(IL3)
               Iv3 = Lv(IL3)
               I3 = IL3 - MinC + 1
!
               DO IL4 = MinD, MaxD
!                  KLSh = KLSh + 1
                  It4 = Lt(IL4)
                  Iu4 = Lu(IL4)
                  Iv4 = Lv(IL4)
                  I4 = IL4 - MinD + 1
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
                  N0 = I4 + LenD*( (I3-1) + LenC*( (I2-1) + LenB*(I1-1) ) )
                  N1 = I4 + LenD*( (I3-1) + LenC*( (I1-1) + LenB*(I2-1) ) )
                  ERI_Array(N0) = ERI
                  IF (IEqJ) ERI_Array(N1) = ERI

               END DO
!
            END DO
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
      MinA = LtuvMin_Sph(IAnglA)
      MinB = LtuvMin_Sph(IAnglB)
      MinC = LtuvMin_Sph(IAnglC)
      MinD = LtuvMin_Sph(IAnglD)
      MaxA = LtuvMax_Sph(IAnglA)
      MaxB = LtuvMax_Sph(IAnglB)
      MaxC = LtuvMax_Sph(IAnglC)
      MaxD = LtuvMax_Sph(IAnglD)
      LenA = MaxA - MinA + 1
      LenB = MaxB - MinB + 1
      LenC = MaxC - MinC + 1
      LenD = MaxD - MinD + 1

      DO IL1 = MinA, MaxA
         I1 = IL1 - MinA + 1
         Labelp = I1 + KLoc_Sph(ISh)
!
         DO IL2 = MinB, MaxB
            I2 = IL2 - MinB + 1
            Labelq = I2 + KLoc_Sph(JSh)
!
            IF (IEqJ .AND. (IL2 > IL1)) CYCLE

            DO IL3 = MinC, MaxC
               I3 = IL3 - MinC + 1
               NR = I3
!
               DO IL4 = MinD, MaxD
                  I4 = IL4 - MinD + 1
!
                  N_IJKL = I4 + LenD*( (I3-1) + LenC*( (I2-1) + LenB*(I1-1) ) )
                  ERI = ERI_Array(N_IJKL)
               END DO
!
!              o Copy ERIs into ERI matrix
!
               ERI2e(Labelp, labelq, NR) = ERI
               ERI2e(Labelq, labelp, NR) = ERI
!
            END DO
         END DO
      END DO
!
      END SUBROUTINE
