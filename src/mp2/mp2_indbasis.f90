      SUBROUTINE MP2_IndBasis
!
      USE MP2_Basis_Module, ONLY : Lt, Lu, Lv, IndLtuv, NLtuv_Car, NLtuv_Sph, LtuvMin_Car, LtuvMin_Sph, LtuvMax_Car, LtuvMax_Sph
      USE MP2_Parameter_Module, ONLY : Maxtuv
!
      IMPLICIT NONE
!
      INTEGER :: IL, It, Iu, Iv, Ituv
!
      NLtuv_Car(0) = 1
      LtuvMin_Car(0) = 1
      LtuvMax_Car(0) = 1
      NLtuv_Sph(0) = 1
      LtuvMin_Sph(0) = 1
      LtuvMax_Sph(0) = 1
      DO IL = 1, Maxtuv
         NLtuv_Car(IL) = ((IL + 1) * (IL + 2)) / 2
         LtuvMin_Car(IL) = LtuvMin_Car(IL-1) + NLtuv_Car(IL-1)
         LtuvMax_Car(IL) = LtuvMax_Car(IL-1) + NLtuv_Car(IL)
         NLtuv_Sph(IL) = IL + IL + 1
         LtuvMin_Sph(IL) = LtuvMin_Sph(IL-1) + NLtuv_Sph(IL-1)
         LtuvMax_Sph(IL) = LtuvMax_Sph(IL-1) + NLtuv_Sph(IL)
      END DO
!
      Ituv = 0
      DO IL = 0, Maxtuv
         DO It = Maxtuv, 0, -1
            DO Iu = Maxtuv, 0, -1
               DO Iv = Maxtuv, 0, -1
                  IF ((It + Iu + Iv) == IL) THEN
                     Ituv = Ituv + 1
                     Lt(Ituv) = It
                     Lu(Ituv) = Iu
                     Lv(Ituv) = Iv
                     IndLtuv(It,Iu,Iv) = Ituv - 1
                  END IF
               END DO
            END DO
         END DO
      END DO
!
      END SUBROUTINE
