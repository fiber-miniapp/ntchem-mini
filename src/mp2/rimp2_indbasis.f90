      SUBROUTINE RIMP2_IndBasis
!
      USE RIMP2_Basis_Module, ONLY : NLtuv_RI_Car, NLtuv_RI_Sph
      USE RIMP2_Parameter_Module, ONLY : Maxtuv_RI
!
      IMPLICIT NONE
!
      INTEGER :: IL
!
      NLtuv_RI_Car(0) = 1
      NLtuv_RI_Sph(0) = 1
      DO IL = 1, Maxtuv_RI
         NLtuv_RI_Car(IL) = ((IL + 1) * (IL + 2)) / 2
         NLtuv_RI_Sph(IL) = IL + IL + 1
      END DO
!
      END SUBROUTINE
