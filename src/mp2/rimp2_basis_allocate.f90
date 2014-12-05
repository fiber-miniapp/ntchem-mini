      SUBROUTINE RIMP2_Basis_Allocate
!
      USE RIMP2_Basis_Module, ONLY : Expnt_RI, CCoef_RI, KStart_RI, KAtom_RI, KType_RI, KontG_RI, KLoc_RI_Car, KLoc_RI_Sph, &
     &   KMin_RI, KMax_RI, NLtuv_RI_Car, NLtuv_RI_Sph
      USE RIMP2_Parameter_Module, ONLY : MaxPrim_RI, MaxShel_RI, Maxtuv_RI
!
      IMPLICIT NONE
!
      ALLOCATE(Expnt_RI(MaxPrim_RI))
      ALLOCATE(CCoef_RI(MaxPrim_RI))
      ALLOCATE(KStart_RI(MaxShel_RI))
      ALLOCATE(KAtom_RI(MaxShel_RI))
      ALLOCATE(KType_RI(MaxShel_RI))
      ALLOCATE(KontG_RI(MaxShel_RI))
      ALLOCATE(KLoc_RI_Car(MaxShel_RI))
      ALLOCATE(KLoc_RI_Sph(MaxShel_RI))
      ALLOCATE(KMin_RI(MaxShel_RI))   ! For DRK
      ALLOCATE(KMax_RI(MaxShel_RI))   ! For DRK
      ALLOCATE(NLtuv_RI_Car(0:Maxtuv_RI))
      ALLOCATE(NLtuv_RI_Sph(0:Maxtuv_RI))
!
      END SUBROUTINE
