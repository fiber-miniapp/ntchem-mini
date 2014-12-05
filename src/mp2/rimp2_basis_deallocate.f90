      SUBROUTINE RIMP2_Basis_Deallocate
!
      USE RIMP2_Basis_Module, ONLY : Expnt_RI, CCoef_RI, KStart_RI, KAtom_RI, KType_RI, KontG_RI, KLoc_RI_Car, KLoc_RI_Sph, &
     &   KMin_RI, KMax_RI, NLtuv_RI_Car, NLtuv_RI_Sph
!
      IMPLICIT NONE
!
      DEALLOCATE(Expnt_RI)
      DEALLOCATE(CCoef_RI)
      DEALLOCATE(KStart_RI)
      DEALLOCATE(KAtom_RI)
      DEALLOCATE(KType_RI)
      DEALLOCATE(KontG_RI)
      DEALLOCATE(KLoc_RI_Car)
      DEALLOCATE(KLoc_RI_Sph)
      DEALLOCATE(KMin_RI)   ! For DRK
      DEALLOCATE(KMax_RI)   ! For DRK
      DEALLOCATE(NLtuv_RI_Car)
      DEALLOCATE(NLtuv_RI_Sph)
!
      END SUBROUTINE
