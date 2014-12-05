      SUBROUTINE Int2_Basis_Deallocate
!
      USE MP2_Basis_Module, ONLY : Expnt_HRR, CCoef_HRR, KStart_HRR, KAtom_HRR, KType_HRR, KontG_HRR, &
     &   KLoc_HRR_Car, KLoc_HRR_Sph, MapShel
!
      IMPLICIT NONE
!
      DEALLOCATE(Expnt_HRR)
      DEALLOCATE(CCoef_HRR)
      DEALLOCATE(KStart_HRR)
      DEALLOCATE(KAtom_HRR)
      DEALLOCATE(KType_HRR)
      DEALLOCATE(KontG_HRR)
      DEALLOCATE(KLoc_HRR_Car)
      DEALLOCATE(KLoc_HRR_Sph)
      DEALLOCATE(MapShel)
!
      END SUBROUTINE
