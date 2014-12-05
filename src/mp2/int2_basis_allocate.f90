      SUBROUTINE Int2_Basis_Allocate
!
      USE MP2_Basis_Module, ONLY : Expnt_HRR, CCoef_HRR, KStart_HRR, KAtom_HRR, KType_HRR, KontG_HRR, &
     &   KLoc_HRR_Car, KLoc_HRR_Sph, MapShel
      USE MP2_Parameter_Module, ONLY : MaxPrim, MaxShel
!
      IMPLICIT NONE
!
      ALLOCATE(Expnt_HRR(MaxPrim))
      ALLOCATE(CCoef_HRR(MaxPrim))
      ALLOCATE(KStart_HRR(MaxShel))
      ALLOCATE(KAtom_HRR(MaxShel))
      ALLOCATE(KType_HRR(MaxShel))
      ALLOCATE(KontG_HRR(MaxShel))
      ALLOCATE(KLoc_HRR_Car(MaxShel))
      ALLOCATE(KLoc_HRR_Sph(MaxShel))
      ALLOCATE(MapShel(MaxShel))
!
      END SUBROUTINE
