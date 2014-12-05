      SUBROUTINE MP2_Basis_Allocate
!
      USE MP2_Basis_Module, ONLY : Centr, Expnt, CCoef, KStart, KAtom, KType, KontG, KLoc_Car, KMin, KMax, &
     &   Lt, Lu, Lv, IndLtuv, NLtuv_Car, LtuvMin_Car, LtuvMax_Car, KLoc_Sph, NLtuv_Sph, LtuvMin_Sph, LtuvMax_Sph, FacDen, &
     &   LLShl, ULShl, LLAO_Car, ULAO_Car, LLAO, ULAO, LLChi_Car, ULChi_Car, SphCoef, Spherical
      USE MP2_Parameter_Module, ONLY : MaxAtom, MaxPrim, MaxShel, MaxAngl, Maxtuv, MaxLtuv
!
      IMPLICIT NONE
!
      INTEGER :: LMax
!
      ALLOCATE(Centr(3,MaxAtom))
      ALLOCATE(Expnt(MaxPrim))
      ALLOCATE(CCoef(MaxPrim))
      ALLOCATE(KStart(MaxShel))
      ALLOCATE(KAtom(MaxShel))
      ALLOCATE(KType(MaxShel))
      ALLOCATE(KontG(MaxShel))
      ALLOCATE(KLoc_Car(MaxShel))
      ALLOCATE(KLoc_Sph(MaxShel))
      ALLOCATE(KMin(MaxShel))   ! For DRK
      ALLOCATE(KMax(MaxShel))   ! For DRK
      ALLOCATE(Lt(MaxLtuv))
      ALLOCATE(Lu(MaxLtuv))
      ALLOCATE(Lv(MaxLtuv))
      ALLOCATE(IndLtuv(0:Maxtuv,0:Maxtuv,0:Maxtuv))
      ALLOCATE(NLtuv_Car(0:Maxtuv))
      ALLOCATE(NLtuv_Sph(0:Maxtuv))
      ALLOCATE(LtuvMin_Car(0:Maxtuv))
      ALLOCATE(LtuvMin_Sph(0:Maxtuv))
      ALLOCATE(LtuvMax_Car(0:Maxtuv))
      ALLOCATE(LtuvMax_Sph(0:Maxtuv))
      ALLOCATE(FacDen(0:MaxAngl,0:MaxAngl,0:MaxAngl))
      ALLOCATE(LLShl(MaxAtom))
      ALLOCATE(ULShl(MaxAtom))
      ALLOCATE(LLAO_Car(MaxAtom))
      ALLOCATE(LLAO(MaxAtom))   ! For Sparse
      ALLOCATE(ULAO_Car(MaxAtom))
      ALLOCATE(ULAO(MaxAtom))   ! For Sparse
      ALLOCATE(LLChi_Car(MaxShel))
      ALLOCATE(ULChi_Car(MaxShel))
      IF (Spherical) THEN
         LMax = MAX(MaxAngl, 4)
         ALLOCATE(SphCoef(LMax*2+1,(LMax+1)*(LMax+2)/2,0:LMax))
      END IF
!
      END SUBROUTINE
