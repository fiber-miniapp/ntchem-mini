      SUBROUTINE MP2_Basis_Deallocate
!
      USE MP2_Basis_Module, ONLY : Centr, Expnt, CCoef, KStart, KAtom, KType, KontG, KLoc_Car, KMin, KMax, &
     &   Lt, Lu, Lv, IndLtuv, NLtuv_Car, LtuvMin_Car, LtuvMax_Car, KLoc_Sph, NLtuv_Sph, LtuvMin_Sph, LtuvMax_Sph, FacDen, &
     &   LLShl, ULShl, LLAO_Car, ULAO_Car, LLAO, ULAO, LLChi_Car, ULChi_Car, SphCoef, Spherical
!
      IMPLICIT NONE
!
      DEALLOCATE(Centr)
      DEALLOCATE(Expnt)
      DEALLOCATE(CCoef)
      DEALLOCATE(KStart)
      DEALLOCATE(KAtom)
      DEALLOCATE(KType)
      DEALLOCATE(KontG)
      DEALLOCATE(KLoc_Car)
      DEALLOCATE(KLoc_Sph)
      DEALLOCATE(KMin)   ! For DRK
      DEALLOCATE(KMax)   ! For DRK
      DEALLOCATE(Lt)
      DEALLOCATE(Lu)
      DEALLOCATE(Lv)
      DEALLOCATE(IndLtuv)
      DEALLOCATE(NLtuv_Car)
      DEALLOCATE(NLtuv_Sph)
      DEALLOCATE(LtuvMin_Car)
      DEALLOCATE(LtuvMin_Sph)
      DEALLOCATE(LtuvMax_Car)
      DEALLOCATE(LtuvMax_Sph)
      DEALLOCATE(FacDen)
      DEALLOCATE(LLShl)
      DEALLOCATE(ULShl)
      DEALLOCATE(LLAO_Car)
      DEALLOCATE(LLAO)
      DEALLOCATE(ULAO_Car)
      DEALLOCATE(ULAO)
      DEALLOCATE(LLChi_Car)
      DEALLOCATE(ULChi_Car)
      IF (Spherical) THEN
         DEALLOCATE(SphCoef)
      END IF
!
      END SUBROUTINE
