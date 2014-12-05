      SUBROUTINE MP2_FacDens
!
      USE MP2_Basis_Module, ONLY : Lt, Lu, Lv, LtuvMin_Car, LtuvMax_Car, FacDen
      USE MP2_Parameter_Module, ONLY : MaxAngl
!
!     o Make density factors
!
      IMPLICIT NONE
!
      INTEGER :: IAngl, IL, It, Iu, Iv, NFactor
      REAL(8) :: FacD
      INTEGER, ALLOCATABLE :: NFacOdd(:)
!
      ALLOCATE(NFacOdd(0:MaxAngl))
!
      NFacOdd(0) = 1
      NFactor = -1
      DO IAngl = 1, MaxAngl
         NFactor = NFactor + 2
         NFacOdd(IAngl) = NFacOdd(IAngl-1) * NFactor
      END DO
!
      DO IAngl = 0, MaxAngl
         DO IL = LtuvMin_Car(IAngl), LtuvMax_Car(IAngl)
            It = Lt(IL)
            Iu = Lu(IL)
            Iv = Lv(IL)
            FacD = DBLE(NFacOdd(IAngl)) / DBLE(NFacOdd(It) * NFacOdd(Iu) * NFacOdd(Iv))
            FacDen(It,Iu,Iv) = SQRT(FacD)
         END DO
      END DO
!
      DEALLOCATE(NFacOdd)
!
      END SUBROUTINE
