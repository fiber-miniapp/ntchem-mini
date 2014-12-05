      SUBROUTINE Int2_CCont(CCont, CCont12, ISh, JSh, IAngl, JAngl, NPrim)
!
      USE MP2_Basis_Module, ONLY : Lt, Lu, Lv, LtuvMin_Car, LtuvMax_Car, FacDen
      USE MP2_Parameter_Module, ONLY : MaxAngl, MaxSgm2
!
      IMPLICIT NONE
!
      INTEGER :: ISh, JSh, IAngl, JAngl, NPrim
      REAL(8) :: CCont(MaxSgm2)
      REAL(8) :: CCont12(MaxSgm2,0:MaxAngl,0:MaxAngl,0:MaxAngl,0:MaxAngl,0:MaxAngl,0:MaxAngl)
!
      LOGICAL :: IEqJ
      INTEGER :: IJ, IL1, IL2, IL2Temp, It1, Iu1, Iv1, It2, Iu2, Iv2
      REAL(8) :: FacDen12
!
      IEqJ = (ISh == JSh)
!
      DO IL1 = LtuvMin_Car(IAngl), LtuvMax_Car(IAngl)
         It1 = Lt(IL1)
         Iu1 = Lu(IL1)
         Iv1 = Lv(IL1)
         IL2Temp = LtuvMax_Car(JAngl)
         IF (IEqJ) IL2Temp = IL1
         DO IL2 = LtuvMin_Car(JAngl), IL2Temp
            It2 = Lt(IL2)
            Iu2 = Lu(IL2)
            Iv2 = Lv(IL2)
            FacDen12 = FacDen(It1,Iu1,Iv1) * FacDen(It2,Iu2,Iv2)
            DO IJ = 1, NPrim
               CCont12(IJ,It1,Iu1,Iv1,It2,Iu2,Iv2) = CCont(IJ) * FacDen12
            END DO
         END DO
      END DO
!
      END SUBROUTINE
