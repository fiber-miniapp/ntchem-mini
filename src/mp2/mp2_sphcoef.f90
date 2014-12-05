      SUBROUTINE MP2_SphCoef
!
      USE MP2_Basis_Module, ONLY : Spherical, SphCoef
      USE MP2_Parameter_Module, ONLY : MaxAngl
      USE MP2_Constant_Module, ONLY : Zero, One
!
      IMPLICIT NONE
!
      INTEGER :: LMax, L, LL
      INTEGER :: I, M
      INTEGER :: IPX1, IPY1, IPZ1, IPY2, IPZ2, IPX2
      INTEGER :: IX, IY, IZ
      REAL(8) :: Denom, PreFac, Arg
      REAL(8) :: Temp1, Temp2
!
      LMax = 4
      IF (Spherical) LMax = MAX(MaxAngl, LMax)
!
!     o Initialization
!
      SphCoef(:,:,:) = Zero
!
!     o s coefficients
!
      SphCoef(1,1,0) = One
      IF (LMax == 0) RETURN
!
!     o p coefficients
!
      SphCoef(1,2,1) = One
      SphCoef(2,3,1) = One
      SphCoef(3,1,1) = One
      IF (LMax == 1) GO TO 100
!     
!     o d coefficients
!
!     SphCoef(1,2,2) =  Sqrt3
!     SphCoef(2,5,2) =  Sqrt3
!     SphCoef(3,1,2) = -Half
!     SphCoef(3,4,2) = -Half
!     SphCoef(3,6,2) =  One
!     SphCoef(4,3,2) =  Sqrt3
!     SphCoef(5,1,2) =  Half * Sqrt3
!     SphCoef(5,4,2) = -Half * Sqrt3
!
!     o Higher than d 
!
      DO LL = 2, LMax
         L = LL - 1
         Arg = DBLE(L + L + 1) / DBLE(L + L + 2)
         PreFac = SQRT(Arg)
!
         I = 0
         DO IX = L, 0, -1
            DO IY = (L - IX), 0, -1
               IZ = L - IX - IY
               I = I + 1
               IPX1 = I
               IPY1 = I + (L - IX + 1)   ! Increase 1 Y
               IPZ1 = I + (L - IX + 2)   ! Increase 1 Z
               Temp1 = PreFac * SphCoef(L+L+1,I,L)
               Temp2 = PreFac * SphCoef(1,I,L)
!
!              --- Diagonal recurrence ---
!
               SphCoef(LL+LL+1,IPX1,LL) = SphCoef(LL+LL+1,IPX1,LL) + Temp1
               SphCoef(LL+LL+1,IPY1,LL) = SphCoef(LL+LL+1,IPY1,LL) - Temp2
!
               SphCoef(1,IPY1,LL) = SphCoef(1,IPY1,LL) + Temp1
               SphCoef(1,IPX1,LL) = SphCoef(1,IPX1,LL) + Temp2
!
!              --- Vertical recurrence (1) ---
!
               DO M = (-LL + 1), (LL - 1)
                  Denom = SQRT(DBLE((L + M + 1) * (L - M + 1)))
                  SphCoef(LL+M+1,IPZ1,LL) = SphCoef(LL+M+1,IPZ1,LL) + (DBLE(L + L + 1) / Denom) * SphCoef(L+M+1,I,L)
               END DO
!
            END DO   ! IY
         END DO   ! IX
!
!        --- Vertical recurrence (2) ---
!
         I = 0
         DO IX = (L - 1), 0, -1
            DO IY = (L - 1 - IX), 0, -1
               IZ = L - 1 - IX - IY
               I = I + 1
               IPX2 = I
               IPY2 = I + (L - 1 - IX + 1) + (L - 1 - IX + 2)   ! Increase 2 Y
               IPZ2 = I + (L - 1 - IX + 2) + (L - 1 - IX + 3)   ! Increase 2 Z
!
               DO M = (-LL + 1), (LL - 1)
                  IF ((L + M) /= 0) THEN
                     Arg = DBLE((L + M) * (L - M)) / DBLE((L + M + 1) * (L - M + 1))
                     PreFac = SQRT(Arg)
                     Temp1 = PreFac * SphCoef(L+M,I,L-1)
                     SphCoef(LL+M+1,IPX2,LL) = SphCoef(LL+M+1,IPX2,LL) - Temp1
                     SphCoef(LL+M+1,IPY2,LL) = SphCoef(LL+M+1,IPY2,LL) - Temp1
                     SphCoef(LL+M+1,IPZ2,LL) = SphCoef(LL+M+1,IPZ2,LL) - Temp1
                  END IF
               END DO   ! M
!
            END DO   ! IY
         END DO   ! IX
!
      END DO   ! LL
!
  100 CONTINUE
!
!     o p-1, p0, p+1 -> px, py, pz
!
      SphCoef(:,:,1) = Zero
      SphCoef(1,1,1) = One
      SphCoef(2,2,1) = One
      SphCoef(3,3,1) = One
!
      END SUBROUTINE
