      SUBROUTINE Int2_NormBasis(IAngl, K1, K2)
!
      USE MP2_Basis_Module, ONLY : CCoef_HRR, Expnt_HRR, NormP, NormF
      USE MP2_Constant_Module, ONLY : Zero, Pi32
!
!     o Normalization of basis sets
!
      IMPLICIT NONE
!
      INTEGER :: IAngl, K1, K2
!
      REAL(8), PARAMETER :: TM10 = 1.0D-10
      INTEGER :: Ig, Jg, Kg, NFactor
      REAL(8) :: EE, EE2, Dum, FacS, FacN
!
!DebugIF (NormP) WRITE(*, 9090)
!DebugIF (NormF) WRITE(*, 9100)
!
!     IF (NormP) ... Unnormalization of the primitive functions.
!     If contraction coefficients are given in terms of normalized primitive functions,
!     change them to go with unnormalized primitives.  For D shells, the input coefficients CD
!     must be the coefficients corresponding to the normalized primitive X**2 *Exp(-A*R**2).
!
!     FacN = (2N-1)!!/(2EE)**N * (Pi/EE)**(3/2)
!
      IF (.NOT. NormP) GO TO 100
      DO Ig = K1, K2
         EE = Expnt_HRR(Ig) + Expnt_HRR(Ig)
         FacS = Pi32 / (EE * SQRT(EE))
         FacN = FacS
         IF (IAngl >= 1) THEN
            EE2 = EE + EE
            NFactor = -1
            DO Kg = 1, IAngl
               NFactor = NFactor + 2
               FacN = FacN * DBLE(NFactor) / EE2
            END DO
         END IF
         CCoef_HRR(Ig) = CCoef_HRR(Ig) / SQRT(FacN)
      END DO  
!
!     IF(NormF) Normalize the contracted basis functions.
!
  100 CONTINUE
      IF (.NOT. NormF) GO TO 200
!
      FacS = Zero
      DO Ig = K1, K2
         DO Jg = K1, Ig
            EE = Expnt_HRR(Ig) + Expnt_HRR(Jg)
            FacN = EE * SQRT(EE)
            Dum = CCoef_HRR(Ig) * CCoef_HRR(Jg) / FacN
            IF (IAngl >= 1) THEN
               EE2 = EE + EE
               NFactor = -1  
               DO Kg = 1, IAngl
                  NFactor = NFactor + 2
                  Dum = Dum * DBLE(NFactor) / EE2
               END DO
            END IF
            IF (Ig /= Jg) Dum = Dum + Dum
            FacS = FacS + Dum
         END DO
      END DO
!
      DO Ig = K1, K2
         IF (FacS < TM10) THEN
            CCoef_HRR(Ig) = Zero
         ELSE
            CCoef_HRR(Ig) = CCoef_HRR(Ig) / SQRT(FacS * Pi32)
         END IF
      END DO
!
  200 CONTINUE
!
!9090 FORMAT(1X, 'The contracted primitive functions have been unnormalized')
!9100 FORMAT(1X, 'The contracted basis functions are now normalized to unity')
!
      END SUBROUTINE
