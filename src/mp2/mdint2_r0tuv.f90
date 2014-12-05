      SUBROUTINE MDInt2_R0tuv(NPIJ, NPKL)
!
      USE MP2_Basis_Module, ONLY : Lt, Lu, Lv, IndLtuv, LtuvMin_Car, LtuvMax_Car
      USE MP2_Constant_Module, ONLY : Zero, One
      USE MP2_Parameter_Module, ONLY : MaxLtuv, MaxSgm2
      USE Int2_Module, ONLY : IAnglA, IAnglB, IAnglC, IAnglD, ThrInt
      USE Int2_Int2e_Module, ONLY : R0tuv, ExpntPQ1, ExpntPQ2, PQX, PQY, PQZ
      USE Int2_Gamma_Module, ONLY : FF
!
!     o Generate R-integrals (Hermite integrals) with recurence relation
!
      IMPLICIT NONE
!
      INTEGER :: NPIJ, NPKL
!
      INTEGER :: JJ, IJ, KL, Jtuv, IL, It, Iu, Iv, Mxtuv, MxJtuv, Indtuv0, Indtuv1, Indtuv2
      REAL(8) :: R2PQ, WValue, ExpPQ2, ExpPQj
      REAL(8), ALLOCATABLE :: RJtuv(:,:,:,:)
!
      Mxtuv = IAnglA + IAnglB + IAnglC + IAnglD
      ALLOCATE(RJtuv(NPKL,NPIJ,0:Mxtuv,0:MaxLtuv-1))
!
!     o Calculate basic auxiliary integrals RJ[000]
!
      DO IJ = 1, NPIJ
         DO KL = 1, NPKL
            IF (ABS(ExpntPQ2(KL,IJ)) >= ThrInt) THEN
               R2PQ = PQX(KL,IJ) * PQX(KL,IJ) + PQY(KL,IJ) * PQY(KL,IJ) + PQZ(KL,IJ) * PQZ(KL,IJ)
               WValue = ExpntPQ1(KL,IJ) * R2PQ
               CALL Int2_FFunc(WValue, Mxtuv+1)
               ExpPQ2 = - ExpntPQ1(KL,IJ)
               ExpPQ2 = ExpPQ2 + ExpPQ2
               ExpPQJ = One
               DO JJ = 0, Mxtuv
                  RJtuv(KL,IJ,JJ,0) = ExpPQJ * FF(JJ+1) * ExpntPQ2(KL,IJ)
                  ExpPQJ = ExpPQJ * ExpPQ2
               END DO
            ELSE
               DO JJ = 0, Mxtuv
                  RJtuv(KL,IJ,JJ,0) = Zero
               END DO
            END IF
         END DO
      END DO
!
!     o Reculsive scheme for R-integrals
!
      DO Jtuv = 1, Mxtuv
         MxJtuv = Mxtuv - Jtuv
         DO IL = LtuvMin_Car(Jtuv), LtuvMax_Car(Jtuv)
            It = Lt(IL)
            Iu = Lu(IL)
            Iv = Lv(IL)
            Indtuv0 = IndLtuv(It,Iu,Iv)
!
            IF (It == 1) THEN
               Indtuv1 = IndLtuv(It-1,Iu,Iv)
               DO JJ = 0, MxJtuv
                  DO IJ = 1, NPIJ
                     DO KL = 1, NPKL
                        RJtuv(KL,IJ,JJ,Indtuv0) = PQX(KL,IJ) * RJtuv(KL,IJ,JJ+1,Indtuv1)
                     END DO
                  END DO
               END DO
               GO TO 100
            ELSE IF (It >= 2) THEN
               Indtuv1 = IndLtuv(It-1,Iu,Iv)
               Indtuv2 = IndLtuv(It-2,Iu,Iv)
               DO JJ = 0, MxJtuv
                  DO IJ = 1, NPIJ
                     DO KL = 1, NPKL
                        RJtuv(KL,IJ,JJ,Indtuv0) = PQX(KL,IJ) * RJtuv(KL,IJ,JJ+1,Indtuv1) &
     &                                          + (It-1)     * RJtuv(KL,IJ,JJ+1,Indtuv2)
                     END DO
                  END DO
               END DO
               GO TO 100
            END IF
!
            IF (Iu == 1) THEN
               Indtuv1 = IndLtuv(It,Iu-1,Iv)
               DO JJ = 0, MxJtuv
                  DO IJ = 1, NPIJ
                     DO KL = 1, NPKL
                        RJtuv(KL,IJ,JJ,Indtuv0) = PQY(KL,IJ) * RJtuv(KL,IJ,JJ+1,Indtuv1)
                     END DO
                  END DO
               END DO
               GO TO 100
            ELSE IF (Iu >= 2) THEN
               Indtuv1 = IndLtuv(It,Iu-1,Iv)
               Indtuv2 = IndLtuv(It,Iu-2,Iv)
               DO JJ = 0, MxJtuv
                  DO IJ = 1, NPIJ
                     DO KL = 1, NPKL
                        RJtuv(KL,IJ,JJ,Indtuv0) = PQY(KL,IJ) * RJtuv(KL,IJ,JJ+1,Indtuv1) &
     &                                          + (Iu-1)     * RJtuv(KL,IJ,JJ+1,Indtuv2)
                     END DO
                  END DO
               END DO
               GO TO 100
            END IF
!
            IF (Iv == 1) THEN
               Indtuv1 = IndLtuv(It,Iu,Iv-1)
               DO JJ = 0, MxJtuv
                  DO IJ = 1, NPIJ
                     DO KL = 1, NPKL
                        RJtuv(KL,IJ,JJ,Indtuv0) = PQZ(KL,IJ) * RJtuv(KL,IJ,JJ+1,Indtuv1)
                     END DO
                  END DO
               END DO
               GO TO 100
            ELSE IF (Iv >= 2) THEN
               Indtuv1 = IndLtuv(It,Iu,Iv-1)
               Indtuv2 = IndLtuv(It,Iu,Iv-2)
               DO JJ = 0, MxJtuv
                  DO IJ = 1, NPIJ
                     DO KL = 1, NPKL
                        RJtuv(KL,IJ,JJ,Indtuv0) = PQZ(KL,IJ) * RJtuv(KL,IJ,JJ+1,Indtuv1) &
     &                                          + (Iv-1)     * RJtuv(KL,IJ,JJ+1,Indtuv2)
                     END DO
                  END DO
               END DO
               GO TO 100
            END IF
!
            STOP 'Error: Stop in MDInt2_R0tuv'
!
  100       CONTINUE
         END DO
      END DO
!
!     o Store R0-integrals R0[tuv]
!
      DO Jtuv = 0, Mxtuv
         DO IL = LtuvMin_Car(Jtuv), LtuvMax_Car(Jtuv)
            It = Lt(IL)
            Iu = Lu(IL)
            Iv = Lv(IL)
            Indtuv0 = IndLtuv(It,Iu,Iv)
            DO IJ = 1, NPIJ
               DO KL = 1, NPKL
                  R0tuv(KL,IJ,Indtuv0) = RJtuv(KL,IJ,0,Indtuv0)
               END DO
            END DO
         END DO
      END DO
!
      DEALLOCATE(RJtuv)
!
      END SUBROUTINE

!
! ************************************************************
!

      SUBROUTINE MDInt2_R0tuv_test(NPIJ, NPKL)
!
      USE MP2_Basis_Module, ONLY : Lt, Lu, Lv, IndLtuv, LtuvMin_Car, LtuvMax_Car
      USE MP2_Constant_Module, ONLY : Zero, One, Pi252
      USE MP2_Parameter_Module, ONLY : MaxLtuv, MaxSgm2
      USE Int2_Module, ONLY : IAnglA, IAnglB, IAnglC, IAnglD, ThrInt, &
           ExpntP, ExpntQ, PX, PY, PZ, QX, QY, QZ
      USE Int2_Int2e_Module, ONLY : R0tuv, ExpntPQ1, ExpntPQ2, PQX, PQY, PQZ
      USE Int2_Gamma_Module, ONLY : FF
!
!     o Generate R-integrals (Hermite integrals) with recurence relation
!
      IMPLICIT NONE
!
      INTEGER :: NPIJ, NPKL
!
      INTEGER :: JJ, IJ, KL, Jtuv, IL, It, Iu, Iv, Mxtuv, MxJtuv, Indtuv0, Indtuv1, Indtuv2
      REAL(8) :: R2PQ, WValue, ExpPQ2, ExpPQj
      REAL(8) :: ExpnPQ_IJKL, ExpntPQ1_IJKL, ExpntPQ2_IJKL
      REAL(8) :: PQX_IJKL, PQY_IJKL, PQZ_IJKL
      REAL(8), ALLOCATABLE :: RJtuv(:,:,:,:)
!
      Mxtuv = IAnglA + IAnglB + IAnglC + IAnglD
      ! ALLOCATE(RJtuv(NPKL,NPIJ,0:Mxtuv,0:MaxLtuv-1))
      ALLOCATE(RJtuv(0:Mxtuv, NPKL, NPIJ, 0:MaxLtuv-1))
!
!     o Calculate basic auxiliary integrals RJ[000]
!
      DO IJ = 1, NPIJ
         DO KL = 1, NPKL

            ExpnPQ_IJKL   = ExpntP(IJ) + ExpntQ(KL)
            ExpntPQ1_IJKL = ExpntP(IJ) * ExpntQ(KL) / ExpnPQ_IJKL
            ExpnPQ_IJKL   = ExpntP(IJ) * ExpntQ(KL) * SQRT(ExpnPQ_IJKL)
            ExpntPQ2_IJKL = Pi252 / ExpnPQ_IJKL   ! Adsorption

            PQX_IJKL = PX(IJ) - QX(KL)
            PQY_IJKL = PY(IJ) - QY(KL)
            PQZ_IJKL = PZ(IJ) - QZ(KL)

            IF (ABS(ExpntPQ2_IJKL) >= ThrInt) THEN
               R2PQ = (PQX_IJKL * PQX_IJKL) + (PQY_IJKL * PQY_IJKL) + (PQZ_IJKL * PQZ_IJKL)
               WValue = ExpntPQ1_IJKL * R2PQ
               CALL Int2_FFunc(WValue, Mxtuv+1)
               ExpPQ2 = - ExpntPQ1_IJKL
               ExpPQ2 = ExpPQ2 + ExpPQ2
               ExpPQJ = One
               DO JJ = 0, Mxtuv
                  RJtuv(JJ,KL,IJ,0) = ExpPQJ * FF(JJ+1) * ExpntPQ2_IJKL
                  ExpPQJ = ExpPQJ * ExpPQ2
               END DO
            ELSE
               DO JJ = 0, Mxtuv
                  RJtuv(JJ,KL,IJ,0) = Zero
               END DO
            END IF
         END DO
      END DO
!
!     o Reculsive scheme for R-integrals
!
      DO Jtuv = 1, Mxtuv
         MxJtuv = Mxtuv - Jtuv
         DO IL = LtuvMin_Car(Jtuv), LtuvMax_Car(Jtuv)
            It = Lt(IL)
            Iu = Lu(IL)
            Iv = Lv(IL)
            Indtuv0 = IndLtuv(It,Iu,Iv)
!
            IF (It == 1) THEN
               Indtuv1 = IndLtuv(It-1,Iu,Iv)
               DO IJ = 1, NPIJ
                  DO KL = 1, NPKL
                     PQX_IJKL = PX(IJ) - QX(KL) !
                     DO JJ = 0, MxJtuv
                        RJtuv(JJ,KL,IJ,Indtuv0) = PQX_IJKL * RJtuv(JJ+1,KL,IJ,Indtuv1)
                     END DO
                  END DO
               END DO
               GO TO 100
            ELSE IF (It >= 2) THEN
               Indtuv1 = IndLtuv(It-1,Iu,Iv)
               Indtuv2 = IndLtuv(It-2,Iu,Iv)
               DO IJ = 1, NPIJ
                  DO KL = 1, NPKL
                     PQX_IJKL = PX(IJ) - QX(KL) !
                     DO JJ = 0, MxJtuv
                        RJtuv(JJ,KL,IJ,Indtuv0) = PQX_IJKL * RJtuv(JJ+1,KL,IJ,Indtuv1) &
     &                                          + (It-1)     * RJtuv(JJ+1,KL,IJ,Indtuv2)
                     END DO
                  END DO
               END DO
               GO TO 100
            END IF
!
            IF (Iu == 1) THEN
               Indtuv1 = IndLtuv(It,Iu-1,Iv)
               DO IJ = 1, NPIJ
                  DO KL = 1, NPKL
                     PQY_IJKL = PY(IJ) - QY(KL) !
                     DO JJ = 0, MxJtuv
                        RJtuv(JJ,KL,IJ,Indtuv0) = PQY_IJKL * RJtuv(JJ+1,KL,IJ,Indtuv1)
                     END DO
                  END DO
               END DO
               GO TO 100
            ELSE IF (Iu >= 2) THEN
               Indtuv1 = IndLtuv(It,Iu-1,Iv)
               Indtuv2 = IndLtuv(It,Iu-2,Iv)
               DO IJ = 1, NPIJ
                  DO KL = 1, NPKL
                     PQY_IJKL = PY(IJ) - QY(KL) !
                     DO JJ = 0, MxJtuv
                        RJtuv(JJ,KL,IJ,Indtuv0) = PQY_IJKL * RJtuv(JJ+1,KL,IJ,Indtuv1) &
     &                                          + (Iu-1)     * RJtuv(JJ+1,KL,IJ,Indtuv2)
                     END DO
                  END DO
               END DO
               GO TO 100
            END IF
!
            IF (Iv == 1) THEN
               Indtuv1 = IndLtuv(It,Iu,Iv-1)
               DO IJ = 1, NPIJ
                  DO KL = 1, NPKL
                     PQZ_IJKL = PZ(IJ) - QZ(KL) !
                     DO JJ = 0, MxJtuv
                        RJtuv(JJ,KL,IJ,Indtuv0) = PQZ_IJKL * RJtuv(JJ+1,KL,IJ,Indtuv1)
                     END DO
                  END DO
               END DO
               GO TO 100
            ELSE IF (Iv >= 2) THEN
               Indtuv1 = IndLtuv(It,Iu,Iv-1)
               Indtuv2 = IndLtuv(It,Iu,Iv-2)
               DO IJ = 1, NPIJ
                  DO KL = 1, NPKL
                     PQZ_IJKL = PZ(IJ) - QZ(KL) !
                     DO JJ = 0, MxJtuv
                        RJtuv(JJ,KL,IJ,Indtuv0) = PQZ_IJKL * RJtuv(JJ+1,KL,IJ,Indtuv1) &
     &                                          + (Iv-1)     * RJtuv(JJ+1,KL,IJ,Indtuv2)
                     END DO
                  END DO
               END DO
               GO TO 100
            END IF
!
            STOP 'Error: Stop in MDInt2_R0tuv'
!
  100       CONTINUE
         END DO
      END DO
!
!     o Store R0-integrals R0[tuv]
!
      DO Jtuv = 0, Mxtuv
         DO IL = LtuvMin_Car(Jtuv), LtuvMax_Car(Jtuv)
            It = Lt(IL)
            Iu = Lu(IL)
            Iv = Lv(IL)
            Indtuv0 = IndLtuv(It,Iu,Iv)
            DO IJ = 1, NPIJ
               DO KL = 1, NPKL
                  R0tuv(KL,IJ,Indtuv0) = RJtuv(0,KL,IJ,Indtuv0)
               END DO
            END DO
         END DO
      END DO
!
      DEALLOCATE(RJtuv)
!
      END SUBROUTINE
