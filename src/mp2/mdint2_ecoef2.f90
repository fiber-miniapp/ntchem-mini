      SUBROUTINE MDInt2_ECoef2(ECoefX, ECoefY, ECoefZ, PXAX, PYAY, PZAZ, PXBX, PYBY, PZBZ, ExpHalf, PreFact, IAngl, JAngl, NPrim)
!
      USE MP2_Constant_Module, ONLY : One
      USE MP2_Parameter_Module, ONLY : MaxAngl, MaxAng2, MaxSgm2
!
!     o Generate E-coefficients (Hermite expansion coefficients) with the three-term recurence relation
!
      IMPLICIT NONE
!
      INTEGER :: IAngl, JAngl, NPrim
      REAL(8) :: ECoefX(MaxSgm2,0:MaxAngl,0:MaxAngl,0:MaxAng2)
      REAL(8) :: ECoefY(MaxSgm2,0:MaxAngl,0:MaxAngl,0:MaxAng2)
      REAL(8) :: ECoefZ(MaxSgm2,0:MaxAngl,0:MaxAngl,0:MaxAng2)
      REAL(8) :: PXAX(*), PXBX(*), PYAY(*), PYBY(*), PZAZ(*), PZBZ(*)
      REAL(8) :: ExpHalf(*), PreFact(*)
!
      INTEGER :: I, J, It, IJ, JTemp
!
!     o E(NPrim,0,0,0), The prefactor is incorporated.
!
      DO IJ = 1, NPrim
         ECoefX(IJ,0,0,0) = One 
         ECoefY(IJ,0,0,0) = One 
         ECoefZ(IJ,0,0,0) = PreFact(IJ)   ! Adsorption of prefactor
      END DO
!
!     o Recurence relation
!       E(I,J,t,NPrim), 0 <= t <= (I + J), IAngl >= JAngl
!
      IF (JAngl == 0) RETURN
      DO J = 1, JAngl
         JTemp = J
         IF (J == JAngl) JTemp = IAngl
         DO I = 0, JTemp
!
!           * Case of It = 0
!
            It = 0
            IF ((I + J) <= 1) THEN
               DO IJ = 1, NPrim
                  ECoefX(IJ,I,J,It) = PXBX(IJ) * ECoefX(IJ,I,J-1,It)
                  ECoefY(IJ,I,J,It) = PYBY(IJ) * ECoefY(IJ,I,J-1,It)
                  ECoefZ(IJ,I,J,It) = PZBZ(IJ) * ECoefZ(IJ,I,J-1,It)
               END DO
            ELSE
               DO IJ = 1, NPrim
                  ECoefX(IJ,I,J,It) = PXBX(IJ) * ECoefX(IJ,I,J-1,It  )   &
     &                              +            ECoefX(IJ,I,J-1,It+1)
                  ECoefY(IJ,I,J,It) = PYBY(IJ) * ECoefY(IJ,I,J-1,It  )   &
     &                              +            ECoefY(IJ,I,J-1,It+1)
                  ECoefZ(IJ,I,J,It) = PZBZ(IJ) * ECoefZ(IJ,I,J-1,It  )   &
     &                              +            ECoefZ(IJ,I,J-1,It+1)
               END DO
            END IF
!
!           * Case of It >= 1
!
            DO It = 1, (I + J)
               IF (It == (I + J)) THEN
                  DO IJ = 1, NPrim
                     ECoefX(IJ,I,J,It) = ExpHalf(IJ) * ECoefX(IJ,I,J-1,It-1)
                     ECoefY(IJ,I,J,It) = ExpHalf(IJ) * ECoefY(IJ,I,J-1,It-1)
                     ECoefZ(IJ,I,J,It) = ExpHalf(IJ) * ECoefZ(IJ,I,J-1,It-1)
                  END DO
               ELSE IF (It == (I + J - 1)) THEN
                  DO IJ = 1, NPrim
                     ECoefX(IJ,I,J,It) = ExpHalf(IJ) * ECoefX(IJ,I,J-1,It-1)   &
     &                                 + PXBX(IJ)    * ECoefX(IJ,I,J-1,It  )
                     ECoefY(IJ,I,J,It) = ExpHalf(IJ) * ECoefY(IJ,I,J-1,It-1)   &
     &                                 + PYBY(IJ)    * ECoefY(IJ,I,J-1,It  )
                     ECoefZ(IJ,I,J,It) = ExpHalf(IJ) * ECoefZ(IJ,I,J-1,It-1)   &
     &                                 + PZBZ(IJ)    * ECoefZ(IJ,I,J-1,It  )
                  END DO
               ELSE
                  DO IJ = 1, NPrim
                     ECoefX(IJ,I,J,It) = ExpHalf(IJ) * ECoefX(IJ,I,J-1,It-1)   &
     &                                 + PXBX(IJ)    * ECoefX(IJ,I,J-1,It  )   &
     &                                 + (It + 1)    * ECoefX(IJ,I,J-1,It+1)
                     ECoefY(IJ,I,J,It) = ExpHalf(IJ) * ECoefY(IJ,I,J-1,It-1)   &
     &                                 + PYBY(IJ)    * ECoefY(IJ,I,J-1,It  )   &
     &                                 + (It + 1)    * ECoefY(IJ,I,J-1,It+1)
                     ECoefZ(IJ,I,J,It) = ExpHalf(IJ) * ECoefZ(IJ,I,J-1,It-1)   &
     &                                 + PZBZ(IJ)    * ECoefZ(IJ,I,J-1,It  )   &
     &                                 + (It + 1)    * ECoefZ(IJ,I,J-1,It+1)
                  END DO
               END IF
            END DO
!
            IF (I /= J) THEN
!
!              * Case of It = 0
!
               It = 0
               IF ((I + J) <= 1) THEN
                  DO IJ = 1, NPrim
                     ECoefX(IJ,J,I,It) = PXAX(IJ) * ECoefX(IJ,J-1,I,It)
                     ECoefY(IJ,J,I,It) = PYAY(IJ) * ECoefY(IJ,J-1,I,It)
                     ECoefZ(IJ,J,I,It) = PZAZ(IJ) * ECoefZ(IJ,J-1,I,It)
                  END DO
               ELSE
                  DO IJ = 1, NPrim
                     ECoefX(IJ,J,I,It) = PXAX(IJ) * ECoefX(IJ,J-1,I,It  )   &
     &                                 +            ECoefX(IJ,J-1,I,It+1)
                     ECoefY(IJ,J,I,It) = PYAY(IJ) * ECoefY(IJ,J-1,I,It  )   &
     &                                 +            ECoefY(IJ,J-1,I,It+1)
                     ECoefZ(IJ,J,I,It) = PZAZ(IJ) * ECoefZ(IJ,J-1,I,It  )   &
     &                                 +            ECoefZ(IJ,J-1,I,It+1)
                  END DO
               END IF
!
!              * Case of It >= 1
!
               DO It = 1, (I + J)
                  IF (It == (I + J)) THEN
                     DO IJ = 1, NPrim
                        ECoefX(IJ,J,I,It) = ExpHalf(IJ) * ECoefX(IJ,J-1,I,It-1)
                        ECoefY(IJ,J,I,It) = ExpHalf(IJ) * ECoefY(IJ,J-1,I,It-1)
                        ECoefZ(IJ,J,I,It) = ExpHalf(IJ) * ECoefZ(IJ,J-1,I,It-1)
                     END DO
                  ELSE IF (It == (I + J - 1)) THEN
                     DO IJ = 1, NPrim
                        ECoefX(IJ,J,I,It) = ExpHalf(IJ) * ECoefX(IJ,J-1,I,It-1)   &
     &                                    + PXAX(IJ)    * ECoefX(IJ,J-1,I,It  )
                        ECoefY(IJ,J,I,It) = ExpHalf(IJ) * ECoefY(IJ,J-1,I,It-1)   &
     &                                    + PYAY(IJ)    * ECoefY(IJ,J-1,I,It  )
                        ECoefZ(IJ,J,I,It) = ExpHalf(IJ) * ECoefZ(IJ,J-1,I,It-1)   &
     &                                    + PZAZ(IJ)    * ECoefZ(IJ,J-1,I,It  )
                     END DO
                  ELSE
                     DO IJ = 1, NPrim
                        ECoefX(IJ,J,I,It) = ExpHalf(IJ) * ECoefX(IJ,J-1,I,It-1)   &
     &                                    + PXAX(IJ)    * ECoefX(IJ,J-1,I,It  )   &
     &                                    + (It + 1)    * ECoefX(IJ,J-1,I,It+1)
                        ECoefY(IJ,J,I,It) = ExpHalf(IJ) * ECoefY(IJ,J-1,I,It-1)   &
     &                                    + PYAY(IJ)    * ECoefY(IJ,J-1,I,It  )   &
     &                                    + (It + 1)    * ECoefY(IJ,J-1,I,It+1)
                        ECoefZ(IJ,J,I,It) = ExpHalf(IJ) * ECoefZ(IJ,J-1,I,It-1)   &
     &                                    + PZAZ(IJ)    * ECoefZ(IJ,J-1,I,It  )   &
     &                                    + (It + 1)    * ECoefZ(IJ,J-1,I,It+1)
                     END DO
                  END IF
               END DO
!
            END IF
!
         END DO
      END DO
!
      END SUBROUTINE
